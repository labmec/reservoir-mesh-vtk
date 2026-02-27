#include <filesystem>
#include <iostream>
#include <string>
#include <thread>

#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZAnalyticSolution.h"
#include "TPZGeoMeshTools.h"
#include "TPZGmshReader.h"
#include "TPZRefPatternDataBase.h"
#include "pzstepsolver.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "pzlog.h"

// ================
// Global variables
// ================

int gthreads = 0;

// Permeability
REAL gperm = 0.1;

// Material IDs for domain and boundaries
enum EnumMatIds { 
  EDomain = 1, 
  EFarfield = 2, 
  ECylinder = 3, 
  ETampa = 4,
  EMarkedPyramide = 200, 
  ENone = -1,
  EError = -2
};

enum EnumBCType {
  EDirichlet = 0, // Dirichlet boundary condition
  ENeumann = 1,   // Neumann boundary condition
  EMixed = 2      // Mixed boundary condition
};

// ===================
// Function prototypes
// ===================

// Creates a geometric mesh using TPZGenGrid2D
TPZGeoMesh *createGeoMesh(std::string file);

// Check consistency of geometric mesh (optional, for debugging)
bool checkGeoMesh(TPZGeoMesh *gmesh);

// Marks pyramid elements in the geometric mesh
void MarkPyramids(TPZGeoMesh *gmesh);

// Divides pyramid elements in the geometric mesh into tetrahedra
void DividePyramids(TPZGeoMesh *gmesh);

// Marks tetrahedral elements in the geometric mesh
void MarkETetraedros(TPZGeoMesh *gmesh);

// Computational meshes (H1 and Mixed)
TPZCompMesh* createCompMesh(TPZGeoMesh* gmesh);
TPZMultiphysicsCompMesh* createCompMeshMixed(TPZGeoMesh* gmesh, int pOrder);

// =============
// Main function
// =============

int main(int argc, char *const argv[]) {
  
  #ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
  #endif
  gRefDBase.InitializeRefPatterns();
  
  std::cout << "Reading mesh and creating geometric mesh..." << std::endl;

  // Read .msh file to create geometric mesh
  TPZGeoMesh *gmesh = createGeoMesh("wellPlusReservoir.msh");

  // Check geometric mesh for connectivity issues (optional)
  bool errorFlag = checkGeoMesh(gmesh);

  {
    std::ofstream out("gmeshOG.vtk");
    std::ofstream out2("gmeshOG.txt");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    gmesh->Print(out2);
  }

  if (errorFlag) {
    std::cout << "Something is wrong withe the imported mesh!" << std::endl;
    return 1;
  }

  // Divide pyramids into tetrahedra to avoid issues with H(div) spaces
  // MarkPyramids(gmesh);
  DividePyramids(gmesh);

  {
    std::ofstream out("gmeshAfterDivision.vtk");
    std::ofstream out2("gmshAfterDivision.txt");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    gmesh->Print(out2);
  }

  // H1 solver
  std::cout << "\nCreating computational mesh H1..." << std::endl;
  TPZCompMesh *cmesh = createCompMesh(gmesh);
  TPZLinearAnalysis an(cmesh);
  TPZSSpStructMatrix<STATE, TPZStructMatrixOR <STATE >> matsp(cmesh);
  matsp.SetNumThreads(gthreads);
  an.SetStructuralMatrix(matsp);
  TPZStepSolver<STATE> step;
  step.SetDirect(ECholesky);
  an.SetSolver(step);
  an.Run();

  // Mixed solver (H(div)-L2)
  std::cout << "\nCreating computational mesh mixed..." << std::endl;
  TPZMultiphysicsCompMesh *cmeshMixed = createCompMeshMixed(gmesh, 1);
  TPZLinearAnalysis anMixed(cmeshMixed);
  TPZSSpStructMatrix<STATE> matMixed(cmeshMixed);
  matMixed.SetNumThreads(gthreads);
  anMixed.SetStructuralMatrix(matMixed);
  TPZStepSolver<STATE> stepMixed;
  stepMixed.SetDirect(ELDLt);
  anMixed.SetSolver(stepMixed);
  anMixed.Run();

  // Print computational meshes
  {
    std::ofstream out("cmeshH1.txt");
    std::ofstream out2("cmeshHdiv.txt");
    cmesh->Print(out);
    cmeshMixed->Print(out2);
  }

  // Generate VTK files for visualization
  {
    constexpr int vtkRes{0}; 
    const std::string plotfile = "H1_solution"; 
    TPZManVector<std::string, 2> fields = {"Pressure", "Flux"};
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);  
    vtk.Do();
  }

  {
    constexpr int vtkRes{0}; 
    const std::string plotfile = "Mixed_solution"; 
    TPZManVector<std::string, 2> fields = {"Pressure", "Flux"};
    auto vtk = TPZVTKGenerator(cmeshMixed, fields, plotfile, vtkRes);  
    vtk.Do();
  }

  return 0;
}

// =========
// Functions
// =========

TPZGeoMesh *createGeoMesh(std::string file) {

  std::string currentPath = std::filesystem::current_path();
  std::string fatherPath = std::filesystem::path(currentPath).parent_path();
  std::string path(fatherPath + "/gmsh/" + file);
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  {
    TPZGmshReader reader;
    TPZManVector<std::map<std::string, int>, 4> stringtoint(4);
    stringtoint[3]["volume_nearwell"] = EDomain;
    stringtoint[3]["volume_reservoir"] = EDomain;

    stringtoint[2]["surface_wellbore_cylinder"] = ECylinder;
    stringtoint[2]["surface_wellbore_heel"] = ETampa;
    stringtoint[2]["surface_wellbore_toe"] = ETampa;
    stringtoint[2]["surface_farfield"] = EFarfield;

    // Non-used physical groups (for reference)
    stringtoint[0]["point_toe"] = ENone;
    stringtoint[0]["point_wellbore"] = ENone;
    stringtoint[1]["curve_wellbore"] = ENone;
    stringtoint[1]["curve_toe"] = ENone;
    stringtoint[1]["curve_heel"] = ENone;
    stringtoint[2]["surface_cap_rock"] = ENone;


    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(path, gmesh);
  }

  // Remove gmsh boundary elements and create GeoElBC so normals are consistent
  // int64_t nel = gmesh->NElements();
  // for (int64_t el = 0; el < nel; el++) {
  //   TPZGeoEl *gel = gmesh->Element(el);
  //   if (!gel || gel->Dimension() != gmesh->Dimension() - 1)
  //     continue;
  //   TPZGeoElSide gelside(gel);
  //   TPZGeoElSide neigh = gelside.Neighbour();
  //   gel->RemoveConnectivities();
  //   int matid = gel->MaterialId();
  //   delete gel;
  //   TPZGeoElBC gbc(neigh, matid);
  // }

  // Remove unused elements
  // Don't think it is necessary, but just in case
  // int64_t nel = gmesh->NElements();
  // for (int64_t el = 0; el < nel; el++) {
  //   TPZGeoEl *gel = gmesh->Element(el);
  //   if (!gel || gel->MaterialId() <= 4) continue;
  //   gel->RemoveConnectivities();
  //   delete gel;
  // }

  return gmesh;
}

bool checkGeoMesh(TPZGeoMesh *gmesh) {
  bool errorFlag = false;
  int64_t NElements = gmesh->NElements();
  for (int64_t el = 0; el < NElements; el++) {
    TPZGeoEl *geoel = gmesh->Element(el);
    if (!geoel || geoel->MaterialId() != EDomain)
      continue;

    // Check number of neighbours
    int numNeighbours = 0;
    for (int side = 0; side < geoel->NSides(); side++) {
      TPZGeoElSide gelside(geoel, side);
      if (gelside.Dimension() != 2)
        continue; // Only check faces
      TPZGeoElSide neigh = gelside.Neighbour();
      while (neigh != gelside) {
        numNeighbours++;
        neigh = neigh.Neighbour();
      }

      if (numNeighbours < 1) {
        gelside.Element()->SetMaterialId(EError); // Mark element with error material id
        errorFlag = true;
      }
    }
  }
  return errorFlag;
}

void MarkPyramids(TPZGeoMesh * gmesh) {
  auto refpatter = gRefDBase.FindRefPattern("PyrTwoTets");
  if (!refpatter) {
    std::cout << "Refinement pattern for pyramids not found!" << std::endl;
    DebugStop();
  }
  int64_t nelements = gmesh->NElements();
  for (int64_t el = 0; el < nelements; el++) {
    TPZGeoEl *geoel = gmesh->Element(el);
    if (!geoel)
      continue;
    if (geoel->Type() == EPiramide) {
      geoel->SetMaterialId(EMarkedPyramide);
      std::cout << "Pyramidal element found with indexmaterial " << geoel->MaterialId() << " marked as " << EMarkedPyramide << std::endl;
    }
  }
}

void DividePyramids(TPZGeoMesh *gmesh) {
  auto refpatter = gRefDBase.FindRefPattern("PyrTwoTets");

  if (!refpatter) {
    std::cout << "Refinement pattern for pyramids not found!" << std::endl;
    DebugStop();
  }
  int64_t nelements = gmesh->NElements();
  for (int64_t el = 0; el < nelements; el++) {
    TPZGeoEl *geoel = gmesh->Element(el);
    if (!geoel)
      continue;
    if (geoel->Type() == EPiramide) {
      geoel->SetMaterialId(EMarkedPyramide); // Mark elements to help debugging
      geoel->SetRefPattern(refpatter);
      TPZManVector<TPZGeoEl *> el(0);
      geoel->Divide(el);
    }
  }
}

TPZCompMesh* createCompMesh(TPZGeoMesh* gmesh) {
  TPZCompMesh* cmesh = new TPZCompMesh(gmesh);

  cmesh-> SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(1);
  cmesh->SetAllCreateFunctionsContinuous();

  TPZDarcyFlow *mat = new TPZDarcyFlow(EDomain, gmesh-> Dimension());
  TPZDarcyFlow *matPyramid = new TPZDarcyFlow(EMarkedPyramide, gmesh-> Dimension());
  mat -> SetConstantPermeability(gperm);
  matPyramid -> SetConstantPermeability(gperm);
  cmesh -> InsertMaterialObject(mat);
  cmesh -> InsertMaterialObject(matPyramid);
  
  TPZFMatrix<REAL> val1(1,1,0.);
  TPZManVector<REAL,1> val2(1, 0.);
  
  val2[0] = 1.;
  cmesh ->InsertMaterialObject(
  mat ->CreateBC(mat, EFarfield, EDirichlet, val1, val2));

  val2[0] = 1.;
  cmesh -> InsertMaterialObject(
  mat ->CreateBC(mat, ECylinder, EDirichlet, val1, val2));

  val2[0] = 1.;
  cmesh -> InsertMaterialObject(
  mat ->CreateBC(mat, ETampa, EDirichlet, val1, val2));
  
  cmesh->AutoBuild();
  cmesh->CleanUpUnconnectedNodes();
  
  return cmesh;
}

TPZMultiphysicsCompMesh* createCompMeshMixed(TPZGeoMesh *gmesh, int order) {

  // ------ Flux atomic cmesh -------

  TPZCompMesh *cmeshFlux = new TPZCompMesh(gmesh);
  cmeshFlux->SetDimModel(gmesh->Dimension());
  cmeshFlux->SetDefaultOrder(order);
  cmeshFlux->SetAllCreateFunctionsHDiv();

  TPZNullMaterial<STATE> *mat = new TPZNullMaterial(EDomain, gmesh->Dimension()); 
  TPZNullMaterial<STATE> *matPyramid = new TPZNullMaterial(EMarkedPyramide, gmesh->Dimension());   
  cmeshFlux->InsertMaterialObject(mat);
  cmeshFlux->InsertMaterialObject(matPyramid);

  // Create boundary conditions
  TPZManVector<REAL,1> val2(1,0.); // Part that goes to the RHS vector
  TPZFMatrix<REAL> val1(1,1,0.); // Part that goes to the Stiffnes matrix
  TPZBndCondT<REAL> *bcond;

  bcond = mat->CreateBC(mat, ECylinder, 1, val1, val2);
  cmeshFlux->InsertMaterialObject(bcond);

  bcond = mat ->CreateBC(mat, EFarfield, 1, val1, val2);
  cmeshFlux->InsertMaterialObject(bcond);

  bcond = mat ->CreateBC(mat, ETampa, 1, val1, val2);
  cmeshFlux->InsertMaterialObject(bcond);

  cmeshFlux->AutoBuild();

  // ------ Pressure atomic cmesh -------

  TPZCompMesh *cmeshPressure = new TPZCompMesh(gmesh);
  cmeshPressure->SetDimModel(gmesh->Dimension());
  cmeshPressure->SetDefaultOrder(order);
  if (order < 1) {
    cmeshPressure->SetAllCreateFunctionsDiscontinuous();
  } else {
    cmeshPressure->SetAllCreateFunctionsContinuous();
    cmeshPressure->ApproxSpace().CreateDisconnectedElements(true);
  }

  // Add materials (weak formulation)
  cmeshPressure->InsertMaterialObject(mat);
  cmeshPressure->InsertMaterialObject(matPyramid);

  // Set up the computational mesh
  cmeshPressure->AutoBuild();

  int ncon = cmeshPressure->NConnects();
  const int lagLevel = 1; // Lagrange multiplier level
  for(int i=0; i<ncon; i++)
  {
      TPZConnect &newnod = cmeshPressure->ConnectVec()[i]; 
      newnod.SetLagrangeMultiplier(lagLevel);
  }

  // ------ Multiphysics mesh -------

  TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(gmesh);
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(order);
  cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;

  // Add materials (weak formulation)
  TPZMixedDarcyFlow *matDarcy = new TPZMixedDarcyFlow(EDomain, gmesh->Dimension());  
  TPZMixedDarcyFlow *matDarcyPyramid = new TPZMixedDarcyFlow(EMarkedPyramide, gmesh->Dimension());
  matDarcy->SetConstantPermeability(1.0);
  matDarcyPyramid->SetConstantPermeability(1.0);
  cmesh->InsertMaterialObject(matDarcy);
  cmesh->InsertMaterialObject(matDarcyPyramid);

  // Create, set and add boundary conditions
  val2[0] = 3.;
  bcond = matDarcy->CreateBC(matDarcy, EFarfield, EDirichlet, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 1.;
  bcond = matDarcy->CreateBC(matDarcy, ECylinder, EDirichlet, val1, val2);
  cmesh->InsertMaterialObject(bcond);
  
  val2[0] = 1.;
  bcond = matDarcy->CreateBC(matDarcy, ETampa, EDirichlet, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  // Incorporate the atomic meshes into the multiphysics mesh
  TPZManVector<TPZCompMesh *,2> cmeshes(2);
  cmeshes[0] = cmeshFlux;
  cmeshes[1] = cmeshPressure;

  TPZManVector<int> active(cmeshes.size(),1);    
  cmesh->BuildMultiphysicsSpace(active, cmeshes);

  return cmesh;
}