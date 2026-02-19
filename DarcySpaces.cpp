#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZAnalyticSolution.h"
#include "TPZGenGrid2D.h"
#include "TPZGenGrid3D.h"
#include "TPZGeoMeshTools.h"
#include "TPZGmshReader.h"
#include "TPZRefPatternDataBase.h"
#include "pzstepsolver.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "pzlog.h"
#include <filesystem>
#include <iostream>
#include <string>
#include <thread>
#include "TPZRefPatternDataBase.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"

using namespace std;

// ================
// Global variables
// ================

int gthreads = 0;
// gthreads = std::thread::hardware_concurrency(); // Number of threads for parallel execution
// matsp.SetNumThreads(gthreads);

// Exact solution
TLaplaceExample1 gexact;

// Permeability
REAL gperm = 0.1;

// Material IDs for domain and boundaries
enum EnumMatIds { 
  EDomain = 1, 
  EFarfield = 2, 
  ECylinder = 3, 
  ETampa = 4,
  ECurveTampa = 5,
  EBoundary = 6,
  EMarkedPyramide = 100, 
  ENone = -1,
  EMarkedTetraedro = 110
};

// ===================
// Function prototypes
// ===================

// Creates a geometric mesh using TPZGenGrid2D
TPZGeoMesh *createGeoMesh(std::string file);

// Marks pyramid elements in the geometric mesh
void MarkPyramids(TPZGeoMesh *gmesh);
void DividePyramids(TPZGeoMesh *gmesh);

// Marks pyramid elements in the geometric mesh
void MarkETetraedros(TPZGeoMesh *gmesh);

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
    stringtoint[2]["surface_cap_rock"] = EFarfield;
    stringtoint[2]["surface_farfield_reservoir"] = EFarfield;
    stringtoint[2]["nome_do_msh"] = EBoundary;

    stringtoint[1]["curve_wellbore"] = ENone;
    stringtoint[1]["curve_heel"] = ECurveTampa;
    stringtoint[1]["curve_toe"] = ECurveTampa;

    stringtoint[0]["point_heel"] = ENone;
    stringtoint[0]["point_toe"] = ENone;

    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(path, gmesh);
  }

  return gmesh;
}

void MarkPyramids(TPZGeoMesh *gmesh) {
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
      cout << "Pyramidal element found with indexmaterial " << geoel->MaterialId() << " marked as " << EMarkedPyramide << endl;
    }
  }
  gmesh->BuildConnectivity();
  cout << "Pyramidal marker callled ";
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
      geoel->SetRefPattern(refpatter);
      TPZManVector<TPZGeoEl *> el(0);
      geoel->Divide(el);
     
     
    }
  }
  gmesh->BuildConnectivity();
  cout << "Pyramidal marker callled ";
}

void MarkETetraedros(TPZGeoMesh *gmesh) {
  int64_t nelements = gmesh->NElements();
  for (int64_t el = 0; el < nelements; el++) {
    TPZGeoEl *geoel = gmesh->Element(el);
    if (!geoel)
      continue;
    if (geoel->Type() == ETetraedro) {
      // cout << "Tetrahedral element found with indexmaterial " << geoel->MaterialId() << " marked as " << EMarkedTetraedro << endl;
      geoel->SetMaterialId(EMarkedTetraedro);
    }
  }
  gmesh->BuildConnectivity();
  cout << "Tetrahedral marker callled ";
}

// MALHA COMPUTACIONAL

TPZCompMesh* createCompMesh(TPZGeoMesh* gmesh) {
  TPZCompMesh* cmesh = new TPZCompMesh(gmesh);

  cmesh-> SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(1);
  cmesh->SetAllCreateFunctionsContinuous();

  TPZDarcyFlow *material = new TPZDarcyFlow(EDomain, gmesh-> Dimension());
  material -> SetConstantPermeability(gperm);
  cmesh -> InsertMaterialObject(material);
  
  TPZFMatrix<REAL> val1(1,1,0.);
  TPZManVector<REAL,1> val2(1, 9.);
  
  int diritype = 0, neumanntype = 1, robinntype = 2;
  
  
  val2[0] = 3.;
    cmesh ->InsertMaterialObject(
    material ->CreateBC(material, EFarfield, diritype, val1, val2));

  val2[0] = 1.;
  cmesh -> InsertMaterialObject(
    material ->CreateBC(material, ECylinder, diritype, val1, val2));

  val2[0] = 0.;
  cmesh -> InsertMaterialObject(
    material ->CreateBC(material, ETampa, neumanntype, val1, val2));
  
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

  // Create boundary conditions
  TPZManVector<REAL,1> val2(1,3.); // Part that goes to the RHS vector
  TPZFMatrix<REAL> val1(1,1,0.); // Part that goes to the Stiffnes matrix

  TPZNullMaterial<STATE> *mat = new TPZNullMaterial(EDomain, gmesh->Dimension());  
  cmeshFlux->InsertMaterialObject(mat);

  TPZBndCondT<REAL> *bcond = mat->CreateBC(mat, ECylinder, 1, val1, val2);
  cmeshFlux->InsertMaterialObject(bcond);

  val2[0] = 1.;
  bcond = mat ->CreateBC(mat, EFarfield, 0, val1, val2);
  cmeshFlux->InsertMaterialObject(bcond);

  val2[0] = 0.;
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
  cmesh->SetDefaultOrder(1); // Needed? Wasn't it already set in the atomic meshes?
  cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;

  // Add materials (weak formulation)
  TPZMixedDarcyFlow *matDarcy = new TPZMixedDarcyFlow(EDomain, gmesh->Dimension());  
  matDarcy->SetConstantPermeability(1.0);
  matDarcy->SetForcingFunction(gexact.ForceFunc(),4);
  matDarcy->SetExactSol(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(matDarcy);

  // Create, set and add boundary conditions
  val2[0] = 3.;
  bcond = matDarcy->CreateBC(matDarcy, ECylinder, 0, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 1.;
  bcond = matDarcy->CreateBC(matDarcy, EFarfield, 0, val1, val2);
  cmesh->InsertMaterialObject(bcond);
  
  val2[0] = 0.;
  bcond = matDarcy->CreateBC(matDarcy, ETampa, 1, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  // Incorporate the atomic meshes into the multiphysics mesh
  TPZManVector<TPZCompMesh *,2> cmeshes(2);
  cmeshes[0] = cmeshFlux;
  cmeshes[1] = cmeshPressure;

  TPZManVector<int> active(cmeshes.size(),1);    
  cmesh->BuildMultiphysicsSpace(active, cmeshes);

  return cmesh;
}


// =============
// Main function
// =============

int main(int argc, char *const argv[]) {
  
  #ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
  #endif
  gRefDBase.InitializeRefPatterns();
  
  cout << "Reading mesh and creating geometric mesh..." << endl;

  // Read .msh file to create geometric mesh
  TPZGeoMesh *gmesh = createGeoMesh("Moving_Reservoircopy.msh");

  // mark pyramids
  // MarkPyramids(gmesh);
  // DividePyramids(gmesh);

  cout << "Creating computational mesh H1..." << endl;
  TPZCompMesh *cmesh = createCompMesh(gmesh);
  // cmesh->Print(std::cout);
  TPZLinearAnalysis an(cmesh);
  TPZSSpStructMatrix<STATE, TPZStructMatrixOR <STATE >> matsp(cmesh);
  matsp.SetNumThreads(gthreads);
  an.SetStructuralMatrix(matsp);
  TPZStepSolver<STATE> step;
  step.SetDirect(ECholesky);
  

  an.SetSolver(step);
  an.Run();

  cout << "Creating computational mesh mixed..." << endl;
  TPZMultiphysicsCompMesh *cmeshMixed = createCompMeshMixed(gmesh, 1);
  // cmeshMixed->Print(std::cout);
  TPZLinearAnalysis anMixed(cmeshMixed);
  TPZSSpStructMatrix<STATE> matmixed(cmeshMixed);
  matmixed.SetNumThreads(gthreads);
  anMixed.SetStructuralMatrix(matmixed);
  TPZStepSolver<STATE> stepMixed;
  stepMixed.SetDirect(ELDLt);
 
  anMixed.SetSolver(stepMixed);
  anMixed.Run();

  // Plotting
  constexpr int vtkRes{0}; 

  {
    const std::string plotfile = "postproct"; 
    TPZManVector<std::string, 2> fields = {"Pressure", "Flux"};
    auto vtk = TPZVTKGenerator(cmeshMixed, fields, plotfile, vtkRes);  
    vtk.Do();
  }

  {
    const std::string plotfile = "postproct_mixed"; 
    TPZManVector<std::string, 2> fields = {"Pressure", "Flux"};
    auto vtk = TPZVTKGenerator(cmeshMixed, fields, plotfile, vtkRes);  
    vtk.Do();
  }


  // an.Solution().Print("Solution");

  // cout<< "Computational mesh solution: "<<endl;
  
  // cmesh->Solution().Print("CompMesh Solution");
  
 
 
  // TPZRefPatternTools::RefinePyramids(gmesh, EMarkedPyramide, 1);
  
  // Plot gmesh
  // std::ofstream out("geomesh.vtk");
  // TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  return 0;
}