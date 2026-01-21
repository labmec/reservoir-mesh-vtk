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
  ENone = -1 };

// ===================
// Function prototypes
// ===================

// Creates a geometric mesh using TPZGenGrid2D
TPZGeoMesh *createGeoMesh(std::string file);

// Marks pyramid elements in the geometric mesh
void MarkPyramids(TPZGeoMesh *gmesh);

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
  int64_t nelements = gmesh->NElements();
  for (int64_t el = 0; el < nelements; el++) {
    TPZGeoEl *geoel = gmesh->Element(el);
    if (!geoel)
      continue;
    if (geoel->Type() == EPiramide) {
      geoel->SetMaterialId(EMarkedPyramide);
      // cout << "Pyramidal element found with indexmaterial " << geoel->MaterialId() << " marked as " << EMarkedPyramide << endl;
    }
  }
  gmesh->BuildConnectivity();
  cout << "Pyramidal marker callled ";
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

  val2[0] = 0.;
  cmesh -> InsertMaterialObject(
    material ->CreateBC(material, ECylinder, neumanntype, val1, val2));

  val2[0] = 1.;
  cmesh -> InsertMaterialObject(
    material ->CreateBC(material, ETampa, diritype, val1, val2));
  
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  
  return cmesh;


}

// =============
// Main function
// =============

int main(int argc, char *const argv[]) {
  
  #ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
  #endif
  
  cout << "Reading mesh and creating geometric mesh..." << endl;

  // Read .msh file to create geometric mesh
  TPZGeoMesh *gmesh = createGeoMesh("Moving_Reservoir.msh");
 0,

  cout << "Creating computational mesh..." << endl;
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

  const std::string plotfile = "postproct"; 

  constexpr int vtkRes{0}; 

  TPZManVector<std::string, 2> fields = {"Pressure", "Flux"};
  
  auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);  
  vtk.Do();

  
  
  an.Solution().Print("Solution");

  // mark pyramids
  MarkPyramids(gmesh);
 
  // TPZRefPatternTools::RefinePyramids(gmesh, EMarkedPyramide, 1);
  
  // Plot gmesh
  std::ofstream out("geomesh.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  return 0;
}