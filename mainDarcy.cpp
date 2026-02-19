/*
- Where are the refinement patterns being used? Maybe internally
  in the TPZGenGrid2D class?

- How should I setup the logging system?

- TODO: Set a better wat to define the boundary conditions and
  the problem in general.
*/

#include <iostream>
#include "TPZGenGrid2D.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "TPZVTKGenerator.h"
#include "TPZNullMaterial.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZAnalyticSolution.h"
#include "pzlog.h"

// ----------------
// Global variables
// ----------------

TLaplaceExample1 gexact;

enum EnumMatIds {
  EMatId = 1,  // Material ID for the domain
  EBottom = 2, // Material ID bottom boundary
  ERight = 3,  // Material ID right boundary
  ETop = 4,    // Material ID top boundary
  ELeft = 5,    // Material ID left boundary
  // ------
  ELagrange = 6, // Material ID for flux
  EWrap = 10, // Auxiliary material ID?
  EIntLeft = 11, // Interface material ID for left boundary
  EIntRight = 12, // Interface material ID for right boundary
  EFlux = 13 // Material ID for flux boundary condition I put it ####
};

enum EnumBCType {
  EDirichlet = 0, // Dirichlet boundary condition
  ENeumann = 1,   // Neumann boundary condition
  EMixed = 2      // Mixed boundary condition
};

// -------------------
// Function prototypes
// -------------------

// Creates a geometric mesh using TPZGenGrid2D
TPZGeoMesh* createGeoMesh(const TPZManVector<int, 2> &nelDiv, 
  const TPZManVector<REAL, 2> &minX, const TPZManVector<REAL, 2> &maxX);

// Computes the diameter of a geometric element
REAL ElementDiameter(TPZGeoEl* gel);

// Computes the diameter of a mesh  
REAL MeshDiameter(TPZGeoMesh *gmesh);

// Add auxiliary interface elements to the geometric mesh
void AddAuxiliaryGeometricElements(TPZGeoMesh *gmesh);

// Creates a computational mesh for H1 approximation
TPZCompMesh* createCompMeshH1(TPZGeoMesh *gmesh, int order = 1);

// Creates a computational mesh for mixed approximation
TPZMultiphysicsCompMesh* createCompMeshMixed(TPZGeoMesh *gmesh, int order = 1);

// Creates a computational mesh for hybrid approximation
TPZMultiphysicsCompMesh* createCompMeshHybrid(TPZGeoMesh *gmesh, int order = 1);

// Computes and prints the convergence order of relevant errors
void ConvOrder(std::string method, TPZFMatrix<REAL> &errorMat, TPZVec<REAL> &hvector, int nref);

// ----
// Main
// ----

int main (int argc, char * const argv[]) {

  // Initialize the logger (How it works?)
  #ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
  #endif

  // Initializing uniform refinements for reference elements
  // (Used in TPZGenGrid2D?)
  gRefDBase.InitializeUniformRefPattern(EOned);
  gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
  gRefDBase.InitializeUniformRefPattern(ETriangle);

  const int nthreads = 0;
  const int refMax = 5;

  TPZGeoMesh* gmesh = nullptr;
  TPZCompMesh* cmesh = nullptr;

  TPZFMatrix<REAL> errorMat(refMax, 5, 0.); // To storage the errors
  TPZVec<REAL> hvector(refMax, 0.); // To storage the mesh size

  std::string FEmethod = "H1"; // {H1, Mixed, Hybrid};
  int order = 2; // Approximation order

  // Set a problem with analytic solution
  gexact.fExact = TLaplaceExample1::ESinSin; 

  // Initial geometric mesh
  gmesh = createGeoMesh({2, 2}, {0., 0.}, {1., 1.});
  // TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  // gmesh->Print();

  // Refinement loop
  for (int ref = 0; ref < refMax; ++ref) {

    // Mesh refinement
    if (gmesh && ref > 0) {
      TPZCheckGeom checkgeom(gmesh);
      checkgeom.UniformRefine(1);
    }

    // Get mesh diameter and store it
    hvector[ref] = MeshDiameter(gmesh);

    // ------ Create computational mesh ------

    if (cmesh) {
      delete cmesh; // Delete previous computational mesh (needed?)
    }

    if (FEmethod == "Mixed") {
      cmesh = createCompMeshMixed(gmesh, order);
    } else if (FEmethod == "Hybrid") {
      cmesh = createCompMeshH1(gmesh, order);
    } else { // Uses H1 by default
      cmesh = createCompMeshH1(gmesh, order);
    }

    // ------ Assemble and solve ------

    TPZLinearAnalysis an(cmesh); // Analysis object
    TPZSSpStructMatrix<STATE> matsp(cmesh); // Sparse matrix structure for assembly
    matsp.SetNumThreads(nthreads); // Number of threads for assembly
    an.SetStructuralMatrix(matsp);

    // Set direct solver
    TPZStepSolver<STATE> step;
    auto directType = (FEmethod == "H1") ? ECholesky : ELDLt;
    step.SetDirect(directType);
    an.SetSolver(step);
    an.Run();

    an.Solution().Print("Solution");
    an.Rhs().Print("Rhs");

    // ------ Compute errors ------

    TPZVec<REAL> errors(5, 0.);
    an.PostProcessError(errors, false, std::cout);

    int nerrors = (FEmethod == "H1") ? 3 : 5; // Number of computed errors

    // Store errors in the matrix
    for(int j=0; j<nerrors; j++){
        errorMat(ref, j) = errors[j];
    }
  }

  // ------ Convergence orders ------
  ConvOrder(FEmethod, errorMat, hvector, refMax);

  // ------ Plotting ------
  // (Only the most refined solution)

  const std::string plotfile = "darcyplot";  // sem o .vtk no final
  constexpr int vtkRes{0};

  TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
  auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

  vtk.Do();

  // ------ Clean up ------

  delete cmesh;
  delete gmesh;
  
  return 0;
}




// ---------
// Functions
// ---------

TPZGeoMesh* createGeoMesh(const TPZManVector<int, 2> &nelDiv, 
                          const TPZManVector<REAL, 2> &minX, 
                          const TPZManVector<REAL, 2> &maxX) {

  TPZGeoMesh *gmesh = new TPZGeoMesh;
  TPZGenGrid2D generator(nelDiv, minX, maxX);
  generator.SetElementType(MMeshType::EQuadrilateral);
  generator.Read(gmesh, EMatId); 
  generator.SetBC(gmesh, 4, EBottom);  
  generator.SetBC(gmesh, 5, ERight);
  generator.SetBC(gmesh, 6, ETop);
  generator.SetBC(gmesh, 7, ELeft);

  return gmesh;
}

TPZCompMesh* createCompMeshH1(TPZGeoMesh *gmesh, int order) {
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(order); // Polynomial order
  cmesh->SetAllCreateFunctionsContinuous(); // H1 Elements

  // Add materials (weak formulation)
  TPZDarcyFlow *mat = new TPZDarcyFlow(EMatId, gmesh->Dimension());
  mat->SetConstantPermeability(1.0); // Set constant permeability
  mat->SetForcingFunction(gexact.ForceFunc(),4);
  mat->SetExactSol(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(mat);

  // Add boundary conditions
  TPZManVector<REAL,1> val2(1,3.); // Part that goes to the RHS vector
  TPZFMatrix<REAL> val1(1,1,0.); // Part that goes to the Stiffnes matrix

  TPZBndCondT<REAL> *bcond = mat->CreateBC(mat, EBottom, ENeumann, val1, val2);
  bcond->SetForcingFunctionBC(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(bcond);

  bcond = mat->CreateBC(mat, ERight, EDirichlet, val1, val2);
  bcond->SetForcingFunctionBC(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(bcond);
  
  bcond = mat->CreateBC(mat, ETop, ENeumann, val1, val2);
  bcond->SetForcingFunctionBC(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(bcond);
  
  bcond = mat->CreateBC(mat, ELeft, EDirichlet, val1, val2);
  bcond->SetForcingFunctionBC(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(bcond);

  // Set up the computational mesh
  cmesh->AutoBuild();

  return cmesh;
}

TPZMultiphysicsCompMesh* createCompMeshMixed(TPZGeoMesh *gmesh, int order) {

  // ------ Flux atomic cmesh -------

  TPZCompMesh *cmeshFlux = new TPZCompMesh(gmesh);
  cmeshFlux->SetDimModel(gmesh->Dimension());
  cmeshFlux->SetDefaultOrder(order);
  cmeshFlux->SetAllCreateFunctionsHDiv();

  // Add materials (weak formulation)
  TPZNullMaterial<STATE> *mat = new TPZNullMaterial(EMatId, gmesh->Dimension());  
  cmeshFlux->InsertMaterialObject(mat);

  // Create boundary conditions
  TPZManVector<REAL,1> val2(1,3.); // Part that goes to the RHS vector
  TPZFMatrix<REAL> val1(1,1,0.); // Part that goes to the Stiffnes matrix

  TPZBndCondT<REAL> *bcond = mat->CreateBC(mat, EBottom, ENeumann, val1, val2);
  cmeshFlux->InsertMaterialObject(bcond);

  bcond = mat->CreateBC(mat, ERight, EDirichlet, val1, val2);
  cmeshFlux->InsertMaterialObject(bcond);
  
  bcond = mat->CreateBC(mat, ETop, ENeumann, val1, val2);
  cmeshFlux->InsertMaterialObject(bcond);
  
  bcond = mat->CreateBC(mat, ELeft, EDirichlet, val1, val2);
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
  TPZMixedDarcyFlow *matDarcy = new TPZMixedDarcyFlow(EMatId, gmesh->Dimension());  
  matDarcy->SetConstantPermeability(1.0);
  matDarcy->SetForcingFunction(gexact.ForceFunc(),4);
  matDarcy->SetExactSol(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(matDarcy);

  // Create, set and add boundary conditions
  bcond = matDarcy->CreateBC(matDarcy, EBottom, ENeumann, val1, val2);
  bcond->SetForcingFunctionBC(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(bcond);

  bcond = matDarcy->CreateBC(matDarcy, ERight, EDirichlet, val1, val2);
  bcond->SetForcingFunctionBC(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(bcond);
  
  bcond = matDarcy->CreateBC(matDarcy, ETop, ENeumann, val1, val2);
  bcond->SetForcingFunctionBC(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(bcond);
  
  bcond = matDarcy->CreateBC(matDarcy, ELeft, EDirichlet, val1, val2);
  bcond->SetForcingFunctionBC(gexact.ExactSolution(),4);
  cmesh->InsertMaterialObject(bcond);

  // Incorporate the atomic meshes into the multiphysics mesh
  TPZManVector<TPZCompMesh *,2> cmeshes(2);
  cmeshes[0] = cmeshFlux;
  cmeshes[1] = cmeshPressure;

  TPZManVector<int> active(cmeshes.size(),1);    
  cmesh->BuildMultiphysicsSpace(active, cmeshes);

  return cmesh;
}

TPZMultiphysicsCompMesh *CompMesHybrid(TPZGeoMesh *gmesh, int order) {
  
  // Copy of the geometric mesh
  TPZGeoMesh *gmeshlocal = new TPZGeoMesh(*gmesh);

  int dim = gmeshlocal->Dimension();

  // Add interface/flux elements I guess
  AddAuxiliaryGeometricElements(gmeshlocal);

  // TPZCompMesh *h1disc = CreateAtomicH1Hybrid(problem, gmeshlocal);
  // TPZCompMesh *h1flux = CreateHDivFluxes(problem, gmeshlocal);
  // h1disc->ComputeNodElCon();
  // h1flux->ComputeNodElCon();

    // {
    //     std::ofstream out1("h1hybrid.txt");
    //     h1disc->Print(out1);
    //     std::ofstream out2("fluxmesh.txt");
    //     h1flux->Print(out2);
    // }

    //TPZManVector<TPZCompMesh *,2> meshvec = {h1flux,h1disc};
    auto *mphys = new TPZMultiphysicsCompMesh(gmeshlocal);
    //InsertMaterialObjects(problem, mphys);
    // mphys->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    // mphys->BuildMultiphysicsSpace(meshvec);
    // {
    //     TPZLagrangeMultiplierCS<> *matleft = new TPZLagrangeMultiplierCS<>(problem.fLeftInterfaceMatid.first, dim-1);
    //     TPZLagrangeMultiplierCS<> *matright = new TPZLagrangeMultiplierCS<>(problem.fLeftInterfaceMatid.second, dim-1);
    //     matright->SetMultiplier(-1.);
    //     mphys->InsertMaterialObject(matleft);
    //     mphys->InsertMaterialObject(matright);
    // }

    // InsertLagrangeElements(problem, mphys);
    // TPZH1ApproxCreator create(gmeshlocal);
    // create.HybridType() = HybridizationType::EStandard;
    // create.ProbType() = ProblemType::EDarcy;
    // create.GroupAndCondenseElements(mphys);
    // mphys->ComputeNodElCon();
    // mphys->CleanUpUnconnectedNodes();
    // mphys->SaddlePermute();
    // {
    //     std::ofstream out("mphys.txt");
    //     mphys->Print(out);
    // }
    return mphys;
}

void AddAuxiliaryGeometricElements(TPZGeoMesh *gmesh) {
  int64_t nel = gmesh->NElements();
  int dim = gmesh->Dimension();
  for(int64_t el = 0; el<nel; el++) {
    TPZGeoEl *gel = gmesh->Element(el);
    if(!gel || gel->HasSubElement() || gel->Dimension() != dim) continue;
    int firstside = gel->FirstSide(dim-1);
    int lastside = gel->NSides()-1;
    for(int side = firstside; side<lastside; side++) {
      TPZGeoElSide gelside(gel,side);
      if(gelside.HasNeighbour({ETop, EBottom, ELeft, ERight})) continue; // Skip the boundaries?
      TPZGeoElBC gbc(gel,side,EWrap); // Not sure about the material ID
      int orient = gel->NormalOrientation(side);
      int interfacematid = (orient < 0) ? EIntLeft : EIntRight;
      TPZGeoElBC(gbc.CreatedElement(),interfacematid);
    }
  }
  /// create the geometric element of the Lagrange multipliers
  /// The lagrange multiplier is associated with the larger element
  nel = gmesh->NElements();
  for(int64_t el = 0; el<nel; el++) {
    TPZGeoEl *gel = gmesh->Element(el);
    if(!gel || gel->HasSubElement() || gel->MaterialId() != EWrap) continue;
    TPZGeoElSide gelside(gel);
    if(gelside.HasNeighbour(EFlux)) continue;
    if(gelside.HasLowerLevelNeighbour(EWrap)) continue;
    TPZGeoElBC(gelside.Neighbour(), EFlux); // Material ID for the Lagrange multiplier
    }
}

REAL ElementDiameter(TPZGeoEl* gel) {
    REAL maxdist = 0.;
    int nnodes = gel->NNodes();
    for (int i = 0; i < nnodes; ++i) {
        TPZManVector<REAL,3> xi(3,0.), xj(3,0.);
        gel->Node(i).GetCoordinates(xi);
        for (int j = i+1; j < nnodes; ++j) {
            gel->Node(j).GetCoordinates(xj);
            REAL dist = 0.;
            for (int d = 0; d < gel->Dimension(); ++d) {
                dist += (xi[d] - xj[d]) * (xi[d] - xj[d]);
            }
            dist = sqrt(dist);
            if (dist > maxdist) maxdist = dist;
        }
    }
    return maxdist;
}

REAL MeshDiameter(TPZGeoMesh *gmesh) {
  REAL h = 0.;
  int64_t nel = gmesh->NElements();
  for (int64_t el = 0; el < nel; ++el) {
    TPZGeoEl *gel = gmesh->Element(el);
    if (gel && gel->Dimension() == gmesh->Dimension() && !gel->HasSubElement()) {
      REAL elSize = ElementDiameter(gel);
      if (elSize > h) h = elSize;
    }
  }
  return h;
}

void ConvOrder(std::string method, TPZFMatrix<REAL> &errorMat, TPZVec<REAL> &hvector, int nref) {

const char* errnames_h1[3] = {"H1", "L2", "Energy"};
const char* errnames_mixed[5] = {"L2 pressure", "L2 flux", "L2 div flux", "L2 grad pressure", "H(div) flux"};
const char** errnames = nullptr;
int nerrors = 0;

if (method == "Mixed") {
    nerrors = 5;
    errnames = errnames_mixed;
} else { // H1 is the default
    nerrors = 3;
    errnames = errnames_h1;
}

  REAL order;

  std::cout << std::endl;
  for(int j = 0; j < nerrors; j++) {
    std::cout << errnames[j] << " error and order:" << std::endl;
    std::cout << "ref level 0 & " << std::scientific << std::setprecision(3) << errorMat(0, j)
              << " & - " << std::endl; // No order for the first refinement
    for(int i = 1; i < nref; i++) {
      order = log(errorMat(i, j) / errorMat(i - 1, j)) / log(hvector[i] / hvector[i - 1]);
      std::cout << "ref level " << i 
                << " & " << std::scientific << std::setprecision(3) << errorMat(i, j)
                << " & " << std::fixed << std::setprecision(2) << order
                << std::endl;
    }
        std::cout << std::endl;
    }
}