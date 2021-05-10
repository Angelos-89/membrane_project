#define MAIN
#include <string.h>
#include <mpi.h>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <random>
#include "McmemLib.hpp"
#include "fft.h"

double PI = 4.*atan(1.0);

int main(int argc, char* argv[]){

  /*------ Initialize MPI. -------- */

  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  /*----------------- Define strings needed for filenames. -------------------*/

  std::string inputFilename  = "input_"        + std::to_string(rank) + ".txt";
  std::string outputFilename = "timeseries_"   + std::to_string(rank) + ".txt";
  std::string hfieldFilename = "hfield_"       + std::to_string(rank) +  ".h5";
  std::string hsnapsFilename = "snapshots_"    + std::to_string(rank) +  ".h5";
  std::string hspecFilename  = "hspec_"        + std::to_string(rank) + ".txt";
  std::string pinsetFilename = "pinned_sites_" + std::to_string(rank) + ".txt";
  const char* cXtend = hsnapsFilename.c_str(); 
  const char* cField = hfieldFilename.c_str();
  const char* cSpect = hspecFilename.c_str();
  char inputFieldFilename[25] = {};
  char srank[5]; sprintf(srank, "%d", rank);
  strcat(inputFieldFilename, "hfield_eq_");
  strcat(inputFieldFilename, srank);    
  strcat(inputFieldFilename, ".h5");
  
  /*----------------- Input variables declaration --------------*/

  int isSim;         // Check if we start from earlier simulation
  int wSnap;         // If set to 1, it writes snapshots in extendible HDF5 
  int sampleEvery;   // Sample every these accepted moves
  double maxiter;    // Max number of iterations
  double sig;        // Internal tension (k_BT/a_0^2)
  double tau;        // Frame tension    (k_BT/a_0^2)
  double epsilon;    // Maximum height perturbation
  double minChange;  // Min fraction of lattice size change
  double maxChange;  // Max fraction of lattice size change
  double pinRatio;   // Ratio of the total DoFs to be pinned
  double Eactive;    // Membrane activity (k_BT)
  
  /*-------------------------- Variables definition -------------------------*/

  int N = InputDoFs(argc,argv);   // Degrees of freedom (DoFs) per dimension
  int DoF = pow(N,2);             // Total number of DoFs
  int nGhost = 2;                  // Ghost points per boundary point
  int attemptLatticeChange = 5;  // Iterations to attempt a lattice change
  int iter = 0;                    
  double rig = 10.0;               // Bending rigidity (k_BT)
  double potStrength = 14000.;    // Strength of the pinning potential
  double h0 = 0.;                  // Equilibrium position of pinned sites
  double alpha = 1.0;              // Lattice spacing (distance between 2 DoFs)
 
  /* Read the input files. */ 
  
  ReadInput(inputFilename, isSim, wSnap, sampleEvery, maxiter, sig, tau,
	    epsilon, minChange, maxChange, pinRatio, Eactive);
  
  PrintOut(1,rank); // Input file is read

  /* Output the parameters of the run to txt files. */
  
  OutputParams(maxiter, N, DoF, nGhost, rig, sig, tau, epsilon, minChange,//this can be done from python gen.input
	       maxChange, alpha, pinRatio, sampleEvery, rank, Eactive);

  if (wSnap ==1 )
    hid_t file = H5Fcreate(cXtend, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      
  /*---------- Definition of variables needed for Monte Carlo --------------*/

  double prjArea = 0.0;         // Projected membrane area
  double totArea = 0.0;         // Total membrane area
  double tauEnergy = 0.0;       // Energy due to frame tension
  double crvEnergy = 0.0;       // Energy due to curvature
  double sigEnergy = 0.0;       // Energy due to internal tension
  double corEnergy = 0.0;       // Entropy correction energy
  double pinEnergy = 0.0;       // Pinning energy
  double totEnergy = 0.0;       // Total energy
  double localEnergyPre = 0.0;  // Local energy before the attempted move
  double localEnergyAft = 0.0;  // Local energy after the attempted move
  double localAreaPre = 0.0;    // Local area before the attempted move
  double localAreaAft = 0.0;    // Local area before the attempted move
  double dA_local = 0.0;        // Difference in local area
  double dE_local = 0.0;        // Difference in local energy
  double perturb;               // Value of the random perturbation
  int latticeAttempts = 0;      // Number of lattice change attempts
  int latticeChanges = 0;       // Number of accepted lattice moves
  int heightChanges = 0;        // Number of accepted height moves
  int totalMoves = 0;           // Total accepted moves
  int x,y;                      // Coordinates of selected site
  bool latticeAccept = 0;       // Indicates the acceptance of a lattice move
  bool heightAccept = 0;        // Indicates the acceptance of a height move
  bool where = 0;               // Indicator of boundary or bulk point
  bool pin = 0;                 // Indicator of pinned point
  
  /* Different number of neighbors are
     needed for the calculation
     of the change in a local quantity. 
     See GetNeighbors in McmemLib.cpp.    */

  int lenArea = 5;              
  int lenCorr = 4*nGhost+1;
  int lenEner = pow( (2*nGhost+1) ,2);
  Site site;
  Site neighborsArea[lenArea];
  Site neighborsCorr[lenCorr];
  Site neighborsEner[lenEner];
  std::unordered_set<Site> pinnedSites={};
 
  /*------------- Block-pinning ---------------*/

  int blockRadius = 0;
  int blockLength = pow( (2*blockRadius+1) ,2);
  Site neighbors[blockLength];
  
  /*-- Set up spectrum computations --*/

  int specSteps = 0;
  int qdiagMax = floor(sqrt(2)*N)+1;
  double L_mean = (double) N;
  double dk = 2.0*PI/L_mean;
  double* hx  = new double[N*(N+2)]();
  double* s1D = new double[qdiagMax]();
  fftw_complex* hq = (fftw_complex*) hx;
  fft_setup2d(N,hx,hq);
  
  /*-------- Random number generators and activity addition ------------*/

  std::uniform_int_distribution<int> RandInt(0,N-1);  
  std::uniform_real_distribution<double>  RandDouble(-epsilon,epsilon);
  AddShift(Eactive); // Shifts energy in metropolis for the active case

  /*------------------------ START OF ALGORITHM --------------------------*/
  
  /*------- (1) Initialize pinning and the height field hfield(i,j) ------*/
  RectMesh hfield(N,N,nGhost);

  if (isSim == 0 and pinRatio == 0)
    InitSurface(hfield,pinnedSites,-0.1,+0.1);
  
  if (isSim == 0 and pinRatio != 0){
    pinnedSites = InitPinning(hfield, pinRatio, neighbors, blockRadius); 
    WritePinnedSites(pinsetFilename, pinnedSites);
    InitSurface(hfield, pinnedSites, -0.1, +0.1);}
  
  if (isSim == 1 and pinRatio == 0)
    hfield.readH5(inputFieldFilename);
  
  if (isSim == 1 and pinRatio != 0){
    ReadPinnedSites(pinsetFilename, pinnedSites);
    hfield.readH5(inputFieldFilename);}
  
  PrintOut(2,rank); // Height field and pinning ok
  
  /* (2) Calculate energies and areas of the membrane and write the data. */

  CalculateTotal(hfield, rig, sig, tau, totEnergy, tauEnergy, crvEnergy,
		 sigEnergy, corEnergy, pinEnergy, totArea, prjArea, alpha,
   		 pinnedSites, potStrength, h0);
  
  Sample(outputFilename, iter, totalMoves, totEnergy, crvEnergy,
	    corEnergy, pinEnergy, totArea, prjArea, alpha, DoF);
  
  PrintOut(3,rank); // MC-Loop initialized
  if (isSim == 0) PrintOut(4,rank); // Spectrum averaged over non-equilibrium
  
  /*----------------------------------MC Loop---------------------------------*/
  
  for (iter = 1; iter < maxiter+1; iter ++)
    {
      
      /* (3) Randomly choose a lattice site (i,j), check whether it 
   	 belongs to the boundaries or to the bulk, if it is a pinned
   	 site or not, and find and store all its neighbors. */
      
      x = RandInt(mt); y = RandInt(mt); site.set(x,y);
      GetNeighbors(hfield, site, neighborsArea, neighborsCorr, neighborsEner);
      where = WhereIs(site,N,N,nGhost); pin = Ispinned(site, pinnedSites);
      
      /* (4) Calculate the local area and energy of that point. */
      
      localEnergyPre = LocalEnergy(hfield, neighborsArea, neighborsCorr,
				   neighborsEner, alpha, rig, sig, tau,
				   potStrength, h0, pin);
      
      localAreaPre = LocalArea(hfield, neighborsArea, alpha);
      
      /* (5) Randomly perturbate the height of the chosen point. */
      
      perturb = RandDouble(mt);
      hfield(x,y) += perturb;
      if (where == 1) GhostCopy(hfield);
      
      /* (6) Calculate the new local energy and local area. */
      
      localEnergyAft = LocalEnergy(hfield, neighborsArea, neighborsCorr,
				   neighborsEner, alpha, rig, sig, tau,
				   potStrength, h0, pin);
      
      localAreaAft = LocalArea(hfield, neighborsArea, alpha);
      
      /* (7) Calculate the difference in area and energy. */
      
      dA_local = localAreaAft-localAreaPre;
      dE_local = localEnergyAft-localEnergyPre;
      
      /* (8) Metropolis criterion for move acceptance
	 If the move is accepted, update total area, total energy, 
	 and sample. Otherwise, return to previous state. */
      
      heightAccept = UpdateState(hfield, site, where, totArea, totEnergy,
			   dA_local, dE_local, heightChanges, perturb);
      
      if (heightAccept)
	{
	  totalMoves++;
	  if (totalMoves % sampleEvery == 0)
	    {
	  
	      Sample(outputFilename,iter,totalMoves,totEnergy,
		     crvEnergy,corEnergy,pinEnergy,totArea,prjArea,alpha,DoF);
	      
	      specSteps ++;
	      CopyFieldToArray(hfield,hx);
	      fft();
	      onedspec2d(s1D,N,hx,alpha,dk,qdiagMax);
	      hfield.writeH5(cField);
	      if (wSnap == 1)
		WriteToExtendibleH5(cXtend, hfield);
	    } 
	}	  
      
      
      /* (9) After "attempt_lattice_change" iterations, randomly change 
   	 alpha, compute the new projected area and update the 
   	 total energy. */
      
      if (iter % attemptLatticeChange == 0)
	{
	  latticeAccept = UpdateLattice(hfield,minChange,maxChange,rig,sig,tau,
					prjArea, totArea, totEnergy,
					tauEnergy, crvEnergy, sigEnergy,
					corEnergy, pinEnergy, alpha,
					latticeAttempts, latticeChanges,
					pinnedSites, potStrength, h0);
	  if (latticeAccept)
	    {
	      totalMoves++;
	      if ( totalMoves % sampleEvery == 0 )
		{

		  Sample(outputFilename,iter,totalMoves,totEnergy,
			 crvEnergy,corEnergy,pinEnergy,totArea,prjArea,alpha,DoF);

		  specSteps ++;
		  CopyFieldToArray(hfield,hx);
		  fft();
		  onedspec2d(s1D,N,hx,alpha,dk,qdiagMax);
		  hfield.writeH5(cField);
		  if (wSnap == 1) WriteToExtendibleH5(cXtend, hfield);
		}
	    }	 
	}
  
    }// end of MC-loop

  /*------------------------ END OF ALGORITHM -----------------------------*/
  
  /* Attach metadata to extendible HDF5 set */

  if (wSnap == 1)
    {
      hfield_metadata wdata[9]; 
      wdata[0].value = specSteps;  wdata[0].field = "Samples";  
      wdata[1].value = sig;        wdata[1].field = "Sigma";
      wdata[2].value = tau;        wdata[2].field = "Tau";
      wdata[3].value = rig;        wdata[3].field = "Rigidity";
      wdata[4].value = pinRatio;   wdata[4].field = "Fraction of pinning";
      wdata[5].value = Eactive;    wdata[5].field = "Activity";
      wdata[6].value = N;          wdata[6].field = "Rows";
      wdata[7].value = N;          wdata[7].field = "Cols";
      wdata[8].value = nGhost;     wdata[8].field = "Ghost points";      
      WriteMetadataToH5File(cXtend,wdata,9);
    }

  PrintOut(5,rank); // MC-Loop finished successfully
  
  /* Average power spectrum and write it to a file */
  
  WriteSpectrum(hspecFilename, s1D, specSteps, qdiagMax, dk);
  
  /* Write acceptance ratios and number of spectrum calculations. */
  
  WStats(maxiter, heightChanges, latticeAttempts,
	 latticeChanges, specSteps, rank);
  
  MPI_Finalize();
  PrintOut(6,rank); //Program terminated successfully
  return 0;
}
