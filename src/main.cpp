#include <string.h>
#include <mpi.h>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <random>
#include "McmemLib.hpp"
#include "fft.h"

/*------------------------ Hash function for the Site class ------------------*/
namespace std
{
  template <>
  struct hash<Site>
  {
    size_t operator()( const Site& p ) const
    {
      return((53 + std::hash<int>()(p.getx()))*53 + std::hash<int>()(p.gety()));
    }
  }; 
}

/*------ Global -------*/
double PI = 4.*atan(1.0);
std::random_device RD;
std::mt19937 MT(RD());

/*---------------------------- Main Function ---------------------------------*/
int main(int argc, char* argv[])
{

  /*------ Initialize MPI. -------- */
  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  /*----------------- Define strings needed for filenames. -------------------*/
  std::string input_filename  = "input_"        + std::to_string(rank) + ".txt";
  std::string output_filename = "timeseries_"   + std::to_string(rank) + ".txt";
  std::string hfield_filename = "hfield_"       + std::to_string(rank) +  ".h5";
  std::string hsnaps_filename = "snapshots_"    + std::to_string(rank) +  ".h5";
  std::string hspec_filename  = "hspec_"        + std::to_string(rank) + ".txt";
  std::string pinset_filename = "pinned_sites_" + std::to_string(rank) + ".txt";
  const char* cextend = hsnaps_filename.c_str(); 
  const char* cfield = hfield_filename.c_str();
  const char* cspec =  hspec_filename.c_str();
  // The lines below are needed to read equilibrated height fields if is_sim==1.
  char input_field_filename[25] = {};
  char srank[5];
  sprintf(srank, "%d", rank);
  strcat(input_field_filename, "hfield_eq_");
  strcat(input_field_filename, srank);    
  strcat(input_field_filename, ".h5");
  
  /*----------------- Input variables declaration --------------*/
  // These are read from an input file.
  int is_sim;         // Check if we start from earlier simulation
  int sample_every;   // Sample every these accepted moves
  double maxiter;     // Max number of iterations
  double sig;         // Internal tension (k_BT/a_0^2)
  double tau;         // Frame tension    (k_BT/a_0^2)
  double epsilon;     // Maximum height perturbation
  double min_change;  // Min fraction of lattice size change
  double max_change;  // Max fraction of lattice size change
  double pin_ratio;   // Ratio of the total DoFs to be pinned
  double Eactive;     // Membrane activity (k_BT)
  
  //why not read these from a file too and the Output to txt is done in Python?
  /*-------------------------- Variable definition --------------------------*/
  int N = Input_DoFs(argc,argv);   // Degrees of freedom (DoFs) per dimension
  int DoFs = pow(N,2);             // Total number of DoFs
  int Nghost = 2;                  // Ghost points per boundary point
  int attempt_lattice_change = 5;  // Iterations to attempt a lattice change
  int iter = 0;                    
  double rig = 10.0;               // Bending rigidity (k_BT)
  double pot_strength = 14000.;    // Strength of the pinning potential
  double h0 = 0.;                  // Equilibrium position of pinned sites
  double alpha = 1.0;              // Lattice spacing (distance between 2 DoFs)
 
  /* Read the input files. */ 
  
  ReadInput(input_filename, is_sim, sample_every, maxiter, sig, tau, epsilon,
	    min_change, max_change, pin_ratio, Eactive);
  
  PrintOut(1,rank); // Input file is read

  /* Output the parameters of the run to txt files. */
  
  OutputParams(maxiter, N, DoFs, Nghost, rig, sig, tau, epsilon, min_change,
	       max_change, alpha, pin_ratio, sample_every, rank, Eactive);

  if (is_sim == 1)
    {hid_t file = H5Fcreate(cextend, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);}
      
  
  /*---------- Definition of variables needed for Monte Carlo --------------*/
  double prj_area = 0.0;         // Projected membrane area
  double tot_area = 0.0;         // Total membrane area
  double tau_energy = 0.0;       // Energy due to frame tension
  double crv_energy = 0.0;       // Energy due to curvature
  double sig_energy = 0.0;       // Energy due to internal tension
  double cor_energy = 0.0;       // Entropy correction energy
  double pin_energy = 0.0;       // Pinning energy
  double tot_energy = 0.0;       // Total energy
  double local_energy_pre = 0.0; // Local energy before the attempted move
  double local_energy_aft = 0.0; // Local energy after the attempted move
  double local_area_pre = 0.0;   // Local area before the attempted move
  double local_area_aft = 0.0;   // Local area before the attempted move
  double dA_local = 0.0;         // Difference in local area
  double dE_local = 0.0;         // Difference in local energy
  double perturb;                // Value of the random perturbation
  int lattice_attempts = 0;      // Number of lattice change attempts
  int lattice_changes = 0;       // Number of accepted lattice moves
  int height_changes = 0;        // Number of accepted height moves
  int total_moves = 0;           // Total accepted moves
  int x,y;                       // Coordinates of selected site
  bool lattice_accept = 0;       // Indicates the acceptance of a lattice move
  bool accept = 0;               // Indicates the acceptance of a height move
  bool where = 0;                // Indicator of boundary or bulk point
  bool pin = 0;                  // Indicator of pinned point
  
  /* Different number of neighbors are
     needed for the calculation
     of the change in a local quantity. 
     See GetNeighbors in McmemLib.cpp.    */
  int len_area = 5;              
  int len_corr = 4*Nghost+1;
  int len_ener = pow( (2*Nghost+1) ,2);
  Site site;
  Site neighbors_area[len_area];
  Site neighbors_corr[len_corr];
  Site neighbors_ener[len_ener];
  std::unordered_set<Site> pinned_sites={};
 
  /*------------- Block-pinning ---------------*/
  int block_radius = 0;
  int block_length = pow( (2*block_radius+1) ,2);
  Site neighbors[block_length];
  
  /*-- Set up spectrum computations --*/
  int spec_steps = 0;
  int qdiag_max = floor(sqrt(2)*N)+1;
  double L_mean = (double) N;
  double dk = 2.0*PI/L_mean;
  double* hx  = new double[N*(N+2)]();
  double* S1d = new double[qdiag_max]();
  fftw_complex* hq = (fftw_complex*) hx;
  fft_setup2d(N,hx,hq);
  
  /*-------- Random number generators and activity addition ------------*/
  std::uniform_int_distribution<int> RandInt(0 , N-1);  
  std::uniform_real_distribution<double>  RandDouble(-epsilon,epsilon);
  AddShift(Eactive); // Shifts energy in metropolis

  /*------------------------ START OF ALGORITHM --------------------------*/
  
  /*------- (1) Initialize pinning and the height field hfield(i,j) ------*/
  RectMesh hfield(N,N,Nghost);

  if (is_sim == 0 and pin_ratio == 0)
    InitSurface(hfield,pinned_sites,-0.1,+0.1);
  
  if (is_sim == 0 and pin_ratio != 0){
    pinned_sites = InitPinning(hfield, pin_ratio, neighbors, block_radius); 
    WritePinnedSites(pinset_filename,pinned_sites);
    InitSurface(hfield,pinned_sites,-0.1,+0.1);}
  
  if (is_sim == 1 and pin_ratio == 0)
    hfield.readH5(input_field_filename);
  
  if (is_sim == 1 and pin_ratio != 0){
    ReadPinnedSites(pinset_filename, pinned_sites);
    hfield.readH5(input_field_filename);}
  
  PrintOut(2,rank); // Height field and pinning ok
  
  /*----------------------------------------------------------------------*/
  
  /* (2) Calculate energies and areas of the membrane and write the data. */
  CalculateTotal(hfield, rig, sig, tau, tot_energy, tau_energy, crv_energy,
		 sig_energy, cor_energy, pin_energy, tot_area, prj_area, alpha,
   		 pinned_sites, pot_strength, h0);
  
  WriteData(output_filename, iter, total_moves, tot_energy, crv_energy,
	    cor_energy, pin_energy, tot_area, prj_area, alpha, DoFs);
  /*----------------------------------------------------------------------*/
  
  PrintOut(3,rank); // MC-Loop initialized
  if (is_sim==0) PrintOut(4,rank); // Spectrum averaged over non-equilibrium
  
  /*----------------------------------MC Loop---------------------------------*/
  
  for (iter=1; iter<maxiter+1; iter++)
    {
      
      /* (3) Randomly choose a lattice site (i,j), check whether it 
   	 belongs to the boundaries or to the bulk, if it is a pinned
   	 site or not, and find and store all its neighbors. */
      
      x = RandInt(MT); y = RandInt(MT); site.set(x,y);
      GetNeighbors(hfield, site, neighbors_area, neighbors_corr, neighbors_ener);
      where = WhereIs(site,N,N,Nghost); pin = Ispinned(site,pinned_sites);
      
      /* (4) Calculate the local area and energy of that point. */
      
      local_energy_pre = LocalEnergy(hfield, neighbors_area, neighbors_corr,
				     neighbors_ener, alpha, rig, sig, tau,
				     pot_strength, h0, pin);
      
      local_area_pre = LocalArea(hfield, neighbors_area, alpha);
      
      /* (5) Randomly perturbate the height of the chosen point. */
      
      perturb = RandDouble(MT);
      hfield(x,y) += perturb;
      if (where == 1)
	{GhostCopy(hfield);}
      
      /* (6) Calculate the new local energy and local area. */
      
      local_energy_aft = LocalEnergy(hfield, neighbors_area, neighbors_corr,
				     neighbors_ener, alpha, rig, sig, tau,
				     pot_strength, h0, pin);
      
      local_area_aft = LocalArea(hfield, neighbors_area, alpha);
      
      /* (7) Calculate the energy difference and check whether the 
	 move is accepted or not. */
      
      dA_local = local_area_aft-local_area_pre;
      dE_local = local_energy_aft-local_energy_pre;
      accept = Metropolis(dE_local);
      
      /* (8) If the move is accepted, update total area and energy and sample.
	 Otherwise, return to previous state. */
      
      UpdateState(hfield, site, accept, where, tot_area, tot_energy,
		  dA_local, dE_local, height_changes, perturb);
      
      if (accept){
	total_moves++;
	Sample(output_filename,sample_every,iter,total_moves,tot_energy,
	       crv_energy,cor_energy,pin_energy,tot_area,prj_area,alpha,DoFs);}
      
      /* (9) After "attempt_lattice_change" iterations, randomly change 
   	 alpha, compute the new projected area and update the 
   	 total energy. */
      
      if (iter % attempt_lattice_change == 0){
	lattice_accept = ChangeLattice(hfield, min_change, max_change, rig, sig,
				       tau, prj_area, tot_area, tot_energy,
				       tau_energy, crv_energy, sig_energy,
				       cor_energy, pin_energy, alpha,
				       lattice_attempts, lattice_changes,
				       pinned_sites, pot_strength, h0);
	if (lattice_accept){
	  total_moves++;
	  Sample(output_filename,sample_every,iter,total_moves,tot_energy,
		 crv_energy,cor_energy,pin_energy,tot_area,prj_area,alpha,DoFs);} 
      }
      
      /* (10) Compute radial 1D spectrum and write the height field.*/
      
      if (total_moves % sample_every == 0)
	{
	  spec_steps ++;
	  CopyFieldToArray(hfield,hx); fft();
	  onedspec2d(S1d,N,hx,alpha,dk,qdiag_max);
	  hfield.writeH5(cfield);
	  if (is_sim == 1) // Only for equilibrated simulations
	    Write_to_extendible_H5(cextend, hfield);
	}
      
    }// end of MC-loop

  /*------------------------ END OF ALGORITHM -----------------------------*/
  
  /* Attach metadata to extendible HDF5 set */
  if (is_sim == 1)
    {
      hfield_metadata wdata[9];  

      wdata[0].value = spec_steps; wdata[0].field = "Samples";  
      wdata[1].value = sig;        wdata[1].field = "Sigma";
      wdata[2].value = tau;        wdata[2].field = "Tau";
      wdata[3].value = rig;        wdata[3].field = "Rigidity";
      wdata[4].value = pin_ratio;  wdata[4].field = "Fraction of pinning";
      wdata[5].value = Eactive;    wdata[5].field = "Activity";
      wdata[6].value = N;          wdata[6].field = "Rows";
      wdata[7].value = N;          wdata[7].field = "Cols";
      wdata[8].value = Nghost;     wdata[8].field = "Ghost points";

      Write_metadata_to_H5_file(cextend, wdata, 9);
    }
    
  PrintOut(5,rank); // MC-Loop finished successfully
  
  /* Average power spectrum and write it to a file */
  
  WriteSpectrum(hspec_filename, S1d, spec_steps, qdiag_max, dk);
  
  /* Write acceptance ratios and number of spectrum calculations. */
  
  WStats(maxiter, height_changes, lattice_attempts,
	 lattice_changes, spec_steps, rank);
  
  MPI_Finalize();
  
  PrintOut(6,rank); //Program terminated successfully
  return 0;
}



