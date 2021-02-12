#include <mpi.h>
#include <iomanip>
#include <fstream>
#include <string>
#include <random>
#include "McmemLib.hpp"

std::random_device RD;
std::mt19937 MT(RD()); 

int main(int argc, char* argv[])
{

  /* 0) Initialize MPI and read files with the values of internal 
        and frame tensions.                                                 */ 

  int rank;
  std::cout << "Initiating MPI" << std::endl;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  std::cout << "Initiating MPI: DONE" << std::endl;
/* parameters, some of these are overwritten after reading input files */
  int miter = 1e3;
  const int DoF = 6400;              //number of degrees of freedom
  const int N = sqrt(DoF);           //DoF per dimension
  const int nghost = 2;              //ghost points per boundary point
  const double rig = 10.0;           //bending rigidity
  double s = 0.52;
  double t = 0.5;
  const double epsilon = 0.34;       //maximum possible height perturbation
  const double min_change = 0.98;    //min percentage of lattice size change
  const double max_change = 1.02;    //max percentage of lattice size change
  double alpha = 1.0;                //lattice spacing (distance between 2 DoF)
  double Eactive = 0.;                // Active energy, zero by default
/* If you set the Eactive to non-zero values the membrane is no longer is
 * equilibrium. The value can be both positive or negative. See Kumar &
 * Dasgupta PRE 102, 2020 */  
  double Ea = 0.; 
  std::string input_filename = "input" + std::to_string(rank) + ".txt";
  std::string output_filename = "sampling" + std::to_string(rank) + ".txt";
  std::string hfield_filename = "hfield" + std::to_string(rank) + ".h5";
  const char* cc = hfield_filename.c_str();
//  ReadTensions(input_filename,s,t);
  std::vector<double> input;
  ReadInputs(input_filename, input);
  s = input[0];
  t = input[1];
  std::cout << "Simulation " << rank+1 << ":\n";
  if (input.size() > 1){ 
     Ea = input[2]; 
     std::cout << "Eactive = " << std::setprecision(3) << Ea <<"\n"; 
  }
  if (input.size() > 2){ 
     miter = int(input[3]); 
     std::cout << "read miter = " <<  miter <<"\n"; 
  }
  const double sig = s;              //internal tension
  const double tau = t;              //frame tension
  Eactive = Ea;  
  const int maxiter = miter;           //max no of iterations
  std::cout << "Input file read.\ntau = " << std::setprecision(3) << tau
	    << "\nsigma = " << std::setprecision(6) << sig << std::endl;  
  std::cout << "Eactive = " << std::setprecision(3) << Eactive <<"\n"; 
  std::cout << "max iterations = " <<  maxiter <<"\n"; 
  /* 1) Set values for number of degrees of freedom "DoF", bending rigidity 
     "rig" and lattice spacing "alpha". Define and initialize variables 
     needed for the MC code.                                                */

  
  double prj_area = 0.0;
  double tot_area = 0.0;
  double tau_energy = 0.0;
  double crv_energy = 0.0;
  double sig_energy = 0.0;
  double cor_energy = 0.0;
  double tot_energy = 0.0;
  
  double local_energy_pre = 0.0;
  double local_energy_aft = 0.0;
  double local_area_pre = 0.0;
  double local_area_aft = 0.0;
  double dAlocal = 0.0;        //difference in local area
  double dElocal = 0.0;        //difference in local energy
  double perturb;              //value of the random perturbation
  
  bool where;                  //indicator of boundary or bulk point
  bool accept;                 //indicates the acceptance of a Monte-Carlo move
  
  int lattice_changes = 0;
  int accepted_moves = 0;
  int lattice_moves = 0;
  int move_counter = 0;         

  int x,y;
  int len_area = 5; 
  int len_corr = 4*nghost+1;
  int len_ener = pow((2*nghost+1),2); 
  Site neighbors_area[len_area] = {};
  Site neighbors_corr[len_corr] = {};
  Site neighbors_ener[len_ener] = {};
  Site site;
  
  std::uniform_int_distribution<int>      RandInt(0,N-1);  
  std::uniform_real_distribution<double>  RandDouble(-epsilon,epsilon);

  std::cout << "All variables defined.\n";
  std::cout << "----------------------" << std::endl;
  
  /* 2) Initialize the height field hfield(i,j)                           */ 

  RectMesh hfield(N,N,nghost);
  InitSurface(hfield,-0.1,+0.1);
  std::cout << "Simulation " << rank+1 <<":\n";
  std::cout << "Height field initialized.\n";
  
  /* 3) Calculate the projected membrane area "prj_area", the total area 
     "tot_area" and the energies "tau_energy","sig_energy","crv_energy",
     "cor_energy" and "tot_energy".                                       */

  std::cout << "Calculating total quantities of initial configuration.\n"
    "------------------------------------------------------" << std::endl;
  
  CalculateTotal(hfield,DoF,rig,sig,tau,tot_energy,tau_energy,crv_energy,
  		 sig_energy,cor_energy,tot_area,prj_area,alpha);

  std::cout << "Simulation " << rank+1 << ": Calculations complete."
    "Monte Carlo loop initiated...\n";
  
  /*---------------------------------MC-Loop------------------------------*/
  
  for (int iter=0; iter<maxiter; iter++)
    {
      
      /* 4) Randomly choose a lattice site (i,j) and check whether it 
  	 belongs to the boundaries or to the bulk. Also find and store all
  	 its neighbors.                                                   */

      x = RandInt(MT);
      y = RandInt(MT);
      site.set(x,y);
      GetNeighbors(hfield,site,neighbors_area,neighbors_corr,neighbors_ener);
      where = WhereIs(site,N,N,nghost);

      /* 5) Calculate the local area and energy of that point.            */
      local_energy_pre = LocalEnergy(hfield,neighbors_area,
  				     neighbors_corr,
  				     neighbors_ener,
  				     alpha,rig,sig,tau);
      local_area_pre = LocalArea(hfield,neighbors_area,alpha);
      
      /* 6) Randomly perturbate the height of the chosen point.           */
      perturb = RandDouble(MT); 
      hfield(x,y) += perturb;
      move_counter ++;
      if(where==1)
  	GhostCopy(hfield);
      
      /* 7) Calculate the new local energy and local area.                */

      local_energy_aft = LocalEnergy(hfield,neighbors_area,
  				     neighbors_corr,
  				     neighbors_ener,
  				     alpha,rig,sig,tau);
      local_area_aft = LocalArea(hfield,neighbors_area,alpha);
      
      /* 8) Calculate the energy difference and check whether the 
  	 move is accepted or not.                                         */

      dAlocal = local_area_aft - local_area_pre;
      dElocal = local_energy_aft - local_energy_pre;
      dElocal = dElocal + Eactive;
      accept  = Metropolis(dElocal);
      /* If we are solving for an equilibrium membrane the Eactive is zero */
      
      /* 9) If the move is accepted, update total area and energy.
  	 Otherwise return to previous state.                              */

      AcceptOrDecline(hfield,site,accept,where,tot_area,
  		      tot_energy,dAlocal,dElocal,accepted_moves,perturb);
      
      /* 10) After 5 MC steps, randomly change alpha, compute the 
  	 new projected area and update the total energy.                  */

      ChangeLattice(hfield,min_change,max_change,DoF,rig,sig,tau,prj_area,
		    tot_area,tot_energy,tau_energy,crv_energy,sig_energy,
		    cor_energy,alpha,move_counter,
		    lattice_moves,lattice_changes);
      
      /* 11) Sample                                                       */

      Sample(iter,output_filename,tot_energy,tau_energy,crv_energy,sig_energy,
  	     cor_energy,tot_area,prj_area,alpha,DoF);
      if (iter % (int) 1e5 == 0)
  	hfield.writeH5(cc);
    }

  /* 12) Print acceptance ratios                                          */

  std::cout <<"--------------------------------------------" << std::endl;
  std::cout << "Simulation " << rank+1 << " finished sucessfully." << std::endl;
  PrintAcceptance(maxiter,accepted_moves,lattice_moves,lattice_changes);
  MPI_Finalize();
  return 0;
}

