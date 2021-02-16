#include <mpi.h>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <random>
#include "McmemLib.hpp"

namespace std
{
  template <>
  struct hash<Site>
  {
    size_t operator()( const Site& p ) const
    {
      return((53 + std::hash<int>()(p.getx()))*53
	         + std::hash<int>()(p.gety()));
    }
  }; 
}

std::random_device RD;
std::mt19937 MT(RD()); 

int main(int argc, char* argv[])
{

  /* 0) Initialize MPI and read files with the values of internal 
        and frame tensions.                                                 */ 

  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  double s; // to store the internal tension               
  double t; // to store the frame tension
  std::string input_filename  = "input"    + std::to_string(rank) + ".txt";
  std::string output_filename = "sampling" + std::to_string(rank) + ".txt";
  std::string hfield_filename = "hfield"   + std::to_string(rank) + ".h5";
  const char* cc = hfield_filename.c_str();
  ReadTensions(input_filename,s,t);

  /* 1) Set values for number of degrees of freedom "DoF", bending rigidity 
     "rig" and lattice spacing "alpha". Define and initialize variables 
     needed for the MC code.                                                */

  const int maxiter = 1e8;           //max no of iterations
  const int N = 80;                  //DoF per dimension
  const int DoF = N*N;               //Number of degrees of freedom
  const int nghost = 2;              //ghost points per boundary point
  const double rig = 10.0;           //bending rigidity
  const double sig = s;              //internal tension
  const double tau = t;              //frame tension
  const double epsilon = 0.34;       //maximum possible height perturbation
  const double min_change = 0.98;    //min percentage of lattice size change
  const double max_change = 1.02;    //max percentage of lattice size change
  double alpha = 1.0;                //lattice spacing (distance between 2 DoF)
  const double pn_prcn = 0.1;        //percentage of the total DoF to be pinned
  const double pot_strength = 14000; //strength of the pinning potential
  const double h0 = 0;               //equilibrium position of pinned sites
  
  double prj_area = 0.0;
  double tot_area = 0.0;
  double tau_energy = 0.0;
  double crv_energy = 0.0;
  double sig_energy = 0.0;
  double cor_energy = 0.0;
  double pin_energy = 0.0;
  double tot_energy = 0.0;
  
  double local_energy_pre = 0.0;
  double local_energy_aft = 0.0;
  double local_area_pre = 0.0;
  double local_area_aft = 0.0;
  double dAlocal = 0.0;        //difference in local area
  double dElocal = 0.0;        //difference in local energy
  double perturb;              //value of the random perturbation
  
  bool where;                  //indicator of boundary or bulk point
  bool pin;                    //indicator of pinned point
  bool accept;                 //indicates the acceptance of a Monte-Carlo move
  
  int lattice_changes = 0;
  int accepted_moves = 0;
  int lattice_moves = 0;
  int move_counter = 0;         

  Site site;
  int x,y;
  int len_area = 5; 
  int len_corr = 4*nghost+1;
  int len_ener = pow((2*nghost+1),2); 
  Site neighbors_area[len_area] = {};
  Site neighbors_corr[len_corr] = {};
  Site neighbors_ener[len_ener] = {};
  std::unordered_set<Site> pinned_sites;
  
  std::uniform_int_distribution<int>      RandInt(0,N-1);  
  std::uniform_real_distribution<double>  RandDouble(-epsilon,+epsilon);

  /* 2) Initialize pinning and the height field hfield(i,j)                 */

  RectMesh hfield(N,N,nghost);
  pinned_sites = InitPinning(N,pn_prcn);      //store pinned sites to a set
  InitSurface(hfield,-0.1,+0.1,pinned_sites); //initialize a random surface
  
  /* 3) Calculate the projected membrane area "prj_area", the total area 
     "tot_area" and the energies "tau_energy","sig_energy","crv_energy",
     "cor_energy" and "tot_energy".                                       */
  
  CalculateTotal(hfield,DoF,rig,sig,tau,tot_energy,tau_energy,crv_energy,
  		 sig_energy,cor_energy,pin_energy,tot_area,prj_area,alpha,
		 pinned_sites,pot_strength,h0);

  /*---------------------------------MC-Loop------------------------------*/
  
  for (int iter=0; iter<maxiter; iter++)
    {
      
      /* 4) Randomly choose a lattice site (i,j) and check whether it 
  	 belongs to the boundaries or to the bulk, and if it is a pinned
	 site. Also find and store all its neighbors.                     */

      x = RandInt(MT);
      y = RandInt(MT);
      site.set(x,y);
      GetNeighbors(hfield,site,neighbors_area,neighbors_corr,neighbors_ener);
      where = WhereIs(site,N,N,nghost);
      pin = Ispinned(site,pinned_sites);
      
      /* 5) Calculate the local area and energy of that point.            */

      local_energy_pre = LocalEnergy(hfield,neighbors_area,
  				     neighbors_corr,
  				     neighbors_ener,
  				     alpha,rig,sig,tau,pot_strength,h0,pin);
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
  				     alpha,rig,sig,tau,pot_strength,h0,pin);
      local_area_aft = LocalArea(hfield,neighbors_area,alpha);
      
      /* 8) Calculate the energy difference and check whether the 
  	 move is accepted or not.                                         */

      dAlocal = local_area_aft   - local_area_pre;
      dElocal = local_energy_aft - local_energy_pre;
      accept  = Metropolis(dElocal);
      
      /* 9) If the move is accepted, update total area and energy.
  	 Otherwise return to previous state.                              */

      AcceptOrDecline(hfield,site,accept,where,tot_area,
  		      tot_energy,dAlocal,dElocal,accepted_moves,perturb);
      
      /* 10) After 5 MC steps, randomly change alpha, compute the 
  	 new projected area and update the total energy.                  */

      ChangeLattice(hfield,min_change,max_change,DoF,rig,sig,
		    tau,prj_area,tot_area,tot_energy,tau_energy,
		    crv_energy,sig_energy,cor_energy,pin_energy,
		    alpha,move_counter,lattice_moves,lattice_changes,
		    pinned_sites,pot_strength,h0);
      
      /* 11) Sample                                                       */

      Sample(iter,output_filename,tot_energy,tau_energy,crv_energy,
	     sig_energy,cor_energy,pin_energy,tot_area,prj_area,alpha,DoF);

      if (iter % (int) 1e5 == 0)
  	hfield.writeH5(cc);
    }

  /* 12) Print acceptance ratios                                          */

  PrintAcceptance(maxiter,accepted_moves,lattice_moves,lattice_changes);

  MPI_Finalize();
  return 0;
}

