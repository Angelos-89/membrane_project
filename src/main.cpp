#include <mpi.h>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <random>
#include "McmemLib.hpp"

/* Hash function for the Site class */
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
        and frame tensions.                                                   */ 
  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  int acc_samples;
  double maxit,e,s,t,minchange,maxchange,pin_ratio;
  double Eactive = 0;
  /* Eactive is the active energy, zero by default.
     If you set the Eactive to non-zero values the 
     membrane is no longer in equilibrium. The value 
     can be both positive or negative. 
     See Kumar & Dasgupta PRE 102, 2020 */
  
  std::string input_filename  = "input_"      + std::to_string(rank) + ".txt";
  std::string output_filename = "timeseries_" + std::to_string(rank) + ".txt";
  std::string hfield_filename = "hfield_"     + std::to_string(rank) + ".h5";
  const char* cc = hfield_filename.c_str();

  ReadInput(input_filename,maxit,s,t,e,
	    minchange,maxchange,pin_ratio,acc_samples,Eactive);
  
  const int maxiter = maxit;           //max no of iterations
  const int N = 80;                    //DoF per dimension
  const int DoF = N*N;                 //number of degrees of freedom
  const int nghost = 2;                //ghost points per boundary point
  const double rig = 10.0;             //bending rigidity
  const double sig = s;                //internal tension
  const double tau = t;                //frame tension
  const double epsilon = e;            //maximum possible height perturbation
  const double min_change = minchange; //min percentage of lattice size change
  const double max_change = maxchange; //max percentage of lattice size change
  const double pn_prcn = pin_ratio;    //ratio of the total DoF to be pinned
  const double pot_strength = 14000;   //strength of the pinning potential
  const double h0 = 0;                 //equilibrium position of pinned sites
  double alpha = 1.0;                  //lattice spacing(distance between 2 DoF)
  int sample_every = acc_samples;      //sample when acc_samples are accepted
  int attempt_lattice_change = 5;      //iterations to attempt a lattice change
  int iter = 0;
  
  OutputParams(maxiter,N,DoF,nghost,rig,sig,tau,epsilon,
	       min_change,max_change,alpha,pn_prcn,sample_every,rank,Eactive);
  
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
  
  bool where = 0;              //indicator of boundary or bulk point
  bool pin = 0;                //indicator of pinned point
  bool accept = 0;             //indicates the acceptance of a height move
  bool lattice_accept = 0;     //indicates the acceptance of a lattice move

  int lattice_changes = 0;     //number of accepted lattice moves
  int height_changes = 0;      //number of accepted height moves
  int lattice_moves = 0;       //number of lattice change attempts
  int total_moves = 0;         //total accepted moves
  
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

  AddShift(Eactive); //shift energy in metropolis to implement "activity".
  RectMesh OutSpec(N+2,N);
  RectMesh SpecRes(N+2,N);
  int spec_steps = 0;
  int spec_every = 1e3;
  
  /* 2) Initialize pinning and the height field hfield(i,j)                   */

  RectMesh hfield(N,N,nghost);
  pinned_sites = InitPinning(N,pn_prcn);      //store pinned sites to a set
  InitSurface(hfield,pinned_sites,-0.1,+0.1); //initialize a random surface
    
  /* 3) Calculate the projected membrane area "prj_area", the total area 
     "tot_area" and the energies "tau_energy","sig_energy","crv_energy",
     "cor_energy" and "tot_energy" and write the data.                        */
    
  CalculateTotal(hfield,DoF,rig,sig,tau,tot_energy,tau_energy,crv_energy,
  		 sig_energy,cor_energy,pin_energy,tot_area,prj_area,alpha,
		 pinned_sites,pot_strength,h0);
  
  Sample(iter,total_moves,output_filename,tot_energy,crv_energy,
	 cor_energy,pin_energy,tot_area,prj_area,alpha,DoF);
  
  /*----------------------------------MC Loop---------------------------------*/
  
  for (iter=1; iter<maxiter+1; iter++)
    {
      
      /* 4) Randomly choose a lattice site (i,j) and check whether it 
  	 belongs to the boundaries or to the bulk, and if it is a pinned
	 site. Also find and store all its neighbors.                         */

      x = RandInt(MT);
      y = RandInt(MT);
      site.set(x,y);
      GetNeighbors(hfield,site,neighbors_area,neighbors_corr,neighbors_ener);
      where = WhereIs(site,N,N,nghost);
      pin = Ispinned(site,pinned_sites);
      
      /* 5) Calculate the local area and energy of that point.                */

      local_energy_pre = LocalEnergy(hfield,neighbors_area,
  				     neighbors_corr,
  				     neighbors_ener,
  				     alpha,rig,sig,tau,
				     pot_strength,h0,pin);

      local_area_pre = LocalArea(hfield,neighbors_area,alpha);
      
      /* 6) Randomly perturbate the height of the chosen point.               */

      perturb = RandDouble(MT); 
      hfield(x,y) += perturb;
      if (where==1)
	GhostCopy(hfield);
      
      /* 7) Calculate the new local energy and local area.                    */

      local_energy_aft = LocalEnergy(hfield,neighbors_area,
  				     neighbors_corr,
  				     neighbors_ener,
  				     alpha,rig,sig,tau,
				     pot_strength,h0,pin);

      local_area_aft = LocalArea(hfield,neighbors_area,alpha);
      
      /* 8) Calculate the energy difference and check whether the 
  	 move is accepted or not.                                             */

      dAlocal = local_area_aft - local_area_pre;
      dElocal = local_energy_aft - local_energy_pre;
      accept  = Metropolis(dElocal);
      
      /* 9) If the move is accepted, update total area and energy.
  	 Otherwise return to previous state.                                  */
      
      AcceptOrDecline(hfield,site,accept,where,tot_area,
  		      tot_energy,dAlocal,dElocal,height_changes,perturb);
      
      if (accept) //if the move is accepted update total_moves
	{
	  total_moves ++;
	  if (total_moves % sample_every == 0) //write every sample_every moves 
	    Sample(iter,total_moves,output_filename,
		   tot_energy,crv_energy,cor_energy,pin_energy,
		   tot_area,prj_area,alpha,DoF);
	}
      
      /* 10) After "attempt_lattice_change" iterations, randomly change 
	 alpha, compute the new projected area and update the 
	 total energy.                                                        */
      
      if (iter % attempt_lattice_change == 0)
	{
	  lattice_accept = ChangeLattice(hfield,min_change,max_change,
					 DoF,rig,sig,tau,prj_area,tot_area,
					 tot_energy,tau_energy,crv_energy,
					 sig_energy,cor_energy,pin_energy,
					 alpha,lattice_moves,lattice_changes,
					 pinned_sites,pot_strength,h0);
	  
	  if (lattice_accept) //if the move is accepted update total_moves
	    {  
	      total_moves ++;
	      if (total_moves % sample_every == 0)  
		Sample(iter,total_moves,output_filename,
		       tot_energy,crv_energy,cor_energy,pin_energy,
		       tot_area,prj_area,alpha,DoF);
	    }
	}

      /* Compute Spectrum */

      if (total_moves > (int) 1e5)
	{
	  spec_steps ++;
	  Spectrum(hfield,OutSpec);
	  SpecRes += OutSpec;
	}
      SpecRes = SpecRes/spec_steps;

      
      if (total_moves % (int) 1e3 == 0) //write surface every 1e3 accepted moves
	hfield.writeH5(cc);
    }
  
  /* 11) Print acceptance ratios and finish                                   */

  PrintAcceptance(maxiter,height_changes,lattice_moves,lattice_changes,rank);
  std::cout << lattice_changes << std::endl;
  MPI_Finalize();
  return 0;
}

