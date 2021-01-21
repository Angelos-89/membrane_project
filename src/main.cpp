#include <iomanip>
#include <fstream>
#include <chrono>
#include <random>
#include "McmemLib.hpp"

#include "PARAMS"

std::random_device RD;
std::mt19937 MT(RD()); 

int main()
{
    
  /* 1) Set values for number of degrees of freedom "DoF", bending rigidity 
     "rig", internal tension "sig", frame tension "tau" and lattice spacing
     "alpha".Define and initialize variables needed for the MC code.        */
  
  int N = sqrt(DoF);        //degrees of freedom per dimension
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
  
  int x,y;
  int lattice_changes = 0;
  int accepted_moves = 0;
  int lattice_moves = 0;
  int move_counter = 0;         
  
  int len_area = 5; 
  int len_corr = 4*nghost+1;
  int len_ener = pow((2*nghost+1),2); 
  Site neighbors_area[len_area] = {};
  Site neighbors_corr[len_corr] = {};
  Site neighbors_ener[len_ener] = {};
  Site site;
  
  std::uniform_int_distribution<int> RandInt(0,N-1);  
  std::uniform_real_distribution<double>  RandDouble(-epsilon,epsilon);
  
  std::string filename = "sampling.txt";

  /* 2) Initialize the height field hfield(i,j)                           */ 
  RectMesh hfield(N,N,nghost);
  
  /* 3) Calculate the projected membrane area "prj_area", the total area 
     "tot_area" and the energies "tau_energy","sig_energy","crv_energy",
     "cor_energy" and "tot_energy".                                       */

  CalculateTotal(hfield,tot_energy,tau_energy,crv_energy,
  		 sig_energy,cor_energy,tot_area,prj_area,
  		 DoF,alpha,rig,sig,tau);
  
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
      accept  = Metropolis(dElocal);
      
      /* 9) If the move is accepted, update total area and energy.
  	 Otherwise return to previous state.                              */

      AcceptOrDecline(hfield,site,accept,where,tot_area,
		      tot_energy,dAlocal,dElocal,accepted_moves,perturb);
      
      /* 10) After 5 MC steps, randomly change alpha, compute the 
	 new projected area and update the total energy.                  */

      ChangeLattice(hfield,min_change,max_change,
		    move_counter,alpha,prj_area,
      		    tot_area,DoF,tot_energy,
      		    tau_energy,crv_energy,
      		    sig_energy,cor_energy,
      		    rig,sig,tau,lattice_moves,
      		    lattice_changes);
      
      /* 11) Sample                                                       */

      Sample(iter,filename,tot_energy,tau_energy,crv_energy,sig_energy,
	     cor_energy,tot_area,prj_area,alpha,DoF);
    }

  /* 12) Write the final configuration and print acceptance ratios        */
  
  hfield.write("hfield.h5");
  PrintAcceptance(maxiter,accepted_moves,lattice_moves,lattice_changes); 

  return 0;
}

