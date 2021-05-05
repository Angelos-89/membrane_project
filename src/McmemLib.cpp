/* In the following, when paper1 is mentioned, we refer to the  paper 
   "Monte Carlo study of the frame fluctuation and internal tensions of 
   fluctuating membranes with fixed area" by Shiba et al. (2016). */

#include <stdio.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>
#include <string>
#include "Vec3dLib.hpp"
#include "McmemLib.hpp"
using std::strtol;


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

std::random_device rd;
std::mt19937 mt(rd());

static double ShiftInEnergy = 0.;

/*------------------------- OutputParams ---------------------------*/

/* Prints the parameters of the specific run to the 
   console and writes them to a .txt file                           */

void OutputParams(const int maxiter,const int N,const int DoF,
		  const int nghost,const double rig,const double sig,
		  const double tau,const double epsilon,
		  const double min_change,const double max_change,
		  double alpha,double pn_prcn,int sample_every,
		  int rank,double Eactive)
{
  std::stringstream strm; 
  strm << "Run "                              << rank+1
       << " parameters: \n-----------------\n"
       << "Maxiter: "                         << maxiter                  <<"\n"
       << "Grid Size: "                       << N << "x"  << N           <<"\n"
       << "DoF: "                             << DoF                      <<"\n"
       << "Ghost points per boundary point: " << nghost                   <<"\n"
       << "Bending rigidity: "                << rig         << " (k_BT)" <<"\n"
       << "Internal tension: "                << sig         << " (k_BT)" <<"\n"
       << "Frame tension: "                   << tau         << " (k_BT)" <<"\n"
       << "Max height perturbation: "         << epsilon     << " (a_0)"  <<"\n"
       << "Min change in lattice spacing: "   << min_change               <<"\n"
       << "Max change in lattice spacing: "   << max_change               <<"\n"
       << "Initial lattice spacing: "         << alpha       << " (a_0)"  <<"\n"
       << "Pinning percentage: "              << pn_prcn*100 << "%"       <<"\n"
       << "Sample every: "                    << sample_every<< " moves"  <<"\n"
       << "Eactive: "                         << Eactive     << " (k_BT) " 
       << "\n------------------------------------\n\n";
  std::cout << strm.str();

  /* Now write to a file */
  std::ofstream file;
  std::string filename = "PARAMS_" + std::to_string(rank) + ".txt";
  file.open(filename);
  file << strm.str();
  file.close();
}

/*------------------------------ TileOfNeighbors ----------------------*/

/* Returns a set of neighbors of a site, which form a square of a certain 
   length. The neighbors are stored inside neighbors[] and the site
   is included. */

void TileOfNeighbors(Site neighbors[], RectMesh& field, Site site, int radius)
{
  int cols = field.getcols();
  int rows = field.getrows();
  int nghost = field.getnghost();
  int x = site.getx();
  int y = site.gety();
  int k=0;
  for(int j=y-radius; j<=y+radius; j++)
    {
      for(int i=x-radius; i<=x+radius; i++)
	{
	  site.set(i, j);
	  neighbors[k] = site;
	  k++;
	}
    } 
}

/*---------------------------- InitPinning -------------------------*/

/* Given the fraction (0.0 - 1.0) of the total DoFs that will be 
   pinned, it fills a set with all the pinned sites. N = sqrt(DoFs) */

std::unordered_set<Site> InitPinning(RectMesh& hfield,double pin_ratio,
				     Site neighbors[],int radius) 
{
  if (pin_ratio > 1.0 or pin_ratio < 0){
    std::cout << "InitPinning: Pinning fraction must be in the range "
      "0.0 <= pin_fraction <= 1.0. Exiting." << std::endl;
    exit(EXIT_FAILURE);}
  
  int Nx = hfield.getcols();
  int Ny = hfield.getrows();
  int N2 = Nx*Ny;
  int LEN = N2*pin_ratio;

  Site site_a,site_b;
  std::unordered_set<Site> pinned_sites = {}; 
  int length = pow( (2*radius+1) ,2);
  int x,y;
  while (pinned_sites.size() < LEN)
    {
      std::uniform_int_distribution<int> RandIntX(0,Nx-1);
      std::uniform_int_distribution<int> RandIntY(0,Ny-1);
      x = RandIntX(mt); y = RandIntY(mt);
      site_a.set(x,y); pinned_sites.insert(site_a);

      if (pinned_sites.find(site_a) != pinned_sites.end()) //if site_a does not already exist 
	{
	  TileOfNeighbors(neighbors,hfield,site_a,radius); //find the neighbors
	  for (int k=0; k<length; k++) // find the distance between site_a and all neighbors 
	    {
	      site_b = neighbors[k];
	      if (site_a.dist(site_b) <= radius) //if in correct range, insert
		{
		  site_b.set( (site_b.getx()+Nx)%Nx , (site_b.gety()+Ny)%Ny);
		  pinned_sites.insert(site_b);
		}
	    }//end for
	}//end if
    }//end while

  while (pinned_sites.size() > LEN) //erase some sites to get the correct percentage
    pinned_sites.erase(pinned_sites.begin());
  return pinned_sites;  
}

/*------------------------ InitSurface ---------------------*/

/* Initializes a RectMesh object with random values ranging
   from min to max. Furthermore, it sets the sites contained
   inside the vector to zero.                               */

void InitSurface(RectMesh& hfield,std::unordered_set<Site>& pinned_sites,
		 double min,double max)
{
  std::uniform_real_distribution<double> UnifProb(min,max);
  for (int j=0; j<hfield.getrows(); j++)
  {
    for (int i=0; i<hfield.getcols(); i++)
      hfield(i,j) = UnifProb(mt);
  }
  
  // /* Pinning */
  // int x,y;
  // for (auto it = pinned_sites.begin(); it != pinned_sites.end(); ++it)
  //   {
  //     x = (*it).getx();
  //     y = (*it).gety();
  //     hfield(x,y) = h0;
  //   }
  
  GhostCopy(hfield);
}

/*---------------- GhostCopy ------------------*/

/*This function implements periodic boundary
  conditions by copying the first and last nghost 
  rows and columns of a 2D mesh at the ghost 
  zones of the same  mesh. nghost denotes the 
  number of ghost points per boundary point.

  Example:

  A 2D mesh defined as: RectMesh field(3,4,2);
  with 3 columns, 4 rows and 2 ghost points per
  boundary point, is shown below before the 
  execution of GhostCopy.

  0    0    0    0    0    0    0    
  0    0    0    0    0    0    0    
  0    0    0    1    2    0    0    
  0    0    6    7    8    0    0    
  0    0    3    4    5    0    0    
  0    0    0    1    2    0    0    
  0    0    0    0    0    0    0    
  0    0    0    0    0    0    0    

  After executing GhostCopy(field); 

  4    5    3    4    5    3    4    
  1    2    0    1    2    0    1    
  1    2    0    1    2    0    1    
  7    8    6    7    8    6    7    
  4    5    3    4    5    3    4    
  1    2    0    1    2    0    1    
  1    2    0    1    2    0    1    
  7    8    6    7    8    6    7     */

void GhostCopy(RectMesh& mesh)
{
  int nghost = mesh.getnghost();
  int Nx = mesh.getcols();
  int Ny = mesh.getrows();
  for (int j=-nghost; j<Ny+nghost; j++)
    {
      for (int i=0; i<nghost; i++)
        {
	  mesh(Nx+i,j) = mesh(i,j);
	  mesh(-i-1,j) = mesh(Nx-i-1,j);
        }
    }
  for (int i=-nghost; i<Nx+nghost; i++)
    {
      for (int j=0; j<nghost; j++)
        {
	  mesh(i,Ny+j) = mesh(i,j);
	  mesh(i,-j-1) = mesh(i,Ny-j-1);
        }
    }
}

/*-------------------------  Der ------------------------------*/

/* This function calculates the gradient of a scalar field F(i,j) 
   at a given point, where i corresponds to the x-coordinate and 
   j to the y-coordinate. The pointer grad represents an 1D array 
   where grad[0] stores the partial derivative with respect to x 
   and grad[1] the partial derivative with respect to y. 
   Central finite difference schemes are used to compute the 
   derivatives. alpha is the lattice spacing.                  */

void Der(const RectMesh& mesh, Site site, double alpha, double grad[])
{

  int i = site.getx();
  int j = site.gety();
  int order = 2*mesh.getnghost();

  switch (order)
    {
    case 2:
      
      grad[0] = (mesh(i+1,j) - mesh(i-1,j)) / (2*alpha);
      grad[1] = (mesh(i,j+1) - mesh(i,j-1)) / (2*alpha);
      break;
      
    case 4:
      
      grad[0] = ( -mesh(i+2,j) + 8*mesh(i+1,j) - 8*mesh(i-1,j) + mesh(i-2,j) )
	/ (12*alpha);
      grad[1] = ( -mesh(i,j+2) + 8*mesh(i,j+1) - 8*mesh(i,j-1) + mesh(i,j-2) )
	/ (12*alpha);
      break;
      
    default:
      std::cout << "Enter a valid scheme order: nghost=1 (2nd order)"
	"or nghost=2 (4th order)"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
  
}

/*-------------------------- Der2 ------------------------*/

/* This function calculates the second order derivatives of 
   a scalar field F(i,j) at a given point where i corresponds 
   to the x-coordinate and j to the y-coordinate. The pointer 
   hess represents the Hessian matrix stored as an 1D array 
   where hess[0] = Fxx, hess[1] =  Fxy, hess[2] = Fyx 
   and hess[3] = Fyy. Central finite difference schemes 
   are used to compute the derivatives. 
   alpha is the lattice spacing.                          */

void Der2(const RectMesh& mesh, Site site, double alpha, double hess[])
{
  
  int i = site.getx();
  int j = site.gety();
  int order = 2*mesh.getnghost();
  
  switch(order)
    {
    case 2:
      
      hess[0] = (mesh(i+1,j)+mesh(i-1,j)-2*mesh(i,j)) / (alpha*alpha);
      hess[1] = (mesh(i+1,j+1)-mesh(i-1,j+1)-mesh(i+1,j-1)+mesh(i-1,j-1))
	/(4*alpha*alpha);
      hess[2] = hess[1];
      hess[3] = (mesh(i,j+1)+mesh(i,j-1)-2*mesh(i,j)) / (alpha*alpha);
      break;
      
    case 4:
      
      hess[0] = (-mesh(i+2,j)+16*mesh(i+1,j)-30*mesh(i,j)+16*mesh(i-1,j)
		 - mesh(i-2,j)) / (12*alpha*alpha);
      hess[1] = (+8*(mesh(i+1,j-2)+mesh(i+2,j-1)+mesh(i-2,j+1)+mesh(i-1,j+2))
		 -8*(mesh(i-1,j-2)+mesh(i-2,j-1)+mesh(i+1,j+2)+mesh(i+2,j+1))
		 -1*(mesh(i+2,j-2)+mesh(i-2,j+2)-mesh(i-2,j-2)-mesh(i+2,j+2))
		 +64*(mesh(i-1,j-1)+mesh(i+1,j+1)-mesh(i+1,j-1)-mesh(i-1,j+1)))
	/ (144*alpha*alpha);
      hess[2] = hess[1];
      hess[3] = (-mesh(i,j+2)+16*mesh(i,j+1)
		 -30*mesh(i,j)+16*mesh(i,j-1)
		 -mesh(i,j-2)) / (12*alpha*alpha);
      break;
      
    default:
      std::cout << "Enter a valid scheme order: nghost=1 (2nd order)"
	"or nghost=2 (4th order)"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
}

/*------------------------- SNodeArea ----------------------------*/

/* This function calculates the area corresponding to a single node
   represented by (i,j). We use equation (7) from paper1. 
   That is, we  calculate the area of the four triangles that are 
   formed by the value of the field at (i,j) and its four nearest 
   neighbors and then devide by 2. Note however, that we do not 
   multiply by 0.5 to find the area of the triangle, but multiply 
   by 0.25 at the end instead.                                    */

double SNodeArea(const RectMesh& field,Site site,double alpha)
{
  int i = site.getx();
  int j = site.gety();
  
  Vec3d a( alpha,  0    , field(i+1,j) - field(i,j));
  Vec3d b( 0    ,  alpha, field(i,j+1) - field(i,j));
  Vec3d c(-alpha,  0    , field(i-1,j) - field(i,j));
  Vec3d d( 0    ,  alpha, field(i,j-1) - field(i,j));

  double area1 = Norm(Cross(a,b)); 
  double area2 = Norm(Cross(b,c));
  double area3 = Norm(Cross(c,d));
  double area4 = Norm(Cross(d,a));

  return 0.25*(area1+area2+area3+area4);
}

/*------------------------ LocalArea ---------------------------*/

/* This function calculates the area that corresponds to a site
   (i,j) and its four nearest neighbors as well.                */

double LocalArea(const RectMesh& hfield,Site neighbors[],double alpha)
{
  double area = 0;
  for (int n=0; n<5; n++)
    area += SNodeArea(hfield,neighbors[n],alpha);
  return area;  
}

/*--------------------------- TotalArea --------------------------*/

/* This function calculates the total area of the surface that 
   is formed by a 2D scalar field. The scalar field is represented 
   by a RectMesh object and it is assumed to be periodic in both 
   x and y. The step size in each direction is alpha.             */

double TotalArea(const RectMesh& field,double alpha)
{
  // if (field.getnghost() == 0)
  //   {
  //     std::cout << "TotalArea: The number of ghost points per side of the"
  // 	"RectMesh object must be at least 1! Exiting. \n"
  // 		<< std::endl;
  //     exit(EXIT_FAILURE);
  //   }
  double total_area=0;
  Site site;
  for (int j=0; j<field.getrows(); j++)
    {
      for (int i=0; i<field.getcols(); i++)
	{
	  site.set(i,j);
	  total_area += SNodeArea(field,site,alpha);
	}
    }
  return total_area;
}

/*--------------------- SNodeNormalZ ----------------------*/

/* This function calculates the z coordinate of the unit 
   normal vector at a single point of the surface that is 
   represented by the RectMesh object hfield. alpha is the 
   lattice spacing and site represents the point of the 
   mesh at which the unit normal coordinate is calculated. */

double SNodeNormalZ(const RectMesh& hfield,Site site,double alpha)
{
  double grad[2]={};
  Der(hfield,site,alpha,grad);
  double h_x = grad[0];
  double h_y = grad[1];
  return (1. / sqrt(1. + pow(h_x,2) + pow(h_y,2)));
};

/*------------------- NormalZ ------------------*/

/* This function calculates the z coordinates 
   of the unit normal vectors at every point of 
   the surface represented by hfield. alpha is 
   the lattice spacing.                         */

RectMesh NormalZ(const RectMesh& hfield,double& alpha)
{
  Site site;
  RectMesh normalz(hfield.getcols(),hfield.getrows(),0);
  for(int j=0; j<hfield.getrows(); j++)
    {
      for(int i=0; i<hfield.getcols(); i++)
	{
	  site.set(i,j); 
	  normalz(i,j) = SNodeNormalZ(hfield,site,alpha);
	}
    }
  return normalz;
};

/*------------------------ SNodeCurvature ------------------------------*/

/* This function calculates the local mean curvature, i.e., the mean 
   curvature H at a given point (i,j) in the Monge representaion. Here, 
   hfield is the surface represented by a RectMesh object. Note that the 
   GhostCopy function must be called before calling this function for 
   boundary points.                                                     */

double SNodeCurvature(const RectMesh& h,Site site,double alpha)
{
  double grad[2]={};
  Der(h,site,alpha,grad);
  double hess[4] = {};
  Der2(h,site,alpha,hess);
  double h_x = grad[0];
  double h_y = grad[1];
  double h_xx = hess[0];
  double h_xy = hess[1];
  double h_yx = hess[2];
  double h_yy = hess[3];
  double h_x2 = h_x*h_x;
  double h_y2 = h_y*h_y;
  double numerator = (1.+h_x2)*h_yy + (1.+h_y2)*h_xx - 2.*h_x*h_y*h_xy;
  double denominator = pow(1.0+h_x2+h_y2, 3./2.);
  double curvature = numerator/denominator;
  return curvature;
}

/*----------------------- TotalCurvature --------------------------*/

/* This function returns a RectMesh object with the mean curvature 
   of every point of a surface that is represented by the RectMesh 
   object field.                                                   */

// RectMesh TotalCurvature(const RectMesh& field,double& alpha)
// {
//   Site site;
//   RectMesh H(field.getcols(),field.getrows(),0);
//   for (int j=0; j<field.getrows(); j++)
//     {
//       for (int i=0; i<field.getcols(); i++)
// 	{
// 	  site.set(i,j);
// 	  H(i,j) = SNodeCurvature(field,site,alpha);
// 	}
//     }
//   return H;
// }

/*----------------------- SNodeCurvatureEnergy ------------------------*/

/* This function calculates the curvature energy of a single site 
   due to its mean curvature H. This energy is the integrand of 
   equation (3) found in paper1. Here, rig is the bending 
   rigidity and hfield the RectMesh instance representing the membrane 
   in the Monge gauge.                                                 */

double SNodeCurvatureEnergy(const RectMesh& hfield,Site site,
			    double alpha,const double rig)
{
  double dA = SNodeArea(hfield,site,alpha);
  double H = SNodeCurvature(hfield,site,alpha);
  return 0.5*rig*H*H*dA; 
}

/*---------------- LocalCurvatureEnergy -----------------------*/

double LocalCurvatureEnergy(const RectMesh& field,Site neighbors[],
			    double alpha,const double rig)
{
  int len = pow((2*field.getnghost()+1),2);
  double energy = 0;
  for (int n=0; n<len; n++)  
    energy += SNodeCurvatureEnergy(field,neighbors[n],alpha,rig);  
  return energy;
}

/*---------------------- CurvatureEnergyTotal --------------------*/

/* This function calculates the total energy of the membrane due
   to its curvature by summing up the contributions from every 
   site. Again, hfield is the RectMesh object representing the
   membrane height field and rig is the bending rigidity.         */

double CurvatureEnergyTotal(const RectMesh& hfield,
			    double alpha,const double rig)
{
  double energy = 0;
  Site site;
  for(int j=0; j<hfield.getrows(); j++)   
    {
      for(int i=0; i<hfield.getcols(); i++){
	site.set(i,j);
	energy += SNodeCurvatureEnergy(hfield,site,alpha,rig);
      }
    }
  return energy;
}

/*-------------------- SNodeCorrectionEnergy ----------------------*/

/* This function calculates the correction in energy at a single 
   node, due to the Monge gauge representation of the membrane. 
   See paper1.                                                     */

double SNodeCorrectionEnergy(const RectMesh& hfield,Site site,double alpha)
{
  return -log(SNodeNormalZ(hfield,site,alpha));
}

/*-------------------- LocalCorrectionEnergy ----------------------*/

/* This function calculates the correction energy by summing up the
   contributions from every neighbor that affects the energy when
   the site (i,j) is perturbed. This depends on the order of
   the derivative scheme. The neighbors form a cross with (i,j) in
   the center.                                                     */

double LocalCorrectionEnergy(const RectMesh& field,
			     Site neighbors[],double alpha)
{
  int len = 4*field.getnghost()+1;
  double energy = 0;
  for (int n=0; n<len; n++)  
    energy += SNodeCorrectionEnergy(field,neighbors[n],alpha);  
  return energy;
}

/*----------------------------- LocalPinEnergy --------------------------------*/

double LocalPinEnergy(const RectMesh& hfield,Site& site,
		      const double& pot_strength,const double& h0)
{
  int x = site.getx();
  int y = site.gety();
  double h = hfield(x,y);
  return 0.5 * pot_strength * pow((h-h0),2);
}
  
/*-------------------------CorrectionEnergyTotal------------------------------*/

/* This function calculates the total correction in energy due to the Monge
   gauge representation of the membrane. For details see paper1. The correction 
   term for the energy is given by H_corr = -k_b*T*sum_{i,j}{ln(NormalZ(i,j))}, 
   where k_b is Boltzmann's constant, T the temperature and NormalZ the z 
   coordinates of the unit normal vectors along the surface of the membrane. 
   Note that we work in reduced units i.e., k_b*T=1.                          */

double CorrectionEnergyTotal(const RectMesh& hfield,double alpha)
{
  Site site;
  double sum = 0.0;
  for(int j=0; j<hfield.getrows(); j++)
    {
      for(int i=0; i<hfield.getcols(); i++){
	site.set(i,j);
	sum += log(SNodeNormalZ(hfield,site,alpha));}
    }
  return -sum;

  // RectMesh temp = NormalZ(hfield,alpha);
  // temp.ln();
  // return -temp.sum();
}

/*------------------------- PinningEnergyTotal --------------------------*/

/* Calculates the total energy due to the harmonic oscillator potential
   at each pinned site.                                                  */

double PinningEnergyTotal(const RectMesh& hfield,
			  std::unordered_set<Site>& pinned_sites,
			  const double& pot_strength,const double& h0)
{
  double energy = 0;
  int x,y;
  double h;
  for (auto it = pinned_sites.begin(); it != pinned_sites.end(); it++)
    {
      x = (*it).getx();
      y = (*it).gety();
      h = hfield(x,y);
      energy += pow((h-h0),2);
    }
  energy *= 0.5*pot_strength;
  return energy;
}

/*-------------------------- LocalEnergy -------------------------------------*/

/* This function calculates the local energy by summing up the contributions
   of a block of nearest neighbors. The size of the block depends on the order
   of the finite difference scheme. Contributions from all different energies
   are taken into account.                                                    */

double LocalEnergy(const RectMesh& hfield,
		   Site neighbors_area[],
		   Site neighbors_corr[],
		   Site neighbors_ener[],
		   double alpha,const double rig,
		   const double sig,const double tau,
		   const double& pot_strength,const double& h0,bool pin)
{
  double tau_ener = -tau*alpha*alpha; 
  double crv_ener = LocalCurvatureEnergy(hfield,neighbors_ener,alpha,rig);
  double sig_ener = sig*LocalArea(hfield,neighbors_area,alpha);
  double cor_ener = LocalCorrectionEnergy(hfield,neighbors_corr,alpha);
  double pin_energy = 0;
  if (pin){
      Site site; site = neighbors_area[0];
      pin_energy = LocalPinEnergy(hfield,site,pot_strength,h0);}
  return tau_ener + crv_ener + sig_ener + cor_ener + pin_energy;
}

/*-------------------------------- CalculateTotal --------------------------*/

/* This function updates the total energy, total area and total projected
   area of the membrane. It also updates the different energies involved.   */ 

void CalculateTotal(const RectMesh& hfield,const double& rig,
		    const double& sig,const double& tau,double& tot_energy,
		    double& tau_energy,double& crv_energy,double& sig_energy,
		    double& cor_energy,double& pin_energy,double& tot_area,
		    double& prj_area,double& alpha,
		    std::unordered_set<Site>& pinned_sites,
		    const double& pot_strength,const double& h0)
{
  int DoFs = hfield.getcols()*hfield.getrows();
  prj_area = (double)DoFs*alpha*alpha;
  tot_area = TotalArea(hfield,alpha);
  tau_energy = -tau*prj_area;
  crv_energy = CurvatureEnergyTotal(hfield,alpha,rig);
  sig_energy = sig*tot_area;
  cor_energy = CorrectionEnergyTotal(hfield,alpha);
  pin_energy = PinningEnergyTotal(hfield,pinned_sites,pot_strength,h0);
  tot_energy = tau_energy + crv_energy + sig_energy + cor_energy + pin_energy;
}

/*---------------------- WhereIs ---------------------*/

/* This function takes as arguments a site(i,j), the 
   columns, rows and ghost points of a RectMesh object 
   and returns an integer value that depends on 
   whether it belongs to the boundary or the bulk.    */

bool WhereIs(Site site, int cols, int rows, int nghost)
{
  int i = site.getx();
  int j = site.gety();
  // if (i>=cols || i<0 || j>=rows || j<0)
  // {
  //   std::cout << "WhereIs: It must hold that 0<=i<cols and 0<=j<rows."
  //     "Exiting." << std::endl;
  //   exit(EXIT_FAILURE);
  // }
  if (i<nghost || i>=cols-nghost || j<nghost || j>=rows-nghost)
    return 1; //boundary point
  else return 0; //bulk 
}

/*------------------------- Ispinned ------------------------*/

bool Ispinned(Site& site,std::unordered_set<Site>& pinned_sites)
{
  auto search = pinned_sites.find(site);
  if (search != pinned_sites.end()) 
    return 1;
  else 
    return 0;
}

/*------------------ GetNeighbors ---------------*/

/* This function fills up three lists that 
   contain the neighbors of a point (i,j) and
   the point as well. These lists are needed
   for the calculation of local area, local 
   correcrion energy and local curvature energy. */

void GetNeighbors(const RectMesh& field,Site site,
		  Site neighbors_area[],
		  Site neighbors_corr[],
		  Site neighbors_ener[])
{
  Site init = site;
  int cols = field.getcols();
  int rows = field.getrows();
  int nghost = field.getnghost();
  int x = site.getx();
  int y = site.gety();
  
  /* Fill the list of neighbors needed for 
     the local area calculation. */
  
  neighbors_area[0] = site;
  
  site.set((x+1+cols)%cols, y);
  neighbors_area[1] = site;
  
  site.set(x, (y+1+rows)%rows);
  neighbors_area[2] = site;
  
  site.set((x-1+cols)%cols, y);
  neighbors_area[3] = site;

  site.set(x, (y-1+rows)%rows);
  neighbors_area[4] = site;

  /* Fill the list needed for 
     the local correction energy.
     The neighbors are located in
     a cross with (x,y) 
     in the center. */
  
  site = init;
  neighbors_corr[0] = site;
  int k=1;
  for(int i=-nghost; i<0; i++)
    {
      site.set((i+x+cols)%cols, y);
      neighbors_corr[k] = site;
      k++;
    }
  for(int i=1; i<=nghost; i++)
    {
      site.set((i+x+cols)%cols, y);
      neighbors_corr[k] = site;
      k++;
    }
  for(int j=-nghost; j<0; j++)
    {
      site.set(x, (j+y+rows)%rows);
      neighbors_corr[k] = site;
      k++;
    }
  for(int j=1; j<=nghost; j++)
    {
      site.set(x, (j+y+rows)%rows);
      neighbors_corr[k] = site;
      k++;
    }
  
  /* Fill the list needed for 
     the local curvature energy.
     The neighbors form a block (tile) 
     consisting of (2*nghost+1)^2
     sites. */ 

  k=0;
  for(int j=y-nghost; j<=y+nghost; j++)
    {
      for(int i=x-nghost; i<=x+nghost; i++)
	{
	  site.set((cols+i)%cols, (rows+j)%rows);
	  neighbors_ener[k] = site;
	  k++;
	}
    }
}

/*------------ PrintNeighbors -------------*/

void PrintNeighbors(Site neighbors[],int len)
{
  for(int i=0; i<len; i++)
    neighbors[i].print();
}

/*------------- AddShift ----------------*/

void AddShift(double& dE)
{
  ShiftInEnergy = dE;
}

/*------------------- Metropolis ----------------------*/

/* It returns 1 if the move is accepted and 
   0 otherwise according to the Boltzmann criterion.   */

bool Metropolis(double& dElocal)
{
  double dE = dElocal + ShiftInEnergy;
  //dElocal += ShiftInEnergy;
  if(dE < 0) return 1;
  else
    {
      std::uniform_real_distribution<double> UnifProb(0,1);
      double r = UnifProb(mt);
      if (r < exp(-dE)) return 1;
      else return 0;
    }
}

/*------------------------- UpdateState -----------------------------*/

/* This function updates the total energy and total area of the membrane
   if a height trial move is accepted. Otherwise it returns the membrane
   to its previous state.                                                */

void UpdateState(RectMesh& hfield,Site site,bool accept,bool where,
		     double& tot_area,double& tot_energy,double& dAlocal,
		     double& dElocal,int& height_changes,double& perturb)
{  
  if (accept == 1)
    {
      height_changes ++;
      tot_energy += dElocal;
      tot_area += dAlocal;
    }
  else
    {
      int i = site.getx(), j = site.gety();
      hfield(i,j) -= perturb;
      if (where == 1) GhostCopy(hfield);
    }
}

/*-------------------------- Wstats ------------------------*/

/* It prints the acceptance rations of both height and lattice size 
   trial moves, the number of times the spectrum was averaged, and 
   also writes them to file. */

void WStats(const int maxiter, int height_changes,
	    int lattice_moves, int lattice_changes,
	    int spec_steps,int rank)
{
  double accept_ratio = (double) height_changes/maxiter;
  double lattice_moves_ratio = (double) lattice_changes/lattice_moves; 
  std::stringstream stream; 
  stream << "Simulation " << rank+1
	 << " is finished: \n-------------------------\n"
	 << "Height move ratio: " << accept_ratio << "\n"
	 << "Lattice spacing move ratio: "
	 << lattice_moves_ratio << "\n"
	 << "Power spectrum averaged " << spec_steps << " times"<< "\n"
	 << "--------------------------------------------\n";
  std::cout << stream.str();

  /*Now damp to a file*/ 
  std::ofstream file;
  std::string filename = "stats_" + std::to_string(rank) + ".txt";
  file.open(filename);
  file << stream.str();
  file.close();  
}

/*------------------------- ChangeLattice ---------------------------*/

/* It attemps a lattice size trial move and checks whether it is 
   accepted or not. If it is, it updates everything, if not, it returns 
   the membrane to its previous state.                                  */

bool ChangeLattice(const RectMesh& hfield,const double& min_change,
		   const double& max_change,const double& rig,
		   const double& sig, const double& tau,double& prj_area,
		   double& tot_area,double& tot_energy,double& tau_energy,
		   double& crv_energy,double& sig_energy,double& cor_energy,
		   double& pin_energy,double& alpha,
		   int& lattice_moves,int& lattice_changes,
		   std::unordered_set<Site>& pinned_sites,
		   const double& pot_strength,const double& h0)
{
  lattice_moves ++;
  std::uniform_real_distribution<double> RandDouble(min_change,max_change); 
  double old_prj_area = prj_area;
  double old_tot_area = tot_area;
  double old_energy = tot_energy;
  double old_alpha = alpha;
  double fraction = RandDouble(mt);
  alpha *= fraction;
  
  CalculateTotal(hfield, rig, sig, tau, tot_energy, tau_energy, crv_energy,
		 sig_energy, cor_energy, pin_energy, tot_area, prj_area, alpha,
		 pinned_sites, pot_strength, h0);

  double dE = tot_energy-old_energy; 
  bool accept = Metropolis(dE);

  if (accept == 1){
    lattice_changes ++;
    return 1;}
  else{
    tot_energy = old_energy;
    prj_area = old_prj_area;
    tot_area = old_tot_area;
    alpha = old_alpha;
    return 0;}
}

/*--------------------------- Sample -------------------------*/

/* Stores the data in a txt file. */

int first_call = 1;
void WriteData(std::string filename,int& iter,int& total_moves,
	    double& tot_energy,double& crv_energy,
	    double& cor_energy,double& pin_energy,double& tot_area,
	    double& prj_area,double& alpha,int& DoFs)
{
  std::ofstream file;
  file.open(filename, std::ios::app);
  if (first_call == 1)
    {
      file << "%iter"                                     << "\t"
	   << std::right<<std::setw(12)<<"total_moves"    << "\t"
	   << std::right<<std::setw(12)<<"total_area"     << "\t"
	   << std::right<<std::setw(12)<<"prj_area"       << "\t"
	   << std::right<<std::setw(12)<<"alpha"          << "\t"
	   << std::right<<std::setw(15)<<"curv_energy"    << "\t"
	   << std::right<<std::setw(15)<<"entropic_corr"  << "\t"
	   << std::right<<std::setw(15)<<"pinning_energy" << "\t"
	   << std::right<<std::setw(15)<<"tot_energy"     << "\n";
      first_call = 0;
    }
  double DOF = (double) DoFs;
  file << iter                                                           << "\t"
       <<std::right<<std::setw(12)<<total_moves                          << "\t"
       <<std::right<<std::setw(12)<<std::setprecision(6)<< tot_area/DOF  << "\t"
       <<std::right<<std::setw(12)<<std::setprecision(6)<< prj_area/DOF  << "\t"
       <<std::right<<std::setw(12)<<std::setprecision(6)<< alpha         << "\t"
       <<std::right<<std::setw(15)<<std::setprecision(6)<< crv_energy/DOF<< "\t"
       <<std::right<<std::setw(15)<<std::setprecision(6)<< cor_energy/DOF<< "\t"
       <<std::right<<std::setw(15)<<std::setprecision(6)<< pin_energy/DOF<< "\t"
       <<std::right<<std::setw(15)<<std::setprecision(6)<< tot_energy/DOF<< "\n";
  file.close();
}

void Sample(std::string output_filename, int& sample_every, int& iter,
	    int& total_moves, double& tot_energy, double& crv_energy,
	    double& cor_energy, double& pin_energy, double& tot_area,
	    double& prj_area, double& alpha, int& DoFs)
{
  if (total_moves % sample_every == 0)  
    WriteData(output_filename, iter, total_moves, tot_energy, crv_energy,
	   cor_energy, pin_energy, tot_area, prj_area, alpha, DoFs);      
}




/*-------------------------------- ReadInput -----------------------------*/

void ReadInput(std::string filename,int& sim, int& acc_samples, double& maxiter,
	       double& sig,double& tau,double& epsilon,double& min_change,
	       double& max_change,double& pin_ratio, double& Ea)
{
  std::ifstream infile;
  infile.open(filename);
  if(!infile.is_open())
    {
      std::cout << "File "<< filename <<  " is not open. Exiting." << std::endl;
      exit(EXIT_FAILURE);
    }
  while(!infile.eof())
    infile >> sim >> acc_samples >> maxiter >> sig >> tau >> epsilon >> min_change
	   >> max_change >> pin_ratio >> Ea;

  infile.close();
}

/*-------------------write_to_extendible_H5----------------------*/

void Write_to_extendible_H5(const char* FILENAME, RectMesh& hfield)
{
  hsize_t ndims = 2;
  hsize_t nrows = hfield.getrows() + 2*hfield.getnghost();
  hsize_t ncols = hfield.getcols() + 2*hfield.getnghost();

  /* Create a memory dataspace to indicate the 
     size of our buffer to be written in memory. 
     The dimensions of the buffer do not change 
     during code execution. */

  hsize_t mbuff_dims[ndims];
  mbuff_dims[0] = nrows;
  mbuff_dims[1] = ncols;
  hid_t mem_space = H5Screate_simple(ndims, mbuff_dims, NULL);
  
  /* Open the file. */
  
  hid_t file = H5Fopen(FILENAME, H5F_ACC_RDWR, H5P_DEFAULT);

  if (file<0)
    {
      std::cout << "HDF5 file to contain extendible dataset does not exist. Exiting. " << std::endl;
      exit(EXIT_FAILURE);
    }

  else
    {
  
      /* Check if there is a dataset. */
  
      if ( !H5Lexists(file,"extendibleDset",H5P_DEFAULT) )
	{
	  /* Dataset does not exist. Create it and
	     write the first buffer. */
	  
	  // Create a 2D dataspace.
	  
	  hsize_t dims[ndims];
	  dims[0]=nrows; dims[1]=ncols;
	  hsize_t max_dims[ndims];
	  max_dims[0]=H5S_UNLIMITED; max_dims[1]=ncols;
	  hid_t file_space = H5Screate_simple(ndims, dims, max_dims);
	  
	  // Then create a dataset creation property list.  
	  
	  hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
	  H5Pset_layout(plist, H5D_CHUNKED);
	  hsize_t chunk_dims[ndims];
	  chunk_dims[0]=nrows; chunk_dims[1]=ncols;
	  H5Pset_chunk(plist, ndims, chunk_dims);
	  
	  // Create the dataset.
	  
	  hid_t dset = H5Dcreate(file, "extendibleDset", H5T_NATIVE_DOUBLE,
				 file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
	  
	  /* Close resources. */
	  
	  H5Pclose(plist);
	  	  
	  /* Write the first buffer */
	  
	  // Select hyperslab on file dataset.
	  
	  file_space = H5Dget_space(dset);
	  hsize_t start[2] = {0, 0};
	  hsize_t count[2] = {nrows, ncols};
	  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start,
			      NULL, count, NULL);
	  
	  /* Write buffer to dataset. */
	  
	  H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, file_space,
		   H5P_DEFAULT, hfield.getmemory());
	  
	  /*We can now release resources. */
	  
	  H5Sclose(mem_space);
	  H5Dclose(dset);
	  H5Sclose(file_space);
	  H5Fclose(file);
	}
      else
	{
	  /* Dataset already exists. Extend it and write 
	     the next buffer. */
	  
	  // Open the dataset and get the dimensions of the existing dataset.
	  
	  hid_t dset = H5Dopen(file,"extendibleDset",H5P_DEFAULT);
	  hid_t file_space = H5Dget_space(dset);
	  hsize_t dims[ndims];
	  H5Sget_simple_extent_dims(file_space, dims, NULL);
	  
	  // Extend the dimensions.
	  
	  dims[0] += nrows;      
	  dims[1] = ncols;
	  H5Dset_extent(dset, dims);
	  file_space = H5Dget_space(dset);
	  
	  // Select hyperslab
	  
	  hsize_t start[2] = {dims[0]-nrows, 0};
	  hsize_t count[2] = {nrows, ncols};
	  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start,
			      NULL, count, NULL);
	  
	  // Write buffer
	  
	  H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, file_space,
		   H5P_DEFAULT, hfield.getmemory());
	  
	  /*We can now release resources. */

	  H5Sclose(mem_space);
	  H5Dclose(dset);
	  H5Sclose(file_space);
	  H5Fclose(file);
	}
    }    
}

/*-----------------------------------------------------------------------*/

void Write_metadata_to_H5_file(const char* FILENAME, hfield_metadata* wdata, hsize_t DIM0) 
{
  
  hid_t   file, filetype, memtype, strtype, space, dset;
  herr_t  status;
  hsize_t  dims[1] = {DIM0};
  
  /* Open file using the default properties. */

  file = H5Fopen(FILENAME, H5F_ACC_RDWR, H5P_DEFAULT);

  /* Create variable-length string datatype. */

  strtype = H5Tcopy (H5T_C_S1);
  status = H5Tset_size (strtype, H5T_VARIABLE);
  
  /* Create the compound datatype for memory.*/
  
  memtype = H5Tcreate (H5T_COMPOUND, sizeof (hfield_metadata));
  status  = H5Tinsert (memtype, "Value",
		       HOFFSET (hfield_metadata, value), H5T_NATIVE_DOUBLE);
  status  = H5Tinsert (memtype, "Field", HOFFSET (hfield_metadata, field),
		       strtype);
  
  /* Create the compound datatype for the file. Because the standard
     types we are using for the file may have different sizes than
     the corresponding native types, we must manually calculate the
     offset of each member. */

  filetype = H5Tcreate (H5T_COMPOUND, 8 + sizeof (hvl_t) + 8 + 8);
  status = H5Tinsert (filetype, "Value", 0, H5T_NATIVE_DOUBLE);
  status = H5Tinsert (filetype, "Field", 8, strtype);
  
  /* Create dataspace.  Setting maximum size to NULL sets the maximum
     size to be the current size. */

  space = H5Screate_simple (1, dims, NULL);

    /* Create the dataset and write the compound data to it.*/

  dset = H5Dcreate (file, "metadata", filetype, space,
		    H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
  status = H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);

  /* Close and release resources. */

  status = H5Dclose (dset);
  status = H5Sclose (space);
  status = H5Tclose (filetype);
  status = H5Fclose (file);
}

/*---------------------------------------------------------------------*/

int Input_DoFs(int argc, char* argv[])
{  
  if (argc == 1 or argc > 2)
    {
      std::cout << "One command line argument is needed: The number of degrees" 
	"of freedom per dimension. Exiting." << std::endl;
      exit(-1);
    }
  else
    {
      int N = strtol(argv[1], nullptr, 0);
      if( (N>0) and (N%2 == 0) )       
	return N;
      else
	{
	  std::cout << "The number of degrees of freedom per" 
	    " dimension must be a positive and even number. Exiting."
		    << std::endl;
	  exit(-1);
	}
    }
}


/* WritePinnedSites------------- */
void WritePinnedSites(std::string pinset_filename,std::unordered_set<Site>& set)
{
  std::ofstream pinned;
  pinned.open(pinset_filename);
  for (auto it = set.begin(); it != set.end(); it++)
    pinned << (*it).getx() << " " << (*it).gety() << "\n";
  pinned.close();
}
/*------------------------*/

void ReadPinnedSites(std::string pinset_filename,
		      std::unordered_set<Site>& pinned_set)
{
  Site site;
  int x,y;
  std::ifstream pinset(pinset_filename); 
  if (pinset.is_open())
    {
      while (!pinset.eof())
	{
	  pinset >> x >> y;
	  site.set(x,y);
	  pinned_set.insert(site);
	}
    }
  else
    {
      std::cout << "File pinned_sites_rank.txt is not found." << std::endl;
      exit(EXIT_FAILURE);
    }
}

/*---------------------------------------------------------------------*/

void CopyFieldToArray(RectMesh& hfield, double* hx)
{
  int Nx = hfield.getcols();
  int Ny = hfield.getrows();
  for(int j=0; j<Ny; j++){
    for(int i=0; i<Nx; i++)
      hx[ i + (Nx+2)*j ] = hfield(i,j);
  }
}

/*---------------------------------------------------------------------*/

void WriteSpectrum(std::string hspec_filename, double* S1d, int spec_steps,
		   int qdiag_max, double dk)
{
  std::ofstream radSpecFile;
  radSpecFile.open(hspec_filename); 
  for (int i=0; i<qdiag_max/2; i++)
    radSpecFile << i*dk << "\t" << 4.0*S1d[i]/(double)spec_steps << "\n";
  radSpecFile.close();
}

/*---------------------------------------------------------------------*/

void PrintOut(int choose, int rank)
{
  int sim=rank+1;
  std::stringstream strm;
  switch (choose)
    {
    case 1:

      strm << "Simulation " + std::to_string(sim) + ": Input file is read.\n";
      std::cout << strm.str();
      break;
      
    case 2:

      strm << "Simulation " + std::to_string(sim) + ": Height field"
	" initialized.\n";
      std::cout << strm.str();
      break;
      
    case 3:
      
      strm << "Simulation " + std::to_string(sim) + ": MC-loop started.\n";
      std::cout << strm.str();
      break;

    case 4:
      
      strm << "Simulation " + std::to_string(sim) + ": WARNING:"
	" The spectrum is averaged over non equilibrium states.\n";
      std::cout << strm.str();
      break;
      
    case 5:
      
      strm << "Simulation " + std::to_string(sim) + ": Termination of"
	" MC-loop.\n";
      std::cout << strm.str();
      break;

    case 6: 

      strm << "Simulation " + std::to_string(sim) + ": Program finished"
	" successfully.\n";
      std::cout << strm.str();
      break;
      
    default:
      std::cout << "Here is a cookie as PrintOut" << std::endl;
      
    }
}
/*---------------------------------------------------------------------*/
