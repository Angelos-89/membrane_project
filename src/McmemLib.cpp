#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>
#include <string>
#include "Vec3dLib.hpp"
#include "McmemLib.hpp"

std::random_device rd;
std::mt19937 mt(rd());

/*-------------------InitSurface---------------------------*/

void InitSurface(RectMesh& hfield,double min,double max)
{
  std::uniform_real_distribution<double> UnifProb(min,max);
  for (int j=0; j<hfield.getrows(); j++)
  {
    for (int i=0; i<hfield.getcols(); i++)
      hfield(i,j) = UnifProb(mt);
  }
  GhostCopy(hfield);
}

/*----------------GhostCopy----------------*/

/*This function implements periodic boundary
  conditions by copying the first and last n 
  rows and columns of a 2D mesh at the ghost 
  zones of the same  mesh. n denotes the number 
  of ghost points per side.

  Example:

  A 2D mesh defined as: RectMesh field(3,4,2);
  with 3 columns, 4 rows and 2 ghost points 
  is shown below before the GhostCopy.

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

/*------------------------Der----------------------------------*/

/* This function calculates the gradient of a scalar field F(x,y) 
   at a given point, where i corresponds to the x-coordinate and 
   j to the y-coordinate. The pointer grad represents an 1D array 
   where grad[0] stores the partial derivative with respect to x 
   and grad[1] the partial derivative with respect to y. 
   Central finite difference schemes are used to compute the 
   derivatives. alpha is the lattice spacing.                  */

void Der(RectMesh& mesh, Site site, double alpha, double grad[])
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

/*-------------------------Der2-------------------------*/

/* This function calculates the second order derivatives of 
   a scalar field F(x,y) at a given point where j corresponds 
   to the x-coordinate and i to the y-coordinate. The pointer 
   hess represents the Hessian matrix as a 1D array where 
   hess[0] stores Fxx, hess[1] stores Fxy, hess[2] stores Fyx 
   and hess[3] stores Fyy. Central finite difference schemes 
   are used to compute the derivatives. alpha is the lattice spacing.
*/

void Der2(RectMesh& mesh, Site site, double alpha, double hess[])
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

/*---------------------------SNodeArea----------------------------*/

/* This function calculates the local area corresponding to a node
   represented by (i,j). We use equation (7) from the paper "Monte 
   Carlo study of the frame fluctuation and internal tensions of 
   fluctuating membranes with fixed area" by Shiba et al. (2016). 
   That is, we  calculate the area of the four triangles that are 
   formed by the value of the field at (i,j) and its nearest 
   neighbors and then devide by 2. Note however, that we do not 
   multiply by 0.5 to find the area of the triangle, but multiply by 
   0.25 at the end instead.                                       */

double SNodeArea(RectMesh& field,Site site,double alpha)
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

/*-------------------------LocalArea----------------------------*/

double LocalArea(RectMesh& hfield, Site neighbors[], double alpha)
{
  double area = 0;
  for (int n=0; n<5; n++)
    area += SNodeArea(hfield,neighbors[n],alpha);
  return area;  
}

/*---------------------------TotalArea------------------------*/

/* This function calculates the total area of the surface that 
   is represented by the values of a scalar field. The scalar 
   field is represented by the class RectMesh and it is assumed 
   to be periodic in both x and y. The step size in each 
   direction is alpha.                                        */

double TotalArea(RectMesh& field,double alpha)
{
  if (field.getnghost() == 0)
    {
    std::cout << "TotalArea: The number of ghost points per side of the"
      "RectMesh object must be at least 1! Exiting. \n"
	      << std::endl;
    exit(EXIT_FAILURE);
  }
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

/*----------------------SNodeNormalZ----------------------*/

/* This function calculates the z coordinate of the unit 
   normal vector at one point of the surface that is 
   represented by the RectMesh object hfield. alpha is the 
   lattice spacing and site represents the point of the mesh 
   at which the unit normal coordinate is calculated.
*/

double SNodeNormalZ(RectMesh& hfield,Site site,double alpha)
{
  double grad[2]={};
  Der(hfield,site,alpha,grad);
  double h_x = grad[0];
  double h_y = grad[1];
  return (1. / sqrt(1. + pow(h_x,2) + pow(h_y,2)));
};

/*-------------------NormalZ-----------------*/

/* This function calculates the z coordinates 
   of the unit normal vectors at each point of 
   a surface represented by hfield. alpha is the 
   lattice spacing. NormalZ calls LocalNormalZ 
   for every point of the surface.           */

RectMesh NormalZ(RectMesh& hfield, double alpha)
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

/*-----------------------SNodeCurvature------------------------------*/

/* This function calculates the local mean curvature, i.e., the mean 
   curvature H at a given point (i,j) in the Monge representaion. Here, 
   hfield is the surface represented by a RectMesh object. Note that the 
   GhostCopy function must be called before calling this function for 
   boundary points.                                                  */

double SNodeCurvature(RectMesh& h,Site site,double alpha)
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
  double denominator = 2.*pow(1.0+h_x2+h_y2, 3./2.);
  double curvature = numerator/denominator;
  return curvature;
}

/*------------------------TotalCurvature------------------------*/

/* This function returns an Eigen Matrix with the mean curvature 
   at each point of a surface that is represented by a RectMesh 
   object.                                                      */

RectMesh TotalCurvature(RectMesh& field,double alpha)
{
  Site site;
  RectMesh H(field.getcols(),field.getrows(),0);
  for (int j=0; j<field.getrows(); j++)
    {
      for (int i=0; i<field.getcols(); i++)
	{
	  site.set(i,j);
	  H(i,j) = SNodeCurvature(field,site,alpha);
	}
    }
  return H;
}

/*----------------------SNodeCurvatureEnergy------------------------*/

/* This function calculates the curvature energy of a single site 
   due to its mean curvature H. This energy is the integrand of 
   equation (3) found in the paper "Monte Carlo study of the frame 
   fluctuation and internal tensions of fluctuating membranes with 
   fixed area" by Shiba et al. (2016). Here, kk is the bending 
   rigidity and hfield the RectMesh instance representing the membrane 
   in the Monge gauge.                                              */

double SNodeCurvatureEnergy(RectMesh& hfield,Site site,double alpha,double rig)
{
  double dA = SNodeArea(hfield,site,alpha);
  double H = SNodeCurvature(hfield,site,alpha);
  return 0.5*rig*H*H*dA; 
}

/*----------------------LocalCurvatureEnergy------------------------*/

double LocalCurvatureEnergy(RectMesh& field,Site neighbors[],
			    double alpha,double rig)
{
  int len = pow((2*field.getnghost()+1),2);
  double energy = 0;
  for (int n=0; n<len; n++)  
    energy += SNodeCurvatureEnergy(field,neighbors[n],alpha,rig);  
  return energy;
}

/*-----------------------CurvatureEnergyTotal---------------------*/

/* This function calculates the total energy of the membrane due
   to its curvature by summing up the contributions from each 
   site. Again, hfield is the RectMesh object representing the
   membrane height field and kk is the bending rigidity.          */

double CurvatureEnergyTotal(RectMesh& hfield,double alpha,double rig)
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

/*---------------------SNodeCorrectionEnergy-----------------------*/

/* This function calculates the correction in energy at a single 
   node,due to the Monge gauge representation of the membrane.     */

double SNodeCorrectionEnergy(RectMesh& hfield,Site site,double alpha)
{
  return -log(SNodeNormalZ(hfield,site,alpha));
}

/*---------------------LocalCorrectionEnergy-----------------------*/

double LocalCorrectionEnergy(RectMesh& field,Site neighbors[],double alpha)
{
  int len = 4*field.getnghost()+1;
  double energy = 0;
  for (int n=0; n<len; n++)  
    energy += SNodeCorrectionEnergy(field,neighbors[n],alpha);  
  return energy;
}

/*-------------------------CorrectionEnergyTotal------------------------------*/

/* This function calculates the total correction in energy due to the Monge
   gauge representation of the membrane. For details see  "Monte Carlo study 
   of the frame fluctuation and internal tensions of fluctuating membranes with
   fixed area" by Shiba et al. (2016). The correction term for the energy is
   given by H_{coor} = -k_b*T*sum_{i,j}{ln(NormalZ(i,j))}, where k_b is
   Boltzmann's constant, T the temperature and NormalZ the z coordinates of the 
   unit normal vectors along the surface of the membrane. Note that we work in 
   reduced units i.e., k_b*T=1.                                               */

double CorrectionEnergyTotal(RectMesh& hfield, double alpha)
{
  RectMesh temp = NormalZ(hfield,alpha);
  temp.ln();
  return -temp.sum();
}

/*---------------------------LocalEnergy----------------------------------*/

double LocalEnergy(RectMesh& hfield,
		   Site neighbors_area[],
		   Site neighbors_corr[],
		   Site neighbors_ener[],
		   double alpha,double rig,
		   double sig, double tau)
{
  double tau_ener = -tau*alpha*alpha; 
  double crv_ener = LocalCurvatureEnergy(hfield,neighbors_ener,alpha,rig);
  double sig_ener = sig*LocalArea(hfield,neighbors_area,alpha);
  double cor_ener = LocalCorrectionEnergy(hfield,neighbors_corr,alpha);
  return tau_ener + crv_ener + sig_ener + cor_ener;
}

/*-----------------------TotalEnergy---------------------------*/

double TotalEnergy(RectMesh& hfield,double& tot_area,
		   double& prj_area,double alpha,
		   double rig, double sig,double tau)
{
  double tau_ener = -tau*prj_area;
  double crv_ener = CurvatureEnergyTotal(hfield,alpha,rig);
  double sig_ener = sig*tot_area;
  double cor_ener = CorrectionEnergyTotal(hfield,alpha);
  return tau_ener + crv_ener + sig_ener + cor_ener;
}

/*--------------------------------CalculateTotal----------------------------*/

void CalculateTotal(RectMesh& hfield, double& tot_energy, double& tau_energy,
		    double& crv_energy, double& sig_energy,
		    double& cor_energy, double& tot_area,
		    double& prj_area, int& DoF, double& alpha,
		    double& rig,double& sig, double& tau)
{
  prj_area = (double)DoF*alpha*alpha;
  tot_area = TotalArea(hfield,alpha);
  tau_energy = -tau*prj_area;
  crv_energy = CurvatureEnergyTotal(hfield,alpha,rig);
  sig_energy = sig*tot_area;
  cor_energy = CorrectionEnergyTotal(hfield,alpha);
  tot_energy = tau_energy + crv_energy + sig_energy + cor_energy;
}

/*---------------------WhereIs----------------------*/

/* This function takes as arguments a site(i,j), the 
   columns, rows and ghost points of a RectMesh object 
   and returns an integer value that depends on 
   whether it belongs to the boundary or the bulk.  */

bool WhereIs(Site site, int cols, int rows, int nghost)
{
  int i = site.getx();
  int j = site.gety();
  if (i>=cols || i<0 || j>=rows || j<0)
  {
    std::cout << "WhereIs: It must hold that 0<i<cols and 0<j<rows."
      "Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (i<nghost || i>=cols-nghost || j<nghost || j>=rows-nghost) return 1;
  else return 0;
}

/*----------------GetNeighbors-------------*/

void GetNeighbors(RectMesh& field,Site site,
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
  
  
  /* fill the list of neighbors needed for 
     the local area calculation */
  
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
     the local correction energy */ 

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
     the local curvature energy */ 

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

/*-------------PrintNeighbors--------------*/
void PrintNeighbors(Site neighbors[],int len)
{
  for(int i=0; i<len; i++)
    neighbors[i].print();
}

/*---------------Metropolis----------------*/

bool Metropolis(double dE)
{
  if(dE < 0) return 1;
  else
  {
    std::uniform_real_distribution<double> UnifProb(0,1);
    double r = UnifProb(mt);
    if (r < exp(-dE)) return 1;
    else return 0;
  }
}

/*--------------------------AcceptOrdecline---------------------------*/

void AcceptOrDecline(RectMesh& hfield,Site site,bool accept,bool where,
		     double& tot_area,double& tot_energy,double dAlocal,
		     double dElocal,int& accepted_moves,double perturb)
{  
  if (accept == 1)
    {
      accepted_moves ++;
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

/*---------------------------PrintAcceptance-------------------------*/

void PrintAcceptance(int maxiter, int accepted_moves,
		     int lattice_moves, int lattice_changes)
{
  double accept_ratio = (double) accepted_moves/maxiter;
  double lattice_moves_ratio = (double) lattice_changes/lattice_moves;
  std::cout << std::endl; 
  std::cout << "Height move acceptance ratio: " << accept_ratio << "\n"
	    << "Lattice spacing move acceptance ratio: "
	    << lattice_moves_ratio << std::endl;
  std::cout << std::endl;
}

/*-----------------------Sample----------------------------------------*/

void Sample(int& iter,std::string filename,double& tot_energy,
	    double& tau_energy,double& crv_energy,
	    double& sig_energy,double& cor_energy,
	    double& tot_area,double& prj_area,
	    double& alpha,int& DoF)
{
  if(iter%1000==0)
    {
      double DOF = (double) DoF;
      std::ofstream file;
      file.open(filename, std::ios::app);
      file << std::setprecision(12) << tot_area/DOF     << "\t"
	   << std::setprecision(12) << prj_area/DOF     << "\t"
	   << std::setprecision(12) << alpha            << "\t"
	   << std::setprecision(12) << tau_energy/DOF   << "\t"
	   << std::setprecision(12) << crv_energy/DOF   << "\t"
	   << std::setprecision(12) << sig_energy/DOF   << "\t"
	   << std::setprecision(12) << cor_energy/DOF   << "\t"
	   << std::setprecision(12) << tot_energy/DOF
	   << std::endl;
      file.close();
    }
}

/*-------------------------ChangeLattice---------------------------*/

void ChangeLattice(RectMesh& hfield,int& move_counter,double& alpha,
		       double& prj_area,double& tot_area,int& DoF,
		       double& tot_energy,double& tau_energy,
		       double& crv_energy,double& sig_energy,
		       double& cor_energy,double& rig,
		       double& sig,double& tau,int& lattice_moves,
		   int& lattice_changes)
{
  if(move_counter !=0 && move_counter % 5 == 0)
    {
      lattice_moves ++;
      std::uniform_real_distribution<double>  RandDouble(0.98,1.02); 
      double old_alpha = alpha;
      double epsilon = RandDouble(mt);
      alpha *= epsilon;
      double old_prj_area = prj_area;
      double old_tot_area = tot_area;
      double old_energy = tot_energy;
      CalculateTotal(hfield,tot_energy,tau_energy,crv_energy,
		     sig_energy,cor_energy,tot_area,prj_area,DoF,
		     alpha,rig,sig,tau);
      double dE = tot_energy - old_energy; 
      bool accept = Metropolis(dE);
      if (accept == 1)
	lattice_changes ++;
      else
	{
	  tot_energy = old_energy;
	  prj_area = old_prj_area;
	  tot_area = old_tot_area;
	  alpha = old_alpha;
	}
      move_counter = 0;
    }
}





// void ReadTxt(const char filename[], std::vector<double> &data)
// {
//   double number;
//   std::ifstream myfile(filename);
//   if (myfile.is_open())
//     {
//       while (myfile >> number)
// 	data.push_back(number);
//       myfile.close();
//     }

//   else
//     {
//       std::cout << "Unable to open file. Exiting." << std::endl; 
//       exit(EXIT_FAILURE);
//     }
// }

// void ReadTxtInt(const char filename[], std::vector<int> &data)
// {
//   int number;
//   std::ifstream myfile(filename);
//   if (myfile.is_open())
//     {
//       while (myfile >> number)
// 	data.push_back(number);
//       myfile.close();
//     }

//   else
//     {
//       std::cout << "Unable to open file. Exiting." << std::endl; 
//       exit(EXIT_FAILURE);
//     }
// }
