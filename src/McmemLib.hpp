/* This header file contains all the functions that are needed for the 
   simulation of the membrane. From the calculation of the derivatives and the
   local curvature of the membrane, to the updating of the energy and area and 
   sampling. */

#ifndef MCMEMLIB_HPP
#define MCMEMLIB_HPP

//#ifdef MAIN
//#else
extern int rank;
extern std::random_device rd;
extern std::mt19937 mt;
//#endif

/*---------------------*/
#include <unordered_set>
#include "Site.hpp"
#include "RectMesh.hpp"
#include <vector>
#include <string>
/*---------------------*/

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

/*------------------*/
typedef struct
{  
  const char*  field;
  double value;
  
}hfield_metadata;
/*------------------*/

/*---------------------------------------------------------------------*/
void OutputParams(const int maxiter,const int N,const int DoF,
		  const int nghost,const double rig,const double sig,
		  const double tau,const double epsilon,
		  const double min_change,const double max_change,
		  double pn_prcn,int blockRadius, double potStrength,
		  int sample_every, int rank, double Eactive);

std::unordered_set<Site> InitPinning(RectMesh& hfield, double pin_ratio,
				    Site neighbors[],int radius);

void InitSurface(RectMesh& hfield,std::unordered_set<Site>&,
		 const double min,const double max);

void GhostCopy(RectMesh& mesh);

void Der(const RectMesh& mesh, Site site, double alpha, double grad[]);

void Der2(const RectMesh& mesh, Site site, double alpha, double hess[]);

double SNodeArea(const RectMesh& field,Site site,double alpha);

double LocalArea(const RectMesh& hfield,Site neighbors[], double alpha);

double TotalArea(const RectMesh& field,double alpha);

double SNodeNormalZ(const RectMesh& field,Site site,double alpha);

RectMesh NormalZ(const RectMesh& field,double& alpha);

double SNodeCurvature(const RectMesh& field,Site site,double alpha);

double SNodeCurvatureEnergy(const RectMesh& field,Site site,
			    double alpha,const double rig);

double LocalCurvatureEnergy(const RectMesh& field,Site neighbors[],
			    double alpha,const double rig);

double CurvatureEnergyTotal(const RectMesh& field,
			    double alpha,const double rig);

double SNodeCorrectionEnergy(const RectMesh& field,Site site,double alpha);

double LocalCorrectionEnergy(const RectMesh& field,
			     const Site neighbors[],double alpha);

double LocalPinEnergy(const RectMesh& hfield,Site& site,
		      const double& pot_strength,const double& h0);

double CorrectionEnergyTotal(const RectMesh& field,double alpha);

double PinningEnergyTotal(const RectMesh& hfield,
			  std::unordered_set<Site>& pinned_sites,
			  double& pot_strength,double h0);

double LocalEnergy(const RectMesh& hfield,
		   Site neighbors_area[],
		   Site neighbors_corr[],
		   Site neighbors_ener[],
		   double alpha,const double rig,
		   const double sig,const double tau,
		   const double& pot_strength,const double& h0,bool pin);

void CalculateTotal(const RectMesh& hfield,const double& rig,
		    const double& sig,const double& tau,double& tot_energy,
		    double& tau_energy,double& crv_energy,double& sig_energy,
		    double& cor_energy,double& pin_energy,double& tot_area,
		    double& prj_area,double& alpha,
		    std::unordered_set<Site>& pinned_sites,
		    const double& pot_strength,const double& h0);

bool WhereIs(Site site,int cols,int rows,int nghost);

bool Ispinned(Site& site,std::unordered_set<Site>& pinned_sites);

void GetNeighbors(const RectMesh& field,Site site,Site neighbors_area[],
		  Site neighbors_corr[],Site neighbors_ener[]);

void PrintNeighbors(Site neighbors[],int len);

bool Metropolis(double& dElocal);

bool UpdateState(RectMesh& hfield,Site site,
		     bool where,double& tot_area,double& tot_energy,
		     double& dAlocal,double& dElocal,
		     int& height_moves,double& perturb);

void WriteStats(const int maxiter, int height_moves, int lattice_moves,
	    int lattice_changes,int spec_steps,int rank);

bool UpdateLattice(const RectMesh& hfield,const double& min_change,
		   const double& max_change,const double& rig,
		   const double& sig, const double& tau,double& prj_area,
		   double& tot_area,double& tot_energy,double& tau_energy,
		   double& crv_energy,double& sig_energy,double& cor_energy,
		   double& pin_energy,double& alpha,
		   int& lattice_moves,int& lattice_changes,
		   std::unordered_set<Site>& pinned_sites,
		   const double& pot_strength,const double& h0);

void Sample(std::string output_filename, int& iter,
	    int& total_moves, double& tot_energy, double& crv_energy,
	    double& cor_energy, double& pin_energy, double& tot_area,
	    double& prj_area, double& alpha, int& DoFs);

void ReadInput(std::string filename,int& sim, int& wSnap, int& acc_samples,
	       int& blockRadius,double& maxiter,double& sig,double& tau,
	       double& epsilon,double& min_change,double& max_change,
	       double& pin_ratio,double& potStrength, double& Ea);

void AddShift(double& dE);

void Spectrum(RectMesh& Input_Field,RectMesh& Output);

void WriteToExtendibleH5(const char* FILENAME, RectMesh& hfield);

void WriteMetadataToH5File(const char* FILENAME,
			       hfield_metadata* wdata, hsize_t DIM0); 

int InputDoFs(int argc, char* argv[]);

void ReadPinnedSites(std::string pinset_filename,
		     std::unordered_set<Site>& pinned_set);

void WritePinnedSites(std::string pinset_filename,std::unordered_set<Site>& set);

void TileOfNeighbors(Site neighbors[], RectMesh& field, Site site, int radius);

void CopyFieldToArray(RectMesh& hfield, double* hx);

void PrintOut(int choose, int rank);

void ReadAlpha(std::string filename, double& alpha);

void WriteAlpha(int rank, double alpha);

/*---------------------------------------------------------------------*/

#endif
