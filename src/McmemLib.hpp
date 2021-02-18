/* This header file contains all the functions that are needed for the 
   simulation of the membrane. From the calculation of the derivatives and the
   local curvature of the membrane, to the updating of the energy and area and 
   sampling. */

#ifndef MCMEMLIB_HPP
#define MCMEMLIB_HPP

#include <unordered_set>
#include "Site.hpp"
#include "RectMesh.hpp"
#include <vector>
#include <string>


void OutputParams(const int maxiter,const int N,const int DoF,
		  const int nghost,const double rig,const double sig,
		  const double tau,const double epsilon,
		  const double min_change,const double max_change,
		  double alpha,double pn_prcn,int sample_every,int rank);
std::unordered_set<Site> InitPinning(int N,double pn_prcn);
void InitSurface(RectMesh& hfield,const double min,const double max,
		 std::unordered_set<Site>& pinned_sites);
void GhostCopy(RectMesh& mesh);
void Der(const RectMesh& mesh, Site site, double alpha, double grad[]);
void Der2(const RectMesh& mesh, Site site, double alpha, double hess[]);
double SNodeArea(const RectMesh& field,Site site,double alpha);
double LocalArea(const RectMesh& hfield,Site neighbors[], double alpha);
double TotalArea(const RectMesh& field,double alpha);
double SNodeNormalZ(const RectMesh& field,Site site,double alpha);
RectMesh NormalZ(const RectMesh& field,double& alpha);
double SNodeCurvature(const RectMesh& field,Site site,double alpha);
//RectMesh TotalCurvature(RectMesh& field,double& alpha);
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
// double TotalEnergy(RectMesh& hfield,double& total_area,
// 		   double& prj_area,double alpha,
// 		   double rig, double sig, double tau);
void CalculateTotal(const RectMesh& hfield,const int& DoF,const double& rig,
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
void AcceptOrDecline(RectMesh& hfield,Site site,bool accept,
		     bool where,double& tot_area,double& tot_energy,
		     double& dAlocal,double& dElocal,
		     int& accepted_moves,double& perturb);
void PrintAcceptance(const int maxiter, int accepted_moves,
		     int lattice_moves, int lattice_changes, int rank);
void ChangeLattice(const RectMesh& hfield,const double& min_change,
		   const double& max_change,const int& DoF,const double& rig,
		   const double& sig, const double& tau,double& prj_area,
		   double& tot_area,double& tot_energy,double& tau_energy,
		   double& crv_energy,double& sig_energy,double& cor_energy,
		   double& pin_energy,double& alpha,int& move_counter,
		   int& lattice_moves,int& lattice_changes,
		   std::unordered_set<Site>& pinned_sites,
		   const double& pot_strength,const double& h0);
void Sample(int& iter,int& sample_every,
	    int& lattice_changes,std::string filename,
	    double& tot_energy,double& tau_energy,
	    double& crv_energy,double& sig_energy,
	    double& cor_energy,double& pin_energy,double& tot_area,
	    double& prj_area,double& alpha,const int& DoF);
void ReadInput(std::string filename,double& maxiter,double& sig,double& tau,
	       double& epsilon,double& min_change,double& max_change,
	       double& pin_ratio,int& acc_samples);

// void ReadTxt(const char filename, std::vector<double> &data);
// void ReadTxtInt(const char filename[], std::vector<int> &data);

#endif
