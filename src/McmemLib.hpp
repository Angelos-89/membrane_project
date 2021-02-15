/* This header file contains all the functions that are needed for the 
   simulation of the membrane. From the calculation of the derivatives and the
   local curvature of the membrane, to the updating of the energy and area and 
   sampling. */

#ifndef MCMEMLIB_HPP
#define MCMEMLIB_HPP

#include "Site.hpp"
#include "RectMesh.hpp"
#include <vector>
#include <string>

void InitSurface(RectMesh& hfield,const double min,const double max);
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
double CorrectionEnergyTotal(const RectMesh& field, double alpha);
double LocalEnergy(const RectMesh& hfield,
		   Site neighbors_area[],
		   Site neighbors_corr[],
		   Site neighbors_ener[],
		   double alpha,const double rig,
		   const double sig,const double tau);
// double TotalEnergy(RectMesh& hfield,double& total_area,
// 		   double& prj_area,double alpha,
// 		   double rig, double sig, double tau);
void CalculateTotal(const RectMesh& field,const int& DoF,const double& rig,
		    const double& sig,const double& tau, double& tot_energy,
		    double& tau_energy,double& crv_energy,double& sig_energy,
		    double& cor_energy,double& tot_area,
		    double& prj_area,double& alpha);
bool WhereIs(Site site,int cols,int rows,int nghost);
void GetNeighbors(const RectMesh& field,Site site,Site neighbors_area[],
		  Site neighbors_corr[],Site neighbors_ener[]);
void PrintNeighbors(Site neighbors[],int len);
bool Metropolis(double& dElocal);
void AcceptOrDecline(RectMesh& hfield,Site site,bool accept,
		     bool where,double& tot_area,double& tot_energy,
		     double& dAlocal,double& dElocal,
		     int& accepted_moves,double& perturb);
void PrintAcceptance(const int maxiter, int accepted_moves,
		     int lattice_moves, int lattice_changes);

void ChangeLattice(const RectMesh& hfield,const double& min_change,
		   const double& max_change,const int& DoF,const double& rig,
		   const double& sig, const double& tau,double& prj_area,
		   double& tot_area,double& tot_energy,double& tau_energy,
		   double& crv_energy,double& sig_energy,double& cor_energy,
		   double& alpha,int& move_counter,int& lattice_moves,
		   int& lattice_changes);
void Sample(int& iter, std::string filename,double& tot_energy,
	    double& tau_energy,double& crv_energy,
	    double& sig_energy,double& cor_energy,
	    double& tot_area,double& prj_area,
	    double& alpha,const int& DoF);
void ReadTensions(std::string filename,double& sig,double& tau);
void ReadInputs(std::string filename,std::vector<double>& input);
void AddShift(double& dE);

// void ReadTxt(const char filename, std::vector<double> &data);
// void ReadTxtInt(const char filename[], std::vector<int> &data);

#endif
