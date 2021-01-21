/* This header file contains all the functions that are needed for the 
   simulation of the membrane. From the calculation of the derivatives and the
   local curvature of the membrane, to the updating of the energy and area and 
   sampling. */

#ifndef MCMEMLIB_HPP
#define MCMEMLIB_HPP

#include "Site.hpp"
#include "RectMesh.hpp"
#include <vector>

void InitSurface(RectMesh& hfield,double min,double max);
void GhostCopy(RectMesh& mesh);
void Der(RectMesh& mesh, Site site, double alpha, double grad[]);
void Der2(RectMesh& mesh, Site site, double alpha, double hess[]);
double SNodeArea(RectMesh& field,Site site,double alpha);
double LocalArea(RectMesh& hfield, Site neighbors[], double alpha);
double TotalArea(RectMesh& field,double alpha);
double SNodeNormalZ(RectMesh& field,Site site,double alpha);
double SNodeCurvature(RectMesh& field,Site site,double alpha);
double SNodeCurvatureEnergy(RectMesh& field,Site site,double alpha,double rig);
double LocalCurvatureEnergy(RectMesh& field,Site neighbors[],
			    double alpha,double rig);
double CurvatureEnergyTotal(RectMesh& field, double alpha, double rig);
double SNodeCorrectionEnergy(RectMesh& field,Site site,double alpha);
double LocalCorrectionEnergy(RectMesh& field,Site neighbors[],double alpha);
double CorrectionEnergyTotal(RectMesh& field, double alpha);
double LocalEnergy(RectMesh& hfield,
		   Site neighbors_area[],
		   Site neighbors_corr[],
		   Site neighbors_ener[],
		   double alpha,double rig,
		   double sig,double tau);
double TotalEnergy(RectMesh& hfield,double& total_area,
		   double& prj_area,double alpha,
		   double rig, double sig, double tau);
RectMesh TotalCurvature(RectMesh& field,double alpha);
RectMesh NormalZ(RectMesh& field,double alpha);
void CalculateTotal(RectMesh& field, double& tot_energy,double& tau_energy,
		    double& crv_energy,double& sig_energy,
		    double& cor_energy,double& tot_area,
		    double& prj_area,int& DoF,double& alpha,
		    double& rig, double& sig,double& tau);
void GetNeighbors(RectMesh& field,Site site, Site neighbors_area[],
		  Site neighbors_corr[],Site neighbors_ener[]);
void PrintNeighbors(Site neighbors[],int len);
void AcceptOrDecline(RectMesh& hfield,Site site,bool accept,
		     bool where,double& tot_area,double& tot_energy,
		     double dAlocal,double dElocal,
		     int& accepted_moves,double perturb);
void PrintAcceptance(int maxiter, int accepted_moves,
		     int lattice_moves, int lattice_changes);
void Sample(int& iter, std::string filename,double& tot_energy,
	    double& tau_energy,double& crv_energy,
	    double& sig_energy,double& cor_energy,
	    double& tot_area,double& prj_area,
	    double& alpha,int& DoF);
void ChangeLattice(RectMesh& hfield,double& min_change,double& max_change,
		   int& move_counter,double& alpha,
		   double& prj_area,double& tot_area,int& DoF,
		   double& tot_energy,double& tau_energy,
		   double& crv_energy,double& sig_energy,
		   double& cor_energy,double& rig,
		   double& sig,double& tau,int& lattice_moves,
		   int& lattice_changes);
bool WhereIs(Site site,int cols,int rows,int nghost);
bool Metropolis(double dE);


// void ReadTxt(const char filename, std::vector<double> &data);
// void ReadTxtInt(const char filename[], std::vector<int> &data);

#endif
