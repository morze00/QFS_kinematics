#ifndef _HEADERS_HH_
#define _HEADERS_HH_

//A structure to store momentum and angles of two
//scattered particles in CM system
struct cm_values
{
	//Internal cluster
	double	e_clust;
	double	p_clust;
	double	theta_clust;
	//Scattered particle
	double	e_scat;
	double	p_scat;
	double	theta_scat;
	//indicates satisfactory kinematics (i.e. energy & momentum conservation)
	bool 	good;
	double T;
};

TVector3 DREHUNG(TVector3 v1,TVector3 v2);
double CINEMA(double,double,double);
cm_values CENMASS(double,double,double,double,bool,bool,bool);
std::pair<double, double> Lorentz(double,double,double,double);
double get_T(double,double);
double get_CosThetaCM_Cugnon(double,bool);
TH1F* Theoretical_Momentum(char* filename);

#endif //_HEADERS_HH_
