#ifndef _INFO_HH_
#define _INFO_HH_

const int MAX_STORY=100000;

//const Double_t TKE_projectile = 4767.4804;//total kinetic energy of the projectile A in inve
//const Double_t TKE_projectile = 1*4767.4804;//total kinetic energy of the projectile A in inverse kinematics case
const Double_t TKE_projectile = 630 * 25;//total kinetic energy of the projectile A in inverse kinematics case
const double ENERGY=400; //in direct-kinematics case - the energy of incident nucleon in MeV

const double NEUTRON_MASS= 939.565;	// MeV
const Double_t PROTON_MASS = 938.27208816;//MeV/c2
const Double_t MASS_11B_MEV_C2 = 10252.548;//MeV/c2
const Double_t MASS_12C_MEV_C2 = 11174.864;//MeV/c2
const Double_t MASS_11C_MEV_C2 = 10254.019;//MeV/c2

const double MASS_A = MASS_12C_MEV_C2; 
const double MASS_B = MASS_11C_MEV_C2;
//const double MASS_B = MASS_11B_MEV_C2;


const double EXE=0.0; //Excitation energy (MeV) of the outgoing heavy fragment
const double MOM_SIGMA= 106.5;	// Internal momentum spread (Gauss distribution)

//Constants
const double PI = TMath::Pi();
const double R2D = 57.29577951;
const double D2R = 0.017453292;
const double HBAR = 197.326960;
const double ALPHA = 0.007297352;
const double LN2 = 0.6931471805599452862;// ln2
const double CC = 299792458.;		// Speed of light
const double UNIT = 931.494013; 	// Atomic mass unit  (MeV/c^2)

const double rad2deg = 180./TMath::Pi();
const double deg2rad = TMath::Pi()/180.;

// Physics parameters for kinematic calculations
const double CUGNON_NORMALIZATION = 25.;   // Normalization for Cugnon parameterization
#endif

