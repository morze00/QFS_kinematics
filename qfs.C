/* Original code was written in FORTRAN by Leonid Chulkov l.chulkov@gsi.de
 *
 * QFS kinematics generator for the reactions in direct and inverse kinematics:
 * 
 * A(p,2p)B, A(p,pn)B, A(n,np)B 
 * 
 * use --direct-kinematics, otherwise default inverse kinematics will be used.
 * 
 * In inverse kinematics the program simulates inverse kinematics of quasi-elastic scattering
 * of a target proton "i" (being at rest) on a cluster/nucleon "a" bound in the projectile nucleus A
 * 
 * A --> i  ==> B + (i -> a) =>  B + 1 + 2
 * 
 * The kinematics of the residual heavy spectator "B", the knocked-out outgoing cluster (2)
 * and the scattered target nucleon (1) is stored in the output tree.
 *
 * Modify info.hh to change beam parameters, masses, internal momentum spread, etc.
 * 
 * Pass --direct-kinematics flag to simulate in direct kinematics the scattering of
 * of a beam-like proton "i" on a cluster/nucleon "a" bound in target nucleus A which is at rest.
 *
 * See output tree variables for the kinematics of the particles
 * 1, 2 and B in the lab and CM frame and other kinematical variables
 * 
 * Valerii Panin October 2026 (v.panin@gsi.de)
 */

#include "libs.hh"
#include "headers.hh"
#include "info.hh"

using namespace std;

//executable function called from main()
void run(Char_t * root_filename,
        int nevents,
        Bool_t is_ppn,
        Bool_t is_isotropic,
        Bool_t is_direct,
        Bool_t is_nnp
        )
{
    bool is_p2p = false;

    // Total kinetic energy (MeV) of the projectile
    double Tkin = TKE_projectile;
    if(is_direct) Tkin = ENERGY;

    const double MA   = MASS_A;  	
    const double MB   = MASS_B + EXE; 

    double Mi = PROTON_MASS; // default
    double Ma = PROTON_MASS;//default

    cout << "\n-- Max events = " << nevents;
    
    if(is_direct){
        cout << "\n-- Direct kinematics";
    }
   else{
        cout << "\n-- Inverse kinematics";
    
    }
    
    if(is_ppn){
        Ma   = NEUTRON_MASS;  
        cout << "\n-- Generating (p,pn) reactions" << endl;
        cout << "-- Neutron Mass: " << Ma << " MeV";
    }
    else if (is_nnp)
    {
        if (is_ppn)
        {
            cout << "\n\nChoose either (p,pn) or (n,np) reactions, not both! Abortin...\n\n";
            return;
        }

        cout << "\n-- Generating (n,np) events" << endl;
        cout << "-- Proton Mass: " << Ma << " MeV";
        Mi = NEUTRON_MASS;
    }
    else
    {
        cout << "\n-- Generating (p,2p) events" << endl;
        cout << "-- Proton Mass: " << Ma << " MeV";
        is_p2p = true;
    }

    if (is_isotropic)
    {
        cout << "\n-- Using isotropic NN cross sections";
    }

    else if (!is_isotropic)
    {
        cout << "\n-- Using parametrized NN cross sections (Cugnon)";
    }

    //------- Output root tree ----------
    TFile file(root_filename, "RECREATE");
    TTree *tree = new TTree("tree", "Tree with simulated QFS kinematics");

#include "tree.hh"

    double PA, EA, bA, gA, Pi, bi, gi, EaL, EBL;
    double S_first, sigma;

    //--------------------- Beam parameters -----------------------------

    if (!is_direct) // inverse kinematics
    {
        PA = sqrt(Tkin * (Tkin + 2 * MA));         // Total 3-momentum of the beam
        EA = sqrt(MA * MA + PA * PA);              // Total energy of the beam
        bA = -PA / EA;                             // Beta of the beam
        gA = 1 / sqrt(1 - bA * bA);                // Gamma of the beam
        S_first = (EA + Mi) * (EA + Mi) - PA * PA; // Invariant mass (Mandelstam S-variable)
        sigma = MOM_SIGMA;                         // Internal momentum spread of the cluster "a" inside "A"
        Ei = Mi;

        cout << "\n\n****** Beam parameters ********";
        cout << "\n-- MA:\t" << MA << " MeV";
        cout << "\n-- MB:\t" << MB << " MeV";
        cout << "\n-- Mi:\t" << Mi << " MeV";
        cout << "\n-- Ma:\t" << Ma << " MeV";
        cout << "\n-- Total momentum of the incident nucleus:\t" << PA << " MeV/c";
        cout << "\n-- Total energy of the incident nucleus:\t" << EA << " MeV";
        cout << "\n-- Beta (beam):\t" << (-bA);
        cout << "\n-- Gamma (beam):\t" << gA;
        cout << "\n-- Initial Mandelstam S:\t" << S_first << " MeV^2\n\n";
        cout << "\n-- Generating " << nevents << " events........\n\n";
    }

    else // direct kinematics
    {
        Pi = sqrt(Tkin * (Tkin + 2 * Mi)); // Total 3-momentum of the projectile
        Ei = sqrt(Mi * Mi + Pi * Pi);      // Total energy of the projectile
        bi = -Pi / Ei;                     // Beta of the projectile
        gi = 1 / sqrt(1 - bi * bi);        // Gamma of the projectile
        // S_first = (Ei+MA)*(Ei+MA) - Pi*Pi;	// Invariant mass (Mandelstam S-variable)
        S_first = (Ei + MA) * (Ei + MA) - Pi * Pi; // Invariant mass (Mandelstam S-variable)
        sigma = MOM_SIGMA;                         // Internal momentum spread of the cluster "a" inside "A"

        cout << "\n\n****** Beam parameters (direct kinematics) ********";
        cout << "\n-- MA:\t" << MA << " MeV";
        cout << "\n-- MB:\t" << MB << " MeV";
        cout << "\n-- Mi:\t" << Mi << " MeV";
        cout << "\n-- Ma:\t" << Ma << " MeV";
        cout << "\n-- Total momentum:\t" << Pi << " MeV/c";
        cout << "\n-- Total energy:\t" << Ei << " MeV";
        cout << "\n-- Beta (beam): \t" << (-bi);
        cout << "\n-- Gamma (beam):\t" << gi;
        cout << "\n-- InitialMandelstam S:\t" << S_first << " MeV^2\n\n";
        cout << "\n-- Generating " << nevents << " events........\n\n";
    }
    TRandom3 r1;
    r1.SetSeed(0);

    TLorentzVector LVa; // Lorentz vector of the cluster
    TLorentzVector LVi; // Lorentz vector of the incident particle

    if (is_direct)
        LVi.SetPxPyPzE(0.0, 0.0, Pi, Ei);
    else
        LVi.SetPxPyPzE(0.0, 0.0, 0.0, Mi);

    TVector3 Pa(1e-9, 0.0, 0.0);
    TVector3 P1cm(1e-9, 0.0, 0.0);
    TVector3 P2cm(1e-9, 0.0, 0.0);
    TVector3 PB;

    int events = 0;          // event counter
    while (events < nevents) // eventloop
    {
        // Momentum distribution of the cluster "a" inside the projectile "A"
        // modify sigma in info.hh or implement proper momentum distributions
        Pa.SetX(r1.Gaus(0, sigma));
        Pa.SetY(r1.Gaus(0, sigma));
        Pa.SetZ(r1.Gaus(0, sigma));

        //------------ Internal momentum of the residual-----------------
        PB.SetX(-Pa.X());
        PB.SetY(-Pa.Y());
        PB.SetZ(-Pa.Z());

        // Filling tree variables for the fragmnent
        PBx = PB.X();
        PBy = PB.Y();
        PBz_rf = PB.Z(); // Longitudinal momentum in the restframe

        // From the energy conservation in the virtual dissociation A->B+a
        double rrtt = MA * MA + MB * MB - 2 * MA * sqrt(MB * MB + Pa.Mag2());
        if (rrtt <= 0)
        {
            cout << "\nERROR off-shell mass!! rrtt=" << rrtt;                             // non-zero and real off-shell mass
            cout << "\nP:" << PBx << "\t" << PBy << "\t" << PBz_rf << "\n"; // non-zero and real off-shell mass
            continue;
        }
        // Off-shell mass of the bound cluster
        double Ma_off = sqrt(rrtt);
        // Total energies of "a" and "B" in the restframe of "A"
        double EaL = sqrt(Ma_off * Ma_off + Pa.Mag2());
        double EBL = sqrt(MB * MB + PB.Mag2());

        if (!is_direct)
        {
            //------- Lorentz transformations into laboratory system ---------
            std::pair<double, double> lora = Lorentz(gA, bA, EaL, Pa.Z());
            EaL = lora.first;     // cluster energy in lab
            Pa.SetZ(lora.second); // cluster Pz in lab

            std::pair<double, double> lorB = Lorentz(gA, bA, EBL, PB.Z());
            EBL = lorB.first;     // energy of the residual B in lab
            PB.SetZ(lorB.second); // Pz of the residual B in lab
            //---------- Generating CM scattering process ----------
            // S = Ma_off*Ma_off + Mi*Mi + 2*Mi*EaL; //Mandelstam invariant
        }

        Ea = EaL - Ma_off; // filling tree variable

        LVa.SetPxPyPzE(Pa.X(), Pa.Y(), Pa.Z(), EaL);
        TLorentzVector LVstart = LVa + LVi;
        double S = LVstart.Mag2(); // Mandelstam invariant
        // cout << "\nS=" <<  LVstart.Mag2();

        Mandelstam_S = S; // filling tree variable
        // Now generate CM scattering kinematics
        cm_values CM = CENMASS(S, Ma_off, Mi, Ma, is_isotropic, is_direct, (is_ppn || is_nnp));
        if (!CM.good)
        {
            // cout << "\nBad CM values!\n";
            continue; // non-physical output
        }
        double phi_rand = r1.Uniform(-PI, PI);

        P2cm.SetMag(CM.p_clust);
        P2cm.SetTheta(CM.theta_clust);
        P2cm.SetPhi(phi_rand);

        P1cm.SetX(-P2cm.X());
        P1cm.SetY(-P2cm.Y());
        P1cm.SetZ(-P2cm.Z());

        //-------  In the direction of the quasi-particle (cluster) in inverse kinematics--------
        double beta_cm;
        beta_cm = 0. - LVstart.Beta();

        double gamma_cm = 1 / sqrt(1 - beta_cm * beta_cm);

        std::pair<double, double> lora1 = Lorentz(gamma_cm, beta_cm, CM.e_scat, P1cm.Z());
        std::pair<double, double> lora2 = Lorentz(gamma_cm, beta_cm, CM.e_clust, P2cm.Z());

        P1cm.SetZ(lora1.second); // boost into lab frame
        P2cm.SetZ(lora2.second); // boost into lab frame
        //-------- Rotating back to the beam direction -----------
        TVector3 direction = LVstart.Vect();
        TVector3 P1L;
        TVector3 P2L;

        direction = direction.Unit();
        P1cm.RotateUz(direction);
        P2cm.RotateUz(direction);
        
        // For (p,2p) reactions, randomly assign indistinguishable particles
        if (is_p2p && r1.Uniform(0., 1.) > 0.5) {
            P1L = P2cm;
            P2L = P1cm;
        } else {
            P1L = P1cm;
            P2L = P2cm;
        }

        //--------- Filling in the ROOTTree variables------------
        theta_1 = P1L.Theta() * rad2deg;
        theta_2 = P2L.Theta() * rad2deg;
        theta_B = PB.Theta() * rad2deg;

        phi_1 = P1L.Phi() * rad2deg;
        phi_2 = P2L.Phi() * rad2deg;
        phi_B = PB.Phi() * rad2deg;

        P1x = P1L.X();
        P1y = P1L.Y();
        P1z = P1L.Z();
        P1_tot = P1L.Mag();

        P2x = P2L.X();
        P2y = P2L.Y();
        P2z = P2L.Z();
        P2_tot = P2L.Mag();

        PBx = PB.X();
        PBy = PB.Y();
        PBz_lab = PB.Z();

        E1 = sqrt(Mi*Mi + P1L.Mag2()) - Mi;
        E2 = sqrt(Ma*Ma + P2L.Mag2()) - Ma;
        EB = sqrt(MB*MB + PB.Mag2()) - MB;

        th1_cm  = CM.theta_scat * rad2deg;
        th2_cm  = CM.theta_clust * rad2deg; 
        P1_cm 	= CM.p_scat;
        P2_cm	= CM.p_clust;
        Moff 	= Ma_off;
        Mandelstam_T = CM.T;

        Opang=P1L.Angle(P2L) * rad2deg;

        Dif_phi = phi_1 - phi_2;
        if (Dif_phi < 0.)
            Dif_phi += 360;

        if (events % 10000 == 0)
            cout << "\r" << events << " of " << nevents << " (" << (float)events / nevents * 100 << "%)" << flush;
        
        tree->Fill();
        events++;
    }

    tree->AutoSave();
    file.Close();

    return;
}

// Kinematical function
double CINEMA(double x, double y, double z)
{
    double lambda = x * x + y * y + z * z - 2 * x * y - 2 * x * z - 2 * y * z;
    return lambda;
}

// Calculate elastic scattering kinematics in CM-system (1-target proton, 2-cluster)
cm_values CENMASS(double s, double m2off, double m1, double m2, bool is_isotropic_scattering, bool is_direct, bool is_neutron)
{
    cm_values output;
    output.good = false;
    double X = s;
    double Y = m2off * m2off;
    double Z = m1 * m1;
    double sqrs = sqrt(s);

    // Kinematics before the scattering process
    // (with one off-shell mass)
    double p2_off = sqrt(CINEMA(X, Y, Z)) / 2 / sqrs;
    double p1_off = p2_off;
    // CM energies
    double e1_off = (s + Z - Y) / 2 / sqrs;
    double e2_off = (s + Y - Z) / 2 / sqrs;

    // Now take the real masses (after scattering)
    Y = m2 * m2;
    Z = m1 * m1;
    // And check whether the kinematical function is ok
    // for this specific kinematical case
    double ERROR_CI = CINEMA(X, Y, Z);
    if (ERROR_CI <= 0.)
    {
        return output;
    }

    // Kinematics after the scattering process
    // (with all real masses)
    double p2 = sqrt(CINEMA(X, Y, Z)) / 2 / sqrs;
    double p1 = p2;
    double e1 = (s + Z - Y) / 2 / sqrs;
    double e2 = (s + Y - Z) / 2 / sqrs;

    // Let's consider momentum transfer from the
    // particle i to the cluster a
    double tmax = 2 * (m1 * m1 - e1_off * e1 - p1_off * p1); // COSINE=(-1)
    double tmin = 2 * (m1 * m1 - e1_off * e1 + p1_off * p1); // COSINE=(1)
    //cout << "\n\n tmax=" << tmax << "\t tmin=" << tmin;    

    double t;
    double COSINE;

    if (!is_isotropic_scattering)
    {
        // t = get_T(s,tmax);
        COSINE = get_CosThetaCM_Cugnon(s, is_neutron);
        
        if (is_direct)
            t = 2 * m1 * m1 - 2 * e1_off * e1 + 2 * p1_off * p1 * COSINE;
        else
            t = m2off * m2off + m2 * m2 - 2 * e2_off * e2 + 2 * p2_off * p2 * COSINE;

        if(fabs(t)>fabs(tmax) || fabs(t)<fabs(tmin))
        {
            cout << "\nNon physical t! t=" << t << "\t tmax=" << tmax << "\t tmin=" << tmin << "\n";
            return output;
        }         
        //double COSINE_check = (t - 2 * m1 * m1 + 2 * e1_off * e1) / (2 * p1_off * p1);
        //if (fabs(COSINE_check) > 1 || fabs(t) > fabs(tmax)) return output;
    }
    else 
    {
        t = gRandom->Uniform(tmax, 0); // Isotropic scattering
        COSINE = (t - 2 * m1 * m1 + 2 * e1_off * e1) / (2 * p1_off * p1);
    }
    if (fabs(COSINE) > 1)
    { // momentum transfer out of range
        // cout << "\nNon physical cos(theta cm)\n";
        // cout << "\nCOSINE=" << COSINE << "\n";
        return output;
    }
    // CM scattering angles
    double theta1, theta2;
    if(is_direct)
    {
        theta1 = acos(COSINE);
        theta2 = PI - theta1;
    }
    else
    {
        theta2 = acos(COSINE);
        theta1 = PI - theta2;
    }
    output.e_clust = e2;
    output.p_clust = p2;
    output.theta_clust = theta2;

    output.e_scat = e1;
    output.p_scat = p1;
    output.theta_scat = theta1;

    output.T = t;
    output.good = true;

    return output;
}

std::pair<double, double> Lorentz(double g, double b, double e, double p)
{
    double eL = g * e - g * b * p;
    double pL = g * p - g * b * e;
    return std::make_pair(eL, pL);
}

// Returns a random value of mandelstam T (in (MeV/c)^2 units)
// distributed according to the customly parameterized proton-proton
// cross section. Pass as a parameter "sm" the Mandelstam variable S)
// and the maximum possible momentum transfer
double get_T(double sm, double max)
{
    double Tmax = max * 0.000001; // convert to GeV� units

    Double_t rr = gRandom->Uniform(-1., 1.); // to randomize wrt 90 degrees
    double mandels = sm * 0.000001;          // in GeV^2

    // Probability function from the parameterization
    TF1 foo("foo", "[0]*exp(x*[1])*(1+0.02*exp((-6)*x))", Tmax / 2, 0);
    double c = 0.;
    if (mandels <= 4.79)
        c = -3283.75 + 3064.11 * mandels - 1068.44 * mandels * mandels + 164.844 * pow(mandels, 3) - 9.48152 * pow(mandels, 4);
    else if (mandels > 4.79)
        c = -776.822 + 586.016 * mandels - 175.347 * mandels * mandels + 26.1823 * pow(mandels, 3) - 1.94889 * pow(mandels, 4) + 0.0578352 * pow(mandels, 5);

    foo.FixParameter(0, 25.); // normalization constant (could be anything)
    foo.FixParameter(1, c);

    double Trand = foo.GetRandom(Tmax / 2, 0.); // from 90 to 0 degrees
    if (rr > 0)
        Trand = Tmax - Trand; // symmetrization relative to 90 degrees

    return (Trand * 1000000); // returning value in MeV^2
}

// Using parameterization dsigma/dt
// from Cugnon et al.  NIM in Physics Research B 111(1996) 215-220
double get_CosThetaCM_Cugnon(double sm, bool pn_flag)
{
    double mandels = sm * 0.000001;  // in GeV^2
    double mp = PROTON_MASS * 0.001; // Gev/c2
    double mt = mp;
    if (pn_flag) mt = NEUTRON_MASS * 0.001; // Gev/c2

    // Get equivalent energy for free proton
    double Tfree = (mandels - pow((mp+mt), 2)) / 2. / mt;
    double Pfree = sqrt(Tfree * (2 * mp + Tfree));

    //double Pcm = (1. / 2. / sqrt(mandels)) * sqrt(mandels * (mandels - 4 * mp * mp));
    double Pcm = (1. / 2. / sqrt(mandels)) * 
            sqrt( ( mandels - pow((mt-mp),2) ) * ( mandels - pow((mp+mt),2) ) );

    double Tmax = (-4) * Pcm * Pcm;
    double cos_theta_cm = 0.;
    double Trand = 0.;
    double Bpp = 0.;

    // Main parameter:
    if (Pfree < 2)
        Bpp = 5.5 * pow(Pfree, 8) / (7.7 + pow(Pfree, 8));
    else
        Bpp = 5.334 + 0.67 * (Pfree - 2);

    if (!pn_flag) //---> Parameterization for pp scattering
    {
        // Probability function from the parameterization
        TF1 foo("foo", "[0]*exp([1]*x)", Tmax, 0);
        foo.FixParameter(0, 1); // normalized to unity
        foo.FixParameter(1, Bpp);
        Trand = foo.GetRandom(Tmax, 0.);
    }
    else //--> Parameterization for pn scattering
    {
        double Bnp = 0.;
        double a=0;
        if (Pfree < 0.225)
        {
            Bnp = 0;
        }
        else if (0.225 <= Pfree && Pfree < 0.6)
        {
            Bnp = 16.53 * (Pfree - 0.225);
        }
        else if (0.6 <= Pfree && Pfree < 1.6)
        {
            Bnp = -1.63 * Pfree + 7.16;
        }
        else
        {
            Bnp = Bpp;
        }
        
        if (Pfree > 0.8)
        {
            a = 0.64 / pow(Pfree, 2);
        }
        else
        {
            a = 1;
        }
        TF1 *foo = new TF1("foo", "exp([0]*x)+[1]*exp([0]*([2]-[3]-x))", Tmax, 0);
        foo->FixParameter(0, Bnp);
        foo->FixParameter(1, a);
        foo->FixParameter(2, 2 * (mt * mt + mp * mp));
        foo->FixParameter(3, mandels);
        Trand = foo->GetRandom(Tmax, 0);
    }
    cos_theta_cm = Trand / 2. / Pcm / Pcm + 1;
    return cos_theta_cm; // returning cos of theta in CM system
}

// Helper function to check if argument matches a flag
bool hasFlag(const Char_t *arg, const Char_t *flag)
{
    return (strcmp(arg, flag) == 0);
}

// Helper function to extract value from --key=value format
bool getKeyValue(const Char_t *arg, const Char_t *key, Char_t *&value)
{
    size_t keylen = strlen(key);
    if (strncmp(arg, key, keylen) == 0 && arg[keylen] == '=')
    {
        value = (Char_t *)(arg + keylen + 1);
        return true;
    }
    return false;
}

int main(Int_t argc, Char_t *argv[])
{
    Char_t *root_filename = 0;

    int nevents = MAX_STORY;

    Bool_t NeedHelp = kTRUE;
    Bool_t is_ppn = kFALSE;
    Bool_t is_nnp = kFALSE;
    Bool_t is_ppa = kFALSE;
    Bool_t is_p2p = kTRUE;
    Bool_t is_isotropic = kTRUE;
    Bool_t is_direct = kFALSE;

    if (argc > 1)
    {
        for (Int_t i = 1; i < argc; i++)
        {
            Char_t *value = nullptr;

            if (getKeyValue(argv[i], "--root", value))
            {
                root_filename = value;
                NeedHelp = kFALSE;
            }
            if (getKeyValue(argv[i], "--max-events", value))
            {
                nevents = atoi(value);
            }
            if (hasFlag(argv[i], "--ppn"))
            {
                is_ppn = kTRUE;
            }
            else if (hasFlag(argv[i], "--nnp"))
            {
                is_nnp = kTRUE;
            }
            else
            {
                is_p2p = kTRUE;
            }
            if (hasFlag(argv[i], "--direct-kinematics"))
            {
                is_direct = kTRUE;
                std::cout << "\n-- Assuming direct kinematics";
            }
            if (hasFlag(argv[i], "--cugnon"))
            {
                is_isotropic = kFALSE;
                std::cout << "\n-- Generate anisotropic nucleon-distributions";
            }
            if (hasFlag(argv[i], "-h") || hasFlag(argv[i], "--help"))
            {
                NeedHelp = kTRUE;
            }
        }
    }

    if (NeedHelp)
    {
        std::cout << "\nQFS Kinematics Generator - Usage:\n";
        std::cout << "  ./qfs --root=<filename> [options]\n\n";
        std::cout << "Required Arguments:\n";
        std::cout << "  --root=<file>        Output ROOT file name\n\n";
        std::cout << "Optional Arguments:\n";
        std::cout << "  --max-events=<N>     Number of events to generate (default: 100000)\n";
        std::cout << "  --direct-kinematics  Use direct kinematics (default: inverse)\n";
        std::cout << "  --ppn                Generate (p,pn) reactions (default)\n";
        std::cout << "  --nnp                Generate (n,np) reactions\n";
        std::cout << "  --cugnon             Use Cugnon parameterization for NN elastic scattering(default: isotropic)\n";
        std::cout << "  -h, --help           Display this help message\n\n";
        std::cout << "Examples:\n";
        std::cout << "  ./qfs --root=output.root --max-events=50000 --cugnon\n";
        std::cout << "  ./qfs --root=output.root --direct-kinematics --ppn\n\n";

        return 0;
    }

    // Main generator function
    if (root_filename != NULL)
    {
        run(root_filename,
            nevents,
            is_ppn,
            is_isotropic,
            is_direct,
            is_nnp);
    }

    else
    {
        std::cout << "\nERROR! Output Root file is not specified. \nType: ./qfs \n\n"
                  << std::endl;
        return 0;
    }

    return 0;
}
