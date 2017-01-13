ance/*
 
This program simulates the propagation of Stokes vectors in a turbid
medium.

The original references upon which this code is based are:

"Monte Carlo simulations of the diffuse backscattering 
 Mueller matrix for highly scattering media"
 Appl. Opt. Vol 39, April 1st 2000, p. 1580
 
"Characterizing mammalian cells and cell phantoms by 
 polarized backscatterng fiber-optic measurements"
 J. Mourant, T. Johnson, J. Freyer
 Appl. Opt. Vol 40, No8 Oct 2001, p. 5114
 
 "Mueller matrix of dense polystyrene latex sphere suspensions:
 measurements and Monte Carlo simulation"
 B. Kaplan, G. Ledanois, B. Drevillon
 Appl. Opt. Vol 40, no 16, June 2001, p. 2769-2777

 which are themselves loosely based on code from Wang and Jacques

"MCML-Monte Carlo modeling of light transport" in Optical-Thermal
 response of laser irradiated tissue, A. J. Welch M. van Gemert
 pp. 73-100 http://omlc.ogi.edu/software/mc/index.html
 
*/

#if HAVE_CONFIG_H
    #include <config.h>
#endif


#include <math.h>
#include <iostream>
#include <fstream>

#if HAVE_STDEXCEPT
    #include <stdexcept>
#endif

#include <unistd.h>
#include <sys/param.h>

#include <string>
#include <cmath>
#include <ctime>
#include <limits>
#include <cstdarg>

#include <signal.h>
#include <vector>

#include "constants.h"
#include "rand.h"
#include "configfiles.h"
#include "cubicspline.h"
#include "mydebug.h"
#include "Photon.h"
#include "RealV.h"
#include "StokesV.h"
#include "MuellerM.h"
#include "UTimer.h"

#if HAVE_MPI_H
    #include "mpi.h"
#endif

#include "MCWorld.h"
#include "MCObject.h"
#include "MCBox.h"
#include "MCEllipsoid.h"
#include "MCInfiniteLayers.h"
#include "MCRandomScatterer.h"

#ifndef MAXPATHLEN
    #define MAXPATHLEN 1000
#endif
   
using namespace std;

enum {kPencilBeam = 0, kDiskBeam};

bool gTryToQuitNicely = false;
bool gSaveCurrentStateToDisk = false;

/* This is a very low-level Unix function to
   interact nicely with user.


   Defined below. */
void catch_signal(int sig_num);


int
main(int argc, char *argv[])
{
    /* Before anything, we redirect clog to a file called polmc.log */
    ofstream myclog("linda.polmc.log",ios_base::app);
    clog.rdbuf(myclog.rdbuf());

    char s[MAXPATHLEN];
    getcwd(s, MAXPATHLEN);
    clog << RightNow() << "Starting program" << endl;
    clog << RightNow() << "Current working directory: " << s << endl;

    try {
        /* We try to catch Unix signals: if the user types Ctrl-C
        or "kills" the program with the kill command on
        Unix (kill -TERM), we save our work before quitting.

        If a user send a "USR1" signal, we save to disk (and keep going).
        */
        signal(SIGINT, catch_signal);
        signal(SIGTERM, catch_signal);
        signal(SIGUSR1, catch_signal);

        /* 0) The variables for this program */

        map<string,string> dict;

        long N_photon;
        time_t time;
        
        string formalism;
        string randomSampling;
        string objectType;
        string inputBeamType;
        long inputBeamTypeID;
        long formalismID;
        long randomSamplingID;
        double inputBeamDiameter;
        double mu_s, mu_a;
        
        vector<StokesV> incidentStokesV;
        vector<string> outputFilenames;

        /* 1) Read parameters necessary to run the simulations.
            After this, dict contains all the parameters for the simulation.
            Below, it will be passed to the scattering object so it can retrieve its
            parameters too.  */

        if (argc > 1) {
            clog << RightNow() << "Reading parameters from file " << argv[1] << endl ;
            ifstream fInput(argv[1], std::ios_base::in);
            ReadParametersFromStream(fInput, dict);
        } else {
            clog << RightNow() << "Reading parameters from stdin\n" ;
            ReadParametersFromStream(cin, dict);
        }

        /* 2) Initialize variables for Monte Carlo calculation  */

        double Io, Qo, Uo, Vo;
        Io=1;
        Qo=1;
        Uo=0;
        Vo=0;
        
        string outputFilename;
        InitVariable(dict, "outputFilename",outputFilename);
        incidentStokesV.push_back( StokesV(Io,Qo,Uo,Vo) );
        outputFilenames.push_back( outputFilename);
 
        InitVariable(dict, "N_photon", N_photon);

        if (VariableExists(dict, "inputBeamType")) {
            InitVariable(dict, "inputBeamType", inputBeamType);

            if (inputBeamType == string("pencil")) {
                inputBeamTypeID = kPencilBeam;
                clog << RightNow() << "Using pencil beam incident at (0,0,0)\n";
            } else if (inputBeamType == string("disk")) {
                inputBeamTypeID = kDiskBeam;
                InitVariable(dict, "inputBeamDiameter", inputBeamDiameter);
                clog << RightNow() << "Using disk beam centered at (0,0,0) with diameter " << inputBeamDiameter << endl;
            }
        } else {
            inputBeamTypeID = kPencilBeam;
            clog << RightNow() << "Using pencil beam at (0,0,0)\n";
        }
        
        InitVariable(dict, "formalism",formalism);
        if (formalism == string("Intensity")) {
            formalismID = kIntensity;
        } else if (formalism == string("Bartel")) {
            formalismID = kBartel;
        } else if (formalism == string("Mourant")) {
            formalismID = kMourant;
        } else if (formalism == string("Cote")) {
            formalismID = kCote;
        } else if (formalism == string("Jaillon")) {
            formalismID = kJaillon;
        } else {
            throw runtime_error("Undefined formalism: "+formalism);
        }
        
        /* 2.0 Random scatterer */ 
        InitVariable(dict, "randomSampling",randomSampling);
        MCRandomScatterer* randomScatterer;
        double anisotropy;
        
        if (randomSampling == string("Bartel")) {
            randomSamplingID = kBartel;
        } else if (randomSampling == string("Henyey-Greenstein") || randomSampling == string("HG")) {
            randomSamplingID = kHenyeyGreenstein;
        } else if (randomSampling == string("Kaplan")) {
            randomSamplingID = kKaplan;
        } else if (randomSampling == string("Jaillon")) {
            randomSamplingID = kJaillon;
        }  else {
            throw runtime_error("Undefined random angle sampling: "+randomSampling);
        }


#if WITH_POLARIZATION
       /* 1.1) The object responsible for giving us random angles and distances */
        switch (randomSamplingID) {
            case kKaplan:
                randomScatterer = new MCRandomScattererKaplan(dict);
                break;
            case kHenyeyGreenstein:
                if (! VariableExists(dict, "g") ) {
                    clog <<  RightNow() << "The variable 'g' (the anisotropy) is not defined.\n";
                    randomScatterer = new MCRandomScattererHenyeyGreenstein(dict);
                } else {
                    double mu_s, mu_a;
                    InitVariable(dict, "g", anisotropy);
                    InitVariable(dict, "mu_s", mu_s);
                    InitVariable(dict, "mu_a", mu_a);
                    randomScatterer = new MCRandomScattererHenyeyGreenstein(anisotropy, mu_s, mu_a);
                }

                break;
            case kJaillon:
                randomScatterer = new MCRandomScattererJaillon(dict);
                break;
        }

#else
        if (randomSamplingID != kHenyeyGreenstein) {
            clog <<  RightNow() << "Program compiled without the SCATMECH library: only Henyey-Greenstein random sampling allowed\n";
            throw runtime_error("must use Henyey-Greenstein distribution function for scattering");
        }
        if (! VariableExists(dict, "g") ) {
            clog <<  RightNow() << "The variable 'g' (the anisotropy factor) must be defined to determine distribution.\n";
            throw runtime_error("anisotropy factor 'g' not defined");
        } else {
            InitVariable(dict, "g", anisotropy);

        }

        /* 1.1) The object responsible for giving us random angles and distances
            Without the scatmech library, one can only use the Henyey Greenstein function */

        randomScatterer = new MCRandomScattererHenyeyGreenstein(anisotropy, 0, 0);

#endif

        double wavelength;

        long Nx, Ny, Nz, acceptanceCosineElements;
        double xwidth, yheight, zdepth;
        long Vol_Nx, Vol_Ny, Vol_Nz;
        double Vol_XMin, Vol_XMax, Vol_YMin, Vol_YMax, Vol_ZMin, Vol_ZMax;
        
        InitVariable(dict, "wavelength", wavelength);
        InitVariable(dict, "Nx", Nx);
        InitVariable(dict, "Ny", Ny);
        InitVariable(dict, "Nz", Nz);
        InitVariable(dict, "zdepth",zdepth);
        InitVariable(dict, "xwidth",xwidth);
        InitVariable(dict, "yheight",yheight);
        InitVariable(dict, "acceptanceCosineElements", acceptanceCosineElements);
        
        InitVariable(dict, "Vol_Nx",Vol_Nx);
        InitVariable(dict, "Vol_Ny",Vol_Ny);
        InitVariable(dict, "Vol_Nz",Vol_Nz);
        InitVariable(dict, "Vol_XMin",Vol_XMin);
        InitVariable(dict, "Vol_YMin",Vol_YMin);
        InitVariable(dict, "Vol_ZMin",Vol_ZMin);
        InitVariable(dict, "Vol_XMax",Vol_XMax);
        InitVariable(dict, "Vol_YMax",Vol_YMax);
        InitVariable(dict, "Vol_ZMax",Vol_ZMax);

        /* 2.1) Creation of the object through which photons will propagate
            All the parameters necessary for constructions are in dict */
        MCWorld world(Vol_XMin, Vol_XMax, Vol_YMin, Vol_YMax, Vol_ZMin, Vol_ZMax, Vol_Nx, Vol_Ny, Vol_Nz);
        auto_ptr<MCObject> sample;

        sample.reset(new MCBox(1,1,1, 10,10,10, 1.33, 1, 0, 10, 0));
        sample->SetRandomScatterer(randomScatterer);
        RealV origin;
        world.PlaceIntoObject(sample.get(), origin);


        if (! world.IsGeometryConsistent() ) {
            throw runtime_error("Geometry inconsistent.");
        }
        
        
        /* 3) The main loop, going through all Stokes vectors
            and all photons */
        
        for (long s = 0; s < incidentStokesV.size(); s++) {
            /* 3.1 ) Reset simple internal variables (arrays) */
            clog << RightNow() << "Detection formalism demanded: " << formalism << endl;
            clog << RightNow() << "Entering Stokes vector loop for S(I,Q,U,V)=" << incidentStokesV.at(s) << endl;

            long l = 10;

            /* 3.2) Loop through all photons */
            UTimer timer;

            for (long photon = 1; photon <= N_photon; photon++) {
                
                if (photon % l == 0)    {
                    clog << RightNow() << "Photon [" << photon << "/" << N_photon
                        << "] in " << timer.DurationInText(timer.GetTicks()) << ", expect end in " << timer.DurationInText(TimeT(timer.GetTicks()*(N_photon/photon-1))) << endl;
                    clog.flush();
                    l *= 10;
                }

                // Deleted automatically after end of scope
                auto_ptr<Photon> thePhoton;

                switch (formalismID) {
                    case kIntensity:
                        thePhoton.reset(new PhotonIntensity(incidentStokesV.at(s).mI));
                        break;
                    case kBartel:
                        thePhoton.reset(new PhotonBartel(incidentStokesV.at(s)));
                        break;
                    case kMourant:
                        thePhoton.reset(new PhotonMourant(incidentStokesV.at(s)));
                        break;
                    case kCote:
                        thePhoton.reset(new PhotonCote(incidentStokesV.at(s)));
                        break;
                    case kJaillon:
                        thePhoton.reset(new PhotonJaillon(incidentStokesV.at(s)));
                        break;
                    default:
                        throw runtime_error("Undefined formalism: "+formalism);
                }

                thePhoton->SetGlobalOrigin(sample->GetGlobalOrigin(), false);
                thePhoton->SetCurrentObject(sample.get());
                RealV pos;
                InitVariable(dict, "originx",pos.x);
                InitVariable(dict, "originy",pos.y);
                InitVariable(dict, "originz",pos.z);
                thePhoton->SetLocalPosition(pos);
                double phi = 2*PI*MCRandomScatterer::RandomFloat();
                double theta = acos(1. - 2. * MCRandomScatterer::RandomFloat());
                RealV d(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
                thePhoton->SetPropagationDirectionInLabFrame(d);

                /* 3.3) This call below is the guts of the calculation.
                    Refer to MCObject::PropagateInObject() for details */

                // Define keepPhotonStats in parameter file for full stats (very slow)
                if (VariableExists(dict, "keepPhotonStats"))
                    thePhoton->KeepStats();

                world.PropagateInObject(thePhoton.get());

                // If keepPhotonStats is defined in parameter file, dump to stdout
                if (VariableExists(dict, "keepPhotonStats")) {
                   thePhoton->DumpStats(cout);
                }
                
                if (gTryToQuitNicely) {
                    // User has tried to quit
                    clog << RightNow() << "User asked to quit.  Quitting nicely." << endl;
                    break;
                }

                if (gSaveCurrentStateToDisk) {
                    // User asked for a safety-save
                    clog << RightNow() << "User asked to save current state." << endl;
                    world.DumpToFile(outputFilenames.at(s)+string(".save"));
                    gSaveCurrentStateToDisk = false;
                }

            } // All photons

            /* 4) Store everything on disk */
            world.DumpToFile(outputFilenames.at(s));

            if (gTryToQuitNicely) {
                // User has tried to quit
                break;
            }
            
        } // All Stokes vectors

    } catch (exception& e) {
        clog << RightNow() << "Exception caught: " << string(e.what()) << endl;
        clog << RightNow() << "Program ended abnormally" << endl;
        return 1;
    }

    if (! gTryToQuitNicely) {
        clog << RightNow() << "Program ended normally" << endl;
    } else {
        clog << RightNow() << "Program was required to quit.  Files were written to disk with the stats as they were at the moment of the user asked to quit." << endl;
    }
    

    return 0;
}



/* This function allows the program to avoid getting killed
   without saving (if possible).  It is also possible
   to ask to program to save its current state. */

 void catch_signal(int sig_num)
 {
     /* Reset signal handler.  Must be done, don't ask. */
     signal(sig_num, catch_signal);

     switch (sig_num) {
         case SIGINT:
         case SIGTERM:
             gTryToQuitNicely = true;
             cerr << "Quitting. Please wait while files are saved." << endl;
             break;
         case SIGUSR1:
             gSaveCurrentStateToDisk = true;
             cerr << "Saving." << endl;
             break;
     }
 }
 

