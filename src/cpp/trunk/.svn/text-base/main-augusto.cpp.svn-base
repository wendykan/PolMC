
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
#include "MCCylinder.h"
#include "MCDetector.h"
#include "MCInfiniteLayers.h"
#include "MCSource.h"
#include "randomscatter.h"

#ifndef MAXPATHLEN
#define MAXPATHLEN 1000
#endif

using namespace std;

enum {kPencilBeam = 0, kDiskBeam};

enum {
	kInitializing,
	kPropagatingPhoton,
	kWritingDataToDisk,
	
};

bool gTryToQuitNicely = false;
bool gSaveCurrentStateToDisk = false;
int	gStatus = kInitializing;

/* This is a very low-level Unix function to
interact nicely with user.


Defined below. */
void catch_signal(int sig_num);


int
main(int argc, char *argv[])
{
    
    int ch;
	char logfile[200];
	
	logfile[0] = 0;
	while ((ch = getopt(argc, argv, "f:")) != -1)
		switch (ch) {
			case 'f':
				strncpy(logfile, optarg, 200);
				break;
			case '?':
			default:
				printf("Usage incorrect: polmc -f logfile parameterfile");
                exit(1);
		}
            argc -= optind;
	argv += optind;
    
    /* Before anything, we redirect clog to a log file if asked by the user */
    
    ofstream myclog(logfile,ios_base::app);
    if (strlen(logfile) != 0) 
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
        long unwrappertimer  = 0;
        
        string formalism;
        string randomSampling;
        string objectType;
        string inputBeamType;
        long inputBeamTypeID;
        long formalismID;
        long randomSamplingID;
        double inputBeamDiameter;
        
        vector<StokesV> incidentStokesV;
        string outputFilename;
        
        /* 1) Read parameters necessary to run the simulations.
            After this, dict contains all the parameters for the simulation.
            Below, it will be passed to the scattering object so it can retrieve its
            parameters too.  */
        
		gStatus = kInitializing;
        
        if (argc == 1) {
            clog << RightNow() << "Reading parameters from file " << argv[argc-1] << endl ;
            ifstream fInput(argv[argc-1], std::ios_base::in);
            ReadParametersFromStream(fInput, dict);
        } else {
            clog << RightNow() << "Reading parameters from stdin\n" ;
            ReadParametersFromStream(cin, dict);
        }
        
        /* 2) Initialize variables for Monte Carlo calculation  */
		
		string debugLevel;
		if (VariableExists(dict, "debugLevel"))   
			InitVariable(dict, "debugLevel", debugLevel);
        
		if (debugLevel == "Normal") {
			gDebugLevel = kNormal;
			clog << RightNow() << "Debug level is normal\n";
		} else if (debugLevel == "Verbose") {
			gDebugLevel = kVerbose;
			clog << RightNow() << "Debug level is verbose\n";
		} else if (debugLevel == "ExtremelyVerbose") {
			gDebugLevel = kExtremelyVerbose;
			clog << RightNow() << "Debug level is extremely verbose\n";
		}
		
        InitVariable(dict, "outputFilename",outputFilename);
        
        InitVariable(dict, "N_photon", N_photon);
        
        formalismID = kIntensity;
        
        /*  The object responsible for giving us random angles and distances
            Without the scatmech library, one can only use the Henyey Greenstein function */
        randomSamplingID = kHenyeyGreenstein;
        
        MCRandomScatterer* randomScatterer;
        double anisotropy, mu_s, mu_a;

        if (! VariableExists(dict, "g") ) {
            clog <<  RightNow() << "The variable 'g' (the anisotropy factor) must be defined to determine distribution.\n";
            throw runtime_error("anisotropy factor 'g' not defined");
        } else {
            InitVariable(dict, "g", anisotropy);
            InitVariable(dict, "mu_s", mu_s);
            InitVariable(dict, "mu_a", mu_a);
        }
        
        randomScatterer = new MCRandomScattererHenyeyGreenstein(anisotropy, mu_s, mu_a);

        
        double wavelength;
        long Vol_Nx, Vol_Ny, Vol_Nz;
        double Vol_XMin, Vol_XMax, Vol_YMin, Vol_YMax, Vol_ZMin, Vol_ZMax;
        
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
        
        MCObject* sample;
        
		/* Most objects use these, so they are declared here */
		RealV origin, size;
		double index_outside, index_medium, rotationPerCmInClearSpace;
		long Nx, Ny, Nz, acceptanceCosineElements;
		InitVariable(dict, "origin", origin);
		InitVariable(dict, "size", size);
		InitVariable(dict, "index_outside", index_outside);
		InitVariable(dict, "index_med", index_medium);
		
		if (VariableExists(dict, "acceptanceCosElements")) {
			InitVariable(dict, "acceptanceCosElements",acceptanceCosineElements);
		} else {
			acceptanceCosineElements = 1;
		}
		
		
        string tableForZ, tableForTheta, tableForPhi;
        fastinterpolate zProb;
        fastinterpolate thetaProb;
        fastinterpolate phiProb;
        
        RealV diffuserPosition;
        RealV diffuserDirection;
        if (VariableExists(dict, "tableForZ")) {
            InitVariable(dict, "tableForZ",tableForZ);
            InitVariable(dict, "tableForTheta",tableForTheta);
            InitVariable(dict, "tableForPhi",tableForPhi);
            InitVariable(dict, "diffuserPosition",diffuserPosition);
            InitVariable(dict, "diffuserDirection",diffuserDirection);
        } else {
            throw runtime_error("Missing info for diffuser");
        }
        
        MCCylindricalDiffuser  diffuser(tableForZ,
                                        tableForTheta,
                                        tableForPhi,
                                        diffuserPosition,
                                        diffuserDirection);
    //Augusto included this line    
	InitVariable(dict, "objectType",objectType);
        
        if (objectType == "boxalone") {
            Nx = 1;
            Ny = 1;
            Nz = 1;
            sample = new MCBox(size.x, size.y, size.z, Nx, Ny, Nz, 
                                   index_medium, index_outside, rotationPerCmInClearSpace,
                                   kIntensityOnly | kLeaving | kEntering, acceptanceCosineElements, false );
            sample->SetRandomScatterer(randomScatterer);
            world.Place(sample, origin);
            diffuser.SetContainerObject(sample);
        } else if (objectType == "boxandcylinder") {
            Nx = 10;
            Ny = 10;
            Nz = 10;
            sample = new MCBox(size.x, size.y, size.z, Nx, Ny, Nz, 
                               index_medium, index_outside, rotationPerCmInClearSpace,
                               kIntensityOnly | kLeaving | kEntering, acceptanceCosineElements, false );
            sample->SetRandomScatterer(randomScatterer);

			double radius, height;
			long Nh, Nr;
            RealV cylinderPosition;
			InitVariable(dict, "Nh",Nh);
			InitVariable(dict, "Nr",Nr);
			InitVariable(dict, "radius",radius);
			InitVariable(dict, "height",height);
			InitVariable(dict, "cylinderPosition",cylinderPosition);
            MCCylinder* cyl;
            cyl = new MCCylinder(radius, height, Nr, Nh, 
                                       index_medium, index_outside, rotationPerCmInClearSpace,
                                    kIntensityOnly | kLeaving, acceptanceCosineElements, false ); 
            
/*            cyl = new MCBox(1.5,1.5,4, 1,1,1, 
                                 index_medium, index_outside, rotationPerCmInClearSpace,
                                 kFromInside, acceptanceCosineElements, false ); */
            cyl->SetOuterObjectTo(sample);
            cyl->SetRandomScatterer(new MCRandomScattererHenyeyGreenstein(0, 0, 0.01));
            
            world.Place(sample, origin);
            world.Place(cyl, cylinderPosition);
            diffuser.SetContainerObject(cyl);
        } else {
            throw runtime_error("Undefined objext type: "+objectType);
        }
        
        
        if (! world.IsGeometryConsistent() ) {
            throw runtime_error("Geometry inconsistent.");
        }
        
        /* 3) The main loop, going through all Stokes vectors
            and all photons */
        
		gStatus = kPropagatingPhoton;
        
        /* 3.1 ) Reset source parameters */
        clog << RightNow() << "Detection formalism demanded: " << formalism << endl;
        
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
            PhotonIntensity thePhoton;
            
            /* 3.3) This call below is the guts of the calculation.
                Refer to MCObject::PropagateInObject() for details */

            // Define keepPhotonStats in parameter file for full stats (very slow)
            if (VariableExists(dict, "keepPhotonStats")) {
                thePhoton.KeepStats();
                diffuser.KeepStats();
            }
            
			diffuser.InitializePhoton(&thePhoton);

            world.PropagateInObject(&thePhoton);
            
            // If keepPhotonStats is defined in parameter file, dump to stdout
            if (VariableExists(dict, "dumpAllPhotonStats")) {
                thePhoton.DumpStats(cout, kRaw);
            }
            
            if (gTryToQuitNicely) {
                // User has tried to quit
                clog << RightNow() << "User asked to quit.  Quitting nicely." << endl;
                break;
            }
            
            if (gSaveCurrentStateToDisk) {
                // User asked for a safety-save
                clog << RightNow() << "User asked to save current state." << endl;
                world.NormalizationFactorForInstensityStatistics(diffuser.GetNumberOfPhotonsLaunched());
                world.DumpToFile(outputFilename+string(".save"));
                cerr << "Done.\n";
                gSaveCurrentStateToDisk = false;
            }
            
        } // All photons

        // If dumpAllDiffuserStats is defined in parameter file, dump to stdout
        if (VariableExists(dict, "dumpAllDiffuserStats")) {
            diffuser.WriteStatsToStream(cout);
        }
        
   		gStatus = kWritingDataToDisk;
        
        /* 4) Store everything on disk */
        world.NormalizationFactorForInstensityStatistics(diffuser.GetNumberOfPhotonsLaunched());
        world.DumpToFile(outputFilename);
        
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
			if (gStatus == kInitializing) {
				cerr << "Giving up before initialization\n";
				exit(0);			
			} else {
	            gTryToQuitNicely = true;
	            cerr << "Trying to quit nicely" << endl;
			}
            break;
        case SIGUSR1:
       		if (gStatus != kWritingDataToDisk) {
                gSaveCurrentStateToDisk = true;
                cerr << "Trying to save." << endl;
			} else {
				cerr << "Already saving data" << endl;
			}
            break;
    }
}


