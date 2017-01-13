/*
 
This program simulates the propagation of Stokes vectors in a turbid
medium.
 
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
#include "MCDetector.h"
#include "MCSource.h"
#include "MCRandomScatterer.h"

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

        double Io, Qo, Uo, Vo;
        Io=1;
        Qo=1;
        Uo=0;
        Vo=0;
        
        string outputFilename;
        InitVariable(dict, "outputFilename",outputFilename);
 
        InitVariable(dict, "N_photon", N_photon);
        
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

        double wavelength;
        double epithelium_thickness, epithelium_scatterer_radius, epithelium_scatterer_index, epithelium_background_index,epithelium_mu_s, epithelium_mu_a;
        double precancer_thickness, precancer_scatterer_radius, precancer_scatterer_index, precancer_background_index,precancer_mu_s, precancer_mu_a;
        double stroma_thickness, stroma_scatterer_radius, stroma_scatterer_index, stroma_background_index,stroma_mu_s, stroma_mu_a;

        long Nx, Ny, acceptanceCosineElements;
        double imgWidth, imgHeight;
        long Vol_Nx, Vol_Ny, Vol_Nz;
        double Vol_XMin, Vol_XMax, Vol_YMin, Vol_YMax, Vol_ZMin, Vol_ZMax;
        
        InitVariable(dict, "wavelength", wavelength);
        InitVariable(dict, "Nx", Nx);
        InitVariable(dict, "Ny", Ny);
        InitVariable(dict, "imgWidth", imgWidth);
        InitVariable(dict, "imgHeight", imgHeight);
        InitVariable(dict, "acceptanceCosineElements", acceptanceCosineElements);
        
        InitVariable(dict, "epithelium_thickness", epithelium_thickness);
        InitVariable(dict, "epithelium_mu_s", epithelium_mu_s);
        InitVariable(dict, "epithelium_mu_a", epithelium_mu_a);
        InitVariable(dict, "epithelium_scatterer_radius", epithelium_scatterer_radius);
        InitVariable(dict, "epithelium_scatterer_index", epithelium_scatterer_index);
        InitVariable(dict, "epithelium_background_index", epithelium_background_index);

        InitVariable(dict, "precancer_thickness", precancer_thickness);
        InitVariable(dict, "precancer_mu_s", precancer_mu_s);
        InitVariable(dict, "precancer_mu_a", precancer_mu_a);
        InitVariable(dict, "precancer_scatterer_radius", precancer_scatterer_radius);
        InitVariable(dict, "precancer_scatterer_index", precancer_scatterer_index);
        InitVariable(dict, "precancer_background_index", precancer_background_index);

        InitVariable(dict, "stroma_thickness", stroma_thickness);
        InitVariable(dict, "stroma_mu_s", stroma_mu_s);
        InitVariable(dict, "stroma_mu_a", stroma_mu_a);
        InitVariable(dict, "stroma_scatterer_radius", stroma_scatterer_radius);
        InitVariable(dict, "stroma_scatterer_index", stroma_scatterer_index);
        InitVariable(dict, "stroma_background_index", stroma_background_index);

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

		RealV detectorPosition;
		long detectorPixels;
		double detectorSize;
		double aroundYInDegrees;
        
		if (VariableExists(dict, "detectorPosition")) {
			InitVariable(dict, "detectorPosition", detectorPosition);
			InitVariable(dict, "detectorSize", detectorSize);
			InitVariable(dict, "detectorPixels", detectorPixels);
			InitVariable(dict, "aroundYInDegrees", aroundYInDegrees);
			
			MCDetector* detector = new MCDetector(detectorSize, detectorPixels, 0, 0);
            detector->RotateObject(0, aroundYInDegrees*PI/180., 0);
			world.PlaceIntoObject(detector, detectorPosition);
		} else {
			throw runtime_error("No external detector defined");
		}
		
        MCRandomScattererJaillon epithelium(epithelium_scatterer_radius, epithelium_scatterer_index, epithelium_background_index, wavelength, 10000, epithelium_mu_s, epithelium_mu_a, 0);
        MCRandomScattererJaillon precancer(precancer_scatterer_radius, precancer_scatterer_index, precancer_background_index, wavelength, 10000, precancer_mu_s, precancer_mu_a, 0);
        MCRandomScattererJaillon stroma(stroma_scatterer_radius, stroma_scatterer_index, stroma_background_index, wavelength, 10000, stroma_mu_s, stroma_mu_a, 0);
        RealV origin;
        MCInfiniteLayer* first,* second, *third;

        first = new MCInfiniteLayer(epithelium_thickness, imgWidth, imgHeight, Nx, Ny, epithelium_background_index, 1, 0, kFromInside, acceptanceCosineElements, 0);
        first->SetRandomScatterer(&epithelium);

        second = new MCInfiniteLayer(precancer_thickness, imgWidth, imgHeight, Nx, Ny, precancer_background_index, 1, 0, kFromInside, acceptanceCosineElements, 0);
        second->SetRandomScatterer(&precancer);

        third = new MCInfiniteLayer(stroma_thickness, imgWidth, imgHeight, Nx, Ny, stroma_background_index, 1, 0, kFromInside, acceptanceCosineElements, 0);
        third->SetRandomScatterer(&stroma);

		
        if (epithelium_thickness == 0 ) {
            origin = RealV(0,0,0);
            world.PlaceIntoObject(second, origin);

            origin = RealV(0,0, precancer_thickness);
            world.PlaceIntoObject(third, origin);
            third->InterfaceTouchesOtherObject(second, kForwardZPlane, kBackwardZPlane);
        } else if (precancer_thickness == 0) {
            origin = RealV(0,0,0);
            world.PlaceIntoObject(first, origin);

            origin = RealV(0,0, epithelium_thickness);
            world.PlaceIntoObject(third, origin);
            third->InterfaceTouchesOtherObject(first, kForwardZPlane, kBackwardZPlane);
        } else {
            origin = RealV(0,0,0);
            world.PlaceIntoObject(first, origin);

            origin = RealV(0,0, epithelium_thickness);
            world.PlaceIntoObject(second, origin);
            second->InterfaceTouchesOtherObject(first, kForwardZPlane, kBackwardZPlane);

            origin = RealV(0,0, epithelium_thickness+precancer_thickness);
            world.PlaceIntoObject(third, origin);
            third->InterfaceTouchesOtherObject(second, kForwardZPlane, kBackwardZPlane);
        }

		/* 2.2 Creation of photon source */
		auto_ptr<MCSource> theSource;
		
		RealV dir, pos;
		
	    if ( VariableExists(dict, "direction") ) {
		    InitVariable(dict, "direction",dir);
			dir.normalize();
		} else {
			dir = RealV(0,0,1);	
		}
		
	    if ( VariableExists(dict, "source") ) {
		    InitVariable(dict, "source",pos);
		} else {
			pos = RealV(0,0,-10);	
		}
		
		
        if (VariableExists(dict, "inputBeamType")) {
            InitVariable(dict, "inputBeamType", inputBeamType);
            if (inputBeamType == string("pencil")) {
                clog << RightNow() << "Using pencil beam\n";
				inputBeamDiameter = 0;
			} else if (inputBeamType == string("disk")) {
                InitVariable(dict, "inputBeamDiameter", inputBeamDiameter);
                clog << RightNow() << "Using disk beam centered with diameter " << inputBeamDiameter << endl;
            }
            theSource.reset(new MCLaserSource(StokesV(Io,Qo,Uo,Vo, RealV(0,1,0), dir), inputBeamDiameter/2., 0, pos, dir));
        } else {
			throw runtime_error("You must define the type of light source you want to use");
        }

        theSource->SetContainerObject(&world);
        
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
            auto_ptr<Photon> thePhoton;
			
            switch (formalismID) {
                case kIntensity:
                    thePhoton.reset(new PhotonIntensity);
                    break;
                case kBartel:
                    thePhoton.reset(new PhotonBartel);
                    break;
                case kMourant:
                    thePhoton.reset(new PhotonMourant);
                    break;
                case kCote:
                    thePhoton.reset(new PhotonCote);
                    break;
                case kJaillon:
                    thePhoton.reset(new PhotonJaillon);
                    break;
                default:
                    throw runtime_error("Undefined formalism: "+formalism);
            }
			
			
			
            /* 3.3) This call below is the guts of the calculation.
                Refer to MCObject::PropagateInObject() for details */
			
            // Define keepPhotonStats in parameter file for full stats (very slow)
            if (VariableExists(dict, "keepPhotonStats"))
                thePhoton->KeepStats();
			
			theSource->InitializePhoton(thePhoton.get());
            world.PropagateInObject(thePhoton.get());
			
            // If keepPhotonStats is defined in parameter file, dump to stdout
            if (VariableExists(dict, "dumpAllPhotonStats")) {
                cout.setf(ios::fixed);
				
                thePhoton->DumpStats(cout, kMathematica);
                cout << ",\n";
                cout.unsetf(ios::fixed);
            }
			
            if (gTryToQuitNicely) {
                // User has tried to quit
                clog << RightNow() << "User asked to quit.  Quitting nicely." << endl;
                break;
            }
			
            if (gSaveCurrentStateToDisk) {
                // User asked for a safety-save
                clog << RightNow() << "User asked to save current state." << endl;
                world.DumpToFile(outputFilename+string(".save"));
                cerr << "Done.\n";
                gSaveCurrentStateToDisk = false;
            }
			
            /*
             
             // Here one can put some functions to analyze the stats and dump
             // them to file.  This is an example I sometimes need:
			 
			 StokesV s;
			 thePhoton->GetStokesVectorInLocalFrame(s);
			 
             RealV pos;
             RealV x(1,0,0);
             RealV y(0,1,0);
             RealV z(0,0,1);
			 
             thePhoton->MeasureStokesVectorInLabFrame(s, y,z,x) ;
             double ipara,iperp;
             thePhoton->IntensityThroughLinearPolarizer(s, y,z,x, ipara,iperp);
             
             pos = thePhoton->GetPosition();
             if (pos.z >= 0.99  && abs(pos.x) < 0.05 && abs(pos.y) < 0.05) {
                 cout << "Contribution: " << ipara << endl;
                 thePhoton->DumpStats(cout);
             }
             */
            
        } // All photons
		
   		gStatus = kWritingDataToDisk;
		
        /* 4) Store everything on disk */
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

