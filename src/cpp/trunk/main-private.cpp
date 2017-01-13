/*
 
This program simulates the propagation of Stokes vectors in a turbid
medium.

The original references upon which this code is based are:

"Monte Carlo simulations ¥ the diffuse backscattering 
 Mueller matrix for highly scattering media"
 Appl. Opt. Vol 39x, April 1st 2000, p. 1580
 
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
#include "MCMaterials.h"

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
#include "MCGenericObject.h"
#include "MCSource.h"
#include "MCRandomScatterer.h"
#include "MCUtils.h"
#include "CFHelper.h"

//#include <Accelerate.h>

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
bool gDumpProximityStats = false;
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
    char buf[255];
	
	map<string,string> dict;

    /* Default values out of range */
    double anisotropy = -2, mu_s = -1, mu_a = -1, wavelength = 0;
    bool forceEvenWithInconsistentGeometry = false;
	string outputFilename("Output.xml");
	
	logfile[0] = 0;
    buf[0] = 0;
	while ((ch = getopt(argc, argv, "o:w:fl:u:a:g:D:")) != -1)
		switch (ch) {
			case 'f':
				forceEvenWithInconsistentGeometry = true;
				break;
			case 'l':
				strncpy(logfile, optarg, 200);
				break;
            case 'u':
                strncpy(buf, optarg, 255);
                sscanf(buf, "%lf", &mu_s);
                break;
            case 'a':
                strncpy(buf, optarg, 255);
                sscanf(buf, "%lf", &mu_a);
                break;
            case 'g':
                strncpy(buf, optarg, 255);
                sscanf(buf, "%lf", &anisotropy);
                break;
            case 'w':
                strncpy(buf, optarg, 255);
                sscanf(buf, "%lf", &wavelength);
                break;
            case 'o':
                outputFilename = string(optarg);
                break;
			case 'D': 
				{
				string temp(optarg);
				vector<string> varDef;
                split(temp, varDef,"=");
				if (varDef.size() != 2)
					clog << "Warning: Improper definition";
				else
					dict[varDef[0]] = varDef[1];
				}
				break;
			case '?':
			default:
				printf("Usage incorrect: polmc -f -l logfile -u mu_s -a mu_a -g g -o outputfile parameterfile\n");
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
        signal(SIGUSR2, catch_signal);


		gStatus = kInitializing;

		
		/* 1) Read parameters necessary to run the simulations.
		 After this, dict contains all the parameters for the simulation. */
		
		CFDictionaryRef propertyList;

        if (argc == 1) {
            clog << RightNow() << "Reading parameters from file " << argv[argc-1] << endl ;

			CFURLRef fileURL;
			CFDataRef propertyListData;

			CFStringRef theStringRef = CFStringCreateWithCString(kCFAllocatorDefault, argv[argc-1],kCFStringEncodingISOLatin1);
			
			fileURL = CFURLCreateWithFileSystemPath( kCFAllocatorDefault, theStringRef,  kCFURLPOSIXPathStyle, false ); 
			CFURLCreateDataAndPropertiesFromResource(kCFAllocatorDefault, fileURL, &propertyListData, NULL, NULL, NULL);
			propertyList = (CFDictionaryRef) CFPropertyListCreateFromXMLData (kCFAllocatorDefault, propertyListData, kCFPropertyListImmutable, NULL);
			
			if ( propertyList == NULL ) {
				cout << "Failed at reading property list\n";
				return -1;
			}
		} else {
            clog << RightNow() << "You must provide a parameter file on the command-line\n" ;
			return -1;
        }

		/* 2) Initialize variables for Monte Carlo calculation 
		 
		 It is done in three major steps:
		 
		 1) calculationParameters: anything that is related to output, statistics, computation in general.
			The object MCWorld contains all the objects
		 2) objects: any object with its properties (type (box, generic, layer, etc...), position, sacttering properties, optical properties, etc...
					 are read here add added to MCWorld
		 3) lightSources: all lightsources with their properties (number of photons, position, type (i.e. laser, isotropic, etc..), wavelength, etc..)
		 */
		CFDictionaryRef calculationDict = (CFDictionaryRef)CFDictionaryGetValue(propertyList, CFSTR("calculationParameters"));
		MCWorld* world = CreateBareWorldFromCFDictionary( (CFDictionaryRef)CFDictionaryGetValue(calculationDict, CFSTR("volumeStats")));
		SetStatisticsForCalculation(world, calculationDict);

		CFArrayRef objects = (CFArrayRef)CFDictionaryGetValue(propertyList, CFSTR("objects"));
		PopulateRootObjectWithObjectsFromCFDictionary(world, world, objects);
		CFArrayRef lightSources = (CFArrayRef)CFDictionaryGetValue(propertyList, CFSTR("lightSources"));
		PopulateWorldWithLightSourcesFromCFDictionary(world, lightSources);
		
		vector<MCSource*> allSources = world->GetSourceList();
		
		gStatus = kPropagatingPhoton;
		
		UTimer timer;
		
		ofstream photonStream;
		if ( world->OutputFilenamePhotonPaths() != string("") ) {
			photonStream.open(world->OutputFilenamePhotonPaths().c_str(), std::ios_base::out);
		}

		for ( unsigned long s = 0; s < allSources.size(); s++ ) {
			MCSource* theSource = allSources[s];
			int l = 1;
			long N_photon = theSource->GetNumberOfPhotonsToLaunch();
			timer.Reset();
			for ( long photon = 1 ; photon <= N_photon; photon++ ) {
				if (photon % l == 0)    {
					clog << RightNow() << "Photon [" << photon << "/" << N_photon
						<< "] in " << timer.DurationInText(timer.GetTicks()) << ", expect end in " << timer.DurationInText(TimeT(timer.GetTicks()*(N_photon/photon-1))) << endl;
					clog.flush();
					l *= 10;
				}

				auto_ptr<Photon> thePhoton(theSource->GetNewPhoton());

				if ( world->KeepPhotonStatistics() ) {
					thePhoton->KeepStats();
				}
				
				/* This is actually where the work gets done */
				world->PropagateInObject(thePhoton.get());
								
				if ( world->OutputFilenamePhotonPaths() != string("")) {
					thePhoton->DumpStats(photonStream, world->OutputFilenamePhotonPathsFormat());
				}

/*				if ( world->DebugOutputPhotonStatistics() ) {
					thePhoton->DumpStats(cout, kRaw);
				}
*/				
				
				if (gTryToQuitNicely) {
					// User has tried to quit
					clog << RightNow() << "User asked to quit.  Quitting nicely." << endl;
					break;
				}

				if (gSaveCurrentStateToDisk) {
					// User asked for a safety-save
					clog << RightNow() << "User asked to save current state." << endl;
					world->NormalizationFactorForInstensityStatistics(theSource->GetNumberOfPhotonsLaunched());
					world->DumpToFile(world->OutputFilename()+string(".save"));
					cerr << "Done.\n";
					gSaveCurrentStateToDisk = false;
				}

				if (gDumpProximityStats) {
					// User asked for a stat dump
					clog << RightNow() << "User asked to dump stats." << endl;
					world->DumpProximityListStatistics(cout);
					cerr << "Done.\n";
					gDumpProximityStats = false;
				}
				
			} // All photons
		}
   		gStatus = kWritingDataToDisk;

        /* 4) Store everything on disk */
		vector<MCSource*> sourceList = world->GetSourceList();
		
		unsigned long totalPhotons = 0;
		for ( unsigned int i = 0; i < sourceList.size(); i++) {
			MCSource* theSource = sourceList[i];
			totalPhotons += theSource->GetNumberOfPhotonsLaunched();
		}
		
		if ( totalPhotons > 0)
			world->NormalizationFactorForInstensityStatistics(totalPhotons);

        world->DumpToFile(world->OutputFilename());

		if ( world->OutputFilenameGeometryWithIntensity() != string("") ) {
			ofstream s(world->OutputFilenameGeometryWithIntensity().c_str(), std::ios_base::out);
			world->DumpGeometryToStream(s, true, world->OutputFilenameGeometryFormat());
		}

		if ( world->OutputFilenameGeometry() != string("") ) {
			ofstream s(world->OutputFilenameGeometry().c_str(), std::ios_base::out);
			world->DumpGeometryToStream(s, true, world->OutputFilenameGeometryFormat());
		}

		

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



		
/* On Unix-based systems, this function allows the program to avoid getting killed
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
				if (gTryToQuitNicely == true) {
					cerr << "User asked again.  Aborting this time. Files might be corrupted." << endl;
					exit(1);
				}
					
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
         case SIGUSR2:
             if (gStatus != kWritingDataToDisk) {
                 gDumpProximityStats = true;
                 cerr << "Trying to dump stats" << endl;
             } else {
                 cerr << "Already saving data" << endl;
             }
             break;
     }
 }
 

