/*
 Pol-MC:  Program to simulate the propagation of light (polarized or not)
 in a three-dimensional turbid medium. 

 Copyright Daniel Cote 2004
 Work done while at Ontario Cancer Institute, University of Toronto
 
 Supported by 
 Natural Science and Engineering Research Council of Canada
 University Health Network
 Canadian Institute for Photonics Innovation

 Main web page where to find the latest version of this code:
 http://www.novajo.ca/ont-canc-inst-biophotonics/
 
 You may use this code for your own research, with no garantee 
 expressed or implied. It has been  tested and partly validated
 (linear depolarization, optical rotation, reflectance, transmittance).
 
 If you do use this program, the appreciated citation is:
 
 Daniel Cote and Alex Vitkin, "Robust concentration determination 
 of optically active molecules in turbid media with validated
 three-dimensional polarization Monte Carlo model"
 Optics Express, Dec 2004
 
 The original references upon which this code is based are:
 
 "Description and time reduction of a Monte Carlo code to simulate
 propagation of polarized light through scattering media", 
 Jaillon et al. , Applied Optics,  42, No 16, p. 3290, (2003)
 
 "Monte Carlo simulations of the diffuse backscattering 
 Mueller matrix for highly scattering media",
 Bartel et al, Appl. Opt. Vol 39, April 1st 2000, p. 1580
 
 "Characterizing mammalian cells and cell phantoms by 
 polarized backscatterng fiber-optic measurements"
 J. Mourant, T. Johnson, J. Freyer
 Appl. Opt. Vol 40, No8 Oct 2001, p. 5114
 
 "Mueller matrix of dense polystyrene latex sphere suspensions:
 measurements and Monte Carlo simulation"
 B. Kaplan, G. Ledanois, B. Drevillon
 Appl. Opt. Vol 40, no 16, June 2001, p. 2769-2777
 
 which are themselves based on code from Wang and Jacques
 
 "MCML-Monte Carlo modeling of light transport" in Optical-Thermal
 response of laser irradiated tissue, A. J. Welch M. van Gemert
 pp. 73-100 http://omlc.ogi.edu/software/mc/index.html
 
 
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include "polmc.h"

#ifndef MAXPATHLEN
#define MAXPATHLEN 1000
#endif

using namespace std;

enum {kPencilBeam = 0, kDiskBeam};

enum {
	kInitializing,
	kPropagatingPhoton,
	kWritingDataToDisk
};

bool gTryToQuitNicely = false;
bool gSaveCurrentStateToDisk = false;
int	gStatus = kInitializing;

/* This is a very low-level Unix function to
interact nicely with user. Defined below. */
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
        
        string formalism;
        string randomSampling;
        string objectType;
        string inputBeamType;
        long formalismID;
        long randomSamplingID;
        double inputBeamDiameter;
        
        vector<StokesV> incidentStokesV;
        string outputFilename;
        
        /*  Read parameters necessary to run the simulations.
            After this, "dict" (which is an associative array, a hash, a dictionary where elements
            are accessed by name) contains all the parameters for the simulation.
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
        
        /* Initialize variables for Monte Carlo calculation  */
		
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
        
        InitVariable(dict, "Io",Io);
        InitVariable(dict, "Qo",Qo);
        InitVariable(dict, "Uo",Uo);
        InitVariable(dict, "Vo",Vo);
        
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
        
#if WITH_POLARIZATION
        clog <<  RightNow() << "Program compiled with the SCATMECH library: full Mueller matrix can be used and polarization is considered.\n";
#else
        if (formalismID != kIntensity) {
            clog <<  RightNow() << "Program compiled without the SCATMECH library: only Henyey-Greenstein random sampling allowed with Intensity formalism.\n";
            throw runtime_error("formalism must be Intensity");
        }
#endif
        

        InitVariable(dict, "randomSampling",randomSampling);
        
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
        
        MCRandomScatterer* randomScatterer;
        double anisotropy, mu_s, mu_a;
        double radius_scat, index_scat, index_med, wavelength;
        long N_scatterpts;
                
#if WITH_POLARIZATION

        switch (randomSamplingID) {
            case kKaplan:
                InitVariable(dict, "mu_s", mu_s);
                InitVariable(dict, "mu_a", mu_a);
                InitVariable(dict, "index_med", index_med);
                InitVariable(dict, "index_scat", index_scat);
                InitVariable(dict, "wavelength", wavelength);
                InitVariable(dict, "radius_scat", radius_scat);
                InitVariable(dict, "N_scatterpts", N_scatterpts);
                    
                randomScatterer = new MCRandomScattererKaplan(radius_scat, index_scat, index_med, wavelength, N_scatterpts, mu_s, mu_a);
                break;
            case kHenyeyGreenstein:
                if (! VariableExists(dict, "g") ) {
                    clog <<  RightNow() << "The variable 'g' (the anisotropy) is not defined.\n";
                    randomScatterer = new MCRandomScattererHenyeyGreenstein(dict);
                } else {
                    
                    InitVariable(dict, "g", anisotropy);
                    InitVariable(dict, "mu_s", mu_s);
                    InitVariable(dict, "mu_a", mu_a);
                    randomScatterer = new MCRandomScattererHenyeyGreenstein(anisotropy, mu_s, mu_a);
                }
                
                break;
            case kJaillon:
                InitVariable(dict, "mu_s", mu_s);
                InitVariable(dict, "mu_a", mu_a);
                InitVariable(dict, "index_med", index_med);
                InitVariable(dict, "index_scat", index_scat);
                InitVariable(dict, "wavelength", wavelength);
                InitVariable(dict, "radius_scat", radius_scat);
                InitVariable(dict, "N_scatterpts", N_scatterpts);
                
                randomScatterer = new MCRandomScattererJaillon(radius_scat, index_scat, index_med, wavelength, N_scatterpts, mu_s, mu_a, 0);
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
            InitVariable(dict, "mu_s", mu_s);
            InitVariable(dict, "mu_a", mu_a);
        }
        
        randomScatterer = new MCRandomScattererHenyeyGreenstein(anisotropy, mu_s, mu_a);
        
#endif
        
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
        
        /* Creation of the object through which photons will propagate                                                                                               
	        All the parameters necessary for constructions are in dict */
        MCWorld world(Vol_XMin, Vol_XMax, Vol_YMin, Vol_YMax, Vol_ZMin, Vol_ZMax, Vol_Nx, Vol_Ny, Vol_Nz);
        
        auto_ptr<MCObject> sample;
        InitVariable(dict, "objectType",objectType);
        
		/* Most objects use these, so they are declared here */
		RealV origin, size;
		double index_outside, index_medium, rotationPerCmInClearSpace;
		long Nx, Ny, Nz, acceptanceCosineElements;
		InitVariable(dict, "origin", origin);
		InitVariable(dict, "index_outside", index_outside);
		InitVariable(dict, "index_med", index_medium);
		
		if (VariableExists(dict, "acceptanceCosElements")) {
			InitVariable(dict, "acceptanceCosElements",acceptanceCosineElements);
		} else {
			acceptanceCosineElements = 1;
		}
		
		if (VariableExists(dict, "rotationPerCmInClearSpace")) {
			clog << RightNow() << "Optical activity is considered\n";
			InitVariable(dict, "rotationPerCmInClearSpace",rotationPerCmInClearSpace);
		} else {
			rotationPerCmInClearSpace = 0.;
		}
		
		
        if (objectType == "box") {
            InitVariable(dict, "size", size);
            InitVariable(dict, "Nx",Nx);
            InitVariable(dict, "Ny",Ny);
            InitVariable(dict, "Nz",Nz);
            sample.reset(new MCBox(size.x, size.y, size.z, Nx, Ny, Nz, 
                                   index_medium, index_outside, rotationPerCmInClearSpace,
                                   kFromInside, acceptanceCosineElements, false ));
            sample->SetRandomScatterer(randomScatterer);
            world.PlaceIntoObject(sample.get(), origin);
        } else if (objectType == "layer") {
            InitVariable(dict, "size", size);
            InitVariable(dict, "Nx",Nx);
            InitVariable(dict, "Ny",Ny);
            sample.reset(new MCInfiniteLayer( size.z, size.x, size.y, Nx, Ny, 
                                              index_medium, index_outside, rotationPerCmInClearSpace,
                                              kFromInside, acceptanceCosineElements, false));
            sample->SetRandomScatterer(randomScatterer);
            world.PlaceIntoObject(sample.get(), origin);
        } else if (objectType == "ellipsoid") {
            InitVariable(dict, "size", size);
            InitVariable(dict, "Nx",Nx);
            InitVariable(dict, "Ny",Ny);
            sample.reset(new MCEllipsoid(size.x, size.y, size.z, Nx, Ny, 
                                         index_medium, index_outside, rotationPerCmInClearSpace,
                                         kFromInside, acceptanceCosineElements, false ));
            sample->SetRandomScatterer(randomScatterer);
            world.PlaceIntoObject(sample.get(), origin);
        } else if (objectType == "cylinder") {
			double radius, height;
			long Nh, Nr;
			InitVariable(dict, "Nh",Nh);
			InitVariable(dict, "Nr",Nr);
			InitVariable(dict, "radius",radius);
			InitVariable(dict, "height",height);
            sample.reset(new MCCylinder(radius, height, Nr, Nh, 
                                        index_medium, index_outside, rotationPerCmInClearSpace,
                                        kFromInside, acceptanceCosineElements, false ));
            sample->SetRandomScatterer(randomScatterer);
            world.PlaceIntoObject(sample.get(), origin);
        } else {
            throw runtime_error("Undefined objext type: "+objectType);
        }
        
		RealV detectorPosition;
		long detectorPixels;
		double detectorSize;
		
		if (VariableExists(dict, "detectorPosition")) {
			InitVariable(dict, "detectorPosition", detectorPosition);
			InitVariable(dict, "detectorSize", detectorSize);
			InitVariable(dict, "detectorPixels", detectorPixels);
			
			MCDetector* detector = new MCDetector(detectorSize, detectorPixels, 0, 0);
			world.PlaceIntoObject(detector, detectorPosition);
            clog << RightNow() << "Detector at " << detectorPosition << "\n";
            
		} else {
			clog << RightNow() << "No external detector defined.\n" ;
		}
		
		/* Creation of photon source */
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
			dir.normalize();
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
            theSource.reset(new MCLaserSource(StokesV(Io,Qo,Uo,Vo, RealV(1,0,0), dir), inputBeamDiameter/2., 0, pos, dir));
        } else {
			throw runtime_error("You must define the type of light source you want to use");
        }
        
		if (VariableExists(dict, "sourceInsideSimpleObject"))
			theSource->SetContainerObject(sample.get());
        else
            theSource->SetContainerObject(&world);
        
        if (! world.IsGeometryConsistent() ) {
            throw runtime_error("Geometry inconsistent.");
        }
        
        /* The main loop, going through all Stokes vectors
            and all photons */
        
		gStatus = kPropagatingPhoton;

        clog << RightNow() << "Detection formalism demanded: " << formalism << endl;
        
        long l = 10;
        
        /* Loop through all photons */
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
            
            
            
            /*  This call below is the guts of the calculation.
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
                clog << RightNow() << "User asked to save current state.  Saving to "<< outputFilename+string(".save") << endl;
                world.NormalizationFactorForInstensityStatistics(theSource->GetNumberOfPhotonsLaunched());
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
        
        /* Store everything on disk */
        world.NormalizationFactorForInstensityStatistics(theSource->GetNumberOfPhotonsLaunched());
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
to ask to program to save its current state.  In Unix, use kill -USR1 processID. */

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


