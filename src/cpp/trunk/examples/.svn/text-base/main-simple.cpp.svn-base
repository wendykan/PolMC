/*
 
This file is a stripped down, simplified version to demonstrate
how to use the computer code in your applications.

See the main.cpp file for extensive comments. 
*/

#include "common.h"

#ifndef MAXPATHLEN
    #define MAXPATHLEN 1000
#endif
   
using namespace std;

int
main(int argc, char *argv[])
{
    try {
		long N_photon = 100000;

#if WITH_POLARIZATION
        clog <<  RightNow() << "Program compiled with the SCATMECH library: full Mueller matrix can be used and polarization is considered.\n";
#else
        clog <<  RightNow() << "Program compiled without the SCATMECH library: only Henyey-Greenstein can be used.\n";
#endif

		double mu_s = 50;
		double mu_a = 1;
		double radiusOfSpheres = 1e-4;
		double indexOfSpheres = 1.59;
		double indexOutside = 1.33;
		double wavelength = 633e-7;
        double anisotropy = 0.9;
        
        MCRandomScatterer *randomScatterer;
        
#if WITH_POLARIZATION
        randomScatterer = new MCRandomScattererJaillon(radiusOfSpheres, indexOfSpheres, indexOutside, wavelength, 10000, mu_s, mu_a, 0);
#else
        randomScatterer = new MCRandomScattererHenyeyGreenstein(anisotropy, mu_s, mu_a);
#endif
 
        MCWorld world;

		RealV origin(0,0,0);
		RealV size(1,1,1);
		double index_outside = 1;
        double index_medium = 1.33;
        double rotationPerCmInClearSpace = 0;
		long Nx = 10, Ny=10, Nz=10;
        long acceptanceCosineElements=10;

        MCObject* sample;
		sample = new MCBox(size.x, size.y, size.z, Nx, Ny, Nz, 
						   index_medium, index_outside, rotationPerCmInClearSpace,
						   kFromInside, acceptanceCosineElements, false );
		sample->SetRandomScatterer(randomScatterer);
		world.PlaceIntoObject(sample, origin);

		
		RealV detectorPosition(0,0,3);
		long detectorPixels = 50;
		double detectorSize = 5;
		
        MCDetector* detector = new MCDetector(detectorSize, detectorPixels, 0, 0);
        world.PlaceIntoObject(detector, detectorPosition);
		
		auto_ptr<MCSource> theSource;

		RealV dir(0,0,1), pos(0,0,-1);
        StokesV inputStokes(1,1,0,0, RealV(1,0,0), dir);

        theSource.reset(new MCLaserSource(inputStokes, 0, 0, pos, dir));
        theSource->SetContainerObject(&world);
        
        if (! world.IsGeometryConsistent() ) {
            throw runtime_error("Geometry inconsistent.");
        }
        
        long l = 10;

        /*Loop through all photons */
        UTimer timer;
		
        for (long photon = 1; photon <= N_photon; photon++) {
            
            if (photon % l == 0)    {
                clog << RightNow() << "Photon [" << photon << "/" << N_photon
                    << "] in " << timer.DurationInText(timer.GetTicks()) << ", expect end in " << timer.DurationInText(TimeT(timer.GetTicks()*(N_photon/photon-1))) << endl;
                clog.flush();
                l *= 10;
            }

            // Deleted automatically after end of scope
            PhotonCote thePhoton;

			theSource->InitializePhoton(&thePhoton);
            world.PropagateInObject(&thePhoton);
            
        } // All photons

        world.NormalizationFactorForInstensityStatistics(theSource->GetNumberOfPhotonsLaunched());
        world.DumpToFile("OutputData");

    } catch (exception& e) {
        clog << RightNow() << "Exception caught: " << string(e.what()) << endl;
        clog << RightNow() << "Program ended abnormally" << endl;
        return 1;
    }

    return 0;
}


