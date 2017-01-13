#include "MCDetector.h"
#include "Photon.h"
#include "MCUtils.h"
#include "MuellerM.h"

MCDetector::MCDetector(double inWidth, long inN, long inDetectionType, long inAcceptanceCosineElements, bool thereIsAPolarizer, int polarizeAxis)
:MCBox(inWidth, inWidth, 0.01, inN, inN, 1, 1, 1, 0, inDetectionType, inAcceptanceCosineElements, false)
{
    mObjectName="detector";
    mBallisticPhotons = false;
	
	for (long i = 0; i < mNumOfInterfaces; i++) {
		mSurfaceElements[i]->TMin = 0e-12;
		mSurfaceElements[i]->TMax = 1000e-12;
		mSurfaceElements[i]->Nt = 100;
		mSurfaceElements[i]->mTimedIntensity = 0;
		mSurfaceElements[i]->init();
	}
	
//		MuellerMSpherical* inM;
//		inM = new MuellerMSpherical(0,0,0,0,0,0);

//	if (thereIsAPolarizer){
//		SetPolarizerAlong(polarizeAxis, inM);
//	}
}


IntersectElement
MCDetector::PropagateInObject(Photon *ioPhoton, double inDist)
{
	RealV pos;
	
	ioPhoton->SetWeight(0);

	IntersectElement ie;
	return ie;
}

void
MCDetector::SetPolarizerAlong(int axis, MuellerM *inM ) //axis: x=1, y=2, z=3
{

	// keep track of oreintation of poalrizer if present 
	double m[4][4]; 
	
	inM->linear_polarizer(&m[0][0], axis);
	
}


void
MCDetector::ScoreOnSurface(Photon *ioPhoton, IntersectElement& inIntersectElement)
{
	
	if (thereIsAPolarizer) {
		StokesV* inS;
		double outIpara, outIperp;
		RealV inEparaLab, inEperpLab, inNormalLab;
		
		switch (polarizeAxis)   
		{
		  case 1:  //polarization allowing x direction pass through
			inEparaLab = RealV(1,0,0);
			inEperpLab = RealV(0,1,0);
			break;
		  case 2:  //polarization allowing y direction pass through
			inEparaLab = RealV(1,0,0);
			inEperpLab = RealV(0,1,0);
			break;
		  default:
			throw logic_error("Warning: Polarizer has polarization direction other than x or y direction");
		}
		
		inNormalLab = RealV(0,0,1);  //set the norm of polarizer to z
		
		ioPhoton->IntensityThroughLinearPolarizer(*inS, inEparaLab, inEperpLab, inNormalLab, outIpara, outIperp);

	}

	if ( ioPhoton->GetNumberOfScatteringEvents() == 0 && ! mBallisticPhotons ) {
		return;
	}
	
	inIntersectElement.surfaceElement->ScoreOnSurface(ioPhoton, inIntersectElement);
}


bool
MCDetector::IsTransmitted(Photon *ioPhoton, IntersectElement& inIntersectElement)
{
	return true;
}
