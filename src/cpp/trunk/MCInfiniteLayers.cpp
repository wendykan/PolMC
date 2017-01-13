#include "MCInfiniteLayers.h"

#include "configfiles.h"
#include "MCObject.h"
#include "Photon.h"
#include "StokesV.h"
#include "MuellerM.h"
#include "MCRandomScatterer.h"
#include "MCUtils.h"


MCInfiniteLayer::MCInfiniteLayer(double inThickness,
                                 double inImgWidth,
                                 double inImgHeight,
                                 long inPtsWidth,
                                 long inPtsHeight,
                                 double inIndexMedium,
                                 double inIndexOutside,
                                 double inRotationPerCmInClearSpace,
								 long inDetectionType,
                                 long inAcceptanceCosineElements,
                                 bool inUseNextEventEstimator):MCObject(inIndexMedium, inIndexOutside, inRotationPerCmInClearSpace, inUseNextEventEstimator)
{
     init(inThickness, inImgWidth, inImgHeight, inPtsWidth, inPtsHeight, inDetectionType, inAcceptanceCosineElements);
}


void
MCInfiniteLayer::init(double inThickness, double inImgWidth, double inImgHeight, long inPtsWidth, long inPtsHeight, long inDetectionType, long inAcceptanceCosineElements)
{
    mObjectName="infinite layer";
    mThickness = inThickness;

    mNumOfInterfaces = 2;

    for (long i = 0; i < mNumOfInterfaces; i++) {
        // We have to allocate individually because
        // new[] must be deleted with delete[], and we can't do that
        // because we give ownership to MCWorld (later on)
        SurfaceElement* se = new SurfaceElement;
        mSurfaceElements.push_back(se);
    }
    
    double hw = inImgWidth/2.;
    double hh = inImgHeight/2.;
	double hthickness = inThickness/2.;
	
    mSurfaceElements[kForwardZPlane]->name = "forward";
    mSurfaceElements[kForwardZPlane]->origin = RealV(-hw,-hh,inThickness);
    mSurfaceElements[kForwardZPlane]->a = RealV(1,0,0);
    mSurfaceElements[kForwardZPlane]->b = RealV(0,1,0);
    mSurfaceElements[kForwardZPlane]->a *= inImgWidth;
    mSurfaceElements[kForwardZPlane]->b *= inImgHeight;
    mSurfaceElements[kForwardZPlane]->Na = inPtsWidth;
    mSurfaceElements[kForwardZPlane]->Nb = inPtsHeight;
    mSurfaceElements[kForwardZPlane]->el = RealV(1,0,0);
    mSurfaceElements[kForwardZPlane]->er = RealV(0,1,0);
    mSurfaceElements[kForwardZPlane]->normal = RealV(0,0,1);
    mSurfaceElements[kForwardZPlane]->surfaceShape = kInfinitePlane;
    mSurfaceElements[kForwardZPlane]->mDetection = inDetectionType;
	mSurfaceElements[kForwardZPlane]->mCosineElements = inAcceptanceCosineElements;
	mSurfaceElements[kForwardZPlane]->init();
    
    mSurfaceElements[kBackwardZPlane]->name = "backward";
	// Wendy says: I changed it to -hw to check signs.
    mSurfaceElements[kBackwardZPlane]->origin = RealV(-hw,-hh,0);
    mSurfaceElements[kBackwardZPlane]->a = RealV(-1,0,0);
    mSurfaceElements[kBackwardZPlane]->b = RealV(0,1,0);
    mSurfaceElements[kBackwardZPlane]->a *= inImgWidth;
    mSurfaceElements[kBackwardZPlane]->b *= inImgHeight;
    mSurfaceElements[kBackwardZPlane]->Na = inPtsWidth;
    mSurfaceElements[kBackwardZPlane]->Nb = inPtsHeight;
    mSurfaceElements[kBackwardZPlane]->el = RealV(-1,0,0);
    mSurfaceElements[kBackwardZPlane]->er = RealV(0,1,0);
    mSurfaceElements[kBackwardZPlane]->normal = RealV(0,0,-1);
    mSurfaceElements[kBackwardZPlane]->surfaceShape = kInfinitePlane;
    mSurfaceElements[kBackwardZPlane]->mDetection = inDetectionType;
	mSurfaceElements[kBackwardZPlane]->mCosineElements = inAcceptanceCosineElements;
	mSurfaceElements[kBackwardZPlane]->init();
    
    FinishCreate();

}

void
MCInfiniteLayer::init(map<string,string>& inDict)
{

    long Nx, Ny, acceptanceCosElements;
    double imagingWidth, imagingHeight;
    
    InitVariable(dict, "thickness", mThickness);
    InitVariable(dict, "Nx",Nx);
    InitVariable(dict, "Ny",Ny);
    InitVariable(dict, "imagingWidth",imagingWidth);
    InitVariable(dict, "imagingHeight",imagingHeight);
    InitVariable(dict, "acceptanceCosElements",acceptanceCosElements);

    init(mThickness, imagingWidth, imagingHeight, Nx, Ny, kFromInside, acceptanceCosElements);
}


bool
MCInfiniteLayer::IsOutsideObject(RealV& inLocalPosition)
{
	// We are now defining an actual boundary region - surface element is not in the interior of the object
	// added Wendy 7/15/10
    if (inLocalPosition.z <= (mThickness + OBJECT_SAFETY_DISTANCE) && inLocalPosition.z >= (0.0 - OBJECT_SAFETY_DISTANCE) ) {
        return false;
    } else {
        return true;
    }

}


//added by Wendy Kan
bool
MCInfiniteLayer::IsInsideObject(RealV& inLocalPosition)
{
	// We are now defining an actual boundary region - surface element is not in the interior of the object
	// added Wendy 7/15/10
    if (inLocalPosition.z <= (mThickness - OBJECT_SAFETY_DISTANCE) && inLocalPosition.z >= (0.0 + OBJECT_SAFETY_DISTANCE) ) {
        return true;
    } else {
        return false;
    }
}


void MCInfiniteLayer::InterfaceTouchesOtherObject(MCObject * inOtherObject, long inOtherInterface, long thisObjectInterface){
	clog << RightNow() << "Running: MCInfiniteLayer::InterfaceTouchesOtherObject, we're touching " << this->GetName() << " and " << inOtherObject->GetName() << endl;

    vector<SurfaceElement*> thisObjectSurfaceElements;
    long thisObjectHowMany;
    thisObjectHowMany = GetSurfaceElements(thisObjectSurfaceElements);

    vector<SurfaceElement*> otherObjectSurfaceElements;
    long otherObjectHowMany;
    otherObjectHowMany =inOtherObject->GetSurfaceElements(otherObjectSurfaceElements);

	// thisObjectInterface and otherSurfaceElement will be either kForwardZPlane or kBackwardZPlane
   SurfaceElement* thisSE = thisObjectSurfaceElements[thisObjectInterface];
   SurfaceElement* otherSE = otherObjectSurfaceElements[inOtherInterface];

	// Sanity checks for programmers
	if ( thisSE == NULL || otherSE == NULL  ) {
		throw logic_error("Programming error: surface elements null in MCInfiniteLayer::InterfaceTouchesOtherObject");
	}
	
	if ( thisSE->objectInside != this || otherSE->objectInside != inOtherObject) {
		throw logic_error("Programming error: objects assigned to incorrect objects in MCInfiniteLayer::InterfaceTouchesOtherObject");	
	}
	
   //
   // Check that thisSurfaceElement and otherSurfaceElement

   // 1) The normals are parallel
   // 2) non NULL
   // 3) Location may be close (?)
   // 5) thisSurfaceElement->objectInside == this  && thisSurfaceElement->objectOutside == NULL
   // 6) otherSurfaceElement->objectInside == inOtherObject && otherSurfaceElement->objectOutside == NULL
	bool parallel = RealV::areParallel(thisSE->normal,otherSE->normal);
	bool areClose = areTheSame( (thisSE->origin).z , (otherSE->origin).z, 3 );
	if ( parallel &&  areClose) {

		if ( otherSE->objectOutside != thisSE->objectOutside ) {
			clog << RightNow() << "Warning: Expected outer objects to be the same: otherSE->objectOutside= "<< otherSE->objectOutside << "thisSE->objectOutside=" << thisSE->objectOutside<< endl;
		}
		
		// If fulfilled, then adjust the SE's objectOutside
		otherSE->objectOutside = this;
		thisSE->objectOutside = inOtherObject;
//		mSurfaceElements[thisObjectInterface] = otherSE;
		clog << RightNow() << "We have just touched " << this->GetName() << " and " << inOtherObject->GetName() << endl;
	}
	else{
		// keep both SE
		clog << RightNow() << "Warning: wanting to touch interface, but not successful. parallel = "<< parallel << ",areClose = " << areClose << ", (thisSE->origin).z = " << (thisSE->origin).z << ", (otherSE->origin).z = " << (otherSE->origin).z << endl;
	}

}

