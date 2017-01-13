#if HAVE_CONFIG_H
#include <config.h>
#endif

#include "configfiles.h"
#include "MCObject.h"
#include "MCWorld.h"
#include "Photon.h"
#include "StokesV.h"
#include "MCRandomScatterer.h"
#include "MCUtils.h"
#include "XMLUtil.h"
#include "MCKDTree.h"
#include <list>

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>




long gDebugInterfaceID = kNoInterface;

long MCObject::gObjectCounter = 1;
double MCObject::gNormalizationIntensity = 1;
long MCObject::gVertexOffset = 0;
long MCObject::gNormalOffset = 0;
long MCObject::gProximitySubdivisions = 7;
long MCObject::gMinPhotonCount = 7;
long MCObject::gMinSurfaceElements = 7;
int MCObject::gInterfaceCrossingAlgorithm = kAllAlgorithm;

MCObject::MCObject()
{
    mTree = NULL;
    
    mGlobalOrigin = RealV(0,0,0);
    mIndexOutside = 1;
    rotationPerCmInClearSpace = 0;
    mNumOfInterfaces = 0;
    mWorld = 0;
    cacheIsValid = false;
    
    mEnergy = NULL;
    mFluence = NULL;
    Vol_Nx = 0;
    Vol_Ny = 0;
    Vol_Nz = 0;

    Vol_XMin = 0.;
    Vol_XMax = 0.;
    Vol_YMin = 0.;
    Vol_YMax = 0.;
    Vol_ZMin = 0.;
    Vol_ZMax = 0.;

    mRandomScatterer = new MCRandomScatterer( 0, 0, 1 );
    useNextEventEstimator = false;
    
    N_photon = 0;
    gNormalizationIntensity = 1;
    mAbsorbedPhotons = 0;
    
	mBallisticPhotons = true;   
	
    mObjectID = gObjectCounter++;

}

MCObject::MCObject(double inIndexMedium, double inIndexOutside, double inRotationPerCmInClearSpace, bool inUseNextEventEstimator)
{
    mObjectID = gObjectCounter++;
    mNumOfInterfaces = 0;

    init(inIndexMedium, inIndexOutside, inRotationPerCmInClearSpace, inUseNextEventEstimator);

}

MCObject::MCObject(map<string,string>& inDict)
{
    mObjectID = gObjectCounter++;
    mNumOfInterfaces = 0;

    init(inDict);

}

void
MCObject::FinishCreate()
{
    /* This function must be called by subclasses when they are done creating themselves
       to ensure MCObject has a chance to finish creating various tables */
    
    if (mVertexList.size() == 0) {
        for (unsigned long i = 0; i < mSurfaceElements.size(); i++) {
            mSurfaceElements[i]->init();
            mSurfaceElements[i]->objectInside = this;
        }
        mNumOfInterfaces = mSurfaceElements.size();
        
        if (InitConnectivityTablesFromSurfaceElements()) 
            clog << RightNow() << "Warning: failed initializing the connectivity table for object "+mObjectName << endl;
    } else {
        if (InitSurfaceElementsFromConnectivityTables()) 
            clog << RightNow() << "Warning: failed initializing the surface elements from connectivity table for object "+mObjectName << endl;

        for (unsigned long i = 0; i < mSurfaceElements.size(); i++) {
            mSurfaceElements[i]->init();
        }
    }
    
}

void
MCObject::init(double inIndexMedium, double inIndexOutside, double inRotationPerCmInClearSpace, bool inUseNextEventEstimator)
{
    mGlobalOrigin = RealV(0,0,0);
    mIndexOutside  = inIndexOutside;
    rotationPerCmInClearSpace = inRotationPerCmInClearSpace;
    useNextEventEstimator = inUseNextEventEstimator;
    mWorld = 0;

	mRandomScatterer = new MCRandomScatterer( 0, 0, 1 );
    useNextEventEstimator = false;

	mBallisticPhotons = true;   

    N_photon = 0;
    gNormalizationIntensity = 1;
    mAbsorbedPhotons = 0;
    cacheIsValid = false;

    mEnergy = NULL;
    mFluence = NULL;
    Vol_Nx = 0;
    Vol_Ny = 0;
    Vol_Nz = 0;

    Vol_XMin = 0.;
    Vol_XMax = 0.;
    Vol_YMin = 0.;
    Vol_YMax = 0.;
    Vol_ZMin = 0.;
    Vol_ZMax = 0.;
    
}


void
MCObject::init(map<string,string>& inDict)
{
    dict = inDict;

    double theIndexMedium;
    double theIndexOutside;
    double theRotationPerCmInClearSpace;
    bool theUseNextEventEstimator;
    
    InitVariable(dict, "index_outside",theIndexOutside);
    InitVariable(dict, "index_med",theIndexMedium);

    if (VariableExists(dict, "rotationPerCmInClearSpace")) {
        clog << RightNow() << "Optical activity is considered\n";
        InitVariable(dict, "rotationPerCmInClearSpace",theRotationPerCmInClearSpace);
    } else {
        theRotationPerCmInClearSpace = 0.;
    }

    if (VariableExists(dict, "Vol_Nx")) {
        InitVariable(dict, "Vol_Nx",Vol_Nx);
        InitVariable(dict, "Vol_Ny",Vol_Ny);
        InitVariable(dict, "Vol_Nz",Vol_Nz);
        InitVariable(dict, "Vol_XMin",Vol_XMin);
        InitVariable(dict, "Vol_YMin",Vol_YMin);
        InitVariable(dict, "Vol_ZMin",Vol_ZMin);
        InitVariable(dict, "Vol_XMax",Vol_XMax);
        InitVariable(dict, "Vol_YMax",Vol_YMax);
        InitVariable(dict, "Vol_ZMax",Vol_ZMax);
    }

    if (VariableExists(dict, "useNextEventEstimator")) {
        InitVariable(dict, "useNextEventEstimator", theUseNextEventEstimator);
        if (useNextEventEstimator) {
            clog << RightNow() << "Using next event estimator"<< endl;
        }
    } else {
        theUseNextEventEstimator = false;
    }

    init(theIndexMedium, theIndexOutside, theRotationPerCmInClearSpace, theUseNextEventEstimator);
    
}

bool
MCObject::InitConnectivityTablesFromSurfaceElements()
{
    /*
     We make two lists: vertices and normals.  The vertex one
     contains all vertices (duplicates are removed), whereas normals
     contains the normals.  Each surfaceElement will then have 
     a list of indicess to vertices that make up its shape as well as an index to
     the normal.
     
     We have to worry about the shape of elements: triangles, parallelogram
     */
    mVertexList.clear();
    mNormalList.clear();
    mFaceList.clear();
    
    for (unsigned long i = 0; i < mSurfaceElements.size() ; ++i) {
        SurfaceElement* se = mSurfaceElements[i];
        Face theFace;
        Vertex theVertex;

        /* Avoid normal duplicates */
        bool dup = false;
        unsigned long k = 0;
        for (k = 0; k < mNormalList.size(); k++) {
            if (se->normal == mNormalList[k]) {
                dup = true;
                break;
            }
        }
        if ( ! dup) {
            mNormalList.push_back(se->normal);
        }
        theVertex.normal = k;

        /* Deal with vertices */
        long Na, Nb;
        
        if (se->surfaceShape == kTriangle) {
            Na = 1;
            Nb = 1;
        } else {
            Na = se->Na;
            Nb = se->Nb;
        }
        
        for (long u = 0; u < Na; u++) {
            
            for (long v = 0; v < Nb; v++) {
                long current = 0;
                bool done = false;

                theFace.clear();

                while (! done) {
                    RealV point;
                    theVertex.vertex = 0;
                    theVertex.texture = 0;
                    theVertex.intensity = 0;
            
                    double du = 1./double(Na);
                    double dv = 1./double(Nb);
                    
                    /* Get vertex in local coordinates: when outputting, add mGlobalOrigin
                       that way, if we move object, we don't need to redo the tables */
                    if (se->surfaceShape == kTriangle) {
                        switch (current) {
                            case 0:
                                point = se->origin + se->a * u * du + se->b * v * dv;
                                break;
                            case 1:
                                point = se->origin + se->a * (u+1) * du + se->b * v * dv;
                                break;
                            case 2:
                                point = se->origin + se->a  * u * du + se->b * (v+1) * dv;
                                done = true;
                                break;
                        }
                    } else if (se->surfaceShape == kParallelogram) {
                        switch (current) {
                            case 0:
                                point = se->origin + se->a * u * du + se->b * v * dv;
                                break;
                            case 1:
                                point = se->origin + se->a * (u+1) * du + se->b * v * dv;
                                break;
                            case 2:
                                point = se->origin + se->a  * (u+1) * du + se->b * (v+1) * dv;
                                break;
                            case 3:
                                point = se->origin + se->a  * u * du + se->b * (v+1) * dv;
                                done = true;
                                break;
                        }
                    } else if (se->surfaceShape == kInfinitePlane) {
                        switch (current) {
                            case 0:
                                point = se->origin + se->a * u * du + se->b * v * dv;
                                break;
                            case 1:
                                point = se->origin + se->a * (u+1) * du + se->b * v * dv;
                                break;
                            case 2:
                                point = se->origin + se->a  * (u+1) * du + se->b * (v+1) * dv;
                                break;
                            case 3:
                                point = se->origin + se->a  * u * du + se->b * (v+1) * dv;
                                done = true;
                                break;
                        }
                    } else {
                        clog << RightNow() << "Shapes other than triangles (with 1 sub element) and parallelograms (infinite plane are cut) are not supported yet\n";
                        return true;
                    }
                    
                    /* Avoid vertex duplicates */
                    bool dup = false;
                    unsigned long j = 0;
                    for (j = 0; j < mVertexList.size(); j++) {
                        if (point == mVertexList[j]) {
                            dup = true;
                            break;
                        }
                    }
                    if ( ! dup) {
                        mVertexList.push_back(point);
                    }

                    MyAssert_(j < mVertexList.size() );

                    theVertex.vertex = j;
                    theVertex.texture = 0;
                    if (se->mDetection & kStokesVector) {
                        // Stokes always use arrays
                        if (Na > 1 && Nb > 1)
                            for(long k = 0; k < se->mCosineElements; k++)
                                theVertex.intensity += se->mStokesV[u][v][k].mI;
                    } else {
                        if (Na > 1 && Nb > 1)
                            theVertex.intensity += se->mIntensityArray[u][v];
                        else {
                            theVertex.intensity += se->mIntensity;
						}
                    }
                    // Add vertex to face (copy)
                    theFace.push_back(theVertex);
                    current++;
                }
                // Add face to object
                MyAssert_(theFace.size() == 3 || theFace.size() == 4);
                mFaceList.push_back(theFace);
                current = 0;
                
            }
        }
    }
    
    return false;
    
}

bool
MCObject::InitSurfaceElementsFromConnectivityTables()
{
	throw runtime_error("Not expecting this call");
    /* we have mVertexList, and mFaceList and we rebuild mNormalList */
    mNormalList.clear();
    
    for (unsigned long i = 0; i < mFaceList.size(); ++i) {
        Face theFace;
        
        theFace = mFaceList[i];

        if (theFace.size() != 3) {
            throw runtime_error("Can only work with triangulated shapes for now.");
        }
        
        SurfaceElement* se = new SurfaceElement();
        
        se->origin = mVertexList[theFace[0].vertex-1];
        se->a = mVertexList[theFace[1].vertex-1] - se->origin;
        se->b = mVertexList[theFace[2].vertex-1] - se->origin;
        RealV normal = RealV::CrossProduct(se->a, se->b);
		se->area = normal.abs();
        normal.normalize();
        mNormalList.push_back(normal);
        se->normal = normal;
        se->Na = 1;
        se->Nb = 1;
//        se->indexIn = 1;
//        se->indexOut = 1;
        se->surfaceShape = kTriangle;
        se->mDetection = kIntensityOnlyFromInside;
        se->mCosineElements = 1;
        se->objectInside = this;
        se->mNormalIndex = mNormalList.size() -1;
        se->mVertexIndices.push_back(theFace[0].vertex-1);
        se->mVertexIndices.push_back(theFace[1].vertex-1);
        se->mVertexIndices.push_back(theFace[2].vertex-1);        

        mSurfaceElements.push_back(se);
    }
    
    return false;
}

MCObject::~MCObject()
{
    delete mRandomScatterer;
	delete mTree;
	
	// We don't own the SurfaceElements.  MCWorld does and will delete them.
}


IntersectElement
MCObject::PropagateInObject(Photon *ioPhoton, double inDist)
{
    RealV src,dest;

    double d, theta, phi, mu_t;
    long interface;
    IntersectElement ie;

    PrintMessageIfLevel_("Entering propagation in object " <<GetName() <<"("<<mObjectID<<")"<< " at position " << ioPhoton->GetGlobalPosition(), kVerbose);
	
    N_photon++;
    
    while ( ioPhoton->IsNotDead() ) {
		if (ioPhoton->GetCurrentObject() == this ) {            
//			PrintMessageIfLevel_("ioPhoton->GetCurrentObject = " << *ioPhoton->GetCurrentObject(), kVerbose);

			if (ie.opticalPathLeftover == 0) {
				d = GetRandomScatteringDistance(ioPhoton);
			} else {
                d = (mu_t = GetTotalExtinctionCoefficient(ioPhoton)) != 0 ? ie.opticalPathLeftover / mu_t : INFINITE_DISTANCE ;
			} 
			
 			ie = IntersectElement();
			
			interface = PathCrossesInterface(ioPhoton, d, ie, kAnyEvent);

			if ( interface > kNoInterface) {
                // The path intersects an interface
                
                // Move to interface, then deal with possible reflections
                // Important: if we don't cross (determined later), we really don't
                // want to have the photon on the other side.  Hence we move in two
                // steps: almost there, and if we cross, we move the rest
                MoveWithinObject(ioPhoton, ie.distanceToMove - PHOTON_SAFETY_DISTANCE);
                PrintMessageIfLevel_("Moved to interface at " << ioPhoton->GetGlobalPosition(), kVerbose);
                PrintMessageIfLevel_("ie.objectInside = " << *ie.objectInside, kVerbose);

                if (ie.objectInside->IsTransmitted(ioPhoton, ie)) {
                    // Photon has been transmitted
                    MoveWithinObject(ioPhoton, 2 * PHOTON_SAFETY_DISTANCE);
                    
                    ie.objectInside->TransmitThroughInterface(ioPhoton, ie);
                    
                    ioPhoton->SetCurrentObject(ie.objectOutside);
                    
					
					// Wendy added this
					PrintMessageIfLevel_("ie.objectOutside = " << ie.objectOutside->GetName()<<"("<<ie.objectOutside->mObjectID<<")", kVerbose);
					PrintMessageIfLevel_("ie.objectInside = " << *ie.objectInside, kVerbose);
					PrintMessageIfLevel_("mSurfaceElements.size = " << mSurfaceElements.size(), kVerbose);
                    return ie;
                
                } else {
                    // Photon has been reflected, returning to the top of this loop will take care of proagating the leftover distance
                    PrintMessageIfLevel_(" Photon is being reflected at interface " << ie.interface << " of object " << ie.objectInside->mObjectID << " at position " << ioPhoton->GetGlobalPosition() << " in direction " << ioPhoton->GetPropagationDirectionInLabFrame() , kVerbose);
                    ie.objectInside->ReflectAtInterface(ioPhoton, ie);
                }

                // At this point, we have reflected/crossed/entered/exited anything in our way, and must keep
                // propagating in the current object. This is accomplished through IntersectElement se:
                // here it still contains the leftover optical path, which will get picked up at the beginning of this loop.
                // Only when we get to the interaction point do we actually decrease the weight.
                
			 } else if (interface == kNoInterface)  {
                 // The photon propagates unimpeded to the interaction point
				 PrintMessageIfLevel_(" Moved  by " << ie.distanceToMove << " through object " << mObjectID << " to position " << ioPhoton->GetGlobalPosition(), kExtremelyVerbose);
				 
				 MoveWithinObject(ioPhoton, ie.distanceToMove );
				 DecreaseWeight(ioPhoton);
				 
				 GetRandomScatteringAngles(ioPhoton, theta, phi);
				 ScatterInObject(ioPhoton, theta, phi);
				 // IntersectElement has opticalPathLeftOver == 0
			 } else if ( interface == kInfinity ) {
				 PrintMessageIfLevel_(" Photon moving to infinity in object " << mObjectID , kVerbose);
				 MoveWithinObject(ioPhoton, INFINITE_DISTANCE );
				 IntersectElement ie;
				 ie.interface = kInfinity;
				 return ie;
			 }
			
			/* This is in case something goes very wrong */
            if (IsOutsideBoundingBox(ioPhoton)) {
				
				if ( IsInsideObject(ioPhoton) ) {
					throw logic_error("Code inconsistent: Outside bounding box but inside object.  This code sucks.");
				}
				
				clog << RightNow() << "Emergency exit: object " << mObjectName << "(" << mObjectID << ") after " << ioPhoton->GetNumberOfScatteringEvents() << " scattering events. Position=" << ioPhoton->GetGlobalPosition() << ", direction= " << ioPhoton->GetPropagationDirectionInLabFrame() << endl;

				if (gDebugLevel >= kVerbose) 
					ioPhoton->DumpStats(clog, kRaw);
                N_photon--;
				ioPhoton->SetWeight(0);
                return IntersectElement();
			}
			
		} 
    }

     return ie;
}

IntersectElement
MCObject::RayTraceInObject(Photon *ioPhoton, double inDist)
{
    RealV src,dest;

    double d, theta, phi;
    long interface;
    IntersectElement ie;
    
    while ( ioPhoton->IsNotDead() ) {
		if (ioPhoton->GetCurrentObject() == this ) {            
			d = INFINITE_DISTANCE ;

 			ie = IntersectElement();
			
			interface = PathCrossesInterface(ioPhoton, d, ie, kAnyEvent);

			if ( interface > kNoInterface) {
                // The path intersects an interface
				throw logic_error("the statements following are unreachable code... FIX THIS! ...");
                return ie;
				
                // Move to interface, then deal with possible reflections
                // Important: if we don't cross (determined later), we really don't
                // want to have the photon on the other side.  Hence we move in two
                // steps: almost there, and if we cross, we move the rest
                MoveWithinObject(ioPhoton, ie.distanceToMove - PHOTON_SAFETY_DISTANCE);
                PrintMessageIfLevel_("Moved to interface at " << ioPhoton->GetGlobalPosition(), kVerbose);
                
                if (ie.objectInside->IsTransmitted(ioPhoton, ie)) {
                    // Photon has been transmitted
                    MoveWithinObject(ioPhoton, 2 * PHOTON_SAFETY_DISTANCE);
                    
                    ie.objectInside->TransmitThroughInterface(ioPhoton, ie);
                    
                    ioPhoton->SetCurrentObject(ie.objectOutside);
					
                    
					return ie;
                
                } else {
                    // Photon has been reflected, returning to the top of this loop will take care of proagating the leftover distance
                    PrintMessageIfLevel_(" Photon is being reflected at interface " << ie.interface << " of object " << ie.objectInside->mObjectID << " at position " << ioPhoton->GetGlobalPosition() << " in direction " << ioPhoton->GetPropagationDirectionInLabFrame() , kVerbose);
                    ie.objectInside->ReflectAtInterface(ioPhoton, ie);
                }

                // At this point, we have reflected/crossed/entered/exited anything in our way, and must keep
                // propagating in the current object. This is accomplished through IntersectElement se:
                // here it still contains the leftover optical path, which will get picked up at the beginning of this loop.
                // Only when we get to the interaction point do we actually decrease the weight.
                
			 } else if (interface == kNoInterface)  {
                 // The photon propagates unimpeded to the interaction point
				 MoveWithinObject(ioPhoton, ie.distanceToMove );
				 DecreaseWeight(ioPhoton);
				 
				 GetRandomScatteringAngles(ioPhoton, theta, phi);
				 ScatterInObject(ioPhoton, theta, phi);
				 // IntersectElement has opticalPathLeftOver == 0
			 } else if ( interface == kInfinity ) {
				 PrintMessageIfLevel_(" Photon moving to infinity in object " << mObjectID , kVerbose);
				 MoveWithinObject(ioPhoton, INFINITE_DISTANCE );
				 IntersectElement ie;
				 ie.interface = kInfinity;
				 return ie;
			 }
			
			/* This is in case something goes very wrong */
            if (IsOutsideBoundingBox(ioPhoton)) {
				clog << RightNow() << "Emergency exit: object " << mObjectName << "(" << mObjectID << ") after " << ioPhoton->GetNumberOfScatteringEvents() << " scattering events. Position=" << ioPhoton->GetGlobalPosition() << ", direction= " << ioPhoton->GetPropagationDirectionInLabFrame() << endl;
				if (gDebugLevel >= kVerbose) 
					ioPhoton->DumpStats(clog, kRaw);
                ioPhoton->SetWeight(0);
                return IntersectElement();
			}
			
		} 
    }

     return ie;
}

long
MCObject::PathCrossesInterface(Photon *ioPhoton, double inDist, IntersectElement& outIntersectElement, long inCrossingEvents)
{
    /* This function is the most critical in the algorithm: most of the time is spent here.
    Any time-saving trick that can be used should, as long as the function remains readable. */
    
    CheckDoubleValue_(inDist);

    /* The strategy is to compute the distance to all interfaces that may be crossed, and then keep the closest one. */

	/* First, we check to see if we cross into any objects contained in this object, then we look for crossing out of the current object. */
    //IntersectElement intersect;

    /* Second, we look for crossing of interfaces out of the current object. */

    RealV dir  = ioPhoton->GetPropagationDirectionInLabFrame();
    RealV d    = dir * inDist;
    RealV src  = ioPhoton->GetLocalPosition();

    vector<IntersectElement> possibleIntersection;
    vector<SurfaceElement*> surfaceElements;
    
#if __MYDEBUG

//	IntersectElement debugIntersect;
    IntersectElement debugOutIntersect;
    
    for ( long debug = 0; debug <= 1; debug++) {
        // We are going to do the calculation twice, with the short list, and long list
        //  If we don't get the same answer, there is a problem
		gDebugInterfaceID = kNoInterface;
        if (debug == 0) {
            GetSurfaceElements(surfaceElements);
			possibleIntersection.clear();
        } else {
            // Intersect and outIntersect must be reinitialized
            // debugIntersect = intersect;
            debugOutIntersect = outIntersectElement;
            gDebugInterfaceID = outIntersectElement.interface;
			
            // intersect = IntersectElement();
            outIntersectElement = IntersectElement();
            surfaceElements.clear();
            GetSurfaceElementsCloseToSegment(surfaceElements, src, src+d);
			      
			possibleIntersection.clear();
        }
#else
		//GetSurfaceElements(surfaceElements);
        GetSurfaceElementsCloseToSegment(surfaceElements, src, src+d);
#endif  
      
        PrintMessageIfLevel_("Elements in surfaceElements " << surfaceElements.size() << "/" << mSurfaceElements.size(), kExtremelyVerbose);

   		size_t size = surfaceElements.size();
        for ( unsigned long i = 0; i < size; i++) {
            SurfaceElement* se = surfaceElements[i];
            // http://www.thepolygoners.com/tutorials/lineplane/lineplane.html
            
			IntersectElement inter;
			if ( se->PathMayCrossSurfaceElement(src, d, inter) ) {
				inter.interface = i;
				possibleIntersection.push_back(inter);
			}
        }
			
//		PrintMessageIfLevel_("possibleIntersection.size() = " << possibleIntersection.size(), kVerbose);

        /* We have a list of possible surface elements that could be intersected.
		   Currently, the vecvtor d intersects the infinite plane in which the triangle lies
		   We need to 1) Find the closest one  and then 2) check that the intersect point is on the triangle. */
        sort(possibleIntersection.begin(), possibleIntersection.end());

        
        vector<IntersectElement>::iterator p = possibleIntersection.begin();

        while (p != possibleIntersection.end()) {
			if ( p->surfaceElement->IsCrossing(ioPhoton, src, d, &(*p)) ) {
				outIntersectElement = *p;
				break;
			} 
			 p++;
        }
      
        if (p == possibleIntersection.end() ) {
            outIntersectElement.interface = kNoInterface;
            outIntersectElement.distanceToMove = inDist;
            outIntersectElement.distanceLeftover = 0;
            outIntersectElement.opticalPathLeftover = 0;
            outIntersectElement.surfaceElement = 0;
        }
	
#ifdef __MYDEBUG
    
    }
    
    if (outIntersectElement.surfaceElement != debugOutIntersect.surfaceElement) {
        // There is a problem for sure
		clog << "Fancy algo: " << outIntersectElement.surfaceElement << " Debug : " << debugOutIntersect.surfaceElement << "\n";
		clog << "outIntersectElement uniqueID: " << outIntersectElement.interface << "\ndebug uniqueID: " << debugOutIntersect.interface << "\n";
        throw logic_error("Problem, intersect not same");         
    }
#endif

/*	if ( outIntersectElement.interface != kNoInterface ) {
		MyAssert_(outIntersectElement.indexFrom >= 0);
		MyAssert_(outIntersectElement.indexTo >= 0);
	}
    MyAssert_(outIntersectElement.distanceToMove >= 0);
  */  
    return outIntersectElement.interface;
}

void 
MCObject::SetProximityProperties(long inSubdivisions, long inMinSurfaceElements, long inMinPhotonCount)
{
    gProximitySubdivisions = inSubdivisions;
	gMinSurfaceElements = inMinSurfaceElements;
	gMinPhotonCount = inMinPhotonCount;
}



void 
MCObject::SetOuterObjectTo(MCObject* inObject)
{
    vector<SurfaceElement*> boundarySurfaceElements;
    long howMany = GetBoundarySurfaceElements(boundarySurfaceElements);
    
	if (howMany != 0) {
		for ( long i = 0; i < howMany; i++) {
            boundarySurfaceElements[i]->objectOutside = inObject;
//            boundarySurfaceElements[i]->indexOut = inObject->GetIndexMedium(NULL);
        }
    } else {
        clog << RightNow() << "Warning: object has no surface elements\n";
    }
    
}

MCWorld*
MCObject::GetWorld()
{
    return mWorld;           
}

void
MCObject::SetWorld(MCWorld* inWorld)
{
    mWorld = inWorld;           
}

void
MCObject::SetGlobalOrigin(RealV inOrigin)
{
    mGlobalOrigin = inOrigin;
}
 
RealV
MCObject::GetGlobalOrigin()
{
    return mGlobalOrigin; 
}

bool
MCObject::IsOutsideBoundingBox(Photon *ioPhoton)
{
    RealV top, bottom;
    
    GetBoundingBox(top, bottom);
    
    RealV pt = ioPhoton->GetGlobalPosition();
    
    if (pt.x < top.x && pt.y < top.y && pt.z < top.z &&
        pt.x > bottom.x && pt.y > bottom.y && pt.z > bottom.z ) {
        return false;
    }
    
    return true;
    
}

bool
MCObject::IsOutsideObject(RealV& inLocalPosition)
{
    throw runtime_error("IsOutsideObject(RealV) not implemented in MCObject");
}
    
bool
MCObject::IsOutsideObject(Photon *ioPhoton)
{
	MyAssert_(ioPhoton != NULL);
	bool isOut = false;
	
    // This is mainly a debugging function.
    // Other objects must implement the specifics.
    RealV pos = ioPhoton->GetLocalPosition();
    
	pos.RotateAroundX(-mXRotation);
	pos.RotateAroundY(-mYRotation);
	pos.RotateAroundZ(-mZRotation);

    isOut = IsOutsideObject(pos);

	pos.RotateAroundX(mXRotation);
	pos.RotateAroundY(mYRotation);
	pos.RotateAroundZ(mZRotation);
	
	return isOut;
}

size_t
MCObject::GetNumberOfIntersectionsWithSurface(Photon *ioPhoton)
{
	size_t intersections = 0;
	long interface;
	IntersectElement ie = IntersectElement();
	
	interface = PathCrossesInterface(ioPhoton, INFINITE_DISTANCE, ie, kAnyEvent);
	
	while(interface > kNoInterface) {
		intersections++;
		//we want to make sure we definitely move the photon to the other side of
		//the surface, because we are merely looking for the number of intersections that
		//this ray has with the surface
		MoveWithinObject(ioPhoton, ie.distanceToMove + PHOTON_SAFETY_DISTANCE);
		
		//now continue searching:
		interface = PathCrossesInterface(ioPhoton, INFINITE_DISTANCE, ie, kAnyEvent);
	}
	
	return intersections;
}


bool
MCObject::IsInsideObject(Photon *ioPhoton)
{
	//We will cast a ray in 6 different directions, and check the number
	//of intersections each ray has with the surface of the object
	RealV propagationDirections[6] = {RealV(0,0,1), RealV(0,0,-1), RealV(0,1,0), 
									  RealV(0,-1,0), RealV(1,0,0), RealV(-1,0,0)};
					
	RealV initialPosition = ioPhoton->GetLocalPosition();
	size_t odds = 0;
	size_t evens = 0;				  
	
	for(unsigned int i = 0; i < 6; i++) {
		ioPhoton->SetPropagationDirectionInLabFrame(propagationDirections[i]);
		ioPhoton->SetLocalPosition(initialPosition);
		size_t intersections = GetNumberOfIntersectionsWithSurface(ioPhoton);
		(intersections % 2 == 0) ? evens++ : odds++;
	} 
	
	if(evens > 4) {
		//at least 5 rays have an even number of intersections, photon is *not*in the object
		return false;
	}
	else if(odds > 4) {	
		//at least 5 rays have an odd number of intersections, photon *is* in the object
		return true;
	}
	else { //we do not have a clear answer...
		std::string errorOutput = "Problem, cannot determine whether the photon is in object: ";
		errorOutput.append((ioPhoton->GetCurrentObject())->GetName());
		throw logic_error(errorOutput);    
	} 
}

bool
MCObject::IsInsideObject(RealV& inPoint)
{
	//We will cast a ray in 6 different directions, and check the number
	//of intersections each ray has with the surface of the object
	RealV propagationDirections[6] = {RealV(0,0,1), RealV(0,0,-1), RealV(0,1,0), 
	RealV(0,-1,0), RealV(1,0,0), RealV(-1,0,0)};
	
	// For now, we take advantage of all the functions that Photon and MCObject have
	// to determine the intersection, so we fake a Photon.
	PhotonIntensity aPhoton;
	
	RealV initialPosition = inPoint;
	size_t odds = 0;
	size_t evens = 0;				  
	
	for(unsigned int i = 0; i < 6; i++) {
		aPhoton.SetPropagationDirectionInLabFrame(propagationDirections[i]);
		aPhoton.SetLocalPosition(initialPosition);
		size_t intersections = GetNumberOfIntersectionsWithSurface(&aPhoton);
		(intersections % 2 == 0) ? evens++ : odds++;
	} 
	
	if(evens > 4) {
		//at least 5 rays have an even number of intersections, photon is *not*in the object
		return false;
	}
	else if(odds > 4) {	
		//at least 5 rays have an odd number of intersections, photon *is* in the object
		return true;
	}
	else { //we do not have a clear answer...
		std::string errorOutput = "Problem, cannot determine whether the point is in object: " + GetName();
		throw logic_error(errorOutput);    
	} 
}


bool
MCObject::IsTransmitted(Photon *ioPhoton, IntersectElement& inIntersectElement)
{
    double prob;

    if (SetFresnelCoefficients( inIntersectElement)) {
        return false;
    }
    
    prob = ioPhoton->GetReflectionProbability(inIntersectElement.normal,
                                              inIntersectElement.Rp,
                                              inIntersectElement.Rs,
                                              inIntersectElement.Tp,
                                              inIntersectElement.Ts);
    MyAssert_(prob <= 1.);
    MyAssert_(mRandomScatterer != NULL);

    double num = MCRandomScatterer::RandomFloat();
    MyAssert_(num <= 1.);

    if (num > prob) {
        PrintMessageIfLevel_("Photon transmitted (prob. of reflection was " <<  prob << ")",kVerbose);
        return true;
    } else {
        PrintMessageIfLevel_("Photon reflected (prob. of reflection was " <<  prob << ")",kVerbose);
        return false;
    }
}


bool
MCObject::SetFresnelCoefficients(  IntersectElement& ioIntersectElement) {

    /* Compute angles */
    double cosThetaFrom = abs(ioIntersectElement.cosine);
    CheckDoubleValue_(cosThetaFrom);
    MyAssert_(cosThetaFrom <= 1);

    double sinThetaFrom = sqrt(1.-cosThetaFrom*cosThetaFrom);
    CheckDoubleValue_(sinThetaFrom);

    double m = ioIntersectElement.indexTo/ioIntersectElement.indexFrom;

    CheckDoubleValue_(ioIntersectElement.indexTo);
    CheckDoubleValue_(ioIntersectElement.indexFrom);
	MyAssert_(ioIntersectElement.indexTo >= 1);
	MyAssert_(ioIntersectElement.indexFrom >= 1);

	if ( m == 1.) {
        ioIntersectElement.Rp = 0.;
        ioIntersectElement.Rs = 0.;
        ioIntersectElement.Tp = 1.;
        ioIntersectElement.Ts = 1.;
        return false;
	}
	
    double sinThetaTo = sinThetaFrom / m;
    CheckDoubleValue_(sinThetaTo);

	
    if (sinThetaTo >= 1.) {
        // Past critical angle, totally reflected
        ioIntersectElement.Rp = 1.;
        ioIntersectElement.Rs = 1.;
        ioIntersectElement.Tp = 0.;
        ioIntersectElement.Ts = 0.;
        return true;
    }

    MyAssert_(sinThetaTo <= 1);
    double cosThetaTo = sqrt(1.-sinThetaTo*sinThetaTo);
    CheckDoubleValue_(cosThetaTo);

    // Fresnel coefficients for fields (not intensities)
    ioIntersectElement.Tp = 2. * cosThetaFrom / (cosThetaTo + m * cosThetaFrom);
    ioIntersectElement.Ts = 2. * cosThetaFrom / (cosThetaFrom + m * cosThetaTo);
    ioIntersectElement.Rp = ioIntersectElement.Tp * m - 1.;
    ioIntersectElement.Rs = ioIntersectElement.Ts - 1.;
	
    CheckDoubleValue_(ioIntersectElement.Tp);
    CheckDoubleValue_(ioIntersectElement.Ts);
    CheckDoubleValue_(ioIntersectElement.Rp);
    CheckDoubleValue_(ioIntersectElement.Rs);

    return false;
}

void
MCObject::NormalizationFactorForInstensityStatistics(double inNorm)
{
    gNormalizationIntensity = inNorm;
}

void
MCObject::SetDetectionStatisticsOnSurfaceElements(unsigned long inDetection)
{
	vector<SurfaceElement*> surfaceElements;
	GetSurfaceElements(surfaceElements);

	unsigned int i;
	for ( i = 0; i < surfaceElements.size(); i++) {
		surfaceElements[i]->ClearMemory();
		surfaceElements[i]->mDetection = inDetection;
		surfaceElements[i]->init();
	}
}

void
MCObject::SetAcceptBallisticPhotons(bool ballistic)
{
	mBallisticPhotons = ballistic;
}

void
MCObject::ScoreOnSurface(Photon *ioPhoton, IntersectElement& inIntersectElement)
{
	if ( ioPhoton->GetNumberOfScatteringEvents() == 0 && ! mBallisticPhotons ) {
		return;
	}
	
	inIntersectElement.surfaceElement->ScoreOnSurface(ioPhoton, inIntersectElement);
}

void
MCObject::TransmitThroughInterface(Photon* ioPhoton, IntersectElement& inIntersectElement)
{

    ioPhoton->TransmitThrough(inIntersectElement.normal,
                              inIntersectElement.indexFrom,
                              inIntersectElement.indexTo,
                              inIntersectElement.Rp,
                              inIntersectElement.Rs,
                              inIntersectElement.Tp,
                              inIntersectElement.Ts);

	ScoreOnSurface(ioPhoton, inIntersectElement);

    
}

void
MCObject::ReflectAtInterface(Photon* ioPhoton, IntersectElement& inIntersectElement)
{
    ioPhoton->ReflectAtInterface(inIntersectElement.normal,
                                 inIntersectElement.indexFrom,
                                 inIntersectElement.indexTo,
                                 inIntersectElement.Rp,
                                 inIntersectElement.Rs,
                                 inIntersectElement.Tp,
                                 inIntersectElement.Ts);
}

void
MCObject::ScoreInVolume(Photon *ioPhoton, double inEnergyDeposited, bool pastHistory)
{
    if (mWorld != NULL) {
        mWorld->ScoreInVolume(ioPhoton, inEnergyDeposited, pastHistory);
    } else {
        if (mEnergy != NULL || mFluence != NULL) {
            RealV pos;
			//double weight;
			
			if ( pastHistory ) {
				;
			} else {
				pos = ioPhoton->GetLocalPosition();

				long xBin = WhichBin(pos.x, Vol_XMin, Vol_XMax, Vol_Nx);
				long yBin = WhichBin(pos.y, Vol_YMin, Vol_YMax, Vol_Ny);
				long zBin = WhichBin(pos.z, Vol_ZMin, Vol_ZMax, Vol_Nz);
				
				if ( xBin < 0 || yBin < 0 || zBin < 0 ) 
					return;
				
				if (mFluence) // Integrated 
					mFluence[To1DIndexFrom4D_(xBin,yBin,zBin,0, Vol_Nx, Vol_Ny, Vol_Nz)] += ioPhoton->GetWeight();
				if (mEnergy && inEnergyDeposited != 0) 
					mEnergy[xBin][yBin][zBin] += inEnergyDeposited;
			}
        }
    }    
}


void
MCObject::DecreaseWeight(Photon *ioPhoton)
{
    MyAssert_(mRandomScatterer != NULL);

	double mu_t;
	
	if ( (mu_t = mRandomScatterer->GetTotalExtinctionCoefficient(ioPhoton)) != 0) {
		double energyDeposited = ioPhoton->GetWeight() * mRandomScatterer->GetAbsorptionCoefficient(ioPhoton) / mu_t;
		ScoreInVolume(ioPhoton, energyDeposited);
		ioPhoton->DecreaseWeightBy(energyDeposited);
        mAbsorbedPhotons += energyDeposited;
	}
	
    mRandomScatterer->Roulette(ioPhoton);

}

void
MCObject::MoveWithinObject(Photon *ioPhoton, double inDist)
{
    CheckDoubleValue_(inDist);
    
	double index = GetIndexMedium(ioPhoton);
	ioPhoton->MoveBy(inDist, index);
	
    if (rotationPerCmInClearSpace != 0.)
        ioPhoton->RotatePolarizationStateBy(rotationPerCmInClearSpace * inDist);

//	double wavelength = ioPhoton->GetWavelength();
//	ioPhoton->IncreasePhaseDueToPopagation( 2. * PI * index / wavelength * inDist);
}


void
MCObject::ScatterInObject(Photon *ioPhoton, double inTheta, double inPhi)
{
    MyAssert_(mRandomScatterer != NULL);

    ioPhoton->ScatterBy(inTheta, inPhi, mRandomScatterer->GetMuellerMatrix());

    ioPhoton->NormalizeStokesV();
}

/*
bool
MCObject::IsGeometryConsistent()
{
    bool ok = true;
    if (mRandomScatterer != NULL) {
        if ( dynamic_cast<MCRandomScattererKaplan*>(mRandomScatterer) != NULL || 
             dynamic_cast<MCRandomScattererJaillon*>(mRandomScatterer) != NULL) {
            // HG scatterer does not keep track of index of refraction.
            if (mRandomScatterer->GetIndexMedium(NULL) != mIndexMedium) {
                clog << RightNow() << "Error: Index of medium (" << mRandomScatterer->GetIndexMedium(NULL) << ") in scatterer and object (" << mIndexMedium <<  ") should match in object " << mObjectID << " (name=" << mObjectName << ")" <<endl;
                ok = false;
            }
        }
    }
    
	long howMany = mSurfaceElements.size();
	
	if (howMany != 0) {
		for ( long i = 0; i < howMany; i++) {
			SurfaceElement* se = mSurfaceElements[i];
			
			if ( ! areTheSame(se->normal.abs(), 1, 6) ) {
				clog << RightNow() << "Normal not unitary (is actually " << se->normal.abs() << ") for surface element " << i  << " of object " << mObjectID << "(" << mObjectName << ") although it should " << endl;
	              	ok = false;
			}
			
			RealV tempNormal = RealV::CrossProduct(se->a, se->b);
			tempNormal.normalize();
			
			if ( ! areTheSame(DotProduct_(tempNormal, se->normal), 1,  6) ) {
				clog << RightNow() << "A cross B is not the normal vector for surface element " << i  << " of object " << mObjectID << "(" << mObjectName << ") although it should " << endl;
	              	ok = false;
			}
        }
	}

    return ok;
}
*/

bool
MCObject::IsGeometryConsistent()
{
    bool ok = true;
			
	vector<SurfaceElement*> listOfSurfaceElements;
	long howMany;
			
	PrintMessage_("Entering check object: " << mObjectID);
	howMany = GetSurfaceElements(listOfSurfaceElements);
	
	double tolerance = 1e-4;
	if (howMany != 0) {
		for ( long i = 0; i < howMany; i++) {
			// We need to rest this first: surfaceelements are shared between objects.
			SurfaceElement &se = *listOfSurfaceElements[i];
			se.oaFound = false;
			se.abFound = false;
			se.boFound = false;
			
			if (se.surfaceShape == kParallelogram) {
				clog << RightNow() << "Skipping geometry check for object " << mObjectName << " since it contains parallelograms." << endl;
				return true;
			}
		}
		
		for ( long i = 0; i < howMany; i++) {
			SurfaceElement &se = *listOfSurfaceElements[i];
			
			if (se.surfaceShape == kTriangle || se.surfaceShape == kParallelogram) { 
				if ( ! areTheSame(se.normal.abs(), 1, 6) ) {
					clog << RightNow() << "Normal not unitary (is actually " << se.normal.abs() << ") for surface element " << i  << " of object " << mObjectID << " although it should " << endl;
						ok = false;
				}
				
				RealV tempNormal = RealV::CrossProduct(se.a, se.b);
				tempNormal.normalize();
				
				if ( ! areTheSame(RealV::DotProduct(tempNormal, se.normal), 1,  6) ) {
					clog << RightNow() << "A cross B is not the normal vector for surface element " << i  << " of object " << mObjectID << " although it should " << endl;
						ok = false;
				}
			}
			// Checking for object closure, each surface element must share each of its edges with one and only one other surface element edge.
			// Will only check those objects with kTriangle surface elements. 
			// Implemented in a "slow and steady" fashion (can be made faster)
			
			if (se.surfaceShape == kTriangle) {
				
				// calculating vertices of current surface element
				RealV origin = se.origin;
				RealV a = se.origin + se.a;
				RealV b = se.origin + se.b;
				
				// loop  through surface elements looking for shared edges
				for ( long j = 0; j < howMany; j++) {
					//clog << "Checking element " << i << " with " << j << endl;
					// check other surface elements
					if (j != i) {
						
						SurfaceElement &seCheck = *listOfSurfaceElements[j];
						
						// true if matching vertex found
						bool aFound = false, bFound = false, oFound = false;
						
						// Calculate vertices of surface element to check
						RealV originCheck = seCheck.origin;
						RealV aCheck = seCheck.origin + seCheck.a;
						RealV bCheck = seCheck.origin + seCheck.b;
						
						// check if origin is shared by elements
						if ( (origin - originCheck).abs() < tolerance || (origin - aCheck).abs() < tolerance || (origin - bCheck).abs() < tolerance )
							oFound = true;
						
						// check if vertex a is shared by elements
						if ( (a - originCheck).abs() < tolerance || (a - aCheck).abs() < tolerance || (a - bCheck).abs() < tolerance )
							aFound = true;
						
						// check if vertex b is shared by elements
						if ( (b - originCheck).abs() < tolerance || (b - aCheck).abs() < tolerance || (b - bCheck).abs() < tolerance )
							bFound = true;
						
						if (oFound && aFound) {
							if (!se.oaFound) {
								// oa neighbour found
								se.oaNeighbour = &seCheck;
								se.oaFound = true;
							} 
							else {
								// duplicate edge neighbours found
							//	clog << RightNow() << " Surface element " << se.uniqueID << " shares edge oa with two neighbours " << (se.oaNeighbour)->uniqueID << " and " << seCheck.uniqueID << endl;
								ok = false;
							}
						}		
						if (aFound && bFound) {
							if (!se.abFound) {
								// ab neighbour found
								se.abNeighbour = &seCheck;
								se.abFound = true;
							}
							else {
								// duplicate edge neighbours found
							//	clog << RightNow() << " Surface element " << se.uniqueID << " shares edge ab with two neighbours " << (se.abNeighbour)->uniqueID << " and " << seCheck.uniqueID << endl;
								ok = false;
							}
						}	
						if (bFound && oFound) {
							if (!se.boFound) {
								// bo neighbour found
								se.boNeighbour = &seCheck;
								se.boFound = true;
							}
							else {
								// duplicate edge neighbours found
							//	clog << RightNow() << " Surface element " << se.uniqueID << " shares edge bo with two neighbours " << (se.boNeighbour)->uniqueID << " and " << seCheck.uniqueID << endl;
								ok = false;
							}
						}
					}
				}
				
				// check to make sure all edges have neighbours
				if (!se.oaFound) {
					clog << RightNow() << mObjectName << ": Surface element " << se.uniqueID << " has no neighbour along edge " << origin << " to " << a << "." << endl;
					ok = false;
				}
				if (!se.abFound) {
					clog << RightNow() << mObjectName << ": Surface element " << se.uniqueID << " has no neighbour along edge " << a << " to " << b << "." << endl;
					ok = false;
				}
				if (!se.boFound) {
					clog << RightNow() << mObjectName << ": Surface element " << se.uniqueID << " has no neighbour along edge " << origin << " to " << b << "." << endl;
					ok = false;
				}
				
				// Checking surface element normals to ensure they are consistent with neighbouring surface elements
				// warning: this will not tell you which surface element's normal is correct
				
				// finding vector for each edge of the surface element
				RealV oa = se.a;
				RealV ab = se.b - se.a;
				RealV bo = se.b;
								
				// do not check edge if neighbour not found
				if (se.oaFound) {
					SurfaceElement* neigh = se.oaNeighbour;
					RealV edge;
					if (se.objectInside == this)
						edge = se.a;
					else
						edge = -se.a;
						
					if ( RealV::DotProduct(neigh->a, edge) >= 0.999  ) {
						clog << RightNow() << " Surface element " << se.uniqueID << " has inconsistent normal with " << neigh->uniqueID << endl;
						ok = false;
					} else if (RealV::DotProduct(neigh->b, edge) <= -0.999) {
						clog << RightNow() << " Surface element " << se.uniqueID << " has inconsistent normal with " << neigh->uniqueID << endl;
						ok = false;
					} else if (RealV::DotProduct(neigh->b - neigh->a, edge) >= 0.999) {
						clog << RightNow() << " Surface element " << se.uniqueID << " has inconsistent normal with " << neigh->uniqueID << endl;
						ok = false;
					}
				}
				
				// do not check edge if neighbour not found
				if (se.abFound) {
					SurfaceElement* neigh = se.abNeighbour;
					RealV edge;
					if (se.objectInside == this)
						edge = se.b - se.a;
					else
						edge = - se.b + se.a;
					
					if ( RealV::DotProduct(neigh->a, edge) >= 0.999  ) {
						clog << RightNow() << " Surface element " << se.uniqueID << " has inconsistent normal with " << neigh->uniqueID << endl;
						ok = false;
					} else if (RealV::DotProduct(neigh->b, edge) <= -0.999) {
						clog << RightNow() << " Surface element " << se.uniqueID << " has inconsistent normal with " << neigh->uniqueID << endl;
						ok = false;
					} else if (RealV::DotProduct(neigh->b - neigh->a, edge) >= 0.999) {
						clog << RightNow() << " Surface element " << se.uniqueID << " has inconsistent normal with " << neigh->uniqueID << endl;
						ok = false;
					}
				}
				
				// do not check edge if neighbour not found
				if (se.boFound) {
					SurfaceElement* neigh = se.boNeighbour;
					RealV edge;
					if (se.objectInside == this)
						edge = -se.b;
					else
						edge = se.b;
					if ( RealV::DotProduct(neigh->a, edge) >= 0.999  ) {
						clog << RightNow() << " Surface element " << se.uniqueID << " has inconsistent normal with " << neigh->uniqueID << endl;
						ok = false;
					} else if (RealV::DotProduct(neigh->b, edge) <= -0.999) {
						clog << RightNow() << " Surface element " << se.uniqueID << " has inconsistent normal with " << neigh->uniqueID << endl;
						ok = false;
					} else if (RealV::DotProduct(neigh->b - neigh->a, edge) >= 0.999) {
						clog << RightNow() << " Surface element " << se.uniqueID << " has inconsistent normal with " << neigh->uniqueID << endl;
						ok = false;
					}
				}
				
			} else if ( se.surfaceShape == kParallelogram ) {
				clog << RightNow() << "Skipping geometry check for object ";
				clog << mObjectName << " since it contains parallelograms (element " <<  i << ")" << endl;
			}
		}
	}
	
	return ok;
}

				
string
MCObject::GetName()const
{
    return mObjectName;
}

void  
MCObject::SetName(string inName)
{
    mObjectName = inName;
}

bool 
MCObject::SetPropertiesByName(string& inMaterialName, double inWavelength)
{
    throw runtime_error("Deprecated SetPropertiesByName()");
}

void
MCObject::SetRandomScatterer(MCRandomScatterer* inRandomScatterer)
{
    MyAssert_(inRandomScatterer != NULL);

    if (inRandomScatterer != NULL) {
        if (mRandomScatterer != NULL)
            delete mRandomScatterer;
        
        mRandomScatterer = inRandomScatterer;
    }
    
}

MCRandomScatterer*
MCObject::GetRandomScatterer()
{
    MyAssert_(mRandomScatterer != NULL);

    return mRandomScatterer;
}

double 
MCObject::GetTotalExtinctionCoefficient(Photon *inPhoton)
{
    MyAssert_(mRandomScatterer != NULL);

    return mRandomScatterer->GetTotalExtinctionCoefficient(inPhoton);
}

double 
MCObject::GetScatteringCoefficient(Photon *inPhoton)
{
    MyAssert_(mRandomScatterer != NULL);

    return mRandomScatterer->GetScatteringCoefficient(inPhoton);
}

double 
MCObject::GetAbsorptionCoefficient(Photon *inPhoton)
{
    MyAssert_(mRandomScatterer != NULL);

    return mRandomScatterer->GetAbsorptionCoefficient(inPhoton); 
}

long
MCObject::GetSurfaceElements(vector<SurfaceElement*>& outSe)
{
    outSe = mSurfaceElements;
    
    return mSurfaceElements.size();
}

long
MCObject::GetSurfaceElementsCloseToSegment(vector<SurfaceElement*>& outSe, RealV inStart, RealV inEnd)
{
	RealV top, bottom;
	GetBoundingBox(top, bottom);
	if ( gInterfaceCrossingAlgorithm ==  kKDTreeAlgorithm ) {
		return GetSurfaceElementsCloseToSegmentInTree(outSe, inStart, inEnd, top, bottom);	
	} else if ( gInterfaceCrossingAlgorithm == kProximityBoxesAlgorithm ) {
		return GetSurfaceElementsCloseToSegmentWithinBox(outSe, inStart, inEnd, top, bottom, mProximityBoxList);
	} else {
		return GetSurfaceElements(outSe);
	}

}

long
MCObject::GetSurfaceElementsCloseToSegmentInTree(vector<SurfaceElement*>& outSe, RealV inStart, RealV inEnd, RealV top, RealV bottom)
{
	outSe.clear();
	
	if(mTree == NULL)
		mTree = new MCKDTree( mSurfaceElements);
	
	return mTree->GetSurfaceElementsCloseToSegment(outSe, inStart, inEnd);
}

long
MCObject::GetSurfaceElementsCloseToSegmentWithinBox(vector<SurfaceElement*>& outSe, RealV inStart, RealV inEnd, RealV top, RealV bottom, map<long,ProximityBox*>& proximityList)
{
    ProximityBox prox;
    prox.count = 0;
        
    long N = gProximitySubdivisions;
    long bit = 0;

    for (long n = N; n != 0; n >>= 1)
        bit++;
    
    unsigned long long index = 0;
    prox.binBottom.i = WhichBin(min(inStart.x,inEnd.x), bottom.x, top.x, N);
    prox.binBottom.i += prox.binBottom.i == -1 ? 1 : 0; 
    index |= prox.binBottom.i;
    prox.binBottom.j = WhichBin(min(inStart.y,inEnd.y), bottom.y, top.y, N);
    prox.binBottom.j += prox.binBottom.j == -1 ? 1 : 0; 
    index <<= bit;
    index |= prox.binBottom.j;
    prox.binBottom.k = WhichBin(min(inStart.z,inEnd.z), bottom.z, top.z, N);
    prox.binBottom.k += prox.binBottom.k == -1 ? 1 : 0; 
    index <<= bit;
    index |= prox.binBottom.k;

    prox.binTop.i = WhichBin(max(inStart.x,inEnd.x), bottom.x, top.x, N);
    prox.binTop.i += prox.binTop.i == -1 ? N : 0; 
    index <<= bit;
    index |= prox.binTop.i;
    prox.binTop.j = WhichBin(max(inStart.y,inEnd.y), bottom.y, top.y, N);
    prox.binTop.j += prox.binTop.j == -1 ? N : 0; 
    index <<= bit;
    index |= prox.binTop.j;
    prox.binTop.k = WhichBin(max(inStart.z,inEnd.z), bottom.z, top.z, N);
    prox.binTop.k += prox.binTop.k == -1 ? N : 0; 
    index <<= bit;
    index |= prox.binTop.k;
    
#ifdef __MYDEBUG
    if (prox.binBottom.i == -1 || prox.binBottom.j == -1 || prox.binBottom.k == -1 ||
        prox.binTop.i == -1 || prox.binTop.j == -1 || prox.binTop.k == -1 ) {
        PrintMessage_("Out of bound : " << bottom << " and " << top << " for points " << inStart << " and " << inEnd);
    }
#endif

    map<long,ProximityBox*>::iterator p = proximityList.find(index);
    
    if (p != proximityList.end()) {
        ProximityBox* theBox = p->second;
		theBox->count++;

		if ( theBox->mProximityBoxList != NULL) {
			int s = GetSurfaceElementsCloseToSegmentWithinBox(outSe, inStart, inEnd, theBox->top, theBox->bottom, *theBox->mProximityBoxList);
		} else if ( theBox->se.size() > gMinSurfaceElements && theBox->count > gMinPhotonCount ) {
			PrintMessage_("Subdividing element " << theBox << " currently with " << theBox->se.size() << " elements and " << theBox->count << " photons");
			theBox->mProximityBoxList = new map<long,ProximityBox*>;
			int s = GetSurfaceElementsCloseToSegmentWithinBox(outSe, inStart, inEnd, theBox->top, theBox->bottom, *theBox->mProximityBoxList);
			//DumpProximityListStatisticsFor(clog, *theBox->mProximityBoxList);
        } else {
			outSe = theBox->se;
		}
        
        return outSe.size();
    }
	
    prox.top.x = MaxBinBoundary(prox.binTop.i, bottom.x, top.x, N) + 0.001;
    prox.top.y = MaxBinBoundary(prox.binTop.j, bottom.y, top.y, N) + 0.001;
    prox.top.z = MaxBinBoundary(prox.binTop.k, bottom.z, top.z, N) + 0.001;
    
    prox.bottom.x = MinBinBoundary(prox.binBottom.i, bottom.x, top.x, N) - 0.001;
    prox.bottom.y = MinBinBoundary(prox.binBottom.j, bottom.y, top.y, N) - 0.001;
    prox.bottom.z = MinBinBoundary(prox.binBottom.k, bottom.z, top.z, N) - 0.001;
    
    GetSurfaceElementsWithinBox(prox.se, prox.bottom, prox.top);

    prox.count++;
    proximityList[index] = new ProximityBox(prox);
    
    outSe = prox.se;
    
    return outSe.size();
   
}

#define FINDMINMAX(x0,x1,x2,min,max) \
min = max = x0;   \
if(x1<min) min=x1;\
if(x1>max) max=x1;\
if(x2<min) min=x2;\
if(x2>max) max=x2;

long 
MCObject::GetSurfaceElementsWithinBox(vector<SurfaceElement*>& outSe, RealV inBottom, RealV inTop)
{
    outSe.clear();
    
    RealV delta = 0.1 * (inTop - inBottom );
    inBottom -= delta;
    inTop    += delta;
    
    RealV center = (inBottom + inTop) / 2.;

    inBottom -= center;
    inTop -= center;

    double min, max;
    
    for ( unsigned long i = 0; i < mSurfaceElements.size(); i++) {
        SurfaceElement* se = mSurfaceElements[i];

        /* We use the "Separating Axes" theorem,
           where we look at the projection of 
           the box and the triangle onto various axes.
           If the shadows don't overlap, then the triangle
           does not intersect the box.  */
        
        if ( se->surfaceShape == kTriangle ) {
            RealV v0 = se->origin - center;
            RealV v1 = se->a + se->origin - center;
            RealV v2 = se->b + se->origin - center;
            RealV c = se->b - se->a;
            
            /* X,Y and Z */
            FINDMINMAX(v0.x, v1.x, v2.x, min, max);
            if (min > inTop.x || max < inBottom.x )
                continue;
            FINDMINMAX(v0.y, v1.y, v2.y, min, max);
            if (min > inTop.y || max < inBottom.y)
                continue;
            FINDMINMAX(v0.z, v1.z, v2.z, min, max);
            if (min > inTop.z || max < inBottom.z)
                continue;
            
            /* 3 x 3 vectors */
            RealV X(1,0,0), Y(0,1,0), Z(0,0,1);
            RealV p;
            double p0,p1,p2;
            double r;
            
            /* X x a */
            p = RealV::CrossProduct(X, se->a);
            p0 = DotProduct_(p, v0);
            p1 = DotProduct_(p, v1);
            p2 = DotProduct_(p, v2);
            FINDMINMAX(p0, p1, p2, min, max);
            r = abs(p.x * inTop.x) + abs(p.y * inTop.y) + abs(p.z * inTop.z);
            if(min > r || max < -r) 
                continue;
            
            /* Y x a */
            p = RealV::CrossProduct(Y, se->a);
            p0 = DotProduct_(p, v0);
            p1 = DotProduct_(p, v1);
            p2 = DotProduct_(p, v2);
            FINDMINMAX(p0, p1, p2, min, max);
            r = abs(p.x * inTop.x) + abs(p.y * inTop.y) + abs(p.z * inTop.z);
            if(min > r || max < -r) 
                continue;
            
            /* Z x a */
            p = RealV::CrossProduct(Z, se->a);
            p0 = DotProduct_(p, v0);
            p1 = DotProduct_(p, v1);
            p2 = DotProduct_(p, v2);
            FINDMINMAX(p0, p1, p2, min, max);
            r = abs(p.x * inTop.x) + abs(p.y * inTop.y) + abs(p.z * inTop.z);
            if(min > r || max < -r) 
                continue;
            
            /* X x b */
            p = RealV::CrossProduct(X, se->b);
            p0 = DotProduct_(p, v0);
            p1 = DotProduct_(p, v1);
            p2 = DotProduct_(p, v2);
            FINDMINMAX(p0, p1, p2, min, max);
            r = abs(p.x * inTop.x) + abs(p.y * inTop.y) + abs(p.z * inTop.z);
            if(min > r || max < -r) 
                continue;
            
            /* Y x b */
            p = RealV::CrossProduct(Y, se->b);
            p0 = DotProduct_(p, v0);
            p1 = DotProduct_(p, v1);
            p2 = DotProduct_(p, v2);
            FINDMINMAX(p0, p1, p2, min, max);
            r = abs(p.x * inTop.x) + abs(p.y * inTop.y) + abs(p.z * inTop.z);
            if(min > r || max < -r) 
                continue;
            
            /* Z x b */
            p = RealV::CrossProduct(Z, se->b);
            p0 = DotProduct_(p, v0);
            p1 = DotProduct_(p, v1);
            p2 = DotProduct_(p, v2);
            FINDMINMAX(p0, p1, p2, min, max);
            r = abs(p.x * inTop.x) + abs(p.y * inTop.y) + abs(p.z * inTop.z);
            if(min > r || max < -r)
                continue;
            
            /* X x c */
            p = RealV::CrossProduct(X, c);
            p0 = DotProduct_(p, v0);
            p1 = DotProduct_(p, v1);
            p2 = DotProduct_(p, v2);
            FINDMINMAX(p0, p1, p2, min, max);
            r = abs(p.x * inTop.x) + abs(p.y * inTop.y) + abs(p.z * inTop.z);
            if(min > r || max < -r) 
                continue;
            
            /* Y x c */
            p = RealV::CrossProduct(Y, c);
            p0 = DotProduct_(p, v0);
            p1 = DotProduct_(p, v1);
            p2 = DotProduct_(p, v2);
            FINDMINMAX(p0, p1, p2, min, max);
            r = abs(p.x * inTop.x) + abs(p.y * inTop.y) + abs(p.z * inTop.z);
            if(min > r || max < -r) 
                continue;
            
            /* Z x c */
            p = RealV::CrossProduct(Z, c);
            p0 = DotProduct_(p, v0);
            p1 = DotProduct_(p, v1);
            p2 = DotProduct_(p, v2);
            FINDMINMAX(p0, p1, p2, min, max);
            r = abs(p.x * inTop.x) + abs(p.y * inTop.y) + abs(p.z * inTop.z);
            if(min > r || max < -r)
                continue;
            
            /* Normal */
            p = se->normal;
            p0 = DotProduct_(p, v0);
            p1 = DotProduct_(p, v1);
            p2 = DotProduct_(p, v2);
            FINDMINMAX(p0, p1, p2, min, max);
            r = abs(p.x * inTop.x) + abs(p.y * inTop.y) + abs(p.z * inTop.z);
            if(min > r || max < -r)
                continue;
        } 
        
        outSe.push_back(se);
    }
    
//    sort(outSe.begin(), outSe.end());
    
    PrintMessageIfLevel_(outSe.size() << " elements added for box (bottom, top): " <<  inBottom << " " << inTop, kExtremelyVerbose);
    return outSe.size();   
}

long
MCObject::GetBoundarySurfaceElements(vector<SurfaceElement*>& outSe)
{
    outSe = mSurfaceElements;
    
    return mSurfaceElements.size();
}

long
MCObject::AddSurfaceElement(SurfaceElement* inSe)
{
    for (unsigned long i = 0; i < mSurfaceElements.size(); i++) {
        if (mSurfaceElements[i] == inSe) {
            // Already in list (not an error).
            return mSurfaceElements.size() - 1;
        }
    }

    mSurfaceElements.push_back(inSe);
    
    return mSurfaceElements.size() - 1;

}

void
MCObject::SetIndexMedium(double inIndex)
{
    
	throw logic_error("Not expecting a call of SetIndexMedium()");
}

double
MCObject::GetIndexMedium(Photon *inPhoton)
{
    MyAssert_(mRandomScatterer != NULL);
    
    return mRandomScatterer->GetIndexMedium(inPhoton);
}

double
MCObject::GetRandomScatteringDistance(Photon* inPhoton)
{
    MyAssert_(mRandomScatterer != NULL);

    return mRandomScatterer->GetRandomScatteringDistance(inPhoton);
}

void
MCObject::GetRandomScatteringAngles(Photon *ioPhoton, double& outTheta, double& outPhi)
{
    MyAssert_(mRandomScatterer != NULL);

    mRandomScatterer->GetRandomScatteringAngles(outTheta, outPhi, ioPhoton);
}

RealV
MCObject::GetRandomPointInsideObjectUniformDistribution()
{
   throw runtime_error("GetRandomPointInsideObjectUniformDistribution not implemented in MCObject"); 
}

RealV
MCObject::GetRandomPointInsideObjectEnergyDistribution()
{
    throw runtime_error("GetRandomPointInsideObjectEnergyDistribution not implemented in MCObject");
}

void
MCObject::GetBoundingBox(RealV& outTop, RealV& outBottom)
{
    
    if (! cacheIsValid ) {
        if (mSurfaceElements.size() > 0) {
            outTop = mSurfaceElements[0]->origin;
            outBottom = mSurfaceElements[0]->origin;
        } else {
            outTop = RealV(0,0,0);
            outBottom = RealV(0,0,0);
        }
        
        for (size_t i = 0; i < mSurfaceElements.size(); i++) {
			if (mSurfaceElements[i]->surfaceShape == kTriangle || mSurfaceElements[i]->surfaceShape == kParallelogram) {
				RealV origin = mSurfaceElements[i]->origin;
				RealV ptA = origin + mSurfaceElements[i]->a;
				RealV ptB = origin + mSurfaceElements[i]->b;
				
				outTop.x = max(origin.x, outTop.x);
				outTop.x = max(ptA.x, outTop.x);
				outTop.x = max(ptB.x, outTop.x);
				outTop.y = max(origin.y, outTop.y);
				outTop.y = max(ptA.y, outTop.y);
				outTop.y = max(ptB.y, outTop.y);
				outTop.z = max(origin.z, outTop.z);
				outTop.z = max(ptA.z, outTop.z);
				outTop.z = max(ptB.z, outTop.z);
				
				outBottom.x = min(origin.x, outBottom.x);
				outBottom.x = min(ptA.x, outBottom.x);
				outBottom.x = min(ptB.x, outBottom.x);
				outBottom.y = min(origin.y, outBottom.y);
				outBottom.y = min(ptA.y, outBottom.y);
				outBottom.y = min(ptB.y, outBottom.y);
				outBottom.z = min(origin.z, outBottom.z);
				outBottom.z = min(ptA.z, outBottom.z);
				outBottom.z = min(ptB.z, outBottom.z);
			} else if (mSurfaceElements[i]->surfaceShape == kConic) {
				RealV origin = mSurfaceElements[i]->origin;
				RealV ptA = origin + mSurfaceElements[i]->a;
				RealV ptB = origin - mSurfaceElements[i]->a;

				outTop.x = max(origin.x, outTop.x);
				outTop.x = max(ptA.x, outTop.x);
				outTop.x = max(ptB.x, outTop.x);
				outTop.y = max(origin.y, outTop.y);
				outTop.y = max(ptA.y, outTop.y);
				outTop.y = max(ptB.y, outTop.y);
				outTop.z = max(origin.z, outTop.z);
				outTop.z = max(ptA.z, outTop.z);
				outTop.z = max(ptB.z, outTop.z);
				
				outBottom.x = min(origin.x, outBottom.x);
				outBottom.x = min(ptA.x, outBottom.x);
				outBottom.x = min(ptB.x, outBottom.x);
				outBottom.y = min(origin.y, outBottom.y);
				outBottom.y = min(ptA.y, outBottom.y);
				outBottom.y = min(ptB.y, outBottom.y);
				outBottom.z = min(origin.z, outBottom.z);
				outBottom.z = min(ptA.z, outBottom.z);
				outBottom.z = min(ptB.z, outBottom.z);
				
			} else if ( mSurfaceElements[i]->surfaceShape == kInfinitePlane ) {
				RealV origin = mSurfaceElements[i]->origin;
				RealV ptA = origin + mSurfaceElements[i]->a * INFINITE_DISTANCE;
				RealV ptB = origin + mSurfaceElements[i]->b * INFINITE_DISTANCE;
				origin = mSurfaceElements[i]->origin - mSurfaceElements[i]->a * INFINITE_DISTANCE - mSurfaceElements[i]->b * INFINITE_DISTANCE ;

				outTop.x = max(origin.x, outTop.x);
				outTop.x = max(ptA.x, outTop.x);
				outTop.x = max(ptB.x, outTop.x);
				outTop.y = max(origin.y, outTop.y);
				outTop.y = max(ptA.y, outTop.y);
				outTop.y = max(ptB.y, outTop.y);
				outTop.z = max(origin.z, outTop.z);
				outTop.z = max(ptA.z, outTop.z);
				outTop.z = max(ptB.z, outTop.z);
				
				outBottom.x = min(origin.x, outBottom.x);
				outBottom.x = min(ptA.x, outBottom.x);
				outBottom.x = min(ptB.x, outBottom.x);
				outBottom.y = min(origin.y, outBottom.y);
				outBottom.y = min(ptA.y, outBottom.y);
				outBottom.y = min(ptB.y, outBottom.y);
				outBottom.z = min(origin.z, outBottom.z);
				outBottom.z = min(ptA.z, outBottom.z);
				outBottom.z = min(ptB.z, outBottom.z);
			} else {
				ThrowRuntimeError("Code incomplete for surface element of type " << mSurfaceElements[i]->surfaceShape );
			}
        }
        
        mBoundingBoxTop = 1.001 * outTop - 0.001 * outBottom;
        mBoundingBoxBottom = 1.001 * outBottom - 0.001 * outTop;
        cacheIsValid = true;
    } else {
        outTop = mBoundingBoxTop;
        outBottom = mBoundingBoxBottom;
    }
    
    return;   
    
}

void
MCObject::RotateObject(double inAlpha, double inBeta, double inGamma)
{

	mXRotation += inAlpha;
	mYRotation += inBeta;
	mZRotation += inGamma;
	
	for ( unsigned long i = 0; i < mSurfaceElements.size(); i++) {
		mSurfaceElements[i]->RotateObject(inAlpha, inBeta, inGamma);
	}
}

void
MCObject::ScaleObject(double inX, double inY, double inZ)
{
    throw runtime_error("ScaleObject not implemented in MCObject");
}

void
MCObject::ScaleObject(double inScale)
{
    throw runtime_error("ScaleObject not implemented in MCObject");
}

void
MCObject::DumpInterfacesToFile(ostream & out)
{
    out << "<object id=\"" << mObjectID << "\" name=\"" << mObjectName << "\">\n";
    out << "<origin coordinate=\"global\">" << mGlobalOrigin<< "</origin>" << endl;

    if ( MCWorld::GetOutputProperties() & kObjectTotalAbsorbance)
		out << "<absorbance>" << mAbsorbedPhotons / gNormalizationIntensity << "</absorbance>" << endl;
    
    double index; 
    try {
        index = mRandomScatterer->GetIndexMedium(NULL);
        out << "<index>" <<  index << "</index>" << endl;
    } catch (...) {
        out << "<index>variable index with wavelength</index>" << endl;
    }

    mRandomScatterer->DumpStatsToStream(out);
    
	if (MCWorld::GetOutputProperties() & ( kObjectInterfaceTotalTransmittance | kObjectInterfaceConstructionParameters | kObjectInterfaceStokesVectorAny | kObjectInterfaceAverageAny | kObjectInterfaceTimeResolvedIntensity ) ) {
		for(unsigned long k = 0; k < mSurfaceElements.size(); k++) {
			out << "<interface id=\"" << k << "\" name=\"" << mSurfaceElements[k]->name << "\">" << endl;
			DumpInterfaceContructionParametersToStream(out, k);
			DumpStokesVToStream(out, k);
			DumpOtherStatsToStream(out, k);
			
			out << "</interface>" << endl;
		}
	}
    out << "</object>\n";
    
}

void
MCObject::DumpInterfaceContructionParametersToStream(ostream & out, long interface)
{
    if ( MCWorld::GetOutputProperties() & kObjectInterfaceConstructionParameters) {
		out << "<normal>"<< mSurfaceElements[interface]->normal << "</normal>" << endl;
		out << "<origin>"  << mSurfaceElements[interface]->origin << "</origin>" << endl;
		out << "<basevector id=\"a\">" << mSurfaceElements[interface]->a << "</basevector>" << endl;
		out << "<basevector id=\"b\">" << mSurfaceElements[interface]->b << "</basevector>" << endl;
		out << "<detectionbasis id=\"el\">" << mSurfaceElements[interface]->el << "</detectionbasis>" << endl;
		out << "<detectionbasis id=\"er\">" << mSurfaceElements[interface]->er << "</detectionbasis>" << endl;
	}
}

void
MCObject::DumpEnergyDepositionToStream(ostream & out )
{
    if ( MCWorld::GetOutputProperties() & kEnergyDistributionBinUnits) {
		out << "<!-- Each row is (x,y,z,energy), where energy is in photon weight.  (x,y,z) is given in bin number (i.e. integers) -->\n";
		out << "<energyBinNumber>\n";
		for (long i = 0; i < Vol_Nx; i++) {
			for (long j = 0; j < Vol_Ny; j++) {
				long k;
				for (k = 0; k < Vol_Nz; k++) {
					out << i << "\t" << j << "\t" << k << "\t" << mEnergy[i][j][k] << "\n";
				}
			}
		}
		out << "</energyBinNumber>\n";
	}
	
    if ( MCWorld::GetOutputProperties() & kEnergyDistributionXYZUnits) {
		out << "<!-- Each row is (x,y,z,energy), where energy is in photon weight. (x,y,z) is given in global coordinates-->\n";
		out << "<energy>\n";
		for (long i = 0; i < Vol_Nx; i++) {
			for (long j = 0; j < Vol_Ny; j++) {
				long k;
				for (k = 0; k < Vol_Nz; k++) {
					double x = BinCenter(i, Vol_XMin, Vol_XMax, Vol_Nx);
					double y = BinCenter(j, Vol_YMin, Vol_YMax, Vol_Ny);
					double z = BinCenter(k, Vol_ZMin, Vol_ZMax, Vol_Nz); 

					out << x << "\t" << y << "\t" << z << "\t" << mEnergy[i][j][k] << "\n";
				}
			}
		}
		out << "</energy>\n";
	}
	
    if ( MCWorld::GetOutputProperties() & kEnergyDistributionRawBinary) {
		out << "<!-- This energyRaw tag is a raw binary stream of " << Vol_Nz << " 16-bit images of dimensions " << Vol_Nx << " by " << Vol_Nx << ".\n";
		unsigned short bigEndian = 1;
		if ( *((unsigned char*)&bigEndian) == 1)
			out << "     The images are stored in little endian format (i.e. 1024 stored as 0x02, then 0x00). This is often referred to as the \"PC format\".\n";
		else
			out << "     The images are stored in big endian format (i.e. 1024 stored as 0x00, then 0x20).  This is often referred to as the \"Mac format\".\n";
		out << " Because XML can't accept binary data, it has been base64 encoded.  On Unix, it is pretty simple to decode using perl:\n";
		out << " xpath thisfilename.xml \"/simulation/energyRaw/text()\" | perl -w -MMIME::Base64 -e \"print decode_base64(<STDIN>);\" > decoded.raw\n";
		out << " -->\n";
		out << "<energyRaw>";
		
		long size = sizeof(unsigned short) * Vol_Nx * Vol_Ny * Vol_Nz;
		unsigned short* raw = (unsigned short*)malloc(size);
		double max = 0;
		for (long k = 0; k < Vol_Nz; k++) {
			for (long j = 0; j < Vol_Ny; j++) {
				long i;
				for (i = 0; i < Vol_Nx; i++) {
					max = mEnergy[i][j][k] > max ? mEnergy[i][j][k] : max;
				}
			}
		}
		
		
		for (long k = 0; k < Vol_Nz; k++) {
			for (long j = 0; j < Vol_Ny; j++) {
				long i;
				for (i = 0; i < Vol_Nx; i++) {
					raw[k*(Vol_Nx*Vol_Ny)+ j*(Vol_Nx) + i ] = (unsigned short) ( mEnergy[i][j][k] / max * 65535);
				}
			}
		}
		char* outBuffer = NULL;
		long outLength = 0;
		if ( XMLUtil::encode_base64((unsigned char*)raw, size, &outBuffer,  &outLength) == kNoErr) {
			for (long i = 0 ; i < outLength; i++) 
				out << outBuffer[i];
			
			free(outBuffer);
		}
		
		free(raw);
		
		out << "</energyRaw>\n";
	}

}

void
MCObject::DumpFluenceToStream(ostream & out )
{
    if ( MCWorld::GetOutputProperties() & kFluenceDistributionBinUnits) {
		out << "<!-- Each row is (x,y,z,fluence), where fluence is in photon weight.  (x,y,z) is given in bin number (i.e. integers) -->\n";
		out << "<fluenceBinNumber>\n";
		for (long i = 0; i < Vol_Nx; i++) {
			for (long j = 0; j < Vol_Ny; j++) {
				long k;
				for (k = 0; k < Vol_Nz; k++) {
					out << i << "\t" << j << "\t" << k << "\t" << mFluence[To1DIndexFrom4D_(i,j,k,0, Vol_Nx, Vol_Ny, Vol_Nz)] << "\n";
				}
			}
		}
		out << "</fluenceBinNumber>\n";
	}
	
    if ( MCWorld::GetOutputProperties() & kFluenceDistributionXYZUnits) {
		out << "<!-- Each row is (x,y,z,fluence), where fluence is in photon weight. (x,y,z) is given in global coordinates-->\n";
		out << "<fluence>\n";
		for (long i = 0; i < Vol_Nx; i++) {
			for (long j = 0; j < Vol_Ny; j++) {
				long k;
				for (k = 0; k < Vol_Nz; k++) {
					double x = BinCenter(i, Vol_XMin, Vol_XMax, Vol_Nx);
					double y = BinCenter(j, Vol_YMin, Vol_YMax, Vol_Ny);
					double z = BinCenter(k, Vol_ZMin, Vol_ZMax, Vol_Nz); 
					
					out << x << "\t" << y << "\t" << z << "\t" << mFluence[To1DIndexFrom4D_(i,j,k,0, Vol_Nx, Vol_Ny, Vol_Nz)] << "\n";
				}
			}
		}
		out << "</fluence>\n";
	}
	
    if ( MCWorld::GetOutputProperties() & kFluenceDistributionRawBinary) {
		out << "<!-- This fluenceRaw tag is a raw binary stream of " << Vol_Nz << " 16-bit images of dimensions " << Vol_Nx << " by " << Vol_Nx << ".\n";
		unsigned short bigEndian = 1;
		if ( *((unsigned char*)&bigEndian) == 1)
			out << "     The images are stored in little endian format (i.e. 1024 stored as 0x02, then 0x00). This is often referred to as the \"PC format\".\n";
		else
			out << "     The images are stored in big endian format (i.e. 1024 stored as 0x00, then 0x20).  This is often referred to as the \"Mac format\".\n";
		out << " Because XML can't accept binary data, it has been base64 encoded.  On Unix, it is pretty simple to decode using perl:\n";
		out << " xpath thisfilename.xml \"/simulation/energyRaw/text()\" | perl -w -MMIME::Base64 -e \"print decode_base64(<STDIN>);\" > decoded.raw\n";
		out << " -->\n";
		out << "<fluenceRaw>";
		
		long size = sizeof(unsigned short) * Vol_Nx * Vol_Ny * Vol_Nz;
		unsigned short* raw = (unsigned short*)malloc(size);
		double max = 0;
		for (long k = 0; k < Vol_Nz; k++) {
			for (long j = 0; j < Vol_Ny; j++) {
				long i;
				for (i = 0; i < Vol_Nx; i++) {
					max = mFluence[To1DIndexFrom4D_(i,j,k,0, Vol_Nx, Vol_Ny, Vol_Nz)] > max ? mFluence[To1DIndexFrom4D_(i,j,k,0, Vol_Nx, Vol_Ny, Vol_Nz)] : max;
				}
			}
		}
		
		
		for (long k = 0; k < Vol_Nz; k++) {
			for (long j = 0; j < Vol_Ny; j++) {
				long i;
				for (i = 0; i < Vol_Nx; i++) {
					raw[To1DIndexFrom4D_(i,j,k,0, Vol_Nx, Vol_Ny, Vol_Nz) ] = (unsigned short) ( mFluence[To1DIndexFrom4D_(i,j,k,0, Vol_Nx, Vol_Ny, Vol_Nz)] / max * 65535);
				}
			}
		}
		char* outBuffer = NULL;
		long outLength = 0;
		if ( XMLUtil::encode_base64((unsigned char*)raw, size, &outBuffer,  &outLength) == kNoErr) {
			for (long i = 0 ; i < outLength; i++) 
				out << outBuffer[i];
			
			free(outBuffer);
		}
		
		free(raw);
		
		out << "</fluenceRaw>\n";
	}
	
}

void
MCObject::DumpStokesVToStream(ostream & out, long interface )
{


	if (mSurfaceElements[interface]->mDetection & kStokesVector ) {
        // We output the values for all acceptance angle up to acceptanceMax, not the individual ranges.
        // We set up the memory here.
        
        long Nx = mSurfaceElements[interface]->Na;
        long Ny = mSurfaceElements[interface]->Nb;
        
        StokesV** summedStokesV = new StokesV*[Nx];
        for (long i = 0; i < Nx; i++) {
            summedStokesV[i] = new StokesV[Ny];
        }
        
        double transmittance = 0;
		// added by Wendy Kan 7/9/10
		double DesiredSignalIntensity = 0;
		
        for (long k = mSurfaceElements[interface]->mCosineElements - 1; k >= 0 ; k--) {
            // We accumulate the values in the table on every loop.
            for (long i = 0; i < Nx; i++) {
                for (long j = 0; j < Ny; j++) {
                    summedStokesV[i][j] += mSurfaceElements[interface]->mStokesV[i][j][k];
                    transmittance += mSurfaceElements[interface]->mStokesV[i][j][k].mI;
					//added by Wendy Kan 7/9/10
					DesiredSignalIntensity += mSurfaceElements[interface]->mDesiredSignalIntensity[i][j][k];
                }
            }
            
            double minAcceptanceCosine;
            double maxAcceptanceCosine;
            if (mSurfaceElements[interface]->mDetection == kFromInside) { 
                minAcceptanceCosine = MinBinBoundary(k, 0, 1, mSurfaceElements[interface]->mCosineElements);
                maxAcceptanceCosine = 1;
            } else if (mSurfaceElements[interface]->mDetection == kFromOutside) {
                minAcceptanceCosine = MinBinBoundary(k, 0, -1, mSurfaceElements[interface]->mCosineElements);
                maxAcceptanceCosine = -1;
            }
            
            out << "<StokesV transmittance=\"" << transmittance/gNormalizationIntensity << "\" acceptanceCosineIndex=\"" << k << "\">\n";

			//added by Wendy Kan 7/9/10
			out << "<DesiredSignalIntensity>" << DesiredSignalIntensity/gNormalizationIntensity << "</DesiredSignalIntensity>";

            out << "<!-- The transmittance is the normalized sum of photon weight that crossed this interface\n";
            out << "within the acceptance angle quoted.  The meaning of the transmittance depends on the interface you are\n";
            out << "looking at (the transmittance of the 'backward' face of a cube is the reflectance for instance)\n";
            out << "All photons with a propagation vector making an angle with the Fresnel normal to the surface\n";
            out << "smaller than " << 180/PI*acos(minAcceptanceCosine) << " degrees are included.  There are " << mSurfaceElements[interface]->mCosineElements << endl;
            out << "elements and the acceptanceCosineIndex represents the bin number.  When acceptanceCosineIndex is 0\n";
            out << "this means all photons are included -->\n";
            out << "<minAcceptanceCosine>" << minAcceptanceCosine << "</minAcceptanceCosine>\n";
            out << "<maxAcceptanceCosine>" << maxAcceptanceCosine << "</maxAcceptanceCosine>\n";

			if ( MCWorld::GetOutputProperties() & ( kObjectInterfaceIntensity | kObjectInterfaceTotalTransmittance) ) {
				out << "<I>\n";
				for (long i = 0; i < Nx; i++) {
					long j;
					for (j = 0; j < Ny-1; j++) {
						out << summedStokesV[i][j].mI << "\t" ;
					}
					out << summedStokesV[i][j].mI << endl;
				}
				out << "</I>\n";
			}
			
			if ( MCWorld::GetOutputProperties() & kObjectInterfaceStokesVectorQ) {
				out << "<Q>\n";
				for (long i = 0; i < Nx; i++) {
					long j;
					for (j = 0; j < Ny-1; j++) {
						out << summedStokesV[i][j].mQ << "\t" ;
					}
					out << summedStokesV[i][j].mQ << endl;
				}
				out << "</Q>\n";
			}
			
			if ( MCWorld::GetOutputProperties() & kObjectInterfaceStokesVectorU) {
				out << "<U>\n";
				for (long i = 0; i < Nx; i++) {
					long j;
					for (j = 0; j < Ny-1; j++) {
						out << summedStokesV[i][j].mU << "\t" ;
					}
					out << summedStokesV[i][j].mU << endl;
				}
				out << "</U>\n";
            }

			if ( MCWorld::GetOutputProperties() & kObjectInterfaceStokesVectorV) {
				out << "<V>\n";
				for (long i = 0; i < Nx; i++) {
					long j;
					for (j = 0; j < Ny-1; j++) {
						out << summedStokesV[i][j].mV << "\t" ;
					}
					out << summedStokesV[i][j].mV << endl;
				}
				out << "</V>\n";
			}
			
			if ( MCWorld::GetOutputProperties() & kObjectInterfaceStokesVectorBetaLinear) {
				out << "<betalin>\n";
				for (long i = 0; i < Nx; i++) {
					long j;
					for (j = 0; j < Ny-1; j++) {
						out << summedStokesV[i][j].degreeOfLinearPolarization()  << "\t" ;
					}
					out << summedStokesV[i][j].degreeOfLinearPolarization()  << endl;
				}
				out << "</betalin>\n";
			}
            
			if ( MCWorld::GetOutputProperties() & kObjectInterfaceStokesVectorBetaCircular) {
				out << "<betacirc>\n";
				for (long i = 0; i < Nx; i++) {
					long j;
					for (j = 0; j < Ny-1; j++) {
						out << summedStokesV[i][j].degreeOfCircularPolarization()  << "\t" ;
					}
					out << summedStokesV[i][j].degreeOfCircularPolarization()  << endl;
				}
				out << "</betacirc>\n";
			}
			
            out << "</StokesV>\n";
        }
    } else if (mSurfaceElements[interface]->mDetection & kIntensityOnly) {
        // Intensity does not discriminate with incidence angle (see ScoreOnSurface).
        long Nx = mSurfaceElements[interface]->Na;
        long Ny = mSurfaceElements[interface]->Nb;
        
        double transmittance = 0;

        if ( ! (Nx == 1 && Ny == 1)) {
            for (long i = 0; i < Nx; i++) {
                for (long j = 0; j < Ny; j++) {
                    transmittance += mSurfaceElements[interface]->mIntensityArray[i][j];
                }
            }
        } else {
            transmittance = mSurfaceElements[interface]->mIntensity;
        }
                        
        out << "<StokesV transmittance=\"" << transmittance/gNormalizationIntensity << "\" acceptanceCosineIndex=\"" << 0 << "\">\n";
        out << "<!-- The transmittance is the normalized sum of photon weight that crossed this interface\n";
        out << "within the acceptance angle quoted.  The meaning of the transmittance depends on the interface you are\n";
        out << "looking at (the transmittance of the 'backward' face of a cube is the reflectance for instance)\n";
        out << "All photons with a propagation vector making an angle with the Fresnel normal to the surface\n";
        out << "smaller than " << 90 << " degrees are included.  There is only " << 1 << endl;
        out << "element only here, which means all photons are included regardless of incidence (notation same as when full Stokes vector";
        out << "are computed -->\n";
        out << "<minAcceptanceCosine>" << 0 << "</minAcceptanceCosine>\n";
        out << "<maxAcceptanceCosine>" << 1 << "</maxAcceptanceCosine>\n";

        
		if ( MCWorld::GetOutputProperties() & ( kObjectInterfaceIntensity | kObjectInterfaceTotalTransmittance ) ) {
			out << "<I>\n";
			if ( ! (Nx == 1 && Ny == 1)) {
				for (long i = 0; i < Nx; i++) {
					long j;
					for (j = 0; j < Ny-1; j++) {
						out << mSurfaceElements[interface]->mIntensityArray[i][j] << "\t" ;
					}
					out << mSurfaceElements[interface]->mIntensityArray[i][j] << endl;
				}
			} else {
				out << mSurfaceElements[interface]->mIntensity << endl;
			}
			out << "</I>\n";
		}
		if ( (MCWorld::GetOutputProperties() & kObjectInterfaceTimeResolvedIntensity) && (mSurfaceElements[interface]->mDetection & kTimeResolved)) {
			out << "<I type=\"time-resolved\" min=\" " << mSurfaceElements[interface]->TMin << "\" max=\"" << mSurfaceElements[interface]->TMax << "\">\n";
			for (long i = 0; i < mSurfaceElements[interface]->Nt; i++) {
				out << mSurfaceElements[interface]->mTimedIntensity[i] << endl;
			}
			out << "</I>\n";
		}
		
        out << "</StokesV>\n";
    }
	
	
}

void
MCObject::DumpOtherStatsToStream(ostream & out, long interface)
{
	
	if (MCWorld::GetOutputProperties() & (kObjectInterfaceAveragePathLength | kObjectInterfaceAveragePolarizedPathLength | kObjectInterfaceAverageNumberOfScatteringEvents )) {
	//    out.setf(ios::scientific);
		// I need to rewrite this for all cases, not just when full Stokes vector are computed
		if (mSurfaceElements[interface]->mDetection & kStokesVector != 0) {
			
			// We output the values for all acceptance angle up to acceptanceMax, not the individual ranges.
			// We set up the memory here.
			
			
			long Nx = mSurfaceElements[interface]->Na;
			long Ny = mSurfaceElements[interface]->Nb;
			
			StokesV** summedStokesV = new StokesV*[Nx];
			double ** summedPath = new double*[Nx];
			double ** summedPathPol = new double*[Nx];
			double ** summedNumScatter = new double*[Nx];
			for (long i = 0; i < Nx; i++) {
				summedStokesV[i] = new StokesV[Ny];
				summedPath[i] = new double[Ny];
				summedPathPol[i] = new double[Ny]; 
				summedNumScatter[i] = new double[Ny]; 
				
				for (long j = 0; j < Ny; j++) {
					summedPath[i][j] = 0;
					summedPathPol[i][j] = 0;
					summedNumScatter[i][j] = 0;
				}
			}
			
			for (long k = mSurfaceElements[interface]->mCosineElements - 1; k >= 0 ; k--) {
				
				// We accumulate the values in the table on every loop.
				for (long i = 0; i < Nx; i++) {
					for (long j = 0; j < Ny; j++) {
						summedStokesV[i][j] += mSurfaceElements[interface]->mStokesV[i][j][k];
						summedPath[i][j] += mSurfaceElements[interface]->mPathStats[i][j][k];
						summedPathPol[i][j] += mSurfaceElements[interface]->mPathPolStats[i][j][k];
						summedNumScatter[i][j] += mSurfaceElements[interface]->mNumScatter[i][j][k];
					}
				}
				
				double minAcceptanceCosine;
				double maxAcceptanceCosine;
				if (mSurfaceElements[interface]->mDetection == kFromInside) { 
					minAcceptanceCosine = MinBinBoundary(k, 0, 1, mSurfaceElements[interface]->mCosineElements);
					maxAcceptanceCosine = 1;
				} else if (mSurfaceElements[interface]->mDetection == kFromOutside) {
					minAcceptanceCosine = MinBinBoundary(k, 0, -1, mSurfaceElements[interface]->mCosineElements);
					maxAcceptanceCosine = -1;
				}
				
				out << "<stats acceptanceCosineIndex=\"" << k << "\">\n";
				out << "<minAcceptanceCosine>" << minAcceptanceCosine << "</minAcceptanceCosine>\n";
				out << "<maxAcceptanceCosine>" << maxAcceptanceCosine << "</maxAcceptanceCosine>\n";

				if (MCWorld::GetOutputProperties() & kObjectInterfaceAveragePathLength) {
					out << "<avgpath>\n";
					
					for (long i = 0; i < Nx; i++) {
						long j;
						for (j = 0; j < Ny-1; j++) {
							double avgDist = summedStokesV[i][j].mI != 0 ? summedPath[i][j] / summedStokesV[i][j].mI: 0;
							out << avgDist  << "\t" ;
						}
						double avgDist = summedStokesV[i][j].mI != 0 ? summedPath[i][j] / summedStokesV[i][j].mI: 0;
						out << avgDist << endl;
					}
					out << "</avgpath>\n";
				}
				
				if (MCWorld::GetOutputProperties() & kObjectInterfaceAveragePolarizedPathLength) {
					out << "<avgpathpol>\n";
					
					for (long i = 0; i < Nx; i++) {
						long j;
						for (j = 0; j < Ny-1; j++) {
							double avgDist = summedStokesV[i][j].mI != 0 ? summedPathPol[i][j] / ((summedStokesV[i][j].mI+summedStokesV[i][j].mQ)/2): 0;
							out << avgDist  << "\t" ;
						}
						double avgDist = summedStokesV[i][j].mI != 0 ? summedPathPol[i][j] / ((summedStokesV[i][j].mI+summedStokesV[i][j].mQ)/2): 0;
						out << avgDist << endl;
					}
					out << "</avgpathpol>\n";
				}
				
				if (MCWorld::GetOutputProperties() & kObjectInterfaceAverageNumberOfScatteringEvents) {
					out << "<numscatter>\n";
					for (long i = 0; i < Nx; i++) {
						long j;
						for (j = 0; j < Ny-1; j++) {
							double numScatter = summedStokesV[i][j].mI != 0 ? summedNumScatter[i][j] / summedStokesV[i][j].mI: 0;
							out << numScatter  << "\t" ;
						}
						double numScatter = summedStokesV[i][j].mI != 0 ? summedNumScatter[i][j] / summedStokesV[i][j].mI: 0;
						out << numScatter << endl;
					}
					out << "</numscatter>\n";
				}
				
				out << "</stats>\n";
			}
			
			
			for (long i = 0; i < Nx; i++) {
				delete [] summedStokesV[i];
				delete [] summedPath[i];
				delete [] summedPathPol[i];
				delete [] summedNumScatter[i];
			}
			
			delete [] summedStokesV;
			delete [] summedPath;
			delete [] summedPathPol;
			delete [] summedNumScatter;
		}
	}
}

void
MCObject::DumpGeometryToStream(ostream & out, bool inAlsoIncludedObjects, long inFormat)
{
    if (inFormat == kMathematicaPolygons) {
        /* Mathematica does not take scientific notation */
        out.setf(ios::fixed);
        out << "{" ;
        
        for (unsigned long i = 0; i < mSurfaceElements.size(); i++) {
            RealV a,b, origin;
			
            /* We want everything in global coordinates */
            origin = mGlobalOrigin + mSurfaceElements[i]->origin;
            a = mSurfaceElements[i]->a;
            b = mSurfaceElements[i]->b;
			
            
            out << "Polygon[ {";
            out << "{ " << origin.x << "," << origin.y << "," << origin.z << " }, ";
            out << "{ " << a.x+origin.x << "," << a.y+origin.y << "," << a.z+origin.z << " }, ";
            if (mSurfaceElements[i]->surfaceShape == kParallelogram || mSurfaceElements[i]->surfaceShape == kInfinitePlane) {
                out << "{ " << a.x+b.x+origin.x << "," << a.y+b.y+origin.y << "," << a.z+b.z+origin.z << " }, ";
				
            }
            out << "{ " << b.x+origin.x << "," << b.y+origin.y << "," << b.z+origin.z << " } } ] ";
			
            if (i != mSurfaceElements.size() -1 ) {
                out << ", \n";
            }
			
        }
        out << "}" ;
        out.unsetf(ios::fixed);
		
    } else if (inFormat == kMatlabPatch) {
			/* Object output as a Matlab patch */
        out << "<object id=\"" << mObjectID << "\">\n" ;
        
        out << "<vertices>\n";
        for (unsigned long i = 0; i < mSurfaceElements.size(); i++) {
            RealV a,b,c, origin;
            
            origin = mSurfaceElements[i]->origin;
            a = mSurfaceElements[i]->a + origin;
            b = mSurfaceElements[i]->b + origin;
            c = mSurfaceElements[i]->a + mSurfaceElements[i]->b + origin;
            
			if (mSurfaceElements[i]->surfaceShape == kTriangle) {
                out << origin.x << "\t" << origin.y << "\t" << origin.z << "\n";
                out << a.x << "\t" << a.y << "\t" << a.z << "\n";
                out << b.x << "\t" << b.y << "\t" << b.z << "\n";
            } else if (mSurfaceElements[i]->surfaceShape == kParallelogram) {
                out << origin.x << "\t" << origin.y << "\t" << origin.z << "\n";
                out << a.x << "\t" << a.y << "\t" << a.z << "\n";
                out << c.x << "\t" << c.y << "\t" << c.z << "\n";
                out << b.x << "\t" << b.y << "\t" << b.z << "\n";
			} else if (mSurfaceElements[i]->surfaceShape == kInfinitePlane) {
                out << origin.x << "\t" << origin.y << "\t" << origin.z << "\n";
                out << a.x << "\t" << a.y << "\t" << a.z << "\n";
                out << c.x << "\t" << c.y << "\t" << c.z << "\n";
                out << b.x << "\t" << b.y << "\t" << b.z << "\n";
            }
        }
        out << "</vertices>\n";

        long c = 0;
        out << "<faces>\n";
        for (unsigned long i = 0; i < mSurfaceElements.size(); i++) {
			if (mSurfaceElements[i]->surfaceShape == kTriangle) {
                out << ++c << " " << ++c<< " " << ++c << "\n";
            } else if (mSurfaceElements[i]->surfaceShape == kParallelogram) {
                out << ++c << " " << ++c<< " " << ++c << " " << ++c << "\n";
			} else if (mSurfaceElements[i]->surfaceShape == kInfinitePlane) {
                out << ++c << " " << ++c<< " " << ++c << " " << ++c << "\n";
            }
        }
        out << "</faces>\n";

        out << "<color>\n";
        for (unsigned long i = 0; i < mSurfaceElements.size(); i++) {
            out << mSurfaceElements[i]->mIntensity << "\n";
        }
        out << "</color>\n";
        
        out << "<command>\npatch('Vertices',vertices,'Faces',faces,'FaceVertexCData',color,'FaceColor','flat','CDataMapping','scaled')\n</command>\n";
        out << "</object>\n" ;  
        
	} else if (inFormat == kGenericXML) {
			/* This is a generic XML format */
        out << "<object id=\"" << mObjectID << "\">\n" ;
		out << "\t<origin>" << mGlobalOrigin << "</origin>\n";
		out << "\t<numInterfaces>" << mSurfaceElements.size() << "</numInterfaces>\n";
		PrintMessage_("Warning: always one cosineElement");
		out << "\t<acceptanceCosineElements>" << 1 << "</acceptanceCosineElements>\n";

        for (unsigned long i = 0; i < mSurfaceElements.size(); i++) {
            RealV a,b, origin;

            origin = mSurfaceElements[i]->origin;
            a = mSurfaceElements[i]->a;
            b = mSurfaceElements[i]->b;

            
            out << "<interface>\n";
			out << "\t<name>" << mSurfaceElements[i]->name << "</name>\n";
			out << "\t<origin>" << mSurfaceElements[i]->origin << "</origin>\n";
			out << "\t<a>" << mSurfaceElements[i]->a << "</a>\n";
			out << "\t<b>" << mSurfaceElements[i]->b << "</b>\n";
			out << "\t<Na>" << mSurfaceElements[i]->Na << "</Na>\n";
			out << "\t<Nb>" << mSurfaceElements[i]->Nb << "</Nb>\n";
			out << "\t<er>" << mSurfaceElements[i]->er << "</er>\n";
			out << "\t<el>" << mSurfaceElements[i]->el << "</el>\n";
			out << "\t<normal>" << mSurfaceElements[i]->normal << "</normal>\n";
//			out << "\t<indexIn>" << mSurfaceElements[i]->indexIn << "</indexIn>\n";
//			out << "\t<indexOut>" << mSurfaceElements[i]->indexOut << "</indexOut>\n";
			
			if (mSurfaceElements[i]->surfaceShape == kTriangle)
				out << "\t<surfaceShape>kTriangle</surfaceShape>\n";
			else if (mSurfaceElements[i]->surfaceShape == kParallelogram)
				out << "\t<surfaceShape>kParallelogram</surfaceShape>\n";
			else if (mSurfaceElements[i]->surfaceShape == kInfinitePlane)
				out << "\t<surfaceShape>kInfinitePlane</surfaceShape>\n";
			
            out << "</interface>\n";

        }
        out << "</object>\n" ;
    } else if (inFormat == kWaveFrontObj) {
        //long howMany = mSurfaceElements.size();
        
        out << "g " << mObjectName << endl;
        for (unsigned long i = 0; i < mVertexList.size(); i++) {
            out << "v " << mGlobalOrigin.x + mVertexList[i].x << " " << mGlobalOrigin.y + mVertexList[i].y << " " << mGlobalOrigin.z + mVertexList[i].z << endl;
        }

        for (unsigned long i = 0; i < mNormalList.size(); i++) {
            out << "vn " << mNormalList[i].x << " " << mNormalList[i].y << " " << mNormalList[i].z << endl;
        }
        
        for (long i = 0; i <= 255; i++) {
            out << "vt 0 " << double(i)/255. << endl;
        }

        double theMax = 1, theLogMax = 1;

        /* Find maximum */
        for (unsigned long i = 0; i < mFaceList.size(); i++) {
			RealV origin, a, b;
			origin = mVertexList[mFaceList[i][0].vertex];
			a = mVertexList[mFaceList[i][1].vertex] - origin;
			b = mVertexList[mFaceList[i][2].vertex] - origin;
			float area = RealV::CrossProduct(a, b).abs();
			
            for (unsigned long j = 0; j < mFaceList[i].size(); j++) {
                theMax = theMax < mFaceList[i][j].intensity / area ? mFaceList[i][j].intensity / area : theMax;
                theLogMax = theLogMax < log10(mFaceList[i][j].intensity / area +1) ? log10(mFaceList[i][j].intensity / area +1) : theLogMax;
            }
        }

		PrintMessage_("Maximum is: " << theMax << " in object " << mObjectName );
        /* Normalize */
        for (unsigned long i = 0; i < mFaceList.size(); i++) {
			RealV origin, a, b;
			origin = mVertexList[mFaceList[i][0].vertex];
			a = mVertexList[mFaceList[i][1].vertex] - origin;
			b = mVertexList[mFaceList[i][2].vertex] - origin;
			float area = RealV::CrossProduct(a, b).abs();
			
            for (unsigned long j = 0; j < mFaceList[i].size(); j++) {
                mFaceList[i][j].texture = long(mFaceList[i][j].intensity / area / theMax * 255.) + 1;
                //mFaceList[i][j].texture = long(log10(mFaceList[i][j].intensity+1) / theLogMax * 255) + 1;
            }
        }
        
        for (unsigned long i = 0; i < mFaceList.size(); i++) {
            out << "f";
            
            for (unsigned long j = 0; j < mFaceList[i].size(); j++) {
                out << " " << mFaceList[i][j].vertex + gVertexOffset << "/" << mFaceList[i][j].texture << "/" << mFaceList[i][j].normal + gNormalOffset;
            }
            out << endl;
        }

        gVertexOffset += mVertexList.size();
        gNormalOffset += mNormalList.size();

            
    } else if ( inFormat == kBinary3DFormat ) {

		int32_t theSize = mSurfaceElements.size();

        /* Find maximum */
		float theMax = 0;
        for (long i = 0; i < theSize; i++) {
			RealV origin, a, b;
			origin = mSurfaceElements[i]->origin;
			a = mSurfaceElements[i]->a;
			b = mSurfaceElements[i]->b;
			float area;
			
			if ( mSurfaceElements[i]->surfaceShape == kTriangle )
				area = RealV::CrossProduct(a, b).abs() / 2.;
			else
				area = RealV::CrossProduct(a, b).abs();

			theMax = theMax < mSurfaceElements[i]->mIntensity / area ? mSurfaceElements[i]->mIntensity / area: theMax;
        }


#define WriteAsFloatToStream(s, x) { float value = x; s.write((char*)&value, sizeof(float)); }
#define WritePointAsFloatToStream(s, p) { WriteAsFloatToStream(s, p.x); WriteAsFloatToStream(s, p.y); WriteAsFloatToStream(s, p.z)  }

		int i;
		for (i = 0; i < theSize; i++) {
			if (  mSurfaceElements[i]->surfaceShape == kTriangle ) {
				int32_t theSize = 1001;
				out.write((char*)&theSize, sizeof(theSize));

				int32_t segments = 3;
				out.write((char*)&segments, sizeof(segments));

				RealV o = mSurfaceElements[i]->origin;
				RealV a = o + mSurfaceElements[i]->a;
				RealV b = o + mSurfaceElements[i]->b;
				RealV n = RealV::NormalizedCrossProduct(a,b);
				float area = RealV::CrossProduct(a, b).abs() / 2.;
				float intensity = mSurfaceElements[i]->mIntensity / area / theMax;
				RealV color(intensity, intensity, intensity);
				
				WritePointAsFloatToStream(out, color);
				WritePointAsFloatToStream(out, o);
				WritePointAsFloatToStream(out, n);
				WritePointAsFloatToStream(out, color);
				WritePointAsFloatToStream(out, a);
				WritePointAsFloatToStream(out, n);
				WritePointAsFloatToStream(out, color);
				WritePointAsFloatToStream(out, b);
				WritePointAsFloatToStream(out, n);
			} else if ( mSurfaceElements[i]->surfaceShape == kParallelogram ) {
				int32_t theSize = 1001;
				out.write((char*)&theSize, sizeof(theSize));

				int32_t segments = 4;
				out.write((char*)&segments, sizeof(segments));

				RealV o = mSurfaceElements[i]->origin;
				RealV a = o + mSurfaceElements[i]->a;
				RealV c = o + mSurfaceElements[i]->a + mSurfaceElements[i]->b;
				RealV b = o + mSurfaceElements[i]->b;
				RealV n = RealV::NormalizedCrossProduct(a,b);
				float area = RealV::CrossProduct(a, b).abs() / 2.;
				float intensity = mSurfaceElements[i]->mIntensity / area / theMax;
				RealV color(intensity, intensity, intensity);

				WritePointAsFloatToStream(out, color);
				WritePointAsFloatToStream(out, o);
				WritePointAsFloatToStream(out, n);
				WritePointAsFloatToStream(out, color);
				WritePointAsFloatToStream(out, a);
				WritePointAsFloatToStream(out, n);
				WritePointAsFloatToStream(out, color);
				WritePointAsFloatToStream(out, c);
				WritePointAsFloatToStream(out, n);
				WritePointAsFloatToStream(out, color);
				WritePointAsFloatToStream(out, b);
				WritePointAsFloatToStream(out, n);			
			} else {
				throw runtime_error("Can't output binary that is not triangles or parallelogram");
			}
		}	
		
	} else {
        throw runtime_error("Unknown geometry format required in MCObject::DumpGeometryToStream()");
    }
}

void
MCObject::DumpProximityListStatistics(ostream & out)
{
    out << mObjectName << endl;

	DumpProximityListStatisticsFor(out, mProximityBoxList);

}

void
MCObject::DumpProximityListStatisticsFor(ostream & out, map<long,ProximityBox*>& proximityList)
{
    
    map<long,ProximityBox*>::iterator the_iterator;
	
	out << " Has " << proximityList.size() << " elements\n";
    for( the_iterator = proximityList.begin(); the_iterator != proximityList.end(); the_iterator++) {
        ProximityBox* prox = the_iterator->second;
        out << "[" << the_iterator->first << "] ";
        out << "Si:(" << prox->binBottom.i << ", " << prox->binBottom.j << ", " << prox->binBottom.k << ") ";
        out << "Ei:(" << prox->binTop.i << ", " << prox->binTop.j << ", " << prox->binTop.k << ") ";
        out << "Sf:" << prox->bottom << " ";
        out << "Ef:" << prox->top << " ";
        out << "SE:" << prox->se.size() << " ";
        out << "C:" << prox->count << endl;
		
		if (prox->mProximityBoxList != NULL) {
			out << "This is subdivided some more into:" << endl;
			DumpProximityListStatisticsFor(out, *prox->mProximityBoxList);
		}
	}
    
}



SurfaceElement::SurfaceElement()
{
	name = "";
	uniqueID = 0;
	origin = RealV(0,0,0); 
	a = RealV(0,0,0);
	b = RealV(0,0,0);
	el = RealV(0,0,0);
	er = RealV(0,0,0);
	normal = RealV(0,0,0);
	Na = 0;
	Nb = 0;
	surfaceShape = kParallelogram;
	objectInside = NULL;
	objectOutside = NULL;
	mDetection = kNone;
	mCosineElements = 0;
	mMemoryInit = false;
	mIntensity = 0;	
	mTimedIntensity = 0;
	TMin = 0;
	TMax = 0;
	Nt = 0;
	oaFound = false;
	abFound = false;
	boFound = false;
}

SurfaceElement::~SurfaceElement()
{
	ClearMemory();
}

void
SurfaceElement::ClearMemory()
{
	if (mIntensityArray) {
		for (long i = 0; i < Na; i++) {
			if (mIntensityArray[i]) {
				delete [] mIntensityArray[i];
			}
		}
		delete[] mIntensityArray;
		mIntensityArray = 0;
	}
	
	if ( mTimedIntensity ) {
		delete[] mTimedIntensity;
		mTimedIntensity = 0;
	}
	
	if ( mStokesV ) {
		for (long i = 0; i < Na; i++) {
			for (long j = 0; j < Nb; j++) {
				delete[] mStokesV[i][j];
				delete[] mPathStats[i][j];
				delete[] mPathPolStats[i][j];
				delete[] mNumScatter[i][j];
				//added Wendy Kan 7/9/10
				delete[] mDesiredSignalIntensity[i][j];
			}
			delete[] mStokesV[i];
			delete[] mPathStats[i];
			delete[] mPathPolStats[i];
			delete[] mNumScatter[i];
			//added Wendy Kan 7/9/10
			delete[] mDesiredSignalIntensity[i];
		}
		delete[] mStokesV;
		delete[] mPathStats;
		delete[] mPathPolStats;
		delete[] mNumScatter;
		//added Wendy Kan 7/9/10
		delete[] mDesiredSignalIntensity;
		
		mStokesV = 0;
		mPathStats = 0;
		mPathPolStats = 0;
		mNumScatter = 0;
		//added Wendy kan 7/9/10
		mDesiredSignalIntensity = 0;
	}
	
	mMemoryInit = false;
}

void
SurfaceElement::init()
{
	if (mMemoryInit)
		ClearMemory();
	
    if (! mMemoryInit) {
        mMemoryInit = true;
        
        if ( mDetection & kIntensityOnly) {
            if ( !(Na == 1 && Nb == 1) ) {
                mIntensityArray = new double*[Na];
                mIntensity = 0;
				
                for (long i = 0; i < Na; i++) {
                    mIntensityArray[i] = new double[Nb];
                    for (long j = 0; j < Nb; j++) {
                        mIntensityArray[i][j] = 0;
                    }
                }
            } else {
                mIntensityArray = 0;
                mIntensity = 0;
            }
			
			if (Nt != 0) {
				mTimedIntensity = new double[Nt];
				for (long i = 0; i < Nt; i++) {
					mTimedIntensity[i] = 0;
				}
			}
        } else {
            mIntensityArray = 0;
            mIntensity = 0;
        }
        
        if ( mDetection & kStokesVector ) {
            mStokesV = new StokesV**[Na];
            mPathStats = new double**[Na];
            mPathPolStats = new double**[Na];
            mNumScatter = new double**[Na];			
			// added Wendy Kan 7/9/10
			mDesiredSignalIntensity = new double**[Na];
            for (long i = 0; i < Na; i++) {
                mStokesV[i] = new StokesV*[Nb];
                mPathStats[i] = new double*[Nb];
                mPathPolStats[i] = new double*[Nb];
                mNumScatter[i] = new double*[Nb];
				//added Wendy Kan 7/9/10
				mDesiredSignalIntensity[i] = new double*[Nb];
				
                for (long j = 0; j < Nb; j++) {
                    mStokesV[i][j] = new StokesV[mCosineElements];
                    mPathStats[i][j] = new double[mCosineElements];
                    mPathPolStats[i][j] = new double[mCosineElements];
                    mNumScatter[i][j] = new double[mCosineElements];
					//added Wendy Kan 7/9/10
					mDesiredSignalIntensity[i][j] = new double[mCosineElements];
					
                    for (long l = 0; l < mCosineElements; l++) {
                        mPathStats[i][j][l] = 0.;
                        mPathPolStats[i][j][l] = 0.;
                        mNumScatter[i][j][l] = 0.;
						//added Wendy Kan 7/9/10
						mDesiredSignalIntensity[i][j][l] = 0.;
                    }
                    
                }
            }
            
        } else {
            mStokesV = 0;
            mPathStats = 0;
            mPathPolStats = 0;
            mNumScatter = 0;
			//added Wendy Kan 7/9/10
			mDesiredSignalIntensity = 0;
        }
    }
}

void
SurfaceElement::Dump(ostream& s)
{
	s << "SurfaceElement: { uniqueID : " << uniqueID << " o,a,b" << origin << a << b << "} \n";
	
}

bool 
SurfaceElement::PathMayCrossSurfaceElement(RealV& src, RealV& d, IntersectElement& intersect)
{
	switch (surfaceShape) {
		case kInfinitePlane:
		{
			// http://www.thepolygoners.com/tutorials/lineplane/lineplane.html
			double t = src.DistanceToPlane(origin, normal, d);
			if ( ! (t < 0) && !(t > 1.)) {
				intersect.distanceToMove = t ;
				intersect.surfaceElement = this;
				intersect.interface = uniqueID;
				return true;
			}
			
		}
			break;
		case kParallelogram:
		case kTriangle:
		{
			// http://www.thepolygoners.com/tutorials/lineplane/lineplane.html
			double t = src.DistanceToPlane(origin, normal, d);
			
			// If path crosses the (infinite) interface, we check
			// to see if we actually are on the (finite) triangle/parallelogram
			// The weird condition computes faster like that
			if ( ! (t < 0) && !(t > 1.)) {
				intersect.distanceToMove = t ;
				intersect.surfaceElement = this;
				intersect.interface = uniqueID;
				return true;
			}
		}
			break;
		case kConic:
		{
			RealV ad = RealV(d.x/a.x, d.y/a.y, d.z/a.z);
			RealV as = RealV(src.x/a.x, src.y/a.y, src.z/a.z);
			RealV ao = RealV(origin.x/a.x, origin.y/a.y, origin.z/a.z);
			double A = ad.norm();
			RealV v = as - ao;
			double B = 2. * RealV::DotProduct(ad, v);
			double C = v.norm() - 1.; 
			
			// throw runtime_error("Conic surfaces not implemented yet");
			double delta = B*B-4.*A*C;
			if (delta < 0) {
				return false;
			}
			
			double t0 = (-B - sqrt(delta) )/2./A;
			if (t0 >= 0 && t0 <= 1) {
				intersect.distanceToMove = t0 ;
				
				intersect.surfaceElement = this;
				intersect.interface = uniqueID;
				return true;
			} else {
				double t1 = (-B + sqrt(delta) )/2./A;
				if (t1 >= 0 && t1 <= 1) {
					intersect.distanceToMove = t1 ;
					intersect.surfaceElement = this;
					intersect.interface = uniqueID;
					return true;
				}				
			}
		}
			break;
	}
	
	return false;
}

bool 
SurfaceElement::IsCrossing(Photon* ioPhoton, RealV& src, RealV& d, IntersectElement* p)
{
	switch (surfaceShape) {
		case kInfinitePlane:
		case kTriangle:
		case kParallelogram:
			{
				// It is possible that the intersect point does not lie on the triangle.
				// (This is done for efficiency purposes)
				// We need to check 
				// http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm
				double t = p->distanceToMove;
				RealV dir = d;
				dir.normalize();
				double inDist = d.abs();
				
				RealV vecPosSrc = src - origin;
				RealV onThePlane = vecPosSrc + t * d;
				
				double pa = DotProduct_(onThePlane,a);
				double pb = DotProduct_(onThePlane,b);
				
				const double ab = DotProduct_(a,b);
				const double aa = DotProduct_(a,a);
				const double bb = DotProduct_(b,b);
				double denominator = ab * ab  - aa * bb;
				
				double u = (ab * pb - bb * pa) / denominator;
				double v = (ab * pa - aa * pb) / denominator;
				
				// Fuzz factor to avoid sneaking between two surface elements
				double tolerance = 1e-11;
			
				if (surfaceShape != kInfinitePlane) {
					if (u >= -tolerance && u <= 1+tolerance && v >= -tolerance && v <= 1+tolerance) {
						if (surfaceShape == kTriangle && u+v > 1+tolerance) {
							return false;
						}
					} else {
						return false;
					}
				}
			
				if ( p != NULL ) {
					// We want some information about the intersection, not just whether or not it crosses
					p->u = u;
					p->v = v;
					p->distanceToMove *= inDist;
					
					// Set the output element
					p->normal   	= normal;
					p->cosine 	= RealV::NormalizedDotProduct(p->normal, dir);
									
					p->i      	= WhichBin(p->u, 0, 1, Na);
					p->j      	= WhichBin(p->v, 0, 1, Nb);
					p->objectInside 	= p->cosine >= 0 ? objectInside : objectOutside;
					p->objectOutside 	= p->cosine >= 0 ? objectOutside : objectInside;
					
					p->distanceLeftover = (inDist - p->distanceToMove);
					MyAssert_(p->distanceLeftover > 0);

					double indexFrom, indexTo;
					
					if ( ioPhoton != NULL ) {
						// We want even more info, including index and remaining paths
						if (p->cosine >= 0) {
							GetIndexMedia(ioPhoton, indexFrom, indexTo);
						} else {
							GetIndexMedia(ioPhoton, indexTo, indexFrom);
						} 
						p->indexFrom = indexFrom;
						p->indexTo = indexTo; 
						p->opticalPathLeftover = p->distanceLeftover * p->objectInside->GetTotalExtinctionCoefficient(ioPhoton);
					} else {
						p->indexFrom  	= -1;
						p->indexTo 	= -1; 
						p->opticalPathLeftover = NAN;
					}
				}
			}
			break;
			
		case kConic:
		{
			// http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm
			double t = p->distanceToMove;
			RealV dir = d;
			dir.normalize();
			double inDist = d.abs();
			
			RealV vecPosSrc = src - origin;
			RealV onTheSphereFromCenter = vecPosSrc + t * d;
			//clog << "On the sphere: " << onTheSphereFromCenter.abs() << " " << onTheSphereFromCenter << endl;
/*			MyAssert_( areTheSame(onTheSphereFromCenter.abs(), 1, 3) );
			if ( ! areTheSame(onTheSphereFromCenter.abs(), 1, 3) ) {
				clog << "Not on sphere " << onTheSphereFromCenter << " " << onTheSphereFromCenter.abs() << " with " << src << d << " t is " << t<< endl;
				return false;
			}*/
			
			if ( p != NULL ) { 
				p->distanceToMove = ( t * d).abs();
				
				// A little bit of clever vector calculus gives the normal easily:
				p->normal = 2. * onTheSphereFromCenter;
				p->normal.x /= a.x*a.x;
				p->normal.y /= a.y*a.y;
				p->normal.z /= a.z*a.z;
				p->normal.normalize();
				
				p->cosine 	= RealV::NormalizedDotProduct(p->normal, dir);
				p->u = atan2(onTheSphereFromCenter.x,onTheSphereFromCenter.y);
				p->v = atan2(sqrt(onTheSphereFromCenter.x * onTheSphereFromCenter.x + onTheSphereFromCenter.y * onTheSphereFromCenter.y),onTheSphereFromCenter.z);
				
				double indexFrom, indexTo;
				
				if (p->cosine >= 0) {
					GetIndexMedia(ioPhoton, indexFrom, indexTo);
				} else {
					GetIndexMedia(ioPhoton, indexTo, indexFrom);
				}
				
				MyAssert_(indexFrom >= 1);
				CheckDoubleValue_(indexFrom);
				MyAssert_(indexTo >= 1);
				CheckDoubleValue_(indexTo);

				p->indexFrom  	= indexFrom;
				p->indexTo 	= indexTo;
				p->i      	= WhichBin(p->u, -PI, PI, Na);
				p->j      	= WhichBin(p->v, -PI, PI, Nb);
				p->objectInside 	= p->cosine >= 0 ? objectInside : objectOutside;
				p->objectOutside 	= p->cosine >= 0 ? objectOutside : objectInside;
				
				p->distanceLeftover = (inDist - p->distanceToMove);
				MyAssert_(p->distanceLeftover > 0);
				p->opticalPathLeftover = p->distanceLeftover * p->objectInside->GetTotalExtinctionCoefficient(ioPhoton);
			}
		}
			break;
	}
	
	return true;
}

void
SurfaceElement::ScoreOnSurface(Photon *ioPhoton, IntersectElement& inIntersectElement)
{

	if (! mMemoryInit) 
		init();

	
	if (mDetection != kNone ) { 
		long cosineBin;
	    if (mDetection & kLeaving ) { 
            if (inIntersectElement.cosine < 0) {
                return;
            } else {
                cosineBin = WhichBin(inIntersectElement.cosine, 0, 1, mCosineElements);
            }
        } else if (mDetection & kEntering ) {
            if (inIntersectElement.cosine > 0) {
                return;
            } else {
                cosineBin = WhichBin(inIntersectElement.cosine, 0, -1, mCosineElements);
            }
        }
        
	    if (cosineBin == -1) {
	        return;
	    }

	    if (inIntersectElement.i < 0 || inIntersectElement.j < 0 ||
	        inIntersectElement.i > Na || inIntersectElement.j > Nb )
	        return;

	    StokesV so;
	    ioPhoton->NormalizeStokesV();
	    ioPhoton->MeasureStokesVectorInLabFrame(so, el, er, normal );
		
		MyAssert_(so.mI != 0.);
		
		if (so.mI != 0. ) {
			MyAssert_(ioPhoton->GetWeight() <= 1);
			CheckDoubleValue_(ioPhoton->GetWeight() / so.mI);
			CheckDoubleValue_(so.mI);
			
			so *= ioPhoton->GetWeight();
			
            if (mStokesV != 0)
                mStokesV[inIntersectElement.i][inIntersectElement.j][cosineBin] +=  so ;
            if (mPathStats != 0)
                mPathStats[inIntersectElement.i][inIntersectElement.j][cosineBin] += (ioPhoton->GetDistanceTraveled() * so.mI);
            if (mPathPolStats != 0)
                mPathPolStats[inIntersectElement.i][inIntersectElement.j][cosineBin] += (ioPhoton->GetDistanceTraveled() * (so.mI+so.mQ)/2.);
            if (mNumScatter != 0)
                mNumScatter[inIntersectElement.i][inIntersectElement.j][cosineBin] += (ioPhoton->GetNumberOfScatteringEvents() * so.mI);
            
			// update "singly scattered depth
			// added by Wendy Kan 7/8/10
			/*if (mDesiredSignalIntensity != 0){
				MCObject* theObject;
				theObject = objectInside;
				MCWorld* theWorld;
				theWorld = objectInside->GetWorld();
				string targetObjectName = theWorld->GetmaximumDepthTargetObjectName();
				MCObject* targetObject;
				targetObject = theWorld->FindObjectByName(targetObjectName);
				RealV lastScatteringPos = ioPhoton->GetMaxDepth();
				bool isInsideTargetObject = targetObject->IsInsideObject(lastScatteringPos);
				
				if (ioPhoton->GetNumberOfScatteringEvents() == 1 & isInsideTargetObject)
					mDesiredSignalIntensity[inIntersectElement.i][inIntersectElement.j][cosineBin] += so.mI;
			// end addition 7/8/10

			}
			*/
            mIntensity += so.mI;
            
            if (mIntensityArray != 0)
                mIntensityArray[inIntersectElement.i][inIntersectElement.j] += so.mI;

			if (mDetection & kTimeResolved) {
				int t = WhichBin(ioPhoton->GetTravelTime(), TMin, TMax, Nt);
				
				if (t != -1)
					mTimedIntensity[t] += so.mI;
			}
		}
	}    
    
}


void
SurfaceElement::RotateObject(double inAlpha, double inBeta, double inGamma)
{
	if (inAlpha != 0) {
		origin.RotateAroundX(inAlpha);
		a.RotateAroundX(inAlpha);
		b.RotateAroundX(inAlpha);
		normal.RotateAroundX(inAlpha);
		el.RotateAroundX(inAlpha);
		er.RotateAroundX(inAlpha);
	}
	
	if (inBeta != 0) {
		origin.RotateAroundY(inBeta);
		a.RotateAroundY(inBeta);
		b.RotateAroundY(inBeta);
		normal.RotateAroundY(inBeta);
		el.RotateAroundY(inBeta);
		er.RotateAroundY(inBeta);
	}

	if (inGamma != 0) {
		origin.RotateAroundZ(inGamma);
		a.RotateAroundZ(inGamma);
		b.RotateAroundZ(inGamma);
		normal.RotateAroundZ(inGamma);
		el.RotateAroundZ(inGamma);
		er.RotateAroundZ(inGamma);
	}	
	
}

bool
SurfaceElement::GetIndexMedia(Photon* inPhoton, double& outInsideMedium, double& outOutsideMedium)
{
    if (objectInside != NULL) {
        outInsideMedium = objectInside->GetIndexMedium(inPhoton);
    } else {
        outInsideMedium = 1;
    }

    if (objectOutside != NULL) {
        outOutsideMedium = objectOutside->GetIndexMedium(inPhoton);
    } else {
        outOutsideMedium = 1;
    }

    return false;
}

/*
IntersectElement::IntersectElement()
{
    interface = kNoInterface;

	distanceLeftover = 0;
    opticalPathLeftover = 0;
	surfaceElement = 0;
    distanceToMove = 0;    
//
	i = -1;
    j = -1;
    normal = RealV(0,0,0);
    u = 0;
    v = 0;
    cosine = 0;
    indexFrom = 0;
    indexTo = 0;
    Rp = 0;
    Rs = 0;
    Ts = 0;
    Tp = 0;
    objectInside = NULL;
    objectOutside = NULL;
    distanceLeftover = 0;
    opticalPathLeftover = 0;
    distanceToMove = 0;    
    surfaceElement = 0;

}
*/
bool 
IntersectElement::operator<(const IntersectElement& inRhs) const
{
    return distanceToMove < inRhs.distanceToMove;
}

ProximityBox::ProximityBox()
{
    count = 0;
	itselfSubdivided = false;
	mProximityBoxList = 0;
}

ProximityBox::~ProximityBox()
{
	delete mProximityBoxList;
}

ProximityBox::ProximityBox(const ProximityBox& inOriginal)
{
    binBottom = inOriginal.binBottom;
    binTop = inOriginal.binTop;
    bottom = inOriginal.bottom;
    top = inOriginal.top;
    se = inOriginal.se;
    count = inOriginal.count;
	itselfSubdivided = inOriginal.itselfSubdivided;
	mProximityBoxList = inOriginal.mProximityBoxList;

}

ProximityBox& 
ProximityBox::operator=(const ProximityBox& inRhs)
{
    binBottom = inRhs.binBottom;
    binTop = inRhs.binTop;
    bottom = inRhs.bottom;
    top = inRhs.top;
    se = inRhs.se;
    count = inRhs.count;
	itselfSubdivided = inRhs.itselfSubdivided;
	mProximityBoxList = inRhs.mProximityBoxList;

    return *this;
}

