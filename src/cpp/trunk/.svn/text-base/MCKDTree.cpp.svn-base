#if HAVE_CONFIG_H
#include <config.h>
#endif

#include "MCKDTree.h"
#include "MCObject.h"
#include "Photon.h"
#include <fstream>

using namespace std;

MCKDTree::~MCKDTree()
{
	FreeMemoryOfChildrenNode(mRootNode);
	delete mRootNode;
}

void
MCKDTree::FreeMemoryOfChildrenNode(SpatialTreeNode* node) 
{
	if ( node->leftNode != NULL ) {
		FreeMemoryOfChildrenNode(node->leftNode);
		delete node->leftNode;
	}
	if ( node->rightNode != NULL ) {
		FreeMemoryOfChildrenNode(node->rightNode);
		delete node->rightNode;
	}
}

BoundingBox 
MCKDTree::BoundingBoxFromBoundingBoxes( BoundingBox box1, BoundingBox box2)
{
	BoundingBox theBox;
	
	theBox.corners[kBottomLeftX] = min(box1.corners[kBottomLeftX], box2.corners[kBottomLeftX]);
	theBox.corners[kBottomLeftY] = min(box1.corners[kBottomLeftY], box2.corners[kBottomLeftY]);
	theBox.corners[kBottomLeftZ] = min(box1.corners[kBottomLeftZ], box2.corners[kBottomLeftZ]);

	theBox.corners[kTopRightX] = max(box1.corners[kTopRightX], box2.corners[kTopRightX]);
	theBox.corners[kTopRightY] = max(box1.corners[kTopRightY], box2.corners[kTopRightY]);
	theBox.corners[kTopRightZ] = max(box1.corners[kTopRightZ], box2.corners[kTopRightZ]);
	theBox.se = NULL;
	
	return theBox;

}

BoundingBox 
MCKDTree::BoundingBoxFromSurfaceElement( SurfaceElement* theSurfaceElement)
{
	BoundingBox theBox;
		
	RealV origin = theSurfaceElement->origin;
	RealV a = theSurfaceElement->a;
	RealV b = theSurfaceElement->b;
	
	theBox.corners[kBottomLeftX] = tripleMin(origin.x, origin.x + a.x, origin.x + b.x);
	theBox.corners[kBottomLeftY] = tripleMin(origin.y, origin.y + a.y, origin.y + b.y);
	theBox.corners[kBottomLeftZ] = tripleMin(origin.z, origin.z + a.z, origin.z + b.z);

	theBox.corners[kTopRightX] = tripleMax(origin.x, origin.x + a.x, origin.x + b.x);
	theBox.corners[kTopRightY] = tripleMax(origin.y, origin.y + a.y, origin.y + b.y);
	theBox.corners[kTopRightZ] = tripleMax(origin.z, origin.z + a.z, origin.z + b.z);
	theBox.se = theSurfaceElement;
	
	return theBox;

}

BoundingBox 
MCKDTree::BoundingBoxFromManySurfaceElements( vector<SurfaceElement*> surfaceElements)
{
	BoundingBox theOverallBox = { {0.,0.,0.,0.,0.,0.}, NULL};
	for ( unsigned long i = 0; i < surfaceElements.size(); i++) {
		BoundingBox theBox = BoundingBoxFromSurfaceElement(surfaceElements[i]);
		theOverallBox = BoundingBoxFromBoundingBoxes(theOverallBox, theBox);
	}
	
	return theOverallBox;

}


MCKDTree::MCKDTree(vector<SurfaceElement*> allSurfaceElements)
{
	
	mMinimumNumberOfElements = kMinimumNumberOfElements;
	
	// We do not own the surface elements, we just keep a reference to the pointer
	BoundingBox theUniverse = BoundingBoxFromManySurfaceElements(allSurfaceElements);
	mRootNode = BuildNode( theUniverse, allSurfaceElements, kXYPlane, NULL );
}

MCKDTree::MCKDTree(vector<SurfaceElement*> allSurfaceElements, int inMinimumNumberOfElements)
{
	mMinimumNumberOfElements = inMinimumNumberOfElements;
	
	// We do not own the surface elements, we just keep a reference to the pointer
	BoundingBox theUniverse = BoundingBoxFromManySurfaceElements(allSurfaceElements);
	mRootNode = BuildNode( theUniverse, allSurfaceElements, kXYPlane, NULL );	
}

void 
MCKDTree::Dump(ostream& s)
{
	DumpNodeAndChildren(s, mRootNode);
}

void 
MCKDTree::DumpNodeAndChildren(ostream& s, SpatialTreeNode* inNode)
{
	if ( inNode->leftNode == NULL &&  inNode->rightNode == NULL  ) {
		BoundingBox bb = inNode->boundingBox;
		s << "Level : " << NodeLevel(inNode) << endl;
		s << "BoundingBox : (b: " << bb.corners[kBottomLeftX] << ", " << bb.corners[kBottomLeftY] << "," << bb.corners[kBottomLeftZ] << ") to (t: " << bb.corners[kTopRightX] << ", " << bb.corners[kTopRightY] << "," << bb.corners[kTopRightZ] << ") \n";
		s << "SurfaceElements : " << inNode->surfaceElements.size() << endl;
	} else {
		DumpNodeAndChildren(s, inNode->leftNode);
		DumpNodeAndChildren(s, inNode->rightNode);
	}
}


int 
MCKDTree::NodeLevel(SpatialTreeNode* inNode)
{
	int level = 0;
	
	while ( 1 ) {
		if ( inNode->parentNode != NULL ) {
			inNode = inNode->parentNode;
			level++;
		} else {
			break;
		}
	}
	
	return level;
}

SpatialTreeNode*
MCKDTree::BuildNode(BoundingBox& inBox, vector<SurfaceElement*> surfaceElements, int inPlane, SpatialTreeNode* inParent)
{
	SpatialTreeNode* theNode = new SpatialTreeNode;
	theNode->leftNode = NULL; 
	theNode->rightNode = NULL; 
	theNode->parentNode = inParent; 
	theNode->dividingPosition = 0;
	theNode->planeNormal = kNoAxis;
	theNode->boundingBox = inBox;
	
	if ( surfaceElements.size() < kMinimumNumberOfElements || NodeLevel(theNode) >= kMaxTreeLevel) {
		// This is a leaf, return a node with elements
		theNode->surfaceElements = surfaceElements;
	} else {
		// This is a node, pick a dividing position along axis normal to plane, and then assign left/right of dividing plane
		BoundingBox leftBox, rightBox;
		leftBox = inBox;
		rightBox = inBox;
		
		theNode->planeNormal = inPlane;
		
		int topCoordinate;
		int bottomCoordinate;
		switch ( theNode->planeNormal ) {
			case kXYPlane:
				topCoordinate = kTopRightZ;
				bottomCoordinate = kBottomLeftZ;
				break;
			case kYZPlane:
				topCoordinate = kTopRightX;
				bottomCoordinate = kBottomLeftX;
				break;
			case kZXPlane:
				topCoordinate = kTopRightY;
				bottomCoordinate = kBottomLeftY;
				break;
		}
		// This is very, very simple and should be optimized
		// We should be more clever and take median
		theNode->dividingPosition = (inBox.corners[bottomCoordinate] + inBox.corners[topCoordinate])/2; 		
		leftBox.corners[topCoordinate] = theNode->dividingPosition;
		rightBox.corners[bottomCoordinate] = theNode->dividingPosition;
		
		vector<SurfaceElement*> leftSurfaceElements;
		vector<SurfaceElement*> rightSurfaceElements;		
		
		for ( unsigned long i = 0; i < surfaceElements.size(); i++) {
			SurfaceElement* theSurfaceElement = surfaceElements[i];
			BoundingBox theSurfaceElementBox = BoundingBoxFromSurfaceElement(theSurfaceElement);
			
			if (theSurfaceElementBox.corners[topCoordinate] < theNode->dividingPosition ) {
				leftSurfaceElements.push_back(theSurfaceElement);
			} else if ( theSurfaceElementBox.corners[bottomCoordinate] > theNode->dividingPosition ) {
				rightSurfaceElements.push_back(theSurfaceElement);
			} else {
				leftSurfaceElements.push_back(theSurfaceElement);
				rightSurfaceElements.push_back(theSurfaceElement);
			}
			
		}
		
		// All elements in either left or right or both
		int nextPlane = ( inPlane + 1) % 3;
		theNode->leftNode = BuildNode(leftBox, leftSurfaceElements,  nextPlane, theNode);		
		theNode->rightNode = BuildNode(rightBox, rightSurfaceElements, nextPlane, theNode);
	}
	
	
	return theNode;
}

size_t
MCKDTree::GetSurfaceElementsCloseToSegment(vector<SurfaceElement*>& outSe, RealV inStart, RealV inEnd)
{
	if ( Search(mRootNode, inStart, inEnd, outSe) == 0 ) {
		return 0;
	} else {
		return outSe.size();
	}
	
}


size_t 
MCKDTree::Search(SpatialTreeNode* inNode,  RealV inStart, RealV inEnd, vector<SurfaceElement*>& outSe)
{
	if ( inNode->leftNode == NULL ) {
		bool found = false;
		
		if (inNode->rightNode != NULL )
			throw runtime_error("NewKDTree is improperly built.  Check constructor.");
		
		// Figure out if crossing any of inNode->surfaceElements

		vector<IntersectElement> possibleIntersection;
		RealV d = inEnd - inStart;
		for ( size_t i = 0; i < inNode->surfaceElements.size(); i++) {
			SurfaceElement* se = inNode->surfaceElements[i];
			
			IntersectElement inter;
			if ( se->PathMayCrossSurfaceElement(inStart, d, inter) ) {
				inter.interface = i;
				possibleIntersection.push_back(inter);
			}
		}
		
		/* We have a list of possible surface elements that could be intersected.
		 Currently, the vector d intersects the infinite plane in which the triangle lies
		 We need to 1) Find the closest one  and then 2) check that the intersection point is on the triangle. */
		sort(possibleIntersection.begin(), possibleIntersection.end());
		
		vector<IntersectElement>::iterator p = possibleIntersection.begin();
		
		while (p != possibleIntersection.end()) {			
			if ( p->surfaceElement->IsCrossing(NULL, inStart, d, &(*p)) ) {
				outSe.push_back(p->surfaceElement);
				break;
			}
			p++;
		}
		
		return outSe.size();
		
	} else {
		RealV pointOnPlane;
		RealV normal;
		double startComponent, endComponent;
		
		switch ( inNode->planeNormal ) {
			case kXYPlane:
				pointOnPlane = RealV(0,0,inNode->dividingPosition);
				normal = RealV(0,0,1);	
				startComponent = inStart.z;
				endComponent = inEnd.z;
				break;
			case kYZPlane:
				pointOnPlane = RealV(inNode->dividingPosition, 0,0);
				normal = RealV(1,0,0);	
				startComponent = inStart.x;
				endComponent = inEnd.x;
				break;
			case kZXPlane:
				pointOnPlane = RealV(0, inNode->dividingPosition,0);
				normal = RealV(0,1,0);	
				startComponent = inStart.y;
				endComponent = inEnd.y;
				break;
		}

		if ( startComponent < inNode->dividingPosition && endComponent < inNode->dividingPosition ) {
			return Search(inNode->leftNode, inStart, inEnd, outSe);
		} else  if ( startComponent > inNode->dividingPosition && endComponent > inNode->dividingPosition ) {
			return Search(inNode->rightNode, inStart, inEnd, outSe);
		} else {
			RealV d = inEnd - inStart;			
			double t = inStart.DistanceToPlane(pointOnPlane, normal, d);
			MyAssert_(t >= 0 && t <= 1);
			RealV v = (inStart - pointOnPlane) + 1.01 * t * d;
			MyAssert_(RealV::arePerpendicular(v, normal));
			RealV midPoint = inStart + t * d;
#ifdef __MYDEBUG
			switch ( inNode->planeNormal ) {
				case kXYPlane:
					MyAssert_(midPoint.z == inNode->dividingPosition);
					break;
				case kYZPlane:
					MyAssert_(midPoint.x == inNode->dividingPosition);
					break;
				case kZXPlane:
					MyAssert_(midPoint.y == inNode->dividingPosition);
					break;
			}
#endif
			
			
			if ( startComponent < inNode->dividingPosition ) {
				if ( Search(inNode->leftNode, inStart, midPoint, outSe) != 0 ) {
					return outSe.size();
				} else {
					return Search(inNode->rightNode, midPoint, inEnd, outSe);
				}
			} else {
				if ( Search(inNode->rightNode, inStart, midPoint, outSe) != 0 ) {
					return outSe.size();
				} else {
					return Search(inNode->leftNode, midPoint, inEnd, outSe);
				}
			}
			
		}
		
		
	}
	
	throw runtime_error("WTF?");
}
