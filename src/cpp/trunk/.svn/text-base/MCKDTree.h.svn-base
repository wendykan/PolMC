#if HAVE_CONFIG_H
#include <config.h>
#endif

#include "MCObject.h"
#include "Photon.h"
#include <fstream>
using namespace std;



class SpatialTreeNode;

enum {
	kNoAxis = -1,
	kXAxis = 0,
	kYZPlane = 0,
	kYAxis = 1,
	kZXPlane = 1,
	kZAxis = 2,
	kXYPlane = 2,
	
};

enum {
	kBottomLeftX = 0,
	kBottomLeftY = 1,
	kBottomLeftZ = 2,
	kTopRightX = 3,
	kTopRightY = 4,
	kTopRightZ = 5
};

typedef struct BoundingBox { float corners[6]; SurfaceElement* se; } BoundingBox;

typedef struct SpatialTreeNode { 
	struct SpatialTreeNode* leftNode; 
	struct SpatialTreeNode* rightNode;
	struct SpatialTreeNode* parentNode;
	
	vector<SurfaceElement*> surfaceElements;
	BoundingBox boundingBox;
	double dividingPosition;
	int planeNormal;
} SpatialTreeNode;

#define tripleMin(x,y,z) min(min(x,y),z)
//( (x) < (y) && (x) < (z) ? (x) : (y) < (x) && (y) < (z) ? (y) : (z) )
#define tripleMax(x,y,z) max(max(x,y),z)
//( (x) > (y) && (x) > (z) ? (x) : (y) > (x) && (y) > (z) ? (y) : (z) )

#define kMinimumNumberOfElements 5
#define kMaxTreeLevel 15

class MCKDTree {
public:
	MCKDTree(vector<SurfaceElement*> allSurfaceElements);
	MCKDTree(vector<SurfaceElement*> allSurfaceElements, int inMinimumNumberOfElements);
	~MCKDTree();
	void FreeMemoryOfChildrenNode(SpatialTreeNode* node);
	
	BoundingBox BoundingBoxFromSurfaceElement( SurfaceElement* theSurfaceElement);
	BoundingBox BoundingBoxFromManySurfaceElements( vector<SurfaceElement*> surfaceElements);
	BoundingBox BoundingBoxFromBoundingBoxes( BoundingBox box1, BoundingBox box2);
	
	SpatialTreeNode* BuildNode(BoundingBox& inBox, vector<SurfaceElement*> surfaceElements, int inPlane, SpatialTreeNode* inParent);
	int GetDividingPosition(BoundingBox& inBox, vector<SurfaceElement*>& surfaceElements, int inPlane);
	int GetNextPlane(BoundingBox& inBox, vector<SurfaceElement*>& surfaceElements, int inPlane);
	size_t GetSurfaceElementsCloseToSegment(vector<SurfaceElement*>& outSe, RealV inStart, RealV inEnd);
	size_t Search(SpatialTreeNode* node,  RealV inStart, RealV inEnd, vector<SurfaceElement*>& outSe);
	void Dump(ostream& s);
	void DumpNodeAndChildren(ostream& s, SpatialTreeNode* inNode);
	int NodeLevel(SpatialTreeNode* inNode);
	
protected:
	SpatialTreeNode* mRootNode;
	int mMinimumNumberOfElements;
	
};
