/*!
    @header MCObject
    @abstract   This is the main class that includes all the propagation code for photons in objects
    @discussion This is a large class with mutliple functions to perform the propagation of photons.
*/


#ifndef __MCOBJECT
#define __MCOBJECT
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "RealV.h"

class MCKDTree;
class Box;

class MCObject;
class MCWorld;
class StokesV;
class Photon;
class MCRandomScatterer;
class MuellerMSpherical;
class IntersectElement;

extern long gDebugInterfaceID;

/*!
    @enum Common surface elements
    
*/

enum {
	kInfinity = -2,
    kNoInterface = -1,
    kForwardZPlane = 0,
    kBackwardZPlane,
    kPositiveXPlane,
    kNegativeXPlane,
    kPositiveYPlane,
    kNegativeYPlane
};

/*!
    @enum Types of surface elements
    
    These are the various types of surface elements used to delimit objects.  Any three dimensional object
    can be described by a set of these surface elements.  The difference between the various types is
    important when calculating whether or not a path crosses an interface.
    
*/

enum {
    kParallelogram,
    kTriangle,
	kConic,
    kInfinitePlane
};

/*
	@enum Types of detection
	
	The surface of an object can be used as a "detector".  The detector can keep track of photons
	that cross the interface either in, out or both, but can also keep track of nothing.
*/

enum {
    kNone = 0,
    kStokesVector = 1,
    kIntensityOnly = 2,
    kLeaving = 4,
    kEntering = 8,
	kTimeResolved = 16,
    
    /* Older names, kept for compatibility */
	kFromInside = kStokesVector | kLeaving,
	kFromOutside = kStokesVector | kEntering,
	kFromOutsideAndInside = kStokesVector | kEntering | kLeaving,
	kIntensityOnlyFromInside = kIntensityOnly | kLeaving,
    kIntensityOnlyFromOutside = kIntensityOnly | kEntering
};

/*
	@enum Types interface crossing events that are considered
	
    To speed up the interface crossing algorithm, we limit the search to 
    some events only.
*/

enum {
    kOnlyEntering = 10,
    kOnlyLeaving,
    kAnyEvent
};


/*!
    @enum Surface output format
    
    The class can output a description of itself in a 3D format.  Currently, we only support
    Mathematica, but that is very simple to change and it is expected that other types will
    be added later.
*/

enum {
    kMathematicaPolygons,
    kGenericXML,
	kMatlabPatch,
    kWaveFrontObj,
	kBinary3DFormat
};

/*!
    @const  PHOTON_SAFETY_DISTANCE
    
    As can be seen in PropagateInObject(), photons are always moved in two steps to avoid round off 
    errors due to moving photons to an interface and crossing it when we are not supposed to (or not crossing
    when we expect the photon to).  This should be "irrelevantly" small.  Currently, it is 10 nm.
*/

// changes made by Wendy Kan (Kort) 7/1/10
const double PHOTON_SAFETY_DISTANCE=1.0e-6;


/*!
    @const  OBJECT_SAFETY_DISTANCE
	
	A separate boundary width (safety distance) for objects; used to detect object relationships.
*/

// changes made by Wendy Kan (Kort) 7/1/10
const double OBJECT_SAFETY_DISTANCE = PHOTON_SAFETY_DISTANCE / 2.0;


/*!
    @class Vertex
 
    This is a small class in development for importing and dealing with arbitrary structures.
    It is a "vertex", identified with a vertex index, a texture index and a normal index
*/

typedef struct { long vertex, texture, normal; double intensity;} Vertex;

/*!
@class Face

This is a small class in development for importing and dealing with arbitrary structures.
It is a list of vertices, which makes up a face.
*/

typedef vector<Vertex> Face;

/*!
    @class SurfaceElement

    This embodies everything one needs to know about surface elements that bound three dimensional objects.
    This is implemented as a class to enforce initialization of all its elements upon creation.  All its elements
    are public and can be modified directly.  The constructor does not accept arguments and will initialize 
    everything to value that will trigger a numerical error.
*/

class  SurfaceElement {
public:
    /*! 
	@function SurfaceElement
	@discussion Default constructor
	*/
	SurfaceElement();

    /*! 
	@function SurfaceElement
	@discussion Default constructor
	*/
	virtual ~SurfaceElement();
    
    /*! 
		@function ClearMemory
		@discussion Clears all allocated memory from init.
		*/
	void ClearMemory();

	/*!
	@function init
	 @discussion All of the statistics memory is allocated here.  The reason this is not part of the constructor is that the subclass must allocate
				 and initialize all the surface elements before this function can be called (since we store stats for each surface element).
				 The user must call it explicitly after having set all parameters.
	 */
	
	void	init();
	
	
    /*!
	@function FreeStatsMemory
	 @discussion All of the statistics memory is freed here.  The reason this is not part of the constructor is that the subclass must allocate
				 and initialize all the surface elements before this function can be called (since we store stats for each surface element).
				 The subclass is responsible for calling it upon destruction.
	 */
	
    void	FreeStatsMemory();

    /*!
	 @function Dump
	 @discussion This function dump to stream s the content of the surface element.
	 */

	void	Dump( ostream& s);

    /*!
     @function		ScoreOnSurface
	 @discussion	All statistics are scored here after interface crossing.
	 @param			inPhoton The photon 
     @param         inIntersectElement The properties of the intersect element.
	 */

    virtual void	ScoreOnSurface(Photon *ioPhoton, IntersectElement& inIntersectElement);
	

    bool            GetIndexMedia(Photon* inPhoton, double& outInsideMedium, double& outOutsideMedium);

	void RotateObject(double inAlpha, double inBeta, double inGamma);
	
	bool PathMayCrossSurfaceElement(RealV& inSrc, RealV& d, IntersectElement& intersect);
	bool IsCrossing(Photon* ioPhoton, RealV& inSrc, RealV& d, IntersectElement* p);

	/*! @var name    This is the name of the surface element.  This is used when putting out to stream ("forward", "backward", etc...). */
   
    string name;
    
    /*! @var uniqueID    This is a unique ID that identifies this interface */
    
    long uniqueID;
    
    /*! @var origin  Origin of the surface element, in local coordinates */
    
    RealV origin;

    /*! @var normal   Normal to the surface, points outwards. a cross b = normal */

    RealV normal;

    /*! @var a   First vector that spans the surface. a cross b = normal. A linear combination of a and b can be any point on this surface. */

    RealV a;

    /*! @var b   Second vector that spans the surface. a cross b = normal. A linear combination of a and b can be any point on this surface. */

    RealV b;

	/*! @var area Area of the surface element */
	
	double area;
	
    /*! @var el  Reference axis for polariztion detection (E_parallel) */

    RealV el;

    /*! @var er  Reference axis for polariztion detection (E_perpendicular) */

    RealV er;

    /*! @var Na  Number of discrete bins along the a vector */

    long Na;

    /*! @var Nb  Number of discrete bins along the b vector */
    
    long Nb;
    
    /*! @var surfaceShape    Type of surface elements (kTriangle, kParalellogram, kInfinitePlane */
    
    short surfaceShape;

    /*! @var objectInside */
        
    MCObject* objectInside;

    /*! @var objectOutside */

    MCObject* objectOutside;
    
	/*! @var mDetection Determines what photons, if any, are kept for detection statistics */
	
	long mDetection;

    /*! @var		mCosineElements		Number of bins for storing the exit angle of photns on surfaces */

    long mCosineElements;
	
	/*! @var mStokesV		Stokes vectors used for accumulating statistics  */

    StokesV ***mStokesV;

	/*! @var mPathStats			Array used for accumulating statistics for average path length (weight is intensity * weight of photon) */

    double ***mPathStats;

	/*! @var mPathPolStats		Array used for accumulating statistics for average polarized path length (weight is (intensity + Stokes q) * weight of photon)  */

    double ***mPathPolStats;

	/*! @var mNumScatter	Array used for accumulating statistics for average number of collisions (weight is intensity * weight of photon)  */

    double ***mNumScatter;

	//added Wendy Kan 7/9/10
	double ***mDesiredSignalIntensity;

	/*! @var mTimedStokesV		Stokes vectors used for accumulating statistics including binning with time of travel */
	
    StokesV ****mTimedStokesV;
		
	/*! @var mIntensity         Total intensity that has crossed this surface element */
    
    double mIntensity;

    /*! @var mIntensityArray         Intensity that has crossed this surface element, subdivided into subelements (NA x Nb) */
    
    double** mIntensityArray;

	/*! @var mTimedIntensity	Array used for accumulating statistics for intensity with time of travel */
	
    double *mTimedIntensity;
	long Nt;
	double TMin;
	double TMax;
	
    /*! @var mVertexIndices     List of indices from a vertex list that make up the face */
    vector<long> mVertexIndices;
    
	
    /*! @var mNormalIndex       Index from a list of normal vector for this normal */
    long mNormalIndex;
    
	/*! @var oaNeighbour		pointer to the surface elements that sharing oa edge, generic object only */
	
	SurfaceElement *oaNeighbour;
	
	/*! @var abNeighbour		pointer to the surface elements that sharing ab edge, generic object only */
	
	SurfaceElement *abNeighbour;
	
	/*! @var boNeighbour		pointer to the surface elements that sharing bo edge, generic object only */
	
	SurfaceElement *boNeighbour;
	
	/*! @var oaFound 		True if neighbour edge to oa has been found, generic object only */
	
	bool  oaFound;
	
	/*! @var abFound 		True if neighbour edge to ab has been found, generic object only */
	
	bool  abFound;
	
	/*! @var boFound 		True if neighbour edge to bo has been found, generic oject only */
	
	bool  boFound;
	
	/* start: Added by Alex */
	/*! @var mUniverse			Coordinate of the volume that bounds the surface element */
	Box *mUniverse;
	/* end: Added by Alex */
	
	private:
	
	/*! @var mMemoryInit	If memory has been initialized, this is true */
	
	bool	mMemoryInit;

} ;

class IntV {
public:
    long i;
    long j;
    long k;
};

/*!
    @class ProximityBox
  
 */

class ProximityBox { 
	public:
        ProximityBox();
        ProximityBox(const ProximityBox& inOriginal);
        ProximityBox& operator=(const ProximityBox& inRhs);
        ~ProximityBox();

        IntV binBottom;
        IntV binTop;
        RealV bottom;
        RealV top;
        vector<SurfaceElement*> se;
		map<long,ProximityBox*>* mProximityBoxList;

		long count;
		bool itselfSubdivided;

};


/*!
    @class IntersectElement

    This class is a very close cousin of SurfaceElement.  It describes anything one needs to know about the intersection
    between a photon path and a surface element.  It is "more" than just a SurfaceElement since we may cross in or out,
    the Fresnel reflection depends on photon polarization, etc... It uses the same numbering scheme as surface elements.
*/

class IntersectElement {
public:
	inline IntersectElement()
	{
		interface = kNoInterface;
		
		distanceLeftover = 0;
		opticalPathLeftover = 0;
		surfaceElement = 0;
		distanceToMove = 0;    
		/*
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
		 */
	}
    
 
    bool operator<(const IntersectElement& inRhs) const;

    /*! @var interface   interface we crossed */
    
    long interface;

    /*! @var normal   Normal to the surface, points outwards */

    RealV normal;

    /*! @var u   Coefficient for "a" vector decribing position of intersect on plane */
    
    double u;
    
    /*! @var v   Coefficient for "b" vector decribing position of intersect on plane */

    double v;
    
    /*! @var cosine  Cosine of the angle between propagation direction and normal.  Negative if going in, positive if going out. */
    
    double cosine;
    
    /*! @var indexFrom  Index of refraction on the side from which the photon originates */
    
    double indexFrom;

    /*! @var indexTo  Index of refraction on the side to which the photon is going */
    
    double indexTo;
    
    /*! @var i   Bin number of position along a (between 0 and Na) */

    long i;

    /*! @var j   Bin number of position along b (between 0 and Nb) */

    long j;

    /*! @var Rp  Fresnel reflection coefficient for field along "p" polarization (electric field in the plane of incidence ) */
    
    double Rp;

    /*! @var Rs  Fresnel reflection coefficient for field along "s" polarization (electric field perp to the plane of incidence) */

    double Rs;

    /*! @var Tp  Fresnel transmission coefficient for field along "p" polarization (electric field in the plane of incidence ) */
    
    double Tp;

    /*! @var Ts  Fresnel transmission coefficient for field along "s" polarization (electric field perp to the plane of incidence) */
    double Ts;

    /*! @var distanceLeftover Distance left over to propagate after interface */
    double distanceLeftover;

    /*! @var opticalPathLeftover Optical path (mu_t times physical distance) left over to propagate after interface */
    double opticalPathLeftover;

    /*! @var distanceToMove Distance to move until interface */
    double distanceToMove;

    /*! @var objectInside  Object that owns to which this intersect/surface element belongs */
    MCObject*	objectInside;

    /*! @var objectOutside  This intersect/surface element touches another object */
    MCObject*	objectOutside;
    
    /*! @var surfaceElement The surface element being crossed */
    SurfaceElement* surfaceElement;
};


// forward declaration: 
std::ostream& operator<<(std::ostream& os, const MCObject& obj); 

/*!
    @class MCObject
	@abstract Class that takes care of propagating photons within objects. This is the top class for any object.
    @discussion The class MCObject (for Monte Carlo object) is the base class for all "scattering/turbid" objects: it takes care of most of the complexity
    involved in propagating a photon in a scattering medium.  One does not instantiate an MCObject: you instantiate a derived class (MCBox, 
    MCInfiniteLayer, etc...).   MCObject moves the photon, decreases its weight, checks to see if it crosses 
    interfaces, calculates the probability of Fresnel reflection/transmission by querying the Photon when required, etc... 
    It does not assume any particular shape for the object, it does not assume any particular formalism for the scattering interaction, it does
    not assume anything about the photon (wavelength, polarization, etc...) and does not assume a particular random number generator or algorithm 
    for arbitrary distribution.
    
    A scattering object may contain other scattering objects, which are "placed" inside this object using PlaceIntoObject() and local coordinates
    for the object.  Each specific object has its own local origin (look at the documentation, but for MCBox, it is the center of one face).
*/

class MCObject {
public:
	/*!
	@function MCObject constructor
	@discussion Default constructor for MCObject.  It will initialize everything to sane values (index in and out is zero) 
	 and assign a MCRandomScatterer with mu_s and mu_a of zero.
	 */
	
    MCObject();

	/*!
	 @function MCObject constructor
	 @discussion Constructor for MCObject that accepts parameters common to any object of any shape.
	 @param      inIndexMedium	Index of refraction of the medium contained in the object
	 @param      inIndexOutside	Index of refraction of the medium outside the object
	 @param		 inRotationPerCmInClearSpace	Optical activity, given in radian per centimer in clear space
	 @param		 inAcceptanceCosineElements		Number of bins used for storing the exit angle of a photon with the normal to an interface
	 @param		 inUseNextEventEstimator
	 */
	MCObject(double inIndexMedium, double inIndexOutside, double inRotationPerCmInClearSpace, bool inUseNextEventEstimator);

	/*!
	 @function MCObject dictionary constructor 
	 @discussion Constructor for MCObject that accepts parameters common to any object of any shape in the form of a string dictionary (or 'hash' table, or
				 associative array).  It is not recommended (use the other constructors) but can be use for simplicity.
	 @param      inDict	Associative array (map<string,string>) that must contain index_outside, index_med and may contain acceptanceCosElements, rotationPerCmInClearSpace
				 Vol_Nx and such, and useNextEventEstimator.
	 */
    
	MCObject(map<string,string>& inDict);

	/*!
	 @function init
	 @discussion Init function for MCObject that accepts parameters common to any object of any shape.
	 @param      inIndexMedium	Index of refraction of the medium contained in the object
	 @param      inIndexOutside	Index of refraction of the medium outside the object
	 @param		 inRotationPerCmInClearSpace	Optical activity, given in radian per centimer in clear space
	 @param		 inAcceptanceCosineElements		Number of bins used for storing the exit angle of a photon with the normal to an interface
	 @param		 inUseNextEventEstimator
	 */

    void init(double inIndexMedium, double inIndexOutside, double inRotationPerCmInClearSpace, bool inUseNextEventEstimator);

	/*!
		@function init  
	 @discussion Init function for MCObject that accepts parameters common to any object of any shape in the form of a string dictionary (or 'hash' table, or
																																		  associative array).  It is not recommended (use the other constructors) but can be use for simplicity.
	 @param      inDict	Associative array (map<string,string>) that must contain index_outside, index_med and may contain acceptanceCosElements, rotationPerCmInClearSpace
				 Vol_Nx and such, and useNextEventEstimator.
	 */
	
    virtual void	init(map<string,string>& inDict);

    /*!
     @function      InitConnectivityTablesFromSurfaceElements
	 @discussion	This function is under development.  It initializes a list of faces, vertices, common
     in 3D graphics programs (Matlab, .obj, etc...)
	 */
    
    bool            InitConnectivityTablesFromSurfaceElements();
    
    /*!
     @function      InitSurfaceElementsFromConnectivityTables
	 @discussion	This function is under development.  It initializes the object from a list of
     faces and vertices.
	 */
    bool            InitSurfaceElementsFromConnectivityTables();
    
	/*!
     @function	MCObject destructor
	 @discussion The destructor is responsible for destroying anything it allocates except the stats memory which must be freed
				 by the sub class (see functions below AllocateStatsMemory() and FreeStatsMemory())
	 */
    virtual ~MCObject();
    
	/*!
     @function	FinishCreate
	 @discussion This function must be called by subclasses when they are done creating themselves
     to ensure MCObject has a chance to finish creating various tables
	 */
    void    FinishCreate();

	/*!
		@function GetGlobalOrigin
	 @discussion Returns the origin of the object in global coordinates as a RealV.  The definition of where this origin is on the object depends on the subclass.
	 @result     the origin of the object in global coordinates
	 */

    virtual RealV GetGlobalOrigin();

	/*!
		@function SetGlobalOrigin
	 @discussion Sets the origin of the object in global coordinates.  This should not be called directly (use PlaceIntoObject()).
	 @param      the new origin of the object in global coordinates
	 */

    virtual void SetGlobalOrigin(RealV inOrigin);
    
	/*!
		@function PropagateInObject
	 @discussion This is the single most important function of this class, since it takes care of propagating
				 a photon until it crosses an interface an exits the object.  If other objects are included in this object
				 the photon will be moved into it upon crossing their surface. If useNextEventEstimator is set to true,
				 it could use this variance reduction technique (see code).
	 @param      ioPhoton	Pointer to a photon the content of which will be modified upon scattering and displacement.
	 @result     The last IntersectElement crossed by the photon
	 */

    virtual IntersectElement	PropagateInObject(Photon *ioPhoton, double inDist = 0);

	virtual 	IntersectElement RayTraceInObject(Photon *ioPhoton, double inDist);

    /*!
     @function GetName
	 @discussion Get name of object
	 @result     The name of the object as a string. This is used for convenience in the output file.
	 */
    
    string  GetName()const;
    
    /*!
     @function SetName
	 @discussion Set name of object.  This is used for convenience in the output file.
	 @param     The name of the object as a string
	 */
    
    void  SetName(string inName);

    /*!
     @function GetIndexMedium
	 @discussion Get the index of the medium
     @param      ioPhoton The photon
	 @result     The index of the medium
	 */

    double  GetIndexMedium(Photon *inPhoton);

    /*!
        @function SetIndexMedium
	 @discussion Set the index of the medium, an internally adjusts all the surface elements to reflect it
	 @result     inIndex The index of the medium
	 */

    void    SetIndexMedium(double inIndex);

    /*!
		@function GetRandomScatteringDistance
	 @discussion Get a random scattering distance for this photon.  This function actually calls the random scatterer class associated 
	 with this object
	 @param      ioPhoton The photon
	 @result     Random scattering distance
	 */

	double	GetRandomScatteringDistance(Photon *ioPhoton);

    /*!
		@function GetRandomScatteringAngles
	 @discussion Get a random scattering angles for this photon.  This function actually calls the random scatterer class associated 
	 with this object.  The polarization of the photon may be used for calculation.
	 @param      ioPhoton The photon
	 @param      outTheta The scattering angle
	 @param      outPhi The angle that the new scattering plane makes with the current scattering plane
	 */

    void	GetRandomScatteringAngles(Photon *ioPhoton, double& outTheta, double& outPhi);

    /*!
		@function GetTotalExtinctionCoefficient
	 @discussion Get the total extinction coefficient for this photon.  This function actually calls the random scatterer class associated with this object. 
	 @param      ioPhoton The photon
	 @result	 the total extinction coefficient. For discrete scaterers, this is mu_s + mu_a. 
	 */

    virtual double GetTotalExtinctionCoefficient(Photon *inPhoton);

    /*!
		@function GetScatteringCoefficient
	 @discussion Get the scattering coefficient for this photon.  This function actually calls the random scatterer class associated with this object. 
	 @param      ioPhoton The photon
	 @result	 the scattering coefficient. For discrete scaterers, this is mu_s. 
	 */
	
	virtual double GetScatteringCoefficient(Photon *inPhoton);

	/*!
		@function GetAbsorptionCoefficient
	 @discussion Get the absorption coefficient for this photon.  This function actually calls the random scatterer class associated with this object. 
	 @param      inPhoton The photon
	 @result	 the absorption coefficient. For discrete scaterers, this is mu_a. 
	 */

	virtual double GetAbsorptionCoefficient(Photon *inPhoton);

    /*!
        @function SetPropertiesByName
	 @discussion Sets the scattering and optical properties by "tissue type".  See code for specific tissues
     that are implemented. Ideally, this would be read from a file.
	 @param inMaterialName This is a string that corresponds to the tissue type.  See MCObject.cpp for specifics.
	 */
	
    virtual bool SetPropertiesByName(string& inMaterialName, double inWavelength);
    
	/*!
     @function SetRandomScatterer
	 @discussion Sets the object responsible for all the random interactions. MCRandomScatterer takes care of all properties and action related 
	 to the scatterers.  This MCObject does not really need to know the details of the scattering and will simply call the appropriate function of MCRandomScatter. 
	 
	 @param inRandomScatterer	The random scatterer for this object. It is owned by MCObject after the call and should not be freed. 
	 */
	
    void	SetRandomScatterer(MCRandomScatterer* inRandomScatterer);

	/*!
		@function GetRandomScatterer
	 @discussion Get the object responsible for all the random interactions. MCRandomScatterer takes care of all properties and action related 
	 to the scatterers.  This MCObject does not really need to know the details of the scattering and will simply call the appropriate function of MCRandomScatter. 
	 
	 @result	The random scatterer for this object. It is owned by MCObject after the call and should not be freed. 
	 */

	MCRandomScatterer*	GetRandomScatterer();

	void    SetProximityProperties(long inSubdivisions, long inMinSurfaceElements, long inMinPhotonCount);
    
    void    SetOuterObjectTo(MCObject* inObject);
	/*!
     @function GetSurfaceElements
	 @discussion Get an array of all surface elements that define this object (including other objects inside)
	 @param      outSe	A reference to an array pointer that will contain the surface elements.
	 @result        the number of elements.  Indexing the array goes from 0 to (numOfElements-1).
	 */
	
    virtual long GetSurfaceElements(vector<SurfaceElement*>& outSe);

	/*!
        @function GetSurfaceElementsCloseToSegment
	 @discussion Get an array of all surface elements close to line segment, using the object's bounding box as the reference
	 @param      outSe	A reference to an array pointer that will contain the surface elements.
	 @param		inStart Beginning of segment
	 @param		inEnd End of segment
	 @result        the number of elements.  Indexing the array goes from 0 to (numOfElements-1).
	 */
	
    virtual long GetSurfaceElementsCloseToSegment(vector<SurfaceElement*>& outSe, RealV inStart, RealV inEnd);

   /*!
   @function       GetSurfaceElementsCloseToSegmentInTree
   @discussion     Get an array of all surface elements close to line segment
   @param          outSe	A reference to an array pointer that will contain the surface elements.
   @param		inStart Beginning of segment
   @param		inEnd End of segment
      
   @result         the number of elements.  Indexing the array goes from 0 to (numOfElements-1).
   */
    
   virtual long GetSurfaceElementsCloseToSegmentInTree(vector<SurfaceElement*>& outSe, RealV inStart, RealV inEnd, RealV top, RealV bottom);
  
	/*!
        @function GetSurfaceElementsCloseToSegmentWithinBox
	 @discussion Get an array of all surface elements close to line segment
	 @param      outSe	A reference to an array pointer that will contain the surface elements.
	 @param		inStart Beginning of segment
	 @param		inEnd End of segment
	 @param		top Top corner of box (most positive corner)
	 @param		bottom Bottom corner of box (most negative corner)
	 @param		level How many levels into the main box this is.  The first one is zero.
	 @result        the number of elements.  Indexing the array goes from 0 to (numOfElements-1).
	 */
		
    virtual long GetSurfaceElementsCloseToSegmentWithinBox(vector<SurfaceElement*>& outSe, RealV inStart, RealV inEnd, RealV top, RealV bottom, map<long,ProximityBox*>& proximityList);
	
    /*!
        @function GetSurfaceElementsWithinBox
	 @discussion Get an array of all surface elements within a given box
	 @param      outSe	A reference to an array pointer that will contain the surface elements.
	 @result        the number of elements.  Indexing the array goes from 0 to (numOfElements-1).
	 */
	
    virtual long GetSurfaceElementsWithinBox(vector<SurfaceElement*>& outSe, RealV inBottom, RealV inTop);
    
	/*!
        @function GetBoundarySurfaceElements
	 @discussion Get an array of all surface elements that define the boundary of this object
	 @param      outSe	A reference to an array pointer that will contain the surface elements.
	 @result        the number of elements.  Indexing the array goes from 0 to (numOfElements-1).
	 */

    virtual long GetBoundarySurfaceElements(vector<SurfaceElement*>& outSe);
        
    long         AddSurfaceElement(SurfaceElement* inSe);

	/*!
	 @function	 PathCrossesInterface
	 @discussion The most critical function (in terms of speed) for the computation.  This function returns whether or not a photon
				 path crosses a surface element and computes the intersection point.  It will call recusively the same function
				 for any contained object if inOnlyInto is true.
	 @param      ioPhoton	The photon
	 @param      inDist		Distance by which the photon could travel
	 @param		 outIntersectElement Reference to an intersect element which will contain all the information about the intersection point
     @param      inCrossingEvents What crossing events are considered
	 @result	 The interface number that was crossed
	 */
	
    virtual long	PathCrossesInterface(Photon *ioPhoton, double inDist, IntersectElement& outIntersectElement, long inCrossingEvents);

	/*!
		@function	 IsTransmitted
	 @discussion The function is called with an IntersectElement to determine whether or not the photon is transmitted.
				 This is determined by calculating the Fresnel coefficients at that interface, at that angle of incidence
				 with that polarization state.  The photon isn't moved or transformed.  The photon TransmitThrough should
				 be called if the result is true, ReflectAtInterface() should be called if it's false.
	 @param      ioPhoton	The photon
	 @param		 inIntersectElement Reference to an intersect element which contains all the information about the intersection point
				 and interface.  This is obtained from PathCrossesInterface()
	 @result	 whether or not the photon should be transmitted
	 */
	
    virtual bool	IsTransmitted(Photon *ioPhoton, IntersectElement& inIntersectElement);

	/*!
		@function	 TransmitThroughInterface
	 @discussion This function actually transforms the photon upon propagation through the given interface. 
				 This means the photon is refracted, and its polarization is transformed according to the Fresnel coefficients
				 the polarization state and the angle of incidence. Note that the object doesn't actually do the transoromation:
				 it calls the photon's function TransmitThrough() with its relevant parameters to do the actul work.
				 After, it scores the photon values for4 this interface by calling ScoreOnSurface().
	 @param      ioPhoton	The photon
	 @param		 inIntersectElement Reference to an intersect element which contains all the information about the intersection point
				 and interface.  This is obtained from PathCrossesInterface()
	 */
    
	virtual void	TransmitThroughInterface(Photon* ioPhoton, IntersectElement& inIntersectElement);

	/*!
		@function	 ReflectAtInterface
	 @discussion This function actually transforms the photon upon reflection at the given interface. 
				 This means the photon is reflected, and its polarization is transformed according to the Fresnel coefficients
				 the polarization state and the angle of incidence. Note that the object doesn't actually do the transoromation:
				 it calls the photon's function ReflectAtInterface() with its relevant parameters to do the actul work.
				 ScoreOnSurface() is NOT called since the photon never crossed the interface.
	 @param      ioPhoton	The photon
	 @param		 inIntersectElement Reference to an intersect element which contains all the information about the intersection point
				 and interface.  This is obtained from PathCrossesInterface()
	 */

	virtual void	ReflectAtInterface(Photon* ioPhoton, IntersectElement& ioIntersectElement);

	/*!
		@function	 SetFresnelCoefficients
	 @discussion This sets all the necessary values for Fresnel reflections in the IntersectElement (angle of refraction,
																									 Fresnel coefficients). 
	 @param		 ioIntersectElement Reference to an intersect element which contains the interface and angle of incidence,
				 as well as the indices of refraction.  
	 @param		 returns true if the photon is beyond the critical angle and is totally internally reflected.
	 */
	
    static bool		SetFresnelCoefficients(  IntersectElement& inIntersectElement);
    
	/*!
	 @function		RotateObject
	 @discussion	This rotates the object by an angle alpha around the local x axis, beta around y and 
					gamma around z.  The sign of rotation is determined by the right-hand rule: with your thumb
					pointing in the direction of the axis, your hand folds towards the positive direction.
					If the object contains other objects, this function will fail (and throw an exception).   
	 @param			inAlpha Amount by which the object is rotated around the x axis. 
	 @param			inBeta Amount by which the object is rotated around the y axis. 
	 @param			inGamma Amount by which the object is rotated around the z axis. 
	 */
	
    virtual void    RotateObject(double inAlpha, double inBeta, double inGamma);
	
	/*!
	 @function		ScaleObject
	 @discussion	This function scales the object along x, y and z (local frame)
	 @param			inX Factor to scale along x axis
	 @param			inY Factor to scale along y axis
	 @param			inZ Factor to scale along z axis
	 */
	
    virtual void    ScaleObject(double inX, double inY, double inZ);

	/*!
	 @function		ScaleObject
	 @discussion	This function scales the object along x, y and z (local frame)
	 @param			inScale Factor to scale along all axis
	 */

	virtual void    ScaleObject(double inScale);
    
    /*!
        @function		NormalizationFactorForInstensityStatistics
	 @discussion	The transmittance at interfaces is normalized with this factor when
     the DumpXXX() functions are called.  Note that the total number of photons
     launched in this object is not equal to the total number of launched photons by the
     source: there are reflections at interfaces (usually).  However, most of the time
     we want transmittance, reflectance and absorbance in reference to the source
     intensity, not the photons that make it inside the object.
	 @param			inNorm Normalization factor.  Usually obtained from MCSource::GetNumberOfPhotonsLaunched()
     
	 */
    
    virtual void	NormalizationFactorForInstensityStatistics(double inNorm);
    
	/*!
     @function		DumpInterfacesToFile
	 @discussion	This is part of the "Dump" function series to save everything to stream (output or file).
                    This will cycle through all the interfaces (or surface elements) of the object and save
                    their information to stream.  This function is called from DumpToFile()
	 @param			out
	 */

    virtual void    DumpInterfacesToFile(ostream & out);
    
	/*!
     @function      DumpGeometryToStream	
	 @discussion	This is part of the "Dump" function series to save everything to stream (output or file).
                    This will save the geometry to stream in the format demanded by the user. The formats
                    are kMathematicaPolygons (for mathematica), kGenericXML (for this program).
	 @param			out
     @param         inAlsoIncludedObjects
     @param         inFormat
	 */
    
    virtual void    DumpGeometryToStream(ostream & out, bool inAlsoIncludedObjects, long inFormat = kMathematicaPolygons);

    /*!
        @function		IsOutsideBoundingBox
	 @discussion	This function serves a debugging purpose only.
     If the photon is outside the bounding box of the object, this function  returns true. 
     The photon is for sure outside of the object at this point.  It does not mean that if the function
     returns false, the photon is outside: the bounding box is larger than the object.  However,
     it is much faster to compute whether or not the photon is inside a box.
	 @param			ioPhoton the photon
     @result        true if photon is outisde bounding box
	 */
    
    bool	IsOutsideBoundingBox(Photon *ioPhoton);

    
	/*!
     @function		IsOutsideObject
	 @discussion	This function serves a debugging purpose only.
                    If the photon is outside the object, this function should return true.  This will actually 
                    call the function IsOutsideObject(RealV& inLocalPosition) with the position of the photon.
                    Objects should not override this function, they should override 
                    IsOutsideObject(RealV& inLocalPosition); 
	 @param			ioPhoton
     @result        true if IsOutsideObject(RealV& inLocalPosition) returns true.
	 */
    
    bool	IsOutsideObject(Photon *ioPhoton);
    
	/*!
     @function		IsOutsideObject
	 @discussion	This function serves a debugging purpose only. If the position is outside, the function should return true.
                    The default implementation of this function throws an exception.  Subclasses must override.
                    
	 @param			inLocalPosition The position to check
     @result        true if object is outside object.
	 */
    
    virtual bool	IsOutsideObject(RealV& inLocalPosition); 
    
	/*!
     @function      IsInsideObject	
	 @discussion	This function intersects the surface of the object with 6 rays starting at
	 				the origin of the input photon. The number of intersections, odd or even, will determine
	 				whether or not the photon is inside the object. If at least 5 of the 6 rays
	 				have an odd number of intersections, the photon is considered to be inside,
	 				and if at least 5 of the 6 rays have an even amount of intersections, the 
	 				photon is considered outside of the object. We allow for one "mismatch" number
	 				of intersection, because a ray might touch the surface, returning an ambiguous/wrong
	 				answer.
	 @param			ioPhoton
	 @result		Boolean wheteher or not the photon is inside	 */
    
    virtual bool	IsInsideObject(Photon *ioPhoton);

	/*!
     @function      IsInsideObject	
	 @discussion	This function determines if the provided point is inside the object.
	 This is performed by casting rays from a fake photon.  The number of intersections, odd or even, will determine
	 whether or not the photon is inside the object. If at least 5 of the 6 rays
	 have an odd number of intersections, the photon is considered to be inside,
	 and if at least 5 of the 6 rays have an even amount of intersections, the 
	 photon is considered outside of the object. We allow for one "mismatch" number
	 of intersection, because a ray might touch the surface, returning an ambiguous/wrong
	 answer.
	 @param			inPoint
	 @result		Boolean wheteher or not the point is inside
	 */

    virtual bool	IsInsideObject(RealV& inPoint);

    /*!
     @function      GetNumberOfIntersectionsWithSurface
	 @discussion	Given the input photon, this function will determine the number of intersections that the ray of
	 				the photon has with the surface of the object.
	 @param			ioPhoton
	 */
    
    virtual size_t		GetNumberOfIntersectionsWithSurface(Photon *ioPhoton);
    
	/*!
     @function		IsGeometryConsistent
	 @discussion	The function makes multiple checks to ensure everything was entered correctly or coded correctly.
                    It verifies the indices at interfaces (inside, outside) and other parameters.
	 @result        If the geometry is consistent, return true.  If it returns false, there might be a problem in the code.
	 */
    
    bool 	IsGeometryConsistent();
        
	/*!
     @function      GetWorld
	 @discussion	Obtain the highest object that contains this object (i.e. world)
	 */
    
    MCWorld* GetWorld();

    /*!
     @function     	SetWorld
	 @discussion	Set the world that contains this object.
	 */
    
    void SetWorld(MCWorld* inWorld);
    
	/*!
     @function      GetRandomPointInsideObjectUniformDistribution	
	 @discussion	This function is under development. It will be used for obtaining a point uniformly distributed inside object.
	 */
    
    virtual RealV   GetRandomPointInsideObjectUniformDistribution();
    
	/*!
     @function      GetRandomPointInsideObjectEnergyDistribution
	 @discussion	This function is under development. It will be used for obtaining a point distributed according
                    to energy distribution inside object.
	 */
    
    virtual RealV   GetRandomPointInsideObjectEnergyDistribution();
                    
	/*!
     @function      GetBoundingBox
	 @discussion	This function is under development. It returns the two opposite corners of a box
                    that contains the object.
	 */
    
    virtual void   GetBoundingBox(RealV& outTopLeft, RealV& outBottomRight);
        
    /*!
     @function		DecreaseWeight
     @discussion	This decreases the weight of the photon after a scattering event
     @param			ioPhoton
     */
        
    virtual void	DecreaseWeight(Photon *ioPhoton);
    
    
	/*!
     @function		MoveWithinObject
	 @discussion	The photon is displaced by a given distance in the object (between scattering events). 
                    If the medium is optically active, the photon polarization is rotated.
	 @param			ioPhoton The photon to move
     @param         inDist The distance by which the photon moves
	 */
    
    virtual void	MoveWithinObject(Photon *ioPhoton, double inDist);
    
	/*!
     @function		ScatterInObject
	 @discussion	This transforms the photon upon scattering: the Stokes vector of the photon is multiplied
                    by the Mueller matrix that represents the scatterer.  The Stokes vector
                    is normalized after the multiplication.
	 @param			ioPhoton The photon
     @param         inTheta The scattering angle 
     @param         inPhi   The transverse angle by which the current e_parallel needs to be rotated
                    in order to be in the scattering plane
	 */
    
    virtual void	ScatterInObject(Photon *ioPhoton, double inTheta, double inPhi);
    

	void			SetDetectionStatisticsOnSurfaceElements(unsigned long inDetection);

	/*!
		@function		SetAcceptBallisticPhotons
	 @discussion	Determines whether or not ballistic photons are counted.  For objects, the default is true, for detectors it's false.
	 @param			ballistic True if the surfaces count ballistic photons 
	 */
    
     void	SetAcceptBallisticPhotons(bool ballistic);

	/*!
		@function		ScoreOnSurface
	 @discussion	All statistics are scored here after interface crossing.
	 @param			inPhoton The photon 
     @param         inIntersectElement The properties of the intersect element.
	 */
    
    virtual void	ScoreOnSurface(Photon *inPhoton, IntersectElement& inIntersectElement);
	
	/*!
     @function		ScoreInVolume
	 @discussion	All statistics are scored here after energy is deposited in the object.
                    Very important: the statstics are stored "globally" (in MCWorld) not in each
                    object.
	 @param			ioPhoton    The photon
     @param         inEnergyDeposited   The amount of energy that was deposited (in units of photons)
	 */
    
//    virtual void	ScoreInVolume(Photon *ioPhoton, double inEnergyDeposited);
	virtual void	ScoreInVolume(Photon *ioPhoton, double inEnergyDeposited, bool pastHistory = false);

	
	/*!
     @function		DumpInterfaceContructionParametersToStream
	 @discussion	This is part of the "Dump" function series to save everything to stream (output or file).
                    All the parameters related to given interface (normal, a, b vectors, etc...) are output
                    to stream
	 @param			out The output stream
     @param			interface The interface
	 */
    
    virtual void	DumpInterfaceContructionParametersToStream(ostream & out, long interface);
    
	/*!
     @function		DumpEnergyDepositionToStream
	 @discussion	This is part of the "Dump" function series to save everything to stream (output or file).
                    Energy deposited is saved to stream.  Note that only MCWorld will call this function.
	 @param			out The output stream
	 */
    
    virtual void    DumpEnergyDepositionToStream(ostream & out );

	/*!
		@function		DumpFluenceToStream
	 @discussion	This is part of the "Dump" function series to save everything to stream (output or file).
	 The fluence is saved to stream.  Note that only MCWorld will call this function.
	 @param			out The output stream
	 */
    
    virtual void    DumpFluenceToStream(ostream & out );
	
	/*!
     @function		DumpStokesVToStream
	 @discussion	This is part of the "Dump" function series to save everything to stream (output or file).
                    It saves the Stokes vector information accumulated for the photons that have crossed that interface.
	 @param			out Output strean
     @param			interface The interface
	 */
    
    virtual void	DumpStokesVToStream(ostream & out, long interface);
    
	/*!
     @function		DumpOtherStatsToStream
	 @discussion	This is part of the "Dump" function series to save everything to stream (output or file).
                    It outputs to stream the other stats for the photons that have crossed that interface and 
                    that are not Stokes-vector related (path length for instance).
	 @param			out The output stream
     @param			interface The interface
	 */
    
    virtual void	DumpOtherStatsToStream(ostream & out, long interface);
    
	/*!
     @function		DumpProximityListStatistics
	 @discussion	Dump the proximity list statistics
	 @param			out The output stream
	 */

    virtual void    DumpProximityListStatistics(ostream & out);


	/*!
		@function		DumpProximityListStatisticsFor
	 @discussion	Dump the proximity list statistics for a given ProximityList
	 @param			out The output stream
	 @param			proximityList The proximity list
	 */
	
	void DumpProximityListStatisticsFor(ostream & out, map<long,ProximityBox*>& proximityList);

	friend std::ostream& operator<<(std::ostream& os, const MCObject& obj)
	{
	  os<<obj.GetName()<<"("<<obj.mObjectID<<")";
	  return os;
	}

	std::string id_str(void)const
	{
	  std::ostringstream oss;
	  oss<<*this;
	  return oss.str();
	}


protected:

	/*! @var		gInterfaceCrossingAlgorithm			Using the enum's above, determine the algorithm used for interface crossing (kd-tree, proximity boxes, all). */
	static int gInterfaceCrossingAlgorithm;

	/*! @var    *mTree           A KD-tree object holding a spatially optimized search structure */
    MCKDTree *mTree; 
    
    /*! @var 	mObjectID	A unique object ID that identifies this object unambiguously */
    long mObjectID;
    /*! @var 	mObjectName	A name for this object (optional).  For readability in output file */
    string mObjectName;
    /*! @var		dict		Deprecated. Dictionary (associative array) that contains the construction parameters (saved so it can be saved to disk). */
    map<string,string> dict;
    /*! @var		rotationPerCmInClearSpace	Optical activity given in rotation per centimeter in clear media */
    double rotationPerCmInClearSpace;
    /*! @var		mIndexOutside				Index of refraction outside the object */
    double mIndexOutside;
    /*! @var		N_photon					Number of photons that have been launched in this object */
    long N_photon;
    /*! @var		gNormalizationIntensity     Normalization factor for transmittance, reflectance, absorbance */
    static double gNormalizationIntensity;
	/*! @var        mAbsorbedPhotons            Amount of energy (in photon weight) that has been absorbed in this object */
    double mAbsorbedPhotons;

    /*! @var		useNextEventEstimator		Whether or not we use the next event estimator in PropagateInObject() */
    bool useNextEventEstimator;
    
	/*! @var 		mXRotation		By how much the object has been rotated around the X axis (lab frame) */
	double mXRotation;
	/*! @var 		mYRotation		By how much the object has been rotated around the Y axis (lab frame) */
	double mYRotation;
	/*! @var 		mZRotation		By how much the object has been rotated around the Z axis (lab frame) */
	double mZRotation;

	/*! @var mEnergy		Three dimensional array pointer for storing the energy that gets deposited in all objects */
    double ***mEnergy;
	/*! @var mFluence		Three dimensional array pointer for storing the fluence that flows through all objects */
    double *mFluence;

	/*! @var Vol_Nx			Number of bins in x coordinates for mEnergy */
    long Vol_Nx;
	/*! @var Vol_Ny			Number of bins in y coordinates for mEnergy */
    long Vol_Ny ;
	/*! @var Vol_Nz			Number of bins in z coordinates for mEnergy */
    long Vol_Nz;

	/*! @var Vol_XMin		Mininum x value (global coordinates) for tabulating statistics for mEnergy */
    double Vol_XMin;
	/*! @var Vol_XMax		Maximum x value (global coordinates) for tabulating statistics for mEnergy */
    double Vol_XMax;
	/*! @var Vol_YMin		Mininum y value (global coordinates) for tabulating statistics for mEnergy */
    double Vol_YMin;
	/*! @var Vol_YMax		Maximum y value (global coordinates) for tabulating statistics for mEnergy */
    double Vol_YMax;
	/*! @var Vol_ZMin		Mininum z value (global coordinates) for tabulating statistics for mEnergy  */
    double Vol_ZMin;
	/*! @var Vol_ZMax		Maximum z value (global coordinates) for tabulating statistics for mEnergy */
    double Vol_ZMax ;
	/*! @var gObjectCounter	Global counter for assigning unique object ID */
    static long gObjectCounter;
	/*! @var mBallisticPhotons	Variable to determine whether of not ballistic photons are scored.
								Typically, objects will consider ballistic photons, but detectors won't. */
	bool mBallisticPhotons;

	/*! @var surfaceElements	Vector of surface elements that define this object, as well as any interface 
                                to any other object */
    vector<SurfaceElement*> mSurfaceElements;

    long mNumOfInterfaces;
    
	/*! @var	mRandomScatterer	Pointer to object responsible for dealing with any scattering-related issue */
    MCRandomScatterer* mRandomScatterer;
    
    /*! @var mGlobalOrigin	Origin of object in global coordinates */
    RealV	mGlobalOrigin;

    /*! @var mWorld   This is a pointer to the top object (an MCWorld). */
    MCWorld* mWorld;

    /*! @var mVertexList    List of all vertices, in global coordinates */
    vector<RealV> mVertexList;

    /*! @var mNormalList    List of all normals, in global coordinates */
    vector<RealV> mNormalList;

    /*! @var mVertexList    List of all vertices, in global coordinates */
    vector<double> mIntensityList;

    /*! @var mFaceList    List of faces, where a face is defined as a structure containing
                            a list of Face structs, see definition in MCObject.h. It contains
                            indexes of vertex, normal and texture (texture is currently unused). */
    vector<Face> mFaceList;

    /*! @var gVertexOffset  Offset representing the number of vertices output so far for .obj format (internal, do not use) */
    static long gVertexOffset;

    /*! @var gNormalOffset  Offset representing the number of vertices output so far for .obj (internal, do not use)  */
    static long gNormalOffset;
	
    /*! @var mProximity */
    bool mProximity;
    
    RealV mBoundingBoxBottom;

    RealV mBoundingBoxTop;
    
    map<long,ProximityBox*> mProximityBoxList;
    
    bool cacheIsValid;
    
	
    static long gProximitySubdivisions;
    static long gMinPhotonCount ;
	static long gMinSurfaceElements;
};


#endif

