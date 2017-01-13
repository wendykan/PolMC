#ifndef __MCWORLD
#define __MCWORLD

#include "MCObject.h"

class MCSource;

/*
	@enum Information that is output to disk
	
	Instead of saving absolutely everything to disk, we offer
	a simple way to limit the size of files by keeping one global variable
	where each bit indicates whether or not we want that information 
	on the output file.  
 
 */

enum {
	kEnergyDistributionRawBinary = 1 << 0,
	kEnergyDistributionXYZUnits = 1 << 1,
	kEnergyDistributionBinUnits = 1 << 2,
	kGeometryWorldWavefrontFormat = 1 << 3,
	kGeometrySeparateWavefrontFormat = 1 << 4,
	kGeometryXMLFormat = 1 << 5,
	kGeometryMatlabFormat = 1 << 6,
	kGeometryMathematicaFormat = 1 << 7,
	kObjectPhotonScatteringStats = 1 << 8,
	kObjectInterfaceConstructionParameters = 1 << 9,
	kObjectTotalAbsorbance = 1 << 10,
	kObjectInterfaceTotalTransmittance = 1 << 11,
	kObjectInterfaceStokesVectorI = 1 << 12,
	kObjectInterfaceIntensity =		kObjectInterfaceStokesVectorI,
	kObjectInterfaceStokesVectorQ = 1 << 13,
	kObjectInterfaceStokesVectorU = 1 << 14,
	kObjectInterfaceStokesVectorV = 1 << 15,
	kObjectInterfaceStokesVectorBetaLinear = 1 << 16,
	kObjectInterfaceStokesVectorBetaCircular = 1 << 17,
	kObjectInterfaceStokesVectorAny = kObjectInterfaceStokesVectorI | kObjectInterfaceStokesVectorQ | kObjectInterfaceStokesVectorU \
	| kObjectInterfaceStokesVectorV | kObjectInterfaceStokesVectorBetaLinear | kObjectInterfaceStokesVectorBetaCircular,
	kObjectInterfaceAveragePathLength = 1 << 18,
	kObjectInterfaceAveragePolarizedPathLength = 1 << 19,
	kObjectInterfaceAverageNumberOfScatteringEvents = 1 << 20,
	kObjectInterfaceAverageAny = kObjectInterfaceAveragePathLength | kObjectInterfaceAveragePolarizedPathLength | kObjectInterfaceAverageNumberOfScatteringEvents,
	kScattererStatistics1D = 1 << 21,
	kScattererStatistics2D = 1 << 22,
	kFluenceDistributionRawBinary = 1 << 23,
	kFluenceDistributionXYZUnits = 1 << 24,
	kFluenceDistributionBinUnits = 1 << 25,
	kObjectInterfaceTimeResolvedIntensity = 1 << 26
	
};

enum {
	kAllAlgorithm = 1,
	kProximityBoxesAlgorithm = 2,
	kKDTreeAlgorithm = 3
};

#define kDefaultOutput  (kObjectTotalAbsorbance | kObjectInterfaceTotalTransmittance | kGeometrySeparateWavefrontFormat )

class MCWorld : public MCObject {

public:
    MCWorld( double inImgXMin,
             double inImgXMax,
             double inImgYMin,
             double inImgYMax,
             double inImgZMin,
             double inImgZMax,
             long inPtsWidth,
             long inPtsHeight,
             long inPtsDepth);
    MCWorld();

    ~MCWorld();
    
    virtual IntersectElement	PropagateInObject(Photon *ioPhoton, double inDist = 0);
    virtual long	PathCrossesInterface(Photon *ioPhoton, double inDist, IntersectElement& outIntersectElement, long inCrossingEvents);
    virtual bool	IsOutsideObject(RealV& inLocalPosition); 

    long ReadWorldFromFile(string inFilename, double inObjectResolution);

    /*!
		@function Place
	 @discussion Places an object into this object.  The surfaces must not cross.  If they do, you will get unexpected results.
                 You must have set the exterior of this current object prior to calling this.
	 @param      inObject	The object being placed into the current one
	 @param		 inLocalPosition The local position where the object is placed.  The origin depends on the subclass.
	 */
	
    virtual void   Place(MCObject* inObject, RealV inLocalPosition);
    

    /*!
		@function PlaceLightSourceInWorld
	 @discussion Places a light source in the world. The function will try to figure out inside what object the source is.
				 If it can't figure it out, it will throw an error.
	 @param		 inSource	The source being placed. It sposition is already set upon creation.
	 */
    virtual void PlaceLightSourceInWorld(MCSource* inSource);

    /*!
		@function GetSourceList
	 @discussion This function returns the list of sources in the world.  Should be considered read-only.
	 @result	 A vector<> of MCSource, should be considered read-only.
	 */
	
	vector<MCSource*> GetSourceList();


    /*!
		@function GetObjectList
	 @discussion This function returns the list of objects in the world.  Should be considered read-only.
	 @result	 A vector<> of MCObject, should be considered read-only.
	 */
	
	vector<MCObject*> GetObjectList();

		/*!
		@function DumpToFile
	 @discussion This function will save to file in an XML format all the information about this object, including all photon stats.
     The function will call other functions (all starting with Dump...()).
	 @param      inFilename	The filename of the file that is produced
	 */
    
    void	DumpToFile(string inFilename);

	/*!
		@function		IsGeometryConsistent
	 @discussion	The function makes multiple checks to ensure everything was entered correctly or coded correctly.
	 It verifies the indices at interfaces (inside, outside) and other parameters. The world implementation
	 simply calls the function for each object inside of it.
	 @result        If the geometry is consistent, return true.  If it returns false, there might be a problem in the code.
	 */
    
    bool 	IsGeometryConsistent();

	/*!
	 @function		IsInsideWhichObject
	 @discussion	The function asks all the objects in MCWorld wheteher or not this point is inside them in order 
					to figure out where the point is.  The function calls IsInsideObject() for all objects.
	 @result        The object in which the point is, only if unambiguous.
	 */
	
	MCObject*	IsInsideWhichObject(RealV inPoint);

    /*!
        @function GetSurfaceElementsCloseToSegment
	 @discussion Get an array of all surface elements close to line segment
	 @param      outSe	A reference to an array pointer that will contain the surface elements.
	 @result        the number of elements.  Indexing the array goes from 0 to (numOfElements-1).
	 */
	
//    virtual long GetSurfaceElementsCloseToSegment(vector<SurfaceElement*>& outSe, RealV inStart, RealV inEnd);
    
    MCObject* FindObjectByName(string inName);

    long AddVertex(RealV inPoint);

    RealV GetVertex(long inWhich);

    long AddSurfaceElement(SurfaceElement* inSurfaceElement);

    SurfaceElement* GetSurfaceElement(long inWhich);

    /*!
        @function		DumpProximityListStatistics
	 @discussion	Dump the proximity list statistics
	 @param			out The output stream
	 */
    
    virtual void    DumpProximityListStatistics(ostream & out);
   
	IntersectElement RayTraceInObject(Photon *ioPhoton, double inDist);
 
    long    LoadMaterialPropertiesTable(string inFilename);

	/*!
	 @function		SetOutputProperties
	 @discussion	Sets the output properties that will be output to disk.
	 @param			inValue The value using the enums (kObjectAbsorbance, etc...)
	 
	 */
	static	void	SetOutputProperties(unsigned long inValue);
	static	unsigned long	GetOutputProperties();

	
	/*!
    @function     SetKeepPhotonStatistics
    @abstract	  For all photons launched, determines whether travel statistics are kept or not.  When activated, this slows down the computation but is sometimes needed for calculation (time-resolved measurement for instance) or debugging.
    @discussion 
	 */

	void SetKeepPhotonStatistics(bool inValue);

	/*!
	 @function     KeepPhotonStatistics
	 @abstract	   Query if travel statistics are kept or not for all photons.
	 @discussion 
	 */
	
	bool KeepPhotonStatistics();

	/*!
	 @function    SetDebugOutputPhotonStatistics
	 @abstract	  Determines whether or not statistics are dumped to screen when a photon dies.  Very verbose and very slow, debug purposes only.
	 @discussion 
	 */
	
	void SetDebugOutputPhotonStatistics(bool inValue);
	
	/*!
	 @function     DebugOutputPhotonStatistics
	 @abstract	   Query if travel statistics are output to screen or not when a photon dies.
	 @discussion 
	 */
	
	bool DebugOutputPhotonStatistics();
	
	/*!
	 @function     OutputFilename
	 @abstract	   The filename of the output file.
	 @discussion 
	 */
	
	string OutputFilename();

	/*!
	 @function     SetOutputFilename
	 @abstract	   Sets the filename of the output file
	 @discussion 
	 */
	
	void SetOutputFilename(string inOutputFilename);

	/*!
	 @function     SetOutputFilenameGeometry
	 @abstract	   Sets the filename of the output file for the geometry
	 @discussion 
	 */
	
	void SetOutputFilenameGeometry(string inOutputFilename);

	/*!
	 @function     OutputFilenameGeometry
	 @abstract	   The filename of the output file for geometry
	 @discussion 
	 */
	
	string OutputFilenameGeometry();

	/*!
	 @function     SetOutputFilenameGeometryFormat
	 @abstract	   Sets the format of the output file for the geometry
	 @discussion 
	 */
	
	void SetOutputFilenameGeometryFormat(long inFormat);
	
	/*!
	 @function     OutputFilenameGeometryFormat
	 @abstract	   The format of the output file for geometry
	 @discussion 
	 */
	
	long OutputFilenameGeometryFormat();
	
	/*!
	 @function     SetOutputFilenameBinaryGeometryWithIntensity
	 @abstract	   Sets the filename of the output file for the binary geometry with intensity
	 @discussion 
	 */
	
	void SetOutputFilenameGeometryWithIntensity(string inOutputFilename);
	
	/*!
	 @function     OutputFilenameBinaryGeometryWithIntensity
	 @abstract	   The filename of the output file for binary geometry with intensity
	 @discussion 
	 */
	
	string OutputFilenameGeometryWithIntensity();
	
	/*!
	 @function     SetOutputFilenamePhotonPaths
	 @abstract	   Sets the filename of the output file for binary photon paths
	 @discussion 
	 */
	
	void SetOutputFilenamePhotonPaths(string inOutputFilename);

	/*!
	 @function     OutputFilenamePhotonPaths
	 @abstract	   The filename of the output file for binary photon paths
	 @discussion 
	 */
	
	string OutputFilenamePhotonPaths();

	//added Wendy Kan
	void SetOutputPhotonPathsFormat(long inOutputFileFormat);
	long OutputFilenamePhotonPathsFormat();
	
	// added by Wendy Kan 7/8/10
	void SetmaximumDepthTargetObjectName(string inObjectName);
	string GetmaximumDepthTargetObjectName();


	
	
protected:
	/*! @var mObjects This is a vector of pointers to objects that are contained within this object */
	vector<MCObject*>  mObjects;

	vector<SurfaceElement*> mAllSurfaceElements;

	vector<MCObject*> mAllObjects;

	vector<MCSource*> mAllSources;

	/*! @var		gOutputProperties			Using the enum's above, we limit the number of parameters output to disk. */
	static unsigned long gOutputProperties;

	/*! @var		keepPhotonStats				When photons are launched, activates statistics tracking */
	bool keepPhotonStats;

	/*! @var		debugOutputPhotonStats		When photons are dead, will dump all travel statistics to the screen.  Very verbose and very slow.  */
	bool debugOutputPhotonStats;

	/*! @var		outputFilename				The filename that could be used for the output file */
	string outputFilename;

	/*! @var		outputFilenameGeometry				The filename that could be used for the output file geometry */
	string outputFilenameGeometry;

	/*! @var		outputFilenameGeometryFormat				The format that will be used for the output file geometry */
	long outputFilenameGeometryFormat;

	/*! @var		outputFilenameGeometryWithIntensity	The filename that could be used for the output file  geometry with intensity*/
	string outputFilenameGeometryWithIntensity;

	/*! @var		outputFilenamePhotonPaths				The filename that could be used for the output file  photon paths */
	string outputFilenamePhotonPaths;
	
	/* added by Wendy Kan*/
	/*! @var		outputPhotonPathsFormat				The format that could be used for the output file  photon paths */	
	long outputPhotonPathsFormat;

	// added by Wendy Kan 7/8/10
	string maximumDepthTargetObjectName;
	
public:
		/*! @var		speedOfLight				Speed of light in appropriate units for calculation */
		static double	speedOfLight;
};

#endif

