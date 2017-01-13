/*
 *  CFHelper.h
 *  Pol3d-MC
 *
 *  Created by Daniel Cote on 25/06/09.
 *  Copyright 2009 Universite Laval. All rights reserved.
 *
 */

/*!
    @header		CFHelper
    @abstract   Helper functions to take advantage of the CoreFoundation functions to build the 3D environment for the computation
    @discussion The CoreFoundation library defines types (CFArray, CFDictionary, CFNumber, CFString, CFBoolean) that can be used
				safely in a program but also can be written to and read from files in a format known as PropertylList (.plist) which
				is an XML format. On OS X, this format is the basis for all configuration files in the system and their
				is a great Property List Editor.app that comes with XCode. 
 
				The hierarchical nature of the PropertyList and its available library made it a natural choice to use for a parameter file
				for 3D objects with polarization.
 
				The CoreFoundation library is available in OS X as part of the system, but also with CF-Lite as an open-source
				project that compiles on Unix platforms.  A Linux version exists.
				More information at: http://cafeine.crulrg.ulaval.ca/users/dccote/weblog/0514e/CoreFoundation_Lite_on_Linux.html

				Possible Windows port: http://karaoke.kjams.com/wiki/CFLite

				Even more information at: http://developer.apple.com/opensource/cflite.html
 */

#include <CoreFoundation/CoreFoundation.h>
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


/*!
    @function
    @abstract   This function creates a bare MCWorld object that will eventually contain the objects, the light sources, the detectors and parameters for the output file
    @discussion The MCWorld object contains everything about the computation: the objects and the various parameters.
 
				Here I will describe all parameters possible. [To be done]
 
    @param      worldProperties	The properties from a CFDictionaryRef
    @result     An empty MCWorld
*/

MCWorld* 
CreateBareWorldFromCFDictionary(CFDictionaryRef worldProperties);

/*!
    @function
    @abstract   This takes a CFArray of lightSources which are themselves CFDictionary with all the properties of light sources.
    @discussion A computation can have more than one light source (laser), sometimes outside objects, sometimes inside objects (diffuser, bioluminescence).
				This function will place the lightsources in MCWorld and determine in what object they lie.
				MCWorld keeps track of all lightsources.  See MCWorld for details.
    @param      world	The MCWorld object for the computation that will contain the inlcuded light sources
	@param		sources	The CFArray containing the light sources, as CFDictionaryRef 
*/

void
PopulateWorldWithLightSourcesFromCFDictionary(MCWorld* world, CFArrayRef sources);

/*!
 @function	
 @abstract   This takes a CFArray of objects which are themselves CFDictionary with all the properties of objects and includes them in rootObject, which is in world.
 @discussion A computation has multiple objects in it, of various shapes.  This function populates the MCWorld with them.
			 MCWorld keeps track of all objects inserted, see MCWorld for details.
 @param      world The MCWorld that contains the objects.
 @param		 rootObject The object that will contain the inserted objects
 @param		 objects The CFArrayRef of CFDictionaryRef objects that will be inserted.  
 */

void
PopulateRootObjectWithObjectsFromCFDictionary(MCWorld* world, MCObject* rootObject, CFArrayRef objects);

/*!
    @function
    @abstract   Sets the properties for the calculation and the many variables of MCWorld from the CFDictionaryRef 
    @discussion This avoid the creation of multiple global variables and keeps everything neat and tidy within MCWorld
	@param      world The MCWorld for which we want to set the variables
	@param      calculationDict The CFDictionary that contains the variables:
				statistics     [ CFArray  ] : an array of Boolean 
				outputFilename [ CFString ] : The variable that will be used for filename
				keepAllPhotonStats [ CFBoolean ] : Whether or not we keep the stats for all photons
				debugOutputPhotonStats [ CFBooleab ] : Whether or not we output all the stats for all photons 
				debugLevel [ CFNumber ] : An integer from 0 (quiet) to 4 (Extremely verbose) used for debugging.
*/

void
SetStatisticsForCalculation(MCWorld* world, CFDictionaryRef calculationDict);

/*!
 @function	
 @abstract   Creates a lightsource from a CFDictionaryRef of properties
 @discussion The software can handle multiple light sources: lasers, isotropic source, diffuser, etc... They all inherit from MCSource.
 @param      lightSourceProperties The light source properties in a CFDictionaryRef
 @result     The light source, of the right type
 */

MCSource* 
CreateLightSourceFromCFDictionary(CFDictionaryRef lightSourceProperties);

/*!
 @function	
 @abstract   Creates an object to be placed in the world from a CFDictionaryRef of properties
 @discussion The software can handle multiple object types: boxes, ellipsoids, detectors, cylinders, layers, generic objects, etc... They all inherit from MCObject.
 @param      objectProperties The object properties in a CFDictionaryRef
 @result     The object, of the right type
 */

MCObject* 
CreateObjectFromCFDictionary(CFDictionaryRef objectProperties);

/*!
 @function	
 @abstract   Creates a class that represents the scattering and optical properties of an object from a CFDictionaryRef
 @discussion An MCObject can be scattering, absorbing, clear, etc...  All the optical and scattering properties
			 of an object are managed by MCRandomScatterer.  This allows the possibility of using different computation models
			 for random angle sampling (e.g., different publications), various scattering types (e.g., intensity-based or polarization)
			 or evem look-up tables for the properties.
 @param      materialProperties The material properties in a CFDictionaryRef
 @result     The random scatterer object
 */

MCRandomScatterer* 
CreateRandomScattererFromCFDictionary(CFDictionaryRef materialProperties);

/*!
 @function	
 @abstract   Obtains a long list of statistical parameters that may or may not be included in the output file
 @discussion See MCWorld.h for details on the value of each key.
 @param      statisticsDict The CFDictionaryRef with all CFBoolean values representing each statistics
 @result     An unsigned long integer where each bit represents one statistics (true or false)
 */

unsigned long 
GetStatisticsBooleanSettingsFromCFDictionary(CFDictionaryRef statisticsDict);

/*!
 @function	
 @abstract   Helper function to read a boolean value from a dictionary with a default value and error checking
 @discussion If the value is not of the right type, this code will throw an error
 @param      dict   The CFDictionaryRef possibly containing the key
 @param		 key	The value of the key as a CFStringRef
 @result     true or false corresponding to the value in the dictionary
 */

bool
GetBooleanInDictWithDefaultFalse(CFDictionaryRef dict, CFStringRef key);

/*!
 @function	
 @abstract   Helper function to read a CFNumber value from a dictionary with a default value and error checking
 @discussion <#(description)#>
 @param      dict   The CFDictionaryRef possibly containing the key
 @param		 key	The value of the key as a CFStringRef
 @param		 type	The type of integer to which CFNumber must be converted.  This must correspond to the C++ type of value
 @param	     value	The integer, cast to void, to which the value is copied
 */

void	
GetNumberInDictWithDefault(CFDictionaryRef dict, CFStringRef key, CFNumberType type, void* outValue);

/*!
 @function	
 @abstract   Helper function to read an integer value storaed as a CFNumber in a dictionary with error checking
 @discussion Helper function to make number conversion easier and simpler to read.
 @param      dict   The CFDictionaryRef possibly containing the key
 @param		 key	The value of the key as a CFStringRef
 @result	 The integer obtained from CFNumber as a 32-bit signed value
 */

long
GetIntegerFromDict(CFDictionaryRef dict, CFStringRef key);

/*!
 @function	
 @abstract   Helper function to read a double value storaed as a CFNumber in a dictionary with error checking
 @discussion Helper function to make number conversion easier and simpler to read.
 @param      dict   The CFDictionaryRef possibly containing the key
 @param		 key	The value of the key as a CFStringRef
 @result	 The value obtained from CFNumber as a double float value
 */

double
GetDoubleFromDict(CFDictionaryRef dict, CFStringRef key);

/*!
    @function
    @abstract   Obtains a string from a disctionary with error checking
    @discussion <#(description)#>
	@param      dict   The CFDictionaryRef possibly containing the key
	@param		 key	The value of the key as a CFStringRef
    @result     A C++ string object
*/

string
GetStringFromDict(CFDictionaryRef dict, CFStringRef key);

/*!
    @function	
    @abstract   Converts a CFString to a RealV
    @discussion This is a helper function to safely convert of 3D vector from a string in CFStringRef format
    @param      theString A CFStringRef in the format (number,number,number)
    @result     A RealV
*/

RealV 
CFStringToRealV(CFStringRef theString);

#define CFStringToCPtr(cf, cs)  { char Str[255]; CFStringGetCString(cf, Str, 256, kCFStringEncodingISOLatin1); c = Str; }		

string
CFDescriptionToString(CFTypeRef object);

void
SetGlobalCFStringEncoding(long encoding);


string
CFStringToString(CFStringRef cf);

/*!
    @function
    @abstract   Sets the global string encoding for converting CFStrings to char*
    @discussion This simplifies the code by using a global variable to set the encoding. Typically, you would use kCFStringEncodingISOLatin1
    @param      encoding The encoding such as kCFStringEncodingISOLatin1, kCFStringEncodingWindowsLatin1, etc... as defined in CFString.h
*/

void
SetGlobalCFStringEncoding(long encoding);


extern long gStringEncoding;
