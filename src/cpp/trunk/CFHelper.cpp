/*
 *  CFHelper.cpp
 *  Pol3d-MC
 *
 *  Created by Daniel Cote on 25/06/09.
 *  Copyright 2009 Universite Laval. All rights reserved.
 *
 */

#include "CFHelper.h"

long gCFStringEncoding = kCFStringEncodingISOLatin1;

#pragma mark Helper CF functions 

string
CFDescriptionToString(CFTypeRef object) {
    UInt8 buffer[10000];
    CFIndex got;
    CFStringRef description = CFCopyDescription(object);
    CFStringGetBytes(description,
					 CFRangeMake(0, CFStringGetLength(description)),
					 gCFStringEncoding, '?', TRUE, buffer, 10000, &got);
    buffer[got] = (char)0;
    CFRelease(description);

	return string((char*)buffer);
}

void
SetGlobalCFStringEncoding(long encoding)
{
	gCFStringEncoding = encoding;
}

long
GetIntegerFromDict(CFDictionaryRef dict, CFStringRef key)
{
	CFNumberRef theNumber;
	long value;
	
	if ( CFDictionaryGetValueIfPresent(dict, key, (const void**)&theNumber) ) {
		if ( CFGetTypeID(theNumber) == CFNumberGetTypeID() ) {
			CFNumberGetValue(theNumber, kCFNumberLongType, (void*)&value);
		} else if ( CFGetTypeID(theNumber) == CFStringGetTypeID()  ) {
			clog << RightNow() << "Warning: the key " << key << " in dictionary is of type String and should be Number\n";
			value = CFStringGetIntValue((CFStringRef)theNumber);
		} else {
			ThrowRuntimeError("Number expected for key " << CFStringToString(key) << " in dictionary:\n" << CFDescriptionToString(dict) );
		}
	} else {
		ThrowRuntimeError("Number key " << CFStringToString(key) << " is not present in dictionary:\n" << CFDescriptionToString(dict));
	}
	
	return value;
}

double
GetDoubleFromDict(CFDictionaryRef dict, CFStringRef key)
{
	CFNumberRef theNumber;
	double value = 0;
	
	if ( CFDictionaryGetValueIfPresent(dict, key, (const void**)&theNumber) ) {
		if ( CFGetTypeID(theNumber) == CFNumberGetTypeID() ) {
			CFNumberGetValue(theNumber, kCFNumberDoubleType, (void*)&value);
		} else if ( CFGetTypeID(theNumber) == CFStringGetTypeID()  ) {
			clog << RightNow() << "Warning: the key " << key << " in dictionary is of type String and should be Number\n";
			value = CFStringGetDoubleValue((CFStringRef)theNumber);
		} else {
			ThrowRuntimeError("Number expected for key " << CFStringToString(key) << " in dictionary:\n" << CFDescriptionToString(dict) );
		}
	} else {
		ThrowRuntimeError("Number key " << CFStringToString(key) << " is not present in dictionary:\n" << CFDescriptionToString(dict));
	}
	
	return value;
}

double
GetDoubleFromDictWithDefault(CFDictionaryRef dict, CFStringRef key, double inDefault)
{
	double value = inDefault;

	try { 
		return GetDoubleFromDict(dict, key);
	} catch (...) {
		;
	}
	
	return value;
}

void
GetNumberInDict(CFDictionaryRef dict, CFStringRef key, CFNumberType type, void* value)
{
	CFNumberRef theNumber;
	
	if ( CFDictionaryGetValueIfPresent(dict, key, (const void**)&theNumber) ) {
		if ( CFGetTypeID(theNumber) == CFNumberGetTypeID() ) {
			CFNumberGetValue(theNumber, type, value);
		} else {
			ThrowRuntimeError("Number expected for key " << CFStringToString(key) << " in dictionary:\n" << CFDescriptionToString(dict) );
		}
	} else {
		ThrowRuntimeError("Number key " << CFStringToString(key) << " is not present in dictionary:\n" << CFDescriptionToString(dict));
	}
}

string
GetStringFromDict(CFDictionaryRef dict, CFStringRef key)
{
	CFStringRef theString;
	
	if ( CFDictionaryGetValueIfPresent(dict, key, (const void**)&theString) ) {
		return CFStringToString(theString);
	}
	
	return string("");
}

void
GetNumberInDictWithDefault(CFDictionaryRef dict, CFStringRef key, CFNumberType type, void* value)
{
	try { 
		GetNumberInDict(dict, key, type, value);
	} catch (...) {
		;
	}
}

bool
GetBooleanInDictWithDefaultFalse(CFDictionaryRef dict, CFStringRef key)
{
	CFBooleanRef aBoolean;
	if ( CFDictionaryGetValueIfPresent(dict, key, (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			return true;
		} else {
			return false;
		}
	}
	
	return false;
}


#pragma mark Objects reading
MCWorld* 
CreateBareWorldFromCFDictionary(CFDictionaryRef worldProperties)
{
	
	long Nx, Ny, Nz;
	double xMin, xMax, yMin, yMax, zMin, zMax;
	
	// If any of the properties are missing for the size, we simply return a bare world.
	
	try {

		xMin = GetDoubleFromDict(worldProperties, CFSTR("xMin"));
		xMax = GetDoubleFromDict(worldProperties, CFSTR("xMax"));
		yMin = GetDoubleFromDict(worldProperties, CFSTR("yMin"));
		yMax = GetDoubleFromDict(worldProperties, CFSTR("yMax"));
		zMin = GetDoubleFromDict(worldProperties, CFSTR("zMin"));
		zMax = GetDoubleFromDict(worldProperties, CFSTR("zMax"));
		
		Nx = GetIntegerFromDict(worldProperties, CFSTR("Nx"));
		Ny = GetIntegerFromDict(worldProperties, CFSTR("Ny"));
		Nz = GetIntegerFromDict(worldProperties, CFSTR("Nz"));
	} catch (...) {
		return new MCWorld();
	}
	
	MCWorld* theWorld = new MCWorld(xMin, xMax, yMin, yMax, zMin, zMax, Nx, Ny, Nz);
	
	return theWorld;
}

string 
PropertiesForWorld()
{
	char *theHelp = "Properties for MCWorld:\n\n\
	For volume stats (energy and fluence deposition) in tissue:\n\
	xMin [double] : The minimum value on x axis for which statistics will be tabulated for energy and fluence deposition\n\
	xMax [double] : The maximum value on x axis for which statistics will be tabulated for energy and fluence deposition\n\
	yMin [double] : The minimum value on y axis for which statistics will be tabulated for energy and fluence deposition\n\
	yMax [double] : The maximum value on y axis for which statistics will be tabulated for energy and fluence deposition\n\
	zMin [double] : The minimum value on z axis for which statistics will be tabulated for energy and fluence deposition\n\
	zMax [double] : The minimum value on z axis for which statistics will be tabulated for energy and fluence deposition\n\
	Nx   [int]    : The number of bins between [xMin,xMax] into which the statistics will be tabulated. The bin size is [xMax-xMin]/Nx\n\
	Ny   [int]    : The number of bins between [yMin,yMax] into which the statistics will be tabulated. The bin size is [yMax-yMin]/Ny\n\
	Nz   [int]    : The number of bins between [zMin,zMax] into which the statistics will be tabulated. The bin size is [zMax-zMin]/Nz\n";
	
	return string(theHelp);
}

MCObject* 
CreateObjectFromCFDictionary(CFDictionaryRef objectProperties)
{
	CFDictionaryRef aDict;
	CFStringRef aString, typeString;
	
	MCObject* theObject;
	if ( CFDictionaryGetValueIfPresent(objectProperties, CFSTR("type"), (const void**)&typeString) ) {
		
		// Begin common elements
		string name;
		long statisticsType = kNone, acceptanceCosineElements = 1;
		bool acceptsBallisticPhotons = FALSE;
		
		if ( CFDictionaryGetValueIfPresent(objectProperties, CFSTR("name"), (const void**)&aString) ) {
			name = CFStringToString(aString);
		} else {
			clog << RightNow() << "Warning: No name specified for object\n" ;
		}
		
		double tiltX = GetDoubleFromDictWithDefault(objectProperties, CFSTR("tiltX"), 0);
		double tiltY = GetDoubleFromDictWithDefault(objectProperties, CFSTR("tiltY"), 0);
		double tiltZ = GetDoubleFromDictWithDefault(objectProperties, CFSTR("tiltZ"), 0);
		
		if ( CFStringCompare (typeString, CFSTR("detector"), kCFCompareCaseInsensitive) == kCFCompareEqualTo ) {
			// Detector is unique: no optical properties
			double width;
			long Nx;
			
			GetNumberInDict(objectProperties, CFSTR("width"), kCFNumberDoubleType, &width);
			GetNumberInDict(objectProperties, CFSTR("Nx"), kCFNumberLongType, &Nx);
			bool thereIsAPolarizer = GetBooleanInDictWithDefaultFalse(objectProperties, CFSTR("polarizer"));
			int polarizeAxis = GetIntegerFromDict(objectProperties, CFSTR("polarizeAxis"));
			
			CFDictionaryRef statisticsDict;
			statisticsType = kNone;
			
			if ( CFDictionaryGetValueIfPresent(objectProperties, CFSTR("statistics"), (const void**)&statisticsDict) ) {
				statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kStokesVector")) ? kStokesVector : 0;
				statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kIntensityOnly")) ? kIntensityOnly : 0;
				statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kEntering")) ? kEntering : 0;
				statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kLeaving")) ? kLeaving : 0;
				statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kTimeResolved")) ? kTimeResolved : 0;
				acceptsBallisticPhotons = GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kAcceptsBallisticPhotons"));
				
				acceptanceCosineElements = 1;
				GetNumberInDictWithDefault(statisticsDict,CFSTR("acceptanceCosineElements"), kCFNumberLongType, &acceptanceCosineElements);
			}
			
			
			theObject = new MCDetector(width, Nx, statisticsType, acceptanceCosineElements,thereIsAPolarizer, polarizeAxis );
			theObject->SetName(name);
			theObject->SetAcceptBallisticPhotons(acceptsBallisticPhotons);
		} else {
			// Object with optical properties
			MCRandomScatterer* scatterer;
			double opticalActivity = 0;
			
			CFDictionaryRef propertiesDict;
			if ( CFDictionaryGetValueIfPresent(objectProperties, CFSTR("properties"), (const void**)&propertiesDict) ) {
				scatterer = CreateRandomScattererFromCFDictionary(propertiesDict);
				if ( CFDictionaryGetValueIfPresent(propertiesDict, CFSTR("optical"), (const void**)&aDict) ) {
					GetNumberInDictWithDefault(aDict, CFSTR("opticalActivity"), kCFNumberDoubleType, &opticalActivity);
				}
			} else {
				scatterer = NULL;
			}
			
			CFDictionaryRef statisticsDict;
			if ( CFDictionaryGetValueIfPresent(objectProperties, CFSTR("statistics"), (const void**)&statisticsDict) ) {
				statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kStokesVector")) ? kStokesVector : 0;
				statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kIntensity")) ? kIntensityOnly : 0;
				statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kEntering")) ? kEntering : 0;
				statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kLeaving")) ? kLeaving : 0;
				statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kTimeResolved")) ? kTimeResolved : 0;
				acceptsBallisticPhotons = GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kAcceptsBallisticPhotons"));
				
				acceptanceCosineElements = 1;
				GetNumberInDictWithDefault(statisticsDict,CFSTR("acceptanceCosineElements"), kCFNumberLongType, &acceptanceCosineElements);
			}
			// End common
			
			if ( CFStringCompare (typeString, CFSTR("box"), kCFCompareCaseInsensitive) == kCFCompareEqualTo ) {
				double width, height, depth;
				long Nx, Ny, Nz;
				
				GetNumberInDict(objectProperties, CFSTR("width"), kCFNumberDoubleType, &width);
				GetNumberInDict(objectProperties, CFSTR("height"), kCFNumberDoubleType, &height);
				GetNumberInDict(objectProperties, CFSTR("depth"), kCFNumberDoubleType, &depth);
				
				GetNumberInDict(objectProperties, CFSTR("Nx"), kCFNumberLongType, &Nx);
				GetNumberInDict(objectProperties, CFSTR("Ny"), kCFNumberLongType, &Ny);
				GetNumberInDict(objectProperties, CFSTR("Nz"), kCFNumberLongType, &Nz);
				
				theObject = new MCBox(width, height, depth, Nx, Ny, Nz, 
									  1., 1., opticalActivity,
									  statisticsType, acceptanceCosineElements, false );
				theObject->SetRandomScatterer(scatterer);
				theObject->SetName(name);
				theObject->SetAcceptBallisticPhotons(acceptsBallisticPhotons);
				
			} else if ( CFStringCompare (typeString, CFSTR("layer"), kCFCompareCaseInsensitive) == kCFCompareEqualTo ) {
				double width, height, thickness;
				long Nx, Ny;
				
				GetNumberInDict(objectProperties, CFSTR("width"), kCFNumberDoubleType, &width);
				GetNumberInDict(objectProperties, CFSTR("height"), kCFNumberDoubleType, &height);
				GetNumberInDict(objectProperties, CFSTR("thickness"), kCFNumberDoubleType, &thickness);
				
				GetNumberInDict(objectProperties, CFSTR("Nx"), kCFNumberLongType, &Nx);
				GetNumberInDict(objectProperties, CFSTR("Ny"), kCFNumberLongType, &Ny);
				
				theObject = new MCInfiniteLayer(thickness, width, height, Nx, Ny, 
												1., 1., opticalActivity,
												statisticsType, acceptanceCosineElements, false );
				theObject->SetRandomScatterer(scatterer);
				theObject->SetName(name);
				theObject->SetAcceptBallisticPhotons(acceptsBallisticPhotons);
				
				
				
			} else if ( CFStringCompare (typeString, CFSTR("ellipsoid"), kCFCompareCaseInsensitive) == kCFCompareEqualTo ) {
				double width, height, depth;
				long Nx, Ny;
				
				GetNumberInDict(objectProperties, CFSTR("width"), kCFNumberDoubleType, &width);
				GetNumberInDict(objectProperties, CFSTR("height"), kCFNumberDoubleType, &height);
				GetNumberInDict(objectProperties, CFSTR("depth"), kCFNumberDoubleType, &depth);
				
				GetNumberInDict(objectProperties, CFSTR("Nx"), kCFNumberLongType, &Nx);
				GetNumberInDict(objectProperties, CFSTR("Ny"), kCFNumberLongType, &Ny);
				
				theObject = new MCEllipsoid(width, height, depth, Nx, Ny, 
											1., 1., opticalActivity,
											statisticsType, acceptanceCosineElements, false );
				theObject->SetRandomScatterer(scatterer);
				theObject->SetName(name);
				theObject->SetAcceptBallisticPhotons(acceptsBallisticPhotons);
			} else if ( CFStringCompare (typeString, CFSTR("cylinder"), kCFCompareCaseInsensitive) == kCFCompareEqualTo ) {
				double radius, height;
				long Nr, Nh;
				
				GetNumberInDict(objectProperties, CFSTR("radius"), kCFNumberDoubleType, &radius);
				GetNumberInDict(objectProperties, CFSTR("height"), kCFNumberDoubleType, &height);
				
				GetNumberInDict(objectProperties, CFSTR("Nr"), kCFNumberLongType, &Nr);
				GetNumberInDict(objectProperties, CFSTR("Nh"), kCFNumberLongType, &Nh);
				
				theObject = new MCCylinder(radius, height, Nr, Nh, 
										   1., 1., opticalActivity,
										   statisticsType, acceptanceCosineElements, false );
				theObject->SetRandomScatterer(scatterer);
				theObject->SetName(name);
				theObject->SetAcceptBallisticPhotons(acceptsBallisticPhotons);
			} else {
				clog << RightNow() << "Invalid object type\n" ;
				CFShow(objectProperties);
				theObject = NULL;
			}
			
			theObject->RotateObject(tiltX, tiltY, tiltZ);
		}
	} else {
		clog << RightNow() << "Undefined object type\n" ;
	}
	
	return theObject;
	
}

MCRandomScatterer* 
CreateRandomScattererFromCFDictionary(CFDictionaryRef materialProperties)
{
	double mu_s, mu_a, g, index, radius, indexScatterer, wavelength = 0;
	int ptsInTable = 10000;
	CFDictionaryRef aDict;
	CFStringRef aModel;
	
	if ( CFDictionaryGetValueIfPresent(materialProperties, CFSTR("scattering"), (const void**)&aDict) ) {
		
		if ( CFDictionaryContainsKey(aDict, CFSTR("model"))  ) {
			CFDictionaryGetValueIfPresent(aDict, CFSTR("model"), (const void**)&aModel);
/*			//by wendy
			clog << RightNow() << "Wendy debug: The scattering model is "<<CFStringToString(aModel)<<"\n";
*/		
		} else {
			aModel = CFSTR("HG");
		}
		
		mu_s = 0;
		GetNumberInDictWithDefault(aDict, CFSTR("mu_s"), kCFNumberDoubleType, &mu_s);
		mu_a = 0;
		GetNumberInDictWithDefault(aDict, CFSTR("mu_a"), kCFNumberDoubleType, &mu_a);
		
		if ( CFDictionaryContainsKey(aDict, CFSTR("g")) ) {
			GetNumberInDict(aDict, CFSTR("g"), kCFNumberDoubleType, &g);
		}
		
		if ( CFDictionaryContainsKey(aDict, CFSTR("radius") ) ) {
			GetNumberInDict(aDict, CFSTR("radius"), kCFNumberDoubleType, &radius);
		}
		
		if ( CFDictionaryContainsKey(aDict, CFSTR("indexScatterer") ) ) {
			GetNumberInDict(aDict, CFSTR("indexScatterer"), kCFNumberDoubleType, &indexScatterer);
		}
		if ( CFDictionaryContainsKey(aDict, CFSTR("wavelength") ) ) {
			GetNumberInDict(aDict, CFSTR("wavelength"), kCFNumberDoubleType, &wavelength);
		}
		
		if ( CFDictionaryContainsKey(aDict, CFSTR("ptsInTable") ) ) {
			GetNumberInDict(aDict, CFSTR("ptsInTable"), kCFNumberDoubleType, &ptsInTable);
		} else {
			clog << RightNow() << "Using default 1000 points for lookup table in Mie scatterer\n";
		}
		
	}
	
	
	if ( CFDictionaryGetValueIfPresent(materialProperties, CFSTR("optical"), (const void**)&aDict) ) {
		index = 1;
		GetNumberInDictWithDefault(aDict, CFSTR("index"), kCFNumberDoubleType, &index);
	}
	
	if ( CFStringCompare (aModel, CFSTR("HG"), kCFCompareCaseInsensitive) == kCFCompareEqualTo ) {
		return new MCRandomScattererHenyeyGreenstein(g,mu_s, mu_a, index, 0);
	} else if ( CFStringCompare (aModel, CFSTR("Kaplan"), kCFCompareCaseInsensitive) == kCFCompareEqualTo ) {
		return new MCRandomScattererKaplan(radius, indexScatterer, index, wavelength, ptsInTable, mu_s, mu_a, 0);
	} else if ( CFStringCompare (aModel, CFSTR("Jaillon"), kCFCompareCaseInsensitive) == kCFCompareEqualTo ) {
		return new MCRandomScattererJaillon(radius, indexScatterer, index, wavelength, ptsInTable, mu_s, mu_a, 0);
	} else {
		ThrowRuntimeError("No model for scattering provided in dictionary " << CFDescriptionToString(materialProperties));
	}
	
	return NULL;
}

void
SetStatisticsForCalculation(MCWorld* world, CFDictionaryRef calculationDict)
{
	
	CFDictionaryRef statisticsDict;
	unsigned long statistics = 0;
	
	if ( CFDictionaryGetValueIfPresent(calculationDict, CFSTR("statistics"), (const void**)&statisticsDict) ) {
		statistics = GetStatisticsBooleanSettingsFromCFDictionary(statisticsDict);
		MCWorld::SetOutputProperties(statistics);
	}
	
	CFStringRef aString;
	if ( CFDictionaryGetValueIfPresent(calculationDict, CFSTR("outputFilename"), (const void**)&aString) ) {
		world->SetOutputFilename( CFStringToString(aString) );
	}
	
	if ( CFDictionaryGetValueIfPresent(calculationDict, CFSTR("outputFilenameGeometry"), (const void**)&aString) ) {
		world->SetOutputFilenameGeometry( CFStringToString(aString) );

		if ( CFDictionaryGetValueIfPresent(calculationDict, CFSTR("outputFilenameGeometryFormat"), (const void**)&aString) ) {
			string theFormat = CFStringToString(aString);
			long format;
			if ( theFormat == "kBinary3DFormat" ) {
				format = kBinary3DFormat;
			} else if ( theFormat == "kMathematicaPolygons" ) {
				format = kMathematicaPolygons;
			} else if ( theFormat == "kGenericXML" ) {
				format = kGenericXML;
			} else if ( theFormat == "kMatlabPatch" ) {
				format = kMatlabPatch;
			} else if ( theFormat == "kWaveFrontObj" ) {
				format = kWaveFrontObj;
			} else {
				ThrowRuntimeError("Format not recognized for outputFilenameGeometryFormat" << theFormat << ".  Valid options are: kMathematicaPolygons, kGenericXML,kMatlabPatch,kWaveFrontObj,kBinary3DFormat\n");
			}
			world->SetOutputFilenameGeometryFormat(format);
		}		
	}
	
	if ( CFDictionaryGetValueIfPresent(calculationDict, CFSTR("outputFilenameGeometryWithIntensity"), (const void**)&aString) ) {
		world->SetOutputFilenameGeometryWithIntensity( CFStringToString(aString) );
	}
	
	
	if ( CFDictionaryGetValueIfPresent(calculationDict, CFSTR("outputFilenamePhotonPaths"), (const void**)&aString) ) {
		world->SetOutputFilenamePhotonPaths( CFStringToString(aString) );
	
		if ( CFDictionaryGetValueIfPresent(calculationDict, CFSTR("outputFilenamePhotonPathsFormat"), (const void**)&aString) ) {
			string theFormat = CFStringToString(aString);
			long format;
			if ( theFormat == "kRaw" ) {
				format = kRaw;
			} else if ( theFormat == "kMathematica" ) {
				format = kMathematica;
			} else if ( theFormat == "kBinaryFormat" ) {
				format = kBinaryFormat;
			} else {
				ThrowRuntimeError("Format not recognized for outputFilenamePhotonPathsFormat" << theFormat << ".  Valid options are: kMathematica, kRaw, kBinaryFormat\n");
			}
			
			
			
			world->SetOutputPhotonPathsFormat( format );
		}
		
	}

	// added Wendy Kan 7/8/10
	if ( CFDictionaryGetValueIfPresent(calculationDict, CFSTR("maximumDepthTargetObjectName"), (const void**)&aString) ) {
		world->SetmaximumDepthTargetObjectName( CFStringToString(aString) );
	}
	
	world->SetKeepPhotonStatistics( GetBooleanInDictWithDefaultFalse(calculationDict, CFSTR("keepAllPhotonStats")) );
	
	world->SetDebugOutputPhotonStatistics( GetBooleanInDictWithDefaultFalse(calculationDict, CFSTR("debugOutputPhotonStats")) );

	GetNumberInDictWithDefault(calculationDict,CFSTR("debugLevel"), kCFNumberLongType, &gDebugLevel);
	
}

unsigned long 
GetStatisticsBooleanSettingsFromCFDictionary(CFDictionaryRef statisticsDict)
{
	unsigned long statisticsType = 0;
	CFBooleanRef aBoolean;
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kEnergyDistributionBinUnits"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kEnergyDistributionBinUnits;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kEnergyDistributionRawBinary"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kEnergyDistributionRawBinary;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kEnergyDistributionXYZUnits"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kEnergyDistributionXYZUnits;
		}
	}
	
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kFluenceDistributionBinUnits"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kFluenceDistributionBinUnits;
		}
	}
	
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kFluenceDistributionRawBinary"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kFluenceDistributionRawBinary;
		}
	}
	
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kFluenceDistributionXYZUnits"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kFluenceDistributionXYZUnits;
		}
	}
	
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kGeometryMathematicaFormat"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kGeometryMathematicaFormat;
		}
	}
	
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kGeometrySeparateWavefrontFormat"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kGeometrySeparateWavefrontFormat;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kGeometryWorldWavefrontFormat"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kGeometryWorldWavefrontFormat;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kGeometryXMLFormat"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kGeometryXMLFormat;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kObjectInterfaceAverageNumberOfScatteringEvents"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kObjectInterfaceAverageNumberOfScatteringEvents;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kObjectInterfaceAveragePathLength"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kObjectInterfaceAveragePathLength;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kObjectInterfaceAveragePolarizedPathLength"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kObjectInterfaceAveragePolarizedPathLength;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kObjectInterfaceConstructionParameters"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kObjectInterfaceConstructionParameters;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kObjectInterfaceStokesVectorBetaCircular"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kObjectInterfaceStokesVectorBetaCircular;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kObjectInterfaceStokesVectorBetaLinear"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kObjectInterfaceStokesVectorBetaLinear;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kObjectInterfaceStokesVectorI"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kObjectInterfaceStokesVectorI;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kObjectInterfaceStokesVectorQ"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kObjectInterfaceStokesVectorQ;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kObjectInterfaceStokesVectorU"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kObjectInterfaceStokesVectorU;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kObjectInterfaceStokesVectorV"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kObjectInterfaceStokesVectorV;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kObjectInterfaceTimeResolvedIntensity"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kObjectInterfaceTimeResolvedIntensity;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kObjectInterfaceTotalTransmittance"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kObjectInterfaceTotalTransmittance;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kObjectPhotonScatteringStats"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kObjectPhotonScatteringStats;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kObjectTotalAbsorbance"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kObjectTotalAbsorbance;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kScattererStatistics1D"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kScattererStatistics1D;
		}
	}
	
	if ( CFDictionaryGetValueIfPresent(statisticsDict, CFSTR("kScattererStatistics2D"), (const void**)&aBoolean) ) {
		if ( CFBooleanGetValue(aBoolean) ) {
			statisticsType |= kScattererStatistics2D;
		}
	}
	
	
	return statisticsType;
}

MCSource* 
CreateLightSourceFromCFDictionary(CFDictionaryRef lightSourceProperties)
{
	CFStringRef aString;
	CFNumberRef theNumber;
	
	MCSource* theSource;
	
	if ( CFDictionaryGetValueIfPresent(lightSourceProperties, CFSTR("type"), (const void**)&aString) ) {
		/*
		 Common
		 */
		long nPhotons;
		if ( CFDictionaryGetValueIfPresent(lightSourceProperties, CFSTR("numberOfPhotons"), (const void**)&theNumber) ) {
			CFNumberGetValue(theNumber, kCFNumberLongType, &nPhotons);
		} else {
			nPhotons = 0;
			ThrowRuntimeError("Light source must have a non-zero numberOfPhotons");
		}
		
		if ( CFStringCompare (aString, CFSTR("isotropic"), kCFCompareCaseInsensitive) == kCFCompareEqualTo ) {
			RealV pos;
			if ( CFDictionaryGetValueIfPresent(lightSourceProperties, CFSTR("position"), (const void**)&aString) ) {
				pos = CFStringToRealV(aString);
			} else {
				clog << RightNow() << "No position specified for point source\n" ;
				return NULL;
			}
			double wavelength;
			if ( CFDictionaryGetValueIfPresent(lightSourceProperties, CFSTR("wavelength"), (const void**)&theNumber) ) {
				CFNumberGetValue(theNumber, kCFNumberDoubleType, &wavelength);
			} else {
				wavelength = 0;
			}
			theSource =  new MCIsotropicPointSource(pos, wavelength);
			PhotonIntensity aPhoton;
			theSource->SetTemplatePhoton(&aPhoton);
			theSource->SetNumberOfPhotonsToLaunch(nPhotons);
			
		} else if ( CFStringCompare (aString, CFSTR("laser"), kCFCompareCaseInsensitive) == kCFCompareEqualTo ) {
			RealV dir, pos;
			
			if ( CFDictionaryGetValueIfPresent(lightSourceProperties, CFSTR("direction"), (const void**)&aString) ) {
				dir = CFStringToRealV(aString);
				dir.normalize();
			} else {
				clog << RightNow() << "No direction specified for light source\n" ;
				return NULL;
			}		
			
			if ( CFDictionaryGetValueIfPresent(lightSourceProperties, CFSTR("position"), (const void**)&aString) ) {
				pos = CFStringToRealV(aString);
			} else {
				clog << RightNow() << "No position specified for light source\n" ;
				return NULL;
			}
			
			double Io, Qo, Uo, Vo;
			RealV eParallel(1,0,0);
			CFDictionaryRef aDict;
			if ( CFDictionaryGetValueIfPresent(lightSourceProperties, CFSTR("polarization"), (const void**)&aDict) ) {
				// Polarization is present
				if ( CFDictionaryGetValueIfPresent(aDict, CFSTR("I"), (const void**)&theNumber) ) {
					CFNumberGetValue(theNumber, kCFNumberDoubleType, &Io);
				} else {
					clog << RightNow() << "No Stokes intensity I specified for light source\n" ;
					return NULL;
				}
				if ( CFDictionaryGetValueIfPresent(aDict, CFSTR("Q"), (const void**)&theNumber) ) {
					CFNumberGetValue(theNumber, kCFNumberDoubleType, &Qo);
				} else {
					clog << RightNow() << "No Stokes intensity Q specified for light source\n" ;
					return NULL;
				}
				if ( CFDictionaryGetValueIfPresent(aDict, CFSTR("U"), (const void**)&theNumber) ) {
					CFNumberGetValue(theNumber, kCFNumberDoubleType, &Uo);
				} else {
					clog << RightNow() << "No Stokes intensity U specified for light source\n" ;
					return NULL;
				}
				if ( CFDictionaryGetValueIfPresent(aDict, CFSTR("V"), (const void**)&theNumber) ) {
					CFNumberGetValue(theNumber, kCFNumberDoubleType, &Vo);
				} else {
					clog << RightNow() << "No Stokes intensity V specified for light source\n" ;
					return NULL;
				}
				if ( CFDictionaryGetValueIfPresent(aDict, CFSTR("eParallel"), (const void**)&aString) ) {
					eParallel = CFStringToRealV(aString);
				} else {
					clog << RightNow() << "No Stokes eParallel specified for light source\n" ;
					return NULL;
				}
			} else {
				Io = 1;
				Qo = 0;
				Uo = 0;
				Vo = 0;
			}
			
			double beamDiameter;
			if ( CFDictionaryGetValueIfPresent(lightSourceProperties, CFSTR("beamDiameter"), (const void**)&theNumber) ) {
				CFNumberGetValue(theNumber, kCFNumberDoubleType, &beamDiameter);
			} else {
				beamDiameter = 0;
			}
			
			double wavefrontCurvature;
			if ( CFDictionaryGetValueIfPresent(lightSourceProperties, CFSTR("wavefrontCurvature"), (const void**)&theNumber) ) {
				CFNumberGetValue(theNumber, kCFNumberDoubleType, &wavefrontCurvature);
			} else {
				wavefrontCurvature = 0;
			}
			
			double wavelength;
			if ( CFDictionaryGetValueIfPresent(lightSourceProperties, CFSTR("wavelength"), (const void**)&theNumber) ) {
				CFNumberGetValue(theNumber, kCFNumberDoubleType, &wavelength);
			} else {
				wavelength = 0;
			}
			
			theSource = new MCLaserSource(StokesV(Io,Qo,Uo,Vo, eParallel, dir), beamDiameter/2., wavefrontCurvature, 0, pos, dir, wavelength);
			
			CFStringRef theFormalism;
			if ( CFDictionaryGetValueIfPresent(lightSourceProperties, CFSTR("formalism"), (const void**)&theFormalism) ) {
				
			} else {
				theFormalism = CFSTR("intensity");
				clog << RightNow() << "No formalism stated for light source, therefore using intensity only.\n" ;
			}
			
			if (CFStringCompare (theFormalism, CFSTR("intensity"), kCFCompareCaseInsensitive) == kCFCompareEqualTo) {
				PhotonIntensity aPhoton;
				theSource->SetTemplatePhoton(&aPhoton);
				theSource->SetNumberOfPhotonsToLaunch(nPhotons);
			} else if (CFStringCompare (theFormalism, CFSTR("Cote"), kCFCompareCaseInsensitive) == kCFCompareEqualTo) {
				PhotonCote aPhoton;
				theSource->SetTemplatePhoton(&aPhoton);
				theSource->SetNumberOfPhotonsToLaunch(nPhotons);
			} else if (CFStringCompare (theFormalism, CFSTR("Jaillon"), kCFCompareCaseInsensitive) == kCFCompareEqualTo) {
				PhotonJaillon aPhoton;
				theSource->SetTemplatePhoton(&aPhoton);
				theSource->SetNumberOfPhotonsToLaunch(nPhotons);
			}
			
		}
	} else {
		clog << RightNow() << "Undefined light source type\n" ;
	}
	
	return theSource;
	
}


void 
PopulateRootObjectWithObjectsFromCFDictionary(MCWorld* world, MCObject* rootObject, CFArrayRef objects)
{
	for ( long i = 0; i < CFArrayGetCount(objects) ; i++ ) {
		CFDictionaryRef anObject = (CFDictionaryRef)CFArrayGetValueAtIndex (objects, i) ;
		
		CFBooleanRef enabled;
		if (CFDictionaryGetValueIfPresent(anObject, CFSTR("enabled"), (const void**)&enabled)) {
			if ( ! CFBooleanGetValue(enabled) ) {
				continue;
			}
		}
		
		CFStringRef aString;
		
		// If object is generic file, we load from the world (may contain many objects)
		if ( CFDictionaryGetValueIfPresent(anObject, CFSTR("type"), (const void**)&aString) ) {
			if ( CFStringCompare (aString, CFSTR("generic"), kCFCompareCaseInsensitive) == kCFCompareEqualTo ) {
				string filename;
				if ( CFDictionaryGetValueIfPresent(anObject, CFSTR("filename"), (const void**)&aString) ) {
					filename = CFStringGetCStringPtr(aString, gCFStringEncoding);
				} else {
					ThrowRuntimeError("No filename specified for generic object") ;
				}
				//				throw runtime_error("Generic object only possible at top level, representing entire world");
				world->ReadWorldFromFile(filename, 0);
				
				// read and overwrite properties if present
				string objectName;
				if ( CFDictionaryGetValueIfPresent(anObject, CFSTR("name"), (const void**)&aString) ) {
					objectName = CFStringGetCStringPtr(aString, gCFStringEncoding);
					
					long statisticsType = kNone, acceptanceCosineElements = 1;
					bool acceptsBallisticPhotons = FALSE;
					MCRandomScatterer* scatterer;
					double opticalActivity = 0;
					
					CFDictionaryRef propertiesDict;
					if ( CFDictionaryGetValueIfPresent(anObject, CFSTR("properties"), (const void**)&propertiesDict) ) {
						scatterer = CreateRandomScattererFromCFDictionary(propertiesDict);
						CFDictionaryRef aDict;
						if ( CFDictionaryGetValueIfPresent(propertiesDict, CFSTR("optical"), (const void**)&aDict) ) {
							GetNumberInDictWithDefault(aDict, CFSTR("opticalActivity"), kCFNumberDoubleType, &opticalActivity);
						}
					} else {
						scatterer = NULL;
					}
					
					CFDictionaryRef statisticsDict;
					if ( CFDictionaryGetValueIfPresent(anObject, CFSTR("statistics"), (const void**)&statisticsDict) ) {
						statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kStokesVector")) ? kStokesVector : 0;
						statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kIntensity")) ? kIntensityOnly : 0;
						statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kEntering")) ? kEntering : 0;
						statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kLeaving")) ? kLeaving : 0;
						statisticsType |= GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kTimeResolved")) ? kTimeResolved : 0;
						acceptsBallisticPhotons = GetBooleanInDictWithDefaultFalse(statisticsDict, CFSTR("kAcceptsBallisticPhotons"));
						
						acceptanceCosineElements = 1;
						GetNumberInDictWithDefault(statisticsDict,CFSTR("acceptanceCosineElements"), kCFNumberLongType, &acceptanceCosineElements);
					}
					
					MCObject* theObject = world->FindObjectByName(objectName);
					if ( scatterer != NULL && theObject != NULL) {
						//	if ( theObject->GetRandomScatterer() == NULL ) {
						theObject->SetRandomScatterer(scatterer);
						//	} else {
						;//throw runtime_error("Object has properties both in .obj file and in property file .plist");
						//	}
						theObject->SetAcceptBallisticPhotons(acceptsBallisticPhotons);
						theObject->SetDetectionStatisticsOnSurfaceElements(statisticsType);
					}
				}
				
			} else if ( CFStringCompare (aString, CFSTR("stack"), kCFCompareCaseInsensitive) == kCFCompareEqualTo ) {
				// This is different from one layer: it is an array of stacked layers and we need to explicitly
				// state that the layers are touching
				CFArrayRef layers = (CFArrayRef)CFDictionaryGetValue(anObject, CFSTR("layers"));
				
				
				for ( i = 0; i < CFArrayGetCount(layers) ; i++ ) {
					CFDictionaryRef anObject = (CFDictionaryRef)CFArrayGetValueAtIndex (layers, i) ;
					
					if ( CFDictionaryGetValueIfPresent(anObject, CFSTR("type"), (const void**)&aString)  ) {
						if ( CFStringCompare (aString, CFSTR("layer"), kCFCompareCaseInsensitive) != kCFCompareEqualTo ) {
							throw runtime_error("Can only stack layers");
						}
					}
					
					
					MCInfiniteLayer* object = (MCInfiniteLayer*) CreateObjectFromCFDictionary( anObject );
										
					RealV pos;
					if ( CFDictionaryGetValueIfPresent(anObject, CFSTR("position"), (const void**)&aString) ) {
						pos = CFStringToRealV(aString);
					} else {
						clog << RightNow() << "No position specified for object: " << object->GetName();
						return;
					}
					object->SetOuterObjectTo(rootObject);
					
					clog << RightNow() << "Hi: we are placing the object: " << object->GetName() << endl;
					world->Place(object,pos);
					//throw runtime_error("Stack building incomplete");
					
					// added by Wendy Kan 7/13/10
					vector<MCObject*> listOfCreatedObjects(world->GetObjectList());
					
					for (int j=0; j<i; j++){
												
						MCObject* prevObject = listOfCreatedObjects[j];
						
						clog<< "going to touch interface, prevObject = " << prevObject->GetName() << " and thisObject:" << object->GetName() << endl;
						
						object->InterfaceTouchesOtherObject(prevObject, kForwardZPlane, kBackwardZPlane);
						object->InterfaceTouchesOtherObject(prevObject, kBackwardZPlane, kForwardZPlane);
						object->InterfaceTouchesOtherObject(prevObject, kForwardZPlane, kForwardZPlane);
						object->InterfaceTouchesOtherObject(prevObject, kBackwardZPlane, kBackwardZPlane);
						
					}
						
				}
				
/*				// this part added by wendy, modified 7/11/10
				for ( i = 0; i < CFArrayGetCount(layers) ; i++ ) {
					// for each layer, check all the other SE's to see if they're touching
						CFDictionaryRef thisObjectRef = (CFDictionaryRef)CFArrayGetValueAtIndex (layers, i) ;
						MCInfiniteLayer* thisObject = (MCInfiniteLayer*) CreateObjectFromCFDictionary( thisObjectRef );
					for (int j=0; j<i; j++){
						CFDictionaryRef otherObject = (CFDictionaryRef)CFArrayGetValueAtIndex (layers, j) ;
						MCObject* prevObject = CreateObjectFromCFDictionary( otherObject );
						
						thisObject->InterfaceTouchesOtherObject(prevObject, kForwardZPlane, kBackwardZPlane);
					}
				}
*/			
				
			} else {
				MCObject* object = CreateObjectFromCFDictionary( anObject );
				
				RealV pos;
				if ( CFDictionaryGetValueIfPresent(anObject, CFSTR("position"), (const void**)&aString) ) {
					pos = CFStringToRealV(aString);
				} else {
					clog << RightNow() << "No position specified for object: " << object->GetName();
					return;
				}
				object->SetOuterObjectTo(rootObject);
				world->Place(object,pos);
				
				// object->DumpGeometryToStream(clog, false, kWaveFrontObj);
				
				CFArrayRef anArrayOfObjects;
				if ( CFDictionaryGetValueIfPresent(anObject, CFSTR("objects"), (const void**)&anArrayOfObjects) ) {
					PopulateRootObjectWithObjectsFromCFDictionary(world, object, anArrayOfObjects);
				}
				
			}
		} else {
			throw runtime_error("Object has no type");
		}
	}
	
}

void 
PopulateWorldWithLightSourcesFromCFDictionary(MCWorld* world, CFArrayRef sources)
{
	for ( long i = 0; i < CFArrayGetCount(sources) ; i++ ) {
		CFDictionaryRef aLightSource = (CFDictionaryRef)CFArrayGetValueAtIndex (sources, i) ;
		MCSource* source = CreateLightSourceFromCFDictionary( aLightSource );
		
		world->PlaceLightSourceInWorld(source);
	}
}

string
CFStringToString(CFStringRef cf)
{ 
	char Str[255];
	CFStringGetCString(cf, Str, 256, gCFStringEncoding); 
	return string(Str); 
}		

RealV CFStringToRealV(CFStringRef theString)
{
	char Str[255];
	string theCPPString;
	if ( CFStringGetCString(theString,Str, 256,gCFStringEncoding) ) {
		theCPPString = string(Str);
	}
	
	istringstream s(theCPPString);
	RealV aVector;
	s >> aVector; 
	return aVector;
}
