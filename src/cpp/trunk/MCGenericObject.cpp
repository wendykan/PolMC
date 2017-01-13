#include "MCGenericObject.h"

#include "configfiles.h"
#include "MCObject.h"
#include "Photon.h"
#include "StokesV.h"
#include "MuellerM.h"
#include "MCRandomScatterer.h"
#include "MCUtils.h"
#include <sstream>
#include "XMLUtil.h"

void
MCGenericObject::init(map<string,string>& inDict)
{
    dict = inDict;

#if HAVE_LIBXML2
    cerr << "Nothing implemented yet." << endl;
#else
    throw logic_error("You must have libxml-2.0 installed to use generic objects (get it at http://www.xmlsoft.org/)");
#endif
    
}

bool
MCGenericObject::IsOutsideObject(RealV& inLocalPosition)
{
    return false;
}

#if ! HAVE_LIBXML2
MCObject*
MCGenericObject::BuildObjectFromFile(string inFilename)
{
    throw logic_error("You must have libxml-2.0 installed to use generic objects (get it at http://www.xmlsoft.org/)");
}
#else

MCObject*
MCGenericObject::BuildObjectFromFile(string inFilename)
{
    return BuildObjectFromXMLFile(inFilename);
}

MCObject*
MCGenericObject::BuildObjectFromXMLFile(string inFilename)
{
	
	xmlDocPtr doc;
	
	doc = (xmlDocPtr)xmlParseFile(inFilename.c_str());
	
	char outString[1000];
	istringstream input;
	double x,y,z;
	long n, l = 100, acceptanceCosElements;
	
    XMLUtil::extractSingleXPath(doc, BAD_CAST("/object/numInterfaces/text()"), outString);
	input.clear();
	input.str(outString);
	input >> n;
	if ( ! input ) {
		clog << RightNow() << "Invalid number of interfaces in XML file : \"" << outString << "\""<< endl;
		throw runtime_error("Invalid number of interfaces in XML file");
	} 
	
	XMLUtil::extractSingleXPath(doc, BAD_CAST("/object/acceptanceCosineElements/text()"), outString);
	input.clear();
	input.str(outString);
	input >> n;
	if ( ! input ) {
		clog << RightNow() << "Invalid acceptanceCosineElements in XML file : \"" << outString << "\""<< endl;
		throw runtime_error("Invalid acceptanceCosineElements in XML file");
	} else {
		n = 1;
	}
	acceptanceCosElements = n;

	long mNumOfInterfaces = n;

    SurfaceElement* surfaceElements = new SurfaceElement[mNumOfInterfaces];

	clog << RightNow() << "Reading " <<  mNumOfInterfaces <<  " surfaceElements " ;
	for (long i = 0; i < mNumOfInterfaces; i++) {
		string theXPath;
		ostringstream which;
  		which << (i+1);;
  
		string baseXPath = "/object/interface[position()=" + which.str() + "]";
		
				
		theXPath = baseXPath + "/name/text()";
		XMLUtil::extractSingleXPath(doc, BAD_CAST(theXPath.c_str()), outString);
		surfaceElements[i].name = string(outString);

		theXPath = baseXPath + "/origin/text()";
		XMLUtil::extractSingleXPath(doc, BAD_CAST(theXPath.c_str()), outString);
		if (ReadRealVFromString(outString, surfaceElements[i].origin)) {
			clog << RightNow() << "Invalid origin in XML file : \"" << outString << "\", at interface element " << i+1 << " of " << mNumOfInterfaces << endl;
			throw runtime_error("Invalid origin in XML file");
		}
		
		theXPath  = baseXPath + "/a/text()";
		XMLUtil::extractSingleXPath(doc, BAD_CAST(theXPath.c_str()), outString);
		if (ReadRealVFromString(outString, surfaceElements[i].a)) {
			clog << RightNow() << "Invalid a-vector in XML file : \"" << outString << "\", at interface element " << i+1 << " of " << mNumOfInterfaces << endl;
			throw runtime_error("Invalid a-vector in XML file");
		}

		theXPath  = baseXPath + "/b/text()";
		XMLUtil::extractSingleXPath(doc, BAD_CAST(theXPath.c_str()), outString);
		if (ReadRealVFromString(outString, surfaceElements[i].b)) {
			clog << RightNow() << "Invalid b-vector in XML file : \"" << outString << "\", at interface element " << i+1 << " of " << mNumOfInterfaces << endl;
			throw runtime_error("Invalid b-vector in XML file");
		}

		theXPath  = baseXPath + "/Na/text()";
		XMLUtil::extractSingleXPath(doc, BAD_CAST(theXPath.c_str()), outString);
		input.clear();
		input.str(outString);
		input >> n;
/*		if ( ! input ) {
			clog << RightNow() << "Invalid Na in XML file : \"" << outString << "\", at interface element " << i+1 << " of " << NumOfInterface << endl;
			throw runtime_error("Invalid Na in XML file");
		} */
		surfaceElements[i].Na = n;

		theXPath  = baseXPath + "/Nb/text()";
		XMLUtil::extractSingleXPath(doc, BAD_CAST(theXPath.c_str()), outString);
		input.clear();
		input.str(outString);
		input >> n;
/*		if ( ! input ) {
			clog << RightNow() << "Invalid Nb in XML file : \"" << outString << "\", at interface element " << i+1 << " of " << NumOfInterface << endl;
			throw runtime_error("Invalid Nb in XML file");
		} */
		surfaceElements[i].Nb = n;

#if WITH_POLARIZATION
		theXPath  = baseXPath + "/er/text()";
		XMLUtil::extractSingleXPath(doc, BAD_CAST(theXPath.c_str()), outString);
		if (ReadRealVFromString(outString, surfaceElements[i].er)) {
			clog << RightNow() << "Invalid er-vector in XML file : \"" << outString << "\", at interface element " << i+1 << " of " << mNumOfInterfaces << endl;
			throw runtime_error("Invalid er-vector in XML file");
		}

		theXPath  = baseXPath + "/el/text()";
		XMLUtil::extractSingleXPath(doc, BAD_CAST(theXPath.c_str()), outString);
		if (ReadRealVFromString(outString, surfaceElements[i].el)) {
			clog << RightNow() << "Invalid el-vector in XML file : \"" << outString << "\", at interface element " << i+1 << " of " << mNumOfInterfaces << endl;
			throw runtime_error("Invalid el-vector in XML file");
		}
#endif

        /*
		theXPath  = baseXPath + "/normal/text()";
		extractSingleXPath(doc, BAD_CAST(theXPath.c_str()), outString);
		if (ReadRealVFromString(outString, surfaceElements[i].normal)) {
			clog << RightNow() << "Invalid normal-vector in XML file : \"" << outString << "\", at interface element " << i+1 << " of " << mNumOfInterfaces << endl;
			throw runtime_error("Invalid normal-vector in XML file");
		}
		*/
        
		surfaceElements[i].normal = RealV::CrossProduct(surfaceElements[i].a, surfaceElements[i].b);
		surfaceElements[i].normal.normalize();

/*	
		if (RealV::NormalizedCrossProduct(surfaceElements[i].a, surfaceElements[i].b) != surfaceElements[i].normal) {
			clog << RightNow() << "Invalid normal-vector in XML file (not equals to a x b) : \"" << outString << "\", at interface element " << i+1 << " of " << NumOfInterface << endl;
			throw runtime_error("Invalid normal-vector in XML file");
 		}
 */
		theXPath  = baseXPath + "/indexIn/text()";
		XMLUtil::extractSingleXPath(doc, BAD_CAST(theXPath.c_str()), outString);
		input.clear();
		input.str(outString);
		input >> x;
/*		if ( input.fail()) {
			clog << RightNow() << "Invalid indexIn in XML file : \"" << outString << "\", at interface element " << i+1 << " of " << NumOfInterface << endl;
			throw runtime_error("Invalid indexIn in XML file");
		} */
//		surfaceElements[i].indexIn = x;

        /* Index is set when placing this object into another.  Will get index from parent. */
//		surfaceElements[i].indexOut = 1;

		theXPath  = baseXPath + "/surfaceShape/text()";
		XMLUtil::extractSingleXPath(doc, BAD_CAST(theXPath.c_str()), outString);
		
		if (outString == string("kTriangle")) 
			surfaceElements[i].surfaceShape = kTriangle;
		else if (outString == string("kParallelogram")) 
			surfaceElements[i].surfaceShape = kParallelogram;
		else if (outString == string("kInfinitePlane")) 
			surfaceElements[i].surfaceShape = kInfinitePlane;
		else {
			clog << RightNow() << "Invalid surfaceShape in XML file : \"" << outString << "\", at interface element " << i+1 << " of " << mNumOfInterfaces << endl;
			throw runtime_error("Invalid surfaceShape in XML file");
		}
		
		surfaceElements[i].mDetection = kIntensity;
		surfaceElements[i].mCosineElements = 1;
		surfaceElements[i].init();
		
	    if (i % 10 == 0)    {
	        clog << "." ;
	        clog.flush();
	    }

        mSurfaceElements.push_back(surfaceElements + i);
	}
	clog << endl;

    FinishCreate();
    
	return NULL;
}

long
MCGenericObject::BuildObjectFromOBJStream(istream& input)
{
//	ifstream input(inFilename.c_str());
    string s;
    const long max_size=1000;
    s.reserve(max_size);

    vector<string> facesStrings;

    mVertexList.clear();
    mNormalList.clear();
    mFaceList.clear();
    
    RealV v;
    
	double x=NAN,y=NAN,z=NAN;
	char c;

    bool hasSeenGroup = false;
    
	while (! input.eof() ) {
        input >> s;
        
        if (s[0] == '#') {
            // Comment
            getline(input, s);
            if (! input.good()) {
                return -1;
            } else {
                cout << "Read comment: " << s << endl;
            }
        } else if (s[0] == 'g') {
            if (! hasSeenGroup) {
                // Group name
                getline(input, s);
                if (! input.good()) {
                    return -1;
                } else {
                    cout << "Read group name: " << s << endl;
                    hasSeenGroup = true;
                }
                mObjectName=s;
            } else {
                // We should put the character back into the stream:
                // this will be a new object
                input.putback(c);
                return 1;
            }
            
        } else if (s == string("vn")) {
            input >> v;
        } else if (s == string("v")) {
            input >> v;
            mVertexList.push_back(v);
        } else if (s == string("f")) {
            getline(input, s);

            vector<string> vertices;
            split(s, vertices," ");

            Face theFace;
            for (long i = 0; i < vertices.size(); i++) {
                vector<string> indices;
                split(vertices[i], indices,"/");
                MyAssert_(indices.size() == 3);
                Vertex vertex;
                vertex.vertex = atoi(indices[0].c_str());
                vertex.texture = atoi(indices[1].c_str());
                vertex.normal = atoi(indices[2].c_str());
                
                theFace.push_back(vertex);
            }
            
            mFaceList.push_back(theFace);
            
        } else {
            getline(input, s);
        }
    }

    FinishCreate();
 
    return 0;
}

bool
MCGenericObject::ReadRealVFromString(string inString,  RealV& outVector)
{
	// Take a look at C++ Third Edition section 21.3.5 for example
	istringstream input;
	input.str(inString);

	double x=NAN,y=NAN,z=NAN;
	char c;
	
	input >> c;
	
	if (c == '(') {
		// We expect the form (x,y,z)
		input >> x >> c;
		if (c == ',') {
			input >> y >> c;
			if (c == ',') {
				input >> z >> c;
				if (c != ')')
					return true;
			} else {
				return true;
			}
		} else {
			return true;
		}
	} else {
		// We expect the form x y z
		input.putback(c);
		input >> x ;
		if (! input.good()) {
			return true;
		}

		input >> y;
		if (! input.good()) {
			return true;
		}

		input >> z;
		if ( ! input.eof()) {
			return true;
		}
	}

    outVector = RealV(x,y,z);

	return false;

}

void 
MCGenericObject::split(const string& str,
              vector<string>& tokens,
              const string& delimiters)
{
    // Skip delimiters at beginning.
    string::size_type start = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type current     = str.find_first_of(delimiters, start);
    
    if (string::npos == start) 
        return;
    
    while (string::npos != current || string::npos != start)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(start, current - start));
        start = current != string::npos ? current+1 : string::npos;
        // Find next "non-delimiter"
        current = str.find_first_of(delimiters, start);
    }
}

#endif

