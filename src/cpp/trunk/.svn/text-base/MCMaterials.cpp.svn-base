#if HAVE_LIBXML2
#include <libxml2/parser.h>
#include <libxml2/xpath.h>
#endif


#include "XMLUtil.h"
#include "MCUtils.h"
#include "MCMaterials.h"
#include "mydebug.h"

#include <sstream>


MCMaterialsList::MCMaterialsList(string inFilename)
{
    mMaterials.reserve(100);
    
    initFromXMLDocument(inFilename);
}

MCMaterialsList::MCMaterialsList(string inFilename, string inDefaultsFilename)
{
    mMaterials.reserve(100);
    
    initFromXMLDocument(inFilename);
    initDefaultsFromXMLDocument(inDefaultsFilename);
}

MCMaterialsList::~MCMaterialsList()
{
    for (unsigned long i = 0; i < mMaterials.size(); i++)
        delete mMaterials[i];
}

bool
MCMaterialsList::initFromXMLDocument(string inFilename)
{
    xmlDocPtr doc = NULL;
    xmlXPathContextPtr ctx = NULL, ctx2 = NULL;
    xmlXPathObjectPtr pobj = NULL;
    xmlNodeSetPtr nset = NULL;
    int i;
    
    if ((doc = xmlParseFile(inFilename.c_str())) == NULL ) {
        printf("parsing of properties.xml failed.  Check file: %s", inFilename.c_str());
        throw ;
    }
    
    xmlXPathInit();
    ctx = xmlXPathNewContext(doc);
    
    pobj = xmlXPathEvalExpression(BAD_CAST("//material"), ctx);
    nset = pobj->nodesetval;
    for (i = 0; i < nset->nodeNr; i++) {
        ctx2 = xmlXPathNewContext(doc);
        ctx2->node = xmlXPathNodeSetItem(nset, i);
        
        MCMaterial* m = new MCMaterial();
        m->initFromXMLDocument(ctx2);

        mMaterials.push_back(m);
        
        xmlXPathFreeContext(ctx2);
    }
    xmlXPathFreeObject(pobj);
    
    xmlFreeDoc(doc);

    return false;
}

bool
MCMaterialsList::initDefaultsFromXMLDocument(string inFilename)
{
    xmlDocPtr doc = NULL;
    xmlXPathContextPtr ctx = NULL, ctx2 = NULL;
    xmlXPathObjectPtr pobj = NULL;
    xmlNodeSetPtr nset = NULL;
    int i;
    
    if ((doc = xmlParseFile(inFilename.c_str())) == NULL ) {
        printf("parsing of defaults failed.  Check file.");
        throw ;
    }
    
    xmlXPathInit();
    ctx = xmlXPathNewContext(doc);
    
    pobj = xmlXPathEvalExpression(BAD_CAST("//material"), ctx);
    nset = pobj->nodesetval;
    for (i = 0; i < nset->nodeNr; i++) {
        ctx2 = xmlXPathNewContext(doc);
        ctx2->node = xmlXPathNodeSetItem(nset, i);
        
		string key;
		bool err = XMLUtil::InitVariable(ctx2, (char*) "./key/text()", key);
		trim(key);
		MCMaterial* m = GetMaterialByKey(key);

		if (m != NULL)
			m->initDefaultsFromXMLDocument(ctx2);

        xmlXPathFreeContext(ctx2);
    }
    xmlXPathFreeObject(pobj);
    
    xmlFreeDoc(doc);

    return false;
}

MCMaterial* 
MCMaterialsList::GetMaterialByName(string inName)
{
    for (unsigned long i = 0; i < mMaterials.size(); i++) {
        if (mMaterials[i]->GetName() == inName) {
            return mMaterials[i];
        }
    }
    
    return NULL;
}

MCMaterial* 
MCMaterialsList::GetMaterialByKey(string inKey)
{
    for (unsigned long i = 0; i < mMaterials.size(); i++) {
        if (mMaterials[i]->GetKey() == inKey) {
            return mMaterials[i];
        }
    }
    
    return NULL;
}

vector<string>
MCMaterialsList::GetMaterialsListByName()
{
    vector<string> list;
    
    for (unsigned long i = 0; i < mMaterials.size(); i++) {
        list.push_back(mMaterials[i]->GetName());
    }
    
    return list;
}

vector<string>
MCMaterialsList::GetMaterialsListByKey()
{
    vector<string> list;
    
    for (unsigned long i = 0; i < mMaterials.size(); i++) {
        list.push_back(mMaterials[i]->GetKey());
    }
    
    return list;
}

MCMaterial::MCMaterial()
{
    mConcentration = 1;
    mDefaultBackgroundIndex = 0;
    mDefaultScatteringCoefficient = 0;
    mDefaultAbsorptionCoefficient = 0;
    mDefaultAnisotropyCoefficient = 0;
}

MCMaterial::~MCMaterial()
{
    
}

MCMaterial& 
MCMaterial::operator=(const MCMaterial& inMaterial)
{
    mName = inMaterial.mName;
    mKey = inMaterial.mKey;
    mConcentration = inMaterial.mConcentration;
    mDescription = inMaterial.mDescription;
    
    for (long i = 0; i < inMaterial.mBackgroundIndexTables.size(); ++i) {
        DataSource* dataSource = new DataSource;
        *dataSource = *(inMaterial.mBackgroundIndexTables[i]);
        mBackgroundIndexTables.push_back(dataSource);
        
        if (inMaterial.mDefaultBackgroundIndex == inMaterial.mBackgroundIndexTables[i])
            mDefaultBackgroundIndex = dataSource;
    }

    for (long i = 0; i < inMaterial.mScatteringCoefficientTables.size(); ++i) {
        DataSource* dataSource = new DataSource;
        *dataSource = *(inMaterial.mScatteringCoefficientTables[i]);
        mScatteringCoefficientTables.push_back(dataSource);

        if (inMaterial.mDefaultScatteringCoefficient == inMaterial.mScatteringCoefficientTables[i])
            mDefaultScatteringCoefficient = dataSource;
    }

    for (long i = 0; i < inMaterial.mAbsorptionCoefficientTables.size(); ++i) {
        DataSource* dataSource = new DataSource;
        *dataSource = *(inMaterial.mAbsorptionCoefficientTables[i]);
        mAbsorptionCoefficientTables.push_back(dataSource);

        if (inMaterial.mDefaultAbsorptionCoefficient == inMaterial.mAbsorptionCoefficientTables[i])
            mDefaultAbsorptionCoefficient = dataSource;
    }

    for (long i = 0; i < inMaterial.mAnisotropyCoefficientTables.size(); ++i) {
        DataSource* dataSource = new DataSource;
        *dataSource = *(inMaterial.mAnisotropyCoefficientTables[i]);
        mAnisotropyCoefficientTables.push_back(dataSource);

        if (inMaterial.mDefaultAnisotropyCoefficient == inMaterial.mAnisotropyCoefficientTables[i])
            mDefaultAnisotropyCoefficient = dataSource;
    }
    
    
}

bool
MCMaterial::initFromXMLDocument(xmlXPathContextPtr inContext)
{
    bool err, err2;
    
    err = XMLUtil::InitVariable(inContext, (char*) "./name/text()", mName);
	trim(mName);
    err2 = XMLUtil::InitVariable(inContext, (char*) "./key/text()", mKey);
	trim(mKey);

    if (err && err2) {
        clog << RightNow() << " Warning : no key or name";
    }
    
    // Those can be empty, so we don't check errors.
//    XMLUtil::InitVariable(inContext, "./concentration/text()", mConcentration);

    initTablesFromNodeSet(inContext, (char*) "./properties/optical/indexBackground", mBackgroundIndexTables);
    initTablesFromNodeSet(inContext, (char*) "./properties/scattering/mu_s", mScatteringCoefficientTables);
    initTablesFromNodeSet(inContext, (char*) "./properties/scattering/mu_a", mAbsorptionCoefficientTables);
    initTablesFromNodeSet(inContext, (char*) "./properties/scattering/anisotropy", mAnisotropyCoefficientTables);

    return false;
}

bool
MCMaterial::initDefaultsFromXMLDocument(xmlXPathContextPtr inContext)
{
    string defaultKey;
	
    // Those can be empty, so we don't check errors.
    if ( ! XMLUtil::InitVariable(inContext, (char*) "./properties/optical/indexBackground/defaultKey/text()", defaultKey) )
		SetDefaultKeyForBackgroundIndex(defaultKey);
    if ( ! XMLUtil::InitVariable(inContext, (char*) "./properties/scattering/mu_s/defaultKey/text()", defaultKey) )
		SetDefaultKeyForScatteringCoefficient(defaultKey);
    if ( ! XMLUtil::InitVariable(inContext, (char*) "./properties/scattering/mu_a/defaultKey/text()", defaultKey) )
		SetDefaultKeyForAbsorptionCoefficient(defaultKey);
    if ( ! XMLUtil::InitVariable(inContext, (char*) "./properties/scattering/anisotropy/defaultKey/text()", defaultKey) )
		SetDefaultKeyForAnisotropyCoefficient(defaultKey);
	
    return false;
}

bool
MCMaterial::initTablesFromNodeSet(xmlXPathContextPtr inContext, char* inXPath, vector<DataSource*>& inDataSource)
{
    
    xmlXPathObjectPtr pobj = xmlXPathEvalExpression(BAD_CAST(inXPath), inContext);
    
    if (pobj == NULL) {
        return true;
    }
    
    xmlNodeSetPtr nodeSet = pobj->nodesetval;

    if (nodeSet == NULL) {
        return true;
    }
    xmlNodePtr saveRoot = inContext->node;

    for (long i = 0; i < nodeSet->nodeNr; i++) {
        inContext->node = xmlXPathNodeSetItem(nodeSet, i);
        DataSource* data = new DataSource;
        
        if (XMLUtil::InitVariable(inContext, (char*) "./key/text()", data->key)) {
            clog << RightNow() << " Warning: no key for table.\n";
        } else {
			trim(data->key);
		}
		
		XMLUtil::InitVariable(inContext, (char*) "./table/text()", data->table);
        XMLUtil::InitVariable(inContext, (char*) "./description/text()", data->description);
        inDataSource.push_back(data);
    }

    inContext->node = saveRoot;

    xmlXPathFreeObject(pobj);

    return false;
}

string
MCMaterial::GetName()
{
    return mName;
}

string
MCMaterial::GetKey()
{
    return mKey;
}

void
MCMaterial::SetDefaultKey(string inKey)
{
    SetDefaultKeyForBackgroundIndex(inKey);
    SetDefaultKeyForScatteringCoefficient(inKey);
    SetDefaultKeyForAbsorptionCoefficient(inKey);
    SetDefaultKeyForAnisotropyCoefficient(inKey);
}

void 
MCMaterial::SetDefaultKeyForBackgroundIndex(string inKey)
{   
	trim(inKey);
    for (unsigned long i = 0; i < mBackgroundIndexTables.size(); i++) {
        if (mBackgroundIndexTables[i]->key == inKey) {
            mDefaultBackgroundIndex = mBackgroundIndexTables[i];
        }
    }
}

void
MCMaterial::SetDefaultKeyForScatteringCoefficient(string inKey)
{
	trim(inKey);
    for (unsigned long i = 0; i < mScatteringCoefficientTables.size(); i++) {
        if (mScatteringCoefficientTables[i]->key == inKey) {
            mDefaultScatteringCoefficient = mScatteringCoefficientTables[i];
        }
    }
}

void
MCMaterial::SetDefaultKeyForAbsorptionCoefficient(string inKey)
{
	trim(inKey);
    for (unsigned long i = 0; i < mAbsorptionCoefficientTables.size(); i++) {
        if (mAbsorptionCoefficientTables[i]->key == inKey) {
            mDefaultAbsorptionCoefficient = mAbsorptionCoefficientTables[i];
        }
    }
}

void
MCMaterial::SetDefaultKeyForAnisotropyCoefficient(string inKey)
{
	trim(inKey);
    for (unsigned long i = 0; i < mAnisotropyCoefficientTables.size(); i++) {
        if (mAnisotropyCoefficientTables[i]->key == inKey) {
            mDefaultAnisotropyCoefficient = mAnisotropyCoefficientTables[i];
        }
    }
}
/*
void 
MCMaterial::SetConstantPropertyValueOverride(long inProperty, double inWavelength)
{
    return NAN;
}
*/

void 
MCMaterial::SetRelativeConcentration(double inConcentration)
{
    mConcentration = inConcentration;
}

double 
MCMaterial::GetPropertyValue(long inProperty, double inWavelength)
{
    return NAN;
}

fastinterpolate 
MCMaterial::GetPropertyDefaultTable(long inProperty)
{
    if (inProperty == kBackgroundIndex) {
        if (mDefaultBackgroundIndex) {
            return mDefaultBackgroundIndex->table;
        } else if (mBackgroundIndexTables.size() == 1) {
            mDefaultBackgroundIndex = mBackgroundIndexTables[0];
            return mBackgroundIndexTables[0]->table;
        } else {
            throw runtime_error("No default table for background index");
        }
    } else if (inProperty == kScatteringCoefficient) {
        if (mDefaultScatteringCoefficient)
            return mDefaultScatteringCoefficient->table;
        else if (mScatteringCoefficientTables.size() == 1) {
            mDefaultScatteringCoefficient = mScatteringCoefficientTables[0];
            return mScatteringCoefficientTables[0]->table;
        } else {
            throw runtime_error("No default table for scattering coefficient");
        }
    } else if (inProperty == kAbsorptionCoefficient) {
        if (mDefaultAbsorptionCoefficient)
            return mDefaultAbsorptionCoefficient->table;
        else if (mAbsorptionCoefficientTables.size() == 1) {
            mDefaultAbsorptionCoefficient = mAbsorptionCoefficientTables[0];
            return mAbsorptionCoefficientTables[0]->table;
        } else {
            throw runtime_error("No default table for absorption coefficient");
        }
    }
    else if (inProperty == kAnisotropyCoefficient) {
        if (mDefaultAnisotropyCoefficient)
            return mDefaultAnisotropyCoefficient->table;
        else if (mAnisotropyCoefficientTables.size() == 1) {
            mDefaultAnisotropyCoefficient = mAnisotropyCoefficientTables[0];
            return mAnisotropyCoefficientTables[0]->table;
        } else {
            throw runtime_error("No default table for anisotropy coefficient");
        }
    }

    throw runtime_error("Bad property or no default table");
}

fastinterpolate 
MCMaterial::GetPropertyTable(long inProperty, string inKey)
{
    if (inProperty == kBackgroundIndex) {
        for (unsigned long i = 0; i < mBackgroundIndexTables.size(); i++) {
            if (mBackgroundIndexTables[i]->key == inKey)
                return mBackgroundIndexTables[0]->table;
        }
        throw runtime_error("No table for background index with that key");
    } else if (inProperty == kScatteringCoefficient) {
        for (unsigned long i = 0; i < mScatteringCoefficientTables.size(); i++) {
            if (mScatteringCoefficientTables[i]->key == inKey)
                return mScatteringCoefficientTables[0]->table;
        }
        throw runtime_error("No table for scattering coefficient with that key");
    } else if (inProperty == kAbsorptionCoefficient) {
        for (unsigned long i = 0; i < mAbsorptionCoefficientTables.size(); i++) {
            if (mAbsorptionCoefficientTables[i]->key == inKey)
                return mAbsorptionCoefficientTables[0]->table;
        }
        throw runtime_error("No table for absorption coefficient with that key");
    }
    else if (inProperty == kAnisotropyCoefficient) {
        for (unsigned long i = 0; i < mAnisotropyCoefficientTables.size(); i++) {
            if (mAnisotropyCoefficientTables[i]->key == inKey)
                return mAnisotropyCoefficientTables[0]->table;
        }
        throw runtime_error("No table for anisoptropy coefficient with that key");
    }

}

vector<DataSource*>
MCMaterial::GetPropertyDataSource(long inProperty)
{
    vector<DataSource*> source;
    switch (inProperty) {
        case kBackgroundIndex:
            source = mBackgroundIndexTables;
            break;
        case kScatteringCoefficient:
            source = mScatteringCoefficientTables;
            break;
        case kAbsorptionCoefficient:
            source = mAbsorptionCoefficientTables;
            break;
        case kAnisotropyCoefficient:
            source = mAnisotropyCoefficientTables;
            break;
    }
    
    return source;
}

double 
MCMaterial::GetBackgroundIndex(double inWavelength)
{
    if (mDefaultBackgroundIndex) {
        return mDefaultBackgroundIndex->table.y(inWavelength);
    }
    
    if (mBackgroundIndexTables.size() == 1) {
        mDefaultBackgroundIndex = mBackgroundIndexTables[0];
        return mBackgroundIndexTables[0]->table.y(inWavelength);
    } else if (mBackgroundIndexTables.size() == 0) {
        return 1;
	} else if (mBackgroundIndexTables.size() > 1) {
        clog << RightNow() << " There are more than one table and default key for mBackgroundIndexTables not set. Possibilities are: ";
        for (unsigned long i = 0; i < mBackgroundIndexTables.size(); i++) {
            clog << mBackgroundIndexTables[i]->key << " ";
        }
        clog << "\n";
    } 
    
    throw runtime_error("No tables for background index of "+mName);
    
    return NAN;
}

double 
MCMaterial::GetScatteringCoefficient(double inWavelength)
{
    if (mDefaultScatteringCoefficient) {
        return mConcentration * mDefaultScatteringCoefficient->table.y(inWavelength);
    }
    
    if (mScatteringCoefficientTables.size() == 1) {
        mDefaultScatteringCoefficient = mScatteringCoefficientTables[0];
        return mConcentration * mScatteringCoefficientTables[0]->table.y(inWavelength);
    } else if (mScatteringCoefficientTables.size() == 0) {
        return 0;
    } else if (mScatteringCoefficientTables.size() > 1) {
        clog << RightNow() << " There are more than one table and default key for mScatteringCoefficientTables not set. Possibilities are: ";
        for (unsigned long i = 0; i < mScatteringCoefficientTables.size(); i++) {
            clog << mScatteringCoefficientTables[i]->key << " ";
        }
        clog << "\n";
    }

    throw runtime_error("No tables for scattering coefficient of "+mName);

    return NAN;
}

double 
MCMaterial::GetAbsorptionCoefficient(double inWavelength)
{
    if (mDefaultAbsorptionCoefficient) {
        return mConcentration * mDefaultAbsorptionCoefficient->table.y(inWavelength);
    }

    if (mAbsorptionCoefficientTables.size() == 1) {
        mDefaultAbsorptionCoefficient = mAbsorptionCoefficientTables[0];
        return mConcentration * mAbsorptionCoefficientTables[0]->table.y(inWavelength);
    } else if (mAbsorptionCoefficientTables.size() == 0) {
        return 0;
    } else if (mAbsorptionCoefficientTables.size() > 1) {
        clog << RightNow() << " There are more than one table and default key for mAbsorptionCoefficientTables not set. Possibilities are: ";

        for (unsigned long i = 0; i < mAbsorptionCoefficientTables.size(); i++) {
            clog << mAbsorptionCoefficientTables[i]->key << " ";
        }
        clog << "\n";
    }
    
    throw runtime_error("No tables for absorption coefficient of "+mName);

    return NAN;
}

double 
MCMaterial::GetAnisotropyCoefficient(double inWavelength)
{
    if (mDefaultAnisotropyCoefficient) {
        return mDefaultAnisotropyCoefficient->table.y(inWavelength);
    }
    
    if (mAnisotropyCoefficientTables.size() == 1) {
        mDefaultAnisotropyCoefficient = mAnisotropyCoefficientTables[0];
        return mAnisotropyCoefficientTables[0]->table.y(inWavelength);
    } else if (mAnisotropyCoefficientTables.size() == 0) {
        return 0;
    } else if (mAnisotropyCoefficientTables.size() > 1) {
        clog << RightNow() << " There are more than one table and default key for mAnisotropyCoefficientTables not set. Possibilities are: ";
        for (unsigned long i = 0; i < mAnisotropyCoefficientTables.size(); i++) {
            clog << mAnisotropyCoefficientTables[i]->key << " ";
        }
        clog << "\n";
    }
    
    throw runtime_error("No tables for anisotropy coefficient of "+mName);

    return NAN;
}

#ifdef STANDALONE_MATERIALS

#include <fstream>
#include <readline/readline.h>
#include <readline/history.h>
#include "MCUtils.h"

int
main()
{
    try {
        using_history();

        MCMaterialsList materialList("properties.xml","properties.defaults.xml");

        ostream* outStream = &cout;

        while ( 1 ) {
            
            try {
                auto_ptr<char> buf;
                buf.reset(readline ("Material [? for help]: "));

                string text(buf.get());
                trim(text);
                
                vector<string> tokens;

                split(text,tokens, " ");

                if (tokens.size() > 0) {
                    add_history(text.c_str());
                    if (tokens[0] == "?" || tokens[0] == "help") {
                        cout << " Commands are:\n";
                        cout << "   info  [mat key]       : by itself, lists all material keys available.\n"; 
						cout << "                           With a material key, gives info about that material (tables and their keys).\n";
                        cout << "   mu_s  [mat key] [key] : return scattering coefficient for material.\n";
                        cout << "                           If more than one table is available, you must provide reference key.\n";
						cout << "   mu_a  [mat key] [key] : return absorption coefficient for material.\n";
                        cout << "                           If more than one table is available, you must provide reference key.\n";
						cout << "   index [mat key] [key] : return background index for material.\n";
                        cout << "                           If more than one table is available, you must provide reference key.\n";
						cout << "   g     [mat key] [key] : return anisotropy coefficient for material.\n";
                        cout << "                           If more than one table is available, you must provide reference key.\n";
						cout << "   out   [filename]      : set output filename for tables.  By itself, closes the file and return output to screen.\n";
                    } else if (tokens[0] == "out") {
                        if (outStream != &cout)
                            delete outStream;
                        
                        if (tokens.size() == 2) {
                            outStream = new ofstream(tokens[1].c_str(),ios_base::app);
                        } else {
                            outStream = &cout;
                        }
                    } else if (tokens[0] == "mu_a") {
                        if (tokens.size() == 2) {
                            MCMaterial* mat = materialList.GetMaterialByKey(tokens[1]);

                            if (mat != NULL) {
                                fastinterpolate table = mat->GetPropertyDefaultTable(kAbsorptionCoefficient);
                                table.tostream(*outStream);
                            }
                        } else if (tokens.size() == 3) {
                            MCMaterial* mat = materialList.GetMaterialByKey(tokens[1]);
                            if (mat != NULL) {
                                fastinterpolate table = mat->GetPropertyTable(kAbsorptionCoefficient, tokens[2]);
                                table.tostream(*outStream);
                            }
                        }
                    } else if (tokens[0] == "mu_s") {
                        if (tokens.size() == 2) {
                            MCMaterial* mat = materialList.GetMaterialByKey(tokens[1]);
                            
                            if (mat != NULL) {
                                fastinterpolate table = mat->GetPropertyDefaultTable(kScatteringCoefficient);
                                table.tostream(*outStream);
                            }
                        } else if (tokens.size() == 3) {
                            MCMaterial* mat = materialList.GetMaterialByKey(tokens[1]);
                            
                            if (mat != NULL) {
                                fastinterpolate table = mat->GetPropertyTable(kScatteringCoefficient, tokens[2]);
                                table.tostream(*outStream);
                            }
                        }
                    } else if (tokens[0] == "index") {
                        if (tokens.size() == 2) {
                            MCMaterial* mat = materialList.GetMaterialByKey(tokens[1]);
                            
                            if (mat != NULL) {
                                fastinterpolate table = mat->GetPropertyDefaultTable(kBackgroundIndex);
                                table.tostream(*outStream);
                            }
                        } else if (tokens.size() == 3) {
                            MCMaterial* mat = materialList.GetMaterialByKey(tokens[1]);

                            if (mat != NULL) {
                                fastinterpolate table = mat->GetPropertyTable(kBackgroundIndex, tokens[2]);
                                table.tostream(*outStream);
                            }
                        }
                    } else if (tokens[0] == "g") {
                        if (tokens.size() == 2) {
                            MCMaterial* mat = materialList.GetMaterialByKey(tokens[1]);
                            
                            if (mat != NULL) {
                                fastinterpolate table = mat->GetPropertyDefaultTable(kAnisotropyCoefficient);
                                table.tostream(*outStream);
                            }
                        } else if (tokens.size() == 3) {
                            MCMaterial* mat = materialList.GetMaterialByKey(tokens[1]);
                            
                            if (mat != NULL) {
                                fastinterpolate table = mat->GetPropertyTable(kAnisotropyCoefficient, tokens[2]);
                                table.tostream(*outStream);
                            }
                        }
                    } else if (tokens[0] == "info") {
                        if (tokens.size() == 1) {
                            vector<string> list = materialList.GetMaterialsListByKey();
                            
                            cout << "\nAvailable materials:\n";
                            for (unsigned long i = 0; i < list.size(); i++) {
                                cout << "  " << list[i] << endl;
                            }
                            
                        } else if (tokens.size() >= 2) {
                            MCMaterial* mat = materialList.GetMaterialByKey(tokens[1]);
                            
                            if (mat != NULL) {
                                cout << "Stats:\n";
                                vector<DataSource*> source = mat->GetPropertyDataSource(kBackgroundIndex);
                                cout << "  Index: " << source.size() << " tables (";
								string t = ", ";
                                for (long i = 0; i < source.size(); i++) {
									if (i == source.size() - 1)
										t = "";
									cout << source[i]->key << t;
								}
                                cout << ")\n";
                                
                                source = mat->GetPropertyDataSource(kScatteringCoefficient);
                                cout << "  Scattering: " << source.size() << " tables (";
								t = ", ";
                                for (long i = 0; i < source.size(); i++) {
									if (i == source.size() - 1)
										t = "";
									cout << source[i]->key << t;
								}
                                cout << ")\n";
                                source = mat->GetPropertyDataSource(kAbsorptionCoefficient);
                                cout << "  Absorption: " << source.size() << " tables (";
								t = ", ";
                                for (long i = 0; i < source.size(); i++) {
									if (i == source.size() - 1)
										t = "";
									cout << source[i]->key << t;
								}
                                cout << ")\n";
                                source = mat->GetPropertyDataSource(kAnisotropyCoefficient);
                                cout << "  Anisotropy: " << source.size() << " tables (";
								t = ", ";
                                for (long i = 0; i < source.size(); i++) {
									if (i == source.size() - 1)
										t = "";
									cout << source[i]->key << t;
								}
                                cout << ")\n";
                            }
                        }
                    } else if (tokens[0] == "quit") {
                        if (outStream != &cout)
                            delete outStream;
                        
                        return 1;
                    }
                }

            } catch (exception& e) {
                clog << RightNow() << " Error: " << string(e.what()) << endl;
            }
        
        }
        
    } catch (exception& e) {
        clog << RightNow() << "Exception caught: " << string(e.what()) << endl;
        clog << RightNow() << "Program ended abnormally" << endl;
        return 1;
    }
    
    return 0;
}
#endif
