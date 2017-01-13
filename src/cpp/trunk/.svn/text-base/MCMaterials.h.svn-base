
#ifndef MATERIALS_H
#define MATERIALS_H

#include "fastinterpolate.h"

#if HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/xpath.h>
#endif

typedef struct { fastinterpolate table; string description, key; } DataSource;

enum {
    kBackgroundIndex,
    kScatteringCoefficient,
    kAbsorptionCoefficient,
    kAnisotropyCoefficient};
    
class MCMaterial {
public:
     MCMaterial();
    virtual ~MCMaterial();
    MCMaterial& operator=(const MCMaterial& inMaterial);
    
    bool initFromXMLDocument(xmlXPathContextPtr inXPath);
    bool initDefaultsFromXMLDocument(xmlXPathContextPtr inXPath);
    bool initTablesFromNodeSet(xmlXPathContextPtr inContext, char* inXPath, vector<DataSource*>& inDataSource);

    string GetName();
    string GetKey();
    double GetPropertyValue(long inProperty, double inWavelength);
    fastinterpolate GetPropertyDefaultTable(long inProperty);
    fastinterpolate GetPropertyTable(long inProperty, string inKey);

    vector<DataSource*>  GetPropertyDataSource(long inProperty);

    double GetBackgroundIndex(double inWavelength);
    double GetScatteringCoefficient(double inWavelength);
    double GetAbsorptionCoefficient(double inWavelength);
    double GetAnisotropyCoefficient(double inWavelength);
    void SetRelativeConcentration(double inConcentration);
    void SetDefaultKey(string inKey);
    void SetDefaultKeyForScatteringCoefficient(string inKey);
    void SetDefaultKeyForAbsorptionCoefficient(string inKey);
    void SetDefaultKeyForBackgroundIndex(string inKey);
    void SetDefaultKeyForAnisotropyCoefficient(string inKey);
        
protected:
    string mName;
    string mKey;
    string mDescription;
    double mConcentration;
    
    DataSource* mDefaultBackgroundIndex;
    DataSource* mDefaultScatteringCoefficient;
    DataSource* mDefaultAbsorptionCoefficient;
    DataSource* mDefaultAnisotropyCoefficient;
    
    vector<DataSource*> mBackgroundIndexTables;
    vector<DataSource*> mScatteringCoefficientTables;
    vector<DataSource*> mAbsorptionCoefficientTables;
    vector<DataSource*> mAnisotropyCoefficientTables;
    
};

class MCMaterialsList {
    
public:
    MCMaterialsList(string inFilename);
	MCMaterialsList(string inFilename, string inDefaultsFilename);

    virtual ~MCMaterialsList();

    MCMaterial* GetMaterialByName(string inName);
    MCMaterial* GetMaterialByKey(string inKey);

    vector<string> GetMaterialsListByKey();
    vector<string> GetMaterialsListByName();

    bool initFromXMLDocument(string inFilename);
    bool initDefaultsFromXMLDocument(string inFilename);

protected:
        vector <MCMaterial*> mMaterials;
};

#endif
