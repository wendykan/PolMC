#ifndef __MCELLIPSOID
#define __MCELLIPSOID

#include "MCObject.h"

class MCEllipsoid : public MCObject {

public:
    MCEllipsoid(double inWidth,
                double inHeight,
                double inDepth,
                long inPtsWidth,
                long inPtsHeight,
                double inIndexMedium,
                double inIndexOutside,
                double inRotationPerCmInClearSpace,
				long inDetectionType,
                long inAcceptanceCosineElements,
                bool inUseNextEventEstimator);
    
    MCEllipsoid(map<string,string>& inDict):MCObject(inDict) { init(inDict); }

    void init(double inWidth,
              double inHeight,
              double inDepth,
              long inPtsWidth,
              long inPtsHeight,
			  long inDetectionType,
			  long inAcceptanceCosineElements);

        void init(map<string,string>& inDict);
    virtual bool	IsOutsideObject(RealV& inLocalPosition); 

    virtual RealV    GetRandomPointInsideObjectUniformDistribution();

protected:
    double mWidth;
    double mHeight;
    double mDepth;
    
    double minNormForComplicatedCalculation;

};

#endif

