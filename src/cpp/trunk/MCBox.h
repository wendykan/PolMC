#ifndef __MCBOX
#define __MCBOX

#include "MCObject.h"

class MCBox : public MCObject {

public:
    MCBox(double inWidth,
          double inHeight,
          double inDepth,
          long inPtsWidth,
          long inPtsHeight,
          long inPtsDepth,
          double inIndexMedium,
          double inIndexOutside,
          double inRotationPerCmInClearSpace,
		  long inDetectionType,
          long inAcceptanceCosineElements,
          bool inUseNextEventEstimator);
    
    void init(double inWidth,
            double inHeight,
            double inDepth,
            long inPtsWidth,
            long inPtsHeight,
			long inPtsDepth,
			long inDetectionType,
			long inAcceptanceCosineElements);

    MCBox(map<string,string>& inDict);

    void init(map<string,string>& inDict);
    virtual bool	IsOutsideObject(RealV& inLocalPosition); 
	virtual bool	IsInsideObject(RealV& inLocalPosition);

protected:
    double mWidth;
    double mHeight;
    double mDepth;
};

#endif
