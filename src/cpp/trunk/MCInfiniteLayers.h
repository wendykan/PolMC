#ifndef __MCINFINITELAYERS
#define __MCINFINITELAYERS

#include "MCObject.h"

class MCInfiniteLayer : public MCObject {

public:

    MCInfiniteLayer(double inThickness,
                    double inImgWidth,
                    double inImgHeight,
                    long inPtsWidth,
                    long inPtsHeight,
                    double inIndexMedium,
                    double inIndexOutside,
                    double inRotationPerCmInClearSpace,
					long inDetectionType,
                    long inAcceptanceCosineElements,
                    bool inUseNextEventEstimator);
    MCInfiniteLayer(map<string,string>& inDict):MCObject(inDict) { init(inDict); }

    void init(double inThickness, double inImgWidth, double inImgHeight, long inPtsWidth, long inPtsHeight, long inDetectionType, long inAcceptanceCosineElements);
    void init(map<string,string>& inDict);
    virtual bool IsOutsideObject(RealV& inLocalPosition);
	virtual bool IsInsideObject(RealV& inLocalPosition);

    void	InterfaceTouchesOtherObject(MCObject * inOtherObject, long inOtherInterface, long thisObjectInterface);
    
protected:
    double mThickness;
    
};

#endif





