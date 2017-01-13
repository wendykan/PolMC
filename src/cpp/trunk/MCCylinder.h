#ifndef MCCYLINDER_H
#define MCCYLINDER_H

#include "MCObject.h"

class MCCylinder : public MCObject {
public:
	MCCylinder(double inRadius,
			   double inHeight,
			   long inPtsRadius,
			   long inPtsHeight,
			   double inIndexMedium,
			   double inIndexOutside,
			   double inRotationPerCmInClearSpace,
			   long inDetectionType,
			   long inAcceptanceCosineElements,
			   bool inUseNextEventEstimator);
	
  void init(double inRadius,
			double inHeight,
			long inPtsRadius,
			long inPtsHeight,
			long inDetectionType,
			long inAcceptanceCosineElements);

  void init(map<string,string>& inDict);

  bool IsOutsideObject(RealV& inLocalPosition);

protected:
	double mRadius;
	double mHeight;
};

#endif
