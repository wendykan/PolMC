
#include "Photon.h"
#include "RealV.h"

#ifndef MCSOURCE_H
#define MCSOURCE_H

class MCSource {
public:
	
	MCSource(RealV inGlobalOrigin, RealV inDirection, double inWavelength);
	virtual ~MCSource();
	virtual void	InitializePhoton(Photon* outPhoton);
	void SetTemplatePhoton(Photon* inPhoton);
	Photon* GetTemplatePhoton();
	Photon* GetNewPhoton();
	void	SetMainAxis(RealV inDirection);
	void	SetGlobalOrigin(RealV inGlobalOrigin);
	RealV	GetGlobalOrigin();
	void	GetRandomPointOnDisk(double inRadius, double& outX, double& outY);
	void	GetRandomPointOn2DGaussian(double inRadius, double& outX, double& outY);
	double  GetRandomPolarAngle();
	double  GetRandomAzimuthalAngle();
	void	SetContainerObject(MCObject* inContainerObject);
	long	GetNumberOfPhotonsToLaunch();
	void	SetNumberOfPhotonsToLaunch(long nPhoton);

    long    GetNumberOfPhotonsLaunched();
    void    AddToStats(RealV& inPosition, RealV& inDirection);
	void    WriteStatsToStream(ostream& out);
    void    KeepStats();
    void    DontKeepStats();
	virtual void	DumpToStream(ostream & out);

protected:
	MCWorld*	mWorld;
	Photon*		templatePhoton;
    MCObject*   mContainerObject;
    double      mWavelength;
	RealV		mDirection;
	RealV		mGlobalOrigin;
    long        mNumberOfPhotonsToLaunch;
    long        mNumberOfPhotonsLaunched;
    list<RealV> mPosStats;
    list<RealV> mDirStats;
    bool        mKeepStats;
};


class MCLaserSource : public MCSource {
public:
	MCLaserSource(StokesV inStokes, double inRadius, double inRadiusOfCurvature, double inNA, RealV inGlobalOrigin, RealV inDirection, double inWavelength);
	virtual ~MCLaserSource();
    
	void	SetStokesVector(StokesV inStokes);
	StokesV GetStokesVector();
	void	InitializePhoton(Photon* outPhoton);
protected:
		StokesV mStokes;
	double mBeamRadius;
	double mRadiusOfCurvature;
	double mNA;
};

class MCIsotropicPointSource : public MCSource {
public:
	MCIsotropicPointSource(RealV inGlobalOrigin, double inWavelength);
    virtual ~MCIsotropicPointSource();
	void	InitializePhoton(Photon* outPhoton);

};

class MCCylindricalDiffuser : public MCSource {
public: 
	MCCylindricalDiffuser(double inLength,
						  RealV inGlobalOrigin, 
						  RealV inDirection,
                          double inWavelength);
	MCCylindricalDiffuser(string inZFilename,
                          string inThetaFilename,
                          string inPhiFilename,
						  RealV inGlobalOrigin, 
						  RealV inDirection,
                          double inWavelength);    
	MCCylindricalDiffuser(double* inZArray,
                          double* inNormalizedZDistribution,
                          long inZPts,
                          double* inThetaArray,
                          double* inNormalizedThetaDistribution, 
                          long inThetaPts,
                          double* inPhiArray,
                          double* inNormalizedPhiDistribution,
                          long inPhiPts,
                          RealV inGlobalOrigin,
                          RealV inDirection,
                          double inWavelength);
	
	virtual ~MCCylindricalDiffuser();
	void	InitializePhoton(Photon* outPhoton);
protected:
    double mLength;
	long mZPts;
	long mThetaPts;
	long mPhiPts;
	double *mNormalizedZDistribution;
	double *mNormalizedThetaDistribution;
	double *mNormalizedPhiDistribution;
    
    fastinterpolate mZDistribution;
    fastinterpolate mThetaDistribution;
    fastinterpolate mPhiDistribution;
    
	
};
#endif
