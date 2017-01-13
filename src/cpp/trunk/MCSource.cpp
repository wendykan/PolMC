#include "MCSource.h"
#include "MCObject.h"
#include "Photon.h"
#include "RealV.h"
#include "MCRandomScatterer.h"

MCSource::MCSource(RealV inGlobalOrigin, RealV inDirection, double inWavelength)
{

	SetMainAxis(inDirection);
	SetGlobalOrigin(inGlobalOrigin);
    mWavelength = inWavelength;
	mContainerObject = 0;
	mWorld = 0;
    mNumberOfPhotonsLaunched = 0;
    templatePhoton = 0;
    mKeepStats = false;
}

MCSource::~MCSource()
{
    
}

void
MCSource::SetTemplatePhoton(Photon* inPhoton)
{
	templatePhoton = inPhoton->MakeCopy();
}

Photon*
MCSource::GetTemplatePhoton()
{
	return templatePhoton;
}

Photon*
MCSource::GetNewPhoton()
{
	Photon* newPhoton = NULL;
	if ( templatePhoton != NULL ) {
		newPhoton = templatePhoton->MakeCopy();
		InitializePhoton(newPhoton);
	} else
		throw runtime_error("Photon type not set in light source (template is NULL)");
	
	return newPhoton;
}


void
MCSource::InitializePhoton(Photon* outPhoton)
{
	throw logic_error("MCSource cannot be instantiated: use subclasses only");
}

long
MCSource::GetNumberOfPhotonsToLaunch()
{
    return mNumberOfPhotonsToLaunch;
}

void
MCSource::SetNumberOfPhotonsToLaunch(long nPhoton)
{
    mNumberOfPhotonsToLaunch = nPhoton;
}

long
MCSource::GetNumberOfPhotonsLaunched()
{
    return mNumberOfPhotonsLaunched;
}

void
MCSource::SetMainAxis(RealV inDirection)
{
	mDirection = inDirection;
	mDirection.normalize();
}

void
MCSource::SetGlobalOrigin(RealV inGlobalOrigin)
{
	mGlobalOrigin = inGlobalOrigin;
}

RealV
MCSource::GetGlobalOrigin()
{
    return mGlobalOrigin;
}


void 
MCSource::GetRandomPointOnDisk(double inRadius, double& outX, double& outY)
{
	double r = inRadius;
	do {
		outX = (MCRandomScatterer::RandomFloat()-0.5) * 2 * r;
		outY = (MCRandomScatterer::RandomFloat()-0.5) * 2 * r;
	} while (outY * outY + outX * outX > r*r);

}

void 
MCSource::GetRandomPointOn2DGaussian(double inRadius, double& outX, double& outY)
{
	throw runtime_error("Gaussian distribution not implemented yet in sources");
}

double 
MCSource::GetRandomPolarAngle()
{
	return acos(1. - 2. * MCRandomScatterer::RandomFloat());
}

double 
MCSource::GetRandomAzimuthalAngle()
{
	return 2. * PI * MCRandomScatterer::RandomFloat();
}

void
MCSource::SetContainerObject(MCObject* inContainerObject)
{
	MyAssert_(inContainerObject != NULL);

	mContainerObject = inContainerObject;
	mWorld = mContainerObject->GetWorld();
	
}

void
MCSource::AddToStats(RealV& inPosition, RealV& inDirection)
{
    mPosStats.push_back(inPosition);
    mDirStats.push_back(inDirection);
}

void
MCSource::WriteStatsToStream(ostream& out)
{
    RealV pos, dir;

    out << "pos\tdir" << endl;

    while ( mPosStats.size() != 0 ) {
        pos = mPosStats.front();
        dir = mDirStats.front();
        
        mPosStats.pop_front();
        mDirStats.pop_front();

        out << pos << "\t" << dir << "\t" << endl;
    }
    
}

void
MCSource::KeepStats()
{
    mKeepStats = true;
}

void
MCSource::DontKeepStats()
{
    mKeepStats = false;
}


void
MCSource::DumpToStream(ostream & out)
{
	out << "<source>\n";
    out << "<origin coordinate=\"global\">" << mGlobalOrigin<< "</origin>" << endl;
    out << "<direction>" << mDirection << "</direction>" << endl;
	
	if ( mContainerObject != NULL ) {
		out << "<insideObject>" << mContainerObject->GetName()<< "</insideObject>" << endl;
	} else {
		out << "<insideObject>NULL</insideObject>" << endl;
	}
    out << "<wavelength>" << mWavelength<< "</wavelength>" << endl;

	if (mKeepStats) {
		out << "<stats>";
		WriteStatsToStream(out);
		out << "</stats>";
	}
	out << "</source>\n";
	
}

MCLaserSource::MCLaserSource(StokesV inStokes, double inBeamRadius, double inRadiusOfCurvature, double inNA, RealV inGlobalOrigin, RealV inDirection, double inWavelength)
	:MCSource(inGlobalOrigin, inDirection, inWavelength)
{
	mStokes = inStokes;
    MyAssert_(inDirection == inStokes.GetPropagationDirection());

	mBeamRadius = inBeamRadius;
	mRadiusOfCurvature = inRadiusOfCurvature;
	mNA = inNA;

    if (mBeamRadius != 0 && (inDirection != RealV(0,0,1) && inDirection != RealV(0,0,-1)))
        throw logic_error("Unimplemented: the laser source does not accept disk beams if the direction is not along z");        

    if (mNA != 0) 
        throw logic_error("Unimplemented: the laser source does not accept laser beams with NA different from 0.  Modify MCLaserSource::InitializePhoton() ");        
}

MCLaserSource::~MCLaserSource()
{
    
}

void
MCLaserSource::SetStokesVector(StokesV inStokes)
{
	MyAssert_(mDirection == inStokes.GetPropagationDirection());
	mStokes = inStokes;
}

StokesV
MCLaserSource::GetStokesVector()
{
	return mStokes;
}

void
MCLaserSource::InitializePhoton(Photon* outPhoton)
{
    RealV pos, dir;
    
    if (mContainerObject)
        outPhoton->SetCurrentObject(mContainerObject);
    else
        throw logic_error("Object not set properly in MCSource.  You must call SetContainerObject() for your source to explicitly indicate in what object it is contained. ");
    
    double x = 0,y = 0;
    if (mBeamRadius != 0 ) {
		if ( mStokes.GetPropagationDirection() == RealV(0,0,1) ) {
			GetRandomPointOnDisk(mBeamRadius, x, y);
			pos = mGlobalOrigin + RealV(x, y, 0);
		} else {
			throw runtime_error("Beam radius not zero: direction must be (0,0,1) or you must fix the code");
		}
    } else {
        pos = mGlobalOrigin;
    }

    outPhoton->SetLocalPosition(pos);
    outPhoton->SetWavelength(mWavelength);

	
	if (mRadiusOfCurvature != 0) {
		// SetStokesVector has direction in it.
		StokesV copy = mStokes;
		dir = RealV(0, 0, mRadiusOfCurvature ) - RealV(x, y, 0);
		dir.normalize();
		copy.SetPropagationDirectionInLabFrame(dir);
		outPhoton->SetStokesVector(copy);
	} else {
		dir = mStokes.GetPropagationDirection();
		outPhoton->SetStokesVector(mStokes);
	}
	
    if (mKeepStats)
        AddToStats(pos, dir);

    mNumberOfPhotonsLaunched++;
}

MCIsotropicPointSource::MCIsotropicPointSource(RealV inGlobalOrigin, double inWavelength)
:MCSource(inGlobalOrigin, RealV(0,0,0), inWavelength)
{
    
}

MCIsotropicPointSource::~MCIsotropicPointSource()
{
    
}

void
MCIsotropicPointSource::InitializePhoton(Photon* outPhoton)
{
    RealV dir, pos;
    
    if (mContainerObject)
        outPhoton->SetCurrentObject(mContainerObject);
    else
        throw logic_error("Object not set properly in MCSource.  You must call SetContainerObject() for your source to explicitly indicate in what object it is contained. ");

    double phi, theta;
    phi = GetRandomAzimuthalAngle();
    theta = GetRandomPolarAngle();
    
    dir = RealV(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
    pos = mGlobalOrigin - mContainerObject->GetGlobalOrigin();
    
    outPhoton->SetCurrentObject(mContainerObject);
    outPhoton->SetLocalPosition(pos);
	StokesV s(1,dir);
    outPhoton->SetStokesVector(s);
//    outPhoton->SetPropagationDirectionInLabFrame(dir); // not necessary
    outPhoton->SetWavelength(mWavelength);
	
    if (mKeepStats)
        AddToStats(pos, dir);
    
    mNumberOfPhotonsLaunched++;
}

MCCylindricalDiffuser::MCCylindricalDiffuser(double inLength,
									RealV inGlobalOrigin, 
									RealV inDirection, double inWavelength)
									: MCSource(inGlobalOrigin, inDirection, inWavelength)
{
	mLength = inLength;

	mZPts = 0;
	mThetaPts = 0;
	mPhiPts = 0;
	
	mNormalizedZDistribution = 0;
	mNormalizedThetaDistribution = 0;
	mNormalizedPhiDistribution = 0;

}

MCCylindricalDiffuser::MCCylindricalDiffuser(string inZFilename,
                      string inThetaFilename,
                      string inPhiFilename,
                      RealV inGlobalOrigin, 
                      RealV inDirection, double inWavelength):
MCSource(inGlobalOrigin, inDirection, inWavelength),
mZDistribution(inZFilename),mThetaDistribution(inThetaFilename),mPhiDistribution(inPhiFilename)
{
    mLength = 0;   
}


MCCylindricalDiffuser::MCCylindricalDiffuser(double* inZArray,
                                             double* inNormalizedZDistribution,
                                             long inZPts,
                                             double* inThetaArray,
                                             double* inNormalizedThetaDistribution, 
                                             long inThetaPts,
                                             double* inPhiArray,
                                             double* inNormalizedPhiDistribution,
                                             long inPhiPts,
                                             RealV inGlobalOrigin,
                                             RealV inDirection, double inWavelength)
: MCSource(inGlobalOrigin, inDirection, inWavelength)
{
    throw logic_error("Reimplement: insert points into fastinterpolate arrays");
}

MCCylindricalDiffuser::~MCCylindricalDiffuser()
{
	
}

void
MCCylindricalDiffuser::InitializePhoton(Photon* outPhoton)
{
    RealV dir,pos;
    
	if (mLength != 0) {
		// Perfect diffuser
		double theta = GetRandomPolarAngle();
		double phi = GetRandomAzimuthalAngle();
		
		dir = RealV(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
		pos = mGlobalOrigin + mDirection * MCRandomScatterer::RandomFloat() * mLength ;
		
	} else {
		// Arbitrary diffuser
		double z,theta,phi;
		float y;
		//distribution for theta, rejection method
		bool done=false;
        long i = 1000; // Attempts
		while (--i && !done) {
			theta = 180. * MCRandomScatterer::RandomFloat();
			y = MCRandomScatterer::RandomFloat() * mThetaDistribution.yMaximum();
            
			if( y < mThetaDistribution.y(theta)) {
				theta *= PI / 180 ;
				done=true;
            }
        }
		if (i == 0) {
            throw runtime_error("Could not determine theta position of photon on diffuser");
        }
        
        
		phi = GetRandomAzimuthalAngle();
		//distribution for z, rejection method
	    done=false;
	    i = 1000;
		while (--i && !done) {
			z = mZDistribution.xMaximum() * MCRandomScatterer::RandomFloat();
			y = MCRandomScatterer::RandomFloat() * mZDistribution.yMaximum();

			if( y <= mZDistribution.y(z) ) {
				done=true;
            }
        }

        if (i == 0) {
            throw runtime_error("Could not determine z position of photon on diffuser");
        }
        
      	dir = RealV(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
		pos = mGlobalOrigin + mDirection * z;
	}
 
    outPhoton->SetCurrentObject(mContainerObject);
    outPhoton->SetGlobalPosition(pos);
    outPhoton->SetPropagationDirectionInLabFrame(dir);
    outPhoton->SetWavelength(mWavelength);

    if (mContainerObject->IsOutsideBoundingBox(outPhoton)) {
        RealV top, bottom;
        mContainerObject->GetBoundingBox(top, bottom);
        
        clog << "Initial position of photon appears to be outside of container object: check 'z' file and size of container object.";
        clog << "Bounding box (top,bottom): " << top << " and " << bottom << ". Photon global position: " << pos << endl;
    }

    if (mKeepStats)
        AddToStats(pos, dir);

    mNumberOfPhotonsLaunched++;
    
    
}

