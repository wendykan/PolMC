#include "StokesV.h"
#include "MCRandomScatterer.h"
#include "cubicspline.h"

bool StokesV::coherentSum = false;

fastinterpolate gCos;
fastinterpolate gSin;

#define sin_(x) (sin( x ) )
#define cos_(x) (cos( x ) )

StokesV::StokesV():
mEpara(-1,0,0),mEperp(0,1,0),mEprop(0,0,1)
{   
    mI = 0;
    mQ = 0;
    mU = 0;
    mV = 0;
	
	mTotalPropagationPhase = 0;
	
	
	if (gCos.size() == 0) {
		for (double x = -4 * PI; x < 4 * PI; x += 8 * PI / 100.) {
			gCos.insertpoint(x, cos(x));
			gSin.insertpoint(x, sin(x));
		}
		gCos.CheckEquallySpaced();
		gSin.CheckEquallySpaced();
		
		//gCos.initSplineTables();
		//gSin.initSplineTables();

	}
	
}


StokesV::StokesV(double inI, double inQ, double inU, double inV, 
                 RealV inEpara , RealV inEprop)
:mEpara(inEpara),mEperp(RealV::CrossProduct(inEpara, inEprop)),mEprop(inEprop)
{   
    mI = inI;
    mQ = inQ;
    mU = inU;
    mV = inV;

	mTotalPropagationPhase = 0;
}


StokesV::StokesV(double inI, RealV inEprop)
{   
    mI = inI;
    mQ = 0;
    mU = 0;
    mV = 0;

	mTotalPropagationPhase = 0;

    SetPropagationDirectionInLabFrame(inEprop);
}


StokesV::StokesV(Complex inEl, Complex inEr)
:mEpara(-1,0,0),mEperp(0,1,0),mEprop(0,0,1)
{   
    mI = norm(inEl) + norm(inEr);
    mQ = norm(inEl) - norm(inEr);
    
    MyAssert_(imag(conj(inEl)*inEr + inEl*conj(inEr)) < 1e-3);
    MyAssert_(imag(I*(conj(inEl)*inEr - inEl*conj(inEr))) < 1e-3);
    
    mU = real(inEl * conj(inEr)  + conj(inEl) * inEr);
    mV = real(I*(inEl * conj(inEr)  + conj(inEl) * inEr));

	mTotalPropagationPhase = 0;

}


StokesV::StokesV(const StokesV& inS)
{       
    mI = inS.mI;
    mQ = inS.mQ;
    mU = inS.mU;
    mV = inS.mV;

	mTotalPropagationPhase = inS.mTotalPropagationPhase;

    mEpara = inS.mEpara;
    mEperp = inS.mEperp;
    mEprop = inS.mEprop;
}

StokesV::~StokesV()
{
	
}

void 
StokesV::SetCoherentSum(bool inCoherent)
{
	coherentSum = inCoherent;
}

bool 
StokesV::GetCoherentSum()
{
	return coherentSum;
}


void
StokesV::Normalize()
{   
    if (mI != 0) {
        mQ /= mI;
        mU /= mI;
        mV /= mI;           
        mI = 1.;
    } else {
        throw runtime_error("Intensity is zero during normalization");
    }
}

StokesV&
StokesV::operator=(const StokesV& inRhs)
{   
    mI = inRhs.mI;
    mQ = inRhs.mQ;
    mU = inRhs.mU;
    mV = inRhs.mV;
    
	mTotalPropagationPhase = inRhs.mTotalPropagationPhase;

    mEpara = inRhs.mEpara;
    mEperp = inRhs.mEperp;
    mEprop = inRhs.mEprop;
    
    CheckDoubleValue_(mI);
    CheckDoubleValue_(mQ);
    CheckDoubleValue_(mU);
    CheckDoubleValue_(mV);
    
    return *this;
}

StokesV&
StokesV::operator+=(StokesV inRhs)
{             
	if ( ! coherentSum ) {
		mI += inRhs.mI;
		mQ += inRhs.mQ;
		mU += inRhs.mU;
		mV += inRhs.mV;
	} else {
		Complex El, Er, rhsEr, rhsEl;
		
		GetLocalComplexFields(El, Er);
		inRhs.GetLocalComplexFields(rhsEl, rhsEr);
		
		El += rhsEl;
		Er += rhsEr;

	    mI = norm(El) + norm(Er);
		mQ = norm(El) - norm(Er);
		
		MyAssert_(imag(conj(El)*Er + El*conj(Er)) < 1e-3);
		MyAssert_(imag(I*(conj(El)*Er - El*conj(Er))) < 1e-3);
		
		mU = real(El * conj(Er)  + conj(El) * Er);
		mV = real(I*(El * conj(Er)  + conj(El) * Er));
		
		mTotalPropagationPhase = 0;
		
	}

    CheckDoubleValue_(mI);
    CheckDoubleValue_(mQ);
    CheckDoubleValue_(mU);
    CheckDoubleValue_(mV);
    
    return *this;
}

StokesV&
StokesV::operator*=(double inRhs)
{   
    mI *= inRhs;
    mQ *= inRhs;
    mU *= inRhs;
    mV *= inRhs;
    
    return *this;
}

StokesV&
StokesV::operator/=(double& inRhs)
{       
    mI /= inRhs;
    mQ /= inRhs;
    mU /= inRhs;
    mV /= inRhs;
    
    return *this;
}

bool
StokesV::operator==(StokesV inRhs)
{   
    if (mI != inRhs.mI)
        return false;
    else if (mQ == inRhs.mQ)
        return false;
    else if (mU == inRhs.mU)
        return false;
    else if (mV == inRhs.mV)
        return false;
    
    return true;
}

bool
StokesV::operator!=(StokesV inRhs)
{   
    return !(*this == inRhs);
}

void
StokesV::GetReferenceFrame(RealV& outEperp, RealV& outEpara, RealV& outEprop)
{   
    outEperp = mEperp;
    outEpara = mEpara;
    outEprop = mEprop;
}

void
StokesV::SetReferenceFrame( RealV& inEperp, RealV& inEpara, RealV& inEprop)
{   
    mEperp = inEperp;
    mEpara = inEpara;
    mEprop = inEprop;
    
    CheckTriad_(mEperp, mEpara, mEprop);
    
}

RealV
StokesV::GetEparaAxis()
{   
    return mEpara;
}

RealV
StokesV::GetEperpAxis()
{
    return mEperp;
}

RealV
StokesV::GetPropagationDirection()
{
    return mEprop;
}

void
StokesV::SetPropagationDirectionInLabFrame(RealV inDirection)
{
    if (inDirection.z != 0)
        mEperp = RealV(1,0, -inDirection.x/inDirection.z);
    else if (inDirection.y != 0)
        mEperp = RealV(1,-inDirection.x/inDirection.y, 0);
    else
        mEperp = RealV(0,1,0);
    
    mEperp.normalize();
    
    mEpara = RealV::NormalizedCrossProduct(inDirection, mEperp);
    
    mEprop = inDirection;
}

void
StokesV::RotateReferenceFrameAroundPropagationDirectionBy(double inPhi)
{
    CheckDoubleValue_(inPhi);
    
    RealV el(mEpara), er(mEperp),ez(mEprop);
    
    double cos_phi = cos_(inPhi);
    double sin_phi = sin_(inPhi);
    
    CheckTriad_(mEperp, mEpara, mEprop);
    
    mEperp.x = er.x * cos_phi + el.x * sin_phi;
    mEperp.y = er.y * cos_phi + el.y * sin_phi;
    mEperp.z = er.z * cos_phi + el.z * sin_phi;
    
    mEpara.x = - er.x * sin_phi + el.x * cos_phi;
    mEpara.y = - er.y * sin_phi + el.y * cos_phi;
    mEpara.z = - er.z * sin_phi + el.z * cos_phi;
    
    CheckTriad_(mEperp, mEpara, mEprop);
    
    // The new Stokes vector could be recalculated
    // from scratch but this is too long (and stupid):
    // the new Stokes vector is simply a rotated version
    // of the previous one that has followed 
    RotatePolarizationStateBy(-inPhi);
    
}

void
StokesV::RotatePolarizationStateBy(double inPhi)
{
    // MyAssert_(inPhi != 0);
    CheckDoubleValue_(inPhi);
    
    double cos_2phi = cos_(2. * inPhi);
    double sin_2phi = sin_(2. * inPhi);
    
    RotatePolarizationStateBy(cos_2phi, sin_2phi);
}

void
StokesV::RotatePolarizationStateBy(double inCos_2phi, double inSin_2phi)
{
    double q = mQ;
    
    CheckDoubleValue_(inCos_2phi);
    CheckDoubleValue_(inSin_2phi);
    
    mQ = inCos_2phi * mQ - inSin_2phi * mU ;
    mU = inSin_2phi * q  + inCos_2phi * mU ;
    
    CheckDoubleValue_(mQ);
    CheckDoubleValue_(mU);
    
}

void
StokesV::ChangePropagationDirectionAroundEPerpBy(double inTheta)
{
    RealV el(mEpara), er(mEperp),ez(mEprop);
    
    CheckTriad_(mEperp, mEpara, mEprop);
    
    CheckDoubleValue_(inTheta);
    
    double cos_theta = cos_(inTheta);
    double sin_theta = sin_(inTheta);
    
    mEpara.x = el.x * cos_theta + ez.x * sin_theta;
    mEpara.y = el.y * cos_theta + ez.y * sin_theta;
    mEpara.z = el.z * cos_theta + ez.z * sin_theta;
    
    mEprop.x = - el.x * sin_theta + ez.x * cos_theta;
    mEprop.y = - el.y * sin_theta + ez.y * cos_theta;
    mEprop.z = - el.z * sin_theta + ez.z * cos_theta;
    
    CheckTriad_(mEperp, mEpara, mEprop);
    
}

void
StokesV::GetLocalComplexFields(Complex& outEl, Complex& outEr)
{
    CheckDoubleValue_(mI);
    CheckDoubleValue_(mQ);
    CheckDoubleValue_(mU);
    CheckDoubleValue_(mV);
        
    outEl = sqrt(abs(mI+mQ)/2.);
    outEr = sqrt(abs(mI-mQ)/2.);
    CheckValue_(outEl);
    CheckValue_(outEr);
    
	Complex delta_theta;
	
	delta_theta = atan2(mV, mU);
	CheckValue_(delta_theta);
	CheckValue_(outEl);

	outEl *= exp(-I*delta_theta/2.);
	outEr *= exp(I*delta_theta/2.);

	if ( coherentSum ) {
		// Only relevant if we need the coherent sum of the fields
		if ( areTheSame(degreeOfPolarization(),1,4) ) {
			outEl *= exp(I*(mTotalPropagationPhase));
			outEr *= exp(I*(mTotalPropagationPhase));
		} else {
			PrintMessageIfLevel_( "Degree of polarization not 1 although coherent sum was requested (" << *this << ")", kVerbose);
			outEl = 0;
			outEr = 0;
		}
	}
#ifdef __MYDEBUG
    StokesV debug(outEl,outEr);
    
    //       MyAssert_(debug == *this);
    if (debug != *this) {
        //clog << RightNow() <<  "Conversion error in StokesV GetLocalComplexFields():\n" << debug << "!= " << *this << endl;
    }
#endif
    
}

void 
StokesV::TransmitThrough(RealV& inNormal, double rp, double rs, double tp, double ts)
{
	
}

void 
StokesV::ReflectAtInterface(RealV& inNormal, double rp, double rs, double tp, double ts)
{
    CheckDoubleValue_(rp);
    CheckDoubleValue_(rs);
    CheckDoubleValue_(tp);
    CheckDoubleValue_(ts);
    MyAssert_(rp * rp + rs * rs != 0);
    
    // First we modify the Stokes vector, then we modifiy the propagation direction
    // Need to get the reference frame the same
    double tempI;
	
    tempI = rp * rp / 2. * ( mI + mQ) + rs * rs / 2. * ( mI - mQ);
    mQ = rp * rp / 2. * ( mI + mQ) + rs * rs / 2. * ( - mI + mQ);
    mI = tempI;
	mU = rp * rs * mU ;
    mV = rp * rs * mV ;
    
	if (rp < 0 && rs < 0) {
		mTotalPropagationPhase += PI;
	}
	
    Normalize();
    
    CheckDoubleValue_(mI);
    CheckDoubleValue_(mQ);
    CheckDoubleValue_(mU);
    CheckDoubleValue_(mV);
    
    /* Compute angles */
    RealV el,er,e3;
    GetReferenceFrame(er,el,e3);
    PrintMessageIfLevel_("Reflection: Normal is " << inNormal << " Er (and s) is " << er << " E3 is " << e3, kExtremelyVerbose);
    /* The normal can be into or out of the plane.  We want the angle of incidence (must take abs()) */
    double thetaFrom = acos(abs(RealV::NormalizedDotProduct(e3, inNormal)));
    CheckDoubleValue_(thetaFrom);
    
    ChangePropagationDirectionAroundEPerpBy( - PI + 2. * thetaFrom);
//    PrintMessage_("Reflection angle is: " << thetaFrom);
    MyAssert_(mI != 0);
	
}

void 
StokesV::IncreasePhaseDueToPopagation(double inPhi)
{
	mTotalPropagationPhase += inPhi;
}

double
StokesV::Orientation()
{
    CheckDoubleValue_(mQ);
    CheckDoubleValue_(mU);
    
    return atan2(mU,mQ);
}

double
StokesV::degreeOfLinearPolarization()
{
    CheckDoubleValue_(mI);
    CheckDoubleValue_(mQ);
    CheckDoubleValue_(mU);
    
    return mI != 0 ? sqrt(mQ*mQ + mU*mU) / mI : 0.;
}

double
StokesV::degreeOfCircularPolarization()
{
    CheckDoubleValue_(mI);
    CheckDoubleValue_(mV);
    
    return mI != 0 ? abs(mV) / mI : 0.;
}

double
StokesV::degreeOfPolarization()
{
    CheckDoubleValue_(mI);
    CheckDoubleValue_(mQ);
    CheckDoubleValue_(mU);
    CheckDoubleValue_(mV);
    
    return mI != 0 ? sqrt(mQ*mQ + mU*mU + mV*mV) / mI : 0.;
}


std::ostream& operator<<(std::ostream& os, const StokesV& s) 
{
    return os << "(" << s.mI <<","<< s.mQ<<"," << s.mU<<"," << s.mV << ")";
}

