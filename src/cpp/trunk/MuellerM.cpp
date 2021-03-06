#include "MuellerM.h"


MuellerM::~MuellerM()
{
}

double
MuellerM::Normalization()
{
    return mNormalization; 
}

double
MuellerM::MaxThetaPhi()
{
    return mMaxThetaPhi; 
}

double
MuellerM::GetAnisotropyCoefficient()
{
    return mG;
}


//linear polarizer added by Wendy Kan
void
MuellerM::linear_polarizer(double *inM, int axis)
{                    
	if (axis == 1) { // polarize along x axis
    // To get [row][col] right, check Kernighan Ritchie, p. 113
    // [row][col] == 4*row+col
    // first row
		inM[0] = 0.5;
		inM[1] = -0.5;
		inM[2] = 0.;
		inM[3] = 0.;
		
		// second row
		inM[4] = -0.5 ;
		inM[5] = 0.5 ;
		inM[6] = 0. ;
		inM[7] = 0. ;
		
		// third row
		inM[8] = 0. ;
		inM[9] = 0. ;
		inM[10] =  0. ;
		inM[11] =  0. ;
		
		// fourth row
		inM[12] = 0.;
		inM[13] = 0.;
		inM[14] = 0.;
		inM[15] = 0.;
	}
	else if (axis == 2){ //polarize along y axis
		inM[0] = 0.5;
		inM[1] = 0.5;
		inM[2] = 0.;
		inM[3] = 0.;
		
		// second row
		inM[4] = 0.5 ;
		inM[5] = 0.5 ;
		inM[6] = 0. ;
		inM[7] = 0. ;
		
		// third row
		inM[8] = 0. ;
		inM[9] = 0. ;
		inM[10] =  0. ;
		inM[11] =  0. ;
		
		// fourth row
		inM[12] = 0.;
		inM[13] = 0.;
		inM[14] = 0.;
		inM[15] = 0.;
	
	
	}
	else{
		throw runtime_error("not valid, polarizer should not be passing z direction");

	}
	
}



MuellerMRayleigh::~MuellerMRayleigh()
{
}

double
MuellerMRayleigh::Normalization()
{
    throw; 
}

#ifndef __NOOPTIMIZE
void
MuellerMRayleigh::GetMatrixForTheta(double *inM, double inTheta)
{                    double a,b,d,e;
    
    /*
     This is Rayleigh scattering
     */
    
    a = 3./16./PI * (1. + cos(inTheta)*cos(inTheta));
    b = 3./16./PI * (-1. + cos(inTheta)*cos(inTheta));
    d = 3./8./PI * cos(inTheta) ;
    e = 0.;
    
    // To get [row][col] right, check Kernighan Ritchie, p. 113
    // [row][col] == 4*row+col
    // first row
    inM[0] = a;
    inM[1] = b;
    inM[2] = 0.;
    inM[3] = 0.;
    
    // second row
    inM[4] = b ;
    inM[5] = a ;
    inM[6] = d ;
    inM[7] = e ;
    
    // third row
    inM[8] = - b ;
    inM[9] = - a ;
    inM[10] =  d ;
    inM[11] =  e ;
    
    // fourth row
    inM[12] = 0.;
    inM[13] = 0.;
    inM[14] = -e;
    inM[15] = d;
}
#else
void
MuellerMRayleigh::GetMatrixForTheta(double inM[4][4], double inTheta)
{
    double a,b,d,e;
    
    /*  
    This is Rayleigh scattering
    */
    
    a = 3./16./PI * (1. + cos(inTheta)*cos(inTheta));
    b = 3./16./PI * (-1. + cos(inTheta)*cos(inTheta));
    d = 3./8./PI * cos(inTheta) ;
    e = 0.;
    
    double cos_2phi = cos(2. * inPhi);
    double sin_2phi = sin(2. * inPhi);
    
    // first row 
    
    inM[0][0] = a;
    inM[0][1] = b;
    inM[0][2] = 0.;
    inM[0][3] = 0.;
    
    // second row 
    inM[1][0] = b * cos_2phi;
    inM[1][1] = a * cos_2phi;
    inM[1][2] = d * sin_2phi;
    inM[1][3] = e * sin_2phi;
    
    // third row 
    inM[2][0] = - b * sin_2phi;
    inM[2][1] = - a * sin_2phi;
    inM[2][2] =   d * cos_2phi;
    inM[2][3] =   e * cos_2phi;
    
    // fourth row 
    inM[3][0] = 0.;
    inM[3][1] = 0.;
    inM[3][2] = -e;
    inM[3][3] = d;
}
#endif

MuellerMSpherical::MuellerMSpherical(string S11filename, string S12filename,string S33filename,string S34filename)
:S11(S11filename),S12(S12filename),S33(S33filename),S34(S34filename)
{
    init();
}

MuellerMSpherical::MuellerMSpherical(double* inX, double* inS11, double* inS12,double* inS33,double* inS34, long inN)
:S11(inX, inS11, inN),S12(inX, inS12, inN),S33(inX, inS33, inN),S34(inX, inS34, inN)
{
    init();
}

MuellerMSpherical::~MuellerMSpherical()
{
    init();
}

void
MuellerMSpherical::init()
{
    // Normalization function does not use I(theta,phi), just I(theta).
    mNormalization = 0.;
    double dx = PI/1000.;
    for (double x = 0; x <= PI; x += dx) {
        mNormalization += S11.y(x) * sin(x);
    }
    mNormalization *= dx * 2. * PI;
    mMaxThetaPhi = S11.y(0.);
    
    clog <<  RightNow() << "WARNING: Normalization factor for Mueller matrix uses S11(theta) only" << endl;
    clog <<  RightNow() << "Normalization factor for Mueller matrix = " << mNormalization <<  endl;
    clog <<  RightNow() << "Maximum for I(0,0) = " << mMaxThetaPhi <<  endl;
    
    dx = PI/10000.;
    double theSum=0;
    double theNorm=0;
    for (double x = 0; x <= PI; x += dx) {
        theNorm += S11.y(x) * sin(x);
    }
    
    theSum = 0;
    for (double x = 0; x <= PI; x += dx) {
        theSum += S11.y(x) * cos(x) * sin(x);
    }
    mG = theSum/theNorm;
    clog <<  RightNow() << "Calculated g from matrix = " << mG <<  endl;
    
}
#ifndef __NOOPTIMIZE
void
MuellerMSpherical::GetMatrixForTheta(double *inM, double inTheta)
{        double a,b,d,e;
    
    /*
     This is Mie scattering
     Obtains data from a table and interpolate
     */
    
    a = S11.y(inTheta);
    b = S12.y(inTheta);
    d = S33.y(inTheta);
    e = S34.y(inTheta);
    
    // To get [row][col] right, check Kernighan Ritchie, p. 113
    // [row][col] == 4*row+col
    // first row
    inM[0] = a;
    inM[1] = b;
    inM[2] = 0.;
    inM[3] = 0.;
    
    // second row
    inM[4] = b ;
    inM[5] = a ;
    inM[6] = 0 ;
    inM[7] = 0 ;
    
    // third row
    inM[8] = 0 ;
    inM[9] = 0 ;
    inM[10] =  d ;
    inM[11] =  e ;
    
    // fourth row
    inM[12] = 0.;
    inM[13] = 0.;
    inM[14] = -e;
    inM[15] = d;
}
#else
void
MuellerMSpherical::GetMatrixForTheta(double inM[4][4], double inTheta)
{
    double a,b,d,e;
    
    /*
     This is Mie scattering
     Obtains data from a table and interpolate
     */
    
    a = S11.y(inTheta);
    b = S12.y(inTheta);
    d = S33.y(inTheta);
    e = S34.y(inTheta);
    
    // first row
    
    inM[0][0] = a;
    inM[0][1] = b;
    inM[0][2] = 0.;
    inM[0][3] = 0.;
    
    // second row
    inM[1][0] = b ;
    inM[1][1] = a ;
    inM[1][2] = 0 ;
    inM[1][3] = 0 ;
    
    // third row
    inM[2][0] = 0 ;
    inM[2][1] = 0 ;
    inM[2][2] = d ;
    inM[2][3] = e ;
    
    // fourth row
    inM[3][0] = 0.;
    inM[3][1] = 0.;
    inM[3][2] = -e;
    inM[3][3] = d;
}
#endif

void
MuellerMSpherical::GetMatrixTopTwoElementsForTheta(double& outM11, double& outM12, double inTheta)
{ 
    /*
     This is Mie scattering
     Obtains data from a table and interpolate
     */
    
    outM11 = S11.y(inTheta);
    outM12 = S12.y(inTheta);
    
}

void
MuellerMSpherical::GetCopyOfInterpolationTables(fastinterpolate& outS11, fastinterpolate& outS12, fastinterpolate& outS33, fastinterpolate& outS34)
{        outS11 = S11;
    outS12 = S12;
    outS33 = S33;
    outS34 = S34;
}

