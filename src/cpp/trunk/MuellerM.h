#ifndef __MUELLERM__H
#define __MUELLERM__H

/*!
@class MuellerM
The class MuellerM is an abstract class used to obtain the matrix elements of a Mueller matrix.
The goal of this class is to be efficient enough for Monte Carlo calculations.
Hence, the descendants of this class can implement the various calculations using
lookup tables.  The necessary accesor functions are declared in MuellerM.

When this file is compiled with the __NOOPTIMIZE #defined, then "slower" functions are used.

*/
#include "fastinterpolate.h"
#include "constants.h"
#include "mydebug.h"

class MuellerM {
public:
	/*!
	@function ~MuellerM
     @discussion Destructor: because we have virtual functions, we must define a vritual destructor
     */
	
	virtual ~MuellerM();
	
    /*! 
    @function GetMatrixForTheta
    @discussion Function to obtain the Mueller matrix of a Rayleigh scatterer for a given scattering angle theta.
    
    @param inM This parameter contains the Mueller matrix on exit.  The mueller matrix is returned through this parameter 
	(which is a 
	double[16] array (if the code is compiled without __NOOPTIMIZE, the default) or double[4][4] array 
	(if the code is compiled with __NOOPTIMIZE). The caller must have allocated the array. A few things to know:
	matrices: m[row][col], hence using a single index one must do: [4*row+col] == [row][col]
	@param inTheta Scattering angle requested (in the scattering plane).
*/
#ifndef __NOOPTIMIZE
    virtual void GetMatrixForTheta(double *inM, double inTheta) = 0;
#else
    virtual void GetMatrixForTheta(double inM[4][4], double inTheta) = 0;
#endif
    virtual void GetMatrixTopTwoElementsForTheta(double& outM11, double& outM12, double inTheta) = 0 ;
    
    /*! 
        @function Normalization
        @discussion Returns the normalization factor corresponding to an integral of the intensity distribution
        over 2 pi.  The normalization should be calculated in the constructor of the matrix.
        @result normalization factor is : 2 pi Integrate[Intensity[theta,0],{x,0, 2pi}]
        */
    virtual double Normalization();
    /*! 
        @function MaxThetaPhi
        @discussion Returns the maximum value of Intensity[theta,0].  The maximum value should be calculated
        in the constructor of the matrix.
        @result Maximum factor is : Max[Intensity[theta,0]]
        */
    virtual double MaxThetaPhi() ;
    
    /*! 
        @function GetAnisotropyCoefficient
        @discussion Returns the anisotropy coefficient as calculated for unpolarized incident light.
        @result 
        */
    
    virtual double GetAnisotropyCoefficient() ;    
	
	virtual void linear_polarizer(double *inM, int axis);


protected:
        /*! @var mNormalization            Contains normalization factor */
        double mNormalization;
    /*! @var mMaxThetaPhi            Contains maximum intensity */
    double mMaxThetaPhi;
    /*! @var mG            Holds the value of "g" (i.e. avg(cosine)) for fast look-up */
    
    double mG;
    
};

/*!
@class MuellerMRayleigh
Provides the necessary encapsulation for a Rayleigh scatterer.

*/


class MuellerMRayleigh : public MuellerM {
    
	/*!
	@function ~MuellerMRayleigh
     @discussion Destructor: because we have virtual functions, we must define a vritual destructor
     */
	
	virtual ~MuellerMRayleigh();
	
    /*
     @function Normalization
     @discussion Returns the normalization factor corresponding to an integral of the intensity distribution
     over 2 pi.  The normalization is not implemented yet for Rayleigh matrices.  
     Avoid problems by just dying here with a throw.
     
     */
    
    virtual double Normalization() ;
    
    /*! 
    @function GetMatrixForTheta
    @discussion Function to obtain the Mueller matrix of a Rayleigh scatterer for a given scattering angle theta.  
    The matrix elements are calculated analytically.
    
    @param inM This parameter contains the Mueller matrix on exit.  The mueller matrix is returned through this parameter (which is a 
               double[16] array (if the code is compiled without __NOOPTIMIZE, the default) or double[4][4] array 
               (if the code is compiled with __NOOPTIMIZE). The caller must have allocated the array. A few things to know:
               matrices: m[row][col], hence using a single index one must do: [4*row+col] == [row][col]
    @param inTheta Scattering angle requested (in the scattering plane).
    @param inPhi This parameter *would* allow a rotation of the preferred axis (not implemented now). Always zero.
               */
    
#ifndef __NOOPTIMIZE
    virtual void GetMatrixForTheta(double *inM, double inTheta) ;
#else
    virtual void GetMatrixForTheta(double inM[4][4], double inTheta);
#endif
    
		
    
};

/*!
@class MuellerMSpherical
Provides the necessary encapsulation for a spherical scatterer and calculates the scattering elements
from four interpolating tables (provided by the user).

*/


class MuellerMSpherical : public MuellerM {
public:
    /*
     @function MuellerMSpherical
     @discussion This constructor allows a quick calculation of the Mueller matrix elements of a spherical scatterer
     by using lookup tables.  The matrix is of the form:
     First row of the matrix  : (S11 , S12,  0,  0)
     Second row of the matrix : (S12 , S11,  0, 0)
     Third row of the matrix  : (0 ,  0,  S33, S34)
     Fourth row of the matrix : (0 ,  0, -S34, S33)
     
     @param S11filename File containing two columns with angles in degree [0,pi] (first) and S11 (second)
     @param S12filename File containing two columns with angles in degree [0,pi] (first) and S12 (second)
     @param S33filename File containing two columns with angles in degree [0,pi] (first) and S33 (second)
     @param S34filename File containing two columns with angles in degree [0,pi] (first) and S34 (second)
     */
    
    MuellerMSpherical(string S11filename, string S12filename,string S33filename,string S34filename);    
    /*
     @function MuellerMSpherical
     @discussion This constructor allows a quick calculation of the Mueller matrix elements of a spherical scatterer
     by using lookup tables.  The matrix is of the form:
     First row of the matrix  : (S11 , S12,  0,  0)
     Second row of the matrix : (S12 , S11,  0, 0)
     Third row of the matrix  : (0 ,  0,  S33, S34)
     Fourth row of the matrix : (0 ,  0, -S34, S33)
     
     @param inX contains a pointer to an array of x values
     @param inS11 contains a pointer to an array of S11 values
     @param inS12 contains a pointer to an array of S12 values
     @param inS33 contains a pointer to an array of S33 values
     @param inS34 contains a pointer to an array of S34 values
     @param inN contains the number of points in array
     
     */
    
    MuellerMSpherical(double* inX, double* inS11, double* inS12,double* inS33,double* inS34, long inN);

    /*! 
        @function ~MuellerMSpherical
        @discussion Empty virtual destructor.
        
        */
    
    virtual ~MuellerMSpherical();

    /*! 
    @function init
    @discussion Function to pre-calculate various parameters for the matrix.
    
    */
    
    
    void init() ;
    /*! 
    @function GetMatrixForTheta
    @discussion Function to obtain the Mueller matrix of a spherical scatterer for a given scattering angle theta.  
    The matrix elements are obtained from an interpolating table provided when the class is instantiated.
    
    @param inM This parameter contains the Mueller matrix on exit.  The mueller matrix is returned through this parameter (which is a 
                                                                                                                           double[16] array (if the code is compiled without __NOOPTIMIZE, the default) or double[4][4] array 
                                                                                                                           (if the code is compiled with __NOOPTIMIZE). The caller must have allocated the array. A few things to know:
                                                                                                                           matrices: m[row][col], hence using a single index one must do: [4*row+col] == [row][col]
                                                                                                                           @param inTheta Scattering angle requested (in the scattering plane).
                                                                                                                           */
    
#ifndef __NOOPTIMIZE
    virtual void GetMatrixForTheta(double *inM, double inTheta);
#else
    virtual void GetMatrixForTheta(double inM[4][4], double inTheta);
#endif
    
    virtual void GetMatrixTopTwoElementsForTheta(double& outM11, double& outM12, double inTheta);
    
    virtual void GetCopyOfInterpolationTables(fastinterpolate& outS11, fastinterpolate& outS12, fastinterpolate& outS33, fastinterpolate& outS34) ;
    

    
protected:
    /*! @var S11	   	        Contains interpolating table for S11 element */
    fastinterpolate S11;
    /*! @var S12    		Contains interpolating table for S12 element */
    fastinterpolate S12;
    /*! @var S33            	Contains interpolating table for S33 element */
    fastinterpolate S33;
    /*! @var S34            	Contains interpolating table for S34 element */
    fastinterpolate S34;
    
};


#endif
