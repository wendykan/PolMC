
#ifndef __STOKESV__H
#define __STOKESV__H

#include "mydebug.h"
#include "constants.h"
#include <stdexcept>
#include "RealV.h"

#define CheckStokesV_(s) MyAssert_(s.degreeOfPolarization() <= 1.001)

using namespace std;

class StokesV;



/*!
@class StokesV
The class StokesV provides the necessary encapsulation and manipulation for 
keeping track of the four Stokes parameters for a photon as well as some
functions for manipulating reference frames.  The reference frames must be 
handled separately.

For a reference on the definition of the parameters, see Bohren and Huffman p. 52, eq. 2.84.   We use:

The I parameter (first Stokes parameter) is E_parallel*conjugate(E_parallel) + E_perp*conjugate(E_perp)
The Q parameter (second Stokes parameter) is E_parallel*conjugate(E_parallel) - E_perp*conjugate(E_perp)
The U parameter (third Stokes parameter) is E_perp*conjugate(E_parallel) + E_parallel*conjugate(E_perp)
The V parameter (fourth Stokes parameter) is sqrt(-1) * (E_perp*conjugate(E_parallel) - E_parallel*conjugate(E_perp))

The triad E_perpendicular, E_parallel, E_propagation is right handed 
(i.e. (E_perpendicular cross E_parallel ) dot E_propagation == 1). Since most books take E_parallel
to be horizontal, we need to have either E_parallel = (-1,0,0) or E_perpendicular (0,-1,0).
I chose E_parallel = (-1,0,0) to match figure 2.14 and 2.15 in Bohren and Huffman,

*/

class StokesV {
    
public:
    
    /*!
    @function StokesV
     @discussion Default constructor (no arguments) is made such that the intensity is 0 so no confusion possible.
     */
    StokesV() ;
    /*! 
    @function StokesV
    @discussion Constructor that uses four arguments (the four Stokes parameters)
    which are all real numbers.  
    
    @param inI The intensity (first Stokes parameter) 
    @param inQ The Q parameter (second Stokes parameter) 
    @param inU The U parameter (third Stokes parameter)
    @param inV The V parameter (fourth Stokes parameter)        
    */
    StokesV(double inI, double inQ, double inU, double inV, 
            RealV inEpara , RealV inEprop);
    /*! 
        @function StokesV
        @discussion Constructor that uses four arguments (the four Stokes parameters)
        which are all real numbers.  
        
        @param inI The intensity (first Stokes parameter) 
        @param inEprop Direction of propagation
        */
    StokesV(double inI, RealV inEprop);
    /*! 
        @function StokesV
        @discussion Constructor that uses two arguments (the two electric fields)
        which are both complex numbers.
        
        @param inEl The parallel electric field E_parallel
        @param inEr The perpendicular electric field E_perp
        */
    
    StokesV(Complex inEl, Complex inEr);
    /*!
        @function StokesV
     @discussion Constructor that uses another StokesV as an argument
     
     @param inS The Stokes parameter
     */
    StokesV(const StokesV& inS);

	/*!
        @function ~StokesV
     @discussion Destructor: because we have virtual functions, we must define a vritual destructor
     */
    virtual ~StokesV();
	
	static void SetCoherentSum(bool inCoherent);
	static bool GetCoherentSum();
	
    /*! 
        @function Normalize
        @discussion This function normalizes the intensity of the Stokes vector to 1.
        (i.e. it divides all the Stokes parameters by the current intensity)
        If the intensity is zero, it throws an error.
        
        */
    void   Normalize() ;        
    /*! 
        @function operator=
        @discussion This function makes a copy of a StokesV object.  This function is called
        when a StokesV is passed to or returned from a function.   
        @param inRhs This is an implicit parameter corresponding to the right hand side of the
        assignment.
        */
    
    StokesV& operator=(const StokesV& inRhs) ;
    /*! 
        @function operator+=
        @discussion This function adds a StokesV to another (cumulative) StokesV.  The purpose of
        this function is to allow easy statistics to be compiled when running Monte Carlo calculations.
        Statistics can be obtained simply by adding the current StokesV to a cumulative StokesV.
        
        @param inRhs This is an implicit parameter corresponding to the right hand side of the
        assignment.
        */
    StokesV& operator+=(StokesV inRhs) ;
    /*! 
        @function operator*=
        @discussion This function multiplies a StokesV by a real constant. 
        
        @param inRhs This is an implicit parameter corresponding to the right hand side of the
        assignment and represents the constant by which the StokesV is multiplied.
        */
    
    StokesV& operator*=(double inRhs);
    
    /*! 
        @function operator/=
        @discussion This function divides a StokesV by a real constant. 
        
        @param inRhs This is an implicit parameter corresponding to the right hand side of the
        assignment and represents the constant by which the StokesV is divided.
        */
    
    StokesV& operator/=(double& inRhs);
    /*!
        @function operator==
     @discussion This function compares two StokesV objects and returns true is they are identical.
     
     @param inRhs This is an implicit parameter corresponding to the right hand side of the
     comparison.
     */
    
    bool operator==(StokesV inRhs);
    /*!
        @function operator!=
     @discussion This function compares two StokesV objects and returns true is they are different.
     
     @param inRhs This is an implicit parameter corresponding to the right hand side of the
     comparison.
     */
    
    bool operator!=(StokesV inRhs);    
    
    /*!
        @function GetReferenceFrame
     @discussion This function returns reference frame mEpara and mEperp.
     
     @param outEperp The value of E_perp is returned through this parameter
     @param outEpara The value of E_parallel is returned through this parameter
     @param outEperp The value of E_prop is returned through this parameter
     
     */
    void GetReferenceFrame(RealV& outEperp, RealV& outEpara, RealV& outEprop) ;
    /*!
        @function SetReferenceFrame
     @discussion This function sets the reference frame mEpara and mEperp.
     
     @param inEperp The new value of E_perp
     @param inEpara The new value of E_parallel
     @param inEperp The new value of E_prop
     
     */
    void SetReferenceFrame( RealV& inEperp, RealV& inEpara, RealV& inEprop) ;    
    /*!
        @function GetEparaAxis
     @discussion This function returns reference axis mEpara.
     
     @result The value of E_parallel is returned
     
     */
    RealV GetEparaAxis() ;
    /*!
        @function GetEperpAxis
     @discussion This function returns reference axis mEperp.
     
     @result The value of E_perp is returned
     
     */
    RealV GetEperpAxis() ;
    /*!
        @function GetPropagationDirection
     @discussion This function returns reference axis mEprop.
     
     @result The direction of propagation is returned
     
     */
    RealV GetPropagationDirection() ;    
	/*!
        @function SetPropagationDirectionInLabFrame
	 @discussion Sets the propagation direction of the photon.  
	 mEperp and mEpara are valid, but arbitrary after the function call.  The Stokes vector
	 coefficient are not adjusted.
	 */
    virtual void SetPropagationDirectionInLabFrame(RealV inDirection);
    
    /*!
        @function RotateReferenceFrameAroundPropagationDirectionBy
     @discussion This function rotates the reference frame mEpara and mEperp by an angle
     phi around mEprop.  This does not change the polarization state of the StokesV: it only
     changes the reference frame.
     
     @param inPhi Angle by which the frame of reference is rotated.  Can be positivie or negative.
     Positive is counterclockwise.
     
     The transformation follows this reasoning:
     
     In a frame of reference {er, el, ez}, those vectors are (1,0,0), (0,1,0) and (0,0,1) respectively.
     A rotation around ez is therefore simple and is given by the following matrix:
     
                /  Cos phi -Sin phi 0 \
     Rz(phi)  = |  Sin phi Cos phi 0 |
                \     0       0    1 /
     
     The matrix for changing basis from {e} to {xyz} is
     made of {el,er,ez} in column, that is:
     
                             / er.x el.x ez.x \
     ChgBasisFrom_e_to_xyz = | er.y el.y ez.y |
                             \ er.z el.z ez.z /
     
     this gives the transformed vector in xyz coordinates. For instance, 
     
     er = ChgBasisFrom_e_to_xyz . Rz(phi) . (1,0,0)
     el = ChgBasisFrom_e_to_xyz . Rz(phi) . (0,1,0)
     ez = unchanged
     
     */
    
    void RotateReferenceFrameAroundPropagationDirectionBy(double inPhi);
    /*! 
        @function RotatePolarizationStateBy
        @discussion Rotates the Stokes vector by an angle phi.  This function calls another
        overloaded function that does the actual work.
        
        @param inPhi Angle by which the polarization vector is rotated.  Can be positivie or negative.
        
        */
    
    void RotatePolarizationStateBy(double inPhi) ;
    
    /*! 
        @function RotatePolarizationStateBy
        @discussion Rotates the polarization state by an angle phi around the direction of propagation.
        This applies the rotation matrix to the Stokes vector:
        
        First row of the matrix  : (1 ,  0,          0,           0)
        Second row of the matrix : (0 ,  inCos_2phi, -inSin_2phi , 0)
        Third row of the matrix  : (0 ,  inSin_2phi, inCos_2phi , 0)
        Fourth row of the matrix : (0 ,  0,          0,           1)
        
        @param inCos_2phi Cosine of the angle by which the Stokes vector is rotated.  Can be positive or negative.
        @param inSin_2phi Sine of the angle by which the Stokes vector is rotated.  Can be positive or negative.
        
        */
    void RotatePolarizationStateBy(double inCos_2phi, double inSin_2phi);
    /*!
        @function ChangePropagationDirectionAroundEPerpBy
     @discussion This function changes the direction of propagation while keeping the polarization
     unchanged in the local frame.
     
     @param inTheta Angle by which the Stokes vector is rotated.  Can be positivie or negative.
     Positive is counterclockwise.
     
     The transformation follows this reasoning:
     
     In a frame of reference {er,el,ez}, those vectors are (1,0,0), (0,1,0) and (0,0,1) respectively.
     A rotation around er is therefore simple and is given by the following matrix:
     
                /  1        0            0    \
     Rx(phi)  = |   0    Cos theta  -Sin theta  |
                \   0    Sin theta   Cos theta /
     
     The matrix for changing basis from {e} to {xyz} is
     made of {el,er,ez} in column, that is:
     
                             / er.x el.x ez.x \
     ChgBasisFrom_e_to_xyz = | er.y el.y ez.y |
                             \ er.z el.z ez.z /
     
     this gives the transformed vector in xyz coordinates. For instance,
     
     er = unchanged
     el = ChgBasisFrom_e_to_xyz . Rx(theta) . (0,1,0)
     ez = ChgBasisFrom_e_to_xyz . Rx(theta) . (0,0,1)
     
     */
    
    void ChangePropagationDirectionAroundEPerpBy(double inTheta);
    
	/*!
		@function TransmitThrough
		@discussion
	 */
	
	virtual void TransmitThrough(RealV& inNormal, double rp, double rs, double tp, double ts) ;

	/*!
		@function ReflectAtInterface
		@discussion
	 */
	
	virtual void ReflectAtInterface(RealV& inNormal, double rp, double rs, double tp, double ts);
	
    /*!
	 @function GetLocalComplexFields
     @discussion From the Stokes parameters, this function returns the complex fields
     E_parallel and E_perp. Don't forget that absolute phase is meaningless: one can multiply both electric
     fields by -1 without changing the physical meaning.
     
     @param outEl The value of E_parallel is returned through this parameter
     @param outEr The value of E_perp is returned through this parameter
     
     */
	
    void GetLocalComplexFields(Complex& outEl, Complex& outEr) ;

	/*!
	@function IncreasePhaseDueToPopagation
	 @discussion Adds a constant phase to both fields.
	 
	 @param inPhi Can be positive or negative.
	 
	 */
	virtual void IncreasePhaseDueToPopagation(double inPhi);
	
	/*!
        @function Orientation
     @discussion Calculates the orientation of the polarized component.
     @result A double representing the angle (in radians) of polarization with the E parallel direction.
     */

    
    double Orientation();
    /*!
        @function degreeOfLinearPolarization
     @discussion Calculates the degree of linear polarization of the Stokes vector.
     @result A double containing the degree of linear polarization. If intensity is zero, returns 0.
     */
    
    double degreeOfLinearPolarization() ;
    /*! 
        @function degreeOfCircularPolarization
        @discussion Calculates the degree of circular polarization of the Stokes vector.
        @result A double containing the degree of circular polarization. If intensity is zero, returns 0.
        */
    double degreeOfCircularPolarization() ;
    /*! 
        @function degreeOfPolarization
        @discussion Calculates the degree of total polarization of the Stokes vector.
        @result A double containing the degree of total polarization. If intensity is zero, returns 0.
        */
    double degreeOfPolarization() ;
    
public:
    /*! @var mI            First Stokes parameter E_parallel*conjugate(E_parallel) + E_perp*conjugate(E_perp) */
    double mI;
    /*! @var mQ            Second Stokes parameter E_parallel*conjugate(E_parallel) - E_perp*conjugate(E_perp) */
    double mQ;
    /*! @var mU            Third Stokes parameter E_perp*conjugate(E_parallel) + E_parallel*conjugate(E_perp)*/
    double mU;
    /*! @var mV            Fourth Stokes parameter sqrt(-1) * (E_perp*conjugate(E_parallel) - E_parallel*conjugate(E_perp))*/
    double mV;
    
    /*! @var mEpara        3D vector representing the E_parallel direction.   E_perp cross E_parallel = EProp */
    RealV mEpara;
    /*! @var mEperp        3D vector representing the E_perp direction.  E_perp cross E_parallel = EProp  */
    RealV mEperp;
    /*! @var mEprop        3D vector representing the propagation direction.  E_perp cross E_parallel = EProp. */
    RealV mEprop;

	/* @var mTotalPropagationPhase Total propagation phase */
	double mTotalPropagationPhase;
	
	/* @var coherentSum		Boolean value that determines if Stokes vectors are added coherently or incoherently */
	static bool coherentSum;
};

/*!
@function operator<<
 @discussion Outputs the Stokes vector to a stream (standard out or a file for instance).
 
 @param os This is an implicit parameter corresponding to the stream
 @param s  This is an implicit parameter corresponding to the StokesV
 */

ostream& operator<<(ostream& os, const StokesV& s);


#endif 
