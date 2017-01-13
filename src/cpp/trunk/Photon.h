/*! @header Photon.h
    This file contains the core of the photon propagation code.
    Various formalisms have been presented in the literature.
    The original references upon which this code is based are:

"Monte Carlo simulations of the diffuse backscattering 
 Mueller matrix for highly scattering media", Bartel & Hielscher
 Appl. Opt. Vol 39, April 1st 2000, p. 1580
 
"Characterizing mammalian cells and cell phantoms by 
 polarized backscatterng fiber-optic measurements"
 J. Mourant, T. Johnson, J. Freyer
 Appl. Opt. Vol 40, No8 Oct 2001, p. 5114

"Description and time reduction of a Monte Carlo code ...",
F. Jaillon, H. Saint Jalmes
Appl. Opt. Vol 42, no 16, p.3290

 which are themselves loosely based on code from Wang and Jacques

"MCML-Monte Carlo modeling of light transport" in Optical-Thermal
 response of laser irradiated tissue, A. J. Welch M. van Gemert 
 pp. 73-100 http://omlc.ogi.edu/software/mc/index.html


This C/C++ file contains both formalisms and it is possible to use one
or the other at run time. This allows you to compare the results from
either formalism. This is done by defining  distinct Photon classes
for Bartel, Cote, Jaillon and Mourant.  Although the variables (i.e. reference frames)
are conceptually the same, the notation used in both papers is
different.

The goal of this code is to be correct and expandable, not necessarily fast.
In debug mode, there are tons of consistency checks that are done with various
functions that end with a "_". These functions do their check if
__MYDEBUG is defined and do absolutely nothing if __MYDEBUG is not
defined.  Hence: it does not hurt to have them everywhere since they will not do
anything in  the final version.
The idea is therefore to check that the code is right by
running a few test simulations with __MYDEBUG defined, then do "the real
stuff" without __MYDEBUG.
The ./configure script provides a --enable-debug option which does just that.

If one wanted to make a class called  PhotonFast, the variable functions
could be rewritten and optimized.
*/

#ifndef __PHOTON_H
#define __PHOTON_H

#include <list>
#include "mydebug.h"
#include "StokesV.h"
#include "RealV.h"
#include "MuellerM.h"
#include "MCObject.h"

enum { kBartel=1, kMourant, kCote, kEvans, kKaplan, kHenyeyGreenstein, kJaillon, kIntensity};

enum { kRaw, kMathematica, kBinaryFormat };

class MCObject;


/*!
    @class Photon

    The class photon provides the necessary encapsulation for manipulating photon wavepackets
    with their polarization.  It keeps track of direction, distance traveled, position
    and weight, where the weight is initially 1 and is reduced due to partial absorption.
    The class Photon is an abstract class.  Other classes (PhotonBartel, PhotonCote, and PhotonMourant)
    implement some of the abstract functions following published implementations.
*/

class Photon {
    public:
         /*! 
        @function Photon
        @discussion Default constructor
		*/
    Photon() ;

         /*! 
        @function Photon
        @discussion Constructor of a photon from the four Stokes parameters.
        Position is set to (0,0,0), weight to 1.
                    
        @param inI The intensity (first Stokes parameter) 
        @param inEprop The direction of propagation
	*/
    Photon(double inI, RealV inEprop) ;
         /*! 
        @function Photon
        @discussion Constructor of a photon from the four Stokes parameters.
        Position is set to (0,0,0), weight to 1.
                    
        @param inI The intensity (first Stokes parameter) 
        @param inQ The Q parameter (second Stokes parameter) 
        @param inU The U parameter (third Stokes parameter)
        @param inV The V parameter (fourth Stokes parameter)        
        @param inEpara A vector in the parallel plane (can be Epara, but doesn't have to be):
					If Eprop and Epara are not perpendicular, Epara will be adjusted.
        @param inEprop The direction of propagation
        */
    Photon(double inI, double inQ, double inU, double inV, RealV inEpara , RealV inEprop);
    /*! 
        @function Photon
        @discussion Constructor of a photon from the a Stokes vector.
        Position is set to (0,0,0), weight to 1.
                    
        @param inStokesVector The initial Stokes vector
        */

    Photon(StokesV& inStokesVector );
		/*
		@function ClearAll
		@discussion Clears all variable
		*/
    virtual void ClearAll();		
        /*!
            @function MakeCopy
         @discussion Make a copy of a photon from the another Photon.

         @param inPhoton The other Photon
         */

        virtual Photon* MakeCopy() = 0;        

        /*!
         @function ~Photon
         @discussion Destructor

         */

        virtual ~Photon();        
        /*!
            @function init
            @discussion Initialize anything that is left to initialize.
        */
        virtual void init();

         /*! 
        @function IsDead
        @discussion If the weight is zero, the photon is dead.

        @result True if photon is dead.
        */

        bool IsDead() ;
         /*! 
        @function IsNotDead
        @discussion If the weight is zero, the photon is dead.

        @result True if photon is not dead.
        */

        bool IsNotDead() ;
         /*! 
        @function GetWeight
        @discussion Returns the weight of a photon.  The weight should always be between 0 and 1.

        @result The weight of the photon.
        */

        virtual double GetWeight() ;
         /*! 
        @function SetWeight
        @discussion Set the weight of a photon

        @param inWeight The weight of the photon.
        */
        virtual void SetWeight(double inWeight);
        /*! 
        @function MultiplyWeightBy
        @discussion Multiply the weight of a photon by a constant

        @param inFactor The factor by which the photon is multiplied
        @result 0 if weight is still within legal range [0,1], 1 if not.
        */
        virtual bool MultiplyWeightBy(double inFactor);
        
        /*! 
        @function DecreaseWeightBy
        @discussion Decrease the weight of a photon by a constant

        @param inDeltaWeight The amount by which the photon is decreased
        @result 0 if weight is still within legal range [0,1], 1 if not.
        */
        virtual bool DecreaseWeightBy(double inDeltaWeight) ; 
        /*! 
        @function GetStokesVectorInLocalFrame
        @discussion Returns the Stokes vector in the local frame of reference of
        the photon.

        @param outS The Stokes vector is returned in this parameter.
        */
        virtual void GetStokesVectorInLocalFrame(StokesV& outS) ;
        /*!
            @function SetStokesVector
         @discussion Sets the Stokes vector in the local frame of reference of
         the photon.  Important: this also sets direction of propagation and reference
         frame for polarization.
         */
        virtual void SetStokesVector(StokesV& inS) ;
        
        /*!
            @function RotateReferenceFrameInFresnelPlane
         @discussion Rotates the reference plane of the Stokes vector in the Fresnel s and p plane
         @param inNormal
         */
        
        virtual void RotateReferenceFrameInFresnelPlane(RealV& inNormal)   ;
        /*!
            @function GetReflectionProbability
         @discussion Provides the Fresnel reflection and transmission coefficients
         @param inNormal
         @param inIndexFrom
         @param inIndexTO
         */

        virtual double GetReflectionProbability(RealV& inNormal, double rp, double rs, double tp, double ts);        
                                            
 
        /*!
            @function TransmitThroughInterface
         */

        virtual void TransmitThrough(RealV& inNormal, double inIndexFrom, double inIndexTo,
                                     double rp, double rs, double tp, double ts) ;
        /*!
            @function ReflectAtInterface
         */

        virtual void ReflectAtInterface(RealV& inNormal, double inIndexFrom, double inIndexTo,
                                        double rp, double rs, double tp, double ts);
        
        virtual void NormalizeStokesV();
        /*! 
        @function MultiplyStokesV
        @discussion Multiply the Stokes vector of a photon by a constant

        @param inFactor The factor by which the Stokes vector is multiplied
        */
        virtual  void MultiplyStokesV(double inFactor);
        /*!
            @function ChangePropagationDirectionAroundEPerpBy
         @discussion

         @param inTheta The angle by which the frame of reference is rotated
         */
        virtual void ChangePropagationDirectionAroundEPerpBy(double inTheta);        
        /*!
          @function This function rotates the reference frame mEpara and mEperp by an angle
          phi around mEprop.  This does not change the polarization state of the StokesV: it only
          changes the reference frame.
          
        @param inPhi The angle by which the frame of reference is rotated
        */
        virtual void RotateReferenceFrameAroundPropagationDirectionBy(double inPhi) ;
        
        /*!
         @function RotatePolarizationStateBy
         @discussion Rotates the polarization state by an angle phi.
          
         @param inPhi Angle by which the polarization state is rotated.  Can be positivie or negative.

         */
        virtual void RotatePolarizationStateBy(double inPhi);

        /*!
			@function IncreasePhaseDueToPopagation
         @discussion Adds a constant phase to both fields.
		 
         @param inPhi Can be positive or negative.
		 
         */
        virtual void IncreasePhaseDueToPopagation(double inPhi);
		
        double GetWavelength();
        
        void SetWavelength(double inWavelength);
        
        void SetCurrentObject(MCObject* inObject) ;
        
        MCObject* GetCurrentObject() ;
        
        /*! 
        @function GetLocalPosition
        @discussion Returns the position of the photon in local coordinates.
        @result The vector representing the position (x,y,z).
        */
        virtual RealV GetLocalPosition() ;


        /*!
            @function SetLocalPosition
         @discussion Sets the position of the photon in local coordinates.
         */
        virtual void SetLocalPosition(RealV inPos ) ;
        
        /*!
            @function SetGlobalPosition
         @discussion Sets the position of the photon in global coordinates.
         */
        virtual void SetGlobalPosition(RealV inPos ) ;
        
         /*! 
        @function GetGlobalOrigin
        @discussion Returns the position of the origin in global coordinates.
        @result The vector representing the position (x,y,z).
        */
        virtual RealV GetGlobalOrigin() ;

        /*!
            @function SetGlobalOrigin
         @discussion Sets the position of the origin in global coordinates, and adjust
         the mLocalPos to reflect it if demanded.
         */
        virtual void SetGlobalOrigin(RealV inNewOrigin, bool inTransformLocalCoordinates ) ;
        /*!
            @function GetGlobalPosition
         @discussion Returns the position of the photon in global coordinates.
         @result The vector representing the position (x,y,z).
         */
        virtual RealV GetGlobalPosition() ;

		/*! 
			@function GetTravelTime
			@discussion Returns the travel time of the photon.
			@result The time since last Distance resset.
			*/
        virtual double GetTravelTime() ;
		
        /*! 
        @function GetDistanceTraveled
        @discussion Returns the distance traveled by the photon.
        @result The distance traveled
        */
        virtual double GetDistanceTraveled() ;

        /*!
            @function ResetDistanceTraveled
         @discussion Resets distance traveled.
         */
        virtual void ResetDistanceTraveled() ;
        
        /*!
            @function GetNumberOfScatteringevents
         @discussion Returns the number of scattering events suffered by the photon.
         @result The number of scattering events
         */
        virtual long GetNumberOfScatteringEvents();        
        /*!
            @function GetPropagationDirectionInLabFrame
         @discussion Gives the propagation direction of the photon (mE3)
         @result The propagation direction mE3 of the photon
         */
        virtual RealV GetPropagationDirectionInLabFrame() ;

        /*!
        @function SetPropagationDirectionInLabFrame
         @discussion Sets the propagation direction of the photon (mE3).  
         ePerp and ePara are valid, but arbitrary after
         function call.
         */
        virtual void SetPropagationDirectionInLabFrame(RealV& inDirection);
         /*!
            @function IntensityThroughLinearPolarizer
         @discussion Gives the intensity after photon has propagated through linear polarizer
         @result The intensity
         */
        virtual void IntensityThroughLinearPolarizer(StokesV& inS, RealV& inEparaLab, RealV& inEperpLab, RealV& inNormalLab, double& outIpara, double& outIperp) ;
      
        /*!
		 @function IntensityThroughCircularPolarizer
         @discussion Gives the intensity after photon has propagated through a circular polarizer
         @result The intensity
         */
        virtual void IntensityThroughCircularPolarizer(StokesV& inS, RealV& inEparaLab, RealV& inEperpLab, RealV& inNormalLab, double& outIRightCirc, double& outILeftCirc);
  
        /*!
            @function MoveBy
         @discussion Displaces the photon by a specified distance, in the direction of the current
         propagation direction.
         @param      inDz distance by which the photon travels
		 @param		 index index of refraction, used to calculate optical path (i.e. time) and phase
         */

        virtual void MoveBy(double inDz, double index)  ;
        /*!
            @function ScatterBy
         @discussion This function scatters is the core of the class.  It computes the
         scattered photon at scattering angle  and azimuthal angle  given
         the Mueller matrix.   Upon return, the direction and reference frame
         of the photon will have changed.

         @param  inTheta Scattering angle [0:pi] (not in degrees, in radians) measured in the same frame of reference as the Mueller matrix
         @param  inPhi Azimuthal angle [0:2 pi] (not in degrees, in radians) measured in the same frame of reference as the Mueller matrix
         @param inM Mueller matrix of the scatterer, given in standard frame of reference (parallel and perpendicular)
         */

        virtual void  ScatterBy(double inTheta, double inPhi, MuellerM* inM);        

        /* These functions are defined below in PhotonBartel and
           PhotonMourant using the notation of the respective publication */
        
        virtual void MeasureStokesVectorInLabFrame(StokesV& outS, RealV& inEparaLab, RealV& inEperpLab, RealV& inNormalLab) = 0;
        
        virtual void ElectricFieldInLabFrame(StokesV& inS,
														   Complex& outFieldPara,
														   Complex& outFieldPerp , 
														   RealV& inEparaLab, 
														   RealV& inEperpLab, 
														   RealV& inNormalLab);
/*        virtual void ElectricFieldInLabFrameWithFullPhase(StokesV& inS, 
														  Complex& outFieldPara, 
														  Complex& outFieldPerp , 
														  RealV& inEparaLab, 
														  RealV& inEperpLab, 
														  RealV& inNormalLab);
*/		
		bool GetStats(int interaction, RealV& outPosition, double& outWeight);

		void StoreStats();
		
		void KeepStats();
        
        void DontKeepStats();

        virtual void DumpStats(ostream& s, long inFormat);
		
		RealV GetMaxDepth();

		
		
    protected:

        StokesV mS;
        double mWavelength;
        
        RealV mLocalPos;
        RealV mGlobalOrigin;
        MCObject* mCurrentlyInObject;
        
        double mDistanceTraveled;
        double mOpticalPathTraveled;
        double mWeight;

        list<StokesV> mStokesStats;
        list<RealV> mPosStats;
        list<RealV> mDirStats;
        list<double> mWeightStats;
        list<double> mTimeStats;
        bool mKeepStats;
        long mNumEvent;
        long mStatsMove;

        bool mWarningStatsGiven;

};

/*!
@class PhotonIntensity
The class PhotonIntensity implements the various actions performed on a photon wavepacket
without any consideration to polarization.

The goal of this implementation is to be fast, hence it overrides the polarization
function of Photon to short circuit them and make everything faster.  Since the Stokes vector
also contains the reference frame (scattering plane, direction of propagation), we still need
to compute the changes to those variables.


*/
class PhotonIntensity : public Photon {

public:
     /*!
     @function PhotonIntensity
     @discussion Default constructor
	*/

    PhotonIntensity( );
   /*!
    @function PhotonIntensity
     @discussion The Stokes vector is set to (intensity,0,0,0) which is used for Fresnel coefficients
     (and avoids having to reimplement everything).
     @param inI The intensity (first Stokes parameter)
     */

    PhotonIntensity( double inIntensity, RealV inProp);
    /*!
    @function init
     @discussion Initialize anything that is left to initialize, which for now means
     spitting out a warning to stdlog to make sure the user knows what he is using.
     */
    virtual void init() ;
    /*!
    @function MakeCopy
     @discussion Make a copy of a photon from the another Photon.

     @param inPhoton The other Photon
     */

    virtual Photon* MakeCopy();

    /*!
    @function MeasureStokesVectorInLabFrame
     @discussion The StokesV of this photon is always fully unpolarized.
     @param outS The Stokes vector is returned via this parameter (fully unpolarized photon)
     */

    virtual void MeasureStokesVectorInLabFrame(StokesV& outS, RealV& inEparaLab, RealV& inEperpLab, RealV& inNormalLab) ;

    /*!
    @function ScatterBy
     @discussion This function scatters is the core of the class.  It computes the
     scattered photon at scattering angle  and azimuthal angle.  Since we only consider intensity,
     we don't need to calculate the effect on polarization.  The Mueller matrix is not used
     (and should be set to zero).
     
     @param  inTheta Scattering angle [0:pi] (not in degrees, in radians) measured in the same frame of reference as the Mueller matrix
     @param  inPhi Azimuthal angle [0:2 pi] (not in degrees, in radians) measured in the same frame of reference as the Mueller matrix
     @param inM Mueller matrix of the scatterer, not used
     */

    virtual void  ScatterBy(double inTheta, double inPhi, MuellerM* inM) ;
    /*!
    @function RotatePolarizationStateBy
     @discussion Since we only consider intensity, we short circuit this calculation (for speed).

     @param inPhi Angle, not used.

     */
    virtual void RotatePolarizationStateBy(double inPhi);
    
};


/*!
    @class PhotonBartel
    The class PhotonBartel implements the various actions done on a photon wavepacket according to
    
    "Monte Carlo simulations of the diffuse backscattering 
    Mueller matrix for highly scattering media", Bartel & Hielscher
    Appl. Opt. Vol 39, April 1st 2000, p. 1580.
    
    The notation used in the implementation of this class is similar to the one used
    in the article to avoid confusion.  The goal of this implementation is to be readable, not necessarily
    fast.

*/
class PhotonBartel : public Photon {
 
 public:
    /*!
        @function PhotonBartel
        @discussion Default constructor
		
    */

    PhotonBartel() ;
    /*!
        @function PhotonBartel
        @discussion This constructor makes a photon with polarization state given by 
        the four Stokes parameters passed as arguments.  The photon is assumed to travel along the (0,0,1)
        direction and have the reference axes for the Stokes parameter to be e_parallel = (1,0,0) and
        e_perp = (0,1,0)
        @param inI The intensity (first Stokes parameter) 
        @param inQ The Q parameter (second Stokes parameter) 
        @param inU The U parameter (third Stokes parameter)
        @param inV The V parameter (fourth Stokes parameter)        
    */

    PhotonBartel(double inI, double inQ, double inU, double inV, RealV inEpara , RealV inEprop); 
    /*!
        @function PhotonBartel
        @discussion This constructor makes a photon with polarization state given by 
        the Stokes vector passed as argument.  The photon is assumed to travel along the (0,0,1)
        direction and have the reference axes for the Stokes parameter to be e_parallel = (1,0,0) and
        e_perp = (0,1,0)
        @param inStokesVector The Stokes vector 
    */

     
    PhotonBartel(StokesV& inStokesVector);
    /*!
        @function init
        @discussion Initialize anything that is left to initialize, which for now means
        spitting out a warning to stdlog to make sure the user knows what he is using.
    */
    virtual void init() ;
    /*!
    @function MakeCopy
     @discussion Make a copy of a photon from the another Photon.

     @param inPhoton The other Photon
     */

    virtual Photon* MakeCopy();    

    /*!
        @function MeasureStokesVectorInLabFrame
        @discussion Gives the Stokes vector with respect to the lab frame, by calculating
        it from the local frame of reference of the photon.
        @param      outS    The Stokes vector is returned via this parameter
    */

    virtual void MeasureStokesVectorInLabFrame(StokesV& outS, RealV& inEparaLab, RealV& inEperpLab, RealV& inNormalLab);


};


/*!
    @class PhotonCote
    The class PhotonCote implements the various actions done on a photon wavepacket according to
    
    "Monte Carlo simulations of the diffuse backscattering 
    Mueller matrix for highly scattering media", Bartel & Hielscher
    Appl. Opt. Vol 39, April 1st 2000, p. 1580.
    
    with a few modifications in the way the Photon is detected.
*/
class PhotonCote : public Photon {
 
 public:
    /*!
        @function PhotonCote
        @discussion Default constructor
		
    */

    PhotonCote()  ;
    
    /*!
        @function PhotonCote
        @discussion This constructor makes a photon with polarization state given by 
        the four Stokes parameters passed as arguments.  The photon is assumed to travel along the (0,0,1)
        direction and have the reference axes for the Stokes parameter to be e_parallel = (1,0,0) and
        e_perp = (0,1,0)
        @param inI The intensity (first Stokes parameter) 
        @param inQ The Q parameter (second Stokes parameter) 
        @param inU The U parameter (third Stokes parameter)
        @param inV The V parameter (fourth Stokes parameter)        
    */

    PhotonCote(double inI, double inQ, double inU, double inV, RealV inEpara , RealV inEprop); 
    /*!
        @function PhotonCote
        @discussion This constructor makes a photon with polarization state given by 
        the Stokes vector passed as argument.  The photon is assumed to travel along the (0,0,1)
        direction and have the reference axes for the Stokes parameter to be e_parallel = (1,0,0) and
        e_perp = (0,1,0)
        @param inStokesVector The Stokes vector 
    */

     
    PhotonCote(StokesV& inStokesVector);
    /*!
        @function init
        @discussion Initialize anything that is left to initialize, which for now means
        spitting out a warning to stdlog to make sure the user knows what he is using.
    */
    virtual void init() ;
    /*!
    @function MakeCopy
     @discussion Make a copy of a photon from the another Photon.

     @param inPhoton The other Photon
     */

    virtual Photon* MakeCopy();


    /*!
        @function MeasureStokesVectorInLabFrame
        @discussion Gives the Stokes vector with respect to the lab frame, by calculating
        it from the local frame of reference of the photon.
        @param      outS    The Stokes vector is returned via this parameter
    */

    virtual void MeasureStokesVectorInLabFrame(StokesV& outS, RealV& inEparaLab, RealV& inEperpLab, RealV& inNormalLab) ;
};

/*!
@class PhotonEvans
The class PhotonEvans implements the various actions done on a photon wavepacket according to

"Three Monte Carlo programs of polarized light 
transport into scattering media: part I"
13 June 2005 / Vol. 13,  No. 12 / OPTICS EXPRESS  4420

and is based on the work of Evans.  It simply sums the Stokes vector, with no rotation in a common reference
frame.  I don't see how this would be implemented in practice (i.e an experiment), but that's not relevant
for testing.
*/
class PhotonEvans : public Photon {
	
public:
    /*!
	@function PhotonEvans
	 @discussion Default constructor
	 
	 */
	
    PhotonEvans()  ;
    
    /*!
	@function PhotonEvans
	 @discussion This constructor makes a photon with polarization state given by 
	 the four Stokes parameters passed as arguments.  The photon is assumed to travel along the (0,0,1)
	 direction and have the reference axes for the Stokes parameter to be e_parallel = (1,0,0) and
	 e_perp = (0,1,0)
	 @param inI The intensity (first Stokes parameter) 
	 @param inQ The Q parameter (second Stokes parameter) 
	 @param inU The U parameter (third Stokes parameter)
	 @param inV The V parameter (fourth Stokes parameter)        
	 */
	
    PhotonEvans(double inI, double inQ, double inU, double inV, RealV inEpara , RealV inEprop); 
    /*!
        @function PhotonCote
	 @discussion This constructor makes a photon with polarization state given by 
	 the Stokes vector passed as argument.  The photon is assumed to travel along the (0,0,1)
	 direction and have the reference axes for the Stokes parameter to be e_parallel = (1,0,0) and
	 e_perp = (0,1,0)
	 @param inStokesVector The Stokes vector 
	 */
	
	
    PhotonEvans(StokesV& inStokesVector);
    /*!
        @function init
	 @discussion Initialize anything that is left to initialize, which for now means
	 spitting out a warning to stdlog to make sure the user knows what he is using.
	 */
    virtual void init() ;
    /*!
		@function MakeCopy
     @discussion Make a copy of a photon from the another Photon.
	 
     @param inPhoton The other Photon
     */
	
    virtual Photon* MakeCopy();
	
	
    /*!
        @function MeasureStokesVectorInLabFrame
	 @discussion This class does not perform any rotation of the reference frame.  Everything is simply added
				 up as discussed by J. Ramella Roman and Evans.
	 @param      outS    The Stokes vector is returned via this parameter
	 */
	
    virtual void MeasureStokesVectorInLabFrame(StokesV& outS, RealV& inEparaLab, RealV& inEperpLab, RealV& inNormalLab) ;
};


/*!
    @class PhotonMourant
    The class PhotonMourant implements the various actions done on a photon wavepacket according to
    
    "Characterizing mammalian cells and cell phantoms by 
    polarized backscatterng fiber-optic measurements"
    J. Mourant, T. Johnson, J. Freyer
    Appl. Opt. Vol 40, No8 Oct 2001, p. 5114
    
    The notation used in the implementation of this class is similar to the one used
 in the Bartel article.  The correspondance between Mourant and Bartel is:
    a = prel = el
    c = u = e3
    b = c cross a = -er
    No other comments exist for this class since it has not been tested and I am not using it.
*/

class PhotonMourant : public Photon {

public:
    /*!
    @function PhotonMourant
     @discussion Default constructor
	      */

    PhotonMourant() ;
    /*!
    @function PhotonMourant
     @discussion This constructor makes a photon with polarization state given by
     the four Stokes parameters passed as arguments.  The photon is assumed to travel along the (0,0,1)
     direction and have the reference axes for the Stokes parameter to be e_parallel = (1,0,0) and
     e_perp = (0,1,0)
     @param inI The intensity (first Stokes parameter)
     @param inQ The Q parameter (second Stokes parameter)
     @param inU The U parameter (third Stokes parameter)
     @param inV The V parameter (fourth Stokes parameter)
     */

    PhotonMourant(double inI, double inQ, double inU, double inV, RealV inEpara , RealV inEprop);
    /*!
    @function PhotonMourant
     @discussion This constructor makes a photon with polarization state given by
     the Stokes vector passed as argument.  The photon is assumed to travel along the (0,0,1)
     direction and have the reference axes for the Stokes parameter to be e_parallel = (1,0,0) and
     e_perp = (0,1,0)
     @param inStokesVector The Stokes vector
     */


    PhotonMourant(StokesV& inStokesVector);

    /*!
    @function init
     @discussion Initialize anything that is left to initialize, which for now means
     spitting out a warning to stdlog to make sure the user knows what he is using.
     */
    virtual void init() ;
    /*!
    @function MakeCopy
     @discussion Make a copy of a photon from the another Photon.

     @param inPhoton The other Photon
     */

    virtual Photon* MakeCopy();
    
    virtual void MeasureStokesVectorInLabFrame(StokesV& outS, RealV& inEparaLab, RealV& inEperpLab, RealV& inNormalLab);
    protected:

};



/*!
@class PhotonJaillon
The class PhotonJaillon implements the various actions done on a photon wavepacket according to

"Monte Carlo simulations of the diffuse backscattering
Mueller matrix for highly scattering media", Bartel & Hielscher
Appl. Opt. Vol 39, April 1st 2000, p. 1580.

with a few modifications in the way the Photon is detected following:

"Description and time reduction of a Monte Carlo code ...",
F. Jaillon, H. Saint Jalmes
Appl. Opt. Vol 42, no 16, p.3290

*/
class PhotonJaillon : public Photon {

public:
    /*!
     @function PhotonJaillon
     @discussion Default constructor
	*/

    PhotonJaillon();
    /*!
    @function PhotonJaillon
     @discussion This constructor makes a photon with polarization state given by
     the four Stokes parameters passed as arguments.  The photon is assumed to travel along the (0,0,1)
     direction and have the reference axes for the Stokes parameter to be e_parallel = (1,0,0) and
     e_perp = (0,1,0)
     @param inI The intensity (first Stokes parameter)
     @param inQ The Q parameter (second Stokes parameter)
     @param inU The U parameter (third Stokes parameter)
     @param inV The V parameter (fourth Stokes parameter)
     */

    PhotonJaillon(double inI, double inQ, double inU, double inV, RealV inEpara , RealV inEprop);
    /*!
    @function PhotonJaillon
     @discussion This constructor makes a photon with polarization state given by
     the Stokes vector passed as argument.  The photon is assumed to travel along the (0,0,1)
     direction and have the reference axes for the Stokes parameter to be e_parallel = (1,0,0) and
     e_perp = (0,1,0)
     @param inStokesVector The Stokes vector
     */


    PhotonJaillon(StokesV& inStokesVector);

    /*!
    @function init
     @discussion Initialize anything that is left to initialize, which for now means
     spitting out a warning to stdlog to make sure the user knows what he is using.
     */
    virtual void init() ;
    /*!
    @function MakeCopy
     @discussion Make a copy of a photon from the another Photon.

     @param inPhoton The other Photon
     */

    virtual Photon* MakeCopy();    
    /*!
    @function MeasureStokesVectorInLabFrame
     @discussion Gives the Stokes vector with respect to the lab frame, by calculating
     it from the local frame of reference of the photon.
     @param      outS    The Stokes vector is returned via this parameter
     */

    virtual void MeasureStokesVectorInLabFrame(StokesV& outS, RealV& inEparaLab, RealV& inEperpLab, RealV& inNormalLab) ;
};

#endif


