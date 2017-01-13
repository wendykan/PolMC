/*!
    @header RealV
    @abstract   Simple class that represents a real 3D vector
    @discussion This class represnets a real vector with numerous functions that can be applied onto vectors, such as dot product, cross product, rotation of vectors in 3D etc...
				Utility functions to read and write RealV to screen or to disk are provided using the simple << and >> operators form C++.
*/

#ifndef __REALV__H
#define __REALV__H

#include "mydebug.h"
#include "constants.h"
#include <iostream>

using namespace std;

/*!
@defined CheckTriad_
 Provided three RealV vectors, in order x,y and z, will warn the user
 if does not form an orthonormal set and debug mode is active (see mydebug.h)  
 */
#define CheckTriad_(x, y, z) MyAssert_(! RealV::CheckTriad(x, y, z));

class RealV;

ostream& operator<<(ostream&s, const RealV& v);
istream& operator>>(istream&s, RealV& v);

#define DotProduct_(u, v) (u.x * v.x + u.y * v.y + u.z * v.z)
#define CrossProduct_(u,v) (RealV(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x))

/*!
@class RealV
Provides encapsulation for a real three-dimensional vector.
*/

class RealV {
public:
    /*!
    @function RealV
     @discussion Default constructor.  Vector set to (0,0,0).
     */
    RealV();
    /*!
    @function RealV
     @discussion Constructor with x,y,z arguments.  Vector set to (inX, inY, inZ).
     @param inX  x cartesian coordinate (or projection)
     @param inY  y cartesian coordinate (or projection)
     @param inZ  z cartesian coordinate (or projection)
     */
    RealV(double inX, double inY, double inZ);
    /*!
        @function RealV
     @discussion Constructor from other RealV .
     @param inRhs implicit argument corresponding to right hand side
     */
    RealV(const RealV& inRhs);    
    /*!
        @function operator=
     @discussion Copy operator to assign vector to another
     @param inRhs implicit argument corresponding to right hand side
     @result The current vector, to which the the other vector has been added
     */
    
    RealV& operator=(const RealV& inRhs) ;
    
    /*!
        @function operator==
     @discussion Comparison operator to compare a vector to another
     @param inRhs implicit argument corresponding to right hand side
     @result true is vectors are identical
     */
    
    bool operator==(const RealV& inRhs);
    
    /*!
        @function operator!=
     @discussion Comparison operator to compare a vector to another
     @param inRhs implicit argument corresponding to right hand side
     @result true is vectors are different
     */
    
    bool operator!=(const RealV& inRhs);    
    
    
    /*!
        @function operator+=
     @discussion Addition operator to add two vectors together
     @param inRhs implicit argument corresponding to right hand side
     @result The current vector, to which the the other vector has been added
     */
    
    inline RealV& operator+=(RealV inRhs)
    {
        x += inRhs.x;
        y += inRhs.y;
        z += inRhs.z;
        
        return *this;
    }
    
    /*!
        @function operator-=
     @discussion Subtraction operator to subtract two vectors together
     @param inRhs implicit argument corresponding to right hand side
     @result The current vector, to which the the other vector has been subtracted
     */
    
    inline RealV& operator-=(RealV inRhs)
    {
        x -= inRhs.x;
        y -= inRhs.y;
        z -= inRhs.z;
        
        return *this;
    }
    
    /*!
        @function operator*=
     @discussion Multiplication operator to multiply a vector by a real constant
     @param inRhs implicit argument corresponding to the constant
     @result The current vector, which has been multiplied by the constant
     */
    
    RealV& operator *= (double inRhs) ;
    
    /*!
        @function operator/=
     @discussion Division operator to divide a vector by a real constant
     @param inRhs implicit argument corresponding to the constant
     @result The current vector, which has been divided by the constant
     */
    
    RealV& operator /= (double inRhs) ;    
    
    /*!
        @function abs
     @discussion Returns the magnitude of the vector
     @result The magnitude sqrt(x*x + y*y + z*z)
     */
    
    double abs () ;
    /*!
        @function norm
     @discussion Returns the norm of the vector (the square of the magnitude)
     @result The x*x + y*y + z*z
     */
    
    double norm () ;    
    /*!
        @function normalize
     @discussion Normalizes the vector
     */
    
    void normalize ();    
    /*!
        @function DotProduct
     @discussion Computes the dot product of two real vectors
     @result The dot product of the two vectors
     */

    static double DotProduct(const RealV& u, const RealV& v) ;
    
    /*!
        @function NormalizedDotProduct
     @discussion Computes the dot product of two real vectors, normalized. Will ensure
     the resulting number is always less than one.  This means the results can safely be used
     in acos() and asin().  Without this test, the resulting number may be above one in rare situations (1 in 10^8).
     
     @result The cosine of the angle between the two vectors.
     */
    
    static double NormalizedDotProduct(RealV& u, RealV& v) ;
    
    /*!
        @function CrossProduct
     @discussion Computes the cross product of two real vectors
     @result The cross product of the two vectors
     */
    
    static RealV CrossProduct(RealV& u, RealV& v);
    
    /*!
        @function NormalizedCrossProduct
     @discussion Computes the cross product of two real vectors and normalizes the results.  Will ensure
     the magnitude of the resulting vector is always less than one.  This means the results can safely be used
     in acos() and asin().  Without this test, the resulting vector may be above one in rare situations (1 in 10^8).
     @result The cross product of the two vectors.  The amplitude is the sine of the angle between the two.
     */
    
    static RealV NormalizedCrossProduct(RealV& u, RealV& v) ;
    
    /*!
        @function TripleProduct
     @discussion Computes the triple product of three real vectors (a x b) . c
     @result The triple product of the two vectors
     */
    
    static double TripleProduct(RealV& u, RealV& v, RealV& w);
    
    /*!
        @function OrientedAngleBetween
     @discussion Computes the angle between two vectors u and v by which vector u needs to be
     rotated by around a third vector w to give the second vector v.
     @result The angle.  The sign is important.
     */
    
    static double OrientedAngleBetween(RealV& u, RealV& v, RealV& w) ;
    
    /*!
        @function areParallel
     @discussion Returns whether or not two vectors are parallel
     @result True if vectors are parallel
     */
    
    static bool areParallel(RealV& u, RealV& v) ;
    
    /*!
        @function arePerpendicular
     @discussion Returns whether or not two vectors are perpendicular
     @result True if vectors are perpendicular
     */
    
    static bool arePerpendicular(RealV& u, RealV& v);
    
    /*!
        @function CheckTriad
     @discussion Given three vectors x, y, z, checks whether or not
     they form an orthonormal set
     @result Non zero if triad is not orthonormal.
     */
    
    static bool CheckTriad(RealV& x, RealV& y, RealV& z) ;
    
	void RotateAroundX(double inPhi);
    
	void RotateAroundY(double inPhi);
	
	/*! @function RotateAroundZ
        @discussion Rotation around the z axis of the lab frame.
        
                   /  Cos phi -Sin phi 0 \
        Rz(phi)  = |  Sin phi Cos phi 0 |
                   \     0       0    1 /
        */
	
    void RotateAroundZ(double inPhi);
    
	void RotateAroundVector();
	
	static RealV RotateAroundVector(double inPhi, RealV vector, RealV axis);

	double DistanceToPlane(const RealV& origin, const RealV& normal, const RealV& d );

public:
	/*! @var x  x Cartesian coordinate of vector */
	double x;
    /*! @var y  y Cartesian coordinate of vector */
    double y;
    /*! @var z  z Cartesian coordinate of vector */
    double z;
};

/*
RealV operator+(const RealV& inLhs, const RealV& inRhs);
RealV operator-(const RealV& inLhs, const RealV& inRhs);
RealV operator+(RealV inRhs);
RealV operator-(RealV inRhs);

RealV operator*(RealV inVec, double inFloat);
RealV operator*(double inFloat, RealV inVec);
RealV operator/(RealV inVec, double inFloat);
*/

inline RealV operator+(const RealV& inLhs, const RealV& inRhs)
{
    return RealV(inLhs.x+inRhs.x,inLhs.y+inRhs.y,inLhs.z+inRhs.z);
}


inline RealV operator-(const RealV& inLhs, const RealV& inRhs)
{
    return RealV(inLhs.x-inRhs.x,inLhs.y-inRhs.y,inLhs.z-inRhs.z);
}


inline RealV operator-(RealV inRhs)
{
    return (inRhs *= -1);
}


inline RealV operator*(RealV inVec, double inFloat)
{
    return (inVec *= inFloat);
}


inline RealV operator*(double inFloat, RealV inVec)
{
    return (inVec *= inFloat);
}


inline RealV operator/(RealV inVec, double inFloat)
{
    return (inVec /= inFloat);
}

#endif
