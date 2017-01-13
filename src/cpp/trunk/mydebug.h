// If the following line is commented out, debug is deactivated.

/*!
    @header mydebug
    @abstract   Utility functions for debugging calculations
    @discussion Multiple utility functions are defined here to check the validity of numbers, complex numbers, check for difference, etc...
				
				Utility functions are defined to help display meaningful messages to the user.
 
				If debug is activated (MYDEBUG_H defined), the the macros ending in _() will perform their checks. If MYDEBUG_H is not defined, then
				the macros will not generate any code and will not slow down the code.
 
				
*/

//#define __MYDEBUG

#ifndef __DEBUG_H
#define __DEBUG_H

#include <iostream>
#include <complex>
#include <valarray>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#if HAVE_CONFIG_H
    #include <config.h>
#endif

#define ThrowRuntimeError(message)  { ostringstream s; s << message ; throw runtime_error(s.str()); }

// Use these macros (ending with _, which will call the real function if debug is activated

#ifdef __MYDEBUG
    #define MyAssert_(a)            if ( !(a) ) std::cerr << "Assertion failed at line " << __LINE__ << ", file "<<  __FILE__ << std::endl;
    #define CheckValue_(a)          CheckValue(a, long(__LINE__), __FILE__)
    #define CheckDoubleValue_(a)    CheckDoubleValue(a, long(__LINE__), __FILE__)
    #define CheckArray_(a)          CheckArray(a, long(__LINE__), __FILE__)
    #define CheckRealArray_(a)      CheckRealArray(a, long(__LINE__), __FILE__)
    #define CheckIfArrayReal_(a)    CheckIfArrayReal(a, long(__LINE__), __FILE__)
    #define CheckIfValidXForFFT_(a) CheckIfValidXForFFT(a, long(__LINE__), __FILE__)
    #define CheckIfHermitian_(a,b)  CheckIfHermitian(a,b,long(__LINE__), __FILE__)
	#define PrintMessage_(toPrint)  std::cerr << RightNow() << toPrint << " on line " << __LINE__ << ", file "<<  __FILE__ << std::endl; 
	#define PrintMessageIfLevel_(toPrint, level)  { if (level <= gDebugLevel ) { std::clog << RightNow() << toPrint << " on line " << __LINE__ << ", file "<<  __FILE__ << std::endl; } }
    #define IfThenPrint_(cond, toPrint)  if ( cond ) { std::log << RightNow() << toPrint << " on line " << __LINE__ << ", file "<<  __FILE__ << std::endl; }
    #define IfThen_(cond, do)        if ( cond) { do; }
    #define ThrowIf_(cond, message)  if ( cond ) { ostringstream s; s << "exception thrown from line " << __LINE__ << ", file "<<  __FILE__; \
                                                   throw runtime_error(s.str()); }
#else
    #define MyAssert_(a)  ;
    #define CheckValue_(a)          ;
    #define CheckDoubleValue_(a)    ;
    #define CheckArray_(a)          ;
    #define CheckRealArray_(a)      ;
    #define CheckIfArrayReal_(a)    ;
    #define CheckIfValidXForFFT_(a) ;
    #define CheckIfHermitian_(a,b)  ;
	#define PrintMessage_(toPrint)  ;
	#define PrintMessageIfLevel_(toPrint, level)  ;
    #define IfThenPrint_(cond, do)  ;
    #define IfThen_(cond, do)       ;
    #define ThrowIf_(cond, message) ;
#endif


extern long gDebugLevel;

enum {
	kQuiet,
	kNormal,
	kVerbose,
	kExtremelyVerbose,
	kInsane
};

/*!
    @typedef Double
    @discussion Type double
*/
typedef double Double;
/*!
    @typedef complex_float
    @discussion Complex number, single precision
*/
typedef std::complex<float> complex_float;
/*!
    @typedef complex_double
    @discussion Complex number, double precision
*/
typedef std::complex<double> complex_double;
/*!
    @typedef Complex
    @discussion Complex number, double precision
*/
typedef std::complex<Double> Complex;


#ifndef NAN
    /*!
     @defined NAN
     @discussion "Not a number", can be used as a number.
     */

    #define NAN numeric_limits<double>::signaling_NaN()
#endif


#ifndef isnan
    /*!
     @defined isnan
     @discussion If isnan(x) does not exist, defines a function
     that always returns false so that we can use isnan() on platforms that have it.
     */

    #define isnan(x) 0
#endif

#ifndef isfinite 
    /*!
     @defined isfinite
     @discussion If isfinite(x) does not exist, fake one by comparing value to a very large value.
     */

    #define isfinite(x) (x < 1e300)  
#endif


/*!
    @function RightNow
    @discussion Returns current time and date, for use typically to output to stderr or stdlog
    @result     Pointer to a static buffer with the time and date.
*/
const char* RightNow();

/*!
    @function CheckArray
    @discussion This function should not be called directly.  See the macro CheckArray_().
    Check to see if complex array contains elements that are infinite or NAN.  Will output a warning
    to stdlog if any, and will state the file and line number where the function was called.
    @param  a    A valarray of complex numbers
    @param  l Line number from which function was called
    @param  s Filename from which function was called
    @result     Non zero if array contains elements that are infinite or NAN.
*/

bool CheckArray( const std::valarray<Complex>& a, long l, char* s);

/*!
    @function CheckIfArrayReal
    @discussion This function should not be called directly.  See the macro CheckIfArrayReal_().
    Check to see if complex array contains elements that have a non zero imaginary part.  Will output a warning
    to stdlog if any, and will state the file and line number where the function was called.
    @param  a    A valarray of complex numbers
    @param  l Line number from which function was called
    @param  s Filename from which function was called
    @result     Non zero if array contains elements that have a non zero imaginary part.
*/

bool CheckIfArrayReal( const std::valarray<Complex>& a, long l, char* s);

/*!
    @function CheckRealArray
    @discussion This function should not be called directly.  See the macro CheckIfArrayReal_().
    Check to see if real array contains elements that are infinite or NAN.  Will output a warning
    to stdlog if any, and will state the file and line number where the function was called.
    @param  a   A valarray of real numbers
    @param  l Line number from which function was called
    @param  s Filename from which function was called
    @result     Non zero if array contains elements that are infinite or NAN.
*/

bool CheckRealArray( const std::valarray<Double>& a, long l, char* s);

/*!
    @function CheckValue
    @discussion This function should not be called directly.  See the macro CheckValue_().
    Check to see if complex number is infinite or NAN.  Will output a warning
    to stdlog if any, and will state the file and line number where the function was called.
    @param  a    A complex number
    @param  l Line number from which function was called
    @param  s Filename from which function was called
    @result     Non zero if number is infinite or NAN.
*/

bool CheckValue( const Complex& a, long l, char* s);

/*!
    @function CheckDoubleValue
    @discussion This function should not be called directly.  See the macro CheckDoubleValue_().
    Check to see if real number is infinite or NAN.  Will output a warning
    to stdlog if any, and will state the file and line number where the function was called.
    @param  a    A real number
    @param  l Line number from which function was called
    @param  s Filename from which function was called
    @result     Non zero if number is infinite or NAN.
*/

bool CheckDoubleValue( const Double& a, long l, char* s);

/*!
    @function CheckIfValidXForFFT
    @discussion This function should not be called directly.  See the macro CheckIfValidXForFFT_().
    Check to see if a array of real numbers is appropriate for use with the FFT algorithm, that is, if
    the first element is zero, and if these elements are symmetric with respect to zero (i.e. a[i] == -a[a.size() - i]).
    
    @param  a An array of real number
    @param  l Line number from which function was called
    @param  s Filename from which function was called
    @result     Non zero if invalid.
*/

bool CheckIfValidXForFFT( const std::valarray<Double>& a, long l, char* s);

/*!
    @function CheckIfHermitian
    @discussion This function should not be called directly.  See the macro CheckIfHermitian_().
    Check to see if array of complex numbers y is hermitian, assuming pairs (x,y). 
    Will output a warning to stdlog if any, and will state the file and line number where the function was called.
	@param  x An array of complex numbers
	@param  x A second array of complex numbers
    @param  l Line number from which function was called
    @param  s Filename from which function was called
    @result     Non zero if not hermitian.
*/

bool CheckIfHermitian( const std::valarray<Double>& x, const std::valarray<Complex>& y, long l, char* s);

/*!
    @function ArgDiff
    @discussion Computes the difference between the argument of two complex numbers (the angle in Euler representation).
    @param      a First complex number
    @param      b Second complex number
    @result     Difference between the angles
*/

double ArgDiff(Complex a, Complex b);

/*!
 @function AbsDiff
 @discussion COmputes the difference between the modulus of two complex numbers (the angle in Euler representation).
 @param      a First complex number
 @param      b Second complex number
 @result     Difference between the absoluate values
 */

double AbsDiff(Complex a, Complex b);


#define areTheSame(x,y, order)	( (x) == (y) || ( (x)*(y) >= 0 && -log10( std::abs( ((x)/(y))-1. ) )  >= order ) )
#define isZero(x, order)	( std::abs( log10( std::abs(x) ) ) >= order )



#endif
