
#if HAVE_CONFIG_H
    #include <config.h>
#endif

#include <iostream>

#if HAVE_UNISTD_H
    #include <unistd.h>
#endif

#include <ctime>
#include <limits>
#include <cmath>


#include "mydebug.h"

#define SBUF 255

using namespace std;

long gDebugLevel = kNormal;

const char* RightNow()
{
        static char s1[SBUF];
        static char s2[SBUF];
        struct tm *date;
        time_t now = time(NULL);
        date = localtime(&now);

        strftime(s1, SBUF, "%x %H:%M:%S",date);
#if HAVE_UNISTD_H
        sprintf(s2," [%d] ", getpid());
        return strcat(s1,s2);
#else
        return s1;
#endif

}


bool CheckArray( const std::valarray<Complex>& a, long l, char* s) 
{
    bool err = false;
    
        for (long i = 0; i < a.size(); i++) {
            if ( ! isfinite(a[i].real()) ) {
                std::cerr << "Value not finite (" << a[i].real() << ") in array at index " << i << " check line " << l << " in file " << s << std::endl; 
                err = true;
            }
            
            if  ( ! isfinite(a[i].imag())) {
                std::cerr << "Value not finite (" << a[i].imag() << ") in array at index " << i << " check line " << l << " in file " << s << std::endl; 
                err = true;
            } 
    
            if ( isnan(a[i].real()) ) {
                std::cerr << "Value NaN (" << a[i].real() << ") in array at index " << i << " check line " << l <<" in file " << s <<  std::endl; 
                err = true;
            }   
                
            if  ( isnan(a[i].imag())) {
                std::cerr << "Value NaN (" << a[i].imag() << ") in array at index " << i << " check line " << l << " in file " << s << std::endl; 
                err = true;
            } 
        }       
    
    return err;
}

bool CheckIfArrayReal( const std::valarray<Complex>& a, long l, char* s) 
{
    bool err = false;
    
        for (long i = 0; i < a.size(); i++) {
            if ( 1e-7 < abs( Double(a[i].imag()/a[i].real())) ) {
                std::cerr << "Imaginary part significant (" << a[i] << ") in array at index " << i << " check line " << l << " in file " << s << std::endl; 
                err = true;
            }
        }       
    
    return err;
}

bool CheckRealArray( const std::valarray<Double>& a, long l, char* s) 
{
    bool err = false;
    
        for (long i = 0; i < a.size(); i++) {
            if ( ! isfinite(a[i]) ) {
                std::cerr << "Value not finite (" << a[i] << ") in array at index " << i << " check line " << l << " in file " << s << std::endl; 
                err = true;
            }
            
            if ( isnan(a[i]) ) {
                std::cerr << "Value NaN (" << a[i] << ") in array at index " << i << " check line " << l << " in file " << s << std::endl; 
                err = true;
            }   
                
        }       
    
    return err;
}

bool CheckIfValidXForFFT( const std::valarray<Double>& a, long l, char* s) 
{
    bool err = false;
    
    if (a[0] != 0) {
        std::cerr << "First value not zero. Called from line " << l << " in file " << s << std::endl; 
        err = true;
    }
    
    for (long i = 1; i < a.size() / 2; ++i) {
        if ( a[i] != -a[a.size() - i]) {
            std::cerr << "Array not symmetric. Called from line " << l << " in file " << s << std::endl;
            err = true;
        }
        
    }
    
    return err;
}

bool CheckIfHermitian( const std::valarray<Double>& x, const std::valarray<Complex>& y, long l, char* s) 
{
    bool err = false;
    
    std::valarray<short> checked(x.size());
    
    for (long i = 0; i < x.size(); ++i) {
        
        if (x[i] == 0.) {
            if ( abs(Double(y[i].imag()/y[i].real())) > 1e-5) {
                std::cerr << "Warning: for value x = 0, y is complex " << y[i] << ". Called from line " << l << " in file " << s << std::endl;
                checked[i] = true;
                err = true;
            } else {
                checked[i] = true;
                continue;
            }
        }

        if (x[i] == x.max()) {
            checked[i] = true;
            continue;
        }
        
        long j;
        for (j = 0; j < x.size(); ++j) {
            if (j == i )
                continue;


            if (x[j] == - x[i]) {
                if ( ArgDiff(y[i],conj(y[j])) > 3 && AbsDiff(y[i],conj(y[j])) > 3)  {
                    std::cerr << "Array not hermitian at x[i] = " << x[i] << " << : ["<< i << "]["<< j <<"] :" << y[i] << " != " << y[j] << "\nCalled from line " << l << " in file " << s << std::endl;
                    std::cerr << x[i] << " and " << x[j] << std::endl;
                    err = true;
                } else {
                    checked[i] = true;
                    checked[j] = true;

                    break;
                }
            }
        }
    }

    if (checked.sum() != checked.size()) {
        std::cerr << checked.sum() << " != " << checked.size() << std::endl;
        std::cerr << "Some values of x do not have a corresponding -x. Called from line " << l << " in file " << s << std::endl;
        err = true;
    }
    
    return err;
}


double ArgDiff(Complex a, Complex b)
{
    double aa = arg(a);
    double ab = arg(b);

    if (aa >= ab && ab != 0) {
        return log10(aa / ab) ;
    } else if (aa <= ab && aa != 0) {
        return log10(ab / aa ) ;
    } else if (aa == 0 && ab != 0) {
        return log10(ab);
    } else if (ab == 0 && aa != 0) {
        return log10(aa);
    } else {
        return 0;
    }
}

double AbsDiff(Complex a, Complex b)
{
    double aa = abs(a);
    double ab = abs(b);

    if (aa >= ab && ab != 0) {
        return log10(aa / ab) ;
    } else if (aa <= ab && aa != 0) {
        return log10(ab / aa ) ;
    } else if (aa == 0 && ab != 0) {
        return log10(ab);
    } else if (ab == 0 && aa != 0) {
        return log10(aa);
    } else {
        return 0;
    }
}

bool CheckValue( const Complex & a, long l, char* s) 
{
    bool err = false;

        if ( ! isfinite(a.real()) ) {
            std::cerr << "Value not finite (" << a.real() << ") check line " << l << " in file " << s << std::endl; 
            err = true;
        }
        
        if  ( ! isfinite(a.imag())) {
            std::cerr << "Value not finite (" << a.imag() << ") check line " << l << " in file " << s << std::endl; 
            err = true;
        } 

        if ( isnan(a.real()) ) {
            std::cerr << "Value NaN (" << a.real() << ") check line " << l <<" in file " << s <<  std::endl; 
            err = true;
        }
        
        if  ( isnan(a.imag())) {
            std::cerr << "Value NaN (" << a.imag() << ") check line " << l << " in file " << s << std::endl; 
            err = true;
        } 
    
    return err;
    
}

bool CheckDoubleValue( const Double & a, long l, char* s) 
{
    bool err = false;

        if ( ! isfinite(a) ) {
            std::cerr << "Value not finite (" << a << ") check line " << l << " in file " << s << std::endl; 
            err = true;
        }
        
        if ( isnan(a) ) {
            std::cerr << "Value NaN (" << a << ") check line " << l << " in file " << s << std::endl; 
            err = true;
        }

    return err;
    
}
