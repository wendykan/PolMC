#include "config.h"
#include "Photon.h"
#include "MCRandomScatterer.h"
#include "StokesV.h"
#include "MCObject.h"
#include "constants.h"
#include "mydebug.h"
#include "fastinterpolate.h"
#include <iostream>
#include <fstream>

#include "mtRand.h"
#define RandomNum() gRand.gendrand()

extern mtRand gRand;

int
main()
{
    cout.setf(ios::scientific);
    
    long failed = 0;

    cout << " === Test Fresnel === " << endl;
    double c;
    
    for (c = 1; c >= 0; c -= 0.01) {
        IntersectElement element;
        element.normal = RealV(1,0,0);
        element.cosine = c;
        element.indexFrom = 1.54;
        element.indexTo = 1.33;
        
        MCObject::SetFresnelCoefficients(element);
        
        cout << acos(c) << "\t" << element.Rs << "\t" <<  element.Rp << "\t" <<element.Ts << "\t" <<element.Tp << endl;

    }
    
    cout << " === Test displacement === " << endl;

    for (long i = 0; i < 100000; i++) {
        PhotonJaillon p(1,0,0,0);

        double phi = 2*PI*RandomNum();
        double theta = PI*RandomNum();
        double d = 1;
        p.RotateReferenceFrameAroundPropagationDirectionBy(phi);
        p.ChangePropagationDirectionAroundEPerpBy(theta);
        p.MoveBy(d);

        RealV position = p.GetLocalPosition();

        if (! areTheSame(position.abs(),d, -6) ) {
            cout << position.abs() << " != " <<  d  << endl;
        }

    }

    cout << " === Test vectors === " << endl;
    RealV x(1,0,0),y(0,1,0),z(0,0,1);
    RealV xpi4(1./sqrt(2.),1./sqrt(2.),0);
    RealV x3pi4(-1./sqrt(2.),1./sqrt(2.),0);
    RealV x5pi4(-1./sqrt(2.),-1./sqrt(2.),0);
    RealV x7pi4(1./sqrt(2.),-1./sqrt(2.),0);

    RealV t;
    t = RealV::CrossProduct(x,y);
    double r;
    
    if (t != z) {
        cout << "Failure: Cross product x * y != z " << t << " != " << z << endl ;
    } else {
        cout << "Success: Cross product x * y == z " << endl ;
    }
        
    if ( !areTheSame(r = RealV::OrientedAngleBetween(x,y,z), PI/2, 6)) {
        cout << "Angle " << r << " not " << PI/2 << endl;
        failed++;
    }
    if ( !areTheSame(r = RealV::OrientedAngleBetween(y,z,x),PI/2., 6)){
        cout << "Angle " << r << " not " << PI/2 << endl;
        failed++;
    }
    if ( !areTheSame(r = RealV::OrientedAngleBetween(z,x,y), PI/2., 6)){
        cout << "Angle " << r << " not " << PI/2 << endl;
        failed++;
    }
    if ( !areTheSame(r = RealV::OrientedAngleBetween(x,z,y),-PI/2., 6)){
        cout << "Angle " << r << " not " << PI/2 << endl;
        failed++;
    }
    if ( !areTheSame(r = RealV::OrientedAngleBetween(z,y,x), -PI/2., 6)){
        cout << "Angle " << r << " not " << PI/2 << endl;
        failed++;
    }
    if ( !areTheSame(r = RealV::OrientedAngleBetween(y,x,z),-PI/2., 6)){
        cout << "Angle " << r << " not " << PI/2 << endl;
        failed++;
    }

    if (failed)
        cout << "Failures of oriented angle calculation: " << failed << endl;
    else
        cout << "Success of oriented angle calculation." << endl;
    failed = 0;
    
    
    for (long i = 0; i < 100000; i++) {
        RealV v1(2*RandomNum()-1,2*RandomNum()-1,2*RandomNum()-1);
        RealV v2(2*RandomNum()-1,2*RandomNum()-1,2*RandomNum()-1);
        RealV v3,v4,v5,v6;
        
        RealV t1,t2;
        t1 = RealV::CrossProduct(v1,v2);
        t2 = RealV::NormalizedCrossProduct(v1,v2);

        if (!RealV::arePerpendicular(v1,t1) || !RealV::arePerpendicular(v2,t1) ) {
            failed++;
        }
        

        v3 = t1;
        v4 = -v3;
        v5 = t2;

        if ( !areTheSame(RealV::NormalizedDotProduct(v3,v4),-1, 6) ) {
            cout << "Mirror inversion failed " << -v3 << " != " << v4 << endl;
        }
        
        if ( !areTheSame(RealV::OrientedAngleBetween(v1,v2,v3), RealV::OrientedAngleBetween(v1,v2,v5), 6)) {
            cout << "Oriented angle failed: " << RealV::OrientedAngleBetween(v1,v2,v3) << " != " << RealV::OrientedAngleBetween(v1,v2,v4) <<  endl ;
            failed++;
        }

        if ( !areTheSame(RealV::OrientedAngleBetween(v1,v2,v3), -RealV::OrientedAngleBetween(v1,v2,v4), 6)) {
            cout << "Oriented angle failed: " << RealV::OrientedAngleBetween(v1,v2,v3) << " != " << -RealV::OrientedAngleBetween(v1,v2,v4) <<  endl ;
            failed++;
        }
        
    }

    if (failed)
        cout << "Failures of random cross product and normalized cross-product: " << failed << endl;
    else
        cout << "Success of random cross product and normalized cross-product." << endl;



    
    cout << " === Test frame of reference === " << endl;
    
    for (double theta = 0.; theta < 2*PI; theta+= 0.1) {
        StokesV S(1,1);
        S.RotateReferenceFrameAroundPropagationDirectionBy(theta);
        Complex El,Er;
        RealV el,er,ez;
        
        S.GetLocalComplexFields(El,Er);
        S.GetReferenceFrame(er,el,ez);
        CheckTriad_(el,er,ez);

        cout << S << endl;
        cout << " Ref: el= " << el << " and ";
        cout << "      er= " << er << endl;
        cout << " field along x = " << real(El)*el.x + real(Er)*er.x << endl;
        cout << " field along y = " << real(El)*el.y + real(Er)*er.y << endl;
        
    }


    cout << " === Test rotation of Stokes vector === " << endl;

    for (double theta = 0.; theta < 2*PI; theta+= 0.1) {
        StokesV S(1,0);
        S.RotatePolarizationStateBy(theta);
        Complex El,Er;
        RealV el,er,ez;

        S.GetLocalComplexFields(El,Er);
        S.GetReferenceFrame(er,el,ez);
        CheckTriad_(el,er,ez);

        cout << " Ref: el= " << el << " and ";
        cout << "      er= " << er << endl;
        cout << " field along x = " << real(El)*el.x + real(Er)*er.x << endl;
        cout << " field along y = " << real(El)*el.y + real(Er)*er.y << endl;

    }

    cout << " === Test change of direction of propagation === " << endl;

    for (double theta = 0.; theta < 2*PI; theta+= 0.1) {
        StokesV S(1,0);
        S.ChangePropagationDirectionAroundEPerpBy(theta);
        Complex El,Er;
        RealV el,er,ez;

        S.GetLocalComplexFields(El,Er);
        S.GetReferenceFrame(er,el,ez);
        CheckTriad_(el,er,ez);

        cout << " Ref: el= " << el << " and ";
        cout << "      er= " << er << " and ";
        cout << "      e3= " << ez << endl;
        cout << " field along x = " << real(El)*el.x + real(Er)*er.x << endl;
        cout << " field along y = " << real(El)*el.y + real(Er)*er.y << endl;
        cout << " field along z = " << real(El)*el.z + real(Er)*er.z << endl;

    }


    cout << " === Test change of direction of propagation === " << endl;

    for (double theta = 0.; theta < 2*PI; theta+= 0.1) {
        StokesV S(1,0);
        S.RotateReferenceFrameAroundPropagationDirectionBy(PI/4.);
        S.ChangePropagationDirectionAroundEPerpBy(theta);
        Complex El,Er;
        RealV el,er,ez;

        S.GetLocalComplexFields(El,Er);
        S.GetReferenceFrame(er,el,ez);
        CheckTriad_(el,er,ez);

        cout << " Ref: el= " << el << " and ";
        cout << "      er= " << er << " and ";
        cout << "      e3= " << ez << endl;
        cout << " field along x = " << real(El)*el.x + real(Er)*er.x << endl;
        cout << " field along y = " << real(El)*el.y + real(Er)*er.y << endl;
        cout << " field along z = " << real(El)*el.z + real(Er)*er.z << endl;

    }
    

}
