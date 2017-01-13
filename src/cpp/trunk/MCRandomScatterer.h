#include <map>
#include <string>
#include "MCMaterials.h"

#ifndef MCRANDOMSCATTERER_H
#define MCRANDOMSCATTERER_H
#define PI                      3.1415926535897932385
#define CHANCE 0.1
#define WeightThreshold 1e-6

class MuellerM;
class MuellerMSpherical;
class Photon;
class MCMaterial;

const long N_Points_stats = 90;
const double INFINITE_DISTANCE = 1000;

/*!
    @class MCRandomScatterer

    The class MCRandomScatterer provides the necessary encapsulation for manipulating random scatterers.
	It isolates the user from the details of how the random angles and random distances are calculated.
	The random scatterer can be a any kind of discrete scatterer (multiple implementations) or something else (like random
	dielectric fluctuations (currently not implemented, but "easily" added)).  One would also add any 
	fluorescence extension here (but some modifications need to be done to base class (i.e. wavelength of photon 
	must be considered).
	
	MCRandomScatterer provides the basis for several other classes, which are different sampling algorithms of a Mueller matrix.
	What is common to all of them (including their interface to the outside world) is declared here.
	
*/

class MCRandomScatterer {
public:
    MCRandomScatterer(double inMu_s,
                      double inMu_a,
                      double inIndexMedium,
                      unsigned long inSeed = 0);
    
    MCRandomScatterer(double inScattererRadius,
                      double inIndexScatterer,
                      double inIndexMedium,
                      double inWavelength,
                      long inPointsTable,
                      double inMu_s,
                      double inMu_a,
                      unsigned long inSeed = 0);
    MCRandomScatterer(map<string,string> inDict);

    virtual ~MCRandomScatterer();

    MCRandomScatterer& operator=(const MCRandomScatterer& inRhs);

    static double RandomFloat();

    void Roulette(Photon* ioPhoton);

    MuellerM* GetMuellerMatrix();
    MuellerMSpherical* GetMuellerMatrixSpherical();

    void initScattererProperties(double inScattererRadius,
                                 double inIndexScatterer,
                                 double inIndexMedium,
                                 double inWavelength,
                                 long inPointsTable = 10000);

    void initScattererDistributions(double inMu_s,
                                    double inMu_a);
    
    void initStats();

    virtual double GetRandomScatteringDistance(Photon* inPhoton);

    virtual void SetScatteringCoefficient(double inMu_s);
    virtual void SetAbsorptionCoefficient(double inMu_a);

    virtual double GetIndexMedium(Photon* inPhoton);
    virtual double GetScatteringCoefficient(Photon* inPhoton);
    virtual double GetAbsorptionCoefficient(Photon* inPhoton);
    virtual double GetTotalExtinctionCoefficient(Photon* inPhoton);

    virtual void GetRandomScatteringAngles(  double&    outTheta,
                                             double&    outPhi,
                                             Photon*   inPhoton);

    void	SetSeed(unsigned long inSeed);
    unsigned long GetSeed();
    
    void KeepStats();
    void DontKeepStats();
    void AddToHistograms(double inTheta, double inPhi);
    double GetSampledAnisotropy();
    virtual void DumpStatsToStream(ostream& out);
    
protected:
     /*! @var mMuellerM          	Pointer to a Mueller matrix*/
	MuellerM* mMuellerM;
     /*! @var mMuellerMSpherical	Pointer to a Mueller matrix for a spherical scatterer (assumes certain symmetry) */
    MuellerMSpherical* mMuellerMSpherical;
     /*! @var mMu_s	                Scattering coefficient */
    double mMu_s;
     /*! @var mMu_a	                Absorption coefficient */
    double mMu_a;
     /*! @var mMu_t	                Total extinction coefficient */
    double mMu_t;
     /*! @var mScattererRadius      Radius of the scattering particle */
    double mScattererRadius;
     /*! @var mIndexScatterer       Index of scattering particle */
    double mIndexScatterer;
    /*! @var mIndexMedium          Index of surrounding medium */
    double mIndexMedium;
    /*! @var mWavelength          Wavelength for which the Mueller matrix is calculated */
    double mWavelength;
    /* @var mPointsTable 	  Number of points for interpolation tables */
    long mPointsTable;
    
     /*! @var mUserWarned			Flag for remembering whether or not we have told the user who we are and what we do */
    bool mUserWarned;

	 /*! @var mSeed					Seed for the random number generator.  Can be set manually so we get the same sequence of numbers
									Useful for debugging. */
    unsigned long mSeed;
    
     /*! @var mKeepStats			Variable telling if we keep the stats on the generated angles.  Very useful for sanity checks. */
    bool mKeepStats;
     /*! @var mAvgCosine			Variable for the running anisotropy coefficient, from sampled angles.  Very useful for sanity checks. */
    double mAvgCosine;
     /*! @var mTotalCount			Variable for the number of generated angles.  Very useful for sanity checks. */
    long mTotalCount;
    
     /*! @var mHistTheta			Histogram for the scattering angle theta (in the plane of scattering). */
    long mHistTheta[N_Points_stats];
     /*! @var mHistPhi				Histogram for the angle phi (around propagation direction). */
    long mHistPhi[N_Points_stats];

    long m2DHistThetaPhi[N_Points_stats][N_Points_stats];
	
	/*! @var hasBeenSeeded			This is a class variable that is set to TRUE when the random number generator has been seeded */
	static bool hasBeenSeeded;

};


/*!
    @class MCRandomScattererKaplan

	Kaplan et al. "Mueller matrix of dense polystyrene latex sphere suspensions: 
	 measurements and Monte Carlo simulation", Applied Optics,  40, No 16, p. 2769, (2001)
		
*/

class MCRandomScattererKaplan :   public MCRandomScatterer {

public:
    MCRandomScattererKaplan(double inScattererRadius,
                      double inIndexScatterer,
                      double inIndexMedium,
                      double inWavelength,
                      long inPointsTable,
                      double inMu_s,
                      double inMu_a,
                      unsigned long inSeed=0);
    

    MCRandomScattererKaplan(map<string,string> inDict);

    ~MCRandomScattererKaplan();
    
    void init();

    virtual void DumpStats();
	virtual void DumpStatsToStream(ostream& out);

    void GetRandomScatteringAngles(  double&    outTheta,
                                     double&    outPhi,
                                     Photon*   inPhoton);

protected:

    /*! @var mS1_sqr          	Contains interpolating table for Mie scattering function S1_sqr */
    fastinterpolate mS1_sqr;
    /*! @var mS2_sqr           	Contains interpolating table for Mie scattering function S2_sqr */
    fastinterpolate mS2_sqr;
    /*! @var mT1			Contains integral of S1_sqr * sin(theta) from 0 to pi */
    double mT1;
    /*! @var mT2			Contains integral of S2_sqr * sin(theta) from 0 to pi */
    double mT2;

    /*! @var mInvF1_phi          	Contains interpolating table for function mInvF1_phi */
    fastinterpolate mInvF1_phi;
    /*! @var mInvF2_phi           	Contains interpolating table for function mInvF2_phi */
    fastinterpolate mInvF2_phi;

    /*! @var mInvCumulProbS1_sqr_SinTheta           Contains interpolating table for inverse cumulative probability distribution function of S1 */
    fastinterpolate mInvCumulProbS1_sqr_SinTheta;
    /*! @var mInvCumulProbS2_sqr_SinTheta           Contains interpolating table for inverse cumulative probability distribution function of S2 */
    fastinterpolate mInvCumulProbS2_sqr_SinTheta;


};

/*!
    @class MCRandomScattererHenyeyGreenstein

    The class MCRandomScattererHenyeyGreenstein is a random angle sampling routine based on the Henyey-Greenstein distribution.
	The constructor can be initialized either with a g parameter or from a Mueller matrix, from which is extracted the parameter g.
*/

class MCRandomScattererHenyeyGreenstein :  public MCRandomScatterer {

public:
    MCRandomScattererHenyeyGreenstein(map<string,string> inDict);

    MCRandomScattererHenyeyGreenstein(double inAnisotropy,
                                      double inMu_s,
                                      double inMu_a,
                                      double inIndexMedium,
                                      unsigned long inSeed = 0);

    virtual ~MCRandomScattererHenyeyGreenstein();
    
    void init();
    void init(double inAnisotropy);
    
    void GetRandomScatteringAngles(  double&    outTheta,
                                     double&    outPhi,
                                     Photon*   inPhoton);
    virtual void DumpStatsToStream(ostream& out);

protected:

    /*! @var mG          	Contains the computed anisotropy parameter*/
    double mG;
};

/*!
    @class MCRandomScattererJaillon

    The class MCRandomScattererJaillon is a random angle sampling routine based on the Jaillon
	Jaillon et al. "Description and time reduction of a Monte Carlo code to simulate 
	propagation of polarized light through scattering media", Applied Optics,  42, No 16, p. 3290, (2003)

*/

class MCRandomScattererJaillon :  public MCRandomScatterer {

public:
    MCRandomScattererJaillon(double inScattererRadius,
                            double inIndexScatterer,
                            double inIndexMedium,
                            double inWavelength,
                            long inPointsTable,
                            double inMu_s,
                            double inMu_a,
                             unsigned long inSeed);
    
    MCRandomScattererJaillon(map<string,string> inDict);

    void init();

    void GetRandomScatteringAngles(  double&    outTheta,
                                     double&    outPhi,
                                     Photon*   inPhoton);
    virtual void DumpStatsToStream(ostream& out);

protected:
   /*! @var mInvCumulComparison          	Contains the interpolating table for the inverse cumulative probability distribution */
      fastinterpolate   mInvCumulComparison;


};

/*!
@class MCRandomScattererMaterial

General class that obtains material properties from a MCMaterial.


*/

class MCRandomScattererMaterial :  public MCRandomScattererHenyeyGreenstein {
    
public:
    MCRandomScattererMaterial(MCMaterial& inMaterial);
    virtual ~MCRandomScattererMaterial();
    
    void SetRelativeConcentration(double inConcentration);
    virtual double GetIndexMedium(Photon* inPhoton);
    virtual double GetTotalExtinctionCoefficient(Photon* inPhoton);
    virtual double GetScatteringCoefficient(Photon* inPhoton);
    virtual double GetAbsorptionCoefficient(Photon* inPhoton);
    virtual double GetAnisotropyCoefficient(Photon* inPhoton);
    void GetRandomScatteringAngles(  double&    outTheta,
                                     double&    outPhi,
                                     Photon*   inPhoton);
	virtual void DumpStatsToStream(ostream& out);

protected:
        /*! @var mMaterial  Material property object. */
        MCMaterial mMaterial;
    
        double cWavelength;
        double cBackgroundIndex;
        double cScatteringCoefficient;
        double cAbsorptionCoefficient;
        double cAnisotropyCoefficient;
        double cExtinctionCoefficient;
 
        bool cCacheScat;
        bool cCacheAbs;
        bool cCacheAnis;
        bool cCacheIndex;
        bool cCacheExt;
        
    
};
#endif