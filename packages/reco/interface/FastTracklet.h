/*

Supporting structure of fast tracking, include SignedHit and Tracklet

Author: Kun Liu, liuk@fnal.gov
Created: 06-09-2013

*/

#ifndef _FASTTRACKLET_H
#define _FASTTRACKLET_H

#include <GlobalConsts.h>

#include <list>
#include <vector>

#include <TVector3.h>

#include <geom_svc/GeomSvc.h>
#include <phool/PHObject.h>

#include "SRawEvent.h"
#include "SRecEvent.h"

class SignedHit : public PHObject
{
public:
    SignedHit();
    explicit SignedHit(int detectorID);
    SignedHit(Hit hit_input, int sign_input);

    //PHObject virtual overloads
    void identify(std::ostream& os = std::cout) const;
    void Reset() { hit = Hit(); sign = 0; }
    int  isValid() const { return hit.index > 0; }
    SignedHit* Clone() const { return (new SignedHit(hit, sign)); }

    //comparision operators for sorting
    bool operator<(const SignedHit elem) const { return hit.detectorID < elem.hit.detectorID; }
    bool operator==(const SignedHit elem) const { return hit.index == elem.hit.index; }

    //Get the real hit position
    double pos() { return hit.pos + sign*hit.driftDistance; }
    double pos(int sign_input) { return hit.pos + sign_input*hit.driftDistance; }

    //Data members
    Hit hit;
    int sign;

    ClassDef(SignedHit, 1)
};

class PropSegment : public PHObject
{
public:
    PropSegment();

    //PHObject virtual overloads
    void identify(std::ostream& os = std::cout) const { print(os); };
    void Reset() { init(); }
    int  isValid() const;
    PropSegment* Clone() const { return (new PropSegment(*this)); }

    //init -- temporary, only used for tests
    void init();

    //Quality cut
    //bool isValid();

    //Debugging output
    void print(std::ostream& os = std::cout) const;

#ifndef __CINT__
    //Get expected position at a given Z
    double getExpPosition(double z) const { return a*z + b; }

    //Get the closest approach to a given space position/proptube wire
    double getClosestApproach(double z, double pos);

    //Get reference pos at first two planes
    double getPosRef(double default_val = -9999.);

    //Number of hits
    int getNHits() const;

    //Number of planes
    int getNPlanes() const;

    //Fit the segment -- naive linear fit
    void fit();   // external call
    void fit_2hits();
    void fit_34hits();

    //linear chisq fitter
    void linearFit_simple();
    void linearFit_iterative();

    //resolve left/right
    void resolveLR();
    void resolveLR(int setting);
#endif

    //track slope the interception
    double a;
    double b;
    double err_a;
    double err_b;

    //chisq of the segment
    double chisq;

    //Auxilary hodoscope hit list, cannot be possibly more than 4 hodoscope hits
    int nHodoHits;
    Hit hodoHits[4];

    //Hit list -- only 4 hits at most
    SignedHit hits[4];

    ClassDef(PropSegment, 4)
};

class Tracklet : public PHObject
{
public:
    Tracklet();

    //PHObject virtual overloads
    void identify(std::ostream& os = std::cout) const { os << "Tracklet @sID=" << stationID << std::endl;}
    void Reset() { *this = Tracklet(); }
    int  isValid() const;
    Tracklet* Clone() const { return (new Tracklet(*this)); }

    //Basic quality cut
    //bool isValid();

    //Debuggin output
    void print(std::ostream& os = std::cout);

#ifndef __CINT__
    //Sort hit list
    void sortHits() { hits.sort(); }

    //Update/get number of real hits
    int getNHits() const { return nXHits + nUHits + nVHits; }

    //Number of all hits (even excluded)
    int getNAllHits() { return hits.size(); }

    //Get the probabilities
    double getProb() const;

    //Get the momentum probabilities
    double getMomProb() const;

    //Get the chi square
    double getChisq() const { return chisq; }

    //Get x and y positions at a given z
    double getExpPositionX(double z) const;
    double getExpPosErrorX(double z) const;
    double getExpPositionY(double z) const;
    double getExpPosErrorY(double z) const;
    double getExpPositionW(int detectorID) const;
    int    getExpElementID(int detectorID) const;

    //Get momentum upstream/downstream
    TVector3 getMomentumSt1() const;
    TVector3 getMomentumSt3() const;
    TVector3 getExpMomentum(double z) const;

    //Get the i-th signed hit
    SignedHit getSignedHit(int index);

    //Kernal function to calculate chi square for minimizer
    double Eval(const double* par);
    double Eval4(const double* par);
    double calcChisq();
    double calcChisq_noDrift();

    //Add dummy hits
    void addDummyHits();

    //Momentum estimation using back partial
    double getMomentum() const;

    //Decide charge by KMag bending direction
    //int getCharge() const { return x0*KMAGSTR > tx ? 1 : -1; }
    int getCharge() const;

    //Get the slope and intersection in station 1
    void getXZInfoInSt1(double& tx_st1, double& x0_st1) const;
    void getXZErrorInSt1(double& err_tx_st1, double& err_x0_st1) const;

    //For sorting tracklet list
    bool operator<(const Tracklet& elem) const;

    //For reducing similar tracklets
    bool similarity(const Tracklet& elem) const;

    //Merge the hit list from two tracklets
    Tracklet merge(Tracklet& elem);

    //For adding two tracklets together to form a back partial track
    Tracklet operator+(const Tracklet& elem) const;

    //For adding two tracklets together to form a global track
    Tracklet operator*(const Tracklet& elem) const;

    //Convert to a SRecTrack
    SRecTrack getSRecTrack(bool hyptest = true);
#endif

    //Station ID, ranging from 1 to nStation, nStation-1 means back partial track, nStation means global track
    int stationID;

    //Number of hits
    mutable int nXHits;
    mutable int nUHits;
    mutable int nVHits;

    //Chi square
    double chisq;

    //chisq at vertex
    double chisq_vtx;

    //List of signed hits
    std::list<SignedHit> hits;
  SignedHit getHit(int _i){
    std::list<SignedHit>::iterator it = hits.begin();
    for(int i=0; i<_i; i++){
      ++it;
    }
    return *it;
  }
  
    //Corresponding prop. tube segments
    PropSegment seg_x;
    PropSegment seg_y;

    //Slope, intersection, momentum and their errors
    double tx;
    double ty;
    double x0;
    double y0;
    double invP;

    //This is an admittedly messy way of keeping track of various bits of needed information to describe the possible particle trajectories in a single-station tracklet
    struct linedef {
      double slopeX;
      double slopeY;
      double initialX;
      double initialY;
      double initialZ;
      double slopeU;
      double initialU;
      double slopeV;
      double initialV;

      double wireHit1Pos;
      double wireHit2Pos;
      double wireHit1PosX;
      double wireHit2PosX;
      double wireHit1PosY;
      double wireHit2PosY;
      double wireHit1PosZ;
      double wireHit2PosZ;
      double wire1Slope;
      double wire2Slope;
      double wireIntercept1;
      double wireIntercept2;
      
      void print(){
	std::cout<<"slopeX: "<<slopeX<<" slopeY: "<<slopeY<<" initialX: "<<initialX<<" initialY: "<<initialY<<" initialZ: "<<initialZ<<" slopeU: "<<slopeU<<" initialU: "<<initialU<<" slopeV: "<<slopeV<<" initialV: "<<initialV<<std::endl;
      }

    } ;

  double st2X;
  double st3X;
  double st2Xsl;
  double st3Xsl;
  double st2U;
  double st3U;
  double st3V;
  double st2V;
  double st2Z;
  double st3Z;
  double st2Y;
  double st3Y;
  double st2Usl;
  double st2Vsl;
  double st3Usl;
  double st3Vsl;
  double st2UZ;
  double st2VZ;

  linedef acceptedXLine2;
  linedef acceptedULine2;
  linedef acceptedVLine2;
  linedef acceptedXLine3;
  linedef acceptedULine3;
  linedef acceptedVLine3;
  
    std::vector<linedef> possibleXLines;
    std::vector<linedef> possibleULines;
    std::vector<linedef> possibleVLines;
  
    void getSlopesX(Hit hit1, Hit hit2);
    void getSlopesU(Hit hit1, Hit hit2);
    void getSlopesV(Hit hit1, Hit hit2);

    double err_tx;
    double err_ty;
    double err_x0;
    double err_y0;
    double err_invP;

    //Residuals of all pos
    double residual[nChamberPlanes];

    bool _chargeSet = false;
    int _charge = 0;
    void setCharge(int chrg);

    ClassDef(Tracklet, 4)
};

class TrackletVector : public PHObject
{
public:
    TrackletVector();
    virtual ~TrackletVector();

    void identify(std::ostream& os = std::cout) const;
    void Reset();
    int  isValid() const { return 1; };
    TrackletVector* Clone() const { return (new TrackletVector(*this)); }

    bool empty() const { return trackletVec.empty(); }
    size_t size() const { return trackletVec.size(); }
    void clear() { Reset(); }

    const Tracklet* at(const size_t index) const;
    Tracklet* at(const size_t index);
    void push_back(const Tracklet* tracklet);
    size_t erase(const size_t index);

    std::vector<Tracklet*>::const_iterator begin() const { return trackletVec.begin(); }
    std::vector<Tracklet*>::const_iterator end()   const { return trackletVec.end(); }  

    std::vector<Tracklet*>::iterator begin() { return trackletVec.begin(); }
    std::vector<Tracklet*>::iterator end()   { return trackletVec.end(); }

private:
    std::vector<Tracklet*> trackletVec;

    ClassDef(TrackletVector, 1)
};


#endif
