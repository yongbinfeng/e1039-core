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

    // PHObject virtual overloads
    void identify(std::ostream &os = std::cout) const;
    void Reset()
    {
        hit = Hit();
        sign = 0;
    }
    int isValid() const { return hit.index > 0; }
    SignedHit *Clone() const { return (new SignedHit(hit, sign)); }

    // comparision operators for sorting
    bool operator<(const SignedHit elem) const { return hit.detectorID < elem.hit.detectorID; }
    bool operator==(const SignedHit elem) const { return hit.index == elem.hit.index; }

    // Get the real hit position
    double pos() { return hit.pos + sign * hit.driftDistance; }
    double pos(int sign_input) { return hit.pos + sign_input * hit.driftDistance; }

    // Data members
    Hit hit;
    int sign;

    ClassDef(SignedHit, 1)
};

class PropSegment : public PHObject
{
public:
    PropSegment();

    // PHObject virtual overloads
    void identify(std::ostream &os = std::cout) const { print(os); };
    void Reset() { init(); }
    int isValid() const;
    PropSegment *Clone() const { return (new PropSegment(*this)); }

    // init -- temporary, only used for tests
    void init();

    // Quality cut
    // bool isValid();

    // Debugging output
    void print(std::ostream &os = std::cout) const;

#ifndef __CINT__
    // Get expected position at a given Z
    double getExpPosition(double z) const { return a * z + b; }

    // Get the closest approach to a given space position/proptube wire
    double getClosestApproach(double z, double pos);

    // Get reference pos at first two planes
    double getPosRef(double default_val = -9999.);

    // Number of hits
    int getNHits() const;

    // Number of planes
    int getNPlanes() const;

    // Fit the segment -- naive linear fit
    void fit(); // external call
    void fit_2hits();
    void fit_34hits();

    // linear chisq fitter
    void linearFit_simple();
    void linearFit_iterative();

    // resolve left/right
    void resolveLR();
    void resolveLR(int setting);
#endif

    // track slope the interception
    double a;
    double b;
    double err_a;
    double err_b;

    // chisq of the segment
    double chisq;

    // Auxilary hodoscope hit list, cannot be possibly more than 4 hodoscope hits
    int nHodoHits;
    Hit hodoHits[4];

    // Hit list -- only 4 hits at most
    SignedHit hits[4];

    ClassDef(PropSegment, 4)
};

class Tracklet : public PHObject
{
public:
    Tracklet();

    // PHObject virtual overloads
    void identify(std::ostream &os = std::cout) const { os << "Tracklet @sID=" << stationID << std::endl; }
    void Reset() { *this = Tracklet(); }
    int isValid();
    Tracklet *Clone() const { return (new Tracklet(*this)); }

    // Basic quality cut
    // bool isValid();

    // Debuggin output
    void print(std::ostream &os = std::cout);

#ifndef __CINT__
    // Sort hit list
    void sortHits() { hits.sort(); }

    // Update/get number of real hits
    int getNHits() const { return nHits[0] + nHits[1] + nHits[2]; }

    // Number of all hits (even excluded)
    int getNAllHits() { return hits.size(); }

    // Get the probabilities
    double getProb() const;

    // Get the momentum probabilities
    double getMomProb() const;

    // Get the chi square
    double getChisq() const { return chisq; }

    // Get x and y positions at a given z
    double getExpPositionX(double z) const;
    double getExpPosErrorX(double z) const;
    double getExpPositionY(double z) const;
    double getExpPosErrorY(double z) const;
    double getExpPositionW(int detectorID) const;
    int getExpElementID(int detectorID) const;

    // Get momentum upstream/downstream
    TVector3 getMomentumSt1() const;
    TVector3 getMomentumSt3() const;
    TVector3 getExpMomentum(double z) const;

    // Get the i-th signed hit
    SignedHit getSignedHit(int index);

    // Kernal function to calculate chi square for minimizer
    double Eval(const double *par);
    double Eval4(const double *par);
    double calcChisq();
    double calcChisq_noDrift();
    double calcChisq_verbose();
    double calcChisq_st1_squares();
    double calcChisq_st1_sum();
    double calcChisq_st1_absSum();

    // Add dummy hits
    void addDummyHits();

    // Momentum estimation using back partial
    double getMomentum() const;

    // Decide charge by KMag bending direction
    // int getCharge() const { return x0*KMAGSTR > tx ? 1 : -1; }
    int getCharge() const;

    // Get the slope and intersection in station 1
    void getXZInfoInSt1(double &tx_st1, double &x0_st1) const;
    void getXZErrorInSt1(double &err_tx_st1, double &err_x0_st1) const;

    // For sorting tracklet list
    bool operator<(const Tracklet &elem) const;
    bool operator==(const Tracklet &elem) const;

    // For reducing similar tracklets
    bool similarity(const Tracklet &elem) const;
    double similarity_st1(const Tracklet &elem) const;
    bool similarityAllowed(const Tracklet &elem) const;
    bool elementSimilarity(const Tracklet &elem) const;

    // Merge the hit list from two tracklets
    Tracklet merge(Tracklet &elem);

    // For adding two tracklets together to form a back partial track
    Tracklet operator+(const Tracklet &elem);

    // For adding two tracklets together to form a global track
    Tracklet operator*(const Tracklet &elem);

    // Convert to a SRecTrack
    SRecTrack getSRecTrack(bool hyptest = true);
    SignedHit getHit(int _i);

    void getSlopes(Hit hit1, Hit hit2, int plane_type);
    void setCharge(int chrg) { _charge = chrg; } // This is needed for displaced tracking, where hte above formula does not work
    static std::pair<int, int> resolvePair(int a) { return std::make_pair(((4 - a) / 2) * 2 - 1, (a % 2) * 2 - 1); }
#endif

    // Station ID, ranging from 1 to nStation, nStation-1 means back partial track, nStation means global track
    int stationID;

    int nHits[3]; /// 0 for X plane, 1 for U plane, 2 for Z plane

    // Chi square
    double chisq;

    // chisq at vertex
    double chisq_vtx;

    // List of signed hits
    std::list<SignedHit> hits;

    // Corresponding prop. tube segments
    PropSegment seg_x;
    PropSegment seg_y;

    // Slope, intersection, momentum and their errors
    double tx;
    double ty;
    double x0;
    double y0;
    double invP;

    double err_tx;
    double err_ty;
    double err_x0;
    double err_y0;
    double err_invP;

    // Charge
    int _charge;

    // Residuals of all pos
    double residual[nChamberPlanes];

    std::vector<double> vtxHypos;

    /// XinL: Below is the information only used in building the "trackletslim",
    /// which only function between building the "actual tracklet". I don't know
    /// if it is a good idea to build a new class

    // for station 3 tracklets keep track of whether the particle is in D3m or D3p
    bool isM;

    // Normal of the plane spaned by tracklet and wire
    TVector3 planeNorm;

    // Temporary information for 4 sign combinations in a station
    struct linedef
    {
        double slope;
        double initial_pos;
        double initialZ;
    };

    std::vector<linedef> possibleLines;
    linedef acceptedLine2;
    linedef acceptedLine3;

    // Temporary st2 and st3 information
    double pos_st2;
    double z_st2;
    double slope_st2;

    double pos_st3;
    double z_st3;
    double slope_st3;

    // Hodo hits information used to match tracklets in st2 and st3
    std::vector<std::pair<int, int>> allowedHodos;
    std::vector<std::pair<int, int>> matching_combos;

    struct UXCombo
    {
        int trackletUIndex;
        std::vector<std::pair<int, int>> hodoMatches;
        double tx;
        double ty;
        double y_st2;
    };
    std::vector<UXCombo> allowedUXCombos;

    ClassDef(Tracklet, 4)
};

class TrackletVector : public PHObject
{
public:
    TrackletVector();
    virtual ~TrackletVector();

    void identify(std::ostream &os = std::cout) const;
    void Reset();
    int isValid() const { return 1; };
    TrackletVector *Clone() const { return (new TrackletVector(*this)); }

    bool empty() const { return trackletVec.empty(); }
    size_t size() const { return trackletVec.size(); }
    void clear() { Reset(); }

    const Tracklet *at(const size_t index) const;
    Tracklet *at(const size_t index);
    void push_back(const Tracklet *tracklet);
    size_t erase(const size_t index);

    std::vector<Tracklet *>::const_iterator begin() const { return trackletVec.begin(); }
    std::vector<Tracklet *>::const_iterator end() const { return trackletVec.end(); }

    std::vector<Tracklet *>::iterator begin() { return trackletVec.begin(); }
    std::vector<Tracklet *>::iterator end() { return trackletVec.end(); }

    std::vector<Tracklet *> trackletVec;

    ClassDef(TrackletVector, 1)
};

#endif
