/*
KalmanFastTracking.cxx

Implementation of class Tracklet, KalmanFastTracking

Author: Kun Liu, liuk@fnal.gov
Created: 05-28-2013
*/

#include <phfield/PHField.h>
#include <phool/PHTimer.h>
#include <phool/recoConsts.h>
#include <fun4all/Fun4AllBase.h>

#include <iostream>
#include <algorithm>
#include <cmath>

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TBox.h>
#include <TMatrixD.h>

#include "KalmanFastTracking.h"
#include "TriggerRoad.h"

#define _DEBUG_ON

namespace 
{
    //static flag to indicate the initialized has been done
    static bool inited = false;

    //Event acceptance cut
    static int MaxHitsDC0;
    static int MaxHitsDC1;
    static int MaxHitsDC2;
    static int MaxHitsDC3p;
    static int MaxHitsDC3m;

    //Sagitta ratio
    static double SAGITTA_DUMP_CENTER;
    static double SAGITTA_DUMP_WIDTH;
    static double SAGITTA_TARGET_CENTER;
    static double SAGITTA_TARGET_WIDTH;
    static double Z_TARGET;
    static double Z_DUMP;

    //Track quality cuts
    static double TX_MAX;
    static double TY_MAX;
    static double X0_MAX;
    static double Y0_MAX;
    static double INVP_MAX;
    static double INVP_MIN;
    static double Z_KMAG_BEND;

    //MuID cuts 
    static double MUID_REJECTION;
    static double MUID_Z_REF;
    static double MUID_R_CUT;
    static double MUID_THE_P0;
    static double MUID_EMP_P0;
    static double MUID_EMP_P1;
    static double MUID_EMP_P2;
    static int    MUID_MINHITS;

    //Track merging threshold
    static double MERGE_THRES;

    //static flag of kmag on/off
	static bool KMAG_ON;

    //running mode
    static bool MC_MODE;
    static bool COSMIC_MODE;
    static bool COARSE_MODE;

    //if displaced, skip fit to the target/vertex
    static bool NOT_DISPLACED;
  static bool TRACK_ELECTRONS; //please see comment in framework/phool/recoConsts.cc
  static bool TRACK_DISPLACED; //please see comment in framework/phool/recoConsts.cc
  static bool OLD_TRACKING; //please see comment in framework/phool/recoConsts.cc

    //initialize global variables
    void initGlobalVariables()
    {
        if(!inited) 
        {
            inited = true;

            recoConsts* rc = recoConsts::instance();
            MC_MODE = rc->get_BoolFlag("MC_MODE");
            KMAG_ON = rc->get_BoolFlag("KMAG_ON");
            COSMIC_MODE = rc->get_BoolFlag("COSMIC_MODE");
            COARSE_MODE = rc->get_BoolFlag("COARSE_MODE");

            NOT_DISPLACED = rc->get_BoolFlag("NOT_DISPLACED");
            TRACK_ELECTRONS = rc->get_BoolFlag("TRACK_ELECTRONS");
            TRACK_DISPLACED = rc->get_BoolFlag("TRACK_DISPLACED");
            OLD_TRACKING = rc->get_BoolFlag("OLD_TRACKING");

            MaxHitsDC0 = rc->get_IntFlag("MaxHitsDC0");
            MaxHitsDC1 = rc->get_IntFlag("MaxHitsDC1");
            MaxHitsDC2 = rc->get_IntFlag("MaxHitsDC2");
            MaxHitsDC3p = rc->get_IntFlag("MaxHitsDC3p");
            MaxHitsDC3m = rc->get_IntFlag("MaxHitsDC3m");
            
            TX_MAX = rc->get_DoubleFlag("TX_MAX");
            TY_MAX = rc->get_DoubleFlag("TY_MAX");
            X0_MAX = rc->get_DoubleFlag("X0_MAX");
            Y0_MAX = rc->get_DoubleFlag("Y0_MAX");
            INVP_MAX = rc->get_DoubleFlag("INVP_MAX");
            INVP_MIN = rc->get_DoubleFlag("INVP_MIN");
            Z_KMAG_BEND = rc->get_DoubleFlag("Z_KMAG_BEND");

            SAGITTA_TARGET_CENTER = rc->get_DoubleFlag("SAGITTA_TARGET_CENTER");
            SAGITTA_TARGET_WIDTH = rc->get_DoubleFlag("SAGITTA_TARGET_WIDTH");
            SAGITTA_DUMP_CENTER = rc->get_DoubleFlag("SAGITTA_DUMP_CENTER");
            SAGITTA_DUMP_WIDTH = rc->get_DoubleFlag("SAGITTA_DUMP_WIDTH");
            Z_TARGET = rc->get_DoubleFlag("Z_TARGET");
            Z_DUMP = rc->get_DoubleFlag("Z_DUMP");

            MUID_REJECTION = rc->get_DoubleFlag("MUID_REJECTION");
            MUID_Z_REF = rc->get_DoubleFlag("MUID_Z_REF");
            MUID_R_CUT = rc->get_DoubleFlag("MUID_R_CUT");
            MUID_THE_P0 = rc->get_DoubleFlag("MUID_THE_P0");
            MUID_EMP_P0 = rc->get_DoubleFlag("MUID_EMP_P0");
            MUID_EMP_P1 = rc->get_DoubleFlag("MUID_EMP_P1");
            MUID_EMP_P2 = rc->get_DoubleFlag("MUID_EMP_P2");
            MUID_MINHITS = rc->get_IntFlag("MUID_MINHITS");
        }
    }
}

KalmanFastTracking::KalmanFastTracking(const PHField* field, const TGeoManager* geom, bool flag): verbosity(0), enable_KF(flag), outputListIdx(4)
{
    using namespace std;
    initGlobalVariables();

#ifdef _DEBUG_ON
    cout << "Initialization of KalmanFastTracking ..." << endl;
    cout << "========================================" << endl;
#endif

    _timers.insert(std::make_pair<std::string, PHTimer*>("st2", new PHTimer("st2")));
    _timers.insert(std::make_pair<std::string, PHTimer*>("st3", new PHTimer("st3")));
    _timers.insert(std::make_pair<std::string, PHTimer*>("st23", new PHTimer("st23")));
    _timers.insert(std::make_pair<std::string, PHTimer*>("global", new PHTimer("global")));
    _timers.insert(std::make_pair<std::string, PHTimer*>("global_st1", new PHTimer("global_st1")));
    _timers.insert(std::make_pair<std::string, PHTimer*>("global_link", new PHTimer("global_link")));
    _timers.insert(std::make_pair<std::string, PHTimer*>("global_kalman", new PHTimer("global_kalman")));
    _timers.insert(std::make_pair<std::string, PHTimer*>("kalman", new PHTimer("kalman")));

    //Initialize Kalman fitter
    if(enable_KF)
    {
        kmfitter = new KalmanFitter(field, geom);
        kmfitter->setControlParameter(50, 0.001);
    }

    //Initialize minuit minimizer
    minimizer[0] = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
    minimizer[1] = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
    fcn = ROOT::Math::Functor(&tracklet_curr, &Tracklet::Eval, KMAG_ON ? 5 : 4);
    for(int i = 0; i < 2; ++i)
    {
        minimizer[i]->SetMaxFunctionCalls(1000000);
        minimizer[i]->SetMaxIterations(100);
        minimizer[i]->SetTolerance(1E-2);
        minimizer[i]->SetFunction(fcn);
        minimizer[i]->SetPrintLevel(0);
    }

    //Minimize ROOT output
    extern Int_t gErrorIgnoreLevel;
    gErrorIgnoreLevel = 9999;

    //Initialize geometry service
    p_geomSvc = GeomSvc::instance();
#ifdef _DEBUG_ON
    p_geomSvc->printTable();
    p_geomSvc->printWirePosition();
    p_geomSvc->printAlignPar();
#endif

    //Initialize plane angles for all planes
    for(int i = 1; i <= nChamberPlanes; ++i)
    {
        costheta_plane[i] = p_geomSvc->getCostheta(i);
        sintheta_plane[i] = p_geomSvc->getSintheta(i);
    }

    //Initialize hodoscope IDs
    detectorIDs_mask[0] = p_geomSvc->getDetectorIDs("H1");
    detectorIDs_mask[1] = p_geomSvc->getDetectorIDs("H2");
    detectorIDs_mask[2] = p_geomSvc->getDetectorIDs("H3");
    detectorIDs_mask[3] = p_geomSvc->getDetectorIDs("H4");
    detectorIDs_maskX[0] = p_geomSvc->getDetectorIDs("H1[TB]");
    detectorIDs_maskX[1] = p_geomSvc->getDetectorIDs("H2[TB]");
    detectorIDs_maskX[2] = p_geomSvc->getDetectorIDs("H3[TB]");
    detectorIDs_maskX[3] = p_geomSvc->getDetectorIDs("H4[TB]");
    detectorIDs_maskY[0] = p_geomSvc->getDetectorIDs("H1[LR]");
    detectorIDs_maskY[1] = p_geomSvc->getDetectorIDs("H2[LR]");
    detectorIDs_maskY[2] = p_geomSvc->getDetectorIDs("H4Y1[LR]");
    detectorIDs_maskY[3] = p_geomSvc->getDetectorIDs("H4Y2[LR]");
    detectorIDs_muidHodoAid[0] = p_geomSvc->getDetectorIDs("H4[TB]");
    detectorIDs_muidHodoAid[1] = p_geomSvc->getDetectorIDs("H4Y");

    //Register masking stations for tracklets in station-0/1, 2, 3+/-
    stationIDs_mask[0].push_back(1);
    stationIDs_mask[1].push_back(1);
    stationIDs_mask[2].push_back(2);
    stationIDs_mask[3].push_back(3);
    stationIDs_mask[4].push_back(3);

    //Masking stations for back partial
    stationIDs_mask[5].push_back(2);
    stationIDs_mask[5].push_back(3);
    stationIDs_mask[5].push_back(4);

    //Masking stations for global track
    stationIDs_mask[6].push_back(1);
    stationIDs_mask[6].push_back(2);
    stationIDs_mask[6].push_back(3);
    stationIDs_mask[6].push_back(4);

    //prop. tube IDs for mu id
    detectorIDs_muid[0][0] = p_geomSvc->getDetectorID("P1X1");
    detectorIDs_muid[0][1] = p_geomSvc->getDetectorID("P1X2");
    detectorIDs_muid[0][2] = p_geomSvc->getDetectorID("P2X1");
    detectorIDs_muid[0][3] = p_geomSvc->getDetectorID("P2X2");
    detectorIDs_muid[1][0] = p_geomSvc->getDetectorID("P1Y1");
    detectorIDs_muid[1][1] = p_geomSvc->getDetectorID("P1Y2");
    detectorIDs_muid[1][2] = p_geomSvc->getDetectorID("P2Y1");
    detectorIDs_muid[1][3] = p_geomSvc->getDetectorID("P2Y2");

    //Reference z_ref for mu id
    z_ref_muid[0][0] = MUID_Z_REF;
    z_ref_muid[0][1] = MUID_Z_REF;
    z_ref_muid[0][2] = 0.5*(p_geomSvc->getPlanePosition(detectorIDs_muid[0][0]) + p_geomSvc->getPlanePosition(detectorIDs_muid[0][1]));
    z_ref_muid[0][3] = z_ref_muid[0][2];

    z_ref_muid[1][0] = MUID_Z_REF;
    z_ref_muid[1][1] = MUID_Z_REF;
    z_ref_muid[1][2] = 0.5*(p_geomSvc->getPlanePosition(detectorIDs_muid[1][0]) + p_geomSvc->getPlanePosition(detectorIDs_muid[1][1]));
    z_ref_muid[1][3] = z_ref_muid[1][2];

    //Initialize masking window sizes, with optimized contingency
    for(int i = nChamberPlanes+1; i <= nChamberPlanes+nHodoPlanes+nPropPlanes; i++)
    {
        double factor = 0.;
        if(i > nChamberPlanes    && i <= nChamberPlanes+4)           factor = 0.25;    //for station-1 hodo
        if(i > nChamberPlanes+4  && i <= nChamberPlanes+8)           factor = 0.2;     //for station-2 hodo
        if(i > nChamberPlanes+8  && i <= nChamberPlanes+10)          factor = 0.15;    //for station-3 hodo
        if(i > nChamberPlanes+10 && i <= nChamberPlanes+nHodoPlanes) factor = 0.;      //for station-4 hodo
        if(i > nChamberPlanes+nHodoPlanes)                           factor = 0.15;    //for station-4 proptube

        z_mask[i-nChamberPlanes-1] = p_geomSvc->getPlanePosition(i);
        for(int j = 1; j <= p_geomSvc->getPlaneNElements(i); j++)
        {
            double x_min, x_max, y_min, y_max;
            p_geomSvc->get2DBoxSize(i, j, x_min, x_max, y_min, y_max);

            x_min -= (factor*(x_max - x_min));
            x_max += (factor*(x_max - x_min));
            y_min -= (factor*(y_max - y_min));
            y_max += (factor*(y_max - y_min));

            x_mask_min[i-nChamberPlanes-1][j-1] = x_min;
            x_mask_max[i-nChamberPlanes-1][j-1] = x_max;
            y_mask_min[i-nChamberPlanes-1][j-1] = y_min;
            y_mask_max[i-nChamberPlanes-1][j-1] = y_max;
        }
    }

#ifdef _DEBUG_ON
    cout << "========================" << endl;
    cout << "Hodo. masking settings: " << endl;
    for(int i = 0; i < 4; i++)
    {
        cout << "For station " << i+1 << endl;
        for(std::vector<int>::iterator iter = detectorIDs_mask[i].begin();  iter != detectorIDs_mask[i].end();  ++iter) cout << "All: " << *iter << endl;
        for(std::vector<int>::iterator iter = detectorIDs_maskX[i].begin(); iter != detectorIDs_maskX[i].end(); ++iter) cout << "X:   " << *iter << endl;
        for(std::vector<int>::iterator iter = detectorIDs_maskY[i].begin(); iter != detectorIDs_maskY[i].end(); ++iter) cout << "Y:   " << *iter << endl;
    }

    for(int i = 0; i < nStations; ++i)
    {
        std::cout << "Masking stations for tracklets with stationID = " << i + 1 << ": " << std::endl;
        for(std::vector<int>::iterator iter = stationIDs_mask[i].begin(); iter != stationIDs_mask[i].end(); ++iter)
        {
            std::cout << *iter << "  ";
        }
        std::cout << std::endl;
    }
#endif

    //Initialize super stationIDs
    for(int i = 0; i < nChamberPlanes/6+2; i++) superIDs[i].clear();
    superIDs[0].push_back((p_geomSvc->getDetectorIDs("D0X")[0] + 1)/2);
    superIDs[0].push_back((p_geomSvc->getDetectorIDs("D0U")[0] + 1)/2);
    superIDs[0].push_back((p_geomSvc->getDetectorIDs("D0V")[0] + 1)/2);
    superIDs[1].push_back((p_geomSvc->getDetectorIDs("D1X")[0] + 1)/2);
    superIDs[1].push_back((p_geomSvc->getDetectorIDs("D1U")[0] + 1)/2);
    superIDs[1].push_back((p_geomSvc->getDetectorIDs("D1V")[0] + 1)/2);
    superIDs[2].push_back((p_geomSvc->getDetectorIDs("D2X")[0] + 1)/2);
    superIDs[2].push_back((p_geomSvc->getDetectorIDs("D2U")[0] + 1)/2);
    superIDs[2].push_back((p_geomSvc->getDetectorIDs("D2V")[0] + 1)/2);
    superIDs[3].push_back((p_geomSvc->getDetectorIDs("D3pX")[0] + 1)/2);
    superIDs[3].push_back((p_geomSvc->getDetectorIDs("D3pU")[0] + 1)/2);
    superIDs[3].push_back((p_geomSvc->getDetectorIDs("D3pV")[0] + 1)/2);
    superIDs[4].push_back((p_geomSvc->getDetectorIDs("D3mX")[0] + 1)/2);
    superIDs[4].push_back((p_geomSvc->getDetectorIDs("D3mU")[0] + 1)/2);
    superIDs[4].push_back((p_geomSvc->getDetectorIDs("D3mV")[0] + 1)/2);

    superIDs[5].push_back((p_geomSvc->getDetectorIDs("P1X")[0] + 1)/2);
    superIDs[5].push_back((p_geomSvc->getDetectorIDs("P2X")[0] + 1)/2);
    superIDs[6].push_back((p_geomSvc->getDetectorIDs("P1Y")[0] + 1)/2);
    superIDs[6].push_back((p_geomSvc->getDetectorIDs("P2Y")[0] + 1)/2);

#ifdef _DEBUG_ON
    cout << "=============" << endl;
    cout << "Chamber IDs: " << endl;
    TString stereoNames[3] = {"X", "U", "V"};
    for(int i = 0; i < nChamberPlanes/6; i++)
    {
        for(int j = 0; j < 3; j++) cout << i << "  " << stereoNames[j].Data() << ": " << superIDs[i][j] << endl;
    }

    cout << "Proptube IDs: " << endl;
    for(int i = nChamberPlanes/6; i < nChamberPlanes/6+2; i++)
    {
        for(int j = 0; j < 2; j++) cout << i << "  " << j << ": " << superIDs[i][j] << endl;
    }

    //Initialize widow sizes for X-U matching and z positions of all chambers
    cout << "======================" << endl;
    cout << "U plane window sizes: " << endl;
#endif

    double u_factor[] = {5., 5., 5., 15., 15.};
    for(int i = 0; i < nChamberPlanes/6; i++)
    {
        int xID = 2*superIDs[i][0] - 1;
        int uID = 2*superIDs[i][1] - 1;
        int vID = 2*superIDs[i][2] - 1;
        double spacing = p_geomSvc->getPlaneSpacing(uID);
        double x_span = p_geomSvc->getPlaneScaleY(uID);

        z_plane_x[i] = 0.5*(p_geomSvc->getPlanePosition(xID) + p_geomSvc->getPlanePosition(xID+1));
        z_plane_u[i] = 0.5*(p_geomSvc->getPlanePosition(uID) + p_geomSvc->getPlanePosition(uID+1));
        z_plane_v[i] = 0.5*(p_geomSvc->getPlanePosition(vID) + p_geomSvc->getPlanePosition(vID+1));

        u_costheta[i] = costheta_plane[uID];
        u_sintheta[i] = sintheta_plane[uID];

        //u_win[i] = fabs(0.5*x_span/(spacing/sintheta_plane[uID])) + 2.*spacing + u_factor[i];
        u_win[i] = fabs(0.5*x_span*sintheta_plane[uID]) + TX_MAX*fabs((z_plane_u[i] - z_plane_x[i])*u_costheta[i]) + TY_MAX*fabs((z_plane_u[i] - z_plane_x[i])*u_sintheta[i]) + 2.*spacing + u_factor[i];

#ifdef _DEBUG_ON
        cout << "Station " << i << ": " << xID << "  " << uID << "  " << vID << "  " << u_win[i] << endl;
#endif
    }

    //Initialize Z positions and maximum parameters of all planes
    for(int i = 1; i <= nChamberPlanes; i++)
    {
        z_plane[i] = p_geomSvc->getPlanePosition(i);
        slope_max[i] = costheta_plane[i]*TX_MAX + sintheta_plane[i]*TY_MAX;
        intersection_max[i] = costheta_plane[i]*X0_MAX + sintheta_plane[i]*Y0_MAX;

#ifdef COARSE_MODE
        resol_plane[i] = 3.*p_geomSvc->getPlaneSpacing(i)/sqrt(12.);
#else
        if(i <= 6)
        {
            resol_plane[i] = recoConsts::instance()->get_DoubleFlag("RejectWinDC0");
        }
        else if(i <= 12)
        {
            resol_plane[i] = recoConsts::instance()->get_DoubleFlag("RejectWinDC1");
        }
        else if(i <= 18)
        {
            resol_plane[i] = recoConsts::instance()->get_DoubleFlag("RejectWinDC2");
        }
        else if(i <= 24)
        {
            resol_plane[i] = recoConsts::instance()->get_DoubleFlag("RejectWinDC3p");
        }
        else
        {
            resol_plane[i] = recoConsts::instance()->get_DoubleFlag("RejectWinDC3m");
        }
#endif
        spacing_plane[i] = p_geomSvc->getPlaneSpacing(i);
    }

#ifdef _DEBUG_ON
    cout << "======================================" << endl;
    cout << "Maximum local slope and intersection: " << endl;
#endif
    for(int i = 1; i <= nChamberPlanes/2; i++)
    {
        double d_slope = (p_geomSvc->getPlaneResolution(2*i - 1) + p_geomSvc->getPlaneResolution(2*i))/(z_plane[2*i] - z_plane[2*i-1]);
        double d_intersection = d_slope*z_plane[2*i];

        slope_max[2*i-1] += d_slope;
        intersection_max[2*i-1] += d_intersection;
        slope_max[2*i] += d_slope;
        intersection_max[2*i] += d_intersection;

#ifdef _DEBUG_ON
        cout << "Super plane " << i << ": " << slope_max[2*i-1] << "  " << intersection_max[2*i-1] << endl;
#endif
    }

    //Initialize sagitta ratios, index 0, 1, 2 are for X, U, V, this is the incrementing order of plane type
    s_detectorID[0] = p_geomSvc->getDetectorID("D2X");
    s_detectorID[1] = p_geomSvc->getDetectorID("D2Up");
    s_detectorID[2] = p_geomSvc->getDetectorID("D2Vp");
}

KalmanFastTracking::~KalmanFastTracking()
{
    if(enable_KF) delete kmfitter;
    delete minimizer[0];
    delete minimizer[1];
}

void KalmanFastTracking::setRawEventDebug(SRawEvent* event_input)
{
    rawEvent = event_input;
    hitAll = event_input->getAllHits();
}

int KalmanFastTracking::setRawEvent(SRawEvent* event_input)
{
	//reset timer
    for(auto iter=_timers.begin(); iter != _timers.end(); ++iter) 
    {
        iter->second->reset();
    }

    //Initialize tracklet lists
    for(int i = 0; i < 5; i++){
      trackletsInSt[i].clear();
      trackletsInStSlimX[i].clear();
      trackletsInStSlimU[i].clear();
      trackletsInStSlimV[i].clear();
    }
    stracks.clear();

    //pre-tracking cuts
    rawEvent = event_input;
    if(!acceptEvent(rawEvent)) return TFEXIT_FAIL_MULTIPLICITY;
    hitAll = event_input->getAllHits();
#ifdef _DEBUG_ON
    for(std::vector<Hit>::iterator iter = hitAll.begin(); iter != hitAll.end(); ++iter) iter->print();
#endif

    //Initialize hodo and masking IDs
    for(int i = 0; i < 4; i++)
    {
        //std::cout << "For station " << i << std::endl;
        hitIDs_mask[i].clear();
        if(MC_MODE || COSMIC_MODE || rawEvent->isFPGATriggered())
        {
            hitIDs_mask[i] = rawEvent->getHitsIndexInDetectors(detectorIDs_maskX[i]);
        }
        else
        {
            hitIDs_mask[i] = rawEvent->getHitsIndexInDetectors(detectorIDs_maskY[i]);
        }

        //for(std::list<int>::iterator iter = hitIDs_mask[i].begin(); iter != hitIDs_mask[i].end(); ++iter) std::cout << *iter << " " << hitAll[*iter].detectorID << " === ";
        //std::cout << std::endl;
    }

    //Initialize prop. tube IDs
    for(int i = 0; i < 2; ++i)
    {
        for(int j = 0; j < 4; ++j)
        {
            hitIDs_muid[i][j].clear();
            hitIDs_muid[i][j] = rawEvent->getHitsIndexInDetector(detectorIDs_muid[i][j]);
        }
        hitIDs_muidHodoAid[i].clear();
        hitIDs_muidHodoAid[i] = rawEvent->getHitsIndexInDetectors(detectorIDs_muidHodoAid[i]);
    }

    if(!COARSE_MODE)
    {
        buildPropSegments();
        if(propSegs[0].empty() || propSegs[1].empty())
        {
#ifdef _DEBUG_ON
            LogInfo("Failed in prop tube segment building: " << propSegs[0].size() << ", " << propSegs[1].size());
#endif
            //return TFEXIT_FAIL_ROUGH_MUONID;
        }
    }

    //Build tracklets in station 2, 3+, 3-
    _timers["st2"]->restart();
    //buildTrackletsInStation(3, 1);   //3 for station-2, 1 for list position 1
    buildTrackletsInStationSlim(3, 1);   //3 for station-2, 1 for list position 1
    buildTrackletsInStationSlimU(3, 1);   //3 for station-2, 1 for list position 1
    buildTrackletsInStationSlimV(3, 1);   //3 for station-2, 1 for list position 1 
    /*if(!TRACK_DISPLACED){
      buildTrackletsInStation(3, 1);   //3 for station-2, 1 for list position 1
    } else{
      for(int pxs = -99; pxs < 100; pxs = pxs+2){ //WPM
	double pos_exp[3], window[3];
	pos_exp[0] = pxs;
	pos_exp[1] = p_geomSvc->getCostheta(1) * pos_exp[0] + p_geomSvc->getSintheta(1) * 0;
	pos_exp[2] = p_geomSvc->getCostheta(5) * pos_exp[0] + p_geomSvc->getSintheta(5) * 0;
	std::cout<<"The station 2 window I'm using is centered at "<<pxs<<", "<<p_geomSvc->getCostheta(1) * pos_exp[0] + p_geomSvc->getSintheta(1) * 0<<", "<<p_geomSvc->getCostheta(5) * pos_exp[0] + p_geomSvc->getSintheta(5) * 0<<std::endl; //WPM
	window[0] = 1;
	window[1] = 10;
	window[2] = 10;
	buildTrackletsInStation(3, 1, pos_exp, window);
      }
      }*/
    _timers["st2"]->stop();
    if(verbosity >= 2) LogInfo("NTracklets in St2: " << trackletsInSt[1].size());
    
    if(outputListIdx > 1 && trackletsInSt[1].empty())
    {
#ifdef _DEBUG_ON
        LogInfo("Failed in tracklet build at station 2");
#endif
        //return TFEXIT_FAIL_ST2_TRACKLET;
    }

    _timers["st3"]->restart();
    //buildTrackletsInStation(4, 2);   //4 for station-3+
    //buildTrackletsInStation(5, 2);   //5 for station-3-
    buildTrackletsInStationSlim(4, 2);   //4 for station-3+
    buildTrackletsInStationSlim(5, 2);   //5 for station-3-
    buildTrackletsInStationSlimU(4, 2);   //4 for station-3+
    buildTrackletsInStationSlimU(5, 2);   //5 for station-3-
    buildTrackletsInStationSlimV(4, 2);   //4 for station-3+
    buildTrackletsInStationSlimV(5, 2);   //5 for station-3-
    _timers["st3"]->stop();
    if(verbosity >= 2) LogInfo("NTracklets in St3: " << trackletsInSt[2].size());
    
    if(outputListIdx > 2 && trackletsInSt[2].empty())
    {
#ifdef _DEBUG_ON
        LogInfo("Failed in tracklet build at station 3");
#endif
        //return TFEXIT_FAIL_ST3_TRACKLET;
    }

    //Build back partial tracks in station 2, 3+ and 3-
    _timers["st23"]->restart();
    //buildBackPartialTracks();
    buildBackPartialTracksSlimX();
    buildBackPartialTracksSlimU();
    buildBackPartialTracksSlimV();
    buildBackPartialTracksSlim();
    //for(std::list<Tracklet>::iterator tracklet23 = trackletsInStSlimX[3].begin(); tracklet23 != trackletsInStSlimX[3].end(); ++tracklet23){
    //  buildTrackletsInStationWithUV(6, 3, *tracklet23);
    //}
    _timers["st23"]->stop();
    if(verbosity >= 2) LogInfo("NTracklets St2+St3: " << trackletsInSt[3].size());


    if(outputListIdx > 3 && trackletsInSt[3].empty())
    {
#ifdef _DEBUG_ON
        LogInfo("Failed in connecting station-2 and 3 tracks");
#endif
        return TFEXIT_FAIL_BACKPARTIAL;
    }

    //Connect tracklets in station 2/3 and station 1 to form global tracks
    _timers["global"]->restart();
    if(!TRACK_DISPLACED){
      buildGlobalTracks();
    } else{
      buildGlobalTracksDisplaced();
    }
    _timers["global"]->stop();
    if(verbosity >= 2) LogInfo("NTracklets Global: " << trackletsInSt[4].size());
    
#ifdef _DEBUG_ON
    for(int i = 0; i < 2; ++i)
    {
        std::cout << "=======================================================================================" << std::endl;
        LogInfo("Prop tube segments in " << (i == 0 ? "X-Z" : "Y-Z"));
        for(std::list<PropSegment>::iterator seg = propSegs[i].begin(); seg != propSegs[i].end(); ++seg)
        {
            seg->print();
        }
        std::cout << "=======================================================================================" << std::endl;
    }

    for(int i = 0; i <= 4; i++)
    {
        std::cout << "=======================================================================================" << std::endl;
        LogInfo("Final tracklets in station: " << i+1 << " is " << trackletsInSt[i].size());
        for(std::list<Tracklet>::iterator tracklet = trackletsInSt[i].begin(); tracklet != trackletsInSt[i].end(); ++tracklet)
        {
            tracklet->print();
        }
        std::cout << "=======================================================================================" << std::endl;
    }
#endif

    if(outputListIdx == 4 && trackletsInSt[4].empty()) return TFEXIT_FAIL_GLOABL;
    if(!enable_KF) return TFEXIT_SUCCESS;

    //Build kalman tracks
    _timers["kalman"]->restart();
    for(std::list<Tracklet>::iterator tracklet = trackletsInSt[outputListIdx].begin(); tracklet != trackletsInSt[outputListIdx].end(); ++tracklet)
    {
        SRecTrack strack = processOneTracklet(*tracklet);
        stracks.push_back(strack);
    }
    _timers["kalman"]->stop();

#ifdef _DEBUG_ON
    LogInfo(stracks.size() << " final tracks:");
    for(std::list<SRecTrack>::iterator strack = stracks.begin(); strack != stracks.end(); ++strack)
    {
        strack->print();
    }
#endif

    return TFEXIT_SUCCESS;
}

bool KalmanFastTracking::acceptEvent(SRawEvent* rawEvent)
{
    if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT) 
    {
        LogInfo("D0: " << rawEvent->getNHitsInD0());
        LogInfo("D1: " << rawEvent->getNHitsInD1());
        LogInfo("D2: " << rawEvent->getNHitsInD2());
        LogInfo("D3p: " << rawEvent->getNHitsInD3p());
        LogInfo("D3m: " << rawEvent->getNHitsInD3m());
        LogInfo("H1: " << rawEvent->getNHitsInDetectors(detectorIDs_maskX[0]));
        LogInfo("H2: " << rawEvent->getNHitsInDetectors(detectorIDs_maskX[1]));
        LogInfo("H3: " << rawEvent->getNHitsInDetectors(detectorIDs_maskX[2]));
        LogInfo("H4: " << rawEvent->getNHitsInDetectors(detectorIDs_maskX[3]));
        LogInfo("Prop:" << rawEvent->getNPropHitsAll());
        LogInfo("NRoadsPos: " << rawEvent->getNRoadsPos());
        LogInfo("NRoadsNeg: " << rawEvent->getNRoadsNeg());
    }

    if(rawEvent->getNHitsInD0() > MaxHitsDC0) return false;
    if(rawEvent->getNHitsInD1() > MaxHitsDC1) return false; 
    if(rawEvent->getNHitsInD2() > MaxHitsDC2) return false; 
    if(rawEvent->getNHitsInD3p() > MaxHitsDC3p) return false;
    if(rawEvent->getNHitsInD3m() > MaxHitsDC3m) return false;

    /*
    if(rawEvent->getNHitsInDetectors(detectorIDs_maskX[0]) > 15) return false;
    if(rawEvent->getNHitsInDetectors(detectorIDs_maskX[1]) > 10) return false;
    if(rawEvent->getNHitsInDetectors(detectorIDs_maskX[2]) > 10) return false;
    if(rawEvent->getNHitsInDetectors(detectorIDs_maskX[3]) > 10) return false;
    if(rawEvent->getNPropHitsAll() > 300) return false;

    if(rawEvent->getNRoadsPos() > 5) return false;
    if(rawEvent->getNRoadsNeg() > 5) return false;
    */

    return true;
}

void KalmanFastTracking::buildBackPartialTracks()
{
    //Temporary container for a simple chisq fit
    int nHitsX2, nHitsX3;
    double z_fit[4], x_fit[4];
    double a, b;

    for(std::list<Tracklet>::iterator tracklet3 = trackletsInSt[2].begin(); tracklet3 != trackletsInSt[2].end(); ++tracklet3)
    {
      //tracklet3->print();
        if(!COARSE_MODE)
        {
            //Extract the X hits only from station-3 tracks
            nHitsX3 = 0;
            for(std::list<SignedHit>::iterator ptr_hit = tracklet3->hits.begin(); ptr_hit != tracklet3->hits.end(); ++ptr_hit)
            {
	      if(ptr_hit->hit.index < 0) continue;
	      if(p_geomSvc->getPlaneType(ptr_hit->hit.detectorID) == 1)
                {
		  z_fit[nHitsX3] = z_plane[ptr_hit->hit.detectorID];
		  x_fit[nHitsX3] = ptr_hit->hit.pos;
		  ++nHitsX3;
                }
	    }
        }

        Tracklet tracklet_best;
        for(std::list<Tracklet>::iterator tracklet2 = trackletsInSt[1].begin(); tracklet2 != trackletsInSt[1].end(); ++tracklet2)
        {
	  //tracklet2->print();
            if(!COARSE_MODE)
            {
	      if(OLD_TRACKING){
		if(fabs(tracklet2->tx - tracklet3->tx) > 0.15 || fabs(tracklet2->ty - tracklet3->ty) > 0.1) continue;
	      }
                //Extract the X hits from station-2 tracke
                nHitsX2 = nHitsX3;
                for(std::list<SignedHit>::iterator ptr_hit = tracklet2->hits.begin(); ptr_hit != tracklet2->hits.end(); ++ptr_hit)
                {
                    if(ptr_hit->hit.index < 0) continue;
                    if(p_geomSvc->getPlaneType(ptr_hit->hit.detectorID) == 1)
                    {
                        z_fit[nHitsX2] = z_plane[ptr_hit->hit.detectorID];
                        x_fit[nHitsX2] = ptr_hit->hit.pos;
                        ++nHitsX2;
                    }
                }

                //Apply a simple linear fit to get rough estimation of X-Z slope and intersection
                chi2fit(nHitsX2, z_fit, x_fit, a, b);
                if(fabs(a) > 2.*TX_MAX || fabs(b) > 2.*X0_MAX) continue;

                //Project to proportional tubes to see if there is enough
                int nPropHits = 0;
                for(int i = 0; i < 4; ++i)
                {
                    double x_exp = a*z_mask[detectorIDs_muid[0][i] - nChamberPlanes - 1] + b;
                    for(std::list<int>::iterator iter = hitIDs_muid[0][i].begin(); iter != hitIDs_muid[0][i].end(); ++iter)
                    {
                        if(fabs(hitAll[*iter].pos - x_exp) < 5.08)
                        {
                            ++nPropHits;
                            break;
                        }
                    }
                    if(nPropHits > 0) break;
                }
                if(!TRACK_ELECTRONS && nPropHits == 0) continue; //Turned off by Patrick for electron tracks
            }

	    Tracklet tracklet_23;
	    if(OLD_TRACKING){
	      tracklet_23 = (*tracklet2) + (*tracklet3);
	    }
	    else{
	      if(compareTracklets(*tracklet2, *tracklet3)){
		tracklet_23 = (*tracklet2) + (*tracklet3);
		tracklet_23.tx = tracklet2->tx; //This is needed to "seed" the tracklet fit that happens below.  This tx and ty information is assigned in compareTracklets below
		tracklet_23.ty = tracklet2->ty;

		tracklet_23.st2Z = tracklet2->st2Z;
		tracklet_23.st2X = tracklet2->st2X;
		tracklet_23.st2Y = tracklet2->st2Y;
		tracklet_23.st3Y = tracklet2->st3Y;
		tracklet_23.st2U = tracklet2->st2U;
		tracklet_23.st2V = tracklet2->st2V;
		tracklet_23.st2Usl = tracklet2->st2Usl;
		tracklet_23.st2Vsl = tracklet2->st2Vsl;
	      }
	      else{
		continue;
	      }
	    }
#ifdef _DEBUG_ON
            LogInfo("Using following two tracklets:");
            tracklet2->print();
            tracklet3->print();
            LogInfo("Yield this combination:");
            tracklet_23.print();
#endif
            fitTracklet(tracklet_23); //This is the fit that needs the seeded tx and ty information. Without the seed information, the fit occasionally finds bad slope and X0 or Y0 values, much in the same way that it does for single-station tracklets.  Note from Patrick: this fit could probably be throw out, as we already know the tx and ty information from compareTracklets.  I would just need to extrapolate back to z = 0 and calculate the chisq by hand
            if(tracklet_23.chisq > 9000.)
            {
#ifdef _DEBUG_ON
                tracklet_23.print();
                LogInfo("Impossible combination!");
#endif
                continue;
            }

            if(!COARSE_MODE && !hodoMask(tracklet_23))
            {
#ifdef _DEBUG_ON
                LogInfo("Hodomasking failed!");
#endif
                continue;
            }
#ifdef _DEBUG_ON
            LogInfo("Hodomasking Scucess!");
#endif

            if(!COARSE_MODE)
            {
	      resolveLeftRight(tracklet_23, 40.); //resolveLeftRight at this stage is, in truth, not needed.  We already know which side of the wire the particle passed by from compareTracklet.  However, getting rid of this would take a bit of work
                resolveLeftRight(tracklet_23, 150.);
            }

            ///Remove bad hits if needed
            removeBadHits(tracklet_23);

#ifdef _DEBUG_ON
            LogInfo("New tracklet: ");
            tracklet_23.print();

            LogInfo("Current best:");
            tracklet_best.print();

            LogInfo("Comparison: " << (tracklet_23 < tracklet_best));
            LogInfo("Quality: " << acceptTracklet(tracklet_23));
#endif

            //If current tracklet is better than the best tracklet up-to-now
            if(acceptTracklet(tracklet_23) && tracklet_23 < tracklet_best)
            {
                tracklet_best = tracklet_23;
            }
#ifdef _DEBUG_ON
            else
            {
                LogInfo("Rejected!!");
            }
#endif
        }

        if(tracklet_best.isValid() > 0) trackletsInSt[3].push_back(tracklet_best);
    }

    reduceTrackletList(trackletsInSt[3]);
    trackletsInSt[3].sort();
}


void KalmanFastTracking::buildBackPartialTracksSlim()
{
#ifdef _DEBUG_ON
  LogInfo("HELLO I'm in buildBackPartialTracksSlim");
#endif
  std::cout<<"trackletsInStSlimV[3].size() = "<<trackletsInStSlimV[3].size()<<std::endl;
  for(std::list<Tracklet>::iterator trackletV = trackletsInStSlimV[3].begin(); trackletV != trackletsInStSlimV[3].end(); ++trackletV){
    std::cout<<"FIFTH TEST trackletV->acceptedVLine2.wireHit1PosZ = "<<(*trackletV).acceptedVLine2.wireHit1PosZ<<std::endl;
    (*trackletV).acceptedVLine2.print();
  }
  
  for(std::list<Tracklet>::iterator trackletX = trackletsInStSlimX[3].begin(); trackletX != trackletsInStSlimX[3].end(); ++trackletX){
    Tracklet tracklet_best;
    for(std::list<Tracklet>::iterator trackletU = trackletsInStSlimU[3].begin(); trackletU != trackletsInStSlimU[3].end(); ++trackletU){
      Tracklet tracklet_best_UX;
      for(std::list<Tracklet>::iterator trackletV = trackletsInStSlimV[3].begin(); trackletV != trackletsInStSlimV[3].end(); ++trackletV){

	double posy2U = ((*trackletU).st2U - p_geomSvc->getCostheta((*trackletX).getHit(0).hit.detectorID+2)*((*trackletX).st2X + (*trackletX).st2Xsl * ((*trackletU).st2Z - (*trackletX).st2Z) ))/p_geomSvc->getSintheta((*trackletX).getHit(0).hit.detectorID+2);
	double posy2V = ((*trackletV).st2V - p_geomSvc->getCostheta((*trackletX).getHit(0).hit.detectorID-2)*((*trackletX).st2X + (*trackletX).st2Xsl * ((*trackletV).st2Z - (*trackletX).st2Z) ))/p_geomSvc->getSintheta((*trackletX).getHit(0).hit.detectorID-2);
	double posy3U = ((*trackletU).st3U - p_geomSvc->getCostheta((*trackletX).getHit(2).hit.detectorID+2)*((*trackletX).st3X + (*trackletX).st3Xsl * ((*trackletU).st3Z - (*trackletX).st3Z) ))/p_geomSvc->getSintheta((*trackletX).getHit(2).hit.detectorID+2);
	double posy3V = ((*trackletV).st3V - p_geomSvc->getCostheta((*trackletX).getHit(2).hit.detectorID-2)*((*trackletX).st3X + (*trackletX).st3Xsl * ((*trackletV).st3Z - (*trackletX).st3Z) ))/p_geomSvc->getSintheta((*trackletX).getHit(2).hit.detectorID-2);

	if( std::abs(posy2U) > 100 || std::abs(posy2V) > 100 || std::abs(posy3U) > 100 || std::abs(posy3V) > 100 ) continue;
	std::cout<<"posy2U="<<posy2U<<"; posy2V="<<posy2V<<"; posy3U="<<posy3U<<"; posy3V="<<posy3V<<std::endl;

	
	int LR1 = 0;
	int LR2 = 0;
	Tracklet tracklet_st2;
	tracklet_st2.stationID = 3;
	
	//resolveLeftRight(*xiter, LR1, LR2);
	tracklet_st2.hits.push_back(SignedHit((*trackletX).getHit(0).hit, LR1));
	tracklet_st2.nXHits++;
	tracklet_st2.hits.push_back(SignedHit((*trackletX).getHit(1).hit, LR1));
	tracklet_st2.nXHits++;
	tracklet_st2.getSlopesX((*trackletX).getHit(0).hit, (*trackletX).getHit(1).hit);
	
	tracklet_st2.hits.push_back(SignedHit((*trackletU).getHit(0).hit, LR1));
	tracklet_st2.nUHits++;
	tracklet_st2.hits.push_back(SignedHit((*trackletU).getHit(1).hit, LR1));
	tracklet_st2.nUHits++;
	tracklet_st2.getSlopesU((*trackletU).getHit(0).hit, (*trackletU).getHit(1).hit);
	
	tracklet_st2.hits.push_back(SignedHit((*trackletV).getHit(0).hit, LR1));
	tracklet_st2.nVHits++;
	tracklet_st2.hits.push_back(SignedHit((*trackletV).getHit(1).hit, LR1));
	tracklet_st2.nVHits++;
	tracklet_st2.getSlopesV((*trackletV).getHit(0).hit, (*trackletV).getHit(1).hit);
	
	tracklet_st2.sortHits();
	if(tracklet_st2.isValid() == 0) //TODO: What IS THIS?
	  {
	    fitTracklet(tracklet_st2); //This is where the original DCA minimization is performed
	  }
	else
	  {
	    continue;
	  }
	
#ifdef _DEBUG_ON
	tracklet_st2.print();
#endif
	if(acceptTracklet(tracklet_st2))
	  {
	    trackletsInSt[1].push_back(tracklet_st2);
	  }
#ifdef _DEBUG_ON
	else
	  {
	    LogInfo("Rejected!!!");
	  }
#endif

	
	Tracklet tracklet_st3;
	if((*trackletX).getHit(2).hit.detectorID > 18 && (*trackletX).getHit(2).hit.detectorID < 25){
	  tracklet_st3.stationID = 4;
	}
	
	if((*trackletX).getHit(2).hit.detectorID > 24 && (*trackletX).getHit(2).hit.detectorID < 31){
	  tracklet_st3.stationID = 5;
	}
	
	//resolveLeftRight(*xiter, LR1, LR2);
	tracklet_st3.hits.push_back(SignedHit((*trackletX).getHit(2).hit, LR1));
	tracklet_st3.nXHits++;
	tracklet_st3.hits.push_back(SignedHit((*trackletX).getHit(3).hit, LR1));
	tracklet_st3.nXHits++;
	tracklet_st3.getSlopesX((*trackletX).getHit(2).hit, (*trackletX).getHit(3).hit);
	
	tracklet_st3.hits.push_back(SignedHit((*trackletU).getHit(2).hit, LR1));
	tracklet_st3.nUHits++;
	tracklet_st3.hits.push_back(SignedHit((*trackletU).getHit(3).hit, LR1));
	tracklet_st3.nUHits++;
	tracklet_st3.getSlopesU((*trackletU).getHit(2).hit, (*trackletU).getHit(3).hit);
	
	tracklet_st3.hits.push_back(SignedHit((*trackletV).getHit(2).hit, LR1));
	tracklet_st3.nVHits++;
	tracklet_st3.hits.push_back(SignedHit((*trackletV).getHit(3).hit, LR1));
	tracklet_st3.nVHits++;
	tracklet_st3.getSlopesV((*trackletV).getHit(2).hit, (*trackletV).getHit(3).hit);

	tracklet_st3.sortHits();
	if(tracklet_st3.isValid() == 0) //TODO: What IS THIS?
	  {
	    fitTracklet(tracklet_st3); //This is where the original DCA minimization is performed
	  }
	else
	  {
	    continue;
	  }
	
#ifdef _DEBUG_ON
	tracklet_st3.print();
#endif
	if(acceptTracklet(tracklet_st3))
	  {
	    trackletsInSt[2].push_back(tracklet_st3);
	  }
#ifdef _DEBUG_ON
	else
	  {
	    LogInfo("Rejected!!!");
	  }
#endif

	
	Tracklet tracklet_23;
	//Tracklet tracklet_23_reverse;
	if(compareTracklets(tracklet_st2, tracklet_st3)){
	  tracklet_23 = (tracklet_st2) + (tracklet_st3);
	  tracklet_23.tx = tracklet_st2.tx; //This is needed to "seed" the tracklet fit that happens below.  This tx and ty information is assigned in compareTracklets below
	  tracklet_23.ty = tracklet_st2.ty;
	  
	  tracklet_23.st2Z = tracklet_st2.st2Z;
	  tracklet_23.st2X = tracklet_st2.st2X;
	  tracklet_23.st2Y = tracklet_st2.st2Y;
	  tracklet_23.st3Y = tracklet_st2.st3Y;
	  tracklet_23.st2U = tracklet_st2.st2U;
	  tracklet_23.st2V = tracklet_st2.st2V;
	  tracklet_23.st2Usl = tracklet_st2.st2Usl;
	  tracklet_23.st2Vsl = tracklet_st2.st2Vsl;

	  //tracklet_23_reverse = tracklet_23;
	  //tracklet_23_reverse.setCharge(-1*tracklet_23_reverse.getCharge());
	}

#ifdef _DEBUG_ON
            LogInfo("Using following two tracklets:");
            tracklet_st2.print();
            tracklet_st3.print();
            LogInfo("Yield this combination:");
            tracklet_23.print();
	    //tracklet_23_reverse.print();
#endif
            fitTracklet(tracklet_23); //This is the fit that needs the seeded tx and ty information. Without the seed information, the fit occasionally finds bad slope and X0 or Y0 values, much in the same way that it does for single-station tracklets.  Note from Patrick: this fit could probably be throw out, as we already know the tx and ty information from compareTracklets.  I would just need to extrapolate back to z = 0 and calculate the chisq by hand
            if(tracklet_23.chisq > 9000.)
            {
#ifdef _DEBUG_ON
                tracklet_23.print();
                LogInfo("Impossible combination!");
#endif
                continue;
            }
	    /*
            if(!COARSE_MODE && !hodoMask(tracklet_23))
            {
#ifdef _DEBUG_ON
                LogInfo("Hodomasking failed!");
#endif
                continue;
            }
#ifdef _DEBUG_ON
            LogInfo("Hodomasking Scucess!");
#endif
	    */
            if(!COARSE_MODE)
            {
	      resolveLeftRight(tracklet_23, 40.); //resolveLeftRight at this stage is, in truth, not needed.  We already know which side of the wire the particle passed by from compareTracklet.  However, getting rid of this would take a bit of work
                resolveLeftRight(tracklet_23, 150.);
            }

            ///Remove bad hits if needed
            removeBadHits(tracklet_23);

#ifdef _DEBUG_ON
            LogInfo("New tracklet: ");
            tracklet_23.print();

            LogInfo("Current best:");
            tracklet_best_UX.print();

            LogInfo("Comparison: " << (tracklet_23 < tracklet_best_UX));
            LogInfo("Quality: " << acceptTracklet(tracklet_23));
#endif

            //If current tracklet is better than the best tracklet up-to-now
            if(acceptTracklet(tracklet_23) && tracklet_23 < tracklet_best_UX)
            {
                tracklet_best_UX = tracklet_23;
            }
#ifdef _DEBUG_ON
            else
            {
                LogInfo("Rejected!!");
            }
#endif

	
		/*	
	std::cout<<"FOURTH TEST (*trackletU).acceptedULine3.wireHit1PosZ = "<<(*trackletU).acceptedULine3.wireHit1PosZ<<"; trackletV->acceptedVLine2.wireHit1PosZ = "<<trackletV->acceptedVLine2.wireHit1PosZ<<std::endl;
	
	double posy2U = ((*trackletU).st2U - p_geomSvc->getCostheta((*trackletX).getHit(0).hit.detectorID+2)*((*trackletX).st2X + (*trackletX).st2Xsl * ((*trackletU).st2Z - (*trackletX).st2Z) ))/p_geomSvc->getSintheta((*trackletX).getHit(0).hit.detectorID+2);
	double posy2V = ((*trackletV).st2V - p_geomSvc->getCostheta((*trackletX).getHit(0).hit.detectorID-2)*((*trackletX).st2X + (*trackletX).st2Xsl * ((*trackletV).st2Z - (*trackletX).st2Z) ))/p_geomSvc->getSintheta((*trackletX).getHit(0).hit.detectorID-2);
	double posy3U = ((*trackletU).st3U - p_geomSvc->getCostheta((*trackletX).getHit(2).hit.detectorID+2)*((*trackletX).st3X + (*trackletX).st3Xsl * ((*trackletU).st3Z - (*trackletX).st3Z) ))/p_geomSvc->getSintheta((*trackletX).getHit(2).hit.detectorID+2);
	double posy3V = ((*trackletV).st3V - p_geomSvc->getCostheta((*trackletX).getHit(2).hit.detectorID-2)*((*trackletX).st3X + (*trackletX).st3Xsl * ((*trackletV).st3Z - (*trackletX).st3Z) ))/p_geomSvc->getSintheta((*trackletX).getHit(2).hit.detectorID-2);

	if( std::abs(posy2U) > 100 || std::abs(posy2V) > 100 || std::abs(posy3U) > 100 || std::abs(posy3V) > 100 ) continue;
	std::cout<<"posy2U="<<posy2U<<"; posy2V="<<posy2V<<"; posy3U="<<posy3U<<"; posy3V="<<posy3V<<std::endl;
	
	Tracklet test_tracklet;
	test_tracklet = (*trackletX) + (*trackletU) + (*trackletV);
	test_tracklet.tx = (*trackletX).tx;
	test_tracklet.ty = (((posy3U+posy3V)/2. - (posy2U+posy2V)/2.))/(z_plane[(*trackletX).getHit(2).hit.detectorID] - z_plane[(*trackletX).getHit(0).hit.detectorID]);

	double tracklet2Ys_D = 0.; //This is the sum of the Y-values of the hits in the U and V planes of the two tracklets.  The sum is taken to be used in an average.  The _D here stands for driftDistance.  I originally did this part of the code without taking the drift distance in the slanted layers into account
	double tracklet3Ys_D = 0.;
	tracklet2Ys_D += (*trackletV).acceptedVLine2.wire1Slope * ((*trackletX).acceptedXLine2.slopeX*((*trackletV).acceptedVLine2.wireHit1PosZ - (*trackletX).acceptedXLine2.initialZ) + (*trackletX).acceptedXLine2.initialX) + (*trackletV).acceptedVLine2.wireIntercept1;
	tracklet2Ys_D += (*trackletV).acceptedVLine2.wire2Slope * ((*trackletX).acceptedXLine2.slopeX*((*trackletV).acceptedVLine2.wireHit2PosZ - (*trackletX).acceptedXLine2.initialZ) + (*trackletX).acceptedXLine2.initialX) + (*trackletV).acceptedVLine2.wireIntercept2;
	tracklet2Ys_D += (*trackletU).acceptedULine2.wire1Slope * ((*trackletX).acceptedXLine2.slopeX*((*trackletU).acceptedULine2.wireHit1PosZ - (*trackletX).acceptedXLine2.initialZ) + (*trackletX).acceptedXLine2.initialX) + (*trackletU).acceptedULine2.wireIntercept1;
	tracklet2Ys_D += (*trackletU).acceptedULine2.wire2Slope * ((*trackletX).acceptedXLine2.slopeX*((*trackletU).acceptedULine2.wireHit2PosZ - (*trackletX).acceptedXLine2.initialZ) + (*trackletX).acceptedXLine2.initialX) + (*trackletU).acceptedULine2.wireIntercept2;

	tracklet3Ys_D += (*trackletV).acceptedVLine3.wire1Slope * ((*trackletX).acceptedXLine3.slopeX*((*trackletV).acceptedVLine3.wireHit1PosZ - (*trackletX).acceptedXLine3.initialZ) + (*trackletX).acceptedXLine3.initialX) + (*trackletV).acceptedVLine3.wireIntercept1;
	tracklet3Ys_D += (*trackletV).acceptedVLine3.wire2Slope * ((*trackletX).acceptedXLine3.slopeX*((*trackletV).acceptedVLine3.wireHit2PosZ - (*trackletX).acceptedXLine3.initialZ) + (*trackletX).acceptedXLine3.initialX) + (*trackletV).acceptedVLine3.wireIntercept2;
	tracklet3Ys_D += (*trackletU).acceptedULine3.wire1Slope * ((*trackletX).acceptedXLine3.slopeX*((*trackletU).acceptedULine3.wireHit1PosZ - (*trackletX).acceptedXLine3.initialZ) + (*trackletX).acceptedXLine3.initialX) + (*trackletU).acceptedULine3.wireIntercept1;
	tracklet3Ys_D += (*trackletU).acceptedULine3.wire2Slope * ((*trackletX).acceptedXLine3.slopeX*((*trackletU).acceptedULine3.wireHit2PosZ - (*trackletX).acceptedXLine3.initialZ) + (*trackletX).acceptedXLine3.initialX) + (*trackletU).acceptedULine3.wireIntercept2;
	
	//Give the station 2 and station 3 tracklets the same tx and ty value.  I could get an X0 and Y0 extrapolation, but that doesn't seem to be strictly necessary.  The X0 and Y0 values are found in the fittracklet function for the combined station 2 + station 3 tracklet
	//std::cout<<"tracklet3Ys_D/4. = "<<tracklet3Ys_D/4.<<"; tracklet2Ys_D/4. = "<<tracklet2Ys_D/4.<<"; (*trackletU).acceptedULine3.wireHit1PosZ = "<<(*trackletU).acceptedULine3.wireHit1PosZ<<"; (*trackletV).acceptedVLine2.wireHit1PosZ = "<<(*trackletV).acceptedVLine2.wireHit1PosZ<<std::endl;
	//test_tracklet.ty = (tracklet3Ys_D/4. - tracklet2Ys_D/4.)/((*trackletU).acceptedULine3.wireHit1PosZ - (*trackletV).acceptedVLine2.wireHit1PosZ); //The y slope is found by taking the average Y position in the station3 and subtracting the average Y position in station2.  This is then divided by the z difference, of course  
	
	std::cout<<"tx is = "<<test_tracklet.tx<<" and ty = "<<test_tracklet.ty<<std::endl;
	
	test_tracklet.st2Z = (*trackletX).st2Z;
	test_tracklet.st2X = (*trackletX).st2X;
	test_tracklet.st2Y = tracklet2Ys_D/4.;
	test_tracklet.st3Y = tracklet3Ys_D/4.;
	test_tracklet.st2U = (*trackletU).st2U;
	test_tracklet.st2V = (*trackletV).st2V;
	test_tracklet.st2Usl = (*trackletU).st2Usl;
	test_tracklet.st2Vsl = (*trackletV).st2Vsl;	

#ifdef _DEBUG_ON
	test_tracklet.print();
#endif
	fitTracklet(test_tracklet); //This is the fit that needs the seeded tx and ty information. Without the seed information, the fit occasionally finds bad slope and X0 or Y0 values, much in the same way that it does for single-station tracklets.  Note from Patrick: this fit could probably be throw out, as we already know the tx and ty information from compareTracklets.  I would just need to extrapolate back to z = 0 and calculate the chisq by hand
	if(test_tracklet.chisq > 9000.)
	  {
#ifdef _DEBUG_ON
	    test_tracklet.print();
	    LogInfo("Impossible combination!");
#endif
	    continue;
	  }
	
	if(!COARSE_MODE && !hodoMask(test_tracklet))
	  {
#ifdef _DEBUG_ON
	    LogInfo("Hodomasking failed!");
#endif
	    continue;
	  }
#ifdef _DEBUG_ON
	LogInfo("Hodomasking Scucess!");
#endif
	
	if(!COARSE_MODE)
	  {
	    resolveLeftRight(test_tracklet, 40.); //resolveLeftRight at this stage is, in truth, not needed.  We already know which side of the wire the particle passed by from compareTracklet.  However, getting rid of this would take a bit of work
	    resolveLeftRight(test_tracklet, 150.);
	  }

            ///Remove bad hits if needed
            removeBadHits(test_tracklet);

#ifdef _DEBUG_ON
            LogInfo("New tracklet: ");
            test_tracklet.print();

            LogInfo("Current best:");
            tracklet_best_UX.print();

            LogInfo("Comparison: " << (test_tracklet < tracklet_best_UX));
            LogInfo("Quality: " << acceptTracklet(test_tracklet));
#endif

            //If current tracklet is better than the best tracklet up-to-now
            if(acceptTracklet(test_tracklet) && test_tracklet < tracklet_best_UX)
            {
                tracklet_best_UX = test_tracklet;
            }
#ifdef _DEBUG_ON
            else
            {
                LogInfo("Rejected!!");
            }
#endif*/
      }
      if(tracklet_best_UX < tracklet_best){
	tracklet_best = tracklet_best_UX;
      }
    }
    if(acceptTracklet(tracklet_best)){
      if(tracklet_best.isValid() > 0) trackletsInSt[3].push_back(tracklet_best);   
    }
  }
  
}

void KalmanFastTracking::buildBackPartialTracksSlimX()
{
    //Temporary container for a simple chisq fit
    int nHitsX2, nHitsX3;
    double z_fit[4], x_fit[4];
    double a, b;

    for(std::list<Tracklet>::iterator tracklet3 = trackletsInStSlimX[2].begin(); tracklet3 != trackletsInStSlimX[2].end(); ++tracklet3)
    {
      //tracklet3->print();
      /*if(!COARSE_MODE)
        {
            //Extract the X hits only from station-3 tracks
            nHitsX3 = 0;
            for(std::list<SignedHit>::iterator ptr_hit = tracklet3->hits.begin(); ptr_hit != tracklet3->hits.end(); ++ptr_hit)
            {
	      if(ptr_hit->hit.index < 0) continue;
	      if(p_geomSvc->getPlaneType(ptr_hit->hit.detectorID) == 1)
                {
		  z_fit[nHitsX3] = z_plane[ptr_hit->hit.detectorID];
		  x_fit[nHitsX3] = ptr_hit->hit.pos;
		  ++nHitsX3;
                }
	    }
        }*/

        Tracklet tracklet_best;
        for(std::list<Tracklet>::iterator tracklet2 = trackletsInStSlimX[1].begin(); tracklet2 != trackletsInStSlimX[1].end(); ++tracklet2)
        {
	  //tracklet2->print();
	  /*if(!COARSE_MODE)
            {
	      if(OLD_TRACKING){
		if(fabs(tracklet2->tx - tracklet3->tx) > 0.15 || fabs(tracklet2->ty - tracklet3->ty) > 0.1) continue;
	      }
                //Extract the X hits from station-2 tracke
                nHitsX2 = nHitsX3;
                for(std::list<SignedHit>::iterator ptr_hit = tracklet2->hits.begin(); ptr_hit != tracklet2->hits.end(); ++ptr_hit)
                {
                    if(ptr_hit->hit.index < 0) continue;
                    if(p_geomSvc->getPlaneType(ptr_hit->hit.detectorID) == 1)
                    {
                        z_fit[nHitsX2] = z_plane[ptr_hit->hit.detectorID];
                        x_fit[nHitsX2] = ptr_hit->hit.pos;
                        ++nHitsX2;
                    }
                }

                //Apply a simple linear fit to get rough estimation of X-Z slope and intersection
                chi2fit(nHitsX2, z_fit, x_fit, a, b);
                if(fabs(a) > 2.*TX_MAX || fabs(b) > 2.*X0_MAX) continue;
		
                //Project to proportional tubes to see if there is enough
                int nPropHits = 0;
                for(int i = 0; i < 4; ++i)
                {
                    double x_exp = a*z_mask[detectorIDs_muid[0][i] - nChamberPlanes - 1] + b;
                    for(std::list<int>::iterator iter = hitIDs_muid[0][i].begin(); iter != hitIDs_muid[0][i].end(); ++iter)
                    {
                        if(fabs(hitAll[*iter].pos - x_exp) < 5.08)
                        {
                            ++nPropHits;
                            break;
                        }
                    }
                    if(nPropHits > 0) break;
                }
                if(!TRACK_ELECTRONS && nPropHits == 0) continue; //Turned off by Patrick for electron tracks
            }*/

	    Tracklet tracklet_23;
	    if(OLD_TRACKING){
	      tracklet_23 = (*tracklet2) + (*tracklet3);
	    }
	    else{
	      if(compareTrackletsSlim(*tracklet2, *tracklet3)){
		tracklet_23 = (*tracklet2) + (*tracklet3);
		tracklet_23.tx = tracklet2->tx; //This is needed to "seed" the tracklet fit that happens below.  This tx and ty information is assigned in compareTracklets below
		tracklet_23.st2X = tracklet2->st2X;
		tracklet_23.st2Z = tracklet2->st2Z;
		tracklet_23.st3X = tracklet3->st2X;
		tracklet_23.st3Z = tracklet3->st2Z;
		tracklet_23.st2Xsl = tracklet2->st2Xsl;
		tracklet_23.st3Xsl = tracklet3->st2Xsl;
		tracklet_23.acceptedXLine2 = tracklet2->acceptedXLine2;
		tracklet_23.acceptedXLine3 = tracklet3->acceptedXLine3;
		std::cout<<"tracklet 2 tx is "<<tracklet2->tx<<" and tracklet 3 tx is "<<tracklet2->tx<<std::endl;
		std::cout<<"tracklet 2 st2X is "<<tracklet2->st2X<<" and tracklet 3 st2X is "<<tracklet3->st2X<<std::endl;

#ifdef _DEBUG_ON
            LogInfo("We had a match using following two tracklets:");
            tracklet2->print();
            tracklet3->print();
            LogInfo("Yield this combination:");
            tracklet_23.print();
#endif
		
		trackletsInStSlimX[3].push_back(tracklet_23);
		/*tracklet_23.ty = tracklet2->ty;

		tracklet_23.st2Z = tracklet2->st2Z;
		tracklet_23.st2X = tracklet2->st2X;
		tracklet_23.st2Y = tracklet2->st2Y;
		tracklet_23.st3Y = tracklet2->st3Y;
		tracklet_23.st2U = tracklet2->st2U;
		tracklet_23.st2V = tracklet2->st2V;
		tracklet_23.st2Usl = tracklet2->st2Usl;
		tracklet_23.st2Vsl = tracklet2->st2Vsl;*/
	      }
	      else{
		continue;
	      }
	    }
	}
    }

    std::cout<<"Number of x pairs for 2+3 tracklets is "<<trackletsInStSlimX[3].size()<<std::endl;
    
	    /*            fitTracklet(tracklet_23); //This is the fit that needs the seeded tx and ty information. Without the seed information, the fit occasionally finds bad slope and X0 or Y0 values, much in the same way that it does for single-station tracklets.  Note from Patrick: this fit could probably be throw out, as we already know the tx and ty information from compareTracklets.  I would just need to extrapolate back to z = 0 and calculate the chisq by hand
            if(tracklet_23.chisq > 9000.)
            {
#ifdef _DEBUG_ON
                tracklet_23.print();
                LogInfo("Impossible combination!");
#endif
                continue;
            }

            if(!COARSE_MODE && !hodoMask(tracklet_23))
            {
#ifdef _DEBUG_ON
                LogInfo("Hodomasking failed!");
#endif
                continue;
            }
#ifdef _DEBUG_ON
            LogInfo("Hodomasking Scucess!");
#endif

            if(!COARSE_MODE)
            {
	      resolveLeftRight(tracklet_23, 40.); //resolveLeftRight at this stage is, in truth, not needed.  We already know which side of the wire the particle passed by from compareTracklet.  However, getting rid of this would take a bit of work
                resolveLeftRight(tracklet_23, 150.);
            }

            ///Remove bad hits if needed
            removeBadHits(tracklet_23);

#ifdef _DEBUG_ON
            LogInfo("New tracklet: ");
            tracklet_23.print();

            LogInfo("Current best:");
            tracklet_best.print();

            LogInfo("Comparison: " << (tracklet_23 < tracklet_best));
            LogInfo("Quality: " << acceptTracklet(tracklet_23));
#endif

            //If current tracklet is better than the best tracklet up-to-now
            if(acceptTracklet(tracklet_23) && tracklet_23 < tracklet_best)
            {
                tracklet_best = tracklet_23;
            }
#ifdef _DEBUG_ON
            else
            {
                LogInfo("Rejected!!");
            }
#endif
        }

        if(tracklet_best.isValid() > 0) trackletsInSt[3].push_back(tracklet_best);
    }

    reduceTrackletList(trackletsInSt[3]);
    trackletsInSt[3].sort();*/
}

void KalmanFastTracking::buildBackPartialTracksSlimU()
{
    //Temporary container for a simple chisq fit
    int nHitsX2, nHitsX3;
    double z_fit[4], x_fit[4];
    double a, b;

    for(std::list<Tracklet>::iterator tracklet3 = trackletsInStSlimU[2].begin(); tracklet3 != trackletsInStSlimU[2].end(); ++tracklet3)
    {
        Tracklet tracklet_best;
        for(std::list<Tracklet>::iterator tracklet2 = trackletsInStSlimU[1].begin(); tracklet2 != trackletsInStSlimU[1].end(); ++tracklet2)
        {
	    Tracklet tracklet_23;
	    if(OLD_TRACKING){
	      tracklet_23 = (*tracklet2) + (*tracklet3);
	    }
	    else{
	      if(compareTrackletsSlimU(*tracklet2, *tracklet3)){
		tracklet_23 = (*tracklet2) + (*tracklet3);
		tracklet_23.tx = tracklet2->tx; //This is needed to "seed" the tracklet fit that happens below.  This tx and ty information is assigned in compareTracklets below
		tracklet_23.st2U = tracklet2->st2U;
		tracklet_23.st2Z = tracklet2->st2Z;
		tracklet_23.st3U = tracklet3->st2U;
		tracklet_23.st3Z = tracklet3->st2Z;
		tracklet_23.st2Usl = tracklet2->st2Usl;
		tracklet_23.st3Usl = tracklet3->st2Usl;
		tracklet_23.acceptedULine2 = tracklet2->acceptedULine2;
		tracklet_23.acceptedULine3 = tracklet3->acceptedULine3;
		std::cout<<"tracklet 2 tx is "<<tracklet2->tx<<" and tracklet 3 tx is "<<tracklet2->tx<<std::endl;
		std::cout<<"tracklet 2 st2U is "<<tracklet2->st2U<<" and tracklet 3 st2U is "<<tracklet3->st2U<<std::endl;

#ifdef _DEBUG_ON
            LogInfo("We had a match using following two tracklets:");
            tracklet2->print();
            tracklet3->print();
            LogInfo("Yield this combination:");
            tracklet_23.print();
#endif
		
		trackletsInStSlimU[3].push_back(tracklet_23);
		/*tracklet_23.ty = tracklet2->ty;

		tracklet_23.st2Z = tracklet2->st2Z;
		tracklet_23.st2X = tracklet2->st2X;
		tracklet_23.st2Y = tracklet2->st2Y;
		tracklet_23.st3Y = tracklet2->st3Y;
		tracklet_23.st2U = tracklet2->st2U;
		tracklet_23.st2V = tracklet2->st2V;
		tracklet_23.st2Usl = tracklet2->st2Usl;
		tracklet_23.st2Vsl = tracklet2->st2Vsl;*/
	      }
	      else{
		continue;
	      }
	    }
	}
    }

    std::cout<<"Number of U pairs for 2+3 tracklets is "<<trackletsInStSlimU[3].size()<<std::endl;

}



void KalmanFastTracking::buildBackPartialTracksSlimV()
{
    //Temporary container for a simple chisq fit
    int nHitsX2, nHitsX3;
    double z_fit[4], x_fit[4];
    double a, b;

    for(std::list<Tracklet>::iterator tracklet3 = trackletsInStSlimV[2].begin(); tracklet3 != trackletsInStSlimV[2].end(); ++tracklet3)
    {
        Tracklet tracklet_best;
        for(std::list<Tracklet>::iterator tracklet2 = trackletsInStSlimV[1].begin(); tracklet2 != trackletsInStSlimV[1].end(); ++tracklet2)
        {
	    Tracklet tracklet_23;
	    if(OLD_TRACKING){
	      tracklet_23 = (*tracklet2) + (*tracklet3);
	    }
	    else{
	      if(compareTrackletsSlimV(*tracklet2, *tracklet3)){
		tracklet_23 = (*tracklet2) + (*tracklet3);
		tracklet_23.tx = tracklet2->tx; //This is needed to "seed" the tracklet fit that happens below.  This tx and ty information is assigned in compareTracklets below
		tracklet_23.st2V = tracklet2->st2V;
		tracklet_23.st2Z = tracklet2->st2Z;
		tracklet_23.st3V = tracklet3->st2V;
		tracklet_23.st3Z = tracklet3->st2Z;
		tracklet_23.st2Vsl = tracklet2->st2Vsl;
		tracklet_23.st3Vsl = tracklet3->st2Vsl;
		std::cout<<"SECOND TEST OF LINE2 "<<tracklet2->acceptedVLine2.wireHit1PosZ<<std::endl;
		tracklet_23.acceptedVLine2 = tracklet2->acceptedVLine2;
		std::cout<<"THIRD TEST OF LINE2 "<<tracklet_23.acceptedVLine2.wireHit1PosZ<<std::endl;
		tracklet_23.acceptedVLine2.print();
		tracklet_23.acceptedVLine3 = tracklet3->acceptedVLine3;
		std::cout<<"tracklet 2 tx is "<<tracklet2->tx<<" and tracklet 3 tx is "<<tracklet2->tx<<std::endl;
		std::cout<<"tracklet 2 st2V is "<<tracklet2->st2V<<" and tracklet 3 st2V is "<<tracklet3->st2V<<std::endl;

#ifdef _DEBUG_ON
            LogInfo("We had a match using following two tracklets:");
            tracklet2->print();
            tracklet3->print();
            LogInfo("Yield this combination:");
            tracklet_23.print();
#endif
		
		trackletsInStSlimV[3].push_back(tracklet_23);
		/*tracklet_23.ty = tracklet2->ty;

		tracklet_23.st2Z = tracklet2->st2Z;
		tracklet_23.st2X = tracklet2->st2X;
		tracklet_23.st2Y = tracklet2->st2Y;
		tracklet_23.st3Y = tracklet2->st3Y;
		tracklet_23.st2U = tracklet2->st2U;
		tracklet_23.st2V = tracklet2->st2V;
		tracklet_23.st2Usl = tracklet2->st2Usl;
		tracklet_23.st2Vsl = tracklet2->st2Vsl;*/
	      }
	      else{
		continue;
	      }
	    }
	}
    }

    std::cout<<"Number of V pairs for 2+3 tracklets is "<<trackletsInStSlimV[3].size()<<std::endl;

}

void KalmanFastTracking::buildGlobalTracks()
{
    double pos_exp[3], window[3];
    for(std::list<Tracklet>::iterator tracklet23 = trackletsInSt[3].begin(); tracklet23 != trackletsInSt[3].end(); ++tracklet23)
    { 
        Tracklet tracklet_best[2];
        for(int i = 0; i < 2; ++i) //for two station-1 chambers
        {
            trackletsInSt[0].clear();
	    if(!TRACK_DISPLACED){
	      //Calculate the window in station 1
	      if(KMAG_ON)
		{
		  getSagittaWindowsInSt1(*tracklet23, pos_exp, window, i+1);
		}
	      else
		{
		  getExtrapoWindowsInSt1(*tracklet23, pos_exp, window, i+1);
		}
	    }
#ifdef _DEBUG_ON
            LogInfo("Using this back partial: ");
            tracklet23->print();
            for(int j = 0; j < 3; j++) LogInfo("Extrapo: " << pos_exp[j] << "  " << window[j]);
#endif

            _timers["global_st1"]->restart();
            if(!TRACK_DISPLACED){
	      buildTrackletsInStation(i+1, 0, pos_exp, window);
	    }
	    if(TRACK_DISPLACED){
	      buildTrackletsInStation(i+1, 0);
	    }

            _timers["global_st1"]->stop();

            _timers["global_link"]->restart();
            Tracklet tracklet_best_prob, tracklet_best_vtx;
            for(std::list<Tracklet>::iterator tracklet1 = trackletsInSt[0].begin(); tracklet1 != trackletsInSt[0].end(); ++tracklet1)
            {
#ifdef _DEBUG_ON
                LogInfo("With this station 1 track:");
                tracklet1->print();
#endif

                Tracklet tracklet_global = (*tracklet23) * (*tracklet1);
                fitTracklet(tracklet_global);
                if(!hodoMask(tracklet_global)) continue;

                ///Resolve the left-right with a tight pull cut, then a loose one, then resolve by single projections
                if(!COARSE_MODE)
                {
                    resolveLeftRight(tracklet_global, 75.);
                    resolveLeftRight(tracklet_global, 150.);
                    resolveSingleLeftRight(tracklet_global);
                }
		if(TRACK_DISPLACED){
		  double firstChiSq = tracklet_global.calcChisq();
		  Tracklet tracklet_global2 = (*tracklet23) * (*tracklet1);
		  tracklet_global2.setCharge(-1*tracklet_global2.getCharge()); /**By default, the value returned by getCharge is based on the x0 of the tracklet.  For a particle produced at the target, this is a valid way to extract charge.  However, for a displaced particle, the x0 does not tell you anything useful about the charge of the particle, which is why we need to check both possible charge values.  getCharge is used later on when extracting certain track quality values, so using the wrong charge leads to tracks getting rejected due to poor quality values*/
		  if(!COARSE_MODE)
		    {
		      resolveLeftRight(tracklet_global2, 75.);
		      resolveLeftRight(tracklet_global2, 150.);
		      resolveSingleLeftRight(tracklet_global2);
		    }
		  double secondChiSq = tracklet_global2.calcChisq();
		  if(secondChiSq < firstChiSq){
		    tracklet_global = tracklet_global2;
		  }
		}
		
                ///Remove bad hits if needed
                removeBadHits(tracklet_global);

                //Most basic cuts
                if(!acceptTracklet(tracklet_global)) continue;

                //Get the tracklets that has the best prob
                if(tracklet_global < tracklet_best_prob) tracklet_best_prob = tracklet_global;

                ///Set vertex information - only applied when KF is enabled
                ///TODO: maybe in the future add a Genfit-based equivalent here, for now leave as is
                if(enable_KF && NOT_DISPLACED)
                {
                    _timers["global_kalman"]->restart();
                    SRecTrack recTrack = processOneTracklet(tracklet_global);
                    _timers["global_kalman"]->stop();
                    tracklet_global.chisq_vtx = recTrack.getChisqVertex();

                    if(recTrack.isValid() && tracklet_global.chisq_vtx < tracklet_best_vtx.chisq_vtx) tracklet_best_vtx = tracklet_global;
                }

#ifdef _DEBUG_ON
                LogInfo("New tracklet: ");
                tracklet_global.print();

                LogInfo("Current best by prob:");
                tracklet_best_prob.print();

                LogInfo("Comparison I: " << (tracklet_global < tracklet_best_prob));
                LogInfo("Quality I   : " << acceptTracklet(tracklet_global));

                if(enable_KF && NOT_DISPLACED)
                {
                    LogInfo("Current best by vtx:");
                    tracklet_best_vtx.print();

                    LogInfo("Comparison II: " << (tracklet_global.chisq_vtx < tracklet_best_vtx.chisq_vtx));
                    //LogInfo("Quality II   : " << recTrack.isValid());
                }
#endif
            }
            _timers["global_link"]->stop();

            //The selection logic is, prefer the tracks with best p-value, as long as it's not low-pz
            if(enable_KF && NOT_DISPLACED && tracklet_best_prob.isValid() > 0 && 1./tracklet_best_prob.invP > 18.)
            {
                tracklet_best[i] = tracklet_best_prob;
            }
            else if(enable_KF && NOT_DISPLACED && tracklet_best_vtx.isValid() > 0) //otherwise select the one with best vertex chisq, TODO: maybe add a z-vtx constraint
            {
                tracklet_best[i] = tracklet_best_vtx;
            }
            else if(tracklet_best_prob.isValid() > 0) //then fall back to the default only choice
            {
                tracklet_best[i] = tracklet_best_prob;
            }
        }

        //Merge the tracklets from two stations if necessary
        Tracklet tracklet_merge;
        if(fabs(tracklet_best[0].getMomentum() - tracklet_best[1].getMomentum())/tracklet_best[0].getMomentum() < MERGE_THRES)
        {
            //Merge the track and re-fit
            tracklet_merge = tracklet_best[0].merge(tracklet_best[1]);
            fitTracklet(tracklet_merge);

#ifdef _DEBUG_ON
            LogInfo("Merging two track candidates with momentum: " << tracklet_best[0].getMomentum() << "  " << tracklet_best[1].getMomentum());
            LogInfo("tracklet_best_1:"); tracklet_best[0].print();
            LogInfo("tracklet_best_2:"); tracklet_best[1].print();
            LogInfo("tracklet_merge:"); tracklet_merge.print();
#endif
        }

        if(tracklet_merge.isValid() > 0 && tracklet_merge < tracklet_best[0] && tracklet_merge < tracklet_best[1])
        {
#ifdef _DEBUG_ON
            LogInfo("Choose merged tracklet");
#endif
            trackletsInSt[4].push_back(tracklet_merge);
        }
        else if(tracklet_best[0].isValid() > 0 && tracklet_best[0] < tracklet_best[1])
        {
#ifdef _DEBUG_ON
            LogInfo("Choose tracklet with station-0");
#endif
            trackletsInSt[4].push_back(tracklet_best[0]);
        }
        else if(tracklet_best[1].isValid() > 0)
        {
#ifdef _DEBUG_ON
            LogInfo("Choose tracklet with station-1");
#endif
            trackletsInSt[4].push_back(tracklet_best[1]);
        }
    }

    trackletsInSt[4].sort();
}


void KalmanFastTracking::buildGlobalTracksDisplaced()
{
    double pos_exp[3], window[3];
    for(std::list<Tracklet>::iterator tracklet23 = trackletsInSt[3].begin(); tracklet23 != trackletsInSt[3].end(); ++tracklet23)
    {
      std::cout<<"I'm testing backwards extrapolation... z_plane thing is "<<z_plane[0]<<" "<<z_plane[1]<<" "<<z_plane[2]<<" "<<z_plane[3]<<" "<<z_plane[4]<<" "<<z_plane[5]<<" "<<z_plane[6]<<" "<<z_plane[7]<<std::endl; //WPM
      getSagittaWindowsInSt1(*tracklet23, pos_exp, window, 1); //WPM as a printing test
      double posx = tracklet23->tx * ( z_plane[3] - tracklet23->st2Z ) + tracklet23->st2X; //WPM pick up here.  do u and v extrapolations
      std::cout<<"test of u extrapo.  st2Usl = "<<tracklet23->st2Usl<<" ( z_plane[(i)*6 + 0] - tracklet23->st2UZ ) = "<<( z_plane[3] - tracklet23->st2UZ )<<" tracklet23->st2U = "<<tracklet23->st2U<<std::endl; //WPM
      double posu = tracklet23->st2Usl * ( z_plane[1] - tracklet23->st2Z ) + tracklet23->st2U; //WPM
      double posv = tracklet23->st2Vsl * ( z_plane[5] - tracklet23->st2Z ) + tracklet23->st2V; //WPM
      double posy = tracklet23->ty * ( z_plane[3] - tracklet23->st2Z ) + tracklet23->st2Y; //WPM
      std::cout<<"checking V extrapolation.  tracklet23->st2Vsl: "<<tracklet23->st2Vsl<<", ( z_plane[5] - tracklet23->st2Z ): "<<( z_plane[5] - tracklet23->st2Z )<<", tracklet23->st2V: "<<tracklet23->st2V<<std::endl; //WPM
      std::cout<<"the extrapolations are x: "<<posx<<" and u: "<<posu<<" and v: "<<posv<<std::endl; //WPM
      std::cout<<"the y-extrapolation is "<<tracklet23->ty * ( z_plane[3] - tracklet23->st2Z ) + tracklet23->st2Y<<std::endl; //WPM
      std::cout<<"quick printout of wire cos theta: "<<p_geomSvc->getCostheta(1)<<", "<<p_geomSvc->getCostheta(2)<<", "<<p_geomSvc->getCostheta(3)<<", "<<p_geomSvc->getCostheta(4)<<", "<<p_geomSvc->getCostheta(5)<<", "<<p_geomSvc->getCostheta(6)<<std::endl;
      std::cout<<"quick printout of wire sin theta: "<<p_geomSvc->getSintheta(1)<<", "<<p_geomSvc->getSintheta(2)<<", "<<p_geomSvc->getSintheta(3)<<", "<<p_geomSvc->getSintheta(4)<<", "<<p_geomSvc->getSintheta(5)<<", "<<p_geomSvc->getSintheta(6)<<std::endl;
      
        Tracklet tracklet_best[2];
        for(int i = 0; i < 1; ++i) //for two station-1 chambers //WPM edited so that it only does one chamber.  I think that's all we're using in our simulation...
        {

	  bool validTrackFound = false; //WPM
	  int pxSlices[6] = {1, 3, 5, 7, 9, 11}; //WPM
	  for(int pxs = 0; pxs < 6; pxs++){ //WPM
	    if(validTrackFound) continue; //WPM potentially controversial.  Trying to get out of px window loop
	    
	    int charges[2] = {-1,1}; //WPM
	    for(int ch = 0; ch < 2; ch++){ //WPM
	      
	      if(validTrackFound) continue; //WPM potentially controversial.  Trying to get out of px window loop
	      std::cout<<"Getting st1 tracklets centered at "<<charges[ch]*pxSlices[pxs]<<std::endl; //WPM
	      
	    trackletsInSt[0].clear();
	    if(!TRACK_DISPLACED){
	      //Calculate the window in station 1
	      if(KMAG_ON)
		{
		  getSagittaWindowsInSt1(*tracklet23, pos_exp, window, i+1);
		}
	      else
		{
		  getExtrapoWindowsInSt1(*tracklet23, pos_exp, window, i+1);
		}
	    }
#ifdef _DEBUG_ON
            LogInfo("Using this back partial: ");
            tracklet23->print();
            for(int j = 0; j < 3; j++) LogInfo("Extrapo: " << pos_exp[j] << "  " << window[j]);
#endif

            _timers["global_st1"]->restart();
            if(!TRACK_DISPLACED){
	      buildTrackletsInStation(i+1, 0, pos_exp, window);
	    }
	    if(TRACK_DISPLACED){
	      (*tracklet23).setCharge(charges[ch]); //WPM
	      tracklet23->print();
	      pos_exp[0] = posx+charges[ch]*pxSlices[pxs];
	      pos_exp[1] = p_geomSvc->getCostheta(1) * pos_exp[0] + p_geomSvc->getSintheta(1) * posy;
	      pos_exp[2] = p_geomSvc->getCostheta(5) * pos_exp[0] + p_geomSvc->getSintheta(5) * posy;
	      std::cout<<"The window I'm using is centered at "<<posx+charges[ch]*pxSlices[pxs]<<", "<<p_geomSvc->getCostheta(1) * pos_exp[0] + p_geomSvc->getSintheta(1) * posy<<", "<<p_geomSvc->getCostheta(5) * pos_exp[0] + p_geomSvc->getSintheta(5) * posy<<std::endl; //WPM
	      window[0] = 1;
	      window[1] = 3;
	      window[2] = 3;
	      buildTrackletsInStation(i+1, 0, pos_exp, window);
	    }

	    std::cout<<"HELLO THERE: number of station 1 tracklets for index "<<i<<" is "<<trackletsInSt[0].size()<<std::endl; //WPM
	    //std::cout<<"the last tracklet is "<<trackletsInSt[0].end()->print()<<std::endl; //WPM
            _timers["global_st1"]->stop();

            _timers["global_link"]->restart();
	    int tracklet_counter = 0; //WPM
            Tracklet tracklet_best_prob, tracklet_best_vtx;
            for(std::list<Tracklet>::iterator tracklet1 = trackletsInSt[0].begin(); tracklet1 != trackletsInSt[0].end(); ++tracklet1)
            {
	      tracklet_counter++; //WPM
	      std::cout<<"try tracklet number "<<tracklet_counter<<std::endl; //WPM
#ifdef _DEBUG_ON
                LogInfo("With this station 1 track:");
                tracklet1->print();
#endif

                Tracklet tracklet_global = (*tracklet23) * (*tracklet1);
                fitTracklet(tracklet_global);
                if(!hodoMask(tracklet_global)) continue;

                ///Resolve the left-right with a tight pull cut, then a loose one, then resolve by single projections
                if(!COARSE_MODE)
                {
                    resolveLeftRight(tracklet_global, 75.);
                    resolveLeftRight(tracklet_global, 150.);
                    resolveSingleLeftRight(tracklet_global);
                }
		/*if(TRACK_DISPLACED){
		  double firstChiSq = tracklet_global.calcChisq();
		  Tracklet tracklet_global2 = (*tracklet23) * (*tracklet1);
		  tracklet_global2.setCharge(-1*tracklet_global2.getCharge()); //By default, the value returned by getCharge is based on the x0 of the tracklet.  For a particle produced at the target, this is a valid way to extract charge.  However, for a displaced particle, the x0 does not tell you anything useful about the charge of the particle, which is why we need to check both possible charge values.  getCharge is used later on when extracting certain track quality values, so using the wrong charge leads to tracks getting rejected due to poor quality values
		  if(!COARSE_MODE)
		    {
		      resolveLeftRight(tracklet_global2, 75.);
		      resolveLeftRight(tracklet_global2, 150.);
		      resolveSingleLeftRight(tracklet_global2);
		    }
		  double secondChiSq = tracklet_global2.calcChisq();
		  if(secondChiSq < firstChiSq){
		    tracklet_global = tracklet_global2;
		  }
		}*/
		
                ///Remove bad hits if needed
                removeBadHits(tracklet_global);

                //Most basic cuts
                if(!acceptTracklet(tracklet_global)) continue;

                //Get the tracklets that has the best prob
                if(tracklet_global < tracklet_best_prob) tracklet_best_prob = tracklet_global;

                ///Set vertex information - only applied when KF is enabled
                ///TODO: maybe in the future add a Genfit-based equivalent here, for now leave as is
                if(enable_KF && NOT_DISPLACED)
                {
                    _timers["global_kalman"]->restart();
                    SRecTrack recTrack = processOneTracklet(tracklet_global);
                    _timers["global_kalman"]->stop();
                    tracklet_global.chisq_vtx = recTrack.getChisqVertex();

                    if(recTrack.isValid() && tracklet_global.chisq_vtx < tracklet_best_vtx.chisq_vtx) tracklet_best_vtx = tracklet_global;
                }

#ifdef _DEBUG_ON
                LogInfo("New tracklet: ");
                tracklet_global.print();

                LogInfo("Current best by prob:");
                tracklet_best_prob.print();

                LogInfo("Comparison I: " << (tracklet_global < tracklet_best_prob));
                LogInfo("Quality I   : " << acceptTracklet(tracklet_global));

                if(enable_KF && NOT_DISPLACED)
                {
                    LogInfo("Current best by vtx:");
                    tracklet_best_vtx.print();

                    LogInfo("Comparison II: " << (tracklet_global.chisq_vtx < tracklet_best_vtx.chisq_vtx));
                    //LogInfo("Quality II   : " << recTrack.isValid());
                }
#endif
            }
            _timers["global_link"]->stop();

            //The selection logic is, prefer the tracks with best p-value, as long as it's not low-pz
            if(enable_KF && NOT_DISPLACED && tracklet_best_prob.isValid() > 0 && 1./tracklet_best_prob.invP > 18.)
            {
                tracklet_best[i] = tracklet_best_prob;
		validTrackFound = true;
            }
            else if(enable_KF && NOT_DISPLACED && tracklet_best_vtx.isValid() > 0) //otherwise select the one with best vertex chisq, TODO: maybe add a z-vtx constraint
            {
                tracklet_best[i] = tracklet_best_vtx;
		validTrackFound = true;
            }
            else if(tracklet_best_prob.isValid() > 0) //then fall back to the default only choice
            {
                tracklet_best[i] = tracklet_best_prob;
		validTrackFound = true;
            }
	    }
	  }
	}

        //Merge the tracklets from two stations if necessary
        Tracklet tracklet_merge;
        if(fabs(tracklet_best[0].getMomentum() - tracklet_best[1].getMomentum())/tracklet_best[0].getMomentum() < MERGE_THRES)
        {
            //Merge the track and re-fit
            tracklet_merge = tracklet_best[0].merge(tracklet_best[1]);
            fitTracklet(tracklet_merge);

#ifdef _DEBUG_ON
            LogInfo("Merging two track candidates with momentum: " << tracklet_best[0].getMomentum() << "  " << tracklet_best[1].getMomentum());
            LogInfo("tracklet_best_1:"); tracklet_best[0].print();
            LogInfo("tracklet_best_2:"); tracklet_best[1].print();
            LogInfo("tracklet_merge:"); tracklet_merge.print();
#endif
        }

        if(tracklet_merge.isValid() > 0 && tracklet_merge < tracklet_best[0] && tracklet_merge < tracklet_best[1])
        {
#ifdef _DEBUG_ON
            LogInfo("Choose merged tracklet");
#endif
            trackletsInSt[4].push_back(tracklet_merge);
        }
        else if(tracklet_best[0].isValid() > 0 && tracklet_best[0] < tracklet_best[1])
        {
#ifdef _DEBUG_ON
            LogInfo("Choose tracklet with station-0");
#endif
            trackletsInSt[4].push_back(tracklet_best[0]);
        }
        else if(tracklet_best[1].isValid() > 0)
        {
#ifdef _DEBUG_ON
            LogInfo("Choose tracklet with station-1");
#endif
            trackletsInSt[4].push_back(tracklet_best[1]);
        }
    }

    trackletsInSt[4].sort();
}

void KalmanFastTracking::resolveLeftRight(Tracklet& tracklet, double threshold)
{
#ifdef _DEBUG_ON
    LogInfo("Left right for this track..");
    tracklet.print();
#endif

    //Check if the track has been updated
    bool isUpdated = false;

    //Four possibilities
    int possibility[4][2] = {{1, 1}, {1, -1}, {-1, 1}, {-1, -1}};

    //Total number of hit pairs in this tracklet
    int nPairs = tracklet.hits.size()/2;

    int nResolved = 0;
    std::list<SignedHit>::iterator hit1 = tracklet.hits.begin();
    std::list<SignedHit>::iterator hit2 = tracklet.hits.begin();
    ++hit2;
    while(true)
    {
#ifdef _DEBUG_ON
        LogInfo(hit1->hit.index << "  " << hit2->sign << " === " << hit2->hit.index << "  " << hit2->sign);
        int detectorID1 = hit1->hit.detectorID;
        int detectorID2 = hit2->hit.detectorID;
        LogInfo("Hit1: " << tracklet.getExpPositionX(z_plane[detectorID1])*costheta_plane[detectorID1] + tracklet.getExpPositionY(z_plane[detectorID1])*sintheta_plane[detectorID1] << "  " << hit1->hit.pos + hit1->hit.driftDistance << "  " << hit1->hit.pos - hit1->hit.driftDistance);
        LogInfo("Hit2: " << tracklet.getExpPositionX(z_plane[detectorID2])*costheta_plane[detectorID2] + tracklet.getExpPositionY(z_plane[detectorID2])*sintheta_plane[detectorID2] << "  " << hit2->hit.pos + hit2->hit.driftDistance << "  " << hit2->hit.pos - hit2->hit.driftDistance);
#endif

        if(hit1->hit.index > 0 && hit2->hit.index > 0 && hit1->sign*hit2->sign == 0)
        {
            int index_min = -1;
            double pull_min = 1E6;
            for(int i = 0; i < 4; i++)
            {
                double slope_local = (hit1->pos(possibility[i][0]) - hit2->pos(possibility[i][1]))/(z_plane[hit1->hit.detectorID] - z_plane[hit2->hit.detectorID]);
                double inter_local = hit1->pos(possibility[i][0]) - slope_local*z_plane[hit1->hit.detectorID];

                if(fabs(slope_local) > slope_max[hit1->hit.detectorID] || fabs(inter_local) > intersection_max[hit1->hit.detectorID]) continue;

                double tx, ty, x0, y0;
                double err_tx, err_ty, err_x0, err_y0;
		if(tracklet.stationID == 7 && hit1->hit.detectorID <= 6)
                {
                    tracklet.getXZInfoInSt1(tx, x0);
                    tracklet.getXZErrorInSt1(err_tx, err_x0);
                }
                else
                {
                    tx = tracklet.tx;
                    x0 = tracklet.x0;
                    err_tx = tracklet.err_tx;
                    err_x0 = tracklet.err_x0;
                }
                ty = tracklet.ty;
                y0 = tracklet.y0;
                err_ty = tracklet.err_ty;
                err_y0 = tracklet.err_y0;
		
                double slope_exp = costheta_plane[hit1->hit.detectorID]*tx + sintheta_plane[hit1->hit.detectorID]*ty;
                double err_slope = fabs(costheta_plane[hit1->hit.detectorID]*err_tx) + fabs(sintheta_plane[hit2->hit.detectorID]*err_ty);
                double inter_exp = costheta_plane[hit1->hit.detectorID]*x0 + sintheta_plane[hit1->hit.detectorID]*y0;
                double err_inter = fabs(costheta_plane[hit1->hit.detectorID]*err_x0) + fabs(sintheta_plane[hit2->hit.detectorID]*err_y0);

                double pull = sqrt((slope_exp - slope_local)*(slope_exp - slope_local)/err_slope/err_slope + (inter_exp - inter_local)*(inter_exp - inter_local)/err_inter/err_inter);
                if(pull < pull_min)
                {
                    index_min = i;
                    pull_min = pull;
                }

#ifdef _DEBUG_ON
                LogInfo(hit1->hit.detectorID << ": " << i << "  " << possibility[i][0] << "  " << possibility[i][1]);
                LogInfo(tx << "  " << x0 << "  " << ty << "  " << y0);
                LogInfo("Slope: " << slope_local << "  " << slope_exp << "  " << err_slope);
                LogInfo("Intersection: " << inter_local << "  " << inter_exp << "  " << err_inter);
                LogInfo("Current: " << pull << "  " << index_min << "  " << pull_min);
#endif
            }

            //LogInfo("Final: " << index_min << "  " << pull_min);
            if(index_min >= 0 && pull_min < threshold)//((tracklet.stationID == 5 && pull_min < 25.) || (tracklet.stationID == 6 && pull_min < 100.)))
            {
                hit1->sign = possibility[index_min][0];
                hit2->sign = possibility[index_min][1];
                isUpdated = true;
            }
        }

        ++nResolved;
        if(nResolved >= nPairs) break;

        ++hit1;
        ++hit1;
        ++hit2;
        ++hit2;
    }

    if(isUpdated) fitTracklet(tracklet);
}

void KalmanFastTracking::resolveSingleLeftRight(Tracklet& tracklet)
{
#ifdef _DEBUG_ON
    LogInfo("Single left right for this track..");
    tracklet.print();
#endif

    //Check if the track has been updated
    bool isUpdated = false;
    for(std::list<SignedHit>::iterator hit_sign = tracklet.hits.begin(); hit_sign != tracklet.hits.end(); ++hit_sign)
    {
        if(hit_sign->hit.index < 0 || hit_sign->sign != 0) continue;

        int detectorID = hit_sign->hit.detectorID;
        double pos_exp = tracklet.getExpPositionX(z_plane[detectorID])*costheta_plane[detectorID] + tracklet.getExpPositionY(z_plane[detectorID])*sintheta_plane[detectorID];
        hit_sign->sign = pos_exp > hit_sign->hit.pos ? 1 : -1;

        isUpdated = true;
    }

    if(isUpdated) fitTracklet(tracklet);
}

void KalmanFastTracking::removeBadHits(Tracklet& tracklet)
{
#ifdef _DEBUG_ON
    LogInfo("Removing hits for this track..");
    tracklet.calcChisq();
    tracklet.print();
#endif

    //Check if the track has beed updated
    int signflipflag[nChamberPlanes];
    for(int i = 0; i < nChamberPlanes; ++i) signflipflag[i] = 0;

    bool isUpdated = true;
    while(isUpdated)
    {
        isUpdated = false;
        tracklet.calcChisq();

        SignedHit* hit_remove = nullptr;
        SignedHit* hit_neighbour = nullptr;
        double res_remove1 = -1.;
        double res_remove2 = -1.;
        for(std::list<SignedHit>::iterator hit_sign = tracklet.hits.begin(); hit_sign != tracklet.hits.end(); ++hit_sign)
        {
            if(hit_sign->hit.index < 0) continue;

            int detectorID = hit_sign->hit.detectorID;
            double res_curr = fabs(tracklet.residual[detectorID-1]);
            if(res_remove1 < res_curr)
            {
                res_remove1 = res_curr;
                res_remove2 = fabs(tracklet.residual[detectorID-1] - 2.*hit_sign->sign*hit_sign->hit.driftDistance);
                hit_remove = &(*hit_sign);

                std::list<SignedHit>::iterator iter = hit_sign;
                hit_neighbour = detectorID % 2 == 0 ? &(*(--iter)) : &(*(++iter));
            }
        }
        if(hit_remove == nullptr) continue;
        if(hit_remove->sign == 0 && tracklet.isValid() > 0) continue;  //if sign is undecided, and chisq is OKay, then pass

        double cut = hit_remove->sign == 0 ? hit_remove->hit.driftDistance + resol_plane[hit_remove->hit.detectorID] : resol_plane[hit_remove->hit.detectorID];
        if(res_remove1 > cut)
        {
#ifdef _DEBUG_ON
            LogInfo("Dropping this hit: " << res_remove1 << "  " << res_remove2 << "   " << signflipflag[hit_remove->hit.detectorID-1] << "  " << cut);
            hit_remove->hit.print();
            hit_neighbour->hit.print();
#endif

            //can only be changed less than twice
            if(res_remove2 < cut && signflipflag[hit_remove->hit.detectorID-1] < 2)
            {
                hit_remove->sign = -hit_remove->sign;
                hit_neighbour->sign = 0;
                ++signflipflag[hit_remove->hit.detectorID-1];
#ifdef _DEBUG_ON
                LogInfo("Only changing the sign.");
#endif
            }
            else
            {
                //Set the index of the hit to be removed to -1 so it's not used anymore
                //also set the sign assignment of the neighbour hit to 0 (i.e. undecided)
                hit_remove->hit.index = -1;
                hit_neighbour->sign = 0;
                int planeType = p_geomSvc->getPlaneType(hit_remove->hit.detectorID);
                if(planeType == 1)
                {
                    --tracklet.nXHits;
                }
                else if(planeType == 2)
                {
                    --tracklet.nUHits;
                }
                else
                {
                    --tracklet.nVHits;
                }

                //If both hit pairs are not included, the track can be rejected
                if(hit_neighbour->hit.index < 0)
                {
#ifdef _DEBUG_ON
                    LogInfo("Both hits in a view are missing! Will exit the bad hit removal...");
#endif
                    return;
                }
            }
            isUpdated = true;
        }

        if(isUpdated)
        {
            fitTracklet(tracklet);
            resolveSingleLeftRight(tracklet);
        }
    }
}

void KalmanFastTracking::resolveLeftRight(SRawEvent::hit_pair hpair, int& LR1, int& LR2)
{
    LR1 = 0;
    LR2 = 0;

    //If either hit is missing, no left-right can be assigned
    if(hpair.first < 0 || hpair.second < 0)
    {
        return;
    }

    int possibility[4][2] = {{1, 1}, {1, -1}, {-1, 1}, {-1, -1}};
    int nResolved = 0;
    for(int i = 0; i < 4; i++)
    {
        if(nResolved > 1) break;

        int hitID1 = hpair.first;
        int hitID2 = hpair.second;
        double slope_local = (hitAll[hitID1].pos + possibility[i][0]*hitAll[hitID1].driftDistance - hitAll[hitID2].pos - possibility[i][1]*hitAll[hitID2].driftDistance)/(z_plane[hitAll[hitID1].detectorID] - z_plane[hitAll[hitID2].detectorID]);
        double intersection_local = hitAll[hitID1].pos + possibility[i][0]*hitAll[hitID1].driftDistance - slope_local*z_plane[hitAll[hitID1].detectorID];

        //LogInfo(i << "  " << nResolved << "  " << slope_local << "  " << intersection_local);
        if(fabs(slope_local) < slope_max[hitAll[hitID1].detectorID] && fabs(intersection_local) < intersection_max[hitAll[hitID1].detectorID])
        {
            nResolved++;
            LR1 = possibility[i][0];
            LR2 = possibility[i][1];
        }
    }

    if(nResolved > 1)
    {
        LR1 = 0;
        LR2 = 0;
    }

    //LogInfo("Final: " << LR1 << "  " << LR2);
}

void KalmanFastTracking::buildTrackletsInStation(int stationID, int listID, double* pos_exp, double* window)
{
#ifdef _DEBUG_ON
    LogInfo("Building tracklets in station " << stationID);
#endif

    //actuall ID of the tracklet lists
    int sID = stationID - 1;

    //Extract the X, U, V hit pairs
    std::list<SRawEvent::hit_pair> pairs_X, pairs_U, pairs_V;
    if(pos_exp == nullptr)
    {
        pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0]);
        pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1]);
        pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2]);
    }
    else
    {
        //Note that in pos_exp[], index 0 stands for X, index 1 stands for U, index 2 stands for V
        pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0], pos_exp[0], window[0]);
        pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1], pos_exp[1], window[1]);
        pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2], pos_exp[2], window[2]);
    }

#ifdef _DEBUG_ON
    LogInfo("Hit pairs in this event: ");
    for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_X.begin(); iter != pairs_X.end(); ++iter) LogInfo("X :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
    for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_U.begin(); iter != pairs_U.end(); ++iter) LogInfo("U :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
    for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_V.begin(); iter != pairs_V.end(); ++iter) LogInfo("V :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
#endif

    if(pairs_X.empty() || pairs_U.empty() || pairs_V.empty())
    {
#ifdef _DEBUG_ON
        LogInfo("Not all view has hits in station " << stationID);
#endif
        return;
    }

    //X-U combination first, then add V pairs
    for(std::list<SRawEvent::hit_pair>::iterator xiter = pairs_X.begin(); xiter != pairs_X.end(); ++xiter)
    {
        //U projections from X plane
        double x_pos = xiter->second >= 0 ? 0.5*(hitAll[xiter->first].pos + hitAll[xiter->second].pos) : hitAll[xiter->first].pos;
        double u_min = x_pos*u_costheta[sID] - u_win[sID];
        double u_max = u_min + 2.*u_win[sID];

#ifdef _DEBUG_ON
        LogInfo("Trying X hits " << xiter->first << "  " << xiter->second << "  " << hitAll[xiter->first].elementID << " at " << x_pos);
        LogInfo("U plane window:" << u_min << "  " << u_max);
#endif
        for(std::list<SRawEvent::hit_pair>::iterator uiter = pairs_U.begin(); uiter != pairs_U.end(); ++uiter)
        {
            double u_pos = uiter->second >= 0 ? 0.5*(hitAll[uiter->first].pos + hitAll[uiter->second].pos) : hitAll[uiter->first].pos;
#ifdef _DEBUG_ON
            LogInfo("Trying U hits " << uiter->first << "  " << uiter->second << "  " << hitAll[uiter->first].elementID << " at " << u_pos);
#endif
            if(u_pos < u_min || u_pos > u_max) continue;

            //V projections from X and U plane
            double z_x = xiter->second >= 0 ? z_plane_x[sID] : z_plane[hitAll[xiter->first].detectorID];
            double z_u = uiter->second >= 0 ? z_plane_u[sID] : z_plane[hitAll[uiter->first].detectorID];
            double z_v = z_plane_v[sID];
            double v_win1 = spacing_plane[hitAll[uiter->first].detectorID]*2.*u_costheta[sID];
            double v_win2 = fabs((z_u + z_v - 2.*z_x)*u_costheta[sID]*TX_MAX);
            double v_win3 = fabs((z_v - z_u)*u_sintheta[sID]*TY_MAX);
            double v_win = v_win1 + v_win2 + v_win3 + 2.*spacing_plane[hitAll[uiter->first].detectorID];
            double v_min = 2*x_pos*u_costheta[sID] - u_pos - v_win;
            double v_max = v_min + 2.*v_win;

#ifdef _DEBUG_ON
            LogInfo("V plane window:" << v_min << "  " << v_max);
#endif
            for(std::list<SRawEvent::hit_pair>::iterator viter = pairs_V.begin(); viter != pairs_V.end(); ++viter)
            {
                double v_pos = viter->second >= 0 ? 0.5*(hitAll[viter->first].pos + hitAll[viter->second].pos) : hitAll[viter->first].pos;
#ifdef _DEBUG_ON
                LogInfo("Trying V hits " << viter->first << "  " << viter->second << "  " << hitAll[viter->first].elementID << " at " << v_pos);
#endif
                if(v_pos < v_min || v_pos > v_max) continue;

                //Now add the tracklet
                int LR1 = 0;
                int LR2 = 0;
                Tracklet tracklet_new;
                tracklet_new.stationID = stationID;

                //resolveLeftRight(*xiter, LR1, LR2);
                if(xiter->first >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[xiter->first], LR1));
                    tracklet_new.nXHits++;
                }
                if(xiter->second >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[xiter->second], LR2));
                    tracklet_new.nXHits++;
                }
		if(!OLD_TRACKING){
		  tracklet_new.getSlopesX(hitAll[xiter->first], hitAll[xiter->second]); //Here, we find the four possible X-Z lines
		}
                //resolveLeftRight(*uiter, LR1, LR2);
                if(uiter->first >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[uiter->first], LR1));
                    tracklet_new.nUHits++;
                }
                if(uiter->second >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[uiter->second], LR2));
                    tracklet_new.nUHits++;
                }
		if(!OLD_TRACKING){
		  tracklet_new.getSlopesU(hitAll[uiter->first], hitAll[uiter->second]); //find the four possible U-Z lines
		}

                //resolveLeftRight(*viter, LR1, LR2);
                if(viter->first >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[viter->first], LR1));
                    tracklet_new.nVHits++;
                }
                if(viter->second >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[viter->second], LR2));
                    tracklet_new.nVHits++;
                }
		if(!OLD_TRACKING){
		  tracklet_new.getSlopesV(hitAll[viter->first], hitAll[viter->second]); //find the four possible V-Z lines
		}

                tracklet_new.sortHits();
                if(tracklet_new.isValid() == 0) //TODO: What IS THIS?
                {
		  fitTracklet(tracklet_new); //This is where the original DCA minimization is performed
                }
                else
                {
                    continue;
                }

#ifdef _DEBUG_ON
                tracklet_new.print();
#endif
                if(acceptTracklet(tracklet_new))
                {
                    trackletsInSt[listID].push_back(tracklet_new);
                }
#ifdef _DEBUG_ON
                else
                {
                    LogInfo("Rejected!!!");
                }
#endif
            }
        }
    }

    //Reduce the tracklet list and add dummy hits
    //reduceTrackletList(trackletsInSt[listID]);
    for(std::list<Tracklet>::iterator iter = trackletsInSt[listID].begin(); iter != trackletsInSt[listID].end(); ++iter)
    {
        iter->addDummyHits();
    }

    //Only retain the best 200 tracklets if exceeded
    std::cout<<"NUMBER of tracklets before resize = "<<trackletsInSt[listID].size()<<std::endl; //WPM
    if(trackletsInSt[listID].size() > 1000)
    {
        trackletsInSt[listID].sort();
        trackletsInSt[listID].resize(1000);
    }
}


void KalmanFastTracking::buildTrackletsInStationSlim(int stationID, int listID, double* pos_exp, double* window)
{
#ifdef _DEBUG_ON
    LogInfo("Building tracklets in station (slim version) " << stationID);
#endif

    //actuall ID of the tracklet lists
    int sID = stationID - 1;

    //Extract the X, U, V hit pairs
    std::list<SRawEvent::hit_pair> pairs_X, pairs_U, pairs_V;
    if(pos_exp == nullptr)
    {
        pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0]);
        //pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1]);
        //pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2]);
    }
    else
    {
        //Note that in pos_exp[], index 0 stands for X, index 1 stands for U, index 2 stands for V
        pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0], pos_exp[0], window[0]);
        //pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1], pos_exp[1], window[1]);
        //pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2], pos_exp[2], window[2]);
    }

#ifdef _DEBUG_ON
    LogInfo("Hit pairs in this event: ");
    for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_X.begin(); iter != pairs_X.end(); ++iter) LogInfo("X :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
    //for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_U.begin(); iter != pairs_U.end(); ++iter) LogInfo("U :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
    //for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_V.begin(); iter != pairs_V.end(); ++iter) LogInfo("V :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
#endif

    //if(pairs_X.empty() || pairs_U.empty() || pairs_V.empty())
    if(pairs_X.empty())
    {
#ifdef _DEBUG_ON
        LogInfo("Not all view has hits in station " << stationID);
#endif
        return;
    }

    //X-U combination first, then add V pairs
    for(std::list<SRawEvent::hit_pair>::iterator xiter = pairs_X.begin(); xiter != pairs_X.end(); ++xiter)
    {

      int LR1 = 0;
      int LR2 = 0;
      Tracklet tracklet_new;
      tracklet_new.stationID = stationID;
      
      //resolveLeftRight(*xiter, LR1, LR2);
      if(xiter->first >= 0)
	{
	  tracklet_new.hits.push_back(SignedHit(hitAll[xiter->first], LR1));
	  tracklet_new.nXHits++;
	}
      if(xiter->second >= 0)
	{
	  tracklet_new.hits.push_back(SignedHit(hitAll[xiter->second], LR2));
	  tracklet_new.nXHits++;
	}
      if(!OLD_TRACKING){
	tracklet_new.getSlopesX(hitAll[xiter->first], hitAll[xiter->second]); //Here, we find the four possible X-Z lines
      }
      
      tracklet_new.sortHits();
#ifdef _DEBUG_ON
      std::cout<<"About to print new tracklet"<<std::endl;
      tracklet_new.print();
      std::cout<<"How many hits does it have?! "<<tracklet_new.hits.size()<<std::endl;
#endif
      
      if(tracklet_new.hits.size() == 2){
	trackletsInStSlimX[listID].push_back(tracklet_new);
      }
      /*     
        //U projections from X plane
        double x_pos = xiter->second >= 0 ? 0.5*(hitAll[xiter->first].pos + hitAll[xiter->second].pos) : hitAll[xiter->first].pos;
        double u_min = x_pos*u_costheta[sID] - u_win[sID];
        double u_max = u_min + 2.*u_win[sID];

#ifdef _DEBUG_ON
        LogInfo("Trying X hits " << xiter->first << "  " << xiter->second << "  " << hitAll[xiter->first].elementID << " at " << x_pos);
        LogInfo("U plane window:" << u_min << "  " << u_max);
#endif
        for(std::list<SRawEvent::hit_pair>::iterator uiter = pairs_U.begin(); uiter != pairs_U.end(); ++uiter)
        {
            double u_pos = uiter->second >= 0 ? 0.5*(hitAll[uiter->first].pos + hitAll[uiter->second].pos) : hitAll[uiter->first].pos;
#ifdef _DEBUG_ON
            LogInfo("Trying U hits " << uiter->first << "  " << uiter->second << "  " << hitAll[uiter->first].elementID << " at " << u_pos);
#endif
            if(u_pos < u_min || u_pos > u_max) continue;

            //V projections from X and U plane
            double z_x = xiter->second >= 0 ? z_plane_x[sID] : z_plane[hitAll[xiter->first].detectorID];
            double z_u = uiter->second >= 0 ? z_plane_u[sID] : z_plane[hitAll[uiter->first].detectorID];
            double z_v = z_plane_v[sID];
            double v_win1 = spacing_plane[hitAll[uiter->first].detectorID]*2.*u_costheta[sID];
            double v_win2 = fabs((z_u + z_v - 2.*z_x)*u_costheta[sID]*TX_MAX);
            double v_win3 = fabs((z_v - z_u)*u_sintheta[sID]*TY_MAX);
            double v_win = v_win1 + v_win2 + v_win3 + 2.*spacing_plane[hitAll[uiter->first].detectorID];
            double v_min = 2*x_pos*u_costheta[sID] - u_pos - v_win;
            double v_max = v_min + 2.*v_win;

#ifdef _DEBUG_ON
            LogInfo("V plane window:" << v_min << "  " << v_max);
#endif
            for(std::list<SRawEvent::hit_pair>::iterator viter = pairs_V.begin(); viter != pairs_V.end(); ++viter)
            {
                double v_pos = viter->second >= 0 ? 0.5*(hitAll[viter->first].pos + hitAll[viter->second].pos) : hitAll[viter->first].pos;
#ifdef _DEBUG_ON
                LogInfo("Trying V hits " << viter->first << "  " << viter->second << "  " << hitAll[viter->first].elementID << " at " << v_pos);
#endif
                if(v_pos < v_min || v_pos > v_max) continue;

                //Now add the tracklet
                int LR1 = 0;
                int LR2 = 0;
                Tracklet tracklet_new;
                tracklet_new.stationID = stationID;

                //resolveLeftRight(*xiter, LR1, LR2);
                if(xiter->first >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[xiter->first], LR1));
                    tracklet_new.nXHits++;
                }
                if(xiter->second >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[xiter->second], LR2));
                    tracklet_new.nXHits++;
                }
		if(!OLD_TRACKING){
		  tracklet_new.getSlopesX(hitAll[xiter->first], hitAll[xiter->second]); //Here, we find the four possible X-Z lines
		}
                //resolveLeftRight(*uiter, LR1, LR2);
                if(uiter->first >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[uiter->first], LR1));
                    tracklet_new.nUHits++;
                }
                if(uiter->second >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[uiter->second], LR2));
                    tracklet_new.nUHits++;
                }
		if(!OLD_TRACKING){
		  tracklet_new.getSlopesU(hitAll[uiter->first], hitAll[uiter->second]); //find the four possible U-Z lines
		}

                //resolveLeftRight(*viter, LR1, LR2);
                if(viter->first >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[viter->first], LR1));
                    tracklet_new.nVHits++;
                }
                if(viter->second >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[viter->second], LR2));
                    tracklet_new.nVHits++;
                }
		if(!OLD_TRACKING){
		  tracklet_new.getSlopesV(hitAll[viter->first], hitAll[viter->second]); //find the four possible V-Z lines
		}

                tracklet_new.sortHits();
                if(tracklet_new.isValid() == 0) //TODO: What IS THIS?
                {
		  fitTracklet(tracklet_new); //This is where the original DCA minimization is performed
                }
                else
                {
                    continue;
                }

#ifdef _DEBUG_ON
                tracklet_new.print();
#endif
                if(acceptTracklet(tracklet_new))
                {
                    trackletsInSt[listID].push_back(tracklet_new);
                }
#ifdef _DEBUG_ON
                else
                {
                    LogInfo("Rejected!!!");
                }
#endif
            }
        }*/
    }

    std::cout<<"Number of x pairs in station "<<listID<<" is "<<trackletsInStSlimX[listID].size()<<std::endl;
    /*    //Reduce the tracklet list and add dummy hits
    //reduceTrackletList(trackletsInSt[listID]);
    for(std::list<Tracklet>::iterator iter = trackletsInSt[listID].begin(); iter != trackletsInSt[listID].end(); ++iter)
    {
        iter->addDummyHits();
    }

    //Only retain the best 200 tracklets if exceeded
    std::cout<<"NUMBER of tracklets before resize = "<<trackletsInSt[listID].size()<<std::endl; //WPM
    if(trackletsInSt[listID].size() > 1000)
    {
        trackletsInSt[listID].sort();
        trackletsInSt[listID].resize(1000);
    }*/
}


void KalmanFastTracking::buildTrackletsInStationSlimU(int stationID, int listID, double* pos_exp, double* window)
{
#ifdef _DEBUG_ON
    LogInfo("Building U tracklets in station (slim version) " << stationID);
#endif

    //actuall ID of the tracklet lists
    int sID = stationID - 1;

    //Extract the X, U, V hit pairs
    std::list<SRawEvent::hit_pair> pairs_X, pairs_U, pairs_V;
    if(pos_exp == nullptr)
    {
      //pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0]);
        pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1]);
        //pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2]);
    }
    else
    {
        //Note that in pos_exp[], index 0 stands for X, index 1 stands for U, index 2 stands for V
        //pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0], pos_exp[0], window[0]);
        pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1], pos_exp[1], window[1]);
        //pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2], pos_exp[2], window[2]);
    }

#ifdef _DEBUG_ON
    LogInfo("Hit pairs in this event: ");
    //for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_X.begin(); iter != pairs_X.end(); ++iter) LogInfo("X :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
    for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_U.begin(); iter != pairs_U.end(); ++iter) LogInfo("U :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
    //for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_V.begin(); iter != pairs_V.end(); ++iter) LogInfo("V :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
#endif

    //if(pairs_X.empty() || pairs_U.empty() || pairs_V.empty())
    if(pairs_U.empty())
    {
#ifdef _DEBUG_ON
        LogInfo("Not all view has hits in station " << stationID);
#endif
        return;
    }

    //X-U combination first, then add V pairs
    for(std::list<SRawEvent::hit_pair>::iterator xiter = pairs_U.begin(); xiter != pairs_U.end(); ++xiter)
    {

      int LR1 = 0;
      int LR2 = 0;
      Tracklet tracklet_new;
      tracklet_new.stationID = stationID;
      
      //resolveLeftRight(*xiter, LR1, LR2);
      if(xiter->first >= 0)
	{
	  tracklet_new.hits.push_back(SignedHit(hitAll[xiter->first], LR1));
	  tracklet_new.nUHits++;
	}
      if(xiter->second >= 0)
	{
	  tracklet_new.hits.push_back(SignedHit(hitAll[xiter->second], LR2));
	  tracklet_new.nUHits++;
	}
      if(!OLD_TRACKING){
	tracklet_new.getSlopesU(hitAll[xiter->first], hitAll[xiter->second]); //Here, we find the four possible X-Z lines
      }
      
      tracklet_new.sortHits();
#ifdef _DEBUG_ON
      std::cout<<"About to print new U tracklet"<<std::endl;
      tracklet_new.print();
      std::cout<<"How many hits does it have?! "<<tracklet_new.hits.size()<<std::endl;
#endif
      
      if(tracklet_new.hits.size() == 2){
	trackletsInStSlimU[listID].push_back(tracklet_new);
      }
    }

    std::cout<<"Number of U pairs in station "<<listID<<" is "<<trackletsInStSlimU[listID].size()<<std::endl;

}



void KalmanFastTracking::buildTrackletsInStationSlimV(int stationID, int listID, double* pos_exp, double* window)
{
#ifdef _DEBUG_ON
    LogInfo("Building V tracklets in station (slim version) " << stationID);
#endif

    //actuall ID of the tracklet lists
    int sID = stationID - 1;

    //Extract the X, V, V hit pairs
    std::list<SRawEvent::hit_pair> pairs_X, pairs_U, pairs_V;
    if(pos_exp == nullptr)
    {
      //pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0]);
      //pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1]);
        pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2]);
    }
    else
    {
        //Note that in pos_exp[], index 0 stands for X, index 1 stands for U, index 2 stands for V
        //pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0], pos_exp[0], window[0]);
        //pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1], pos_exp[1], window[1]);
        pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2], pos_exp[2], window[2]);
    }

#ifdef _DEBUG_ON
    LogInfo("Hit pairs in this event: ");
    //for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_X.begin(); iter != pairs_X.end(); ++iter) LogInfo("X :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
    //for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_U.begin(); iter != pairs_U.end(); ++iter) LogInfo("U :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
    for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_V.begin(); iter != pairs_V.end(); ++iter) LogInfo("V :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
#endif

    //if(pairs_X.empty() || pairs_U.empty() || pairs_V.empty())
    if(pairs_V.empty())
    {
#ifdef _DEBUG_ON
        LogInfo("Not all view has hits in station " << stationID);
#endif
        return;
    }

    //X-V combination first, then add V pairs
    for(std::list<SRawEvent::hit_pair>::iterator xiter = pairs_V.begin(); xiter != pairs_V.end(); ++xiter)
    {

      int LR1 = 0;
      int LR2 = 0;
      Tracklet tracklet_new;
      tracklet_new.stationID = stationID;
      
      //resolveLeftRight(*xiter, LR1, LR2);
      if(xiter->first >= 0)
	{
	  tracklet_new.hits.push_back(SignedHit(hitAll[xiter->first], LR1));
	  tracklet_new.nVHits++;
	}
      if(xiter->second >= 0)
	{
	  tracklet_new.hits.push_back(SignedHit(hitAll[xiter->second], LR2));
	  tracklet_new.nVHits++;
	}
      if(!OLD_TRACKING){
	tracklet_new.getSlopesV(hitAll[xiter->first], hitAll[xiter->second]); //Here, we find the four possible X-Z lines
      }
      
      tracklet_new.sortHits();
#ifdef _DEBUG_ON
      std::cout<<"About to print new V tracklet"<<std::endl;
      tracklet_new.print();
      std::cout<<"How many hits does it have?! "<<tracklet_new.hits.size()<<std::endl;
#endif
      
      if(tracklet_new.hits.size() == 2){
	trackletsInStSlimV[listID].push_back(tracklet_new);
      }
    }

    std::cout<<"Number of V pairs in station "<<listID<<" is "<<trackletsInStSlimV[listID].size()<<std::endl;

}



void KalmanFastTracking::buildTrackletsInStationWithUV(int stationID, int listID, Tracklet& tracklet23, double* pos_exp, double* window)
{

  Tracklet tracklet_best;
  
#ifdef _DEBUG_ON
    LogInfo("Building tracklets with UV in station " << stationID);
#endif

    std::cout<<"tracklet23 first hit detID = "<<tracklet23.getHit(0).hit.detectorID<<std::endl;
    std::cout<<"tracklet23 third hit detID = "<<tracklet23.getHit(2).hit.detectorID<<std::endl;

    //actuall ID of the tracklet lists
    int sID = stationID - 1;
    int sID2 = 2;
    int sID3;
    
    //Extract the X, U, V hit pairs
    std::list<SRawEvent::hit_pair> pairs_U2, pairs_V2, pairs_U3, pairs_V3;
    if(pos_exp == nullptr)
    {
      //pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0]);
        pairs_U2 = rawEvent->getPartialHitPairsInSuperDetector(superIDs[2][1]);
        pairs_V2 = rawEvent->getPartialHitPairsInSuperDetector(superIDs[2][2]);
    }
    else
    {
        //Note that in pos_exp[], index 0 stands for X, index 1 stands for U, index 2 stands for V
        //pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0], pos_exp[0], window[0]);
        pairs_U2 = rawEvent->getPartialHitPairsInSuperDetector(superIDs[2][1], pos_exp[1], window[1]);
        pairs_V2 = rawEvent->getPartialHitPairsInSuperDetector(superIDs[2][2], pos_exp[2], window[2]);
    }

    if(tracklet23.getHit(2).hit.detectorID > 18 && tracklet23.getHit(2).hit.detectorID < 25){
      sID3 = 3;
      if(pos_exp == nullptr)
	{
	  //pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0]);
	  pairs_U3 = rawEvent->getPartialHitPairsInSuperDetector(superIDs[3][1]);
	  pairs_V3 = rawEvent->getPartialHitPairsInSuperDetector(superIDs[3][2]);
	}
      else
	{
	  //Note that in pos_exp[], index 0 stands for X, index 1 stands for U, index 2 stands for V
	  //pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0], pos_exp[0], window[0]);
	  pairs_U3 = rawEvent->getPartialHitPairsInSuperDetector(superIDs[3][1], pos_exp[1], window[1]);
	  pairs_V3 = rawEvent->getPartialHitPairsInSuperDetector(superIDs[3][2], pos_exp[2], window[2]);
	} 
    }
    
    if(tracklet23.getHit(2).hit.detectorID > 24 && tracklet23.getHit(2).hit.detectorID < 31){
      sID3 = 4;
      if(pos_exp == nullptr)
	{
	  //pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0]);
	  pairs_U3 = rawEvent->getPartialHitPairsInSuperDetector(superIDs[4][1]);
	  pairs_V3 = rawEvent->getPartialHitPairsInSuperDetector(superIDs[4][2]);
	}
      else
	{
	  //Note that in pos_exp[], index 0 stands for X, index 1 stands for U, index 2 stands for V
	  //pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0], pos_exp[0], window[0]);
	  pairs_U3 = rawEvent->getPartialHitPairsInSuperDetector(superIDs[4][1], pos_exp[1], window[1]);
	  pairs_V3 = rawEvent->getPartialHitPairsInSuperDetector(superIDs[4][2], pos_exp[2], window[2]);
	} 
    }

#ifdef _DEBUG_ON
    LogInfo("Building withUV tracklets.  Hit pairs in this event: ");
    //for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_X.begin(); iter != pairs_X.end(); ++iter) LogInfo("X :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
    for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_U2.begin(); iter != pairs_U2.end(); ++iter) LogInfo("U2 :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
    for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_V2.begin(); iter != pairs_V2.end(); ++iter) LogInfo("V2 :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
    for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_U3.begin(); iter != pairs_U3.end(); ++iter) LogInfo("U3 :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
    for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_V3.begin(); iter != pairs_V3.end(); ++iter) LogInfo("V3 :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
#endif

    if(pairs_U2.empty() || pairs_V2.empty() || pairs_U3.empty() || pairs_V3.empty())
    {
#ifdef _DEBUG_ON
        LogInfo("Not all view has hits in station " << stationID);
#endif
        return;
    }

    //X-U combination first, then add V pairs
    //U projections from X plane
    double x_pos2 = 0.5*(tracklet23.getHit(0).hit.pos + tracklet23.getHit(1).hit.pos);
    double x_pos3 = 0.5*(tracklet23.getHit(2).hit.pos + tracklet23.getHit(3).hit.pos);
    //double x_pos = xiter->second >= 0 ? 0.5*(hitAll[xiter->first].pos + hitAll[xiter->second].pos) : hitAll[xiter->first].pos;
    double u_min2 = x_pos2*u_costheta[sID2] - u_win[sID2];
    double u_max2 = u_min2 + 2.*u_win[sID2];

#ifdef _DEBUG_ON
    //LogInfo("Trying X hits " << xiter->first << "  " << xiter->second << "  " << hitAll[xiter->first].elementID << " at " << x_pos);
        LogInfo("U2 plane window:" << u_min2 << "  " << u_max2);
#endif
        for(std::list<SRawEvent::hit_pair>::iterator uiter2 = pairs_U2.begin(); uiter2 != pairs_U2.end(); ++uiter2)
        {
            double u_pos2 = uiter2->second >= 0 ? 0.5*(hitAll[uiter2->first].pos + hitAll[uiter2->second].pos) : hitAll[uiter2->first].pos;
#ifdef _DEBUG_ON
            LogInfo("Trying U2 hits " << uiter2->first << "  " << uiter2->second << "  " << hitAll[uiter2->first].elementID << " at " << u_pos2);
#endif
            if(u_pos2 < u_min2 || u_pos2 > u_max2) continue;

            //V projections from X and U plane
            //double z_x2 = xiter->second >= 0 ? z_plane_x[sID] : z_plane[hitAll[xiter->first].detectorID];
	    double z_x2 = z_plane_x[sID2];
            double z_u2 = uiter2->second >= 0 ? z_plane_u[sID2] : z_plane[hitAll[uiter2->first].detectorID];
            double z_v2 = z_plane_v[sID2];
            double v_win1_2 = spacing_plane[hitAll[uiter2->first].detectorID]*2.*u_costheta[sID2];
            double v_win2_2 = fabs((z_u2 + z_v2 - 2.*z_x2)*u_costheta[sID2]*TX_MAX);
            double v_win3_2 = fabs((z_v2 - z_u2)*u_sintheta[sID2]*TY_MAX);
            double v_win_2 = v_win1_2 + v_win2_2 + v_win3_2 + 2.*spacing_plane[hitAll[uiter2->first].detectorID];
            double v_min_2 = 2*x_pos2*u_costheta[sID2] - u_pos2 - v_win_2;
            double v_max_2 = v_min_2 + 2.*v_win_2;

#ifdef _DEBUG_ON
            LogInfo("V2 plane window:" << v_min_2 << "  " << v_max_2);
#endif
            for(std::list<SRawEvent::hit_pair>::iterator viter2 = pairs_V2.begin(); viter2 != pairs_V2.end(); ++viter2)
            {
                double v_pos2 = viter2->second >= 0 ? 0.5*(hitAll[viter2->first].pos + hitAll[viter2->second].pos) : hitAll[viter2->first].pos;
#ifdef _DEBUG_ON
                LogInfo("Trying V2 hits " << viter2->first << "  " << viter2->second << "  " << hitAll[viter2->first].elementID << " at " << v_pos2);
#endif
                if(v_pos2 < v_min_2 || v_pos2 > v_max_2) continue;

                //Now add the tracklet
                int LR1 = 0;
                int LR2 = 0;
                Tracklet tracklet_new2;
                tracklet_new2.stationID = sID2+1;

                //resolveLeftRight(*xiter, LR1, LR2);
		tracklet_new2.hits.push_back(SignedHit(tracklet23.getHit(0).hit, LR1));
		tracklet_new2.hits.push_back(SignedHit(tracklet23.getHit(1).hit, LR2));
		if(!OLD_TRACKING){
		  tracklet_new2.getSlopesX(tracklet23.getHit(0).hit, tracklet23.getHit(1).hit); //Here, we find the four possible X-Z lines
		}
                //resolveLeftRight(*uiter, LR1, LR2);
                if(uiter2->first >= 0)
                {
                    tracklet_new2.hits.push_back(SignedHit(hitAll[uiter2->first], LR1));
                    tracklet_new2.nUHits++;
                }
                if(uiter2->second >= 0)
                {
                    tracklet_new2.hits.push_back(SignedHit(hitAll[uiter2->second], LR2));
                    tracklet_new2.nUHits++;
                }
		if(!OLD_TRACKING){
		  tracklet_new2.getSlopesU(hitAll[uiter2->first], hitAll[uiter2->second]); //find the four possible U-Z lines
		}

                //resolveLeftRight(*viter, LR1, LR2);
                if(viter2->first >= 0)
                {
                    tracklet_new2.hits.push_back(SignedHit(hitAll[viter2->first], LR1));
                    tracklet_new2.nVHits++;
                }
                if(viter2->second >= 0)
                {
                    tracklet_new2.hits.push_back(SignedHit(hitAll[viter2->second], LR2));
                    tracklet_new2.nVHits++;
                }
		if(!OLD_TRACKING){
		  tracklet_new2.getSlopesV(hitAll[viter2->first], hitAll[viter2->second]); //find the four possible V-Z lines
		}

                tracklet_new2.sortHits();
                //if(!(tracklet_new2.isValid() == 0)) //TODO: What IS THIS?
                //{
		//  continue;
		//}
		/*
		{
		  //fitTracklet(tracklet_new); //This is where the original DCA minimization is performed
                }
                else
                {
                    continue;
		    }*/

#ifdef _DEBUG_ON
		std::cout<<"OK HERE'S A STATION 2 TRACKLET:"<<std::endl;
                tracklet_new2.print();
#endif

    double u_min3 = x_pos3*u_costheta[sID3] - u_win[sID3];
    double u_max3 = u_min3 + 2.*u_win[sID3];

#ifdef _DEBUG_ON
    //LogInfo("Trying X hits " << xiter->first << "  " << xiter->second << "  " << hitAll[xiter->first].elementID << " at " << x_pos);
        LogInfo("U3 plane window:" << u_min3 << "  " << u_max3);
#endif
        for(std::list<SRawEvent::hit_pair>::iterator uiter3 = pairs_U3.begin(); uiter3 != pairs_U3.end(); ++uiter3)
        {
            double u_pos3 = uiter3->second >= 0 ? 0.5*(hitAll[uiter3->first].pos + hitAll[uiter3->second].pos) : hitAll[uiter3->first].pos;
#ifdef _DEBUG_ON
            LogInfo("Trying U3 hits " << uiter3->first << "  " << uiter3->second << "  " << hitAll[uiter3->first].elementID << " at " << u_pos3);
#endif
            if(u_pos3 < u_min3 || u_pos3 > u_max3) continue;

            //V projections from X and U plane
            //double z_x3 = xiter->second >= 0 ? z_plane_x[sID] : z_plane[hitAll[xiter->first].detectorID];
	    double z_x3 = z_plane_x[sID3];
            double z_u3 = uiter3->second >= 0 ? z_plane_u[sID3] : z_plane[hitAll[uiter3->first].detectorID];
            double z_v3 = z_plane_v[sID3];
            double v_win1_3 = spacing_plane[hitAll[uiter3->first].detectorID]*2.*u_costheta[sID3];
            double v_win2_3 = fabs((z_u3 + z_v3 - 2.*z_x3)*u_costheta[sID3]*TX_MAX);
            double v_win3_3 = fabs((z_v3 - z_u3)*u_sintheta[sID3]*TY_MAX);
            double v_win_3 = v_win1_3 + v_win2_3 + v_win3_3 + 2.*spacing_plane[hitAll[uiter3->first].detectorID];
            double v_min_3 = 2*x_pos3*u_costheta[sID3] - u_pos3 - v_win_3;
            double v_max_3 = v_min_3 + 2.*v_win_3;

#ifdef _DEBUG_ON
            LogInfo("V3 plane window:" << v_min_3 << "  " << v_max_3);
#endif
            for(std::list<SRawEvent::hit_pair>::iterator viter3 = pairs_V3.begin(); viter3 != pairs_V3.end(); ++viter3)
            {
                double v_pos3 = viter3->second >= 0 ? 0.5*(hitAll[viter3->first].pos + hitAll[viter3->second].pos) : hitAll[viter3->first].pos;
#ifdef _DEBUG_ON
                LogInfo("Trying V3 hits " << viter3->first << "  " << viter3->second << "  " << hitAll[viter3->first].elementID << " at " << v_pos3);
#endif
                if(v_pos3 < v_min_3 || v_pos3 > v_max_3) continue;

                //Now add the tracklet
                //int LR1 = 0;
                //int LR2 = 0;
                Tracklet tracklet_new3;
                tracklet_new3.stationID = sID3+1;

                //resolveLeftRight(*xiter, LR1, LR2);
		tracklet_new3.hits.push_back(SignedHit(tracklet23.getHit(2).hit,LR1));
		tracklet_new3.hits.push_back(SignedHit(tracklet23.getHit(3).hit,LR2));
		if(!OLD_TRACKING){
		  tracklet_new3.getSlopesX(tracklet23.getHit(2).hit, tracklet23.getHit(3).hit); //Here, we find the four possible X-Z lines
		}
                //resolveLeftRight(*uiter, LR1, LR2);
                if(uiter3->first >= 0)
                {
                    tracklet_new3.hits.push_back(SignedHit(hitAll[uiter3->first], LR1));
                    tracklet_new3.nUHits++;
                }
                if(uiter3->second >= 0)
                {
                    tracklet_new3.hits.push_back(SignedHit(hitAll[uiter3->second], LR2));
                    tracklet_new3.nUHits++;
                }
		if(!OLD_TRACKING){
		  tracklet_new3.getSlopesU(hitAll[uiter3->first], hitAll[uiter3->second]); //find the four possible U-Z lines
		}

                //resolveLeftRight(*viter, LR1, LR2);
                if(viter3->first >= 0)
                {
                    tracklet_new3.hits.push_back(SignedHit(hitAll[viter3->first], LR1));
                    tracklet_new3.nVHits++;
                }
                if(viter3->second >= 0)
                {
                    tracklet_new3.hits.push_back(SignedHit(hitAll[viter3->second], LR2));
                    tracklet_new3.nVHits++;
                }
		if(!OLD_TRACKING){
		  tracklet_new3.getSlopesV(hitAll[viter3->first], hitAll[viter3->second]); //find the four possible V-Z lines
		}

                tracklet_new3.sortHits();
                //if(!(tracklet_new3.isValid() == 0)) //TODO: What IS THIS?
		//{
		//  continue;
		//}

#ifdef _DEBUG_ON
		std::cout<<"OK HERE'S A STATION 3 TRACKLET:"<<std::endl;
                tracklet_new3.print();
#endif

		
		Tracklet tracklet_new_23;
	      if(compareTracklets(tracklet_new2, tracklet_new3)){
		std::cout<<"THE COMBINED TRACKLET PASSED"<<std::endl;
		tracklet_new_23 = (tracklet_new2) + (tracklet_new3);
		tracklet_new_23.tx = tracklet_new2.tx; //This is needed to "seed" the tracklet fit that happens below.  This tx and ty information is assigned in compareTracklets below
		tracklet_new_23.ty = tracklet_new2.ty;

		tracklet_new_23.st2Z = tracklet_new2.st2Z;
		tracklet_new_23.st2X = tracklet_new2.st2X;
		tracklet_new_23.st2Y = tracklet_new2.st2Y;
		tracklet_new_23.st3Y = tracklet_new2.st3Y;
		tracklet_new_23.st2U = tracklet_new2.st2U;
		tracklet_new_23.st2V = tracklet_new2.st2V;
		tracklet_new_23.st2Usl = tracklet_new2.st2Usl;
		tracklet_new_23.st2Vsl = tracklet_new2.st2Vsl;
		fitTracklet(tracklet_new_23); //This is the fit that needs the seeded tx and ty information. Without the seed information, the fit occasionally finds bad slope and X0 or Y0 values, much in the same way that it does for single-station tracklets.  Note from Patrick: this fit could probably be throw out, as we already know the tx and ty information from compareTracklets.  I would just need to extrapolate back to z = 0 and calculate the chisq by hand
	      }
	      else{
		continue;
	      }

            if(tracklet_new_23.chisq > 9000.)
            {
#ifdef _DEBUG_ON
                tracklet_new_23.print();
                LogInfo("Impossible combination!");
#endif
                continue;
            }

		
                if(acceptTracklet(tracklet_new_23))
                {
		  if(tracklet_new_23 < tracklet_best){
		    tracklet_best = tracklet_new_23;
		  }
		  //trackletsInSt[listID].push_back(tracklet_new_23);
                }
#ifdef _DEBUG_ON
                else
                {
                    LogInfo("Rejected!!!");
                }
#endif
            }
        }
	    }
	}

	if(acceptTracklet(tracklet_best)){
	  trackletsInSt[listID].push_back(tracklet_best);
	}
    //Reduce the tracklet list and add dummy hits
    //reduceTrackletList(trackletsInSt[listID]);
    for(std::list<Tracklet>::iterator iter = trackletsInSt[listID].begin(); iter != trackletsInSt[listID].end(); ++iter)
    {
        iter->addDummyHits();
    }

    //Only retain the best 200 tracklets if exceeded
    std::cout<<"NUMBER of tracklets before resize = "<<trackletsInSt[listID].size()<<std::endl; //WPM
    if(trackletsInSt[listID].size() > 1000)
    {
        trackletsInSt[listID].sort();
        trackletsInSt[listID].resize(1000);
    }
}


bool KalmanFastTracking::acceptTracklet(Tracklet& tracklet)
{
    //Tracklet itself is okay with enough hits (4-out-of-6) and small chi square
    if(tracklet.isValid() == 0)
    {
#ifdef _DEBUG_ON
        LogInfo("Failed in quality check!");
#endif
        return false;
    }

    if(COARSE_MODE) return true;

    //Hodoscope masking requirement
    if(!hodoMask(tracklet)) return false;

    //For back partials, require projection inside KMAG, and muon id in prop. tubes
    if(tracklet.stationID > nStations-2)
    {
      if(!COSMIC_MODE && !p_geomSvc->isInKMAG(tracklet.getExpPositionX(Z_KMAG_BEND), tracklet.getExpPositionY(Z_KMAG_BEND))) return false;
      if(!TRACK_ELECTRONS && !(muonID_comp(tracklet) || muonID_search(tracklet))) return false; //Muon check for 2-3 connected tracklets.  This needs to be off for electron tracks
      if(TRACK_ELECTRONS && !(muonID_comp(tracklet) || muonID_search(tracklet) || tracklet.stationID > 5)){
	return false;
      }
    }

    //If everything is fine ...
#ifdef _DEBUG_ON
    LogInfo("AcceptTracklet!!!");
#endif
    return true;
}

bool KalmanFastTracking::hodoMask(Tracklet& tracklet)
{
    //LogInfo(tracklet.stationID);
  if(TRACK_ELECTRONS && (tracklet.stationID == 4 || tracklet.stationID == 5)) return true; //Patrick's skip of hodoscope checks for station 3 tracks in the electron-tracking setup.  I could actually probably extrapolate backwards the station 2 hodoscope, now that I get an accurate X-Z slope in station 3
  int nHodoHits = 0;
  if(!OLD_TRACKING){
    //Performing a valid extrapolation in X-Z for station 2 tracklets.  This should be improved.  Currently carries around the old fudge factor
    if(tracklet.stationID == 3){
      for(std::vector<int>::iterator stationID = stationIDs_mask[tracklet.stationID-1].begin(); stationID != stationIDs_mask[tracklet.stationID-1].end(); ++stationID){
	bool masked = false;
	for(std::list<int>::iterator iter = hitIDs_mask[*stationID-1].begin(); iter != hitIDs_mask[*stationID-1].end(); ++iter){
	  int detectorID = hitAll[*iter].detectorID;
	  int elementID = hitAll[*iter].elementID;
	  
	  int idx1 = detectorID - nChamberPlanes - 1;
	  int idx2 = elementID - 1;
	  
	  double factor = tracklet.stationID == nChamberPlanes/6-2 ? 5. : 3.;   //special for station-2, based on real data tuning
	  double xfudge = tracklet.stationID < nStations-1 ? 0.5*(x_mask_max[idx1][idx2] - x_mask_min[idx1][idx2]) : 0.15*(x_mask_max[idx1][idx2] - x_mask_min[idx1][idx2]);
	  double z_hodo = z_mask[idx1];
	  
	  for(unsigned int pl = 0; pl < tracklet.possibleXLines.size(); pl++){
	    double extrapolation = tracklet.possibleXLines.at(pl).slopeX*(z_hodo - tracklet.possibleXLines.at(pl).initialZ) + tracklet.possibleXLines.at(pl).initialX;
	    double err_x = std::abs(factor*extrapolation + xfudge);
	    double x_min = x_mask_min[idx1][idx2] - err_x;
	    double x_max = x_mask_max[idx1][idx2] + err_x;
	    if(extrapolation > x_min && extrapolation < x_max){
	      masked = true;
	      break;
	    }
	  }
	}
	if(!masked) return false;
      }
    }
  }

  if(tracklet.stationID > 5){
    for(std::vector<int>::iterator stationID = stationIDs_mask[tracklet.stationID-1].begin(); stationID != stationIDs_mask[tracklet.stationID-1].end(); ++stationID)
    {
        bool masked = false;
        for(std::list<int>::iterator iter = hitIDs_mask[*stationID-1].begin(); iter != hitIDs_mask[*stationID-1].end(); ++iter)
        {
            int detectorID = hitAll[*iter].detectorID;
            int elementID = hitAll[*iter].elementID;

            int idx1 = detectorID - nChamberPlanes - 1;
            int idx2 = elementID - 1;

            double factor = tracklet.stationID == nChamberPlanes/6-2 ? 5. : 3.;   //special for station-2, based on real data tuning
            double xfudge = tracklet.stationID < nStations-1 ? 0.5*(x_mask_max[idx1][idx2] - x_mask_min[idx1][idx2]) : 0.15*(x_mask_max[idx1][idx2] - x_mask_min[idx1][idx2]);
            double z_hodo = z_mask[idx1];
            double x_hodo = tracklet.getExpPositionX(z_hodo);
            double y_hodo = tracklet.getExpPositionY(z_hodo);
            double err_x = factor*tracklet.getExpPosErrorX(z_hodo) + xfudge;
            double err_y = factor*tracklet.getExpPosErrorY(z_hodo);

            double x_min = x_mask_min[idx1][idx2] - err_x;
            double x_max = x_mask_max[idx1][idx2] + err_x;
            double y_min = y_mask_min[idx1][idx2] - err_y;
            double y_max = y_mask_max[idx1][idx2] + err_y;

#ifdef _DEBUG_ON
            LogInfo(*iter);
            hitAll[*iter].print();
            LogInfo(nHodoHits << "/" << stationIDs_mask[tracklet.stationID-1].size() << ":  " << z_hodo << "  " << x_hodo << " +/- " << err_x << "  " << y_hodo << " +/-" << err_y << " : " << x_min << "  " << x_max << "  " << y_min << "  " << y_max);
#endif
            if(x_hodo > x_min && x_hodo < x_max && y_hodo > y_min && y_hodo < y_max)
            {
                nHodoHits++;
                masked = true;

		if(TRACK_ELECTRONS && tracklet.stationID > 5) return true; //Once the first hodoscope hit is found (at z=1420cm), the combined tracklet passes for electron tracks

                break;
            }
        }

        if(!masked) return false;
    }
  }

#ifdef _DEBUG_ON
    LogInfo(tracklet.stationID << "  " << nHodoHits << "  " << stationIDs_mask[tracklet.stationID-1].size());
#endif
    return true;
}

bool KalmanFastTracking::muonID_search(Tracklet& tracklet)
{
    //Set the cut value on multiple scattering
    //multiple scattering: sigma = 0.0136*sqrt(L/L0)*(1. + 0.038*ln(L/L0))/P, L = 1m, L0 = 1.76cm
    double cut = 0.03;
    if(tracklet.stationID == nStations)
    {
        double cut_the = MUID_THE_P0*tracklet.invP;
        double cut_emp = MUID_EMP_P0 + MUID_EMP_P1/tracklet.invP + MUID_EMP_P2/tracklet.invP/tracklet.invP;
        cut = MUID_REJECTION*(cut_the > cut_emp ? cut_the : cut_emp);
    }

    double slope[2] = {tracklet.tx, tracklet.ty};
    double pos_absorb[2] = {tracklet.getExpPositionX(MUID_Z_REF), tracklet.getExpPositionY(MUID_Z_REF)};
    PropSegment* segs[2] = {&(tracklet.seg_x), &(tracklet.seg_y)};
    for(int i = 0; i < 2; ++i)
    {
        //this shorting circuting can only be done to X-Z, Y-Z needs more complicated thing
        //if(i == 0 && segs[i]->getNHits() > 2 && segs[i]->isValid() > 0 && fabs(slope[i] - segs[i]->a) < cut) continue;

        segs[i]->init();
        for(int j = 0; j < 4; ++j)
        {
            int index = detectorIDs_muid[i][j] - nChamberPlanes - 1;
            double pos_ref = j < 2 ? pos_absorb[i] : segs[i]->getPosRef(pos_absorb[i] + slope[i]*(z_ref_muid[i][j] - MUID_Z_REF));
            double pos_exp = slope[i]*(z_mask[index] - z_ref_muid[i][j]) + pos_ref;

            if(!p_geomSvc->isInPlane(detectorIDs_muid[i][j], tracklet.getExpPositionX(z_mask[index]), tracklet.getExpPositionY(z_mask[index]))) continue;

            double win_tight = cut*(z_mask[index] - z_ref_muid[i][j]);
            win_tight = win_tight > 2.54 ? win_tight : 2.54;
            double win_loose = win_tight*2;
            double dist_min = 1E6;
            for(std::list<int>::iterator iter = hitIDs_muid[i][j].begin(); iter != hitIDs_muid[i][j].end(); ++iter)
            {
                double pos = hitAll[*iter].pos;
                double dist = pos - pos_exp;
                if(dist < -win_loose) continue;
                if(dist > win_loose) break;

                double dist_l = fabs(pos - hitAll[*iter].driftDistance - pos_exp);
                double dist_r = fabs(pos + hitAll[*iter].driftDistance - pos_exp);
                dist = dist_l < dist_r ? dist_l : dist_r;
                if(dist < dist_min)
                {
                    dist_min = dist;
                    if(dist < win_tight)
                    {
                        segs[i]->hits[j].hit = hitAll[*iter];
                        segs[i]->hits[j].sign = fabs(pos - hitAll[*iter].driftDistance - pos_exp) < fabs(pos + hitAll[*iter].driftDistance - pos_exp) ? -1 : 1;
                    }
                }
            }
        }
        segs[i]->fit();

        //this shorting circuting can only be done to X-Z, Y-Z needs more complicated thing
        //if(i == 0 && !(segs[i]->isValid() > 0 && fabs(slope[i] - segs[i]->a) < cut)) return false;
    }

    muonID_hodoAid(tracklet);
    if(segs[0]->getNHits() + segs[1]->getNHits() >= MUID_MINHITS)
    {
        return true;
    }
    else if(segs[1]->getNHits() == 1 || segs[1]->getNPlanes() == 1)
    {
        return segs[1]->nHodoHits >= 2;
    }
    return false;
}

bool KalmanFastTracking::muonID_comp(Tracklet& tracklet)
{
    //Set the cut value on multiple scattering
    //multiple scattering: sigma = 0.0136*sqrt(L/L0)*(1. + 0.038*ln(L/L0))/P, L = 1m, L0 = 1.76cm
    double cut = 0.03;
    if(tracklet.stationID == nStations)
    {
        double cut_the = MUID_THE_P0*tracklet.invP;
        double cut_emp = MUID_EMP_P0 + MUID_EMP_P1/tracklet.invP + MUID_EMP_P2/tracklet.invP/tracklet.invP;
        cut = MUID_REJECTION*(cut_the > cut_emp ? cut_the : cut_emp);
    }
#ifdef _DEBUG_ON
    LogInfo("Muon ID cut is: " << cut << " rad.");
#endif

    double slope[2] = {tracklet.tx, tracklet.ty};
    PropSegment* segs[2] = {&(tracklet.seg_x), &(tracklet.seg_y)};

    for(int i = 0; i < 2; ++i)
    {
#ifdef _DEBUG_ON
        if(i == 0) LogInfo("Working in X-Z:");
        if(i == 1) LogInfo("Working in Y-Z:");
#endif

        double pos_ref = i == 0 ? tracklet.getExpPositionX(MUID_Z_REF) : tracklet.getExpPositionY(MUID_Z_REF);
        if(segs[i]->getNHits() > 2 && segs[i]->isValid() > 0 && fabs(slope[i] - segs[i]->a) < cut && fabs(segs[i]->getExpPosition(MUID_Z_REF) - pos_ref) < MUID_R_CUT)
        {
#ifdef _DEBUG_ON
            LogInfo("Muon ID are already avaiable!");
#endif
            continue;
        }

        for(std::list<PropSegment>::iterator iter = propSegs[i].begin(); iter != propSegs[i].end(); ++iter)
        {
#ifdef _DEBUG_ON
            LogInfo("Testing this prop segment, with ref pos = " << pos_ref << ", slope_ref = " << slope[i]);
            iter->print();
#endif
            if(fabs(iter->a - slope[i]) < cut && fabs(iter->getExpPosition(MUID_Z_REF) - pos_ref) < MUID_R_CUT)
            {
                *(segs[i]) = *iter;
#ifdef _DEBUG_ON
                LogInfo("Accepted!");
#endif
                break;
            }
        }

        if(segs[i]->isValid() == 0) return false;
    }

    if(segs[0]->getNHits() + segs[1]->getNHits() < MUID_MINHITS) return false;
    return true;
}

bool KalmanFastTracking::muonID_hodoAid(Tracklet& tracklet)
{
    double win = 0.03;
    double factor = 5.;
    if(tracklet.stationID == nStations)
    {
        double win_the = MUID_THE_P0*tracklet.invP;
        double win_emp = MUID_EMP_P0 + MUID_EMP_P1/tracklet.invP + MUID_EMP_P2/tracklet.invP/tracklet.invP;
        win = MUID_REJECTION*(win_the > win_emp ? win_the : win_emp);
        factor = 3.;
    }

    PropSegment* segs[2] = {&(tracklet.seg_x), &(tracklet.seg_y)};
    for(int i = 0; i < 2; ++i)
    {
        segs[i]->nHodoHits = 0;
        for(std::list<int>::iterator iter = hitIDs_muidHodoAid[i].begin(); iter != hitIDs_muidHodoAid[i].end(); ++iter)
        {
            int detectorID = hitAll[*iter].detectorID;
            int elementID = hitAll[*iter].elementID;

            int idx1 = detectorID - nChamberPlanes - 1;
            int idx2 = elementID - 1;

            double z_hodo = z_mask[idx1];
            double x_hodo = tracklet.getExpPositionX(z_hodo);
            double y_hodo = tracklet.getExpPositionY(z_hodo);
            double err_x = factor*tracklet.getExpPosErrorX(z_hodo) + win*(z_hodo - MUID_Z_REF);
            double err_y = factor*tracklet.getExpPosErrorY(z_hodo) + win*(z_hodo - MUID_Z_REF);

            err_x = err_x/(x_mask_max[idx1][idx2] - x_mask_min[idx1][idx2]) > 0.25 ? 0.25*err_x/(x_mask_max[idx1][idx2] - x_mask_min[idx1][idx2]) : err_x;
            err_y = err_y/(y_mask_max[idx1][idx2] - y_mask_min[idx1][idx2]) > 0.25 ? 0.25*err_y/(y_mask_max[idx1][idx2] - y_mask_min[idx1][idx2]) : err_y;

            double x_min = x_mask_min[idx1][idx2] - err_x;
            double x_max = x_mask_max[idx1][idx2] + err_x;
            double y_min = y_mask_min[idx1][idx2] - err_y;
            double y_max = y_mask_max[idx1][idx2] + err_y;

            if(x_hodo > x_min && x_hodo < x_max && y_hodo > y_min && y_hodo < y_max)
            {
                segs[i]->hodoHits[segs[i]->nHodoHits++] = hitAll[*iter];
                if(segs[i]->nHodoHits > 4) break;
            }
        }
    }

    return true;
}

void KalmanFastTracking::buildPropSegments()
{
#ifdef _DEBUG_ON
    LogInfo("Building prop. tube segments");
#endif

    for(int i = 0; i < 2; ++i)
    {
        propSegs[i].clear();

        //note for prop tubes superID index starts from 4
        std::list<SRawEvent::hit_pair> pairs_forward  = rawEvent->getPartialHitPairsInSuperDetector(superIDs[i+5][0]);
        std::list<SRawEvent::hit_pair> pairs_backward = rawEvent->getPartialHitPairsInSuperDetector(superIDs[i+5][1]);

#ifdef _DEBUG_ON
        std::cout << "superID: " << superIDs[i+5][0] << ", " << superIDs[i+5][1] << std::endl;
        for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_forward.begin(); iter != pairs_forward.end(); ++iter)
        	LogInfo("Forward: " << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << "  " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
        for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_backward.begin(); iter != pairs_backward.end(); ++iter)
        	LogInfo("Backward: " << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << "  " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
#endif

        for(std::list<SRawEvent::hit_pair>::iterator fiter = pairs_forward.begin(); fiter != pairs_forward.end(); ++fiter)
        {
#ifdef _DEBUG_ON
            LogInfo("Trying forward pair " << fiter->first << "  " << fiter->second);
#endif
            for(std::list<SRawEvent::hit_pair>::iterator biter = pairs_backward.begin(); biter != pairs_backward.end(); ++biter)
            {
#ifdef _DEBUG_ON
                LogInfo("Trying backward pair " << biter->first << "  " << biter->second);
#endif

                PropSegment seg;

                //Note that the backward plane comes as the first in pair
                if(fiter->first >= 0) seg.hits[1] = SignedHit(hitAll[fiter->first], 0);
                if(fiter->second >= 0) seg.hits[0] = SignedHit(hitAll[fiter->second], 0);
                if(biter->first >= 0) seg.hits[3] = SignedHit(hitAll[biter->first], 0);
                if(biter->second >= 0) seg.hits[2] = SignedHit(hitAll[biter->second], 0);

#ifdef _DEBUG_ON
                seg.print();
#endif
                seg.fit();
#ifdef _DEBUG_ON
                seg.print();
#endif

                if(seg.isValid() > 0)
                {
                    propSegs[i].push_back(seg);
                }
#ifdef _DEBUG_ON
                else
                {
                    LogInfo("Rejected!");
                }
#endif
            }
        }
    }
}


int KalmanFastTracking::fitTracklet(Tracklet& tracklet)
{
    tracklet_curr = tracklet;

    //idx = 0, using simplex; idx = 1 using migrad
    int idx = 1;
#ifdef _ENABLE_MULTI_MINI
    if(tracklet.stationID < nStations-1) idx = 0;
#endif

    minimizer[idx]->SetLimitedVariable(0, "tx", tracklet.tx, 0.001, -TX_MAX, TX_MAX);
    minimizer[idx]->SetLimitedVariable(1, "ty", tracklet.ty, 0.001, -TY_MAX, TY_MAX);
    minimizer[idx]->SetLimitedVariable(2, "x0", tracklet.x0, 0.1, -X0_MAX, X0_MAX);
    minimizer[idx]->SetLimitedVariable(3, "y0", tracklet.y0, 0.1, -Y0_MAX, Y0_MAX);
    if(KMAG_ON)
    {
        minimizer[idx]->SetLimitedVariable(4, "invP", tracklet.invP, 0.001*tracklet.invP, INVP_MIN, INVP_MAX);
    }
    minimizer[idx]->Minimize();

    tracklet.tx = minimizer[idx]->X()[0];
    tracklet.ty = minimizer[idx]->X()[1];
    tracklet.x0 = minimizer[idx]->X()[2];
    tracklet.y0 = minimizer[idx]->X()[3];

    tracklet.err_tx = minimizer[idx]->Errors()[0];
    tracklet.err_ty = minimizer[idx]->Errors()[1];
    tracklet.err_x0 = minimizer[idx]->Errors()[2];
    tracklet.err_y0 = minimizer[idx]->Errors()[3];

    if(KMAG_ON && tracklet.stationID == nStations)
    {
        tracklet.invP = minimizer[idx]->X()[4];
        tracklet.err_invP = minimizer[idx]->Errors()[4];
    }

    tracklet.chisq = minimizer[idx]->MinValue();

    int status = minimizer[idx]->Status();
    return status;
}

int KalmanFastTracking::reduceTrackletList(std::list<Tracklet>& tracklets)
{
    std::list<Tracklet> targetList;

    tracklets.sort();
    while(!tracklets.empty())
    {
        targetList.push_back(tracklets.front());
        tracklets.pop_front();

#ifdef _DEBUG_ON_LEVEL_2
        LogInfo("Current best tracklet in reduce");
        targetList.back().print();
#endif

        for(std::list<Tracklet>::iterator iter = tracklets.begin(); iter != tracklets.end(); )
        {
            if(iter->similarity(targetList.back()))
            {
#ifdef _DEBUG_ON_LEVEL_2
                LogInfo("Removing this tracklet: ");
                iter->print();
#endif
                iter = tracklets.erase(iter);
                continue;
            }
            else
            {
                ++iter;
            }
        }
    }

    tracklets.assign(targetList.begin(), targetList.end());
    return 0;
}

void KalmanFastTracking::getExtrapoWindowsInSt1(Tracklet& tracklet, double* pos_exp, double* window, int st1ID)
{
    if(tracklet.stationID != nStations-1)
    {
        for(int i = 0; i < 3; i++)
        {
            pos_exp[i] = 9999.;
            window[i] = 0.;
        }
        return;
    }

    for(int i = 0; i < 3; i++)
    {
        int detectorID = (st1ID-1)*6 + 2*i + 2;
        int idx = p_geomSvc->getPlaneType(detectorID) - 1;

        double z_st1 = z_plane[detectorID];
        double x_st1 = tracklet.getExpPositionX(z_st1);
        double y_st1 = tracklet.getExpPositionY(z_st1);
        double err_x = tracklet.getExpPosErrorX(z_st1);
        double err_y = tracklet.getExpPosErrorY(z_st1);

        pos_exp[idx] = p_geomSvc->getUinStereoPlane(detectorID, x_st1, y_st1);
        window[idx]  = 5.*(fabs(costheta_plane[detectorID]*err_x) + fabs(sintheta_plane[detectorID]*err_y));
    }
}

void KalmanFastTracking::getSagittaWindowsInSt1(Tracklet& tracklet, double* pos_exp, double* window, int st1ID)
{
    if(tracklet.stationID != nStations-1)
    {
        for(int i = 0; i < 3; i++)
        {
            pos_exp[i] = 9999.;
            window[i] = 0.;
        }
        return;
    }

    double z_st3 = z_plane[tracklet.hits.back().hit.detectorID];
    double x_st3 = tracklet.getExpPositionX(z_st3);
    double y_st3 = tracklet.getExpPositionY(z_st3);

    //For U, X, and V planes
    for(int i = 0; i < 3; i++)
    {
        int detectorID = (st1ID-1)*6 + 2*i + 2;
        int idx = p_geomSvc->getPlaneType(detectorID) - 1;

        if(!(idx >= 0 && idx <3)) continue;

        double pos_st3 = p_geomSvc->getUinStereoPlane(s_detectorID[idx], x_st3, y_st3);

        double z_st1 = z_plane[detectorID];
        double z_st2 = z_plane[s_detectorID[idx]];
        double x_st2 = tracklet.getExpPositionX(z_st2);
        double y_st2 = tracklet.getExpPositionY(z_st2);
        double pos_st2 = p_geomSvc->getUinStereoPlane(s_detectorID[idx], x_st2, y_st2);

        double s2_target = pos_st2 - pos_st3*(z_st2 - Z_TARGET)/(z_st3 - Z_TARGET);
        double s2_dump   = pos_st2 - pos_st3*(z_st2 - Z_DUMP)/(z_st3 - Z_DUMP);

        double pos_exp_target = SAGITTA_TARGET_CENTER*s2_target + pos_st3*(z_st1 - Z_TARGET)/(z_st3 - Z_TARGET);
        double pos_exp_dump   = SAGITTA_DUMP_CENTER*s2_dump + pos_st3*(z_st1 - Z_DUMP)/(z_st3 - Z_DUMP);
        double win_target = fabs(s2_target*SAGITTA_TARGET_WIDTH);
        double win_dump   = fabs(s2_dump*SAGITTA_DUMP_WIDTH);

        double p_min = std::min(pos_exp_target - win_target, pos_exp_dump - win_dump);
        double p_max = std::max(pos_exp_target + win_target, pos_exp_dump + win_dump);

	std::cout<<"Test sagitta window thing.  The detectorID is "<<detectorID<<" the idx is "<<idx<<" the z_st1 is "<<z_st1<<std::endl; //WPM
	std::cout<<"More sagitta.  pos_exp is "<<0.5*(p_max + p_min)<<" and window = "<<0.5*(p_max - p_min)<<std::endl; //WPM
        pos_exp[idx] = 0.5*(p_max + p_min);
        window[idx]  = 0.5*(p_max - p_min);
    }
}

void KalmanFastTracking::printAtDetectorBack(int stationID, std::string outputFileName)
{
    TCanvas c1;

    std::vector<double> x, y, dx, dy;
    for(std::list<Tracklet>::iterator iter = trackletsInSt[stationID].begin(); iter != trackletsInSt[stationID].end(); ++iter)
    {
        double z = p_geomSvc->getPlanePosition(iter->stationID*6);
        x.push_back(iter->getExpPositionX(z));
        y.push_back(iter->getExpPositionY(z));
        dx.push_back(iter->getExpPosErrorX(z));
        dy.push_back(iter->getExpPosErrorY(z));
    }

    TGraphErrors gr(x.size(), &x[0], &y[0], &dx[0], &dy[0]);
    gr.SetMarkerStyle(8);

    //Add detector frames
    std::vector<double> x_f, y_f, dx_f, dy_f;
    x_f.push_back(p_geomSvc->getPlaneCenterX(stationID*6 + 6));
    y_f.push_back(p_geomSvc->getPlaneCenterY(stationID*6 + 6));
    dx_f.push_back(p_geomSvc->getPlaneScaleX(stationID*6 + 6)*0.5);
    dy_f.push_back(p_geomSvc->getPlaneScaleY(stationID*6 + 6)*0.5);

    if(stationID == 2)
    {
        x_f.push_back(p_geomSvc->getPlaneCenterX(stationID*6 + 12));
        y_f.push_back(p_geomSvc->getPlaneCenterY(stationID*6 + 12));
        dx_f.push_back(p_geomSvc->getPlaneScaleX(stationID*6 + 12)*0.5);
        dy_f.push_back(p_geomSvc->getPlaneScaleY(stationID*6 + 12)*0.5);
    }

    TGraphErrors gr_frame(x_f.size(), &x_f[0], &y_f[0], &dx_f[0], &dy_f[0]);
    gr_frame.SetLineColor(kRed);
    gr_frame.SetLineWidth(2);
    gr_frame.SetFillColor(15);

    c1.cd();
    gr_frame.Draw("A2[]");
    gr.Draw("Psame");

    c1.SaveAs(outputFileName.c_str());
}

SRecTrack KalmanFastTracking::processOneTracklet(Tracklet& tracklet)
{
    //tracklet.print();
    KalmanTrack kmtrk;
    kmtrk.setTracklet(tracklet);

    /*
    //Set the whole hit and node list
    for(std::list<SignedHit>::iterator iter = tracklet.hits.begin(); iter != tracklet.hits.end(); ++iter)
    {
        if(iter->hit.index < 0) continue;

        Node node_add(*iter);
        kmtrk.getNodeList().push_back(node_add);
        kmtrk.getHitIndexList().push_back(iter->sign*iter->hit.index);
    }

    //Set initial state
    TrkPar trkpar_curr;
    trkpar_curr._z = p_geomSvc->getPlanePosition(kmtrk.getNodeList().back().getHit().detectorID);
    //FIXME Debug Testing: sign reverse
    //TODO seems to be fixed
    trkpar_curr._state_kf[0][0] = tracklet.getCharge()*tracklet.invP/sqrt(1. + tracklet.tx*tracklet.tx + tracklet.ty*tracklet.ty);
    trkpar_curr._state_kf[1][0] = tracklet.tx;
    trkpar_curr._state_kf[2][0] = tracklet.ty;
    trkpar_curr._state_kf[3][0] = tracklet.getExpPositionX(trkpar_curr._z);
    trkpar_curr._state_kf[4][0] = tracklet.getExpPositionY(trkpar_curr._z);

    trkpar_curr._covar_kf.Zero();
    trkpar_curr._covar_kf[0][0] = 0.001;//1E6*tracklet.err_invP*tracklet.err_invP;
    trkpar_curr._covar_kf[1][1] = 0.01;//1E6*tracklet.err_tx*tracklet.err_tx;
    trkpar_curr._covar_kf[2][2] = 0.01;//1E6*tracklet.err_ty*tracklet.err_ty;
    trkpar_curr._covar_kf[3][3] = 100;//1E6*tracklet.getExpPosErrorX(trkpar_curr._z)*tracklet.getExpPosErrorX(trkpar_curr._z);
    trkpar_curr._covar_kf[4][4] = 100;//1E6*tracklet.getExpPosErrorY(trkpar_curr._z)*tracklet.getExpPosErrorY(trkpar_curr._z);

    kmtrk.setCurrTrkpar(trkpar_curr);
    kmtrk.getNodeList().back().getPredicted() = trkpar_curr;
    */

    /*
    //Fit the track first with possibily a few nodes unresolved
    if(!fitTrack(kmtrk))
    {
#ifdef _DEBUG_ON
        LogInfo("!fitTrack(kmtrk) - try flip charge");
#endif
        trkpar_curr._state_kf[0][0] *= -1.;
        kmtrk.setCurrTrkpar(trkpar_curr);
        kmtrk.getNodeList().back().getPredicted() = trkpar_curr;
        if(!fitTrack(kmtrk))
        {
#ifdef _DEBUG_ON
            LogInfo("!fitTrack(kmtrk) - failed flip charge also");
#endif
            SRecTrack strack = tracklet.getSRecTrack();
            strack.setKalmanStatus(-1);
            return strack;
        }
    }

    if(!kmtrk.isValid()) 
    {
#ifdef _DEBUG_ON
        LogInfo("!kmtrk.isValid() Chi2 = " << kmtrk.getChisq() << " - try flip charge");
#endif
        trkpar_curr._state_kf[0][0] *= -1.;
        kmtrk.setCurrTrkpar(trkpar_curr);
        kmtrk.getNodeList().back().getPredicted() = trkpar_curr;
        if(!fitTrack(kmtrk))
        {
#ifdef _DEBUG_ON
            LogInfo("!fitTrack(kmtrk) - failed flip charge also");
#endif
            SRecTrack strack = tracklet.getSRecTrack();
            strack.setKalmanStatus(-1);
            return strack;
        }

#ifdef _DEBUG_ON
        LogInfo("Chi2 after flip charge: " << kmtrk.getChisq());
        if(kmtrk.isValid()) 
        {
            LogInfo("flip charge worked!");
        }
#endif
    }

#ifdef _DEBUG_ON
    LogInfo("kmtrk.print()");
    kmtrk.print();
    LogInfo("kmtrk.printNodes()");
    kmtrk.printNodes();
#endif
    */

    //Resolve left-right based on the current solution, re-fit if anything changed
    //resolveLeftRight(kmtrk);
    if(fitTrack(kmtrk) && kmtrk.isValid())
    {
        SRecTrack strack = kmtrk.getSRecTrack();

        //Set trigger road ID
        TriggerRoad road(tracklet);
        strack.setTriggerRoad(road.getRoadID());

        //Set prop tube slopes
        strack.setNHitsInPT(tracklet.seg_x.getNHits(), tracklet.seg_y.getNHits());
        strack.setPTSlope(tracklet.seg_x.a, tracklet.seg_y.a);

        strack.setKalmanStatus(1);

        return strack;
    }
    else
    {
#ifdef _DEBUG_ON
    	LogInfo("!kmtrk.isValid()");
#endif
        SRecTrack strack = tracklet.getSRecTrack();
        strack.setKalmanStatus(-1);

        return strack;
    }
}

bool KalmanFastTracking::fitTrack(KalmanTrack& kmtrk)
{
    if(kmtrk.getNodeList().empty()) return false;

    if(kmfitter->processOneTrack(kmtrk) == 0)
    {
        return false;
    }
    kmfitter->updateTrack(kmtrk);

    return true;
}

void KalmanFastTracking::resolveLeftRight(KalmanTrack& kmtrk)
{
    bool isUpdated = false;

    std::list<int>::iterator hitID = kmtrk.getHitIndexList().begin();
    for(std::list<Node>::iterator node = kmtrk.getNodeList().begin(); node != kmtrk.getNodeList().end(); )
    {
        if(*hitID == 0)
        {
            double x_hit = node->getSmoothed().get_x();
            double y_hit = node->getSmoothed().get_y();
            double pos_hit = p_geomSvc->getUinStereoPlane(node->getHit().detectorID, x_hit, y_hit);

            int sign = 0;
            if(pos_hit > node->getHit().pos)
            {
                sign = 1;
            }
            else
            {
                sign = -1;
            }

            //update the node list
            TMatrixD m(1, 1), dm(1, 1);
            m[0][0] = node->getHit().pos + sign*node->getHit().driftDistance;
            dm[0][0] = p_geomSvc->getPlaneResolution(node->getHit().detectorID)*p_geomSvc->getPlaneResolution(node->getHit().detectorID);
            node->setMeasurement(m, dm);
            *hitID = sign*node->getHit().index;

            isUpdated = true;
        }

        ++node;
        ++hitID;
    }

    if(isUpdated) fitTrack(kmtrk);
}

void KalmanFastTracking::printTimers() {
	std::cout <<"KalmanFastTracking::printTimers: " << std::endl;
	std::cout <<"================================================================" << std::endl;
	std::cout << "Tracklet St2                "<<_timers["st2"]->get_accumulated_time()/1000. << " sec" <<std::endl;
	std::cout << "Tracklet St3                "<<_timers["st3"]->get_accumulated_time()/1000. << " sec" <<std::endl;
	std::cout << "Tracklet St23               "<<_timers["st23"]->get_accumulated_time()/1000. << " sec" <<std::endl;
	std::cout << "Tracklet Global             "<<_timers["global"]->get_accumulated_time()/1000. << " sec" <<std::endl;
	std::cout << "  Global St1                "<<_timers["global_st1"]->get_accumulated_time()/1000. << " sec" <<std::endl;
	std::cout << "  Global Link               "<<_timers["global_link"]->get_accumulated_time()/1000. << " sec" <<std::endl;
	std::cout << "  Global Kalman             "<<_timers["global_kalman"]->get_accumulated_time()/1000. << " sec" <<std::endl;
	std::cout << "Tracklet Kalman             "<<_timers["kalman"]->get_accumulated_time()/1000. << " sec" <<std::endl;
	std::cout <<"================================================================" << std::endl;
}

void KalmanFastTracking::chi2fit(int n, double x[], double y[], double& a, double& b)
{
    double sum = 0.;
    double sx = 0.;
    double sy = 0.;
    double sxx = 0.;
    double syy = 0.;
    double sxy = 0.;

    for(int i = 0; i < n; ++i)
    {
        ++sum;
        sx += x[i];
        sy += y[i];
        sxx += (x[i]*x[i]);
        syy += (y[i]*y[i]);
        sxy += (x[i]*y[i]);
    }

    double det = sum*sxx - sx*sx;
    if(fabs(det) < 1E-20)
    {
        a = 0.;
        b = 0.;

        return;
    }

    a = (sum*sxy - sx*sy)/det;
    b = (sy*sxx - sxy*sx)/det;
}

bool KalmanFastTracking::compareTracklets(Tracklet& tracklet2, Tracklet& tracklet3)
{
  //Here we will compare the possible X-Z slopes within the station 2 and station 3 tracklets
  Tracklet::linedef line2X;
  Tracklet::linedef line3X;
  Tracklet::linedef line2X_v2;
  Tracklet::linedef line3X_v2;
  
  //It is rare, but sometimes, you will have slopes that match coincidentally.  Therefore, I keep track of best two combinations.  This seems to be sufficienct
  double slopeComp = 1.0;
  double secondSlope = 1.1;
  for(unsigned int t2 = 0; t2 < tracklet2.possibleXLines.size(); t2++){
    for(unsigned int t3 = 0; t3 < tracklet3.possibleXLines.size(); t3++){
      if(std::abs(tracklet3.possibleXLines.at(t3).slopeX - tracklet2.possibleXLines.at(t2).slopeX) < slopeComp){
	
	//if the new combination is the closest so far, then the previous closest becomes the second closest...
	secondSlope = slopeComp;
	line2X_v2 = line2X;
	line3X_v2 = line3X;

	slopeComp = std::abs(tracklet3.possibleXLines.at(t3).slopeX - tracklet2.possibleXLines.at(t2).slopeX);
	line2X = tracklet2.possibleXLines.at(t2);
	line3X = tracklet3.possibleXLines.at(t3);
      }
      else if(std::abs(tracklet3.possibleXLines.at(t3).slopeX - tracklet2.possibleXLines.at(t2).slopeX) < secondSlope){
	//not as close as the closest combination, but closer than the previously existing second combination
	secondSlope = std::abs(tracklet3.possibleXLines.at(t3).slopeX - tracklet2.possibleXLines.at(t2).slopeX);
        line2X_v2 = tracklet2.possibleXLines.at(t2);
        line3X_v2 = tracklet3.possibleXLines.at(t3);
      }
    }
  }

  if(slopeComp > 0.04) return false; //This has not been optimized at all.  I just chose a random value (previous slope comparison allowed for a difference of 0.1)
  double extrapolation = line2X.slopeX*(line3X.initialZ - line2X.initialZ) + line2X.initialX;

  if(std::abs(extrapolation - line3X.initialX) > 5. ){ //allow for a 5 cm difference of the tracklet in station 3 from the station 2 extrapolation.  This also should be optimized
    if(secondSlope > 0.04) return false; //If both the closest and second closest slopes don't match, then this is not a good combination
    double extrapolation_v2 = line2X_v2.slopeX*(line3X_v2.initialZ - line2X_v2.initialZ) + line2X_v2.initialX; //Perform the extrapolation in the rare case that the closest slope combination did not yield a valid extrapolation.  Rare, but necessary
    if(std::abs(extrapolation_v2 - line3X_v2.initialX) > 5. ){ //Same window size!  Could be optimized
      return false;
    } else{
      line2X = line2X_v2; //These are the possible X-Z lines that we actually want, if the closest combination wasn't valid based on the extrapolation
      line3X = line3X_v2;
    }
  }

  Tracklet::linedef line2U;
  Tracklet::linedef line3U;
  double slopeCompU = 1.0;
  for(unsigned int t2 = 0; t2 < tracklet2.possibleULines.size(); t2++){
    for(unsigned int t3 = 0; t3 < tracklet3.possibleULines.size(); t3++){
      if(std::abs(tracklet3.possibleULines.at(t3).slopeU - tracklet2.possibleULines.at(t2).slopeU) < slopeCompU){
	slopeCompU = std::abs(tracklet3.possibleULines.at(t3).slopeU - tracklet2.possibleULines.at(t2).slopeU);
	line2U = tracklet2.possibleULines.at(t2);
	line3U = tracklet3.possibleULines.at(t3);
      }
    }
  } //As of now, I don't keep track of the second-closest combination for the U and V layers

  if(slopeCompU > 0.07) return false; //Larger window here.  From what I can tell, the resolution is worse in this plane, or maybe my slope calculations are somewhat incorrect

  Tracklet::linedef line2V;
  Tracklet::linedef line3V;  
  double slopeCompV = 1.0;
  for(unsigned int t2 = 0; t2 < tracklet2.possibleVLines.size(); t2++){
    for(unsigned int t3 = 0; t3 < tracklet3.possibleVLines.size(); t3++){
      if(std::abs(tracklet3.possibleVLines.at(t3).slopeV - tracklet2.possibleVLines.at(t2).slopeV) < slopeCompV){
	slopeCompV = std::abs(tracklet3.possibleVLines.at(t3).slopeV - tracklet2.possibleVLines.at(t2).slopeV);
	line2V = tracklet2.possibleVLines.at(t2);
        line3V = tracklet3.possibleVLines.at(t3);
      }
    }
  }

  if(slopeCompV > 0.07) return false; //same comment about resolution as for the U layer

  //Now we find the Y-values of the hits in the U and V planes
  double tracklet2Ys_D = 0.; //This is the sum of the Y-values of the hits in the U and V planes of the two tracklets.  The sum is taken to be used in an average.  The _D here stands for driftDistance.  I originally did this part of the code without taking the drift distance in the slanted layers into account
  double tracklet3Ys_D = 0.;
  tracklet2Ys_D += line2V.wire1Slope * (line2X.slopeX*(line2V.wireHit1PosZ - line2X.initialZ) + line2X.initialX) + line2V.wireIntercept1;
  tracklet2Ys_D += line2V.wire2Slope * (line2X.slopeX*(line2V.wireHit2PosZ - line2X.initialZ) + line2X.initialX) + line2V.wireIntercept2;
  tracklet2Ys_D += line2U.wire1Slope * (line2X.slopeX*(line2U.wireHit1PosZ - line2X.initialZ) + line2X.initialX) + line2U.wireIntercept1;
  tracklet2Ys_D += line2U.wire2Slope * (line2X.slopeX*(line2U.wireHit2PosZ - line2X.initialZ) + line2X.initialX) + line2U.wireIntercept2;

  tracklet3Ys_D += line3V.wire1Slope * (line3X.slopeX*(line3V.wireHit1PosZ - line3X.initialZ) + line3X.initialX) + line3V.wireIntercept1;
  tracklet3Ys_D += line3V.wire2Slope * (line3X.slopeX*(line3V.wireHit2PosZ - line3X.initialZ) + line3X.initialX) + line3V.wireIntercept2;
  tracklet3Ys_D += line3U.wire1Slope * (line3X.slopeX*(line3U.wireHit1PosZ - line3X.initialZ) + line3X.initialX) + line3U.wireIntercept1;
  tracklet3Ys_D += line3U.wire2Slope * (line3X.slopeX*(line3U.wireHit2PosZ - line3X.initialZ) + line3X.initialX) + line3U.wireIntercept2;

  //Give the station 2 and station 3 tracklets the same tx and ty value.  I could get an X0 and Y0 extrapolation, but that doesn't seem to be strictly necessary.  The X0 and Y0 values are found in the fittracklet function for the combined station 2 + station 3 tracklet
  tracklet2.tx = (line2X.slopeX + line3X.slopeX)/2;
  tracklet2.ty = (tracklet3Ys_D/4. - tracklet2Ys_D/4.)/(line3U.wireHit1PosZ - line2V.wireHit1PosZ); //The y slope is found by taking the average Y position in the station3 and subtracting the average Y position in station2.  This is then divided by the z difference, of course  
  tracklet3.tx = (line2X.slopeX + line3X.slopeX)/2;
  tracklet3.ty = (tracklet3Ys_D/4. - tracklet2Ys_D/4.)/(line3U.wireHit1PosZ - line2V.wireHit1PosZ);

  tracklet2.st2Z = line2X.initialZ;
  tracklet2.st2X = line2X.initialX;
  tracklet2.st2Y = tracklet2Ys_D/4.;
  tracklet2.st3Y = tracklet3Ys_D/4.;
  tracklet2.st2U = line2U.initialU;
  tracklet2.st2V = line2V.initialV;
  tracklet2.st2Usl = line2U.slopeU;
  tracklet2.st2Vsl = line2V.slopeV;
  tracklet2.st2UZ = line2U.wireHit1PosZ;
  tracklet2.st2VZ = line2V.wireHit1PosZ;
  
  return true;
  
}


bool KalmanFastTracking::compareTrackletsSlim(Tracklet& tracklet2, Tracklet& tracklet3)
{
  //Here we will compare the possible X-Z slopes within the station 2 and station 3 tracklets
  Tracklet::linedef line2X;
  Tracklet::linedef line3X;
  Tracklet::linedef line2X_v2;
  Tracklet::linedef line3X_v2;
  
  //It is rare, but sometimes, you will have slopes that match coincidentally.  Therefore, I keep track of best two combinations.  This seems to be sufficienct
  double slopeComp = 1.0;
  double secondSlope = 1.1;
  for(unsigned int t2 = 0; t2 < tracklet2.possibleXLines.size(); t2++){
    for(unsigned int t3 = 0; t3 < tracklet3.possibleXLines.size(); t3++){
      if(std::abs(tracklet3.possibleXLines.at(t3).slopeX - tracklet2.possibleXLines.at(t2).slopeX) < slopeComp){
	
	//if the new combination is the closest so far, then the previous closest becomes the second closest...
	secondSlope = slopeComp;
	line2X_v2 = line2X;
	line3X_v2 = line3X;

	slopeComp = std::abs(tracklet3.possibleXLines.at(t3).slopeX - tracklet2.possibleXLines.at(t2).slopeX);
	line2X = tracklet2.possibleXLines.at(t2);
	line3X = tracklet3.possibleXLines.at(t3);
      }
      else if(std::abs(tracklet3.possibleXLines.at(t3).slopeX - tracklet2.possibleXLines.at(t2).slopeX) < secondSlope){
	//not as close as the closest combination, but closer than the previously existing second combination
	secondSlope = std::abs(tracklet3.possibleXLines.at(t3).slopeX - tracklet2.possibleXLines.at(t2).slopeX);
        line2X_v2 = tracklet2.possibleXLines.at(t2);
        line3X_v2 = tracklet3.possibleXLines.at(t3);
      }
    }
  }

  if(slopeComp > 0.005) return false; //This has not been optimized at all.  I just chose a random value (previous slope comparison allowed for a difference of 0.1)
  double extrapolation = line2X.slopeX*(line3X.initialZ - line2X.initialZ) + line2X.initialX;
  double extrapolation_toSt1 = line2X.slopeX*(z_plane[3] - line2X.initialZ) + line2X.initialX;
  std::cout<<"extrapolation is "<<extrapolation<<"; based on slope "<<line2X.slopeX<<", diff "<<line3X.initialZ - line2X.initialZ<<", and initialX "<<line2X.initialX<<std::endl;
  std::cout<<"to be compared with st3 position "<<line3X.initialX<<".  diff is "<<std::abs(extrapolation - line3X.initialX)<<std::endl;
  std::cout<<"extrapolation to station 1 is "<<extrapolation_toSt1<<std::endl;
  
  if(std::abs(extrapolation - line3X.initialX) > 5. || std::abs(extrapolation_toSt1) > 100. ){ //allow for a 5 cm difference of the tracklet in station 3 from the station 2 extrapolation.  This also should be optimized
    if(secondSlope > 0.005) return false; //If both the closest and second closest slopes don't match, then this is not a good combination
    double extrapolation_v2 = line2X_v2.slopeX*(line3X_v2.initialZ - line2X_v2.initialZ) + line2X_v2.initialX; //Perform the extrapolation in the rare case that the closest slope combination did not yield a valid extrapolation.  Rare, but necessary
    double extrapolation_toSt1_v2 = line2X_v2.slopeX*(z_plane[3] - line2X_v2.initialZ) + line2X_v2.initialX;
      std::cout<<"V2 extrapolation is "<<extrapolation_v2<<"; based on slope "<<line2X_v2.slopeX<<", diff "<<line3X_v2.initialZ - line2X_v2.initialZ<<", and initialX "<<line2X_v2.initialX<<std::endl;
      std::cout<<"V2 to be compared with st3 position "<<line3X_v2.initialX<<".  diff is "<<std::abs(extrapolation_v2 - line3X_v2.initialX)<<std::endl;
    if(std::abs(extrapolation_v2 - line3X_v2.initialX) > 5. || std::abs(extrapolation_toSt1_v2) > 100. ){ //Same window size!  Could be optimized
      return false;
    } else{
      line2X = line2X_v2; //These are the possible X-Z lines that we actually want, if the closest combination wasn't valid based on the extrapolation
      line3X = line3X_v2;
    }
  }

  tracklet2.acceptedXLine2 = line2X;
  tracklet3.acceptedXLine3 = line3X;
  //Give the station 2 and station 3 tracklets the same tx and ty value.  I could get an X0 and Y0 extrapolation, but that doesn't seem to be strictly necessary.  The X0 and Y0 values are found in the fittracklet function for the combined station 2 + station 3 tracklet
  tracklet2.tx = (line2X.slopeX + line3X.slopeX)/2;
  //tracklet2.ty = (tracklet3Ys_D/4. - tracklet2Ys_D/4.)/(line3U.wireHit1PosZ - line2V.wireHit1PosZ); //The y slope is found by taking the average Y position in the station3 and subtracting the average Y position in station2.  This is then divided by the z difference, of course  
  tracklet3.tx = (line2X.slopeX + line3X.slopeX)/2;
  //tracklet3.ty = (tracklet3Ys_D/4. - tracklet2Ys_D/4.)/(line3U.wireHit1PosZ - line2V.wireHit1PosZ);

  tracklet2.st2Z = line2X.initialZ;
  tracklet2.st2X = line2X.initialX;
  tracklet3.st2Z = line3X.initialZ;
  tracklet3.st2X = line3X.initialX;

  tracklet2.st2Xsl = line2X.slopeX;
  tracklet3.st2Xsl = line3X.slopeX;

  std::cout<<"THERE WAS A TRACKLET MATCH!"<<std::endl;
  std::cout<<"st2 slope = "<<line2X.slopeX<<"; st2 X = "<<line2X.initialX<<std::endl;
  std::cout<<"st3 slope = "<<line3X.slopeX<<"; st3 X = "<<line3X.initialX<<std::endl;
  
  return true;
  
}


bool KalmanFastTracking::compareTrackletsSlimU(Tracklet& tracklet2, Tracklet& tracklet3)
{
  //Here we will compare the possible X-Z slopes within the station 2 and station 3 tracklets
  Tracklet::linedef line2U;
  Tracklet::linedef line3U;
  Tracklet::linedef line2U_v2;
  Tracklet::linedef line3U_v2;
  
  //It is rare, but sometimes, you will have slopes that match coincidentally.  Therefore, I keep track of best two combinations.  This seems to be sufficienct
  double slopeComp = 1.0;
  double secondSlope = 1.1;
  for(unsigned int t2 = 0; t2 < tracklet2.possibleULines.size(); t2++){
    for(unsigned int t3 = 0; t3 < tracklet3.possibleULines.size(); t3++){
      if(std::abs(tracklet3.possibleULines.at(t3).slopeU - tracklet2.possibleULines.at(t2).slopeU) < slopeComp){
	
	//if the new combination is the closest so far, then the previous closest becomes the second closest...
	secondSlope = slopeComp;
	line2U_v2 = line2U;
	line3U_v2 = line3U;

	slopeComp = std::abs(tracklet3.possibleULines.at(t3).slopeU - tracklet2.possibleULines.at(t2).slopeU);
	line2U = tracklet2.possibleULines.at(t2);
	line3U = tracklet3.possibleULines.at(t3);
      }
      else if(std::abs(tracklet3.possibleULines.at(t3).slopeU - tracklet2.possibleULines.at(t2).slopeU) < secondSlope){
	//not as close as the closest combination, but closer than the previously existing second combination
	secondSlope = std::abs(tracklet3.possibleULines.at(t3).slopeU - tracklet2.possibleULines.at(t2).slopeU);
        line2U_v2 = tracklet2.possibleULines.at(t2);
        line3U_v2 = tracklet3.possibleULines.at(t3);
      }
    }
  }

  if(slopeComp > 0.015) return false; //This has not been optimized at all.  I just chose a random value (previous slope comparison allowed for a difference of 0.1)
  std::cout<<"not a match yet"<<std::endl;
  std::cout<<"U st2 slope = "<<line2U.slopeU<<"; st2 U = "<<line2U.initialU<<std::endl;
  std::cout<<"U st3 slope = "<<line3U.slopeU<<"; st3 U = "<<line3U.initialU<<std::endl;
  

  double extrapolation = line2U.slopeU*(line3U.initialZ - line2U.initialZ) + line2U.initialU;
  double extrapolation_toSt1 = line2U.slopeU*(z_plane[3] - line2U.initialZ) + line2U.initialU;

  std::cout<<"U extrapolation is "<<extrapolation<<"; based on slope "<<line2U.slopeU<<", diff "<<line3U.initialZ - line2U.initialZ<<", and initialU "<<line2U.initialU<<std::endl;
  std::cout<<"U to be compared with st3 position "<<line3U.initialU<<".  diff is "<<(extrapolation - line3U.initialU)<<std::endl;
  std::cout<<"U extrapolation to station 1 is "<<extrapolation_toSt1<<std::endl;
  std::cout<<"U slope diff are "<<slopeComp<<" and "<<secondSlope<<std::endl;
  
  if(std::abs(extrapolation - line3U.initialU) > 16. || std::abs(extrapolation_toSt1) > 100. || std::abs(line3U.slopeU) > 0.15 || std::abs(line2U.slopeU) > 0.15 ){ //allow for a 5 cm difference of the tracklet in station 3 from the station 2 extrapolation.  This also should be optimized
    if(secondSlope > 0.015) return false; //If both the closest and second closest slopes don't match, then this is not a good combination
    double extrapolation_v2 = line2U_v2.slopeU*(line3U_v2.initialZ - line2U_v2.initialZ) + line2U_v2.initialU; //Perform the extrapolation in the rare case that the closest slope combination did not yield a valid extrapolation.  Rare, but necessary
    double extrapolation_toSt1_v2 = line2U_v2.slopeU*(z_plane[3] - line2U_v2.initialZ) + line2U_v2.initialU;
      std::cout<<"U V2 extrapolation is "<<extrapolation_v2<<"; based on slope "<<line2U_v2.slopeU<<", diff "<<line3U_v2.initialZ - line2U_v2.initialZ<<", and initialU "<<line2U_v2.initialU<<std::endl;
      std::cout<<"U V2 to be compared with st3 position "<<line3U_v2.initialU<<".  diff is "<<(extrapolation_v2 - line3U_v2.initialU)<<std::endl;
      if(std::abs(extrapolation_v2 - line3U_v2.initialU) > 13. || std::abs(extrapolation_toSt1_v2) > 100. || std::abs(line3U_v2.slopeU) > 0.15 || std::abs(line2U_v2.slopeU) > 0.15 ){ //Same window size!  Could be optimized
      return false;
    } else{
      line2U = line2U_v2; //These are the possible U-Z lines that we actually want, if the closest combination wasn't valid based on the extrapolation
      line3U = line3U_v2;
    }
  }

  tracklet2.acceptedULine2 = line2U;
  tracklet3.acceptedULine3 = line3U;
  
  //Give the station 2 and station 3 tracklets the same tx and ty value.  I could get an X0 and Y0 extrapolation, but that doesn't seem to be strictly necessary.  The X0 and Y0 values are found in the fittracklet function for the combined station 2 + station 3 tracklet
  //tracklet2.tx = (line2X.slopeX + line3X.slopeX)/2;
  //tracklet2.ty = (tracklet3Ys_D/4. - tracklet2Ys_D/4.)/(line3U.wireHit1PosZ - line2V.wireHit1PosZ); //The y slope is found by taking the average Y position in the station3 and subtracting the average Y position in station2.  This is then divided by the z difference, of course  
  //tracklet3.tx = (line2X.slopeX + line3X.slopeX)/2;
  //tracklet3.ty = (tracklet3Ys_D/4. - tracklet2Ys_D/4.)/(line3U.wireHit1PosZ - line2V.wireHit1PosZ);

  tracklet2.st2Z = line2U.initialZ;
  tracklet2.st2U = line2U.initialU;
  tracklet3.st2Z = line3U.initialZ;
  tracklet3.st2U = line3U.initialU;

  tracklet2.st2Usl = line2U.slopeU;
  tracklet3.st2Usl = line3U.slopeU;

  std::cout<<"THERE WAS A U TRACKLET MATCH!"<<std::endl;
  std::cout<<"U st2 slope = "<<line2U.slopeX<<"; st2 U = "<<line2U.initialU<<std::endl;
  std::cout<<"U st3 slope = "<<line3U.slopeX<<"; st3 U = "<<line3U.initialU<<std::endl;
  
  return true;
  
}


bool KalmanFastTracking::compareTrackletsSlimV(Tracklet& tracklet2, Tracklet& tracklet3)
{
  //Here we will compare the possible X-Z slopes within the station 2 and station 3 tracklets
  Tracklet::linedef line2V;
  Tracklet::linedef line3V;
  Tracklet::linedef line2V_v2;
  Tracklet::linedef line3V_v2;
  
  //It is rare, but sometimes, you will have slopes that match coincidentally.  Therefore, I keep track of best two combinations.  This seems to be sufficienct
  double slopeComp = 1.0;
  double secondSlope = 1.1;
  for(unsigned int t2 = 0; t2 < tracklet2.possibleVLines.size(); t2++){
    for(unsigned int t3 = 0; t3 < tracklet3.possibleVLines.size(); t3++){
      if(std::abs(tracklet3.possibleVLines.at(t3).slopeV - tracklet2.possibleVLines.at(t2).slopeV) < slopeComp){
	
	//if the new combination is the closest so far, then the previous closest becomes the second closest...
	secondSlope = slopeComp;
	line2V_v2 = line2V;
	line3V_v2 = line3V;

	slopeComp = std::abs(tracklet3.possibleVLines.at(t3).slopeV - tracklet2.possibleVLines.at(t2).slopeV);
	line2V = tracklet2.possibleVLines.at(t2);
	line3V = tracklet3.possibleVLines.at(t3);
      }
      else if(std::abs(tracklet3.possibleVLines.at(t3).slopeV - tracklet2.possibleVLines.at(t2).slopeV) < secondSlope){
	//not as close as the closest combination, but closer than the previously existing second combination
	secondSlope = std::abs(tracklet3.possibleVLines.at(t3).slopeV - tracklet2.possibleVLines.at(t2).slopeV);
        line2V_v2 = tracklet2.possibleVLines.at(t2);
        line3V_v2 = tracklet3.possibleVLines.at(t3);
      }
    }
  }

  if(slopeComp > 0.015) return false; //This has not been optimized at all.  I just chose a random value (previous slope comparison allowed for a difference of 0.1)
  std::cout<<"not a match yet"<<std::endl;
  std::cout<<"V st2 slope = "<<line2V.slopeV<<"; st2 V = "<<line2V.initialV<<std::endl;
  std::cout<<"V st3 slope = "<<line3V.slopeV<<"; st3 V = "<<line3V.initialV<<std::endl;
  

  double extrapolation = line2V.slopeV*(line3V.initialZ - line2V.initialZ) + line2V.initialV;
  double extrapolation_toSt1 = line2V.slopeV*(z_plane[3] - line2V.initialZ) + line2V.initialV;

  std::cout<<"V extrapolation is "<<extrapolation<<"; based on slope "<<line2V.slopeV<<", diff "<<line3V.initialZ - line2V.initialZ<<", and initialV "<<line2V.initialV<<std::endl;
  std::cout<<"V to be compared with st3 position "<<line3V.initialV<<".  diff is "<<(extrapolation - line3V.initialV)<<std::endl;
  std::cout<<"V extrapolation to station 1 is "<<extrapolation_toSt1<<std::endl;
  std::cout<<"V slope diff are "<<slopeComp<<" and "<<secondSlope<<std::endl;
  
  if(std::abs(extrapolation - line3V.initialV) > 16. || std::abs(extrapolation_toSt1) > 100. || std::abs(line3V.slopeV) > 0.15 || std::abs(line2V.slopeV) > 0.15 ){ //allow for a 5 cm difference of the tracklet in station 3 from the station 2 extrapolation.  This also should be optimized
    if(secondSlope > 0.015) return false; //If both the closest and second closest slopes don't match, then this is not a good combination
    double extrapolation_v2 = line2V_v2.slopeV*(line3V_v2.initialZ - line2V_v2.initialZ) + line2V_v2.initialV; //Perform the extrapolation in the rare case that the closest slope combination did not yield a valid extrapolation.  Rare, but necessary
    double extrapolation_toSt1_v2 = line2V_v2.slopeV*(z_plane[3] - line2V_v2.initialZ) + line2V_v2.initialV;
      std::cout<<"V V2 extrapolation is "<<extrapolation_v2<<"; based on slope "<<line2V_v2.slopeV<<", diff "<<line3V_v2.initialZ - line2V_v2.initialZ<<", and initialV "<<line2V_v2.initialV<<std::endl;
      std::cout<<"V V2 to be compared with st3 position "<<line3V_v2.initialV<<".  diff is "<<(extrapolation_v2 - line3V_v2.initialV)<<std::endl;
      if(std::abs(extrapolation_v2 - line3V_v2.initialV) > 13. || std::abs(extrapolation_toSt1_v2) > 100. || std::abs(line3V_v2.slopeV) > 0.15 || std::abs(line2V_v2.slopeV) > 0.15 ){ //Same window size!  Could be optimized
      return false;
    } else{
      line2V = line2V_v2; //These are the possible V-Z lines that we actually want, if the closest combination wasn't valid based on the extrapolation
      line3V = line3V_v2;
    }
  }

  tracklet2.acceptedVLine2 = line2V;
  std::cout<<"TEST OF LINE2 "<<tracklet2.acceptedVLine2.wireHit1PosZ<<std::endl;
  tracklet3.acceptedVLine3 = line3V;
  
  //Give the station 2 and station 3 tracklets the same tx and ty value.  I could get an X0 and Y0 extrapolation, but that doesn't seem to be strictly necessary.  The X0 and Y0 values are found in the fittracklet function for the combined station 2 + station 3 tracklet
  //tracklet2.tx = (line2X.slopeX + line3X.slopeX)/2;
  //tracklet2.ty = (tracklet3Ys_D/4. - tracklet2Ys_D/4.)/(line3V.wireHit1PosZ - line2V.wireHit1PosZ); //The y slope is found by taking the average Y position in the station3 and subtracting the average Y position in station2.  This is then divided by the z difference, of course  
  //tracklet3.tx = (line2X.slopeX + line3X.slopeX)/2;
  //tracklet3.ty = (tracklet3Ys_D/4. - tracklet2Ys_D/4.)/(line3V.wireHit1PosZ - line2V.wireHit1PosZ);

  tracklet2.st2Z = line2V.initialZ;
  tracklet2.st2V = line2V.initialV;
  tracklet3.st2Z = line3V.initialZ;
  tracklet3.st2V = line3V.initialV;

  tracklet2.st2Vsl = line2V.slopeV;
  tracklet3.st2Vsl = line3V.slopeV;

  std::cout<<"THERE WAS A V TRACKLET MATCH!"<<std::endl;
  std::cout<<"V st2 slope = "<<line2V.slopeX<<"; st2 V = "<<line2V.initialV<<std::endl;
  std::cout<<"V st3 slope = "<<line3V.slopeX<<"; st3 V = "<<line3V.initialV<<std::endl;
  
  return true;
  
}
