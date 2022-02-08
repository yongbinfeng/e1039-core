/*
KalmanFastTracking_NEW_2.cxx

Implementation of class Tracklet, KalmanFastTracking_NEW_2

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

#include "KalmanFastTracking_NEW_2.h"
#include "TriggerRoad.h"

//#define _DEBUG_ON
//#define _DEBUG_PATRICK
//#define _DEBUG_PATRICK_EXTRA

//#define _DEBUG_RES
//#define _DEBUG_RES_EMBED

//#define _DEBUG_FAST

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

  static double KMAGSTR;
  static double PT_KICK_KMAG;
  
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

	    KMAGSTR = rc->get_DoubleFlag("KMAGSTR");
	    PT_KICK_KMAG = rc->get_DoubleFlag("PT_KICK_KMAG")*KMAGSTR;
        }
    }
}

KalmanFastTracking_NEW_2::KalmanFastTracking_NEW_2(const PHField* field, const TGeoManager* geom, bool flag): verbosity(0), enable_KF(flag), outputListIdx(4)
{
    using namespace std;
    initGlobalVariables();

#ifdef _DEBUG_ON
    cout << "Initialization of KalmanFastTracking_NEW_2 ..." << endl;
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

KalmanFastTracking_NEW_2::~KalmanFastTracking_NEW_2()
{
    if(enable_KF) delete kmfitter;
    delete minimizer[0];
    delete minimizer[1];
}

void KalmanFastTracking_NEW_2::setRawEventDebug(SRawEvent* event_input)
{
    rawEvent = event_input;
    hitAll = event_input->getAllHits();
}

int KalmanFastTracking_NEW_2::setRawEvent(SRawEvent* event_input)
{
	//reset timer
    for(auto iter=_timers.begin(); iter != _timers.end(); ++iter) 
    {
        iter->second->reset();
    }

    //Initialize tracklet lists
    for(int i = 0; i < 5; i++){
      trackletsInSt[i].clear();
      for(int b = 0; b < 200; b++){
	trackletsInStSlimX[i][b].clear();
	trackletsInStSlimU[i][b].clear();
	trackletsInStSlimV[i][b].clear();
      }
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

    //Get hit combinations in station 2
    _timers["st2"]->restart();
    buildTrackletsInStationSlim(3, 1);   //3 for station-2, 1 for list position 1
    if(verbosity >= 2) std::cout<<"Station 2 x combos = "<<trackletsInStSlimX[1][0].size()<<std::endl;
    buildTrackletsInStationSlimU(3, 1);   //3 for station-2, 1 for list position 1
    if(verbosity >= 2) std::cout<<"Station 2 u combos = "<<trackletsInStSlimU[1][0].size()<<std::endl;
    buildTrackletsInStationSlimV(3, 1);   //3 for station-2, 1 for list position 1
    if(verbosity >= 2) std::cout<<"Station 2 v combos = "<<trackletsInStSlimV[1][0].size()<<std::endl;
    _timers["st2"]->stop();


    //Get hit combinations in station 3+ and 3-
    _timers["st3"]->restart();
    buildTrackletsInStationSlim(4, 2);   //4 for station-3+
    buildTrackletsInStationSlim(5, 2);   //5 for station-3-
    if(verbosity >= 2) std::cout<<"Station 3 x combos = "<<trackletsInStSlimX[2][0].size()<<std::endl;
    buildTrackletsInStationSlimU(4, 2);   //4 for station-3+
    buildTrackletsInStationSlimU(5, 2);   //5 for station-3-
    if(verbosity >= 2) std::cout<<"Station 3 u combos = "<<trackletsInStSlimU[2][0].size()<<std::endl;
    buildTrackletsInStationSlimV(4, 2);   //4 for station-3+
    buildTrackletsInStationSlimV(5, 2);   //5 for station-3-
    if(verbosity >= 2) std::cout<<"Station 3 v combos = "<<trackletsInStSlimV[2][0].size()<<std::endl;
    _timers["st3"]->stop();

    
    //Matching of station 2 and station 3 hits separately for the three wire tilts
    _timers["st23"]->restart();
    buildBackPartialTracksSlimX(1, m_slopeComparison, m_windowSize);
    buildBackPartialTracksSlimU(1, m_slopeComparison, m_windowSize);
    buildBackPartialTracksSlimV(1, m_slopeComparison, m_windowSize);

    getNum23Combos();
    if(verbosity >= 2){
      std::cout<<"Station 2+3 x combos = "<<num23XCombos<<std::endl;
      std::cout<<"Station 2+3 u combos = "<<num23UCombos<<std::endl;
      std::cout<<"Station 2+3 v combos = "<<num23VCombos<<std::endl;
    }
    
    bool checkedX = false;
    bool checkedU = false;
    bool checkedV = false;

    //Use tighter hit requirements if the number of combinations in a wire-tilt set is large.  The below buildBackPartialTracksSlim calls require that all 4 wires have a hit for the particular wire tilt.  NOTE: alternatively, we could try changing the extrapolation/slope-matching requirements rather than the hit requirements!  I.e. we could use m_slopeComparisonMedium/m_slopeComparisonTight and m_windowSizeMedium/m_windowSizeTight
    if(num23XCombos > 900){
      for(int b = 0; b < 200; b++){
	trackletsInStSlimX[3][b].clear();
      }
      buildBackPartialTracksSlimX(2, m_slopeComparison, m_windowSize);
      checkedX = true;
      getNum23Combos();
      if(verbosity >= 2) std::cout<<"Station 2+3 x combos tight = "<<num23XCombos<<std::endl;
    }
    if(num23UCombos > 900){
      for(int b = 0; b < 200; b++){
	trackletsInStSlimU[3][b].clear();
      }
      buildBackPartialTracksSlimU(2, m_slopeComparison, m_windowSize);
      checkedU = true;
      getNum23Combos();
      if(verbosity >= 2) std::cout<<"Station 2+3 u combos tight = "<<num23UCombos<<std::endl;
    }
    if(num23VCombos > 900){
      for(int b = 0; b < 200; b++){
	trackletsInStSlimV[3][b].clear();
      }
      buildBackPartialTracksSlimV(2, m_slopeComparison, m_windowSize);
      checkedV = true;
      getNum23Combos();
      if(verbosity >= 2) std::cout<<"Station 2+3 v combos tight = "<<num23VCombos<<std::endl;
    }
    
    /* //This commented out section is what I used to use for PU mitigation before transitioning to the binned method
    if(trackletsInStSlimX[3].size() > 300 && !checkedX){
      trackletsInStSlimX[3].clear();
      buildBackPartialTracksSlimX(2, m_slopeComparisonMedium, m_windowSizeMedium);
      std::cout<<"Station 2+3 x combos medium = "<<trackletsInStSlimX[3].size()<<std::endl;
      if(trackletsInStSlimX[3].size() > 150){
	trackletsInStSlimX[3].clear();
	buildBackPartialTracksSlimX(2, m_slopeComparisonTight, m_windowSizeTight);
	std::cout<<"Station 2+3 x combos tight = "<<trackletsInStSlimX[3].size()<<std::endl;
      }
    }
    if(trackletsInStSlimU[3].size() > 300 && !checkedU){
      trackletsInStSlimU[3].clear();
      buildBackPartialTracksSlimU(2, m_slopeComparisonMedium, m_windowSizeMedium);
      std::cout<<"Station 2+3 u combos medium = "<<trackletsInStSlimU[3].size()<<std::endl;
      if(trackletsInStSlimU[3].size() > 150){
	trackletsInStSlimU[3].clear();
	buildBackPartialTracksSlimU(2, m_slopeComparisonTight, m_windowSizeTight);
	std::cout<<"Station 2+3 u combos tight = "<<trackletsInStSlimU[3].size()<<std::endl;
      }
    }
    if(trackletsInStSlimV[3].size() > 300 && !checkedV){
      trackletsInStSlimV[3].clear();
      buildBackPartialTracksSlimV(2, m_slopeComparisonMedium, m_windowSizeMedium);
      std::cout<<"Station 2+3 v combos medium = "<<trackletsInStSlimV[3].size()<<std::endl;
      if(trackletsInStSlimV[3].size() > 150){
	trackletsInStSlimV[3].clear();
	buildBackPartialTracksSlimV(2, m_slopeComparisonTight, m_windowSizeTight);
	std::cout<<"Station 2+3 v combos tight = "<<trackletsInStSlimV[3].size()<<std::endl;
      }
    }
    if( (trackletsInStSlimX[3].size() * trackletsInStSlimU[3].size() * trackletsInStSlimV[3].size()) > 2000000 ){
      trackletsInStSlimX[3].clear();
      trackletsInStSlimU[3].clear();
      trackletsInStSlimV[3].clear();
      buildBackPartialTracksSlimX(2, m_slopeComparisonTight, m_windowSizeTight);
      buildBackPartialTracksSlimU(2, m_slopeComparisonTight, m_windowSizeTight);
      buildBackPartialTracksSlimV(2, m_slopeComparisonTight, m_windowSizeTight);
    }*/
    /*if(trackletsInStSlimV[3].size() > 500){
      trackletsInStSlimV[3].clear();
      buildBackPartialTracksSlimV(2);
      std::cout<<"Station 2+3 v combos = "<<trackletsInStSlimV[3].size()<<std::endl;
      if(trackletsInStSlimV[3].size() > 100){
	return TFEXIT_FAIL_BACKPARTIAL;
      }
    }*/

    
    std::cout<<"num23XCombos*num23UCombos*num23VCombos = "<<num23XCombos*num23UCombos*num23VCombos<<std::endl;
    if(num23XCombos*num23UCombos*num23VCombos < 1000000){
      buildBackPartialTracksSlim_v3(0); //This should be relatively fast given the binned combination method
    } else if( isSlimMiddle() ){ //There were over 1,000,000 possible combinations.  Are there relatively few combinations in the middle 60cm of the detector?
      buildBackPartialTracksSlim_v3(1); //This only looks for tracklets that are in the middle 60 cm of station 2.  This is to avoid PU hits in the edges of the detector, though you will miss signal particles outside of [-30cm, 30cm]
    } else{ //There were a lot of combinations even in the middle of the detector.  Let's tighten our hit requirements and extrapolation/slope-matching requirements

      if(num23XCombos > 300){
	for(int b = 0; b < 200; b++){
	  trackletsInStSlimX[3][b].clear();
	}
	buildBackPartialTracksSlimX(2, m_slopeComparisonMedium, m_windowSizeMedium);
	checkedX = true;
	getNum23Combos();
	if(verbosity >= 2) std::cout<<"Station 2+3 x combos slim tight = "<<num23XCombos<<std::endl;
      }
      if(num23UCombos > 300){
	for(int b = 0; b < 200; b++){
	  trackletsInStSlimU[3][b].clear();
	}
	buildBackPartialTracksSlimU(2, m_slopeComparisonMedium, m_windowSizeMedium);
	checkedU = true;
	getNum23Combos();
	if(verbosity >= 2) std::cout<<"Station 2+3 u combos tight = "<<num23UCombos<<std::endl;
      }
      if(num23VCombos > 300){
	for(int b = 0; b < 200; b++){
	  trackletsInStSlimV[3][b].clear();
	}
	buildBackPartialTracksSlimV(2, m_slopeComparisonMedium, m_windowSizeMedium);
	checkedV = true;
	getNum23Combos();
	if(verbosity >= 2) std::cout<<"Station 2+3 v combos tight = "<<num23VCombos<<std::endl;
      }
      
      if(isSlimMiddle()){ //Now that we've tighted the requirements, are there still too many potential combinations in the middle of the detector?
	buildBackPartialTracksSlim_v3(1); 
      } else{ //Just too many potential combinations at this point.  Giving up.  NOTE: there are probably ways to get around having to give up; just need to improve PU mitigation!
	return TFEXIT_FAIL_BACKPARTIAL;
      }      
      
    }
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

bool KalmanFastTracking_NEW_2::acceptEvent(SRawEvent* rawEvent)
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
    //if(rawEvent->getNHitsInD1() > MaxHitsDC1) return false; 
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




void KalmanFastTracking_NEW_2::buildBackPartialTracksSlim_v3(int cut)
{
#ifdef _DEBUG_PATRICK
  LogInfo("In buildBackPartialTracksSlim_v3");
#endif

  /**NOTE: THERE ARE SOME LOW-HANGING FRUITS HERE.  FOR EXAMPLE, THIS IS EFFECTIVELY A TRIPLE LOOP OVER PLANES.  WE COULD DO SOME STRAIGHFORWARD MATH TO FIGURE OUT THE LINE DEFINED BY THE INTERSECTIONS OF THE VARIOUS PLANES TO EXTRACT TRAJECTORY PARAMETERS AND CHECK COMPATIBILITY.**/
  
  double chiSqCut = 300; //This can be dynamically decreased as a means of PU mitigation
  //if( trackletsInStSlimX[3].size() > 150 || trackletsInStSlimU[3].size() > 150 || trackletsInStSlimV[3].size() > 150 ) chiSqCut = 100; //Example of cut decrease from old PU mitigation scheme
  
  int scanWidth = 142; //In the binned PU mitigation scheme, there are 200 bins that st2+st3 hit combinations fall into.  The bins have width of 2cm, and are based on the st2 position information.
  if(cut == 1) scanWidth = 30; //this will scan the innermost 60 cm of station 2 in x from -30 to +30
  
  for(int binx = 2; binx < scanWidth; binx++){ //Rather than scanning just from left to right in the detector, I scan outwards from the center of station 2 (there's more PU at the edges typically).  The outermost component of this triple loop is the vertical wire combinations
    int bx = 100-(binx % 2 == 0 ? -1*floor((binx-1)/2) : floor(binx/2)); //Again, I'm scanning from the center, so the bin steps are 100, 99, 101, 98, 102, 97, 103, etc.
    if(trackletsInStSlimX[3][bx].size()==0) continue; //Don't need to do anything if there aren't any hit combinations in this bin!
    Tracklet tracklet_best_bin; //keep track of the best st2+st3 full combination per bin.
    
    for(std::list<Tracklet>::iterator trackletX = trackletsInStSlimX[3][bx].begin(); trackletX != trackletsInStSlimX[3][bx].end(); ++trackletX){ //loop over the st2+st3 hit combinations for the vertical wires in this positional bin
      Tracklet tracklet_best;

#ifdef _DEBUG_FAST
      std::cout<<"m_v3 try tracklet X:"<<std::endl;
      trackletX->print();
#endif
      
      for(int bu = std::max(0,bx-15); bu < std::min(200,bx+16); bu++){ //Rather than looking at all possible U combinations for each X combination, let's look at +/- 15 U bins.  Relatively large set of bins because U is slanted relative to X, and we don't know the y position or y-slant of the particle!
	if(trackletsInStSlimU[3][bu].size()==0) continue;
	Tracklet tracklet_best_UX_bin;
	
	for(std::list<Tracklet>::iterator trackletU = trackletsInStSlimU[3][bu].begin(); trackletU != trackletsInStSlimU[3][bu].end(); ++trackletU){ //loop over combinations in the U bin
	  Tracklet tracklet_best_UX;

#ifdef _DEBUG_FAST
	  std::cout<<"m_v3 try tracklet U:"<<std::endl;
	  trackletU->print();
#endif
	  
	  for(int bv = std::max(0, bx-1*(bu-bx)-2); bv < std::min(200,bx-1*(bu-bx)+3); bv++){ //For the bv bins, we combine our X and U position information.  We scan a window of combinations on the opposite side of the st2 X position from the U combination in question.  E.g., if the x bin was 55, and the U bin was 50, we would expect the correct V combination to be in the 60 bin.  However, we still do a bit of a scan to allow for imprecision in those rough measurements
	    if(trackletsInStSlimV[3][bv].size()==0) continue;
	    Tracklet tracklet_best_UV_bin;
	    
	    for(std::list<Tracklet>::iterator trackletV = trackletsInStSlimV[3][bv].begin(); trackletV != trackletsInStSlimV[3][bv].end(); ++trackletV){

#ifdef _DEBUG_FAST
	      std::cout<<"m_v3 try tracklet V:"<<std::endl;
	      trackletV->print();
#endif
	      
#ifdef _DEBUG_FAST
	      //if(trackletsInStSlimX[3][bx].size() * trackletsInStSlimU[3][bu].size() * trackletsInStSlimV[3][bv].size() < 300){
		LogInfo("New combo");
		std::cout<<"trackletX->st2X = "<<trackletX->st2X<<"; trackletU->st2U = "<<trackletU->st2U<<"; trackletV->st2V = "<<trackletV->st2V<<";;; trackletX->st3X = "<<trackletX->st3X<<"; trackletU->st3U = "<<trackletU->st3U<<"; trackletV->st3V = "<<trackletV->st3V<<std::endl;
		std::cout<<"trackletX->st2Xsl = "<<trackletX->st2Xsl<<"; trackletU->st2Usl = "<<trackletU->st2Usl<<"; trackletV->st2Vsl = "<<trackletV->st2Vsl<<";;; trackletX->st3Xsl = "<<trackletX->st3Xsl<<"; trackletU->st3Usl = "<<trackletU->st3Usl<<"; trackletV->st3Vsl = "<<trackletV->st3Vsl<<std::endl;
		trackletX->print();
		trackletU->print();
		trackletV->print();
		//}
#endif	

	      //quick checks on hit combination compatibility in the three wire slants
	      if( std::abs( (trackletX->st2Xsl - trackletU->st2Usl) - -1.*(trackletX->st2Xsl - trackletV->st2Vsl) ) > 0.04 ) continue;
	      if( std::abs( (trackletX->st2X - trackletU->st2U) - -1.*(trackletX->st2X - trackletV->st2V) ) > 7. ) continue;
	      if( std::abs( (trackletX->st3X - trackletU->st3U) - -1.*(trackletX->st3X - trackletV->st3V) ) > 7. ) continue;
	      
#ifdef _DEBUG_FAST
	      LogInfo("OK combo");
#endif
	      
	      double st2UWireSin = TMath::Sin(TMath::ATan(1./trackletU->acceptedULine2.wire1Slope));
	      double st2VWireSin = TMath::Sin(TMath::ATan(1./trackletV->acceptedVLine2.wire1Slope));
	      double st3UWireSin = TMath::Sin(TMath::ATan(1./trackletU->acceptedULine3.wire1Slope));
	      double st3VWireSin = TMath::Sin(TMath::ATan(1./trackletV->acceptedVLine3.wire1Slope));
	      
	      double st2UWireCos = TMath::Cos(TMath::ATan(1./trackletU->acceptedULine2.wire1Slope));
	      double st2VWireCos = TMath::Cos(TMath::ATan(1./trackletV->acceptedVLine2.wire1Slope));
	      double st3UWireCos = TMath::Cos(TMath::ATan(1./trackletU->acceptedULine3.wire1Slope));
	      double st3VWireCos = TMath::Cos(TMath::ATan(1./trackletV->acceptedVLine3.wire1Slope));
	      
	      double testTY2 = ( trackletV->st2Vsl - (st2VWireCos/st2UWireCos)*trackletU->st2Usl )/( (st2VWireCos * st2UWireSin)/(st2UWireCos) - st2VWireSin );
	      double testTY3 = ( trackletV->st3Vsl - (st3VWireCos/st3UWireCos)*trackletU->st3Usl )/( (st3VWireCos * st3UWireSin)/(st3UWireCos) - st3VWireSin );
	      double testTX2 = ( trackletU->st2Usl - (st2UWireSin/st2VWireSin)*trackletV->st2Vsl )/( st2UWireCos - ( st2UWireSin * st2VWireCos )/st2VWireSin );
	      double testTX3 = ( trackletU->st3Usl - (st3UWireSin/st3VWireSin)*trackletV->st3Vsl )/( st3UWireCos - ( st3UWireSin * st3VWireCos )/st3VWireSin );
	      
#ifdef _DEBUG_FAST
	      std::cout<<"testTY2 = "<<testTY2<<", and testTY3 = "<<testTY3<<std::endl;
	      std::cout<<"tracklet2.tx = "<<trackletX->st2Xsl<<", and tracklet3.tx = "<<trackletX->st2Xsl<<", and testTX2 = "<<testTX2<<", and testTX3 = "<<testTX3<<std::endl;
#endif
	      if(std::abs(testTY2 - testTY3) > 0.007) continue; //The y-slope extracted from the U wires should match the y-slope extracted from the V wires

#ifdef _DEBUG_FAST
	      LogInfo("past cont 1");
#endif
	      if(std::abs(trackletX->st2Xsl - testTX2) > 0.007 || std::abs(trackletX->st2Xsl - testTX3) > 0.007) continue; //The x-slopes extracted from the U and V wires should match the slope found from the X wires

#ifdef _DEBUG_FAST
	      LogInfo("past cont 2");
#endif


	      //Let's build a tracklet and assign the x0, tx, and ty parameters (and also the dummy invP parameter).  What we don't know right now is y0.  HOWEVER, IF WE DID SOME MATH TO FIND THE LINE DEFINED BY THE INTERSECTION OF THE PLANES, WE WOULD KNOW THE Y0 VALUE!
	      Tracklet tracklet_TEST_23 = (*trackletX) + (*trackletU) + (*trackletV);
	      tracklet_TEST_23.tx = trackletX->st2Xsl;
	      tracklet_TEST_23.ty = (testTY2 + testTY3)/2.; //take an average of the ty value extracted from the U and V values to hedge your bets
	      tracklet_TEST_23.invP = 1./50.;
	      tracklet_TEST_23.x0 = trackletX->st2Xsl * (0. - trackletX->st2Z) + trackletX->st2X; //simple exrapolation back to z = 0
	      
	      tracklet_TEST_23.sortHits();
	      
	      double testY0 = (-1.*p_geomSvc->getDCA((*trackletU).getHit(0).hit.detectorID, (*trackletU).getHit(0).hit.elementID, tracklet_TEST_23.tx, tracklet_TEST_23.ty, tracklet_TEST_23.x0, 0)/std::abs(st2UWireSin) + p_geomSvc->getDCA((*trackletV).getHit(0).hit.detectorID, (*trackletV).getHit(0).hit.elementID, tracklet_TEST_23.tx, tracklet_TEST_23.ty, tracklet_TEST_23.x0, 0)/std::abs(st2UWireSin))/2.; //This is a method to extract y0 based on what the distance of closest approach would be to the st2 U and V hits if the y0 was 0.  A little hacky.  Again, extracting this information from the line defined by the plane intersections would be much smarter and better
	      tracklet_TEST_23.y0 = testY0;
	      
	      
	      if( tracklet_TEST_23.calcChisq_noDrift() > chiSqCut || isnan(tracklet_TEST_23.calcChisq_noDrift()) ) continue; //check the chisq when drift distances are not accounted for.  (When using drift distance, the expected precision changes, driving up the chisq given that we only have a "rough" extraction of the trajectory parameters at this point)
#ifdef _DEBUG_FAST
	      LogInfo("past cont 4");
#endif	
	      //Assign some parameters that get used later
	      tracklet_TEST_23.st2Z = trackletX->st2Z;
	      tracklet_TEST_23.st2X = trackletX->st2X;
	      tracklet_TEST_23.st2Y = tracklet_TEST_23.y0 + tracklet_TEST_23.ty * tracklet_TEST_23.st2Z;
	      tracklet_TEST_23.st3Y = tracklet_TEST_23.y0 + tracklet_TEST_23.ty * trackletX->st3Z;
	      tracklet_TEST_23.st2U = trackletU->st2U;
	      tracklet_TEST_23.st2V = trackletV->st2V;
	      tracklet_TEST_23.st2Usl = trackletU->st2Usl;
	      tracklet_TEST_23.st2Vsl = trackletV->st2Vsl;

	      
	      fitTracklet(tracklet_TEST_23);
	      
	      //Assign hit sign for unknown hits
	      for(std::list<SignedHit>::iterator hit1 = tracklet_TEST_23.hits.begin(); hit1 != tracklet_TEST_23.hits.end(); ++hit1){
		if(hit1->sign == 0){
		  hit1->sign = 1;
		  fitTracklet(tracklet_TEST_23);
		  double dcaPlus = tracklet_TEST_23.chisq;
		  hit1->sign = -1;
		  fitTracklet(tracklet_TEST_23);
		  double dcaMinus = tracklet_TEST_23.chisq;
		  if(std::abs(dcaPlus) < std::abs(dcaMinus)){
		    hit1->sign = 1;
		  } else{
		    hit1->sign = -1;
		  }
		}
	      }
	      
	      ///Remove bad hits if needed;  Right now this doesn't play nicely with this algorithm
	      //removeBadHits(tracklet_TEST_23);
	      
	      fitTracklet(tracklet_TEST_23); //A final fit now that all hits have been assigned signs
	      
	      
	      if(tracklet_TEST_23.chisq > 9000.)
		{
#ifdef _DEBUG_ON
		  tracklet_TEST_23.print();
		  LogInfo("Impossible combination!");
#endif
		  continue;
		}
	      
#ifdef _DEBUG_PATRICK
	      LogInfo("New tracklet: ");
	      tracklet_TEST_23.print();
	      
	      LogInfo("Current best:");
	      tracklet_best_UV_bin.print();
	      
	      LogInfo("Comparison: " << (tracklet_TEST_23 < tracklet_best_UV_bin));
	      LogInfo("Candidate chisq: " << tracklet_TEST_23.chisq);
	      LogInfo("Quality: " << acceptTracklet(tracklet_TEST_23));
#endif
	      
	      //If current tracklet is better than the best tracklet up-to-now
	      if(acceptTracklet(tracklet_TEST_23) && tracklet_TEST_23 < tracklet_best_UV_bin)
		{
		  tracklet_best_UV_bin = tracklet_TEST_23;
		}
#ifdef _DEBUG_ON
	      else
		{
		  LogInfo("Rejected!!");
		}
#endif
	      
	    }
	    if(tracklet_best_UV_bin < tracklet_best_UX){
	      tracklet_best_UX = tracklet_best_UV_bin;
	    }
	  }
	  if(tracklet_best_UX < tracklet_best_UX_bin){
	    tracklet_best_UX_bin = tracklet_best_UX;
	  }
	}
	if(tracklet_best_UX_bin < tracklet_best){
	  tracklet_best = tracklet_best_UX_bin;
	}
      }
      if(tracklet_best < tracklet_best_bin){
	tracklet_best_bin = tracklet_best;
      }
    }
    if(acceptTracklet(tracklet_best_bin)){
      if(tracklet_best_bin.isValid() > 0){
	if(tracklet_best_bin.chisq < 15.){
	  trackletsInSt[3].push_back(tracklet_best_bin);
#ifdef _DEBUG_ON
	  tracklet_best_bin.print();
#endif
	  
	}
      }
    }
  }
}



void KalmanFastTracking_NEW_2::buildBackPartialTracksSlimX(int pass, double slopeComparison, double windowSize)
{
  
  for(std::list<Tracklet>::iterator tracklet3 = trackletsInStSlimX[2][0].begin(); tracklet3 != trackletsInStSlimX[2][0].end(); ++tracklet3)
    {
      
      Tracklet tracklet_best;
      for(std::list<Tracklet>::iterator tracklet2 = trackletsInStSlimX[1][0].begin(); tracklet2 != trackletsInStSlimX[1][0].end(); ++tracklet2)
        {
#ifdef _DEBUG_RES
	  if(trackletsInStSlimX[1][0].size() < 10 && trackletsInStSlimX[2][0].size() < 10){
	    tracklet2->print();
	    tracklet3->print();
	  }
#endif
	  
	  Tracklet tracklet_23;
	  if(OLD_TRACKING){
	    tracklet_23 = (*tracklet2) + (*tracklet3);
	  }
	  else{
	    //check if the combination is valid.  compareTrackletsSlim checks 2+2 hit combinations, and compareTrackletsSlim_3hits checks 2+1 combinations
	    if(compareTrackletsSlim(*tracklet2, *tracklet3, pass, slopeComparison, windowSize) || (pass == 1 && compareTrackletsSlim_3hits(*tracklet2, *tracklet3, pass, slopeComparison, windowSize)) ){
	      tracklet_23 = (*tracklet2) + (*tracklet3);
	      tracklet_23.tx = tracklet2->tx;
	      tracklet_23.st2X = tracklet2->st2X;
	      tracklet_23.st2Z = tracklet2->st2Z;
	      tracklet_23.st3X = tracklet3->st3X;
	      tracklet_23.st3Z = tracklet3->st2Z;
	      tracklet_23.st2Xsl = tracklet2->st2Xsl;
	      tracklet_23.st3Xsl = tracklet3->st2Xsl;
	      tracklet_23.acceptedXLine2 = tracklet2->acceptedXLine2;
	      tracklet_23.acceptedXLine3 = tracklet3->acceptedXLine3;
	      
#ifdef _DEBUG_ON
	      LogInfo("We had a match using following two tracklets:");
	      tracklet2->print();
	      tracklet3->print();
	      LogInfo("Yield this combination:");
	      tracklet_23.print();
#endif
	      
	      if(tracklet_23.st2X > -200. && tracklet_23.st2X < 200){
		int bin = floor( (tracklet_23.st2X + 200) / 2 ); //Put tracklet into a bin based on st2 position
		trackletsInStSlimX[3][bin].push_back(tracklet_23);
	      } else{
		continue;
	      }
	    }
	    else{
	      continue;
	    }
	  }
	}
    }
}

void KalmanFastTracking_NEW_2::buildBackPartialTracksSlimU(int pass, double slopeComparison, double windowSize)
{

  for(std::list<Tracklet>::iterator tracklet3 = trackletsInStSlimU[2][0].begin(); tracklet3 != trackletsInStSlimU[2][0].end(); ++tracklet3)
    {
      Tracklet tracklet_best;
      for(std::list<Tracklet>::iterator tracklet2 = trackletsInStSlimU[1][0].begin(); tracklet2 != trackletsInStSlimU[1][0].end(); ++tracklet2)
        {
#ifdef _DEBUG_RES
	  if(trackletsInStSlimU[1][0].size() < 10 && trackletsInStSlimU[2][0].size() < 10){
	    tracklet2->print();
	    tracklet3->print();
	  }
#endif
	  
	  Tracklet tracklet_23;
	  if(OLD_TRACKING){
	    tracklet_23 = (*tracklet2) + (*tracklet3);
	  }
	  else{
	    if(compareTrackletsSlimU(*tracklet2, *tracklet3, pass, slopeComparison, windowSize) || (pass == 1 && compareTrackletsSlimU_3hits(*tracklet2, *tracklet3, pass, slopeComparison, windowSize)) ){
	      tracklet_23 = (*tracklet2) + (*tracklet3);
	      tracklet_23.tx = tracklet2->tx;
	      tracklet_23.st2U = tracklet2->st2U;
	      tracklet_23.st2Z = tracklet2->st2Z;
	      tracklet_23.st3U = tracklet3->st3U;
	      tracklet_23.st3Z = tracklet3->st2Z;
	      tracklet_23.st2Usl = tracklet2->st2Usl;
	      tracklet_23.st3Usl = tracklet3->st2Usl;
	      tracklet_23.acceptedULine2 = tracklet2->acceptedULine2;
	      tracklet_23.acceptedULine3 = tracklet3->acceptedULine3;
	      
#ifdef _DEBUG_ON
	      LogInfo("We had a match using following two tracklets:");
	      tracklet2->print();
	      tracklet3->print();
	      LogInfo("Yield this combination:");
	      tracklet_23.print();
#endif
	      
	      if(tracklet_23.st2U > -200. && tracklet_23.st2U < 200){
		int bin = floor( (tracklet_23.st2U + 200) / 2 );
		trackletsInStSlimU[3][bin].push_back(tracklet_23);
	      } else{
		continue;
	      }
	    }
	    else{
	      continue;
	    }
	  }
	}
    }
}



void KalmanFastTracking_NEW_2::buildBackPartialTracksSlimV(int pass, double slopeComparison, double windowSize)
{
  
  for(std::list<Tracklet>::iterator tracklet3 = trackletsInStSlimV[2][0].begin(); tracklet3 != trackletsInStSlimV[2][0].end(); ++tracklet3)
    {
      Tracklet tracklet_best;
      for(std::list<Tracklet>::iterator tracklet2 = trackletsInStSlimV[1][0].begin(); tracklet2 != trackletsInStSlimV[1][0].end(); ++tracklet2)
        {
#ifdef _DEBUG_RES
	  if(trackletsInStSlimV[1][0].size() < 10 && trackletsInStSlimV[2][0].size() < 10){
	    tracklet2->print();
	    tracklet3->print();
	  }
#endif
	  
	  Tracklet tracklet_23;
	  if(OLD_TRACKING){
	    tracklet_23 = (*tracklet2) + (*tracklet3);
	  }
	  else{
	    if(compareTrackletsSlimV(*tracklet2, *tracklet3, pass, slopeComparison, windowSize) || (pass == 1 && compareTrackletsSlimV_3hits(*tracklet2, *tracklet3, pass, slopeComparison, windowSize)) ){
	      tracklet_23 = (*tracklet2) + (*tracklet3);
	      tracklet_23.tx = tracklet2->tx;
	      tracklet_23.st2V = tracklet2->st2V;
	      tracklet_23.st2Z = tracklet2->st2Z;
	      tracklet_23.st3V = tracklet3->st3V;
	      tracklet_23.st3Z = tracklet3->st2Z;
	      tracklet_23.st2Vsl = tracklet2->st2Vsl;
	      tracklet_23.st3Vsl = tracklet3->st2Vsl;
	      tracklet_23.acceptedVLine2 = tracklet2->acceptedVLine2;
	      tracklet_23.acceptedVLine3 = tracklet3->acceptedVLine3;
	      
#ifdef _DEBUG_ON
	      LogInfo("We had a match using following two tracklets:");
	      tracklet2->print();
	      tracklet3->print();
	      LogInfo("Yield this combination:");
	      tracklet_23.print();
#endif
	      
	      if(tracklet_23.st2V > -200. && tracklet_23.st2V < 200){
		int bin = floor( (tracklet_23.st2V + 200) / 2 );
		trackletsInStSlimV[3][bin].push_back(tracklet_23);
	      } else{
		continue;
	      }	    
	    }
	    else{
	      continue;
	    }
	  }
	}
    }
}

void KalmanFastTracking_NEW_2::buildGlobalTracks()
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


void KalmanFastTracking_NEW_2::buildGlobalTracksDisplaced()
{
    double pos_exp[3], window[3];
    for(std::list<Tracklet>::iterator tracklet23 = trackletsInSt[3].begin(); tracklet23 != trackletsInSt[3].end(); ++tracklet23)
    {

      Tracklet tracklet_best_prob, tracklet_best_vtx;

      double posx = tracklet23->tx * ( z_plane[3] - tracklet23->st2Z ) + tracklet23->st2X;
      double posu = tracklet23->st2Usl * ( z_plane[1] - tracklet23->st2Z ) + tracklet23->st2U;
      double posv = tracklet23->st2Vsl * ( z_plane[5] - tracklet23->st2Z ) + tracklet23->st2V;
      double posy = tracklet23->ty * ( z_plane[3] - tracklet23->st2Z ) + tracklet23->st2Y;
	    
      
        Tracklet tracklet_best[2];
        for(int i = 0; i < 1; ++i) //for two station-1 chambers //WPM edited so that it only does one chamber.  I think that's all we're using in our simulation...
        {

	  bool validTrackFound = false; //WPM
	  int pxSlices[8] = {1, 3, 5, 7, 9, 11, 13, 15}; //This is another scan outward from a linear extrapolation from station 2+3.  We don't know the momentum or charge at this point!  The range here does limit the pz range that you can find to some extent, but I can get below 10 GeV with this range
	  for(int pxs = 0; pxs < 8; pxs++){ //WPM
	    if(validTrackFound) continue; //WPM potentially controversial.  Trying to get out of px window loop
	    
	    int charges[2] = {-1,1}; //WPM
	    for(int ch = 0; ch < 2; ch++){ //WPM
	      
	      if(validTrackFound) continue; //WPM potentially controversial.  Trying to get out of px window loop.  Finding higher pz tracks takes less time
	      
	      trackletsInStSlimX[0][0].clear();
	      trackletsInStSlimU[0][0].clear();
	      trackletsInStSlimV[0][0].clear();

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

	      (*tracklet23).setCharge(charges[ch]); //WPM
	      double expXZSlope = ((*tracklet23).tx - 0.002 * charges[ch]*pxSlices[pxs]);
	      
	      _timers["global_st1"]->restart();
	      if(!TRACK_DISPLACED){
		buildTrackletsInStation(i+1, 0, pos_exp, window);
	      }
	      if(TRACK_DISPLACED){
		pos_exp[0] = posx+charges[ch]*pxSlices[pxs];
		//get expected U and V positions in station 1
		double testSt1Upos = p_geomSvc->getCostheta(1) * ( expXZSlope * (z_plane[1] - z_plane[3]) + pos_exp[0]) + p_geomSvc->getSintheta(1) * (tracklet23->ty * z_plane[1] + tracklet23->y0);
		double testSt1Vpos = p_geomSvc->getCostheta(5) * ( expXZSlope * (z_plane[5] - z_plane[3]) + pos_exp[0]) + p_geomSvc->getSintheta(5) * (tracklet23->ty * z_plane[5] + tracklet23->y0);
		pos_exp[1] = testSt1Upos;
		pos_exp[2] = testSt1Vpos;
		window[0] = 1.25;
		window[1] = 1.5;
		window[2] = 1.5;
		buildTrackletsInStationSlim(i+1, 0, pos_exp, window);
		buildTrackletsInStationSlimU(i+1, 0, pos_exp, window);
		buildTrackletsInStationSlimV(i+1, 0, pos_exp, window);
		getNum1Combos();
		bool doTight = false;
		if(num1XCombos*num1UCombos*num1VCombos > 3000) doTight = true; //Some station 1 pileup mitigation!
		if(!(buildTrackletsInStation1_NEW(i+1, 0, expXZSlope, (*tracklet23).ty, (*tracklet23).y0, doTight, pos_exp, window))) continue; //Find station 1 hit combinations in the relevant window
	      }
	      _timers["global_st1"]->stop();
	      
	      if(_timers["global_st1"]->get_accumulated_time()/1000. > 10.) return;
	      
	      _timers["global_link"]->restart();
	      int tracklet_counter = 0; //WPM
	      
#ifdef _DEBUG_PATRICK
	      LogInfo("the size of trackletsInSt[0] is "<<trackletsInSt[0].size());
#endif
	      int trackletCounter = 0;
	      for(std::list<Tracklet>::iterator tracklet1 = trackletsInSt[0].begin(); tracklet1 != trackletsInSt[0].end(); ++tracklet1)
		{ //loop over the potential station 1 tracklets that we found
#ifdef _DEBUG_ON
		  LogInfo("a new station 1 tracklet global tracks displaced.  trackletCounter = "<<trackletCounter);
		  tracklet1->print();
		  trackletCounter++;
#endif
		  
		  tracklet_counter++; //WPM
#ifdef _DEBUG_ON
		  LogInfo("With this station 1 track:");
		  tracklet1->print();
#endif
		  
		  Tracklet tracklet_global = (*tracklet23) * (*tracklet1);
		  tracklet_global.setCharge(charges[ch]); //WPM      
		  fitTracklet(tracklet_global);
		  if(!COARSE_MODE)
		    {
		      //I don't actually know if this is necessary anymore
		      resolveLeftRight(tracklet_global, 75.);
		      resolveLeftRight(tracklet_global, 150.);
		      resolveSingleLeftRight(tracklet_global);
		    }
		  
		  ///Remove bad hits if needed
		  //removeBadHits(tracklet_global);
		  
#ifdef _DEBUG_ON
		  LogInfo("removed bad hits in global tracks displaced");
#endif
		  
		  //Most basic cuts
		  if(!acceptTracklet(tracklet_global)) continue;
		  
#ifdef _DEBUG_ON
		  LogInfo("accepted tracklet global tracks displaced");
#endif
		  
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
		  //if(tracklet_best_prob.hits.size()==6) validTrackFound = true;
		}
	      else if(enable_KF && NOT_DISPLACED && tracklet_best_vtx.isValid() > 0) //otherwise select the one with best vertex chisq, TODO: maybe add a z-vtx constraint
		{
		  tracklet_best[i] = tracklet_best_vtx;
		  //if(tracklet_best_vtx.hits.size()==6) validTrackFound = true;
		}
	      else if(tracklet_best_prob.isValid() > 0) //then fall back to the default only choice
		{
		  tracklet_best[i] = tracklet_best_prob;
		  //if(tracklet_best_prob.hits.size()==6) validTrackFound = true;
		}
	    }
	    if(tracklet_best_prob.chisq < 1.) break; //MIGHT NEED TO RECONSIDER THIS ACTUALLY
	  }
	  if(tracklet_best_prob.chisq < 1.) break; //SIMILARLY THIS MIGHT BE TOO TIGHT
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
    
#ifdef _DEBUG_PATRICK
    LogInfo("end of buildGlobalTracksDisplaced");
#endif
    
    return;
}

void KalmanFastTracking_NEW_2::resolveLeftRight(Tracklet& tracklet, double threshold)
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

void KalmanFastTracking_NEW_2::resolveSingleLeftRight(Tracklet& tracklet)
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

void KalmanFastTracking_NEW_2::removeBadHits(Tracklet& tracklet)
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

void KalmanFastTracking_NEW_2::resolveLeftRight(SRawEvent::hit_pair hpair, int& LR1, int& LR2)
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

void KalmanFastTracking_NEW_2::buildTrackletsInStation(int stationID, int listID, double* pos_exp, double* window)
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
    //std::cout<<"NUMBER of tracklets before resize = "<<trackletsInSt[listID].size()<<std::endl; //WPM
    if(trackletsInSt[listID].size() > 1000)
    {
        trackletsInSt[listID].sort();
        trackletsInSt[listID].resize(1000);
    }
}




bool KalmanFastTracking_NEW_2::buildTrackletsInStation1_NEW(int stationID, int listID, double expXZSlope, double expYSlope, double y0, bool tight, double* pos_exp, double* window)
{
#ifdef _DEBUG_ON
  LogInfo("Building tracklets in station " << stationID);
#endif

#ifdef _DEBUG_RES
  LogInfo("In buildTrackletsInStation1_NEW");
#endif


  double slopeComparisonSt1 = (tight ? 0.05 : m_slopeComparisonSt1);
  
  double testTU = p_geomSvc->getCostheta(1)*expXZSlope + p_geomSvc->getSintheta(1)*expYSlope;
  double testTV = p_geomSvc->getCostheta(5)*expXZSlope + p_geomSvc->getSintheta(5)*expYSlope;
  
  double expX0 = -1*expXZSlope*z_plane[3] + pos_exp[0];
  
  bool st1TrackletFound = false;
  
  //actuall ID of the tracklet lists
  int sID = stationID - 1;

  //We have the hit combinations in station 1 X, U, and V wires separately at this point
  
  for(std::list<Tracklet>::iterator trackletX = trackletsInStSlimX[0][0].begin(); trackletX != trackletsInStSlimX[0][0].end(); ++trackletX){
    
    if(trackletX->hits.size() < 2 && tight) continue;
    
#ifdef _DEBUG_RES
    std::cout<<"expected x slope = "<<expXZSlope<<std::endl;
    trackletX->print();
#endif
    
    if(trackletX->hits.size() == 2){ //if there are 2 hits, we should be able to determine the hit signs and roughly check if the hit combination gives the right expected slopes
      double bestSlopeDiffX = 1.0;
      double slopeDiffX = 1.0;
      int bestTrackletX = 5;
      double trackletXslope = 1.0;
      int nValidXSlopes = 0;
      for(int t3 = 0; t3 < trackletX->possibleXLines.size(); t3++){
	slopeDiffX = trackletX->possibleXLines.at(t3).slopeX - expXZSlope;
#ifdef _DEBUG_RES
	std::cout<<"x slope diff = "<<slopeDiffX<<std::endl;
#endif
	if(std::abs(slopeDiffX)<slopeComparisonSt1) nValidXSlopes++;
	if(std::abs(slopeDiffX)<slopeComparisonSt1 && std::abs(slopeDiffX)<std::abs(bestSlopeDiffX)){
	  bestSlopeDiffX = slopeDiffX;
	  bestTrackletX = t3;
	  trackletXslope = trackletX->possibleXLines.at(t3).slopeX;
	}
      }
      
      if(bestTrackletX > 4){
	continue;
      } else{
	for(std::list<SignedHit>::iterator hit1 = trackletX->hits.begin(); hit1 != trackletX->hits.end(); ++ hit1)
	  {
	    if(hit1->hit.detectorID == 4){
	      if(bestTrackletX == 0){
		hit1->sign = 1;
	      }
	      if(bestTrackletX == 1){
		hit1->sign = -1;
	      }
	      if(bestTrackletX == 2){
		hit1->sign = 1;
	      }
	      if(bestTrackletX == 3){
		hit1->sign = -1;
	      }
	    }
	    if(hit1->hit.detectorID == 3){
	      if(bestTrackletX == 0){
		hit1->sign = 1;
	      }
	      if(bestTrackletX == 1){
		hit1->sign = 1;
	      }
	      if(bestTrackletX == 2){
		hit1->sign = -1;
	      }
	      if(bestTrackletX == 3){
		hit1->sign = -1;
	      }
	    }
	  }
      }
    }
    
    for(std::list<Tracklet>::iterator trackletU = trackletsInStSlimU[0][0].begin(); trackletU != trackletsInStSlimU[0][0].end(); ++trackletU){
      
      if(trackletU->hits.size() < 2 && tight) continue;
      
#ifdef _DEBUG_RES
      std::cout<<"expected u slope = "<<testTU<<std::endl;
      trackletU->print();
#endif
      
      
      if(trackletU->hits.size() == 2){
	double bestSlopeDiffU = 1.0;
	double slopeDiffU = 1.0;
	int bestTrackletU = 5;
	double trackletUslope = 1.0;
	int nValidUSlopes = 0;
	for(int t3 = 0; t3 < trackletU->possibleULines.size(); t3++){
	  slopeDiffU = trackletU->possibleULines.at(t3).slopeU - testTU;
#ifdef _DEBUG_RES
	  std::cout<<"u slope diff = "<<slopeDiffU<<std::endl;
#endif
	  
	  if(std::abs(slopeDiffU)<slopeComparisonSt1) nValidUSlopes++;
	  if(std::abs(slopeDiffU)<slopeComparisonSt1 && std::abs(slopeDiffU)<std::abs(bestSlopeDiffU)){
	    bestSlopeDiffU = slopeDiffU;
	    bestTrackletU = t3;
	    trackletUslope = trackletU->possibleULines.at(t3).slopeU;
	  }
	}
	
	if(bestTrackletU > 4){
	  continue;
	} else{
	  for(std::list<SignedHit>::iterator hit1 = trackletU->hits.begin(); hit1 != trackletU->hits.end(); ++hit1)
	    {
	      if(hit1->hit.detectorID == 1){
		if(bestTrackletU == 0){
		  hit1->sign = 1;
		}
		if(bestTrackletU == 1){
		  hit1->sign = -1;
		}
		if(bestTrackletU == 2){
		  hit1->sign = 1;
		}
		if(bestTrackletU == 3){
		  hit1->sign = -1;
		}
	      }
	      if(hit1->hit.detectorID == 2){
		if(bestTrackletU == 0){
		  hit1->sign = 1;
		}
		if(bestTrackletU == 1){
		  hit1->sign = 1;
		}
		if(bestTrackletU == 2){
		  hit1->sign = -1;
		}
		if(bestTrackletU == 3){
		  hit1->sign = -1;
		}
	      }
	    }
	}
      }
      
      
      for(std::list<Tracklet>::iterator trackletV = trackletsInStSlimV[0][0].begin(); trackletV != trackletsInStSlimV[0][0].end(); ++trackletV){
	
	if(trackletV->hits.size() < 2 && tight) continue;
	
#ifdef _DEBUG_RES
	std::cout<<"expected v slope = "<<testTV<<std::endl;
	trackletV->print();
#endif
	
	
	if(trackletV->hits.size() == 2){
	  double bestSlopeDiffV = 1.0;
	  double slopeDiffV = 1.0;
	  int bestTrackletV = 5;
	  double trackletVslope = 1.0;
	  int nValidVSlopes = 0;
	  for(int t3 = 0; t3 < trackletV->possibleVLines.size(); t3++){
	    slopeDiffV = trackletV->possibleVLines.at(t3).slopeV - testTV;
#ifdef _DEBUG_RES
	    std::cout<<"v slope diff = "<<slopeDiffV<<std::endl;
#endif
	    
	    if(std::abs(slopeDiffV)<slopeComparisonSt1) nValidVSlopes++;
	    if(std::abs(slopeDiffV)<slopeComparisonSt1 && std::abs(slopeDiffV)<std::abs(bestSlopeDiffV)){
	      bestSlopeDiffV = slopeDiffV;
	      bestTrackletV = t3;
	      trackletVslope = trackletV->possibleVLines.at(t3).slopeV;
	    }
	  }
	  
	  if(bestTrackletV > 4){
	    continue;
	  } else{
	    for(std::list<SignedHit>::iterator hit1 = trackletV->hits.begin(); hit1 != trackletV->hits.end(); ++ hit1)
	      {
		if(hit1->hit.detectorID == 5){
		  if(bestTrackletV == 0){
		    hit1->sign = 1;
		  }
		  if(bestTrackletV == 1){
		    hit1->sign = -1;
		  }
		  if(bestTrackletV == 2){
		    hit1->sign = 1;
		  }
		  if(bestTrackletV == 3){
		    hit1->sign = -1;
		  }
		}
		if(hit1->hit.detectorID == 6){
		  if(bestTrackletV == 0){
		    hit1->sign = 1;
		  }
		  if(bestTrackletV == 1){
		    hit1->sign = 1;
		  }
		  if(bestTrackletV == 2){
		    hit1->sign = -1;
		  }
		  if(bestTrackletV == 3){
		    hit1->sign = -1;
		  }
		}
	      }
	  }
	}

	if( trackletV->hits.size() < 2 && trackletU->hits.size() < 2 && trackletX->hits.size() < 2 ) return false; //We require at least 4 out of 6 possible hits
	
	Tracklet tracklet_new_Station1;
	
	tracklet_new_Station1 = (*trackletX) + (*trackletU) + (*trackletV);
	tracklet_new_Station1.stationID = stationID;
	
	tracklet_new_Station1.y0 = y0;
	tracklet_new_Station1.ty = expYSlope;
	
	tracklet_new_Station1.x0 = expX0;
	tracklet_new_Station1.tx = expXZSlope;
	
	
	tracklet_new_Station1.sortHits();
#ifdef _DEBUG_RES
	std::cout<<"station 1 chisq using no drift distance = "<<tracklet_new_Station1.calcChisq_noDrift()<<std::endl;
	tracklet_new_Station1.print();
#endif	
	
	if(tracklet_new_Station1.calcChisq_noDrift() > 300 || isnan(tracklet_new_Station1.calcChisq_noDrift()) ) continue;
	fitTracklet(tracklet_new_Station1);

	//hit sign assignment for station 1
	for(std::list<SignedHit>::iterator hit1 = tracklet_new_Station1.hits.begin(); hit1 != tracklet_new_Station1.hits.end(); ++hit1){
	  if(hit1->sign == 0){
	    hit1->sign = 1;
	    double dcaPlus = tracklet_new_Station1.calcChisq();
	    hit1->sign = -1;
	    double dcaMinus = tracklet_new_Station1.calcChisq();
	    if(std::abs(dcaPlus) < std::abs(dcaMinus)){
	      hit1->sign = 1;
	    } else{
	      hit1->sign = -1;
	    }
	  }
	}

	fitTracklet(tracklet_new_Station1);
	
#ifdef _DEBUG_RES
	LogInfo("FINAL station 1 fitting");
	tracklet_new_Station1.print();
#endif
	
	if(tracklet_new_Station1.chisq < 20.){
	  trackletsInSt[listID].push_back(tracklet_new_Station1);
	}
	
      }
    }
  }
  
  //Reduce the tracklet list and add dummy hits
  //reduceTrackletList(trackletsInSt[listID]);
  for(std::list<Tracklet>::iterator iter = trackletsInSt[listID].begin(); iter != trackletsInSt[listID].end(); ++iter)
    {
      iter->addDummyHits();
    }

  if(trackletsInSt[listID].size() > 0) return true;
}


//This function finds valid X hits combinations
void KalmanFastTracking_NEW_2::buildTrackletsInStationSlim(int stationID, int listID, double* pos_exp, double* window)
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
    }
    else
    {
        //Note that in pos_exp[], index 0 stands for X, index 1 stands for U, index 2 stands for V
        pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0], pos_exp[0], window[0]);
    }

#ifdef _DEBUG_ON
    LogInfo("Hit pairs in this event: ");
    for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_X.begin(); iter != pairs_X.end(); ++iter) LogInfo("X :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
#endif

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

      if(xiter->first >= 0){
	if(stationID == 3){
	  tracklet_new.st2X = hitAll[xiter->first].pos;
	  tracklet_new.st2Z = z_plane[hitAll[xiter->first].detectorID];
	}
	if(stationID == 4 || stationID == 5){
	  tracklet_new.st3X = hitAll[xiter->first].pos;
	  tracklet_new.st3Z = z_plane[hitAll[xiter->first].detectorID];
	}
      } else if(xiter->second >= 0){
	if(stationID == 3){
	  tracklet_new.st2X = hitAll[xiter->second].pos;
	  tracklet_new.st2Z = z_plane[hitAll[xiter->second].detectorID];
	}
	if(stationID == 4 || stationID == 5){
	  tracklet_new.st3X = hitAll[xiter->second].pos;
	  tracklet_new.st3Z = z_plane[hitAll[xiter->second].detectorID];
	}
      }
      
      if(!OLD_TRACKING){
	tracklet_new.getSlopesX(hitAll[xiter->first], hitAll[xiter->second]); //Here, we find the four possible X-Z lines
      }
      
      tracklet_new.sortHits();
#ifdef _DEBUG_ON
      tracklet_new.print();
#endif
      
      trackletsInStSlimX[listID][0].push_back(tracklet_new);
    }

}


void KalmanFastTracking_NEW_2::buildTrackletsInStationSlimU(int stationID, int listID, double* pos_exp, double* window)
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
	if(stationID == 4 || stationID == 5){
	  pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2]);
	}
	else{
	  pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1]);
	}
      }
    else
      {
	if(stationID == 4 || stationID == 5){
	  pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2], pos_exp[2], window[2]);
	}
	else{
	  //Note that in pos_exp[], index 0 stands for X, index 1 stands for U, index 2 stands for V
	  pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1], pos_exp[1], window[1]);
	}
      }
    
#ifdef _DEBUG_ON
    LogInfo("Hit pairs in this event: ");
    for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_U.begin(); iter != pairs_U.end(); ++iter) LogInfo("U :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
#endif

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

      if(xiter->first >= 0){
	if(stationID == 3){
	  tracklet_new.st2U = hitAll[xiter->first].pos;
	  tracklet_new.st2Z = z_plane[hitAll[xiter->first].detectorID];
	}
	if(stationID == 4 || stationID == 5){
	  tracklet_new.st3U = hitAll[xiter->first].pos;
	  tracklet_new.st3Z = z_plane[hitAll[xiter->first].detectorID];
	}
      } else if(xiter->second >= 0){
	if(stationID == 3){
	  tracklet_new.st2U = hitAll[xiter->second].pos;
	  tracklet_new.st2Z = z_plane[hitAll[xiter->second].detectorID];
	}
	if(stationID == 4 || stationID == 5){
	  tracklet_new.st3U = hitAll[xiter->second].pos;
	  tracklet_new.st3Z = z_plane[hitAll[xiter->second].detectorID];
	}
      }
      
      if(!OLD_TRACKING){
	tracklet_new.getSlopesU(hitAll[xiter->first], hitAll[xiter->second]); //Here, we find the four possible X-Z lines
      }
      
      tracklet_new.sortHits();
#ifdef _DEBUG_ON
      tracklet_new.print();
#endif
      
      trackletsInStSlimU[listID][0].push_back(tracklet_new);
    }

}



void KalmanFastTracking_NEW_2::buildTrackletsInStationSlimV(int stationID, int listID, double* pos_exp, double* window)
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
      if(stationID == 4 || stationID == 5){
	pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1]);
      }
      else{
        pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2]);
      }
    }
    else
    {
      if(stationID == 4 || stationID == 5){
	pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1], pos_exp[1], window[1]);
      }
      else{
        //Note that in pos_exp[], index 0 stands for X, index 1 stands for U, index 2 stands for V
        pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2], pos_exp[2], window[2]);
      }
    }

#ifdef _DEBUG_ON
    LogInfo("Hit pairs in this event: ");
    for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_V.begin(); iter != pairs_V.end(); ++iter) LogInfo("V :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
#endif

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

      if(xiter->first >= 0){
	if(stationID == 3){
	  tracklet_new.st2V = hitAll[xiter->first].pos;
	  tracklet_new.st2Z = z_plane[hitAll[xiter->first].detectorID];
	}
	if(stationID == 4 || stationID == 5){
	  tracklet_new.st3V = hitAll[xiter->first].pos;
	  tracklet_new.st3Z = z_plane[hitAll[xiter->first].detectorID];
	}
      } else if(xiter->second >= 0){
	if(stationID == 3){
	  tracklet_new.st2V = hitAll[xiter->second].pos;
	  tracklet_new.st2Z = z_plane[hitAll[xiter->second].detectorID];
	}
	if(stationID == 4 || stationID == 5){
	  tracklet_new.st3V = hitAll[xiter->second].pos;
	  tracklet_new.st3Z = z_plane[hitAll[xiter->second].detectorID];
	}
      }
      
      if(!OLD_TRACKING){
	tracklet_new.getSlopesV(hitAll[xiter->first], hitAll[xiter->second]); //Here, we find the four possible X-Z lines
      }
      
      tracklet_new.sortHits();
#ifdef _DEBUG_ON
      tracklet_new.print();
#endif
      
      trackletsInStSlimV[listID][0].push_back(tracklet_new);
    }

}



bool KalmanFastTracking_NEW_2::acceptTracklet(Tracklet& tracklet)
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
    /* //This algorithm ignores hodomasking requirements for now.  It was only losing efficiency from what I can tell
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
    }*/ //WPM

#ifdef _DEBUG_ON
        LogInfo("Made it through various checks");
#endif
    
    //If everything is fine ...
#ifdef _DEBUG_ON
    LogInfo("AcceptTracklet!!!");
#endif
    return true;
}

bool KalmanFastTracking_NEW_2::hodoMask(Tracklet& tracklet)
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

bool KalmanFastTracking_NEW_2::muonID_search(Tracklet& tracklet)
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

bool KalmanFastTracking_NEW_2::muonID_comp(Tracklet& tracklet)
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

bool KalmanFastTracking_NEW_2::muonID_hodoAid(Tracklet& tracklet)
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

void KalmanFastTracking_NEW_2::buildPropSegments()
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
        //std::cout << "superID: " << superIDs[i+5][0] << ", " << superIDs[i+5][1] << std::endl;
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


int KalmanFastTracking_NEW_2::fitTracklet(Tracklet& tracklet)
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

int KalmanFastTracking_NEW_2::reduceTrackletList(std::list<Tracklet>& tracklets)
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

void KalmanFastTracking_NEW_2::getExtrapoWindowsInSt1(Tracklet& tracklet, double* pos_exp, double* window, int st1ID)
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

void KalmanFastTracking_NEW_2::getSagittaWindowsInSt1(Tracklet& tracklet, double* pos_exp, double* window, int st1ID)
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

	//std::cout<<"Test sagitta window thing.  The detectorID is "<<detectorID<<" the idx is "<<idx<<" the z_st1 is "<<z_st1<<std::endl; //WPM
	//std::cout<<"More sagitta.  pos_exp is "<<0.5*(p_max + p_min)<<" and window = "<<0.5*(p_max - p_min)<<std::endl; //WPM
        pos_exp[idx] = 0.5*(p_max + p_min);
        window[idx]  = 0.5*(p_max - p_min);
    }
}

void KalmanFastTracking_NEW_2::printAtDetectorBack(int stationID, std::string outputFileName)
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

SRecTrack KalmanFastTracking_NEW_2::processOneTracklet(Tracklet& tracklet)
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

bool KalmanFastTracking_NEW_2::fitTrack(KalmanTrack& kmtrk)
{
    if(kmtrk.getNodeList().empty()) return false;

    if(kmfitter->processOneTrack(kmtrk) == 0)
    {
        return false;
    }
    kmfitter->updateTrack(kmtrk);

    return true;
}

void KalmanFastTracking_NEW_2::resolveLeftRight(KalmanTrack& kmtrk)
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

void KalmanFastTracking_NEW_2::printTimers() {
	std::cout <<"KalmanFastTracking_NEW_2::printTimers: " << std::endl;
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

void KalmanFastTracking_NEW_2::chi2fit(int n, double x[], double y[], double& a, double& b)
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



//For the case when both st2 and st3 hit combos have 2 hits
bool KalmanFastTracking_NEW_2::compareTrackletsSlim(Tracklet& tracklet2, Tracklet& tracklet3, int pass, double slopeComparison, double windowSize)
{
  if( !( tracklet2.hits.size() == 2 && tracklet3.hits.size() == 2 ) ) return false;
  
  double slopeWindow = 0.005;
  double extrapoWindow = 5.0;

  if(pass == 2){
    slopeWindow = 0.0025;
    extrapoWindow = 2.5;
  }
  
  //Here we will compare the possible X-Z slopes within the station 2 and station 3 tracklets
  Tracklet::linedef line2X;
  Tracklet::linedef line3X;
  Tracklet::linedef line2X_v2;
  Tracklet::linedef line3X_v2;

  int choiceOfT2 = 5;
  int choiceOfT3 = 5;
  int choiceOfT2_v2 = 5;
  int choiceOfT3_v2 = 5;
  
  double st3pos = (tracklet3.getHit(0).hit.pos + tracklet3.getHit(1).hit.pos)/2;
  double st3posZ = (z_plane[tracklet3.getHit(0).hit.detectorID] + z_plane[tracklet3.getHit(1).hit.detectorID])/2;

  //It is rare, but sometimes, you will have slopes that match coincidentally.  Therefore, I keep track of best two combinations.  This seems to be sufficienct
  double slopeComp = 1.0;
  double bestExtrapComp = 100.;

  double slopeComp2 = 1.0;
  double bestExtrapComp2 = 100.;

  double bestExtrapMatch = 100.;
  
  double secondSlope = 1.1;
  for(unsigned int t2 = 0; t2 < tracklet2.possibleXLines.size(); t2++){

    double st2ExtrapMatch = tracklet2.possibleXLines.at(t2).slopeX*(st3posZ - tracklet2.possibleXLines.at(t2).initialZ) + tracklet2.possibleXLines.at(t2).initialX - st3pos;
    if( std::abs(st2ExtrapMatch) < windowSize && std::abs(st2ExtrapMatch) < bestExtrapMatch ){ //Is this extrapolation valid?

      for(unsigned int t3 = 0; t3 < tracklet3.possibleXLines.size(); t3++){

	if(std::abs(tracklet3.possibleXLines.at(t3).slopeX - tracklet2.possibleXLines.at(t2).slopeX) < slopeComp && std::abs(tracklet3.possibleXLines.at(t3).slopeX - tracklet2.possibleXLines.at(t2).slopeX) < slopeComparison){	

	  slopeComp = std::abs(tracklet3.possibleXLines.at(t3).slopeX - tracklet2.possibleXLines.at(t2).slopeX);
	  bestExtrapComp = st2ExtrapMatch;
	  bestExtrapMatch = std::abs(st2ExtrapMatch);
	  line2X = tracklet2.possibleXLines.at(t2);
	  line3X = tracklet3.possibleXLines.at(t3);
	  choiceOfT2 = t2;
	  choiceOfT3 = t3;
	  
	}
      }
    }
  }    
  
#ifdef _DEBUG_RES
  LogInfo("slopeComp = "<<slopeComp);
  LogInfo("bestExtrapComp = "<<bestExtrapComp);
  LogInfo("extrapolation is "<<line2X.slopeX*(line3X.initialZ - line2X.initialZ) + line2X.initialX<<"; based on slope "<<line2X.slopeX<<", diff "<<line3X.initialZ - line2X.initialZ<<", and initialX "<<line2X.initialX);
  LogInfo("to be compared with st3 position "<<line3X.initialX<<".  diff is "<<std::abs(line2X.slopeX*(line3X.initialZ - line2X.initialZ) + line2X.initialX - line3X.initialX));
#endif

  if(slopeComp > slopeComparison) return false;
  if(bestExtrapComp > windowSize) return false;


  tracklet2.acceptedXLine2 = line2X;
  tracklet3.acceptedXLine3 = line3X;
  //Give the station 2 and station 3 tracklets the same tx and ty value.  I could get an X0 and Y0 extrapolation, but that doesn't seem to be strictly necessary.  The X0 and Y0 values are found in the fittracklet function for the combined station 2 + station 3 tracklet
  tracklet2.tx = (line2X.slopeX + line3X.slopeX)/2;
  tracklet3.tx = (line2X.slopeX + line3X.slopeX)/2;

  tracklet2.st2X = line2X.initialX;
  tracklet3.st2X = line3X.initialX; //yes, I know it looks like a typo, but I set the st2x for tracklet 3 here.  It's used elsewhere
  tracklet3.st3X = line3X.initialX;

  //extract slope measurement using the long lever-arm between station 2 and 3
  double newSlopeX = (tracklet3.possibleXLines.at(choiceOfT3).initialX - tracklet2.possibleXLines.at(choiceOfT2).initialX)/(tracklet3.possibleXLines.at(choiceOfT3).initialZ - tracklet2.possibleXLines.at(choiceOfT2).initialZ);

  tracklet2.st2Xsl = newSlopeX;
  tracklet3.st2Xsl = newSlopeX;


  for(std::list<SignedHit>::iterator hit1 = tracklet2.hits.begin(); hit1 != tracklet2.hits.end(); ++hit1)
    {
      if(hit1->hit.detectorID == 15){
	if(choiceOfT2 == 0){
	  hit1->sign = 1;
	}
	if(choiceOfT2 == 1){
	  hit1->sign = -1;
	}
	if(choiceOfT2 == 2){
	  hit1->sign = 1;
	}
	if(choiceOfT2 == 3){
	  hit1->sign = -1;
	}
      } else{
	if(choiceOfT2 == 0){
	  hit1->sign = 1;
	}
	if(choiceOfT2 == 1){
	  hit1->sign = 1;
	}
	if(choiceOfT2 == 2){
	  hit1->sign = -1;
	}
	if(choiceOfT2 == 3){
	  hit1->sign = -1;
	}
      }
    }

    for(std::list<SignedHit>::iterator hit1 = tracklet3.hits.begin(); hit1 != tracklet3.hits.end(); ++hit1)
    {
      if(hit1->hit.detectorID == 21 || hit1->hit.detectorID == 27){
	if(choiceOfT3 == 0){
	  hit1->sign = 1;
	}
	if(choiceOfT3 == 1){
	  hit1->sign = -1;
	}
	if(choiceOfT3 == 2){
	  hit1->sign = 1;
	}
	if(choiceOfT3 == 3){
	  hit1->sign = -1;
	}
      } else{
	if(choiceOfT3 == 0){
	  hit1->sign = 1;
	}
	if(choiceOfT3 == 1){
	  hit1->sign = 1;
	}
	if(choiceOfT3 == 2){
	  hit1->sign = -1;
	}
	if(choiceOfT3 == 3){
	  hit1->sign = -1;
	}
      }
    }

  
  return true;
  
}


//For the case when st2 X combo has 2 hits and st3 has 1, or vise versa
bool KalmanFastTracking_NEW_2::compareTrackletsSlim_3hits(Tracklet& tracklet2, Tracklet& tracklet3, int pass, double slopeComparison, double windowSize)
{
  
  if( ! ( (tracklet2.hits.size() == 2 && tracklet3.hits.size() == 1 ) || (tracklet2.hits.size() == 1 && tracklet3.hits.size() == 2 ) ) ) return false;
  
  double bestMatch = 100.;
  if(tracklet2.hits.size() == 2){
    
    unsigned int bestT2 = 5;
    for(unsigned int t2 = 0; t2 < tracklet2.possibleXLines.size(); t2++){
      double extrapolation = tracklet2.possibleXLines.at(t2).slopeX*(z_plane[tracklet3.getHit(0).hit.detectorID] - tracklet2.possibleXLines.at(t2).initialZ) + tracklet2.possibleXLines.at(t2).initialX;
      if(std::abs(tracklet3.getHit(0).hit.pos - extrapolation) < bestMatch && std::abs(tracklet3.getHit(0).hit.pos - extrapolation) < windowSize){
	bestMatch = std::abs(tracklet3.getHit(0).hit.pos - extrapolation);
	bestT2 = t2;
      }
    }
    
    if(bestMatch > 99. || bestT2 == 5){
      return false;
    } else{
      
      double newSlopeX = (tracklet3.getHit(0).hit.pos - tracklet2.possibleXLines.at(bestT2).initialX)/(z_plane[tracklet3.getHit(0).hit.detectorID] - tracklet2.possibleXLines.at(bestT2).initialZ);
      
      tracklet2.acceptedXLine2 = tracklet2.possibleXLines.at(bestT2);
      tracklet3.acceptedXLine3 = tracklet2.possibleXLines.at(bestT2);
      tracklet2.tx = newSlopeX;
      tracklet3.tx = newSlopeX;
      
      tracklet2.st2Xsl = newSlopeX;
      tracklet3.st2Xsl = newSlopeX;


      for(std::list<SignedHit>::iterator hit1 = tracklet2.hits.begin(); hit1 != tracklet2.hits.end(); ++hit1)
	{
	  if(hit1->hit.detectorID == 15){
	    if(bestT2 == 0){
	      hit1->sign = 1;
	    }
	    if(bestT2 == 1){
	      hit1->sign = -1;
	    }
	    if(bestT2 == 2){
	      hit1->sign = 1;
	    }
	    if(bestT2 == 3){
	      hit1->sign = -1;
	    }
	  } else{
	    if(bestT2 == 0){
	      hit1->sign = 1;
	    }
	    if(bestT2 == 1){
	      hit1->sign = 1;
	    }
	    if(bestT2 == 2){
	      hit1->sign = -1;
	    }
	    if(bestT2 == 3){
	      hit1->sign = -1;
	    }
	  }
	}


      
      return true;
    }
  }
  
  if(tracklet3.hits.size() == 2){
    
    unsigned int bestT3 = 5;
    for(unsigned int t3 = 0; t3 < tracklet3.possibleXLines.size(); t3++){
      double extrapolation = tracklet3.possibleXLines.at(t3).slopeX*(z_plane[tracklet2.getHit(0).hit.detectorID] - tracklet3.possibleXLines.at(t3).initialZ) + tracklet3.possibleXLines.at(t3).initialX;
      if(std::abs(tracklet2.getHit(0).hit.pos - extrapolation) < bestMatch && std::abs(tracklet2.getHit(0).hit.pos - extrapolation) < windowSize){
	bestMatch = std::abs(tracklet2.getHit(0).hit.pos - extrapolation);
	bestT3 = t3;
      }
    }
    
    if(bestMatch > 99. || bestT3 == 5){
      return false;
    } else{
      
      double newSlopeX = (tracklet2.getHit(0).hit.pos - tracklet3.possibleXLines.at(bestT3).initialX)/(z_plane[tracklet2.getHit(0).hit.detectorID] - tracklet3.possibleXLines.at(bestT3).initialZ);
      
      tracklet2.acceptedXLine2 = tracklet3.possibleXLines.at(bestT3);
      tracklet3.acceptedXLine3 = tracklet3.possibleXLines.at(bestT3);
      tracklet2.tx = newSlopeX;
      tracklet3.tx = newSlopeX;
      
      tracklet2.st2Xsl = newSlopeX;
      tracklet3.st2Xsl = newSlopeX;

      
      for(std::list<SignedHit>::iterator hit1 = tracklet3.hits.begin(); hit1 != tracklet3.hits.end(); ++hit1)
	{
	  if(hit1->hit.detectorID == 21 || hit1->hit.detectorID == 27){
	    if(bestT3 == 0){
	      hit1->sign = 1;
	    }
	    if(bestT3 == 1){
	  hit1->sign = -1;
	    }
	    if(bestT3 == 2){
	      hit1->sign = 1;
	    }
	    if(bestT3 == 3){
	      hit1->sign = -1;
	    }
	  } else{
	    if(bestT3 == 0){
	      hit1->sign = 1;
	    }
	    if(bestT3 == 1){
	      hit1->sign = 1;
	    }
	    if(bestT3 == 2){
	      hit1->sign = -1;
	    }
	    if(bestT3 == 3){
	      hit1->sign = -1;
	    }
	  }
	}
      
      return true;
      
    }
  }
}


bool KalmanFastTracking_NEW_2::compareTrackletsSlimU(Tracklet& tracklet2, Tracklet& tracklet3, int pass, double slopeComparison, double windowSize)
{
  if( !( tracklet2.hits.size() == 2 && tracklet3.hits.size() == 2 ) ) return false;
  
  double slopeWindow = 0.005;
  double extrapoWindow = 5.0;
  
  if(pass == 2){
    slopeWindow = 0.008;
    extrapoWindow = 12.;
  }
  
  //Here we will compare the possible X-Z slopes within the station 2 and station 3 tracklets
  Tracklet::linedef line2U;
  Tracklet::linedef line3U;
  Tracklet::linedef line2U_v2;
  Tracklet::linedef line3U_v2;

  int choiceOfT2 = 5;
  int choiceOfT3 = 5;
  int choiceOfT2_v2 = 5;
  int choiceOfT3_v2 = 5;
  
  double st3pos = (tracklet3.getHit(0).hit.pos + tracklet3.getHit(1).hit.pos)/2;
  double st3posZ = (z_plane[tracklet3.getHit(0).hit.detectorID] + z_plane[tracklet3.getHit(1).hit.detectorID])/2;
  
  double slopeComp = 1.0;
  double bestExtrapComp = 100.;

  double slopeComp2 = 1.0;
  double bestExtrapComp2 = 100.;

  double bestExtrapMatch = 100.;
  
  double secondSlope = 1.1;
  for(unsigned int t2 = 0; t2 < tracklet2.possibleULines.size(); t2++){

    double st2ExtrapMatch = tracklet2.possibleULines.at(t2).slopeU*(st3posZ - tracklet2.possibleULines.at(t2).initialZ) + tracklet2.possibleULines.at(t2).initialU - st3pos;
    if( std::abs(st2ExtrapMatch) < windowSize && std::abs(st2ExtrapMatch) < bestExtrapMatch ){

      for(unsigned int t3 = 0; t3 < tracklet3.possibleULines.size(); t3++){

	if(std::abs(tracklet3.possibleULines.at(t3).slopeU - tracklet2.possibleULines.at(t2).slopeU) < slopeComp && std::abs(tracklet3.possibleULines.at(t3).slopeU - tracklet2.possibleULines.at(t2).slopeU) < slopeComparison){	

	  slopeComp = std::abs(tracklet3.possibleULines.at(t3).slopeU - tracklet2.possibleULines.at(t2).slopeU);
	  bestExtrapComp = st2ExtrapMatch;
	  bestExtrapMatch = std::abs(st2ExtrapMatch);
	  line2U = tracklet2.possibleULines.at(t2);
	  line3U = tracklet3.possibleULines.at(t3);
	  choiceOfT2 = t2;
	  choiceOfT3 = t3;

	}
      }
    }
  }    


  
#ifdef _DEBUG_RES
  LogInfo("slopeComp = "<<slopeComp);
  LogInfo("bestExtrapComp = "<<bestExtrapComp);
  LogInfo("extrapolation is "<<line2U.slopeU*(line3U.initialZ - line2U.initialZ) + line2U.initialU<<"; based on slope "<<line2U.slopeU<<", diff "<<line3U.initialZ - line2U.initialZ<<", and initialU "<<line2U.initialU);
  LogInfo("to be compared with st3 position "<<line3U.initialU<<".  diff is "<<std::abs(line2U.slopeU*(line3U.initialZ - line2U.initialZ) + line2U.initialU - line3U.initialU));
#endif

  if(slopeComp > slopeComparison) return false;
  if(bestExtrapComp > windowSize) return false;
  
  
  tracklet2.acceptedULine2 = line2U;
  tracklet3.acceptedULine3 = line3U;

  tracklet2.st2U = line2U.initialU;
  tracklet3.st2U = line3U.initialU;
  tracklet3.st3U = line3U.initialU;

  double newSlopeU = (tracklet3.possibleULines.at(choiceOfT3).initialU - tracklet2.possibleULines.at(choiceOfT2).initialU)/(tracklet3.possibleULines.at(choiceOfT3).initialZ - tracklet2.possibleULines.at(choiceOfT2).initialZ);
  tracklet2.st2Usl = newSlopeU;
  tracklet3.st2Usl = newSlopeU;


  for(std::list<SignedHit>::iterator hit1 = tracklet2.hits.begin(); hit1 != tracklet2.hits.end(); ++hit1)
    {
      if(hit1->hit.detectorID == 17){
	if(choiceOfT2 == 0){
	  hit1->sign = 1;
	}
	if(choiceOfT2 == 1){
	  hit1->sign = -1;
	}
	if(choiceOfT2 == 2){
	  hit1->sign = 1;
	}
	if(choiceOfT2 == 3){
	  hit1->sign = -1;
	}
      } else{
	if(choiceOfT2 == 0){
	  hit1->sign = 1;
	}
	if(choiceOfT2 == 1){
	  hit1->sign = 1;
	}
	if(choiceOfT2 == 2){
	  hit1->sign = -1;
	}
	if(choiceOfT2 == 3){
	  hit1->sign = -1;
	}
      }
    }

    for(std::list<SignedHit>::iterator hit1 = tracklet3.hits.begin(); hit1 != tracklet3.hits.end(); ++hit1)
    {
      if(hit1->hit.detectorID == 19 || hit1->hit.detectorID == 25){
	if(choiceOfT3 == 0){
	  hit1->sign = 1;
	}
	if(choiceOfT3 == 1){
	  hit1->sign = -1;
	}
	if(choiceOfT3 == 2){
	  hit1->sign = 1;
	}
	if(choiceOfT3 == 3){
	  hit1->sign = -1;
	}
      } else{
	if(choiceOfT3 == 0){
	  hit1->sign = 1;
	}
	if(choiceOfT3 == 1){
	  hit1->sign = 1;
	}
	if(choiceOfT3 == 2){
	  hit1->sign = -1;
	}
	if(choiceOfT3 == 3){
	  hit1->sign = -1;
	}
      }
    }
  
  return true;
  
}

bool KalmanFastTracking_NEW_2::compareTrackletsSlimU_3hits(Tracklet& tracklet2, Tracklet& tracklet3, int pass, double slopeComparison, double windowSize)
{
  
  if( ! ( (tracklet2.hits.size() == 2 && tracklet3.hits.size() == 1 ) || (tracklet2.hits.size() == 1 && tracklet3.hits.size() == 2 ) ) ) return false;
  
  double bestMatch = 100.;
  if(tracklet2.hits.size() == 2){
    
    unsigned int bestT2 = 5;
    for(unsigned int t2 = 0; t2 < tracklet2.possibleULines.size(); t2++){
      double extrapolation = tracklet2.possibleULines.at(t2).slopeU*(z_plane[tracklet3.getHit(0).hit.detectorID] - tracklet2.possibleULines.at(t2).initialZ) + tracklet2.possibleULines.at(t2).initialU;
      if(std::abs(tracklet3.getHit(0).hit.pos - extrapolation) < bestMatch && std::abs(tracklet3.getHit(0).hit.pos - extrapolation) < windowSize){
	bestMatch = std::abs(tracklet3.getHit(0).hit.pos - extrapolation);
	bestT2 = t2;
      }
    }
    
    if(bestMatch > 99. || bestT2 == 5){
      return false;
    } else{
      
      double newSlopeU = (tracklet3.getHit(0).hit.pos - tracklet2.possibleULines.at(bestT2).initialU)/(z_plane[tracklet3.getHit(0).hit.detectorID] - tracklet2.possibleULines.at(bestT2).initialZ);
      
      tracklet2.acceptedULine2 = tracklet2.possibleULines.at(bestT2);
      tracklet3.acceptedULine3 = tracklet2.possibleULines.at(bestT2);
      tracklet2.tx = newSlopeU;
      tracklet3.tx = newSlopeU;
      
      tracklet2.st2Usl = newSlopeU;
      tracklet3.st2Usl = newSlopeU;


      for(std::list<SignedHit>::iterator hit1 = tracklet2.hits.begin(); hit1 != tracklet2.hits.end(); ++hit1)
	{
	  if(hit1->hit.detectorID == 17){
	    if(bestT2 == 0){
	      hit1->sign = 1;
	    }
	    if(bestT2 == 1){
	      hit1->sign = -1;
	    }
	    if(bestT2 == 2){
	      hit1->sign = 1;
	    }
	    if(bestT2 == 3){
	      hit1->sign = -1;
	    }
	  } else{
	    if(bestT2 == 0){
	      hit1->sign = 1;
	    }
	    if(bestT2 == 1){
	      hit1->sign = 1;
	    }
	    if(bestT2 == 2){
	      hit1->sign = -1;
	    }
	    if(bestT2 == 3){
	      hit1->sign = -1;
	    }
	  }
	}

      
      return true;
    }
  }
  
  if(tracklet3.hits.size() == 2){
    
    unsigned int bestT3 = 5;
    for(unsigned int t3 = 0; t3 < tracklet3.possibleULines.size(); t3++){
      double extrapolation = tracklet3.possibleULines.at(t3).slopeU*(z_plane[tracklet2.getHit(0).hit.detectorID] - tracklet3.possibleULines.at(t3).initialZ) + tracklet3.possibleULines.at(t3).initialU;
      if(std::abs(tracklet2.getHit(0).hit.pos - extrapolation) < bestMatch && std::abs(tracklet2.getHit(0).hit.pos - extrapolation) < windowSize){
	bestMatch = std::abs(tracklet2.getHit(0).hit.pos - extrapolation);
	bestT3 = t3;
      }
    }
    
    if(bestMatch > 99. || bestT3 == 5){
      return false;
    } else{
      
      double newSlopeU = (tracklet2.getHit(0).hit.pos - tracklet3.possibleULines.at(bestT3).initialU)/(z_plane[tracklet2.getHit(0).hit.detectorID] - tracklet3.possibleULines.at(bestT3).initialZ);
      
      tracklet2.acceptedULine2 = tracklet3.possibleULines.at(bestT3);
      tracklet3.acceptedULine3 = tracklet3.possibleULines.at(bestT3);
      tracklet2.tx = newSlopeU;
      tracklet3.tx = newSlopeU;
      
      tracklet2.st2Usl = newSlopeU;
      tracklet3.st2Usl = newSlopeU;

      for(std::list<SignedHit>::iterator hit1 = tracklet3.hits.begin(); hit1 != tracklet3.hits.end(); ++hit1)
	{
	  if(hit1->hit.detectorID == 19 || hit1->hit.detectorID == 25){
	    if(bestT3 == 0){
	      hit1->sign = 1;
	    }
	    if(bestT3 == 1){
	  hit1->sign = -1;
	    }
	    if(bestT3 == 2){
	      hit1->sign = 1;
	    }
	    if(bestT3 == 3){
	      hit1->sign = -1;
	    }
	  } else{
	    if(bestT3 == 0){
	      hit1->sign = 1;
	    }
	    if(bestT3 == 1){
	      hit1->sign = 1;
	    }
	    if(bestT3 == 2){
	      hit1->sign = -1;
	    }
	    if(bestT3 == 3){
	      hit1->sign = -1;
	    }
	  }
	}

      
      return true;
      
    }
  }
}



bool KalmanFastTracking_NEW_2::compareTrackletsSlimV(Tracklet& tracklet2, Tracklet& tracklet3, int pass, double slopeComparison, double windowSize)
{
  if( !( tracklet2.hits.size() == 2 && tracklet3.hits.size() == 2 ) ) return false;

  double slopeWindow = 0.005;
  double extrapoWindow = 5.0;
  
  if(pass == 2){
    slopeWindow = 0.008;
    extrapoWindow = 12.;
  }
  
  //Here we will compare the possible X-Z slopes within the station 2 and station 3 tracklets
  Tracklet::linedef line2V;
  Tracklet::linedef line3V;
  Tracklet::linedef line2V_v2;
  Tracklet::linedef line3V_v2;

  int choiceOfT2 = 5;
  int choiceOfT3 = 5;
  int choiceOfT2_v2 = 5;
  int choiceOfT3_v2 = 5;

  double st3pos = (tracklet3.getHit(0).hit.pos + tracklet3.getHit(1).hit.pos)/2;
  double st3posZ = (z_plane[tracklet3.getHit(0).hit.detectorID] + z_plane[tracklet3.getHit(1).hit.detectorID])/2;
  
  double slopeComp = 1.0;
  double bestExtrapComp = 100.;

  double slopeComp2 = 1.0;
  double bestExtrapComp2 = 100.;

  double bestExtrapMatch = 100.;
  
  double secondSlope = 1.1;
  for(unsigned int t2 = 0; t2 < tracklet2.possibleVLines.size(); t2++){

    double st2ExtrapMatch = tracklet2.possibleVLines.at(t2).slopeV*(st3posZ - tracklet2.possibleVLines.at(t2).initialZ) + tracklet2.possibleVLines.at(t2).initialV - st3pos;
    if( std::abs(st2ExtrapMatch) < windowSize && std::abs(st2ExtrapMatch) < bestExtrapMatch ){

      for(unsigned int t3 = 0; t3 < tracklet3.possibleVLines.size(); t3++){

	if(std::abs(tracklet3.possibleVLines.at(t3).slopeV - tracklet2.possibleVLines.at(t2).slopeV) < slopeComp && std::abs(tracklet3.possibleVLines.at(t3).slopeV - tracklet2.possibleVLines.at(t2).slopeV) < slopeComparison){	

	  slopeComp = std::abs(tracklet3.possibleVLines.at(t3).slopeV - tracklet2.possibleVLines.at(t2).slopeV);
	  bestExtrapComp = st2ExtrapMatch;
	  bestExtrapMatch = std::abs(st2ExtrapMatch);
	  line2V = tracklet2.possibleVLines.at(t2);
	  line3V = tracklet3.possibleVLines.at(t3);
	  choiceOfT2 = t2;
	  choiceOfT3 = t3;

	}
      }
    }
  }    



#ifdef _DEBUG_RES
  LogInfo("slopeComp = "<<slopeComp);
  LogInfo("bestExtrapComp = "<<bestExtrapComp);
  LogInfo("extrapolation is "<<line2V.slopeV*(line3V.initialZ - line2V.initialZ) + line2V.initialV<<"; based on slope "<<line2V.slopeV<<", diff "<<line3V.initialZ - line2V.initialZ<<", and initialV "<<line2V.initialV);
  LogInfo("to be compared with st3 position "<<line3V.initialV<<".  diff is "<<std::abs(line2V.slopeV*(line3V.initialZ - line2V.initialZ) + line2V.initialV - line3V.initialV));
#endif

  if(slopeComp > slopeComparison) return false;
  if(bestExtrapComp > windowSize) return false;






  
  tracklet2.acceptedVLine2 = line2V;
  tracklet3.acceptedVLine3 = line3V;
  
  tracklet2.st2V = line2V.initialV;
  tracklet3.st2V = line3V.initialV;
  tracklet3.st3V = line3V.initialV;

  double newSlopeV = (tracklet3.possibleVLines.at(choiceOfT3).initialV - tracklet2.possibleVLines.at(choiceOfT2).initialV)/(tracklet3.possibleVLines.at(choiceOfT3).initialZ - tracklet2.possibleVLines.at(choiceOfT2).initialZ);
    
  tracklet2.st2Vsl = newSlopeV;
  tracklet3.st2Vsl = newSlopeV;

  for(std::list<SignedHit>::iterator hit1 = tracklet2.hits.begin(); hit1 != tracklet2.hits.end(); ++hit1)
    {
      if(hit1->hit.detectorID == 13){
	if(choiceOfT2 == 0){
	  hit1->sign = 1;
	}
	if(choiceOfT2 == 1){
	  hit1->sign = -1;
	}
	if(choiceOfT2 == 2){
	  hit1->sign = 1;
	}
	if(choiceOfT2 == 3){
	  hit1->sign = -1;
	}
      } else{
	if(choiceOfT2 == 0){
	  hit1->sign = 1;
	}
	if(choiceOfT2 == 1){
	  hit1->sign = 1;
	}
	if(choiceOfT2 == 2){
	  hit1->sign = -1;
	}
	if(choiceOfT2 == 3){
	  hit1->sign = -1;
	}
      }
    }

    for(std::list<SignedHit>::iterator hit1 = tracklet3.hits.begin(); hit1 != tracklet3.hits.end(); ++hit1)
    {
      if(hit1->hit.detectorID == 23 || hit1->hit.detectorID == 29){
	if(choiceOfT3 == 0){
	  hit1->sign = 1;
	}
	if(choiceOfT3 == 1){
	  hit1->sign = -1;
	}
	if(choiceOfT3 == 2){
	  hit1->sign = 1;
	}
	if(choiceOfT3 == 3){
	  hit1->sign = -1;
	}
      } else{
	if(choiceOfT3 == 0){
	  hit1->sign = 1;
	}
	if(choiceOfT3 == 1){
	  hit1->sign = 1;
	}
	if(choiceOfT3 == 2){
	  hit1->sign = -1;
	}
	if(choiceOfT3 == 3){
	  hit1->sign = -1;
	}
      }
    }
  
  return true;
  
}


bool KalmanFastTracking_NEW_2::compareTrackletsSlimV_3hits(Tracklet& tracklet2, Tracklet& tracklet3, int pass, double slopeComparison, double windowSize)
{
  
  if( ! ( (tracklet2.hits.size() == 2 && tracklet3.hits.size() == 1 ) || (tracklet2.hits.size() == 1 && tracklet3.hits.size() == 2 ) ) ) return false;
  
  double bestMatch = 100.;
  if(tracklet2.hits.size() == 2){
    
    unsigned int bestT2 = 5;
    for(unsigned int t2 = 0; t2 < tracklet2.possibleVLines.size(); t2++){
      double extrapolation = tracklet2.possibleVLines.at(t2).slopeV*(z_plane[tracklet3.getHit(0).hit.detectorID] - tracklet2.possibleVLines.at(t2).initialZ) + tracklet2.possibleVLines.at(t2).initialV;
      if(std::abs(tracklet3.getHit(0).hit.pos - extrapolation) < bestMatch && std::abs(tracklet3.getHit(0).hit.pos - extrapolation) < windowSize){
	bestMatch = std::abs(tracklet3.getHit(0).hit.pos - extrapolation);
	bestT2 = t2;
      }
    }
    
    if(bestMatch > 99. || bestT2 == 5){
      return false;
    } else{
      
      double newSlopeV = (tracklet3.getHit(0).hit.pos - tracklet2.possibleVLines.at(bestT2).initialV)/(z_plane[tracklet3.getHit(0).hit.detectorID] - tracklet2.possibleVLines.at(bestT2).initialZ);
      
      tracklet2.acceptedVLine2 = tracklet2.possibleVLines.at(bestT2);
      tracklet3.acceptedVLine3 = tracklet2.possibleVLines.at(bestT2);
      tracklet2.tx = newSlopeV;
      tracklet3.tx = newSlopeV;
      
      tracklet2.st2Vsl = newSlopeV;
      tracklet3.st2Vsl = newSlopeV;


      for(std::list<SignedHit>::iterator hit1 = tracklet2.hits.begin(); hit1 != tracklet2.hits.end(); ++hit1)
	{
	  if(hit1->hit.detectorID == 13){
	    if(bestT2 == 0){
	      hit1->sign = 1;
	    }
	    if(bestT2 == 1){
	      hit1->sign = -1;
	    }
	    if(bestT2 == 2){
	      hit1->sign = 1;
	    }
	    if(bestT2 == 3){
	      hit1->sign = -1;
	    }
	  } else{
	    if(bestT2 == 0){
	      hit1->sign = 1;
	    }
	    if(bestT2 == 1){
	      hit1->sign = 1;
	    }
	    if(bestT2 == 2){
	      hit1->sign = -1;
	    }
	    if(bestT2 == 3){
	      hit1->sign = -1;
	    }
	  }
	}

      
      return true;
    }
  }
  
  if(tracklet3.hits.size() == 2){
    
    unsigned int bestT3 = 5;
    for(unsigned int t3 = 0; t3 < tracklet3.possibleVLines.size(); t3++){
      double extrapolation = tracklet3.possibleVLines.at(t3).slopeV*(z_plane[tracklet2.getHit(0).hit.detectorID] - tracklet3.possibleVLines.at(t3).initialZ) + tracklet3.possibleVLines.at(t3).initialV;
      if(std::abs(tracklet2.getHit(0).hit.pos - extrapolation) < bestMatch && std::abs(tracklet2.getHit(0).hit.pos - extrapolation) < windowSize){
	bestMatch = std::abs(tracklet2.getHit(0).hit.pos - extrapolation);
	bestT3 = t3;
      }
    }
    
    if(bestMatch > 99. || bestT3 == 5){
      return false;
    } else{
      
      double newSlopeV = (tracklet2.getHit(0).hit.pos - tracklet3.possibleVLines.at(bestT3).initialV)/(z_plane[tracklet2.getHit(0).hit.detectorID] - tracklet3.possibleVLines.at(bestT3).initialZ);
      
      tracklet2.acceptedVLine2 = tracklet3.possibleVLines.at(bestT3);
      tracklet3.acceptedVLine3 = tracklet3.possibleVLines.at(bestT3);
      tracklet2.tx = newSlopeV;
      tracklet3.tx = newSlopeV;
      
      tracklet2.st2Vsl = newSlopeV;
      tracklet3.st2Vsl = newSlopeV;

      for(std::list<SignedHit>::iterator hit1 = tracklet3.hits.begin(); hit1 != tracklet3.hits.end(); ++hit1)
	{
	  if(hit1->hit.detectorID == 23 || hit1->hit.detectorID == 29){
	    if(bestT3 == 0){
	      hit1->sign = 1;
	    }
	    if(bestT3 == 1){
	  hit1->sign = -1;
	    }
	    if(bestT3 == 2){
	      hit1->sign = 1;
	    }
	    if(bestT3 == 3){
	      hit1->sign = -1;
	    }
	  } else{
	    if(bestT3 == 0){
	      hit1->sign = 1;
	    }
	    if(bestT3 == 1){
	      hit1->sign = 1;
	    }
	    if(bestT3 == 2){
	      hit1->sign = -1;
	    }
	    if(bestT3 == 3){
	      hit1->sign = -1;
	    }
	  }
	}

      
      return true;
      
    }
  }
}

