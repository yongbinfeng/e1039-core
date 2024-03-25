/*
KalmanFastTracking_Displaced.cxx

Implementation of class Tracklet, KalmanFastTracking_Displaced

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

#include "KalmanFastTracking_Displaced.h"
#include "TriggerRoad.h"

// #define _DEBUG_ON

namespace
{
    // static flag to indicate the initialized has been done
    static bool inited = false;

    // Event acceptance cut
    static int MaxHitsDC0;
    static int MaxHitsDC1;
    static int MaxHitsDC2;
    static int MaxHitsDC3p;
    static int MaxHitsDC3m;

    // Sagitta ratio
    static double SAGITTA_DUMP_CENTER;
    static double SAGITTA_DUMP_WIDTH;
    static double SAGITTA_TARGET_CENTER;
    static double SAGITTA_TARGET_WIDTH;
    static double Z_TARGET;
    static double Z_DUMP;

    // Track quality cuts
    static double TX_MAX;
    static double TY_MAX;
    static double X0_MAX;
    static double Y0_MAX;
    static double INVP_MAX;
    static double INVP_MIN;
    static double Z_KMAG_BEND;

    // MuID cuts
    static double MUID_REJECTION;
    static double MUID_Z_REF;
    static double MUID_R_CUT;
    static double MUID_THE_P0;
    static double MUID_EMP_P0;
    static double MUID_EMP_P1;
    static double MUID_EMP_P2;
    static int MUID_MINHITS;

    // Track merging threshold
    static double MERGE_THRES;

    // static flag of kmag on/off
    static bool KMAG_ON;

    // running mode
    static bool MC_MODE;
    static bool COSMIC_MODE;
    static bool COARSE_MODE;

    // if displaced, skip fit to the target/vertex
    static bool TRACK_ELECTRONS; // please see comment in framework/phool/recoConsts.cc
    static bool TRACK_DISPLACED; // please see comment in framework/phool/recoConsts.cc

    static double KMAGSTR;
    static double PT_KICK_KMAG;
    // initialize global variables
    void initGlobalVariables()
    {
        if (!inited)
        {
            inited = true;

            recoConsts *rc = recoConsts::instance();
            MC_MODE = rc->get_BoolFlag("MC_MODE");
            KMAG_ON = rc->get_BoolFlag("KMAG_ON");
            COSMIC_MODE = rc->get_BoolFlag("COSMIC_MODE");
            COARSE_MODE = rc->get_BoolFlag("COARSE_MODE");

            TRACK_ELECTRONS = rc->get_BoolFlag("TRACK_ELECTRONS");
            TRACK_DISPLACED = rc->get_BoolFlag("TRACK_DISPLACED");

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
            PT_KICK_KMAG = rc->get_DoubleFlag("PT_KICK_KMAG") * KMAGSTR;
        }
    }
}

KalmanFastTracking_Displaced::KalmanFastTracking_Displaced(const PHField *field, const TGeoManager *geom, bool flag) : verbosity(0), enable_KF(flag), outputListIdx(4)
{
    using namespace std;
    initGlobalVariables();

#ifdef _DEBUG_ON
    cout << "Initialization of KalmanFastTracking_Displaced ..." << endl;
    cout << "========================================" << endl;
#endif

    _timers.insert(std::make_pair<std::string, PHTimer *>("st2", new PHTimer("st2")));
    _timers.insert(std::make_pair<std::string, PHTimer *>("st3", new PHTimer("st3")));
    _timers.insert(std::make_pair<std::string, PHTimer *>("connect", new PHTimer("connect")));
    _timers.insert(std::make_pair<std::string, PHTimer *>("st23", new PHTimer("st23")));
    _timers.insert(std::make_pair<std::string, PHTimer *>("global", new PHTimer("global")));
    _timers.insert(std::make_pair<std::string, PHTimer *>("global_st1", new PHTimer("global_st1")));
    _timers.insert(std::make_pair<std::string, PHTimer *>("global_link", new PHTimer("global_link")));
    _timers.insert(std::make_pair<std::string, PHTimer *>("global_kalman", new PHTimer("global_kalman")));
    _timers.insert(std::make_pair<std::string, PHTimer *>("kalman", new PHTimer("kalman")));

    // Initialize Kalman fitter
    if (enable_KF)
    {
        kmfitter = new KalmanFitter(field, geom);
        kmfitter->setControlParameter(50, 0.001);
    }

    // Initialize minuit minimizer
    minimizer[0] = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
    minimizer[1] = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
    fcn = ROOT::Math::Functor(&tracklet_curr, &Tracklet::Eval, KMAG_ON ? 5 : 4);
    for (int i = 0; i < 2; ++i)
    {
        minimizer[i]->SetMaxFunctionCalls(1000000);
        minimizer[i]->SetMaxIterations(100);
        minimizer[i]->SetTolerance(1E-2);
        minimizer[i]->SetFunction(fcn);
        minimizer[i]->SetPrintLevel(0);
    }

    // Minimize ROOT output
    extern Int_t gErrorIgnoreLevel;
    gErrorIgnoreLevel = 9999;

    // Initialize geometry service
    p_geomSvc = GeomSvc::instance();
#ifdef _DEBUG_ON
    p_geomSvc->printTable();
    p_geomSvc->printWirePosition();
    p_geomSvc->printAlignPar();
#endif

    // Initialize plane angles for all planes
    for (int i = 1; i <= nChamberPlanes; ++i)
    {
        costheta_plane[i] = p_geomSvc->getCostheta(i);
        sintheta_plane[i] = p_geomSvc->getSintheta(i);
    }

    // Initialize hodoscope IDs
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

    // Register masking stations for tracklets in station-0/1, 2, 3+/-
    stationIDs_mask[0].push_back(1);
    stationIDs_mask[1].push_back(1);
    stationIDs_mask[2].push_back(2);
    stationIDs_mask[3].push_back(3);
    stationIDs_mask[4].push_back(3);

    // Masking stations for back partial
    stationIDs_mask[5].push_back(2);
    stationIDs_mask[5].push_back(3);
    stationIDs_mask[5].push_back(4);

    // Masking stations for global track
    stationIDs_mask[6].push_back(1);
    stationIDs_mask[6].push_back(2);
    stationIDs_mask[6].push_back(3);
    stationIDs_mask[6].push_back(4);

    // prop. tube IDs for mu id
    detectorIDs_muid[0][0] = p_geomSvc->getDetectorID("P1X1");
    detectorIDs_muid[0][1] = p_geomSvc->getDetectorID("P1X2");
    detectorIDs_muid[0][2] = p_geomSvc->getDetectorID("P2X1");
    detectorIDs_muid[0][3] = p_geomSvc->getDetectorID("P2X2");
    detectorIDs_muid[1][0] = p_geomSvc->getDetectorID("P1Y1");
    detectorIDs_muid[1][1] = p_geomSvc->getDetectorID("P1Y2");
    detectorIDs_muid[1][2] = p_geomSvc->getDetectorID("P2Y1");
    detectorIDs_muid[1][3] = p_geomSvc->getDetectorID("P2Y2");

    // Reference z_ref for mu id
    z_ref_muid[0][0] = MUID_Z_REF;
    z_ref_muid[0][1] = MUID_Z_REF;
    z_ref_muid[0][2] = 0.5 * (p_geomSvc->getPlanePosition(detectorIDs_muid[0][0]) + p_geomSvc->getPlanePosition(detectorIDs_muid[0][1]));
    z_ref_muid[0][3] = z_ref_muid[0][2];

    z_ref_muid[1][0] = MUID_Z_REF;
    z_ref_muid[1][1] = MUID_Z_REF;
    z_ref_muid[1][2] = 0.5 * (p_geomSvc->getPlanePosition(detectorIDs_muid[1][0]) + p_geomSvc->getPlanePosition(detectorIDs_muid[1][1]));
    z_ref_muid[1][3] = z_ref_muid[1][2];

    // Initialize masking window sizes, with optimized contingency
    for (int i = nChamberPlanes + 1; i <= nChamberPlanes + nHodoPlanes + nPropPlanes; i++)
    {
        double factor = 0.;
        if (i > nChamberPlanes && i <= nChamberPlanes + 4)
            factor = 0.25; // for station-1 hodo
        if (i > nChamberPlanes + 4 && i <= nChamberPlanes + 8)
            factor = 0.2; // for station-2 hodo
        if (i > nChamberPlanes + 8 && i <= nChamberPlanes + 10)
            factor = 0.15; // for station-3 hodo
        if (i > nChamberPlanes + 10 && i <= nChamberPlanes + nHodoPlanes)
            factor = 0.; // for station-4 hodo
        if (i > nChamberPlanes + nHodoPlanes)
            factor = 0.15; // for station-4 proptube

        z_mask[i - nChamberPlanes - 1] = p_geomSvc->getPlanePosition(i);
        for (int j = 1; j <= p_geomSvc->getPlaneNElements(i); j++)
        {
            double x_min, x_max, y_min, y_max;
            p_geomSvc->get2DBoxSize(i, j, x_min, x_max, y_min, y_max);

            x_min -= (factor * (x_max - x_min));
            x_max += (factor * (x_max - x_min));
            y_min -= (factor * (y_max - y_min));
            y_max += (factor * (y_max - y_min));

            x_mask_min[i - nChamberPlanes - 1][j - 1] = x_min;
            x_mask_max[i - nChamberPlanes - 1][j - 1] = x_max;
            y_mask_min[i - nChamberPlanes - 1][j - 1] = y_min;
            y_mask_max[i - nChamberPlanes - 1][j - 1] = y_max;
        }
    }

#ifdef _DEBUG_ON
    cout << "========================" << endl;
    cout << "Hodo. masking settings: " << endl;
    for (int i = 0; i < 4; i++)
    {
        cout << "For station " << i + 1 << endl;
        for (std::vector<int>::iterator iter = detectorIDs_mask[i].begin(); iter != detectorIDs_mask[i].end(); ++iter)
            cout << "All: " << *iter << endl;
        for (std::vector<int>::iterator iter = detectorIDs_maskX[i].begin(); iter != detectorIDs_maskX[i].end(); ++iter)
            cout << "X:   " << *iter << endl;
        for (std::vector<int>::iterator iter = detectorIDs_maskY[i].begin(); iter != detectorIDs_maskY[i].end(); ++iter)
            cout << "Y:   " << *iter << endl;
    }

    for (int i = 0; i < nStations; ++i)
    {
        std::cout << "Masking stations for tracklets with stationID = " << i + 1 << ": " << std::endl;
        for (std::vector<int>::iterator iter = stationIDs_mask[i].begin(); iter != stationIDs_mask[i].end(); ++iter)
        {
            std::cout << *iter << "  ";
        }
        std::cout << std::endl;
    }
#endif

    // Initialize super stationIDs
    for (int i = 0; i < nChamberPlanes / 6 + 2; i++)
        superIDs[i].clear();
    superIDs[0].push_back((p_geomSvc->getDetectorIDs("D0X")[0] + 1) / 2);
    superIDs[0].push_back((p_geomSvc->getDetectorIDs("D0U")[0] + 1) / 2);
    superIDs[0].push_back((p_geomSvc->getDetectorIDs("D0V")[0] + 1) / 2);
    superIDs[1].push_back((p_geomSvc->getDetectorIDs("D1X")[0] + 1) / 2);
    superIDs[1].push_back((p_geomSvc->getDetectorIDs("D1U")[0] + 1) / 2);
    superIDs[1].push_back((p_geomSvc->getDetectorIDs("D1V")[0] + 1) / 2);
    superIDs[2].push_back((p_geomSvc->getDetectorIDs("D2X")[0] + 1) / 2);
    superIDs[2].push_back((p_geomSvc->getDetectorIDs("D2U")[0] + 1) / 2);
    superIDs[2].push_back((p_geomSvc->getDetectorIDs("D2V")[0] + 1) / 2);
    superIDs[3].push_back((p_geomSvc->getDetectorIDs("D3pX")[0] + 1) / 2);
    superIDs[3].push_back((p_geomSvc->getDetectorIDs("D3pU")[0] + 1) / 2);
    superIDs[3].push_back((p_geomSvc->getDetectorIDs("D3pV")[0] + 1) / 2);
    superIDs[4].push_back((p_geomSvc->getDetectorIDs("D3mX")[0] + 1) / 2);
    superIDs[4].push_back((p_geomSvc->getDetectorIDs("D3mU")[0] + 1) / 2);
    superIDs[4].push_back((p_geomSvc->getDetectorIDs("D3mV")[0] + 1) / 2);

    superIDs[5].push_back((p_geomSvc->getDetectorIDs("P1X")[0] + 1) / 2);
    superIDs[5].push_back((p_geomSvc->getDetectorIDs("P2X")[0] + 1) / 2);
    superIDs[6].push_back((p_geomSvc->getDetectorIDs("P1Y")[0] + 1) / 2);
    superIDs[6].push_back((p_geomSvc->getDetectorIDs("P2Y")[0] + 1) / 2);

#ifdef _DEBUG_ON
    cout << "=============" << endl;
    cout << "Chamber IDs: " << endl;
    TString stereoNames[3] = {"X", "U", "V"};
    for (int i = 0; i < nChamberPlanes / 6; i++)
    {
        for (int j = 0; j < 3; j++)
            cout << i << "  " << stereoNames[j].Data() << ": " << superIDs[i][j] << endl;
    }

    cout << "Proptube IDs: " << endl;
    for (int i = nChamberPlanes / 6; i < nChamberPlanes / 6 + 2; i++)
    {
        for (int j = 0; j < 2; j++)
            cout << i << "  " << j << ": " << superIDs[i][j] << endl;
    }

    // Initialize widow sizes for X-U matching and z positions of all chambers
    cout << "======================" << endl;
    cout << "U plane window sizes: " << endl;
#endif

    double u_factor[] = {5., 5., 5., 15., 15.};
    for (int i = 0; i < nChamberPlanes / 6; i++)
    {
        int xID = 2 * superIDs[i][0] - 1;
        int uID = 2 * superIDs[i][1] - 1;
        int vID = 2 * superIDs[i][2] - 1;
        double spacing = p_geomSvc->getPlaneSpacing(uID);
        double x_span = p_geomSvc->getPlaneScaleY(uID);

        z_plane_x[i] = 0.5 * (p_geomSvc->getPlanePosition(xID) + p_geomSvc->getPlanePosition(xID + 1));
        z_plane_u[i] = 0.5 * (p_geomSvc->getPlanePosition(uID) + p_geomSvc->getPlanePosition(uID + 1));
        z_plane_v[i] = 0.5 * (p_geomSvc->getPlanePosition(vID) + p_geomSvc->getPlanePosition(vID + 1));

        u_costheta[i] = costheta_plane[uID];
        u_sintheta[i] = sintheta_plane[uID];

        // u_win[i] = fabs(0.5*x_span/(spacing/sintheta_plane[uID])) + 2.*spacing + u_factor[i];
        u_win[i] = fabs(0.5 * x_span * sintheta_plane[uID]) + TX_MAX * fabs((z_plane_u[i] - z_plane_x[i]) * u_costheta[i]) + TY_MAX * fabs((z_plane_u[i] - z_plane_x[i]) * u_sintheta[i]) + 2. * spacing + u_factor[i];

#ifdef _DEBUG_ON
        cout << "Station " << i << ": " << xID << "  " << uID << "  " << vID << "  " << u_win[i] << endl;
#endif
    }

    // Initialize Z positions and maximum parameters of all planes
    for (int i = 1; i <= nChamberPlanes; i++)
    {
        z_plane[i] = p_geomSvc->getPlanePosition(i);
        slope_max[i] = costheta_plane[i] * TX_MAX + sintheta_plane[i] * TY_MAX;
        intersection_max[i] = costheta_plane[i] * X0_MAX + sintheta_plane[i] * Y0_MAX;

#ifdef COARSE_MODE
        resol_plane[i] = 3. * p_geomSvc->getPlaneSpacing(i) / sqrt(12.);
#else
        if (i <= 6)
        {
            resol_plane[i] = recoConsts::instance()->get_DoubleFlag("RejectWinDC0");
        }
        else if (i <= 12)
        {
            resol_plane[i] = recoConsts::instance()->get_DoubleFlag("RejectWinDC1");
        }
        else if (i <= 18)
        {
            resol_plane[i] = recoConsts::instance()->get_DoubleFlag("RejectWinDC2");
        }
        else if (i <= 24)
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
    for (int i = 1; i <= nChamberPlanes / 2; i++)
    {
        double d_slope = (p_geomSvc->getPlaneResolution(2 * i - 1) + p_geomSvc->getPlaneResolution(2 * i)) / (z_plane[2 * i] - z_plane[2 * i - 1]);
        double d_intersection = d_slope * z_plane[2 * i];

        slope_max[2 * i - 1] += d_slope;
        intersection_max[2 * i - 1] += d_intersection;
        slope_max[2 * i] += d_slope;
        intersection_max[2 * i] += d_intersection;

#ifdef _DEBUG_ON
        cout << "Super plane " << i << ": " << slope_max[2 * i - 1] << "  " << intersection_max[2 * i - 1] << endl;
#endif
    }

    // Initialize sagitta ratios, index 0, 1, 2 are for X, U, V, this is the incrementing order of plane type
    s_detectorID[0] = p_geomSvc->getDetectorID("D2X");
    s_detectorID[1] = p_geomSvc->getDetectorID("D2Up");
    s_detectorID[2] = p_geomSvc->getDetectorID("D2Vp");
}

KalmanFastTracking_Displaced::~KalmanFastTracking_Displaced()
{
    if (enable_KF)
        delete kmfitter;
    delete minimizer[0];
    delete minimizer[1];
}

void KalmanFastTracking_Displaced::setRawEventDebug(SRawEvent *event_input)
{
    rawEvent = event_input;
    hitAll = event_input->getAllHits();
    std::cout << "test mark: "
              << "Hello world" << std::endl;
}

int KalmanFastTracking_Displaced::setRawEventPrep(SRawEvent *event_input)
{
    // reset timer
    for (auto iter = _timers.begin(); iter != _timers.end(); ++iter)
    {
        iter->second->reset();
    }

    // Initialize tracklet lists (below tracklets in individual stations)
    for (int i = 0; i < 5; i++)
    {
        trackletsInSt[i].clear();
        for (int j = 0; j < 3; j++)
            trackletsInStSlim[i][j].clear();
    }
    stracks.clear();

    // Initialize station 2+3 tracklet lists
    trackletsInSt23Slim[0].clear();
    trackletsInSt23Slim[1].clear();
    trackletsInSt23Slim[2].clear();

    // Initialize global tracklet lists
    globalTracklets.clear();
    globalTracklets_resolveSt1.clear();

    // pre-tracking cuts
    rawEvent = event_input;
    if (!acceptEvent(rawEvent))
        return TFEXIT_FAIL_MULTIPLICITY;

    hitAll = event_input->getAllHits();

#ifdef _DEBUG_ON
    for (std::vector<Hit>::iterator iter = hitAll.begin(); iter != hitAll.end(); ++iter)
        iter->print();
#endif

    // Initialize hodo and masking IDs
    for (int i = 0; i < 4; i++)
    {
        hitIDs_mask[i].clear();
        hitIDs_mask[i] = rawEvent->getHitsIndexInDetectorsNoRepeats(detectorIDs_mask[i]);
        hitIDs_maskX[i].clear();
        hitIDs_maskX[i] = rawEvent->getHitsIndexInDetectorsNoRepeats(detectorIDs_maskX[i]);
        hitIDs_maskY[i].clear();
        hitIDs_maskY[i] = rawEvent->getHitsIndexInDetectorsNoRepeats(detectorIDs_maskY[i]);
    }
    // Initialize prop. tube IDs
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            hitIDs_muid[i][j].clear();
            hitIDs_muid[i][j] = rawEvent->getHitsIndexInDetector(detectorIDs_muid[i][j]);
        }
        hitIDs_muidHodoAid[i].clear();
        hitIDs_muidHodoAid[i] = rawEvent->getHitsIndexInDetectors(detectorIDs_muidHodoAid[i]);
    }

    if (!COARSE_MODE)
    {
        buildPropSegments();
        if (propSegs[0].empty() || propSegs[1].empty())
        {
#ifdef _DEBUG_ON
            LogInfo("Failed in prop tube segment building: " << propSegs[0].size() << ", " << propSegs[1].size());
#endif
            // return TFEXIT_FAIL_ROUGH_MUONID;
        }
    }

    return 0;
}

int KalmanFastTracking_Displaced::setRawEvent(SRawEvent *event_input)
{

    int ret = setRawEventPrep(event_input);
    if (ret != 0)
        return ret;

    totalTime = 0;

    // Get hit combinations in station 2
    _timers["st2"]->restart();
    for (int plane_type = 0; plane_type < 3; plane_type++)
        buildTrackletsInStationSlim(3, 1, plane_type); 
    if (verbosity >= 2)
    {
        std::cout << "Station 2 combos- x : u : v";
        for (int i = 0; i < 3; i++)
            std::cout << trackletsInStSlim[1][i].size() << ": ";
        std::cout << std::endl;
    }
    _timers["st2"]->stop();
    totalTime += _timers["st2"]->get_accumulated_time() / 1000.;

    // Get hit combinations in station 3+ and 3-
    _timers["st3"]->restart();
    for (int plane_type = 0; plane_type < 3; plane_type++)
        for (int station_id = 4; station_id < 6; station_id++)
            buildTrackletsInStationSlim(station_id, 2, plane_type);
    if (verbosity >= 2)
    {
        std::cout << "Station 3 combos- x : u : v";
        for (int i = 0; i < 3; i++)
            std::cout << trackletsInStSlim[2][i].size() << ": ";
        std::cout << std::endl;
    }
    _timers["st3"]->stop();
    totalTime += _timers["st3"]->get_accumulated_time() / 1000.;

    // Matching of station 2 and station 3 hits separately for the three wire tilts
    _timers["connect"]->restart();
    if (buildBackPartialTracksSlim())
        return TFEXIT_FAIL_BACKPARTIAL;
    _timers["connect"]->stop();
    totalTime += _timers["connect"]->get_accumulated_time() / 1000.;

    // Now make the actual station 2+3 tracklets by combining tracklets of different plane types
    _timers["st23"]->restart();
    buildFullBackPartialTracksSlim(); // This should be relatively fast given the binned combination method
    if (verbosity >= 2)
        LogInfo("NTracklets St2+St3: " << trackletsInSt[3].size());
    if (outputListIdx > 3 && trackletsInSt[3].empty())
    {
        return TFEXIT_FAIL_BACKPARTIAL;
    }
    _timers["st23"]->stop();
    totalTime += _timers["st23"]->get_accumulated_time() / 1000.;

    // Connect tracklets in station 2+3 and station 1 to form global tracks
    _timers["global"]->restart();
    if (!TRACK_DISPLACED)
        buildGlobalTracks();
    else
        buildGlobalTracksDisplaced();
    _timers["global"]->stop();
    totalTime += _timers["global"]->get_accumulated_time() / 1000.;
    if (outputListIdx == 4 && trackletsInSt[4].empty())
        return TFEXIT_FAIL_GLOABL;
    if (!enable_KF)
        return TFEXIT_SUCCESS;

    // Build kalman tracks
    _timers["kalman"]->restart();
    for (std::list<Tracklet>::iterator tracklet = trackletsInSt[outputListIdx].begin(); tracklet != trackletsInSt[outputListIdx].end(); ++tracklet)
    {
        SRecTrack strack = processOneTracklet(*tracklet);
        stracks.push_back(strack);
    }
    _timers["kalman"]->stop();
    totalTime += _timers["kalman"]->get_accumulated_time() / 1000.;

    return TFEXIT_SUCCESS;
}

bool KalmanFastTracking_Displaced::acceptEvent(SRawEvent *rawEvent)
{
    if (Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
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
        std::cout << MaxHitsDC0 << " " << MaxHitsDC1 << " " << MaxHitsDC2 << " " << MaxHitsDC3p << " " << MaxHitsDC3m << std::endl;
    }

    if (rawEvent->getNHitsInD0() > MaxHitsDC0)
        return false;
    // if(rawEvent->getNHitsInD1() > MaxHitsDC1) return false; ??????????
    if (rawEvent->getNHitsInD2() > MaxHitsDC2)
        return false;
    if (rawEvent->getNHitsInD3p() > MaxHitsDC3p)
        return false;
    if (rawEvent->getNHitsInD3m() > MaxHitsDC3m)
        return false;

    if (rawEvent->getNHitsInDetectors(detectorIDs_maskX[1]) < 1 || rawEvent->getNUniqueHitsInDetectors(detectorIDs_maskX[1]) > 50)
        return false;
    if (rawEvent->getNHitsInDetectors(detectorIDs_maskX[2]) < 1 || rawEvent->getNUniqueHitsInDetectors(detectorIDs_maskX[2]) > 50)
        return false;

    adjusted = false;

    // high pileup, tight cuts on tracklet matching
    slopeComp = .15;
    windowSize = 15;
    reqHits = 2;
    slopeCompSt1 = .18;
    XWinSt1 = 1.25;
    UVWinSt1 = 1.5;
    hodoXWin = 7;
    hodoUVWin = 30;
    hodoXDIFFWin = 7;
    hodoUVDIFFWin = 7;
    XUVSlopeWin = 0.04;
    XUVPosWin = 9.;
    YSlopesDiff = 0.007;
    XSlopesDiff = 0.007;
    chiSqCut = 100;
    st23ChiSqCut = 15;
    highPU = true;

    if (rawEvent->getNHitsInD0() < 1000 && rawEvent->getNHitsInD2() < 500 && rawEvent->getNHitsInD3p() < 500 && rawEvent->getNHitsInD3m() < 500)
    {
        slopeComp = .25;
        windowSize = 27.5;
        reqHits = 2;
        slopeCompSt1 = .25;
        XWinSt1 = 1.25;
        UVWinSt1 = 1.5;
        hodoXWin = 14;
        hodoUVWin = 40;
        hodoXDIFFWin = 10;
        hodoUVDIFFWin = 13;
        XUVSlopeWin = 0.04;
        XUVPosWin = 9.;
        YSlopesDiff = 0.007;
        XSlopesDiff = 0.007;
        chiSqCut = m_chiSqCut;
        st23ChiSqCut = 30;
        highPU = false;

        if (rawEvent->getNHitsInD2() < 300 && rawEvent->getNHitsInD3p() < 300 && rawEvent->getNHitsInD3m() < 300)
        {
            reqHits = 1;

            if (rawEvent->getNHitsInD0() < 500)
            {
                slopeCompSt1 = m_slopeComparisonSt1;
                UVWinSt1 = 2.25;
                hodoUVDIFFWin = 20;
                st23ChiSqCut = 50;

                if (rawEvent->getNHitsInD0() < 400)
                {
                    slopeComp = m_slopeComparison;
                    windowSize = m_windowSize;
                    slopeCompSt1 = m_slopeComparisonSt1;
                    XWinSt1 = m_XWindowSt1;
                    UVWinSt1 = m_UVWindowSt1;
                    hodoXWin = m_hodoXWindow;
                    hodoUVWin = m_hodoUVWindow;
                    hodoXDIFFWin = m_hodoXDIFFWindow;
                    hodoUVDIFFWin = m_hodoUVDIFFWindow;
                    XUVSlopeWin = m_XUVSlopeWindowCoarse;
                    XUVPosWin = m_XUVPosWindowCoarse;
                    YSlopesDiff = m_YSlopesDiff;
                    XSlopesDiff = m_XSlopesDiff;
                    chiSqCut = m_chiSqCut;
                    st23ChiSqCut = m_st23ChiSqCut;
                }
            }
        }
    }
    return true;
}

// This function finds valid hits combinations for a super detector in a station
void KalmanFastTracking_Displaced::buildTrackletsInStationSlim(int stationID, int listID, int plane_type, double *pos_exp, double *window)
{ /// All this function do is to contruct tracklets only contain a hit pair in a certain superdetector with a set possible track lines based on the 4 possible sign combination
    // actuall ID of the tracklet lists
    int sID = stationID - 1;
    
    // Extract the X, U, V hit pairs
    std::list<SRawEvent::hit_pair> pairs;
    //Notice there is a tricky aspect of the detector: the direction of U(V) plane in station 4 and station 5 are different 
    int extract_number = plane_type;
    if (stationID == 4 || stationID == 5)
        extract_number = (3 - plane_type) % 3; 

    if (pos_exp == nullptr)
        pairs = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][extract_number]);
    else
        pairs = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][extract_number], pos_exp[extract_number], window[extract_number]);

    if (pairs.empty())
        return;

    for (auto pair : pairs)
    {

        int LR1 = 0;
        int LR2 = 0;
        Tracklet tracklet_new;
        tracklet_new.stationID = stationID;

        if (pair.first >= 0)
        {
            tracklet_new.hits.push_back(SignedHit(hitAll[pair.first], LR1));
            tracklet_new.nHits[plane_type]++;
        }
        if (pair.second >= 0)
        {
            tracklet_new.hits.push_back(SignedHit(hitAll[pair.second], LR2));
            tracklet_new.nHits[plane_type]++;
        }

        int hitindex = pair.first >= 0 ? pair.first : pair.second;
        if (hitindex >= 0)
        {
            if (stationID == 3)
            {
                tracklet_new.pos_st2 = hitAll[hitindex].pos;

                tracklet_new.z_st2 = z_plane[hitAll[hitindex].detectorID];
            }
            if (stationID == 4 || stationID == 5)
            {
                tracklet_new.pos_st3 = hitAll[hitindex].pos;
                tracklet_new.z_st3 = z_plane[hitAll[hitindex].detectorID];
                if (hitAll[hitindex].detectorID > 18 && hitAll[hitindex].detectorID < 25)
                {
                    tracklet_new.isM = false;
                }
                else
                {
                    tracklet_new.isM = true;
                }
            }
        }
        if (pair.first >= 0 && pair.second >= 0)
        {
            tracklet_new.getSlopes(hitAll[pair.second], hitAll[pair.first], plane_type); 
            tracklet_new.sortHits();
        }
        trackletsInStSlim[listID][plane_type].push_back(tracklet_new);
    }
}

//This function builds back partial tracklets for all three plane types
int KalmanFastTracking_Displaced::buildBackPartialTracksSlim()
{
	///XL: For low pile up situation, We can simply use the commented code. I will investigate the needed cuts for real data in the future
    /*
    for (int i=0;i<3;i++)buildBackPartialTracksSlim(i, reqHits, slopeComp, windowSize);
    if(verbosity >= 2) {
    std::cout<<"Station 2+3 combos: x: y: z"<<std::endl;
    for (int i=0;i<3;i++)std::cout<<trackletsInSt23Slim[i].size()<<": ";
    std::cout<<std::endl;
    }
    buildXUCombos();
    buildXUVCombos();
    return 0;
    */

    buildBackPartialTracksSlim(0, reqHits, slopeComp, windowSize);
    if (verbosity >= 2)
        std::cout << "Station 2+3 x combos = " << trackletsInSt23Slim[0].size() << std::endl;
    if (trackletsInSt23Slim[0].size() == 0)
        return TFEXIT_FAIL_BACKPARTIAL;

    buildBackPartialTracksSlim(1, reqHits, slopeComp, windowSize);
    buildXUCombos();
    if (verbosity >= 2)
        std::cout << "Station 2+3 u combos = " << trackletsInSt23Slim[1].size() << std::endl;
    if (trackletsInSt23Slim[1].size() == 0)
        return TFEXIT_FAIL_BACKPARTIAL;

    int numHodoValidUXCombos = 0;

    for (unsigned int tx = 0; tx < trackletsInSt23Slim[0].size(); tx++)
    {
        numHodoValidUXCombos += trackletsInSt23Slim[0].at(tx).allowedUXCombos.size();
    }
    if (numHodoValidUXCombos == 0)
        return TFEXIT_FAIL_BACKPARTIAL;

    if (numHodoValidUXCombos > 10000 && !highPU)
    {
        cutAdjuster(numHodoValidUXCombos, 1); 
        adjusted = true;
        for (int i = 0; i < 2; i++)
        {
            trackletsInSt23Slim[i].clear();
            buildBackPartialTracksSlim(i, reqHits, slopeComp, windowSize);
        }
        buildXUCombos();

        numHodoValidUXCombos = 0;
        for (unsigned int tx = 0; tx < trackletsInSt23Slim[0].size(); tx++)
        {
            numHodoValidUXCombos += trackletsInSt23Slim[0].at(tx).allowedUXCombos.size();
        }
    }
    else if (numHodoValidUXCombos > 50000 && highPU)
    {
        return TFEXIT_FAIL_BACKPARTIAL;
    }

    // Checked X+U combos of stations 2+3.  Now make V combos of stations 2+3
    buildBackPartialTracksSlim(2, reqHits, slopeComp, windowSize);
    buildXUVCombos();
    if (verbosity >= 2)
    {
        std::cout << "Station 2+3 v combos = " << trackletsInSt23Slim[2].size() << std::endl;
    }
    if (trackletsInSt23Slim[2].size() == 0)
        return TFEXIT_FAIL_BACKPARTIAL;

    int numHodoValidCombos = 0;
    for (unsigned int tx = 0; tx < trackletsInSt23Slim[0].size(); tx++)
    {
        numHodoValidCombos += trackletsInSt23Slim[0].at(tx).matching_combos.size();
    }
    if (numHodoValidCombos == 0)
        return TFEXIT_FAIL_BACKPARTIAL;

    if (numHodoValidCombos > 50000 && !highPU)
    {
        cutAdjuster(numHodoValidCombos, 2);
        for (int i = 0; i < 3; i++)
        {
            trackletsInSt23Slim[i].clear();
            buildBackPartialTracksSlim(i, reqHits, slopeComp, windowSize);
        }
        buildXUCombos();
        buildXUVCombos();
        numHodoValidCombos = 0;
        for (unsigned int tx = 0; tx < trackletsInSt23Slim[0].size(); tx++)
        {
            numHodoValidCombos += trackletsInSt23Slim[0].at(tx).matching_combos.size();
        }
        if (numHodoValidCombos >= 50000000)
            return TFEXIT_FAIL_BACKPARTIAL;
    }
    else if (highPU && numHodoValidCombos > 500000)
    {
        return TFEXIT_FAIL_BACKPARTIAL;
    }
    return 0;
}

//This function builds back partial tracklets for a perticular plane type
void KalmanFastTracking_Displaced::buildBackPartialTracksSlim(int plane_type, int pass, double slopeComparison, double windowSize)
{
    for (std::list<Tracklet>::iterator tracklet3 = trackletsInStSlim[2][plane_type].begin(); tracklet3 != trackletsInStSlim[2][plane_type].end(); ++tracklet3)
    {
        for (std::list<Tracklet>::iterator tracklet2 = trackletsInStSlim[1][plane_type].begin(); tracklet2 != trackletsInStSlim[1][plane_type].end(); ++tracklet2)
        {
            if ((compareTrackletsSlim(&*tracklet2, &*tracklet3, plane_type, pass, slopeComparison, windowSize) || (pass == 1 && compareTrackletsSlim_3hits(&*tracklet2, &*tracklet3, plane_type, pass, slopeComparison, windowSize))))
            {
                Tracklet tracklet_23 = (*tracklet2) + (*tracklet3);
                tracklet_23.pos_st2 = tracklet2->pos_st2;
                tracklet_23.z_st2 = tracklet2->z_st2;
                tracklet_23.pos_st3 = tracklet3->pos_st3;
                tracklet_23.z_st3 = tracklet3->z_st2;
                tracklet_23.slope_st2 = tracklet2->slope_st2;
                tracklet_23.slope_st3 = tracklet3->slope_st3;
                tracklet_23.isM = tracklet3->isM;
                tracklet_23.acceptedLine2 = tracklet2->acceptedLine2;
                tracklet_23.acceptedLine3 = tracklet3->acceptedLine3;

                if (plane_type == 0)
                {
                    tracklet_23.x0 = tracklet2->slope_st2 * (0. - tracklet2->z_st2) + tracklet2->pos_st2; // simple extrapolation back to z = 0
                    tracklet_23.ty = 0.;
                    tracklet_23.y0 = 0.;
                }
                double pos_0 = tracklet2->slope_st2 * (0. - tracklet2->z_st2) + tracklet2->pos_st2; // simple exrapolation back to z = 0

                if (std::abs(tracklet_23.slope_st2) > 0.15 || std::abs(pos_0) > 150)
                    continue;

                for (auto hit_sign : tracklet_23.hits)
                {
                    if (hit_sign.hit.index >= 0)
                    {
                        tracklet_23.planeNorm = {-costheta_plane[hit_sign.hit.detectorID], -sintheta_plane[hit_sign.hit.detectorID], tracklet_23.slope_st2};
                    }
                }

                std::vector<std::pair<int, int>> allowedHodos;
                std::vector<std::pair<int, double>> hodo2Diffs;
                std::vector<std::pair<int, double>> hodo3Diffs;
                bool passHodo2 = false;
                bool passHodo3 = false;

                int h2index = 0;
                for (const int stationID : stationIDs_mask[tracklet2->stationID - 1])
                {
                    for (const int iter : hitIDs_maskX[stationID - 1])
                    {
                        int idx1 = hitAll[iter].detectorID - nChamberPlanes - 1;
                        int idx2 = hitAll[iter].elementID - 1;
                        double z_hodo2 = z_mask[idx1];
                        double hodo2_pos = tracklet2->slope_st2 * (z_hodo2 - tracklet2->z_st2) + tracklet2->pos_st2;
                        double hodoWin = hodoUVWin;
                        if (plane_type == 0)
                            hodoWin = hodoXWin;
                        if (std::abs(hodo2_pos - hitAll[iter].pos) < hodoWin)
                        {
                            passHodo2 = true;
                            hodo2Diffs.push_back(std::make_pair(h2index, hodo2_pos - hitAll[iter].pos));
                        }

                        h2index++;
                    }
                }

                int h3index = 0;
                for (const int stationID : stationIDs_mask[tracklet3->stationID - 1])
                {
                    for (const int iter : hitIDs_maskX[stationID - 1])
                    {
                        int idx1 = hitAll[iter].detectorID - nChamberPlanes - 1;
                        int idx2 = hitAll[iter].elementID - 1;
                        double z_hodo3 = z_mask[idx1];
                        double hodo3_pos = tracklet2->slope_st2 * (z_hodo3 - tracklet2->z_st2) + tracklet2->pos_st2;
                        double hodoWin = hodoUVWin;
                        if (plane_type == 0)
                            hodoWin = hodoXWin;
                        if (std::abs(hodo3_pos - hitAll[iter].pos) < hodoWin)
                        {
                            passHodo3 = true;
                            hodo3Diffs.push_back(std::make_pair(h3index, hodo3_pos - hitAll[iter].pos));
                        }
                        h3index++;
                    }
                }

                if ((!passHodo2) || (!passHodo3))
                    continue;

                double hodoDIFFWin = hodoUVDIFFWin;
                if (plane_type == 0)
                    hodoDIFFWin = hodoXDIFFWin;
                for (int hd2 = 0; hd2 < hodo2Diffs.size(); hd2++)
                {
                    for (int hd3 = 0; hd3 < hodo3Diffs.size(); hd3++)
                    { /// Why do this strance minus operation? In my opinion there shoule be abs()+abs()
                        if (std::abs(hodo2Diffs.at(hd2).second - hodo3Diffs.at(hd3).second) < hodoDIFFWin)
                        { 
                            tracklet_23.allowedHodos.push_back(std::make_pair(hodo2Diffs.at(hd2).first, hodo3Diffs.at(hd3).first));
                        }
                    }
                }

                if (tracklet_23.allowedHodos.size() > 0)
                {
                    trackletsInSt23Slim[plane_type].push_back(tracklet_23);
                }
            }
        }
    }
}

//This function finds the matching U tracklets to a X tracklet
void KalmanFastTracking_Displaced::buildXUCombos()
{
    for (unsigned int tu23 = 0; tu23 < trackletsInSt23Slim[1].size(); tu23++)
    {
        for (unsigned int tx23 = 0; tx23 < trackletsInSt23Slim[0].size(); tx23++)
        {
            Tracklet &trackletu_23 = trackletsInSt23Slim[1].at(tu23);
            Tracklet &trackletx_23 = trackletsInSt23Slim[0].at(tx23);

            if (!(trackletx_23.isM == trackletu_23.isM))
                continue;
            TVector3 testCross = trackletu_23.planeNorm.Cross(trackletx_23.planeNorm);
            testCross *= 1. / testCross.Z();

            if (std::abs(testCross.X()) < .15 && std::abs(testCross.Y()) < .15)
            {
                Tracklet::UXCombo newCombo;
                newCombo.trackletUIndex = tu23;
                for (int aH = 0; aH < trackletu_23.allowedHodos.size(); aH++)
                {
                    for (unsigned int tx23AH = 0; tx23AH < trackletx_23.allowedHodos.size(); tx23AH++)
                    {
                        if (trackletu_23.allowedHodos.at(aH).first == trackletx_23.allowedHodos.at(tx23AH).first && trackletu_23.allowedHodos.at(aH).second == trackletx_23.allowedHodos.at(tx23AH).second)
                        {
                            newCombo.hodoMatches.push_back(std::make_pair(trackletu_23.allowedHodos.at(aH).first, trackletu_23.allowedHodos.at(aH).second));
                        }
                    }
                }

                if (newCombo.hodoMatches.size() > 0)
                {
                    newCombo.tx = testCross.X();
                    newCombo.ty = testCross.Y();
                    newCombo.y_st2 = costheta_plane[17] / sintheta_plane[17] * (trackletu_23.pos_st2 * costheta_plane[17] - trackletx_23.pos_st2) + trackletu_23.pos_st2 * sintheta_plane[17];
                    trackletx_23.allowedUXCombos.push_back(newCombo);
                }
            }
        }
    }
}

///XL: If we don't need to do cut examination after building U tracklets, these two functions can be merged
//This function finds the matching U and V tracklets to a X tracklet
void KalmanFastTracking_Displaced::buildXUVCombos()
{
    for (unsigned int tv23 = 0; tv23 < trackletsInSt23Slim[2].size(); tv23++)
    {
        for (unsigned int tx23 = 0; tx23 < trackletsInSt23Slim[0].size(); tx23++)
        {
            Tracklet &trackletx_23 = trackletsInSt23Slim[0].at(tx23);
            Tracklet &trackletv_23 = trackletsInSt23Slim[2].at(tv23);
            if (!(trackletx_23.isM == trackletv_23.isM))
                continue;
            TVector3 testCross = trackletv_23.planeNorm.Cross(trackletx_23.planeNorm);
            testCross *= 1. / testCross.Z();

            for (unsigned int tx23UC = 0; tx23UC < trackletx_23.allowedUXCombos.size(); tx23UC++)
            {
                bool match = 0;
                for (int aH = 0; aH < trackletv_23.allowedHodos.size(); aH++)
                {
                    for (unsigned int tx23HM = 0; tx23HM < trackletx_23.allowedUXCombos.at(tx23UC).hodoMatches.size(); tx23HM++)
                    {
                        if (trackletv_23.allowedHodos.at(aH).first == trackletx_23.allowedUXCombos.at(tx23UC).hodoMatches.at(tx23HM).first && trackletv_23.allowedHodos.at(aH).second == trackletx_23.allowedUXCombos.at(tx23UC).hodoMatches.at(tx23HM).second)
                        {
                            if (std::abs(testCross.X() - trackletx_23.allowedUXCombos.at(tx23UC).tx) < .01 && std::abs(testCross.Y() - trackletx_23.allowedUXCombos.at(tx23UC).ty) < 0.02)
                            {
                                if (std::abs(trackletx_23.allowedUXCombos.at(tx23UC).y_st2 - (costheta_plane[13] / sintheta_plane[13] * (trackletv_23.pos_st2 * costheta_plane[13] - trackletx_23.pos_st2) + trackletv_23.pos_st2 * sintheta_plane[13])) < 30.)
                                {
                                    match = 1;
                                    trackletx_23.matching_combos.push_back(std::make_pair(trackletx_23.allowedUXCombos.at(tx23UC).trackletUIndex, tv23));
                                    break;
                                }
                            }
                        }
                    }
                    if (match) break;
                }
            }
        }
    }
}

//This function finds the valid combination of tracklets in st2 and st3, and solve the left-right ambigurity
bool KalmanFastTracking_Displaced::compareTrackletsSlim(Tracklet *tracklet2, Tracklet *tracklet3, int plane_type, int pass, double slope_threshold, double windowSize)
{
    if (!(tracklet2->hits.size() == 2 && tracklet3->hits.size() == 2))
        return false;

    // Here we will compare the possible X-Z slopes within the station 2 and station 3 tracklets
    Tracklet::linedef line2;
    Tracklet::linedef line3;

    int choiceOfT2 = 5;
    int choiceOfT3 = 5;

    double pos_hit_st3 = (tracklet3->getHit(0).hit.pos + tracklet3->getHit(1).hit.pos) / 2;
    double z_hit_st3 = (z_plane[tracklet3->getHit(0).hit.detectorID] + z_plane[tracklet3->getHit(1).hit.detectorID]) / 2;

    // It is rare, but sometimes, you will have slopes that match coincidentally.  Therefore, I keep track of best two combinations.  This seems to be sufficienct
    double bestSlopeComp = 1.;
    double bestExtrapMatch = 100.;

    for (unsigned int t2 = 0; t2 < tracklet2->possibleLines.size(); t2++)
    {
        double st2ExtrapMatch = std::abs(tracklet2->possibleLines.at(t2).slope * (z_hit_st3 - tracklet2->possibleLines.at(t2).initialZ) + tracklet2->possibleLines.at(t2).initial_pos - pos_hit_st3);
        if (st2ExtrapMatch < std::min(windowSize, bestExtrapMatch))
        { // Is this extrapolation valid?
            for (unsigned int t3 = 0; t3 < tracklet3->possibleLines.size(); t3++)
            {
                double slope_comp = std::abs(tracklet3->possibleLines.at(t3).slope - tracklet2->possibleLines.at(t2).slope);
                if (slope_comp < std::min(bestSlopeComp, slope_threshold))
                {
                    bestSlopeComp = slope_comp;
                    bestExtrapMatch = st2ExtrapMatch;
                    line2 = tracklet2->possibleLines.at(t2);
                    line3 = tracklet3->possibleLines.at(t3);
                    choiceOfT2 = t2;
                    choiceOfT3 = t3;
                }
            }
        }
    }
    if (bestSlopeComp > slope_threshold)
        return false;
    if (bestExtrapMatch > windowSize)
        return false;

    tracklet2->acceptedLine2 = line2;
    tracklet3->acceptedLine3 = line3;
    // Give the station 2 and station 3 tracklets the same tx and ty value.  I could get an X0 and Y0 extrapolation, but that doesn't seem to be strictly necessary.  The X0 and Y0 values are found in the fittracklet function for the combined station 2 + station 3 tracklet

    tracklet2->pos_st2 = line2.initial_pos;
    tracklet3->pos_st3 = line3.initial_pos;

    // extract slope measurement using the long lever-arm between station 2 and 3
    double newSlope = (tracklet3->possibleLines.at(choiceOfT3).initial_pos - tracklet2->possibleLines.at(choiceOfT2).initial_pos) / (tracklet3->possibleLines.at(choiceOfT3).initialZ - tracklet2->possibleLines.at(choiceOfT2).initialZ);

    tracklet2->slope_st2 = newSlope;
    tracklet3->slope_st3 = newSlope;

    for (auto &hit : tracklet2->hits)
        assignSign(hit, choiceOfT2);
    for (auto &hit : tracklet3->hits)
        assignSign(hit, choiceOfT3);

    return true;
}

//This functions is a 3 hits version of the previous function
bool KalmanFastTracking_Displaced::compareTrackletsSlim_3hits(Tracklet *tracklet2, Tracklet *tracklet3, int plane_type, int pass, double slopeComparison, double windowSize)
{

    if (!((tracklet2->hits.size() == 2 && tracklet3->hits.size() == 1) || (tracklet2->hits.size() == 1 && tracklet3->hits.size() == 2)))
        return false;
    double bestMatch = 100.;

    Tracklet *tracklet_2hits = new Tracklet;
    Tracklet *tracklet_1hit = new Tracklet;
    if (tracklet2->hits.size() == 2)
    {
        tracklet_2hits = tracklet2;
        tracklet_1hit = tracklet3;
    }
    else
    {
        tracklet_2hits = tracklet3;
        tracklet_1hit = tracklet2;
    }
    unsigned int bestT = 5;
    for (unsigned int t = 0; t < tracklet_2hits->possibleLines.size(); t++)
    {
        double extrapolation = tracklet_2hits->possibleLines.at(t).slope * (z_plane[tracklet_1hit->getHit(0).hit.detectorID] - tracklet_2hits->possibleLines.at(t).initialZ) + tracklet_2hits->possibleLines.at(t).initial_pos;
        if (std::abs(tracklet_1hit->getHit(0).hit.pos - extrapolation) < std::min(bestMatch, windowSize))
        {
            bestMatch = std::abs(tracklet_1hit->getHit(0).hit.pos - extrapolation);
            bestT = t;
        }
    }
    if (bestMatch > 99. || bestT == 5)
    {
        return false;
    }
    else
    {
        double newSlope = (tracklet_1hit->getHit(0).hit.pos - tracklet_2hits->possibleLines.at(bestT).initial_pos) / (z_plane[tracklet_1hit->getHit(0).hit.detectorID] - tracklet_2hits->possibleLines.at(bestT).initialZ);

        tracklet2->acceptedLine2 = tracklet_2hits->possibleLines.at(bestT);
        tracklet3->acceptedLine3 = tracklet_2hits->possibleLines.at(bestT);

        tracklet2->slope_st2 = newSlope;
        tracklet3->slope_st3 = newSlope;

        for (auto &hit : tracklet_2hits->hits)
            assignSign(hit, bestT);

        return true;
    }
}

//This function builds full tracklets in st2 and st3
void KalmanFastTracking_Displaced::buildFullBackPartialTracksSlim()
{
    for (unsigned int tx = 0; tx < trackletsInSt23Slim[0].size(); tx++)
    {

        Tracklet tracklet_best;
        Tracklet trackletX = trackletsInSt23Slim[0].at(tx);

        for (auto combo_pair : trackletsInSt23Slim[0].at(tx).matching_combos)
        {
            Tracklet trackletU = trackletsInSt23Slim[1].at(combo_pair.first);
            Tracklet trackletV = trackletsInSt23Slim[2].at(combo_pair.second);

            // quick checks on hit combination compatibility in the three wire slants
            ///XL: I believe there should be a factor cos(14/360*2*pi) for trackletX. I will check it in the future
            if (std::abs((trackletX.slope_st2 - trackletU.slope_st2) - -1. * (trackletX.slope_st2 - trackletV.slope_st2)) > XUVSlopeWin)
                continue;
            if (std::abs((trackletX.pos_st2 - trackletU.pos_st2) - -1. * (trackletX.pos_st2 - trackletV.pos_st2)) > XUVPosWin)
                continue;
            if (std::abs((trackletX.pos_st3 - trackletU.pos_st3) - -1. * (trackletX.pos_st3 - trackletV.pos_st3)) > XUVPosWin)
                continue;
                
            double u_angle_st2 = (std::asin(p_geomSvc->getSintheta(17)) + std::asin(p_geomSvc->getSintheta(18))) / 2;
            double v_angle_st2 = (std::asin(p_geomSvc->getSintheta(13)) + std::asin(p_geomSvc->getSintheta(14))) / 2;
            double u_angle_st3;
            double v_angle_st3;
            if (trackletX.isM == true)
            {
                u_angle_st3 = (std::asin(p_geomSvc->getSintheta(25)) + std::asin(p_geomSvc->getSintheta(26))) / 2;
                v_angle_st3 = (std::asin(p_geomSvc->getSintheta(29)) + std::asin(p_geomSvc->getSintheta(30))) / 2;
            }
            else
            {
                u_angle_st3 = (std::asin(p_geomSvc->getSintheta(19)) + std::asin(p_geomSvc->getSintheta(20))) / 2;
                v_angle_st3 = (std::asin(p_geomSvc->getSintheta(23)) + std::asin(p_geomSvc->getSintheta(24))) / 2;
            }
            double testTX2 = (std::sin(v_angle_st2) * trackletU.slope_st2 - std::sin(u_angle_st2) * trackletV.slope_st2) / std::sin(v_angle_st2 - u_angle_st2);
            double testTX3 = (std::sin(v_angle_st3) * trackletU.slope_st3 - std::sin(u_angle_st3) * trackletV.slope_st3) / std::sin(v_angle_st3 - u_angle_st3);
            double testTY2 = (-std::cos(v_angle_st2) * trackletU.slope_st2 + std::cos(u_angle_st2) * trackletV.slope_st2) / std::sin(v_angle_st2 - u_angle_st2);
            double testTY3 = (-std::cos(v_angle_st3) * trackletU.slope_st3 + std::cos(u_angle_st3) * trackletV.slope_st3) / std::sin(v_angle_st3 - u_angle_st3);

            if (std::abs(testTY2 - testTY3) > YSlopesDiff) continue;

            if (std::abs(trackletX.slope_st2 - testTX2) > XSlopesDiff || std::abs(trackletX.slope_st2 - testTX3) > XSlopesDiff)
                continue; // The x-slopes extracted from the U and V wires should match the slope found from the X wires

            // Let's build a tracklet and assign the x0, tx, and ty parameters (and also the dummy invP parameter).  What we don't know right now is y0.  HOWEVER, IF WE DID SOME MATH TO FIND THE LINE DEFINED BY THE INTERSECTION OF THE PLANES, WE WOULD KNOW THE Y0 VALUE!
            Tracklet tracklet_23 = (trackletX) + (trackletU) + (trackletV);
            tracklet_23.tx = trackletX.slope_st2;
            tracklet_23.ty = (testTY2 + testTY3) / 2.; // take an average of the ty value extracted from the U and V values to hedge your bets
            tracklet_23.invP = 1. / 50.;
            tracklet_23.x0 = trackletX.slope_st2 * (0. - trackletX.z_st2) + trackletX.pos_st2; // simple exrapolation back to z = 0

            tracklet_23.sortHits();
            tracklet_23.y0 = (-1. * p_geomSvc->getDCA((trackletU).getHit(0).hit.detectorID, (trackletU).getHit(0).hit.elementID, tracklet_23.tx, tracklet_23.ty, tracklet_23.x0, 0) / std::abs(std::sin(u_angle_st2)) + p_geomSvc->getDCA((trackletV).getHit(0).hit.detectorID, (trackletV).getHit(0).hit.elementID, tracklet_23.tx, tracklet_23.ty, tracklet_23.x0, 0) / std::abs(std::sin(u_angle_st2))) / 2.; // This is a method to extract y0 based on what the distance of closest approach would be to the st2 U and V hits if the y0 was 0.  A little hacky.  Again, extracting this information from the line defined by the plane intersections would be much smarter and better

            if (tracklet_23.calcChisq_noDrift() > chiSqCut || isnan(tracklet_23.calcChisq_noDrift()))
                continue; // check the chisq when drift distances are not accounted for.  (When using drift distance, the expected precision changes, driving up the chisq given that we only have a "rough" extraction of the trajectory parameters at this point)

            // Assign some parameters that get used later
            tracklet_23.z_st2 = trackletX.z_st2;
            tracklet_23.pos_st2 = trackletX.pos_st2;

            fitTracklet(tracklet_23);
            for (auto &hit : tracklet_23.hits)
            {
                if (hit.sign == 0)
                {
                    hit.sign = 1;
                    fitTracklet(tracklet_23);
                    double dcaPlus = tracklet_23.chisq;
                    hit.sign = -1;
                    fitTracklet(tracklet_23);
                    double dcaMinus = tracklet_23.chisq;
                    if (std::abs(dcaPlus) < std::abs(dcaMinus))
                    {
                        hit.sign = 1;
                    }
                    else
                    {
                        hit.sign = -1;
                    }
                }
            }

            // Remove bad hits if needed;  Right now this doesn't play nicely with this algorithm
            // removeBadHits(tracklet_23);

            fitTracklet(tracklet_23); // A final fit now that all hits have been assigned signs

            if (tracklet_23.chisq > 9000.)
            {
                continue;
            }

            // If current tracklet is better than the best tracklet up-to-now
            if (acceptTracklet(tracklet_23) && tracklet_23 < tracklet_best)
            {
                tracklet_best = tracklet_23;
            }
        } 

        if (acceptTracklet(tracklet_best) && tracklet_best.chisq < st23ChiSqCut)
            trackletsInSt[3].push_back(tracklet_best);
    }
    for (auto &tracklet : trackletsInSt[3])
        checkQuality(tracklet);           /// XL: If add dummy hits here the result will be considerablly different from the original code. I will check it in the future
    reduceTrackletList(trackletsInSt[3]); 
}

//Not displaced version of building global tracklets
void KalmanFastTracking_Displaced::buildGlobalTracks()
{
    double pos_exp[3], window[3];
    for (std::list<Tracklet>::iterator tracklet23 = trackletsInSt[3].begin(); tracklet23 != trackletsInSt[3].end(); ++tracklet23)
    {
        Tracklet tracklet_best[2];
        for (int i = 0; i < 2; ++i) // for two station-1 chambers
        {
            trackletsInSt[0].clear();
            // Calculate the window in station 1
            if (KMAG_ON)
            {
                getSagittaWindowsInSt1(*tracklet23, pos_exp, window, i + 1);
            }
            else
            {
                getExtrapoWindowsInSt1(*tracklet23, pos_exp, window, i + 1);
            }

#ifdef _DEBUG_ON
            LogInfo("Using this back partial: ");
            tracklet23->print();
            for (int j = 0; j < 3; j++)
                LogInfo("Extrapo: " << pos_exp[j] << "  " << window[j]);
#endif

            _timers["global_st1"]->restart();

            buildTrackletsInStation(i + 1, 0, pos_exp, window);

            _timers["global_st1"]->stop();

            _timers["global_link"]->restart();
            Tracklet tracklet_best_prob, tracklet_best_vtx;
            for (std::list<Tracklet>::iterator tracklet1 = trackletsInSt[0].begin(); tracklet1 != trackletsInSt[0].end(); ++tracklet1)
            {
#ifdef _DEBUG_ON
                LogInfo("With this station 1 track:");
                tracklet1->print();
#endif

                Tracklet tracklet_global = (*tracklet23) * (*tracklet1);
                fitTracklet(tracklet_global);

                if (!hodoMask(tracklet_global))
                    continue;

                /// Resolve the left-right with a tight pull cut, then a loose one, then resolve by single projections
                if (!COARSE_MODE)
                {
                    resolveLeftRight(tracklet_global, 75.);
                    resolveLeftRight(tracklet_global, 150.);
                    resolveSingleLeftRight(tracklet_global);
                }

                /// Remove bad hits if needed
                removeBadHits(tracklet_global);

                // Most basic cuts
                if (!acceptTracklet(tracklet_global))
                    continue;

                // Get the tracklets that has the best prob
                if (tracklet_global < tracklet_best_prob)
                    tracklet_best_prob = tracklet_global;

                /// Set vertex information - only applied when KF is enabled
                /// TODO: maybe in the future add a Genfit-based equivalent here, for now leave as is
                if (enable_KF)
                {
                    _timers["global_kalman"]->restart();
                    SRecTrack recTrack = processOneTracklet(tracklet_global);
                    _timers["global_kalman"]->stop();
                    tracklet_global.chisq_vtx = recTrack.getChisqVertex();

                    if (recTrack.isValid() && tracklet_global.chisq_vtx < tracklet_best_vtx.chisq_vtx)
                        tracklet_best_vtx = tracklet_global;
                }

#ifdef _DEBUG_ON
                LogInfo("New tracklet: ");
                tracklet_global.print();

                LogInfo("Current best by prob:");
                tracklet_best_prob.print();

                LogInfo("Comparison I: " << (tracklet_global < tracklet_best_prob));
                LogInfo("Quality I   : " << acceptTracklet(tracklet_global));

                if (enable_KF)
                {
                    LogInfo("Current best by vtx:");
                    tracklet_best_vtx.print();

                    LogInfo("Comparison II: " << (tracklet_global.chisq_vtx < tracklet_best_vtx.chisq_vtx));
                    // LogInfo("Quality II   : " << recTrack.isValid());
                }
#endif
            }
            _timers["global_link"]->stop();

            // The selection logic is, prefer the tracks with best p-value, as long as it's not low-pz
            if (enable_KF && tracklet_best_prob.isValid() > 0 && 1. / tracklet_best_prob.invP > 18.)
            {
                tracklet_best[i] = tracklet_best_prob;
            }
            else if (enable_KF && tracklet_best_vtx.isValid() > 0) // otherwise select the one with best vertex chisq, TODO: maybe add a z-vtx constraint
            {
                tracklet_best[i] = tracklet_best_vtx;
            }
            else if (tracklet_best_prob.isValid() > 0) // then fall back to the default only choice
            {
                tracklet_best[i] = tracklet_best_prob;
            }
        }

        // Merge the tracklets from two stations if necessary
        Tracklet tracklet_merge;
        if (fabs(tracklet_best[0].getMomentum() - tracklet_best[1].getMomentum()) / tracklet_best[0].getMomentum() < MERGE_THRES)
        {
            // Merge the track and re-fit
            tracklet_merge = tracklet_best[0].merge(tracklet_best[1]);
            fitTracklet(tracklet_merge);

#ifdef _DEBUG_ON
            LogInfo("Merging two track candidates with momentum: " << tracklet_best[0].getMomentum() << "  " << tracklet_best[1].getMomentum());
            LogInfo("tracklet_best_1:");
            tracklet_best[0].print();
            LogInfo("tracklet_best_2:");
            tracklet_best[1].print();
            LogInfo("tracklet_merge:");
            tracklet_merge.print();
#endif
        }

        if (tracklet_merge.isValid() > 0 && tracklet_merge < tracklet_best[0] && tracklet_merge < tracklet_best[1])
        {
#ifdef _DEBUG_ON
            LogInfo("Choose merged tracklet");
#endif
            trackletsInSt[4].push_back(tracklet_merge);
        }
        else if (tracklet_best[0].isValid() > 0 && tracklet_best[0] < tracklet_best[1])
        {
#ifdef _DEBUG_ON
            LogInfo("Choose tracklet with station-0");
#endif
            trackletsInSt[4].push_back(tracklet_best[0]);
        }
        else if (tracklet_best[1].isValid() > 0)
        {
#ifdef _DEBUG_ON
            LogInfo("Choose tracklet with station-1");
#endif
            trackletsInSt[4].push_back(tracklet_best[1]);
        }
    }

    trackletsInSt[4].sort();
    reduceTrackletList(trackletsInSt[4]);
    if (resolveStation1Hits())
    {
        reduceTrackletList(trackletsInSt[4]);
    }
}

//Displaced version of building global trackelets
void KalmanFastTracking_Displaced::buildGlobalTracksDisplaced()
{
    for (std::list<Tracklet>::iterator tracklet23 = trackletsInSt[3].begin(); tracklet23 != trackletsInSt[3].end(); ++tracklet23)
    {
        std::vector<Tracklet> possibleGlobalTracks; 
        Tracklet tracklet_best[2];
        for (int i = 0; i < 1; ++i) // for two station-1 chambers //WPM edited so that it only does one chamber.  I think that's all we're using in our simulation...///XL: added back since there are DC1 for E906
        {
            bool validTrackFound = false;
            double valid_track_p = -1;
            Tracklet tracklet_best_prob, tracklet_best_vtx;
            double inv_ps[23];
            inv_ps[0] = 0.0049; /// XL: I find the results overly sensitive to the parameter here. I will check it in the future. These kind of wierd numbers are derived from Patrick's original code
            for (int j = 1; j < 23; j++)
                inv_ps[j] = inv_ps[j - 1] + 0.0097; 
            for (double inv_p : inv_ps)
            {

                int charges[2] = {-1, 1};
                for (int charge : charges)
                {

                    double expXZSlope = (*tracklet23).tx + charge * PT_KICK_KMAG * inv_p;
                    bool hodoFound = false;
                    for (auto hit_index : hitIDs_maskX[0])
                    {
                        if (std::abs(hitAll[hit_index].pos - (expXZSlope * (p_geomSvc->getPlanePosition(hitAll[hit_index].detectorID) - z_plane[i * 6 + 3]) + tracklet23->tx * (z_plane[i * 6 + 3] - tracklet23->z_st2) + tracklet23->pos_st2 - 500 * charge * PT_KICK_KMAG * inv_p)) < 5.)
                        {
                            // if( std::abs( hitAll[hit_index].pos - (expXZSlope*(p_geomSvc->getPlanePosition(hitAll[hit_index].detectorID) - z_plane[i*6+3]) + tracklet23->tx * z_plane[i*6+3] + tracklet23->x0-500*charge*PT_KICK_KMAG*inv_p) ) < 5. ){
                            hodoFound = true;
                            break;
                        }
                    }
                    if (!hodoFound)
                        continue;

                    for (int j = 0; j < 3; j++)
                        trackletsInStSlim[0][j].clear();
                    trackletsInSt[0].clear();

                    (*tracklet23).setCharge(charge);

                    _timers["global_st1"]->restart();
                    double pos_exp[3], window[3];
                    pos_exp[0] = tracklet23->tx * (z_plane[i * 6 + 3] - tracklet23->z_st2) + tracklet23->pos_st2 - 500 * charge * PT_KICK_KMAG * inv_p;
                    ///XL: I believe the commented one is more reasonable, but it will produce some tiny discrepancies to the original code. To be truthful I don't use it here, but I may use it in the future 
                    // pos_exp[0] = tracklet23->tx * z_plane[i*6+3] + tracklet23->x0-500*charge*PT_KICK_KMAG*inv_p;///XinL: I believe the commented one should be more better but it will modify the original results slightly

                    // get expected U and V positions in station 1
                    pos_exp[1] = p_geomSvc->getCostheta(i * 6 + 1) * (expXZSlope * (z_plane[i * 6 + 1] - z_plane[i * 6 + 3]) + pos_exp[0]) + p_geomSvc->getSintheta(i * 6 + 1) * (tracklet23->ty * z_plane[i * 6 + 1] + tracklet23->y0);
                    pos_exp[2] = p_geomSvc->getCostheta(i * 6 + 5) * (expXZSlope * (z_plane[i * 6 + 5] - z_plane[i * 6 + 3]) + pos_exp[0]) + p_geomSvc->getSintheta(i * 6 + 5) * (tracklet23->ty * z_plane[i * 6 + 5] + tracklet23->y0);
                    window[0] = XWinSt1;
                    window[1] = UVWinSt1;
                    window[2] = UVWinSt1;
                    for (int plane_type = 0; plane_type < 3; plane_type++)
                        buildTrackletsInStationSlim(i + 1, 0, plane_type, pos_exp, window);
                    bool doTight = false;

                    if (trackletsInStSlim[0][0].size() * trackletsInStSlim[0][1].size() * trackletsInStSlim[0][2].size() > 3000)
                        doTight = true; // Some station 1 pileup mitigation!

                    if (!(buildTrackletsInStation1(i + 1, 0, expXZSlope, (*tracklet23).ty, (*tracklet23).y0, doTight, pos_exp, window)))
                        continue; // Find station 1 hit combinations in the relevant window

                    _timers["global_st1"]->stop();
                    if (_timers["global_st1"]->get_accumulated_time() / 1000. > 10.)
                        return;
                    _timers["global_link"]->restart();

                    for (std::list<Tracklet>::iterator tracklet1 = trackletsInSt[0].begin(); tracklet1 != trackletsInSt[0].end(); ++tracklet1)
                    { // loop over the potential station 1 tracklets that we found

                        Tracklet tracklet_global = (*tracklet23) * (*tracklet1);
                        tracklet_global.setCharge(charge); // WPM
                        tracklet_global.y0 = (*tracklet23).y0;
                        tracklet_global.x0 = expXZSlope * (0. - z_plane[3]);
                        tracklet_global.tx = expXZSlope;
                        tracklet_global.ty = (*tracklet23).ty;
                        tracklet_global.invP = inv_p;
                        fitTracklet(tracklet_global);

                        /// Remove bad hits if needed
                        // removeBadHits(tracklet_global);

                        // Most basic cuts
                        if (!acceptTracklet(tracklet_global))
                            continue;

                        if (tracklet_global.chisq < 50. && tracklet_global.isValid())
                        {
                            possibleGlobalTracks.push_back(tracklet_global);
                        }

                        // Get the tracklets that has the best prob
                        if (tracklet_global < tracklet_best_prob)
                            tracklet_best_prob = tracklet_global;
                    }

                    _timers["global_link"]->stop();

                    // The selection logic is, prefer the tracks with best p-value, as long as it's not low-pz
                    if (tracklet_best_prob.isValid() > 0) // then fall back to the default only choice
                    {
                        tracklet_best[i] = tracklet_best_prob;
                    }
                }

                if (tracklet_best_prob.chisq < 50. && tracklet_best_prob.isValid())
                {
                    validTrackFound = true;
                    if (valid_track_p < 0)
                        valid_track_p = inv_p;
                }
                if (validTrackFound && (inv_p - valid_track_p) > 4 * (inv_ps[1] - inv_ps[0]) + 0.000001)///XL: I believe there could be more valid momentum when high pile up. 
                    break;
            }
        }

        // Compare the tracklets from two st1 DCs and the merged one, choose the best one
        Tracklet tracklet_merge;
        if (fabs(tracklet_best[0].getMomentum() - tracklet_best[1].getMomentum()) / tracklet_best[0].getMomentum() < MERGE_THRES)
        {
            // Merge the track and re-fit
            tracklet_merge = tracklet_best[0].merge(tracklet_best[1]);
            fitTracklet(tracklet_merge);
        }
        if (tracklet_merge.isValid() > 0 && tracklet_merge < tracklet_best[0] && tracklet_merge < tracklet_best[1])
        {
            globalTracklets.push_back(tracklet_merge);
        }
        else if (tracklet_best[0].isValid() > 0 && tracklet_best[0] < tracklet_best[1])
        {
            globalTracklets.push_back(tracklet_best[0]);
        }
        else if (tracklet_best[1].isValid() > 0)
        {
            globalTracklets.push_back(tracklet_best[1]);
        }

        if (possibleGlobalTracks.size() > 0)
            globalTracklets_resolveSt1.push_back(possibleGlobalTracks);
    }
    for (auto &tracklet : globalTracklets)
    {
        tracklet.addDummyHits();
        checkQuality(tracklet);
    }
    for (auto tracklet : globalTracklets)
    {
        if (tracklet.chisq < 30)
        {
            trackletsInSt[4].push_back(tracklet);
        }
    }
    reduceTrackletList(trackletsInSt[4]);

    if (resolveStation1Hits())
    {
        reduceTrackletList(trackletsInSt[4]);
    }
}

void KalmanFastTracking_Displaced::resolveLeftRight(Tracklet &tracklet, double threshold)
{
#ifdef _DEBUG_ON
    LogInfo("Left right for this track..");
    tracklet.print();
#endif

    // Check if the track has been updated
    bool isUpdated = false;

    // Four possibilities
    int possibility[4][2] = {{1, 1}, {1, -1}, {-1, 1}, {-1, -1}};

    // Total number of hit pairs in this tracklet
    int nPairs = tracklet.hits.size() / 2;

    int nResolved = 0;
    std::list<SignedHit>::iterator hit1 = tracklet.hits.begin();
    std::list<SignedHit>::iterator hit2 = tracklet.hits.begin();
    ++hit2;
    while (true)
    {
#ifdef _DEBUG_ON
        LogInfo(hit1->hit.index << "  " << hit2->sign << " === " << hit2->hit.index << "  " << hit2->sign);
        int detectorID1 = hit1->hit.detectorID;
        int detectorID2 = hit2->hit.detectorID;
        LogInfo("Hit1: " << tracklet.getExpPositionX(z_plane[detectorID1]) * costheta_plane[detectorID1] + tracklet.getExpPositionY(z_plane[detectorID1]) * sintheta_plane[detectorID1] << "  " << hit1->hit.pos + hit1->hit.driftDistance << "  " << hit1->hit.pos - hit1->hit.driftDistance);
        LogInfo("Hit2: " << tracklet.getExpPositionX(z_plane[detectorID2]) * costheta_plane[detectorID2] + tracklet.getExpPositionY(z_plane[detectorID2]) * sintheta_plane[detectorID2] << "  " << hit2->hit.pos + hit2->hit.driftDistance << "  " << hit2->hit.pos - hit2->hit.driftDistance);
#endif

        if (hit1->hit.index > 0 && hit2->hit.index > 0 && hit1->sign * hit2->sign == 0)
        {
            int index_min = -1;
            double pull_min = 1E6;
            for (int i = 0; i < 4; i++)
            {
                double slope_local = (hit1->pos(possibility[i][0]) - hit2->pos(possibility[i][1])) / (z_plane[hit1->hit.detectorID] - z_plane[hit2->hit.detectorID]);
                double inter_local = hit1->pos(possibility[i][0]) - slope_local * z_plane[hit1->hit.detectorID];

                if (fabs(slope_local) > slope_max[hit1->hit.detectorID] || fabs(inter_local) > intersection_max[hit1->hit.detectorID])
                    continue;

                double tx, ty, x0, y0;
                double err_tx, err_ty, err_x0, err_y0;
                if (tracklet.stationID == 7 && hit1->hit.detectorID <= 6)
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

                double slope_exp = costheta_plane[hit1->hit.detectorID] * tx + sintheta_plane[hit1->hit.detectorID] * ty;
                double err_slope = fabs(costheta_plane[hit1->hit.detectorID] * err_tx) + fabs(sintheta_plane[hit2->hit.detectorID] * err_ty);
                double inter_exp = costheta_plane[hit1->hit.detectorID] * x0 + sintheta_plane[hit1->hit.detectorID] * y0;
                double err_inter = fabs(costheta_plane[hit1->hit.detectorID] * err_x0) + fabs(sintheta_plane[hit2->hit.detectorID] * err_y0);

                double pull = sqrt((slope_exp - slope_local) * (slope_exp - slope_local) / err_slope / err_slope + (inter_exp - inter_local) * (inter_exp - inter_local) / err_inter / err_inter);
                if (pull < pull_min)
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

            // LogInfo("Final: " << index_min << "  " << pull_min);
            if (index_min >= 0 && pull_min < threshold) //((tracklet.stationID == 5 && pull_min < 25.) || (tracklet.stationID == 6 && pull_min < 100.)))
            {
                hit1->sign = possibility[index_min][0];
                hit2->sign = possibility[index_min][1];
                isUpdated = true;
            }
        }

        ++nResolved;
        if (nResolved >= nPairs)
            break;

        ++hit1;
        ++hit1;
        ++hit2;
        ++hit2;
    }

    if (isUpdated)
        fitTracklet(tracklet);
}

void KalmanFastTracking_Displaced::resolveSingleLeftRight(Tracklet &tracklet)
{
#ifdef _DEBUG_ON
    LogInfo("Single left right for this track..");
    tracklet.print();
#endif

    // Check if the track has been updated
    bool isUpdated = false;
    for (std::list<SignedHit>::iterator hit_sign = tracklet.hits.begin(); hit_sign != tracklet.hits.end(); ++hit_sign)
    {
        if (hit_sign->hit.index < 0 || hit_sign->sign != 0)
            continue;

        int detectorID = hit_sign->hit.detectorID;
        double pos_exp = tracklet.getExpPositionX(z_plane[detectorID]) * costheta_plane[detectorID] + tracklet.getExpPositionY(z_plane[detectorID]) * sintheta_plane[detectorID];
        hit_sign->sign = pos_exp > hit_sign->hit.pos ? 1 : -1;

        isUpdated = true;
    }

    if (isUpdated)
        fitTracklet(tracklet);
}

bool KalmanFastTracking_Displaced::removeBadHits(Tracklet &tracklet)
{
#ifdef _DEBUG_ON
    LogInfo("Removing hits for this track..");
    tracklet.calcChisq();
    tracklet.print();
#endif

    // Check if the track has beed updated
    int signflipflag[nChamberPlanes];
    for (int i = 0; i < nChamberPlanes; ++i)
        signflipflag[i] = 0;

    bool isUpdated = true;
    while (isUpdated)
    {
        isUpdated = false;
        tracklet.calcChisq();

        SignedHit *hit_remove = nullptr;
        SignedHit *hit_neighbour = nullptr;
        double res_remove1 = -1.;
        double res_remove2 = -1.;
        for (std::list<SignedHit>::iterator hit_sign = tracklet.hits.begin(); hit_sign != tracklet.hits.end(); ++hit_sign)
        {
            if (hit_sign->hit.index < 0)
                continue;

            int detectorID = hit_sign->hit.detectorID;
            double res_curr = fabs(tracklet.residual[detectorID - 1]);
            if (res_remove1 < res_curr)
            {
                res_remove1 = res_curr;
                res_remove2 = fabs(tracklet.residual[detectorID - 1] - 2. * hit_sign->sign * hit_sign->hit.driftDistance);
                hit_remove = &(*hit_sign);

                std::list<SignedHit>::iterator iter = hit_sign;
                hit_neighbour = detectorID % 2 == 0 ? &(*(--iter)) : &(*(++iter));
            }
        }
        if (hit_remove == nullptr)
            continue;
        if (hit_remove->sign == 0 && tracklet.isValid() > 0)
            continue; // if sign is undecided, and chisq is OKay, then pass

        double cut = hit_remove->sign == 0 ? hit_remove->hit.driftDistance + resol_plane[hit_remove->hit.detectorID] : resol_plane[hit_remove->hit.detectorID];
        if (res_remove1 > cut)
        {
#ifdef _DEBUG_ON
            LogInfo("Dropping this hit: " << res_remove1 << "  " << res_remove2 << "   " << signflipflag[hit_remove->hit.detectorID - 1] << "  " << cut);
            hit_remove->hit.print();
            hit_neighbour->hit.print();
#endif

            // can only be changed less than twice
            if (res_remove2 < cut && signflipflag[hit_remove->hit.detectorID - 1] < 2)
            {
                hit_remove->sign = -hit_remove->sign;
                hit_neighbour->sign = 0;
                ++signflipflag[hit_remove->hit.detectorID - 1];
#ifdef _DEBUG_ON
                LogInfo("Only changing the sign.");
#endif
            }
            else
            {
                // Set the index of the hit to be removed to -1 so it's not used anymore
                // also set the sign assignment of the neighbour hit to 0 (i.e. undecided)
                if ((tracklet.nHits[0] + tracklet.nHits[1] + tracklet.nHits[2]) < 15)
                {
                    return false;
                }
                hit_remove->hit.index = -1;
                hit_neighbour->sign = 0;
                int planeType = p_geomSvc->getPlaneType(hit_remove->hit.detectorID);
                if (planeType == 1)
                {
                    --tracklet.nHits[0];
                }
                else if (planeType == 2)
                {
                    --tracklet.nHits[1];
                }
                else
                {
                    --tracklet.nHits[2];
                }

                // If both hit pairs are not included, the track can be rejected
                if (hit_neighbour->hit.index < 0)
                {
#ifdef _DEBUG_ON
                    LogInfo("Both hits in a view are missing! Will exit the bad hit removal...");
#endif
                    return false;
                }
            }
            isUpdated = true;
        }

        if (isUpdated)
        {
            if ((tracklet.nHits[0] + tracklet.nHits[1] + tracklet.nHits[2]) < 14)
            {
                return false;
            }
            fitTracklet(tracklet);
            resolveSingleLeftRight(tracklet);
        }
    }
    return true;
}

void KalmanFastTracking_Displaced::resolveLeftRight(SRawEvent::hit_pair hpair, int &LR1, int &LR2)
{
    LR1 = 0;
    LR2 = 0;

    // If either hit is missing, no left-right can be assigned
    if (hpair.first < 0 || hpair.second < 0)
    {
        return;
    }

    int possibility[4][2] = {{1, 1}, {1, -1}, {-1, 1}, {-1, -1}};
    int nResolved = 0;
    for (int i = 0; i < 4; i++)
    {
        if (nResolved > 1)
            break;

        int hitID1 = hpair.first;
        int hitID2 = hpair.second;
        double slope_local = (hitAll[hitID1].pos + possibility[i][0] * hitAll[hitID1].driftDistance - hitAll[hitID2].pos - possibility[i][1] * hitAll[hitID2].driftDistance) / (z_plane[hitAll[hitID1].detectorID] - z_plane[hitAll[hitID2].detectorID]);
        double intersection_local = hitAll[hitID1].pos + possibility[i][0] * hitAll[hitID1].driftDistance - slope_local * z_plane[hitAll[hitID1].detectorID];

        // LogInfo(i << "  " << nResolved << "  " << slope_local << "  " << intersection_local);
        if (fabs(slope_local) < slope_max[hitAll[hitID1].detectorID] && fabs(intersection_local) < intersection_max[hitAll[hitID1].detectorID])
        {
            nResolved++;
            LR1 = possibility[i][0];
            LR2 = possibility[i][1];
        }
    }

    if (nResolved > 1)
    {
        LR1 = 0;
        LR2 = 0;
    }

    // LogInfo("Final: " << LR1 << "  " << LR2);
}

void KalmanFastTracking_Displaced::buildTrackletsInStation(int stationID, int listID, double *pos_exp, double *window)
{
#ifdef _DEBUG_ON
    LogInfo("Building tracklets in station " << stationID);
#endif

    // actuall ID of the tracklet lists
    int sID = stationID - 1;

    // Extract the X, U, V hit pairs
    std::list<SRawEvent::hit_pair> pairs_X, pairs_U, pairs_V;
    if (pos_exp == nullptr)
    {
        pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0]);
        pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1]);
        pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2]);
    }
    else
    {
        // Note that in pos_exp[], index 0 stands for X, index 1 stands for U, index 2 stands for V
        pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0], pos_exp[0], window[0]);
        pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1], pos_exp[1], window[1]);
        pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2], pos_exp[2], window[2]);
    }

#ifdef _DEBUG_ON
    LogInfo("Hit pairs in this event: ");
    for (std::list<SRawEvent::hit_pair>::iterator iter = pairs_X.begin(); iter != pairs_X.end(); ++iter)
        LogInfo("X :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
    for (std::list<SRawEvent::hit_pair>::iterator iter = pairs_U.begin(); iter != pairs_U.end(); ++iter)
        LogInfo("U :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
    for (std::list<SRawEvent::hit_pair>::iterator iter = pairs_V.begin(); iter != pairs_V.end(); ++iter)
        LogInfo("V :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
#endif

    if (pairs_X.empty() || pairs_U.empty() || pairs_V.empty())
    {
#ifdef _DEBUG_ON
        LogInfo("Not all view has hits in station " << stationID);
#endif
        return;
    }

    // X-U combination first, then add V pairs
    for (std::list<SRawEvent::hit_pair>::iterator xiter = pairs_X.begin(); xiter != pairs_X.end(); ++xiter)
    {
        // U projections from X plane
        double x_pos = xiter->second >= 0 ? 0.5 * (hitAll[xiter->first].pos + hitAll[xiter->second].pos) : hitAll[xiter->first].pos;
        double u_min = x_pos * u_costheta[sID] - u_win[sID];
        double u_max = u_min + 2. * u_win[sID];

#ifdef _DEBUG_ON
        LogInfo("Trying X hits " << xiter->first << "  " << xiter->second << "  " << hitAll[xiter->first].elementID << " at " << x_pos);
        LogInfo("U plane window:" << u_min << "  " << u_max);
#endif
        for (std::list<SRawEvent::hit_pair>::iterator uiter = pairs_U.begin(); uiter != pairs_U.end(); ++uiter)
        {
            double u_pos = uiter->second >= 0 ? 0.5 * (hitAll[uiter->first].pos + hitAll[uiter->second].pos) : hitAll[uiter->first].pos;
#ifdef _DEBUG_ON
            LogInfo("Trying U hits " << uiter->first << "  " << uiter->second << "  " << hitAll[uiter->first].elementID << " at " << u_pos);
#endif
            if (u_pos < u_min || u_pos > u_max)
                continue;

            // V projections from X and U plane
            double z_x = xiter->second >= 0 ? z_plane_x[sID] : z_plane[hitAll[xiter->first].detectorID];
            double z_u = uiter->second >= 0 ? z_plane_u[sID] : z_plane[hitAll[uiter->first].detectorID];
            double z_v = z_plane_v[sID];
            double v_win1 = spacing_plane[hitAll[uiter->first].detectorID] * 2. * u_costheta[sID];
            double v_win2 = fabs((z_u + z_v - 2. * z_x) * u_costheta[sID] * TX_MAX);
            double v_win3 = fabs((z_v - z_u) * u_sintheta[sID] * TY_MAX);
            double v_win = v_win1 + v_win2 + v_win3 + 2. * spacing_plane[hitAll[uiter->first].detectorID];
            double v_min = 2 * x_pos * u_costheta[sID] - u_pos - v_win;
            double v_max = v_min + 2. * v_win;

#ifdef _DEBUG_ON
            LogInfo("V plane window:" << v_min << "  " << v_max);
#endif
            for (std::list<SRawEvent::hit_pair>::iterator viter = pairs_V.begin(); viter != pairs_V.end(); ++viter)
            {
                double v_pos = viter->second >= 0 ? 0.5 * (hitAll[viter->first].pos + hitAll[viter->second].pos) : hitAll[viter->first].pos;
#ifdef _DEBUG_ON
                LogInfo("Trying V hits " << viter->first << "  " << viter->second << "  " << hitAll[viter->first].elementID << " at " << v_pos);
#endif
                if (v_pos < v_min || v_pos > v_max)
                    continue;

                // Now add the tracklet
                int LR1 = 0;
                int LR2 = 0;
                Tracklet tracklet_new;
                tracklet_new.stationID = stationID;

                // resolveLeftRight(*xiter, LR1, LR2);
                if (xiter->first >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[xiter->first], LR1));
                    tracklet_new.nHits[0]++;
                }
                if (xiter->second >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[xiter->second], LR2));
                    tracklet_new.nHits[0]++;
                }

                tracklet_new.getSlopes(hitAll[xiter->first], hitAll[xiter->second], 0); // Here, we find the four possible X-Z lines

                // resolveLeftRight(*uiter, LR1, LR2);
                if (uiter->first >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[uiter->first], LR1));
                    tracklet_new.nHits[1]++;
                }
                if (uiter->second >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[uiter->second], LR2));
                    tracklet_new.nHits[1]++;
                }

                tracklet_new.getSlopes(hitAll[uiter->first], hitAll[uiter->second], 1); // find the four possible U-Z lines

                // resolveLeftRight(*viter, LR1, LR2);
                if (viter->first >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[viter->first], LR1));
                    tracklet_new.nHits[2]++;
                }
                if (viter->second >= 0)
                {
                    tracklet_new.hits.push_back(SignedHit(hitAll[viter->second], LR2));
                    tracklet_new.nHits[2]++;
                }

                tracklet_new.getSlopes(hitAll[viter->first], hitAll[viter->second], 2); // find the four possible V-Z lines

                tracklet_new.sortHits();
                if (tracklet_new.isValid() == 0) // TODO: What IS THIS?
                {
                    fitTracklet(tracklet_new); // This is where the original DCA minimization is performed
                }
                else
                {
                    continue;
                }

#ifdef _DEBUG_ON
                tracklet_new.print();
#endif
                if (acceptTracklet(tracklet_new))
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

    // Reduce the tracklet list and add dummy hits
    // reduceTrackletList(trackletsInSt[listID]);
    for (std::list<Tracklet>::iterator iter = trackletsInSt[listID].begin(); iter != trackletsInSt[listID].end(); ++iter)
    {
        iter->addDummyHits();
    }

    // Only retain the best 1000 tracklets if exceeded
    if (trackletsInSt[listID].size() > 1000)
    {
        trackletsInSt[listID].sort();
        trackletsInSt[listID].resize(1000);
    }
}

//*//This function builds tracklets in st1
bool KalmanFastTracking_Displaced::buildTrackletsInStation1(int stationID, int listID, double expXZSlope, double expYSlope, double y0, bool tight, double *pos_exp, double *window)
{
#ifdef _DEBUG_ON
    LogInfo("Building tracklets in station " << stationID);
#endif

    double slopeComparisonSt1 = (tight ? 0.05 : slopeCompSt1);
    double slope_exp[3];
    slope_exp[0] = expXZSlope;
    slope_exp[1] = p_geomSvc->getCostheta(1) * expXZSlope + p_geomSvc->getSintheta(1) * expYSlope;
    slope_exp[2] = p_geomSvc->getCostheta(5) * expXZSlope + p_geomSvc->getSintheta(5) * expYSlope;

    std::list<Tracklet> valid_tracklets[3];
    for (int plane_type = 0; plane_type < 3; plane_type++)
    {
        for (Tracklet tracklet : trackletsInStSlim[0][plane_type])
        {
            if (tracklet.hits.size() < 2 && tight)
                continue;
            if (tracklet.hits.size() == 1)
            {
                valid_tracklets[plane_type].push_back(tracklet);
            }
            else
            { // if there are 2 hits, we should be able to determine the hit signs and roughly check if the hit combination gives the right expected slopes
                double bestSlopeDiff = slopeComparisonSt1;
                int bestTracklet = 5;
                for (int t = 0; t < tracklet.possibleLines.size(); t++)
                {
                    double slopeDiff = std::abs(tracklet.possibleLines.at(t).slope - slope_exp[plane_type]);
                    if (slopeDiff < bestSlopeDiff)
                    {
                        bestSlopeDiff = slopeDiff;
                        bestTracklet = t;
                    }
                }
                if (bestTracklet > 4)
                {
                    continue;
                }
                else
                {
                    for (auto &hit : tracklet.hits) assignSign(hit, bestTracklet);
                    valid_tracklets[plane_type].push_back(tracklet);
                }
            }
        }
    }

    for (auto trackletX : valid_tracklets[0])
    {
        for (auto trackletU : valid_tracklets[1])
        {
            for (auto trackletV : valid_tracklets[2])
            {
                if (trackletX.hits.size() < 2 && trackletU.hits.size() < 2 && trackletV.hits.size() < 2)
                    continue; // We require at least 4 out of 6 possible hits
                Tracklet tracklet_new_Station1 = trackletX + trackletU + trackletV;

                tracklet_new_Station1.stationID = stationID;
                tracklet_new_Station1.y0 = y0;
                tracklet_new_Station1.ty = expYSlope;
                tracklet_new_Station1.x0 = -1 * expXZSlope * z_plane[3] + pos_exp[0];
                tracklet_new_Station1.tx = expXZSlope;
                tracklet_new_Station1.sortHits();

                if (tracklet_new_Station1.calcChisq_noDrift() > 300 || isnan(tracklet_new_Station1.calcChisq_noDrift()))
                    continue;
                fitTracklet(tracklet_new_Station1);

                // hit sign assignment for station 1
                ///XL: This part is a bit different from the similar chunk in building back tracklets. I on't know if it is bug and I will check it in the future
                for (auto &hit : tracklet_new_Station1.hits)
                {
                    if (hit.sign == 0)
                    {
                        hit.sign = 1;
                        double dcaPlus = tracklet_new_Station1.calcChisq(); 
                        hit.sign = -1;
                     
                        double dcaMinus = tracklet_new_Station1.calcChisq();
                        if (std::abs(dcaPlus) < std::abs(dcaMinus))
                        {
                            hit.sign = 1;
                        }
                        else
                        {
                            hit.sign = -1;
                        }
                    }
                }

                fitTracklet(tracklet_new_Station1);

                if (tracklet_new_Station1.chisq < 20.)
                {
                    tracklet_new_Station1.addDummyHits();
                    trackletsInSt[listID].push_back(tracklet_new_Station1);
                }
            }
        }
    }
    // reduceTrackletList(trackletsInSt[listID]);

    if (trackletsInSt[listID].size() > 0)
        return true;
}
//*/

bool KalmanFastTracking_Displaced::acceptTracklet(Tracklet &tracklet)
{ /// Maybe add more information here in the future
    return tracklet.isValid();
}

/// XL: This function is problematic when used in buildGlobalTracks. Firstly, the stationID=3
/// statement means it will do nothing at all. Secondly, possibleLines is not assgined to global tracks. 
bool KalmanFastTracking_Displaced::hodoMask(Tracklet &tracklet)
{
    // LogInfo(tracklet.stationID);
    if (TRACK_ELECTRONS && (tracklet.stationID == 4 || tracklet.stationID == 5))
        return true; // Skip of hodoscope checks for station 3 tracks in the electron-tracking setup.  Could actually probably extrapolate backwards the station 2 hodoscope, now that I get an accurate X-Z slope in station 3
    int nHodoHits = 0;
    // Performing a valid extrapolation in X-Z for station 2 tracklets.  This should be improved.  Currently carries around the old fudge factor
    if (tracklet.stationID == 3)
    {
        for (std::vector<int>::iterator stationID = stationIDs_mask[tracklet.stationID - 1].begin(); stationID != stationIDs_mask[tracklet.stationID - 1].end(); ++stationID)
        {
            bool masked = false;
            std::cout << "hey!-3 " << tracklet.possibleLines.size() << std::endl;
            for (std::list<int>::iterator iter = hitIDs_mask[*stationID - 1].begin(); iter != hitIDs_mask[*stationID - 1].end(); ++iter)
            {
                int detectorID = hitAll[*iter].detectorID;
                int elementID = hitAll[*iter].elementID;

                int idx1 = detectorID - nChamberPlanes - 1;
                int idx2 = elementID - 1;

                double factor = tracklet.stationID == nChamberPlanes / 6 - 2 ? 5. : 3.; // special for station-2, based on real data tuning
                double xfudge = tracklet.stationID < nStations - 1 ? 0.5 * (x_mask_max[idx1][idx2] - x_mask_min[idx1][idx2]) : 0.15 * (x_mask_max[idx1][idx2] - x_mask_min[idx1][idx2]);
                double z_hodo = z_mask[idx1];
                std::cout << "hey! " << tracklet.possibleLines.size() << std::endl;
                for (unsigned int pl = 0; pl < tracklet.possibleLines.size(); pl++)
                {
                    double extrapolation = tracklet.possibleLines.at(pl).slope * (z_hodo - tracklet.possibleLines.at(pl).initialZ) + tracklet.possibleLines.at(pl).initial_pos;
                    double err_x = std::abs(factor * extrapolation + xfudge);
                    double x_min = x_mask_min[idx1][idx2] - err_x;
                    double x_max = x_mask_max[idx1][idx2] + err_x;
                    if (extrapolation > x_min && extrapolation < x_max)
                    {
                        masked = true;
                        break;
                    }
                }
            }
            if (!masked)
                return false;
        }
    }

    if (tracklet.stationID > 5)
    {
        for (std::vector<int>::iterator stationID = stationIDs_mask[tracklet.stationID - 1].begin(); stationID != stationIDs_mask[tracklet.stationID - 1].end(); ++stationID)
        {
            bool masked = false;
            for (std::list<int>::iterator iter = hitIDs_mask[*stationID - 1].begin(); iter != hitIDs_mask[*stationID - 1].end(); ++iter)
            {
                int detectorID = hitAll[*iter].detectorID;
                int elementID = hitAll[*iter].elementID;

                int idx1 = detectorID - nChamberPlanes - 1;
                int idx2 = elementID - 1;

                double factor = tracklet.stationID == nChamberPlanes / 6 - 2 ? 5. : 3.; // special for station-2, based on real data tuning
                double xfudge = tracklet.stationID < nStations - 1 ? 0.5 * (x_mask_max[idx1][idx2] - x_mask_min[idx1][idx2]) : 0.15 * (x_mask_max[idx1][idx2] - x_mask_min[idx1][idx2]);
                double z_hodo = z_mask[idx1];
                double x_hodo = tracklet.getExpPositionX(z_hodo);
                double y_hodo = tracklet.getExpPositionY(z_hodo);
                double err_x = factor * tracklet.getExpPosErrorX(z_hodo) + xfudge;
                double err_y = factor * tracklet.getExpPosErrorY(z_hodo);

                double x_min = x_mask_min[idx1][idx2] - err_x;
                double x_max = x_mask_max[idx1][idx2] + err_x;
                double y_min = y_mask_min[idx1][idx2] - err_y;
                double y_max = y_mask_max[idx1][idx2] + err_y;

#ifdef _DEBUG_ON
                LogInfo(*iter);
                hitAll[*iter].print();
                LogInfo(nHodoHits << "/" << stationIDs_mask[tracklet.stationID - 1].size() << ":  " << z_hodo << "  " << x_hodo << " +/- " << err_x << "  " << y_hodo << " +/-" << err_y << " : " << x_min << "  " << x_max << "  " << y_min << "  " << y_max);
#endif
                if (x_hodo > x_min && x_hodo < x_max && y_hodo > y_min && y_hodo < y_max)
                {
                    nHodoHits++;
                    masked = true;

                    if (TRACK_ELECTRONS && tracklet.stationID > 5)
                        return true; // Once the first hodoscope hit is found (at z=1420cm), the combined tracklet passes for electron tracks

                    break;
                }
            }

            if (!masked)
                return false;
        }
    }

#ifdef _DEBUG_ON
    LogInfo(tracklet.stationID << "  " << nHodoHits << "  " << stationIDs_mask[tracklet.stationID - 1].size());
#endif
    return true;
}

//prop tube infomation is not used in this algorithm
bool KalmanFastTracking_Displaced::muonID_search(Tracklet &tracklet)
{
    // Set the cut value on multiple scattering
    // multiple scattering: sigma = 0.0136*sqrt(L/L0)*(1. + 0.038*ln(L/L0))/P, L = 1m, L0 = 1.76cm
    double cut = 0.03;
    if (tracklet.stationID == nStations)
    {
        double cut_the = MUID_THE_P0 * tracklet.invP;
        double cut_emp = MUID_EMP_P0 + MUID_EMP_P1 / tracklet.invP + MUID_EMP_P2 / tracklet.invP / tracklet.invP;
        cut = MUID_REJECTION * (cut_the > cut_emp ? cut_the : cut_emp);
    }

    double slope[2] = {tracklet.tx, tracklet.ty};
    double pos_absorb[2] = {tracklet.getExpPositionX(MUID_Z_REF), tracklet.getExpPositionY(MUID_Z_REF)};
    PropSegment *segs[2] = {&(tracklet.seg_x), &(tracklet.seg_y)};
    for (int i = 0; i < 2; ++i)
    {
        // this shorting circuting can only be done to X-Z, Y-Z needs more complicated thing
        // if(i == 0 && segs[i]->getNHits() > 2 && segs[i]->isValid() > 0 && fabs(slope[i] - segs[i]->a) < cut) continue;

        segs[i]->init();
        for (int j = 0; j < 4; ++j)
        {
            int index = detectorIDs_muid[i][j] - nChamberPlanes - 1;
            double pos_ref = j < 2 ? pos_absorb[i] : segs[i]->getPosRef(pos_absorb[i] + slope[i] * (z_ref_muid[i][j] - MUID_Z_REF));
            double pos_exp = slope[i] * (z_mask[index] - z_ref_muid[i][j]) + pos_ref;

            if (!p_geomSvc->isInPlane(detectorIDs_muid[i][j], tracklet.getExpPositionX(z_mask[index]), tracklet.getExpPositionY(z_mask[index])))
                continue;

            double win_tight = cut * (z_mask[index] - z_ref_muid[i][j]);
            win_tight = win_tight > 2.54 ? win_tight : 2.54;
            double win_loose = win_tight * 2;
            double dist_min = 1E6;
            for (std::list<int>::iterator iter = hitIDs_muid[i][j].begin(); iter != hitIDs_muid[i][j].end(); ++iter)
            {
                double pos = hitAll[*iter].pos;
                double dist = pos - pos_exp;
                if (dist < -win_loose)
                    continue;
                if (dist > win_loose)
                    break;

                double dist_l = fabs(pos - hitAll[*iter].driftDistance - pos_exp);
                double dist_r = fabs(pos + hitAll[*iter].driftDistance - pos_exp);
                dist = dist_l < dist_r ? dist_l : dist_r;
                if (dist < dist_min)
                {
                    dist_min = dist;
                    if (dist < win_tight)
                    {
                        segs[i]->hits[j].hit = hitAll[*iter];
                        segs[i]->hits[j].sign = fabs(pos - hitAll[*iter].driftDistance - pos_exp) < fabs(pos + hitAll[*iter].driftDistance - pos_exp) ? -1 : 1;
                    }
                }
            }
        }
        segs[i]->fit();

        // this shorting circuting can only be done to X-Z, Y-Z needs more complicated thing
        // if(i == 0 && !(segs[i]->isValid() > 0 && fabs(slope[i] - segs[i]->a) < cut)) return false;
    }

    muonID_hodoAid(tracklet);
    if (segs[0]->getNHits() + segs[1]->getNHits() >= MUID_MINHITS)
    {
        return true;
    }
    else if (segs[1]->getNHits() == 1 || segs[1]->getNPlanes() == 1)
    {
        return segs[1]->nHodoHits >= 2;
    }
    return false;
}

//prop tube infomation is not used in this algorithm
bool KalmanFastTracking_Displaced::muonID_comp(Tracklet &tracklet)
{
    // Set the cut value on multiple scattering
    // multiple scattering: sigma = 0.0136*sqrt(L/L0)*(1. + 0.038*ln(L/L0))/P, L = 1m, L0 = 1.76cm
    double cut = 0.03;
    if (tracklet.stationID == nStations)
    {
        double cut_the = MUID_THE_P0 * tracklet.invP;
        double cut_emp = MUID_EMP_P0 + MUID_EMP_P1 / tracklet.invP + MUID_EMP_P2 / tracklet.invP / tracklet.invP;
        cut = MUID_REJECTION * (cut_the > cut_emp ? cut_the : cut_emp);
    }
#ifdef _DEBUG_ON
    LogInfo("Muon ID cut is: " << cut << " rad.");
#endif

    double slope[2] = {tracklet.tx, tracklet.ty};
    PropSegment *segs[2] = {&(tracklet.seg_x), &(tracklet.seg_y)};

    for (int i = 0; i < 2; ++i)
    {
#ifdef _DEBUG_ON
        if (i == 0)
            LogInfo("Working in X-Z:");
        if (i == 1)
            LogInfo("Working in Y-Z:");
#endif

        double pos_ref = i == 0 ? tracklet.getExpPositionX(MUID_Z_REF) : tracklet.getExpPositionY(MUID_Z_REF);
        if (segs[i]->getNHits() > 2 && segs[i]->isValid() > 0 && fabs(slope[i] - segs[i]->a) < cut && fabs(segs[i]->getExpPosition(MUID_Z_REF) - pos_ref) < MUID_R_CUT)
        {
#ifdef _DEBUG_ON
            LogInfo("Muon ID are already avaiable!");
#endif
            continue;
        }

        for (std::list<PropSegment>::iterator iter = propSegs[i].begin(); iter != propSegs[i].end(); ++iter)
        {
#ifdef _DEBUG_ON
            LogInfo("Testing this prop segment, with ref pos = " << pos_ref << ", slope_ref = " << slope[i]);
            iter->print();
#endif
            if (fabs(iter->a - slope[i]) < cut && fabs(iter->getExpPosition(MUID_Z_REF) - pos_ref) < MUID_R_CUT)
            {
                *(segs[i]) = *iter;
#ifdef _DEBUG_ON
                LogInfo("Accepted!");
#endif
                break;
            }
        }

        if (segs[i]->isValid() == 0)
            return false;
    }

    if (segs[0]->getNHits() + segs[1]->getNHits() < MUID_MINHITS)
        return false;
    return true;
}

//prop tube infomation is not used in this algorithm
bool KalmanFastTracking_Displaced::muonID_hodoAid(Tracklet &tracklet)
{
    double win = 0.03;
    double factor = 5.;
    if (tracklet.stationID == nStations)
    {
        double win_the = MUID_THE_P0 * tracklet.invP;
        double win_emp = MUID_EMP_P0 + MUID_EMP_P1 / tracklet.invP + MUID_EMP_P2 / tracklet.invP / tracklet.invP;
        win = MUID_REJECTION * (win_the > win_emp ? win_the : win_emp);
        factor = 3.;
    }

    PropSegment *segs[2] = {&(tracklet.seg_x), &(tracklet.seg_y)};
    for (int i = 0; i < 2; ++i)
    {
        segs[i]->nHodoHits = 0;
        for (std::list<int>::iterator iter = hitIDs_muidHodoAid[i].begin(); iter != hitIDs_muidHodoAid[i].end(); ++iter)
        {
            int detectorID = hitAll[*iter].detectorID;
            int elementID = hitAll[*iter].elementID;

            int idx1 = detectorID - nChamberPlanes - 1;
            int idx2 = elementID - 1;

            double z_hodo = z_mask[idx1];
            double x_hodo = tracklet.getExpPositionX(z_hodo);
            double y_hodo = tracklet.getExpPositionY(z_hodo);
            double err_x = factor * tracklet.getExpPosErrorX(z_hodo) + win * (z_hodo - MUID_Z_REF);
            double err_y = factor * tracklet.getExpPosErrorY(z_hodo) + win * (z_hodo - MUID_Z_REF);

            err_x = err_x / (x_mask_max[idx1][idx2] - x_mask_min[idx1][idx2]) > 0.25 ? 0.25 * err_x / (x_mask_max[idx1][idx2] - x_mask_min[idx1][idx2]) : err_x;
            err_y = err_y / (y_mask_max[idx1][idx2] - y_mask_min[idx1][idx2]) > 0.25 ? 0.25 * err_y / (y_mask_max[idx1][idx2] - y_mask_min[idx1][idx2]) : err_y;

            double x_min = x_mask_min[idx1][idx2] - err_x;
            double x_max = x_mask_max[idx1][idx2] + err_x;
            double y_min = y_mask_min[idx1][idx2] - err_y;
            double y_max = y_mask_max[idx1][idx2] + err_y;

            if (x_hodo > x_min && x_hodo < x_max && y_hodo > y_min && y_hodo < y_max)
            {
                segs[i]->hodoHits[segs[i]->nHodoHits++] = hitAll[*iter];
                if (segs[i]->nHodoHits > 4)
                    break;
            }
        }
    }

    return true;
}

//prop tube infomation is not used in this algorithm
void KalmanFastTracking_Displaced::buildPropSegments()
{
#ifdef _DEBUG_ON
    LogInfo("Building prop. tube segments");
#endif

    for (int i = 0; i < 2; ++i)
    {
        propSegs[i].clear();

        // note for prop tubes superID index starts from 4
        std::list<SRawEvent::hit_pair> pairs_forward = rawEvent->getPartialHitPairsInSuperDetector(superIDs[i + 5][0]);
        std::list<SRawEvent::hit_pair> pairs_backward = rawEvent->getPartialHitPairsInSuperDetector(superIDs[i + 5][1]);

#ifdef _DEBUG_ON
        // std::cout << "superID: " << superIDs[i+5][0] << ", " << superIDs[i+5][1] << std::endl;
        for (std::list<SRawEvent::hit_pair>::iterator iter = pairs_forward.begin(); iter != pairs_forward.end(); ++iter)
            LogInfo("Forward: " << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << "  " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
        for (std::list<SRawEvent::hit_pair>::iterator iter = pairs_backward.begin(); iter != pairs_backward.end(); ++iter)
            LogInfo("Backward: " << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << "  " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
#endif

        for (std::list<SRawEvent::hit_pair>::iterator fiter = pairs_forward.begin(); fiter != pairs_forward.end(); ++fiter)
        {
#ifdef _DEBUG_ON
            LogInfo("Trying forward pair " << fiter->first << "  " << fiter->second);
#endif
            for (std::list<SRawEvent::hit_pair>::iterator biter = pairs_backward.begin(); biter != pairs_backward.end(); ++biter)
            {
#ifdef _DEBUG_ON
                LogInfo("Trying backward pair " << biter->first << "  " << biter->second);
#endif

                PropSegment seg;

                // Note that the backward plane comes as the first in pair
                if (fiter->first >= 0)
                    seg.hits[1] = SignedHit(hitAll[fiter->first], 0);
                if (fiter->second >= 0)
                    seg.hits[0] = SignedHit(hitAll[fiter->second], 0);
                if (biter->first >= 0)
                    seg.hits[3] = SignedHit(hitAll[biter->first], 0);
                if (biter->second >= 0)
                    seg.hits[2] = SignedHit(hitAll[biter->second], 0);

#ifdef _DEBUG_ON
                seg.print();
#endif
                seg.fit();
#ifdef _DEBUG_ON
                seg.print();
#endif

                if (seg.isValid() > 0)
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

///ROOT minimizer to find the best tracklet
int KalmanFastTracking_Displaced::fitTracklet(Tracklet &tracklet)
{
    tracklet_curr = tracklet;

    // idx = 0, using simplex; idx = 1 using migrad
    int idx = 1;
#ifdef _ENABLE_MULTI_MINI
    if (tracklet.stationID < nStations - 1)
        idx = 0;
#endif

    minimizer[idx]->SetLimitedVariable(0, "tx", tracklet.tx, 0.001, -TX_MAX, TX_MAX);
    minimizer[idx]->SetLimitedVariable(1, "ty", tracklet.ty, 0.001, -TY_MAX, TY_MAX);
    minimizer[idx]->SetLimitedVariable(2, "x0", tracklet.x0, 0.1, -X0_MAX, X0_MAX);
    minimizer[idx]->SetLimitedVariable(3, "y0", tracklet.y0, 0.1, -Y0_MAX, Y0_MAX);
    if (KMAG_ON)
    {
        minimizer[idx]->SetLimitedVariable(4, "invP", tracklet.invP, 0.001 * tracklet.invP, INVP_MIN, INVP_MAX);
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

    if (KMAG_ON && tracklet.stationID == nStations)
    {
        tracklet.invP = minimizer[idx]->X()[4];
        tracklet.err_invP = minimizer[idx]->Errors()[4];
    }

    tracklet.chisq = minimizer[idx]->MinValue();

    int status = minimizer[idx]->Status();

    return status;
}

///remove similar tracklets
int KalmanFastTracking_Displaced::reduceTrackletList(std::list<Tracklet> &tracklets)
{
    std::list<Tracklet> targetList;
    tracklets.sort();
    while (!tracklets.empty())
    {
        if ((tracklets.front().nHits[0] + tracklets.front().nHits[1] + tracklets.front().nHits[2]) < 14 && tracklets.front().stationID == nStations)
        {
            tracklets.pop_front();
            continue;
        }
        targetList.push_back(tracklets.front());
        tracklets.pop_front();
        for (std::list<Tracklet>::iterator iter = tracklets.begin(); iter != tracklets.end();)
        {
            if (iter->stationID < nStations)
            {
                /// XL: I don't understand the different condition statement here
                if (iter->similarity(targetList.back()) && std::abs(targetList.back().pos_st2 - iter->pos_st2) < 3. && std::abs(targetList.back().pos_st3 - iter->pos_st3) < 3. && std::abs(targetList.back().tx - iter->tx) < 0.01 && std::abs(targetList.back().ty - iter->ty) < 0.01)
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
            else
            {
                if (iter->similarity(targetList.back()) && std::abs(targetList.back().tx - iter->tx) < 0.01 && std::abs(targetList.back().ty - iter->ty) < 0.01 && std::abs(targetList.back().x0 - iter->x0) < 10 && std::abs(targetList.back().y0 - iter->y0) < 10)
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
    }
    tracklets.assign(targetList.begin(), targetList.end());
    return 0;
}

//Used in not displaced version to find the expected hit position in st1
void KalmanFastTracking_Displaced::getExtrapoWindowsInSt1(Tracklet &tracklet, double *pos_exp, double *window, int st1ID)
{
    if (tracklet.stationID != nStations - 1)
    {
        for (int i = 0; i < 3; i++)
        {
            pos_exp[i] = 9999.;
            window[i] = 0.;
        }
        return;
    }

    for (int i = 0; i < 3; i++)
    {
        int detectorID = (st1ID - 1) * 6 + 2 * i + 2;
        int idx = p_geomSvc->getPlaneType(detectorID) - 1;

        double z_st1 = z_plane[detectorID];
        double x_st1 = tracklet.getExpPositionX(z_st1);
        double y_st1 = tracklet.getExpPositionY(z_st1);
        double err_x = tracklet.getExpPosErrorX(z_st1);
        double err_y = tracklet.getExpPosErrorY(z_st1);

        pos_exp[idx] = p_geomSvc->getUinStereoPlane(detectorID, x_st1, y_st1);
        window[idx] = 5. * (fabs(costheta_plane[detectorID] * err_x) + fabs(sintheta_plane[detectorID] * err_y));
    }
}

//Used in not displaced version to find the expected hit position in st1
void KalmanFastTracking_Displaced::getSagittaWindowsInSt1(Tracklet &tracklet, double *pos_exp, double *window, int st1ID)
{
    if (tracklet.stationID != nStations - 1)
    {
        for (int i = 0; i < 3; i++)
        {
            pos_exp[i] = 9999.;
            window[i] = 0.;
        }
        return;
    }

    double z_st3 = z_plane[tracklet.hits.back().hit.detectorID];
    double x_st3 = tracklet.getExpPositionX(z_st3);
    double y_st3 = tracklet.getExpPositionY(z_st3);

    // For U, X, and V planes
    for (int i = 0; i < 3; i++)
    {
        int detectorID = (st1ID - 1) * 6 + 2 * i + 2;
        int idx = p_geomSvc->getPlaneType(detectorID) - 1;

        if (!(idx >= 0 && idx < 3))
            continue;

        double pos_st3 = p_geomSvc->getUinStereoPlane(s_detectorID[idx], x_st3, y_st3);

        double z_st1 = z_plane[detectorID];
        double z_st2 = z_plane[s_detectorID[idx]];
        double x_st2 = tracklet.getExpPositionX(z_st2);
        double y_st2 = tracklet.getExpPositionY(z_st2);
        double pos_st2 = p_geomSvc->getUinStereoPlane(s_detectorID[idx], x_st2, y_st2);

        double s2_target = pos_st2 - pos_st3 * (z_st2 - Z_TARGET) / (z_st3 - Z_TARGET);
        double s2_dump = pos_st2 - pos_st3 * (z_st2 - Z_DUMP) / (z_st3 - Z_DUMP);

        double pos_exp_target = SAGITTA_TARGET_CENTER * s2_target + pos_st3 * (z_st1 - Z_TARGET) / (z_st3 - Z_TARGET);
        double pos_exp_dump = SAGITTA_DUMP_CENTER * s2_dump + pos_st3 * (z_st1 - Z_DUMP) / (z_st3 - Z_DUMP);
        double win_target = fabs(s2_target * SAGITTA_TARGET_WIDTH);
        double win_dump = fabs(s2_dump * SAGITTA_DUMP_WIDTH);

        double p_min = std::min(pos_exp_target - win_target, pos_exp_dump - win_dump);
        double p_max = std::max(pos_exp_target + win_target, pos_exp_dump + win_dump);

        pos_exp[idx] = 0.5 * (p_max + p_min);
        window[idx] = 0.5 * (p_max - p_min);
    }
}

void KalmanFastTracking_Displaced::printAtDetectorBack(int stationID, std::string outputFileName)
{
    TCanvas c1;

    std::vector<double> x, y, dx, dy;
    for (std::list<Tracklet>::iterator iter = trackletsInSt[stationID].begin(); iter != trackletsInSt[stationID].end(); ++iter)
    {
        double z = p_geomSvc->getPlanePosition(iter->stationID * 6);
        x.push_back(iter->getExpPositionX(z));
        y.push_back(iter->getExpPositionY(z));
        dx.push_back(iter->getExpPosErrorX(z));
        dy.push_back(iter->getExpPosErrorY(z));
    }

    TGraphErrors gr(x.size(), &x[0], &y[0], &dx[0], &dy[0]);
    gr.SetMarkerStyle(8);

    // Add detector frames
    std::vector<double> x_f, y_f, dx_f, dy_f;
    x_f.push_back(p_geomSvc->getPlaneCenterX(stationID * 6 + 6));
    y_f.push_back(p_geomSvc->getPlaneCenterY(stationID * 6 + 6));
    dx_f.push_back(p_geomSvc->getPlaneScaleX(stationID * 6 + 6) * 0.5);
    dy_f.push_back(p_geomSvc->getPlaneScaleY(stationID * 6 + 6) * 0.5);

    if (stationID == 2)
    {
        x_f.push_back(p_geomSvc->getPlaneCenterX(stationID * 6 + 12));
        y_f.push_back(p_geomSvc->getPlaneCenterY(stationID * 6 + 12));
        dx_f.push_back(p_geomSvc->getPlaneScaleX(stationID * 6 + 12) * 0.5);
        dy_f.push_back(p_geomSvc->getPlaneScaleY(stationID * 6 + 12) * 0.5);
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

SRecTrack KalmanFastTracking_Displaced::processOneTracklet(Tracklet &tracklet)
{
    // tracklet.print();
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

    // Resolve left-right based on the current solution, re-fit if anything changed
    // resolveLeftRight(kmtrk);
    if (fitTrack(kmtrk) && kmtrk.isValid())
    {
        SRecTrack strack = kmtrk.getSRecTrack();

        // Set trigger road ID
        TriggerRoad road(tracklet);
        strack.setTriggerRoad(road.getRoadID());

        // Set prop tube slopes
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

bool KalmanFastTracking_Displaced::fitTrack(KalmanTrack &kmtrk)
{
    if (kmtrk.getNodeList().empty())
        return false;

    if (kmfitter->processOneTrack(kmtrk) == 0)
    {
        return false;
    }
    kmfitter->updateTrack(kmtrk);

    return true;
}

void KalmanFastTracking_Displaced::resolveLeftRight(KalmanTrack &kmtrk)
{
    bool isUpdated = false;

    std::list<int>::iterator hitID = kmtrk.getHitIndexList().begin();
    for (std::list<Node>::iterator node = kmtrk.getNodeList().begin(); node != kmtrk.getNodeList().end();)
    {
        if (*hitID == 0)
        {
            double x_hit = node->getSmoothed().get_x();
            double y_hit = node->getSmoothed().get_y();
            double pos_hit = p_geomSvc->getUinStereoPlane(node->getHit().detectorID, x_hit, y_hit);

            int sign = 0;
            if (pos_hit > node->getHit().pos)
            {
                sign = 1;
            }
            else
            {
                sign = -1;
            }

            // update the node list
            TMatrixD m(1, 1), dm(1, 1);
            m[0][0] = node->getHit().pos + sign * node->getHit().driftDistance;
            dm[0][0] = p_geomSvc->getPlaneResolution(node->getHit().detectorID) * p_geomSvc->getPlaneResolution(node->getHit().detectorID);
            node->setMeasurement(m, dm);
            *hitID = sign * node->getHit().index;

            isUpdated = true;
        }

        ++node;
        ++hitID;
    }

    if (isUpdated)
        fitTrack(kmtrk);
}

void KalmanFastTracking_Displaced::printTimers()
{
    std::cout << "KalmanFastTracking_Displaced::printTimers: " << std::endl;
    std::cout << "================================================================" << std::endl;
    std::cout << "Tracklet St2                " << _timers["st2"]->get_accumulated_time() / 1000. << " sec" << std::endl;
    std::cout << "Tracklet St3                " << _timers["st3"]->get_accumulated_time() / 1000. << " sec" << std::endl;
    std::cout << "Connection                  " << _timers["connect"]->get_accumulated_time() / 1000. << " sec" << std::endl;
    std::cout << "Tracklet St23               " << _timers["st23"]->get_accumulated_time() / 1000. << " sec" << std::endl;
    std::cout << "Tracklet Global             " << _timers["global"]->get_accumulated_time() / 1000. << " sec" << std::endl;
    std::cout << "  Global St1                " << _timers["global_st1"]->get_accumulated_time() / 1000. << " sec" << std::endl;
    std::cout << "  Global Link               " << _timers["global_link"]->get_accumulated_time() / 1000. << " sec" << std::endl;
    std::cout << "  Global Kalman             " << _timers["global_kalman"]->get_accumulated_time() / 1000. << " sec" << std::endl;
    std::cout << "Tracklet Kalman             " << _timers["kalman"]->get_accumulated_time() / 1000. << " sec" << std::endl;
    std::cout << "================================================================" << std::endl;
}

void KalmanFastTracking_Displaced::chi2fit(int n, double x[], double y[], double &a, double &b)
{
    double sum = 0.;
    double sx = 0.;
    double sy = 0.;
    double sxx = 0.;
    double syy = 0.;
    double sxy = 0.;

    for (int i = 0; i < n; ++i)
    {
        ++sum;
        sx += x[i];
        sy += y[i];
        sxx += (x[i] * x[i]);
        syy += (y[i] * y[i]);
        sxy += (x[i] * y[i]);
    }

    double det = sum * sxx - sx * sx;
    if (fabs(det) < 1E-20)
    {
        a = 0.;
        b = 0.;

        return;
    }

    a = (sum * sxy - sx * sy) / det;
    b = (sy * sxx - sxy * sx) / det;
}

//*
// For the case when both st2 and st3 hit combos have 2 hits

void KalmanFastTracking_Displaced::cutAdjuster(int numCombos, int pass)
{

    slopeComp = .15;
    windowSize = 15;
    reqHits = 2;
    slopeCompSt1 = .15;
    XWinSt1 = 1.25;
    UVWinSt1 = 1.5;
    hodoXWin = 7;
    hodoUVWin = 30;
    hodoXDIFFWin = 7;
    hodoUVDIFFWin = 7;
    XUVSlopeWin = 0.04;
    XUVPosWin = 9.;
    YSlopesDiff = 0.007;
    XSlopesDiff = 0.007;
    chiSqCut = 100;
    st23ChiSqCut = 15;

    if (!adjusted)
    {

        if (rawEvent->getNHitsInD0() < 400 && rawEvent->getNHitsInD2() < 300 && rawEvent->getNHitsInD3p() < 300 && rawEvent->getNHitsInD3m() < 300)
        {
            slopeComp = .25;
            windowSize = 27.5;
            hodoXWin = 14;
            hodoUVWin = 40;
            hodoXDIFFWin = 10;
            hodoUVDIFFWin = 13;
            chiSqCut = m_chiSqCut;
            st23ChiSqCut = 30;
        }

        if ((numCombos < 5000000 && pass == 2) || (numCombos < 500000 && pass == 1))
        {
            slopeComp = .15;
            windowSize = 15;
            reqHits = 2;
            slopeCompSt1 = .15;
            XWinSt1 = 1.25;
            UVWinSt1 = 1.5;
            hodoXWin = 7;
            hodoUVWin = 30;
            hodoXDIFFWin = 7;
            hodoUVDIFFWin = 7;
            XUVSlopeWin = 0.04;
            XUVPosWin = 9.;
            YSlopesDiff = 0.007;
            XSlopesDiff = 0.007;
            chiSqCut = 100;
            st23ChiSqCut = 15;

            if (rawEvent->getNHitsInD0() < 500 && rawEvent->getNHitsInD2() < 300 && rawEvent->getNHitsInD3p() < 300 && rawEvent->getNHitsInD3m() < 300)
            {
                slopeComp = .25;
                windowSize = 27.5;
                hodoXWin = 14;
                hodoUVWin = 40;
                hodoXDIFFWin = 10;
                hodoUVDIFFWin = 13;
                chiSqCut = m_chiSqCut;
                st23ChiSqCut = 30;

                if (rawEvent->getNHitsInD0() < 400)
                {
                    reqHits = 1;
                    slopeCompSt1 = .25;
                    hodoUVDIFFWin = 17;
                }
            }

            if ((numCombos < 500000 && pass == 2) || (numCombos < 50000 && pass == 1))
            {

                slopeComp = .15;
                windowSize = 15;
                reqHits = 2;
                slopeCompSt1 = .15;
                XWinSt1 = 1.25;
                UVWinSt1 = 1.5;
                hodoXWin = 7;
                hodoUVWin = 30;
                hodoXDIFFWin = 7;
                hodoUVDIFFWin = 7;
                XUVSlopeWin = 0.04;
                XUVPosWin = 9.;
                YSlopesDiff = 0.007;
                XSlopesDiff = 0.007;
                chiSqCut = 100;
                st23ChiSqCut = 15;

                if (rawEvent->getNHitsInD0() < 1000 && rawEvent->getNHitsInD2() < 300 && rawEvent->getNHitsInD3p() < 300 && rawEvent->getNHitsInD3m() < 300)
                {
                    slopeComp = .25;
                    windowSize = 27.5;
                    hodoXWin = 14;
                    hodoUVWin = 40;
                    hodoXDIFFWin = 10;
                    hodoUVDIFFWin = 13;
                    chiSqCut = m_chiSqCut;
                    st23ChiSqCut = 30;

                    if (rawEvent->getNHitsInD0() < 500)
                    {
                        reqHits = 1;
                        slopeCompSt1 = .25;
                        hodoUVDIFFWin = 17;

                        if (rawEvent->getNHitsInD0() < 400 && rawEvent->getNHitsInD2() < 300 && rawEvent->getNHitsInD3p() < 300 && rawEvent->getNHitsInD3m() < 300)
                        {
                            slopeComp = .3;
                            windowSize = 50.;
                            slopeCompSt1 = m_slopeComparisonSt1;
                            UVWinSt1 = 2.25;
                            hodoUVDIFFWin = 20;
                            st23ChiSqCut = 50;
                        }
                    }
                }
            }
        }
    }
    else
    {

        if ((numCombos < 5000000 && pass == 2) || (numCombos < 500000 && pass == 1))
        {

            if (rawEvent->getNHitsInD0() < 400 && rawEvent->getNHitsInD2() < 300 && rawEvent->getNHitsInD3p() < 300 && rawEvent->getNHitsInD3m() < 300)
            {
                slopeComp = .25;
                windowSize = 27.5;
                hodoXWin = 14;
                hodoUVWin = 40;
                hodoXDIFFWin = 10;
                hodoUVDIFFWin = 13;
                chiSqCut = m_chiSqCut;
                st23ChiSqCut = 30;
            }

            if ((numCombos < 500000 && pass == 2) || (numCombos < 50000 && pass == 1))
            {

                if (rawEvent->getNHitsInD0() < 500 && rawEvent->getNHitsInD2() < 300 && rawEvent->getNHitsInD3p() < 300 && rawEvent->getNHitsInD3m() < 300)
                {
                    slopeComp = .25;
                    windowSize = 27.5;
                    hodoXWin = 14;
                    hodoUVWin = 40;
                    hodoXDIFFWin = 10;
                    hodoUVDIFFWin = 13;
                    chiSqCut = m_chiSqCut;
                    st23ChiSqCut = 30;

                    if (rawEvent->getNHitsInD0() < 400)
                    {
                        reqHits = 1;
                        slopeCompSt1 = .25;
                        hodoUVDIFFWin = 17;
                    }
                }
            }
        }
    }
}

///XL: I don't understand what's going on here so I didn't do much modification
bool KalmanFastTracking_Displaced::resolveStation1Hits()
{
    if (trackletsInSt[4].size() < 2)
        return false;
    if (trackletsInSt[4].size() == 2)
    {
        if (globalTracklets_resolveSt1.size() < 2)
            return false;
        double similarity = (*trackletsInSt[4].begin()).similarity_st1(*(++trackletsInSt[4].begin()));
        if (similarity <= 0.2)
            return false;
        // double xint, yint;
        // checkIntercepts(*trackletsInSt[4].begin(), *(++trackletsInSt[4].begin()), xint, yint);

        double bestChiSqCombo = 10000;
        Tracklet best1, best2;
        bool validTracks = false;

        for (unsigned int i = 0; i < globalTracklets_resolveSt1.at(0).size(); i++)
        {
            for (unsigned int j = 0; j < globalTracklets_resolveSt1.at(1).size(); j++)
            {
                double xintTest, yintTest;
                checkIntercepts(globalTracklets_resolveSt1.at(0).at(i), globalTracklets_resolveSt1.at(1).at(j), xintTest, yintTest);
                if (globalTracklets_resolveSt1.at(0).at(i).similarity_st1(globalTracklets_resolveSt1.at(1).at(j)) < 0.2)
                {
                    if (sqrt(globalTracklets_resolveSt1.at(0).at(i).calcChisq_st1_squares() * globalTracklets_resolveSt1.at(0).at(i).calcChisq_st1_squares() + globalTracklets_resolveSt1.at(1).at(j).calcChisq_st1_squares() * globalTracklets_resolveSt1.at(1).at(j).calcChisq_st1_squares()) < bestChiSqCombo && (xintTest < 610 || yintTest < 610) && std::abs(xintTest - yintTest) < 100)
                    {
                        best1 = globalTracklets_resolveSt1.at(0).at(i);
                        best2 = globalTracklets_resolveSt1.at(1).at(j);
                        validTracks = true;
                        bestChiSqCombo = sqrt(globalTracklets_resolveSt1.at(0).at(i).calcChisq_st1_squares() * globalTracklets_resolveSt1.at(0).at(i).calcChisq_st1_squares() + globalTracklets_resolveSt1.at(1).at(j).calcChisq_st1_squares() * globalTracklets_resolveSt1.at(1).at(j).calcChisq_st1_squares());
                    }
                }
            }
        }

        if (validTracks)
        {
            if (!(best1.similarityAllowed(best2)))
            {
                trackletsInSt[4].clear();
                return false;
            }
            best1.addDummyHits();
            best2.addDummyHits();
            checkQuality(best1);
            checkQuality(best2);

            trackletsInSt[4].clear();
            trackletsInSt[4].push_back(best1);
            trackletsInSt[4].push_back(best2);
        }
        else
        {
            trackletsInSt[4].clear();
        }

        return validTracks;
    }

    else
    {
        bool noSimilarities = true;
        std::vector<std::pair<int, int>> similarIndices;

        std::vector<Tracklet> st1_tracklets;
        for (std::list<Tracklet>::iterator tr = trackletsInSt[4].begin(); tr != trackletsInSt[4].end(); ++tr)
        {
            st1_tracklets.push_back((*tr));
        }

        std::vector<Tracklet> goodTracklets;
        std::vector<Tracklet> badTracklets;

        for (int i = 0; i < st1_tracklets.size(); i++)
        {
            for (int j = 0; j < st1_tracklets.size(); j++)
            {
                if (i == j)
                    continue;
                double similarity = st1_tracklets.at(i).similarity_st1(st1_tracklets.at(j));
                if (similarity > 0.2)
                    noSimilarities = false;
                similarIndices.push_back(std::make_pair(i, j));
            }
        }

        if (noSimilarities)
            return false;

        int numTries = 0;
        while (!noSimilarities && numTries < 10)
        {
            noSimilarities = true;
            int ii, jj;
            for (int ii = 0; ii < st1_tracklets.size(); ii++)
            {
                bool unique = true;
                for (int jj = 0; jj < st1_tracklets.size(); jj++)
                {
                    if (ii == jj)
                        continue;
                    double similarity = st1_tracklets.at(ii).similarity_st1(st1_tracklets.at(jj));
                    if (similarity > 0.2)
                    {
                        noSimilarities = false;
                        unique = false;
                        break;
                    }
                }
                if (!unique)
                    break;
            }

            int ind1 = ii;
            int ind2 = jj;

            double bestChiSqCombo = 10000;
            Tracklet best1, best2;
            bool validTracks = false;

            if (globalTracklets_resolveSt1.size() < ind1 || globalTracklets_resolveSt1.size() < ind2)
            {
                // something went wrong here... clear all tracks just in case
                trackletsInSt[4].clear();
                return false;
            }

            for (unsigned int i = 0; i < globalTracklets_resolveSt1.at(ind1).size(); i++)
            {
                for (unsigned int j = 0; j < globalTracklets_resolveSt1.at(ind2).size(); j++)
                {
                    double xintTest, yintTest;
                    checkIntercepts(globalTracklets_resolveSt1.at(ind1).at(i), globalTracklets_resolveSt1.at(ind2).at(j), xintTest, yintTest);
                    if (globalTracklets_resolveSt1.at(ind1).at(i).similarity_st1(globalTracklets_resolveSt1.at(ind2).at(j)) < 0.2)
                    {
                        if (sqrt(globalTracklets_resolveSt1.at(ind1).at(i).calcChisq_st1_squares() * globalTracklets_resolveSt1.at(ind1).at(i).calcChisq_st1_squares() + globalTracklets_resolveSt1.at(ind2).at(j).calcChisq_st1_squares() * globalTracklets_resolveSt1.at(ind2).at(j).calcChisq_st1_squares()) < bestChiSqCombo && (xintTest < 610 || yintTest < 610) && std::abs(xintTest - yintTest) < 100 && globalTracklets_resolveSt1.at(ind1).at(i).similarityAllowed(globalTracklets_resolveSt1.at(ind2).at(j)))
                        {
                            best1 = globalTracklets_resolveSt1.at(ind1).at(i);
                            best2 = globalTracklets_resolveSt1.at(ind2).at(j);
                            validTracks = true;
                            bestChiSqCombo = sqrt(globalTracklets_resolveSt1.at(ind1).at(i).calcChisq_st1_squares() * globalTracklets_resolveSt1.at(ind1).at(i).calcChisq_st1_squares() + globalTracklets_resolveSt1.at(ind2).at(j).calcChisq_st1_squares() * globalTracklets_resolveSt1.at(ind2).at(j).calcChisq_st1_squares());
                        }
                    }
                }
            }

            if (validTracks)
            {
                best1.addDummyHits();
                best2.addDummyHits();

                checkQuality(best1);
                checkQuality(best2);
                st1_tracklets.at(ind1) = best1;
                st1_tracklets.at(ind2) = best2;
            }
            else if (numTries == 9)
            {
                trackletsInSt[4].clear();
                return false;
            }
            else
            {

                std::vector<Tracklet> tempVec;
                for (int i = 0; i < st1_tracklets.size(); i++)
                {
                    if (i != ind2)
                    {
                        tempVec.push_back(st1_tracklets.at(i));
                    }
                }
                st1_tracklets.clear();
                for (int i = 0; i < tempVec.size(); i++)
                {
                    st1_tracklets.push_back(tempVec.at(i));
                }
                tempVec.clear();
            }
            numTries++;
        }

        trackletsInSt[4].clear();
        for (int i = 0; i < st1_tracklets.size(); i++)
        {
            trackletsInSt[4].push_back(st1_tracklets.at(i));
            if (Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
            {
                st1_tracklets.at(i).print();
            }
        }
    }
}

void KalmanFastTracking_Displaced::checkIntercepts(Tracklet tracklet1, Tracklet tracklet2, double &xint, double &yint)
{

    double tx_st1_1, x0_st1_1, ty_st1_1, y0_st1_1;
    tracklet1.getXZInfoInSt1(tx_st1_1, x0_st1_1);
    ty_st1_1 = tracklet1.ty;
    y0_st1_1 = tracklet1.y0;
    double tx_st1_2, x0_st1_2, ty_st1_2, y0_st1_2;
    tracklet2.getXZInfoInSt1(tx_st1_2, x0_st1_2);
    ty_st1_2 = tracklet2.ty;
    y0_st1_2 = tracklet2.y0;

    xint = (x0_st1_1 - x0_st1_2) / (tx_st1_2 - tx_st1_1);
    yint = (y0_st1_1 - y0_st1_2) / (ty_st1_2 - ty_st1_1);

    return;
}

//See if flipping the "side" of a wire that hits are on can improve the quality
void KalmanFastTracking_Displaced::checkQuality(Tracklet &tracklet)
{
    double testChisq = tracklet.chisq;
    bool gettingBetter = true;
    int passes = 0;
    while (testChisq > 20. && gettingBetter && passes < 5)
    {
        checkSigns(tracklet); 
        if (tracklet.chisq < testChisq)
        {
            testChisq = tracklet.chisq;
        }
        else
        {
            gettingBetter = false;
        }
        passes++;
    }
}

void KalmanFastTracking_Displaced::checkSigns(Tracklet &tracklet)
{
    double compChiSq = tracklet.chisq;
    for (std::list<SignedHit>::iterator fliped_hit = tracklet.hits.begin(); fliped_hit != tracklet.hits.end(); ++fliped_hit)
    {
        if (fliped_hit->hit.index < 0)
            continue;

        int detectorID = fliped_hit->hit.detectorID;
        int index = detectorID - 1;
        tracklet.residual[index] = fliped_hit->sign * fabs(fliped_hit->hit.driftDistance) - p_geomSvc->getDCA(detectorID, fliped_hit->hit.elementID, tracklet.tx, tracklet.ty, tracklet.x0, tracklet.y0);

        if (std::abs(-1 * fliped_hit->sign * fabs(fliped_hit->hit.driftDistance) - p_geomSvc->getDCA(detectorID, fliped_hit->hit.elementID, tracklet.tx, tracklet.ty, tracklet.x0, tracklet.y0)) < std::abs(2. * tracklet.residual[index]))
        {
            fliped_hit->sign = -fliped_hit->sign;
            fitTracklet(tracklet);
            if (tracklet.calcChisq() < compChiSq)
            {
                compChiSq = tracklet.chisq;
            }
            else
            {
                fliped_hit->sign = -fliped_hit->sign;
                fitTracklet(tracklet);
            }
        }
    }
}

//Assign sign of the drift distance
void KalmanFastTracking_Displaced::assignSign(SignedHit &hit, int line_index)
{
    int is_odd = hit.hit.detectorID % 2;
    if (is_odd)
        hit.sign = 1 - 2 * (line_index % 2);
    else
        hit.sign = 1 - 2 * (line_index / 2);
}
