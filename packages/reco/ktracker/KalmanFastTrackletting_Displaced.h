#ifndef _KALMAN_FAST_TRACKLETTING_Displaced_H
#define _KALMAN_FAST_TRACKLETTING_Displaced_H
#include "KalmanFastTracking_Displaced.h"

class KalmanFastTrackletting_Displaced : public KalmanFastTracking_Displaced
{
  double TX_MAX;
  double TY_MAX;
  double X0_MAX;
  double Y0_MAX;

public:
    explicit KalmanFastTrackletting_Displaced(const PHField* field, const TGeoManager *geom, bool flag = true);
    virtual ~KalmanFastTrackletting_Displaced();

    virtual int setRawEvent(SRawEvent* event_input);

    virtual void buildTrackletsInStation(int stationID, int listID, double* pos_exp = nullptr, double* window = nullptr);
};

#endif // _KALMAN_FAST_TRACKLETTING_Displaced_H
