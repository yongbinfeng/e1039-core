#ifndef _KALMAN_FAST_TRACKLETTING_NEW_HODO_2_H
#define _KALMAN_FAST_TRACKLETTING_NEW_HODO_2_H
#include "KalmanFastTracking_NEW_HODO_2.h"

class KalmanFastTrackletting_NEW_HODO_2 : public KalmanFastTracking_NEW_HODO_2
{
  double TX_MAX;
  double TY_MAX;
  double X0_MAX;
  double Y0_MAX;

public:
    explicit KalmanFastTrackletting_NEW_HODO_2(const PHField* field, const TGeoManager *geom, bool flag = true);
    virtual ~KalmanFastTrackletting_NEW_HODO_2();

    virtual int setRawEvent(SRawEvent* event_input);

    virtual void buildTrackletsInStation(int stationID, int listID, double* pos_exp = nullptr, double* window = nullptr);
};

#endif // _KALMAN_FAST_TRACKLETTING_NEW_HODO_2_H
