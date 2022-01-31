//#include <iomanip>
#include <phool/phool.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQSpillMap.h>
#include <interface_main/SQHitVector.h>
#include <ktracker/SRawEvent.h>
#include "UtilSRawEvent.h"
using namespace std;

bool UtilSRawEvent::SetEvent(SRawEvent* sraw, const SQEvent* evt, const bool do_assert)
{
  if (!evt) {
    if (do_assert) {
      cout << PHWHERE << ": SQEvent == 0.  Abort." << endl;
      exit(1);
    }
    sraw->setEventInfo(0, 0, 0);
    sraw->setTriggerBits(0);
    return false;
  }
  int run_id   = evt->get_run_id();
  int spill_id = evt->get_spill_id();
  int event_id = evt->get_event_id();
  sraw->setEventInfo(run_id, spill_id, event_id);

  int triggers[10];
  for(int i = SQEvent::NIM1; i <= SQEvent::MATRIX5; ++i) {
    triggers[i] = evt->get_trigger(static_cast<SQEvent::TriggerMask>(i));
  }
  sraw->setTriggerBits(triggers);
  return true;
}

bool UtilSRawEvent::SetSpill(SRawEvent* sraw, const SQSpill* sp, const bool do_assert)
{
  if (!sp) {
    if (do_assert) {
      cout << PHWHERE << ": SQSpill == 0.  Abort." << endl;
      exit(1);
    }
    sraw->setTargetPos(1);
    return false;
  }
  sraw->setTargetPos(sp->get_target_pos());
  return true;
}

bool UtilSRawEvent::SetHit(SRawEvent* sraw, const SQHitVector* hit_vec, std::map<int, size_t>* hitID_idx, const bool do_assert)
{
  sraw->emptyHits();
  if (hitID_idx) hitID_idx->clear();

  if (!hit_vec) {
    if (do_assert) {
      cout << PHWHERE << ": SQHitVector == 0.  Abort." << endl;
      exit(1);
    }
    return false;
  }
  for(size_t idx = 0; idx < hit_vec->size(); ++idx) {
    const SQHit* sq_hit = hit_vec->at(idx);
    if (hitID_idx) (*hitID_idx)[sq_hit->get_hit_id()] = idx;

    Hit h;
    h.index         = sq_hit->get_hit_id();
    h.detectorID    = sq_hit->get_detector_id();
    h.elementID     = sq_hit->get_element_id();
    h.tdcTime       = sq_hit->get_tdc_time();
    h.driftDistance = fabs(sq_hit->get_drift_distance()); //MC L-R info removed here
    h.pos           = sq_hit->get_pos();
    //std::cout<<"in sethit.  index = "<<sq_hit->get_hit_id()<<", detID = "<<sq_hit->get_detector_id()<<", elementID = "<<sq_hit->get_element_id()<<", tdcTime = "<<sq_hit->get_tdc_time()<<", driftDistance = "<<fabs(sq_hit->get_drift_distance())<<", pos = "<<sq_hit->get_pos()<<", is_in_time = "<<sq_hit->is_in_time()<<std::endl; //WPM
    if(sq_hit->is_in_time()) h.setInTime();
    //std::cout<<"and second in time check = "<<sq_hit->is_in_time()<<" and "<<h.isInTime()<<std::endl; //WPM
    sraw->insertHit(h);
  }
  sraw->reIndex(true);
  return true;
}

bool UtilSRawEvent::SetTriggerHit(SRawEvent* sraw, const SQHitVector* hit_vec, std::map<int, size_t>* hitID_idx, const bool do_assert)
{
  sraw->emptyTriggerHits();
  if (hitID_idx) hitID_idx->clear();

  if (!hit_vec) {
    if (do_assert) {
      cout << PHWHERE << ": SQHitVector == 0.  Abort." << endl;
      exit(1);
    }
    return false;
  }
  for(size_t idx = 0; idx < hit_vec->size(); ++idx) {
    const SQHit* sq_hit = hit_vec->at(idx);
    if (hitID_idx) (*hitID_idx)[sq_hit->get_hit_id()] = idx;

    Hit h;
    h.index         = sq_hit->get_hit_id();
    h.detectorID    = sq_hit->get_detector_id();
    h.elementID     = sq_hit->get_element_id();
    h.tdcTime       = sq_hit->get_tdc_time();
    h.driftDistance = fabs(sq_hit->get_drift_distance()); //MC L-R info removed here
    h.pos           = sq_hit->get_pos();
    if(sq_hit->is_in_time()) h.setInTime();
    sraw->insertTriggerHit(h);
  }
  return true;
}
