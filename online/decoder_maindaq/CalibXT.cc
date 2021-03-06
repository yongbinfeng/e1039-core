#include <iomanip>
#include <TGraphErrors.h>
#include <interface_main/SQParamDeco.h>
#include <interface_main/SQRun.h>
#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <geom_svc/CalibParamXT.h>
#include <geom_svc/CalibParamInTimeTaiwan.h>
#include "CalibXT.h"
using namespace std;

CalibXT::CalibXT(const std::string& name) : SubsysReco(name), m_cal_xt(0), m_cal_int(0)
{
  ;
}

CalibXT::~CalibXT()
{
  if (m_cal_xt ) delete m_cal_xt ;
  if (m_cal_int) delete m_cal_int;
}

int CalibXT::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int CalibXT::InitRun(PHCompositeNode* topNode)
{
  SQParamDeco* param_deco = findNode::getClass<SQParamDeco>(topNode, "SQParamDeco");
  SQRun*       run_header = findNode::getClass<SQRun      >(topNode, "SQRun");
  if (!param_deco || !run_header) return Fun4AllReturnCodes::ABORTEVENT;

  if (! m_cal_xt) m_cal_xt = new CalibParamXT();
  m_cal_xt->SetMapIDbyDB(run_header->get_run_id());
  m_cal_xt->ReadFromDB();
  param_deco->set_variable(m_cal_xt->GetParamID(), m_cal_xt->GetMapID());

  if (! m_cal_int) m_cal_int = new CalibParamInTimeTaiwan();
  m_cal_int->SetMapIDbyDB(run_header->get_run_id());
  m_cal_int->ReadFromDB();
  param_deco->set_variable(m_cal_int->GetParamID(), m_cal_int->GetMapID());

  return Fun4AllReturnCodes::EVENT_OK;
}

int CalibXT::process_event(PHCompositeNode* topNode)
{
  SQHitVector* hit_vec = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if (! hit_vec) return Fun4AllReturnCodes::ABORTEVENT;

  for (SQHitVector::Iter it = hit_vec->begin(); it != hit_vec->end(); it++) {
    SQHit* hit = *it;
    TGraphErrors* gr_t2x;
    TGraphErrors* gr_t2dx;
    int det = hit->get_detector_id();
    if (m_cal_xt->Find(det, gr_t2x, gr_t2dx)) {
      int ele = hit->get_element_id();
      double center, width;
      if (! m_cal_int->Find(det, ele, center, width)) {
        cerr << "  WARNING:  Cannot find the in-time parameter for det=" << det << " ele=" << ele << " in CalibXT.\n";
        continue;
        //return Fun4AllReturnCodes::ABORTEVENT;
      }
      float t0 = center + width / 2;
      float drift_time = t0 - hit->get_tdc_time();
      hit->set_drift_distance(gr_t2x->Eval(drift_time));
      /// No field for resolution in SQHit now.
      //cout << "check: " << det << " " << ele << " " << t0 << " " << drift_time << " " << hit->get_drift_distance() << endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CalibXT::End(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
