/// OnlMonTrigNim.C
#include <iomanip>
#include <TStyle.h>
#include <TH2D.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <geom_svc/GeomSvc.h>
#include <UtilAna/UtilHist.h>
#include "OnlMonServer.h"
#include "OnlMonTrigNim.h"
using namespace std;

OnlMonTrigNim::OnlMonTrigNim()
{
  NumCanvases(1);
  Name("OnlMonTrigNim");
  Title("NIM Trigger");
}

int OnlMonTrigNim::InitOnlMon(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int OnlMonTrigNim::InitRunOnlMon(PHCompositeNode* topNode)
{
  h2_count = new TH2D("h2_count", "N of events in which each plane has hit;;", 8, 0.5, 8.5,  10, 0.5, 10.5);
  for (int ii = 1; ii <= 5; ii++) {
    ostringstream oss;
    oss << "NIM " << ii;
    h2_count->GetYaxis()->SetBinLabel(ii, oss.str().c_str());
    oss.str("");
    oss << "FPGA " << ii;
    h2_count->GetYaxis()->SetBinLabel(ii + 5, oss.str().c_str());
  }
  h2_count->GetXaxis()->SetBinLabel(1, "H1X");
  h2_count->GetXaxis()->SetBinLabel(2, "H1Y");
  h2_count->GetXaxis()->SetBinLabel(3, "H2X");
  h2_count->GetXaxis()->SetBinLabel(4, "H2Y");
  h2_count->GetXaxis()->SetBinLabel(5, "H3X");
  h2_count->GetXaxis()->SetBinLabel(6, "H4X");
  h2_count->GetXaxis()->SetBinLabel(7, "H4Y1");
  h2_count->GetXaxis()->SetBinLabel(8, "H4Y2");

  RegisterHist(h2_count);

  return Fun4AllReturnCodes::EVENT_OK;
}

int OnlMonTrigNim::ProcessEventOnlMon(PHCompositeNode* topNode)
{
  SQEvent*     event   = findNode::getClass<SQEvent    >(topNode, "SQEvent");
  SQHitVector* hit_vec = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if (!event || !hit_vec) return Fun4AllReturnCodes::ABORTEVENT;

  /// Sort out all hits.
  HitListMap_t lhm;
  for (SQHitVector::ConstIter it = hit_vec->begin(); it != hit_vec->end(); it++) {
    // We better select good hits here in the fugure
    lhm[(*it)->get_detector_id()].push_back(*it);
  }

  /// Judge if each plane should fire (i.e. has a hit) or not
  GeomSvc* geom = GeomSvc::instance();
  bool hit_h1x  = lhm.find(geom->getDetectorID("H1B"  )) != lhm.end()
               || lhm.find(geom->getDetectorID("H1T"  )) != lhm.end();
  bool hit_h1y  = lhm.find(geom->getDetectorID("H1L"  )) != lhm.end()
               || lhm.find(geom->getDetectorID("H1R"  )) != lhm.end();
  bool hit_h2x  = lhm.find(geom->getDetectorID("H2B"  )) != lhm.end()
               || lhm.find(geom->getDetectorID("H2T"  )) != lhm.end();
  bool hit_h2y  = lhm.find(geom->getDetectorID("H2L"  )) != lhm.end()
               || lhm.find(geom->getDetectorID("H2R"  )) != lhm.end();
  bool hit_h3x  = lhm.find(geom->getDetectorID("H3B"  )) != lhm.end()
               || lhm.find(geom->getDetectorID("H3T"  )) != lhm.end();
  bool hit_h4x  = lhm.find(geom->getDetectorID("H4B"  )) != lhm.end()
               || lhm.find(geom->getDetectorID("H4T"  )) != lhm.end();
  bool hit_h4y1 = lhm.find(geom->getDetectorID("H4Y1L")) != lhm.end()
               || lhm.find(geom->getDetectorID("H4Y1R")) != lhm.end();
  bool hit_h4y2 = lhm.find(geom->getDetectorID("H4Y2L")) != lhm.end()
               || lhm.find(geom->getDetectorID("H4Y2R")) != lhm.end();

  /// Fill the histogram per trigger
  bool trig_flag[10] = {
    event->get_trigger(SQEvent::NIM1),
    event->get_trigger(SQEvent::NIM2),
    event->get_trigger(SQEvent::NIM3),
    event->get_trigger(SQEvent::NIM4),
    event->get_trigger(SQEvent::NIM5),
    event->get_trigger(SQEvent::MATRIX1),
    event->get_trigger(SQEvent::MATRIX2),
    event->get_trigger(SQEvent::MATRIX3),
    event->get_trigger(SQEvent::MATRIX4),
    event->get_trigger(SQEvent::MATRIX5)
  };
  for (int i_trig = 0; i_trig < 10; i_trig++) {
    if (! trig_flag[i_trig]) continue;
    h2_count->Fill(0.0, i_trig+1); // Count N of all events in the underflow bin
    if (hit_h1x ) h2_count->Fill(1, i_trig+1);
    if (hit_h1y ) h2_count->Fill(2, i_trig+1);
    if (hit_h2x ) h2_count->Fill(3, i_trig+1);
    if (hit_h2y ) h2_count->Fill(4, i_trig+1);
    if (hit_h3x ) h2_count->Fill(5, i_trig+1);
    if (hit_h4x ) h2_count->Fill(6, i_trig+1);
    if (hit_h4y1) h2_count->Fill(7, i_trig+1);
    if (hit_h4y2) h2_count->Fill(8, i_trig+1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int OnlMonTrigNim::EndOnlMon(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int OnlMonTrigNim::FindAllMonHist()
{
  h2_count = (TH2*)FindMonHist("h2_count");
  return (h2_count  ?  0  :  1);
}

int OnlMonTrigNim::DrawMonitor()
{
  static TH2* h2_rate = 0;
  if (! h2_rate) delete h2_rate;
  h2_rate = (TH2*)h2_count->Clone("h2_rate");
  h2_rate->SetTitle("Percent of events in which each plane has hit");
  for (int iy = 1; iy <= h2_rate->GetNbinsY(); iy++) {
    double cont0 = h2_rate->GetBinContent(0, iy); // underflow bin
    if (cont0 <= 0) continue;
    for (int ix = 1; ix <= h2_rate->GetNbinsX(); ix++) {
      double cont = h2_rate->GetBinContent(ix, iy);
      h2_rate->SetBinContent(ix, iy, 100 * cont / cont0);
    }
  }
  h2_rate->GetZaxis()->SetRangeUser(0, 100);

  OnlMonCanvas* can0 = GetCanvas(0);
  TPad* pad0 = can0->GetMainPad();
  pad0->SetGrid();
  pad0->Divide(1, 2);
  pad0->cd(1);
  h2_count->Draw("colz");
  pad0->cd(2);
  h2_rate->Draw("colz");
  h2_rate->SetMarkerSize(2.0); // = text size (1.0 by default)
  gStyle->SetPaintTextFormat("3.0f");
  h2_rate->Draw("TEXTsame");

  if (h2_count->Integral() == 0) {
    can0->AddMessage("No event.");
    can0->SetStatus(OnlMonCanvas::WARN);
  } else {
    can0->SetStatus(OnlMonCanvas::OK);
  }

  return 0;
}
