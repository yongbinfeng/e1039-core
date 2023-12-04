#ifndef __CALIB_HIT_ELEMENT_POS_H__
#define __CALIB_HIT_ELEMENT_POS_H__
#include <fun4all/SubsysReco.h>

class SQHitVector;

class SRawEvent;

/// SubsysReco module to set the position of SQHit using GeomSvc.
class CalibHitElementPos: public SubsysReco {
  SQHitVector* m_vec_hit;
  SQHitVector* m_vec_trhit;
  SRawEvent* _rawEvent;
  
 public:
  enum INPUT_TYPE  {E906, E1039};
  
  CalibHitElementPos(const std::string &name = "CalibHitElementPos");
  virtual ~CalibHitElementPos();
  
  void setInputTy(CalibHitElementPos::INPUT_TYPE input_ty) { _input_type = input_ty; }
  
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 protected:
  CalibHitElementPos::INPUT_TYPE  _input_type;
};

#endif // __CALIB_HIT_ELEMENT_POS_H__
