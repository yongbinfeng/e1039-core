/*
 * SQMCHit_v1.h
 *
 *  Created on: Oct 29, 2017
 *      Author: yuhw@nmsu.edu
 */

#ifndef _H_SQMCHit_v1_H_
#define _H_SQMCHit_v1_H_


#include <phool/PHObject.h>
#include <iostream>

#include "SQHit_v1.h"

class SQMCHit_v1 : public SQHit_v1 {

public:

  SQMCHit_v1();
  virtual ~SQMCHit_v1() {}

  // PHObject virtual overloads

  void         identify(std::ostream& os = std::cout) const;
  void         Reset() {*this = SQMCHit_v1();}
  int          isValid() const;
  SQHit*        Clone() const {return (new SQMCHit_v1(*this));}

  virtual int          get_track_id() const                             {return _track_id;}
  virtual void         set_track_id(const int a)                        {_track_id = a;}

private:

  int _track_id;  ///< unique identifier within container

  ClassDef(SQMCHit_v1, 1);
};


#endif /* _H_SQMCHit_v1_H_ */
