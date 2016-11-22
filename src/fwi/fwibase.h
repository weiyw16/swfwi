/*
 * fwibase.h
 *
 *  Created on: Nov 22, 2016
 *      Author: cbw
 */

#ifndef SRC_FWI2D_FWIBASE_H_
#define SRC_FWI2D_FWIBASE_H_

#include "damp4t10d.h"

class FwiBase {
public:
  FwiBase(Damp4t10d &fmMethod, const std::vector<float> &wlt,
                  const std::vector<float> &dobs);

  void writeVel(sf_file file) const;
  float getUpdateObj() const;
  float getInitObj() const;

protected:
  Damp4t10d &fmMethod;
  const std::vector<float> &wlt;  /// wavelet
  const std::vector<float> &dobs; /// actual observed data (nt*ng*ns)

protected: /// propagate from other construction
  int ns;
  int ng;
  int nt;
  int nx;
  int nz;
  float dx;
  float dt;

protected:
  std::vector<float> g0;               /// gradient in previous step
  std::vector<float> updateDirection;
  float updateobj;
  float initobj;
	float obj_val4;
};

#endif /* SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_ */
