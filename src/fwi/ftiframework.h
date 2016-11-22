/*
 * essfwiframework.h
 *
 *  Created on: Mar 10, 2016
 *      Author: rice
 */

#ifndef SRC_FWI2D_FTIFRAMEWORK_H_
#define SRC_FWI2D_FTIFRAMEWORK_H_

#include "damp4t10d.h"
#include "fwiupdatevelop.h"
#include "fwiupdatesteplenop.h"
#include "random-code.h"
#include "fwiframework.h"

class FtiFramework: public FwiFramework {
public:
  FtiFramework(Damp4t10d &fmMethod, const FwiUpdateSteplenOp &updateSteplenOp,
                  const FwiUpdateVelOp &updateVelOp, const std::vector<float> &wlt,
                  const std::vector<float> &dobs);
  void epoch(int iter);

};

#endif /* SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_ */
