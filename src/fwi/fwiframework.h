/*
 * essfwiframework.h
 *
 *  Created on: Nov 22, 2016
 *      Author: cbw
 */

#ifndef SRC_FWI2D_FWIFRAMEWORK_H_
#define SRC_FWI2D_FWIFRAMEWORK_H_

#include "damp4t10d.h"
#include "fwibase.h"
#include "fwiupdatevelop.h"
#include "fwiupdatesteplenop.h"
#include "random-code.h"

class FwiFramework : public FwiBase {
public:
  FwiFramework(Damp4t10d &fmMethod, const FwiUpdateSteplenOp &updateSteplenOp,
                  const FwiUpdateVelOp &updateVelOp, const std::vector<float> &wlt,
                  const std::vector<float> &dobs);
	void epoch(int iter);

protected:
  FwiUpdateSteplenOp updateStenlelOp;
  const FwiUpdateVelOp &updateVelOp;
};

#endif /* SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_ */
