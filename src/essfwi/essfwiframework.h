/*
 * essfwiframework.h
 *
 *  Created on: Mar 10, 2016
 *      Author: rice
 */

#ifndef SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_
#define SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_

#include "damp4t10d.h"
#include "fwibase.h"
#include "updatevelop.h"
#include "updatesteplenop.h"
#include "random-code.h"

class EssFwiFramework : public FwiBase {
public:
  EssFwiFramework(Damp4t10d &fmMethod, const UpdateSteplenOp &updateSteplenOp,
                  const UpdateVelOp &updateVelOp, const std::vector<float> &wlt,
                  const std::vector<float> &dobs);

  void epoch(int iter, float lambdaX = 0, float lambdaZ = 0);

private:
  static const int ESS_SEED = 1;

private:
  UpdateSteplenOp updateStenlelOp;
  const UpdateVelOp &updateVelOp;
  RandomCodes essRandomCodes;
};

#endif /* SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_ */
