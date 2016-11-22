/*
 * ftiframework.h
 *
 *  Created on: Nov 22, 2016
 *      Author: cbw
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
	void calgradient(const Damp4t10d &fmMethod,
    const std::vector<float> &encSrc,
    const std::vector<float> &vsrc,
    std::vector<float> &I,
    std::vector<float> &gd,
    int nt, float dt,
		int shot_id, int rank, int H);
	void image_born(const Damp4t10d &fmMethod, const std::vector<float> &encSrc,
    const std::vector<float> &vsrc,
    std::vector<float> &g0,
    int nt, float dt,
		int shot_id, int rank, int H);
	void cross_correlation(float *src_wave, float *vsrc_wave, float *image, int nx, int nz, float scale, int H);

};

#endif /* SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_ */
