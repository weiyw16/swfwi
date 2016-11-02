/*
 * essfwiframework.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: rice
 */


extern "C"
{
#include <rsf.h>
}

#include <time.h>

#include <cmath>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cstdlib>
#include <functional>
#include <vector>
#include <set>

#include "logger.h"
#include "common.h"
#include "ricker-wavelet.h"
#include "sum.h"
#include "sf-velocity-reader.h"
#include "shotdata-reader.h"
#include "random-code.h"
#include "encoder.h"
#include "velocity.h"
#include "sfutil.h"
#include "parabola-vertex.h"
#include "fwiframework.h"

#include "aux.h"

namespace {

void updateGrad(float *pre_gradient, const float *cur_gradient, float *update_direction,
                           int model_size, int iter) {
  if (iter == 0) {
    std::copy(cur_gradient, cur_gradient + model_size, update_direction);
    std::copy(cur_gradient, cur_gradient + model_size, pre_gradient);
  } else {
    float beta = 0.0f;
    float a = 0.0f;
    float b = 0.0f;
    float c = 0.0f;
    int   i = 0;
    for (i = 0; i < model_size; i ++) {
      a += (cur_gradient[i] * cur_gradient[i]);
      b += (cur_gradient[i] * pre_gradient[i]);
      c += (pre_gradient[i] * pre_gradient[i]);
    }

    beta = (a - b) / c;

    if (beta < 0.0f) {
      beta = 0.0f;
    }

    for (i = 0; i < model_size; i ++) {
      update_direction[i] = cur_gradient[i] + beta * update_direction[i];
    }

    TRACE() << "Save current gradient to pre_gradient for the next iteration's computation";
    std::copy(cur_gradient, cur_gradient + model_size, pre_gradient);
  }
}


void second_order_virtual_source_forth_accuracy(float *vsrc, int num) {
  float *tmp_vsrc = (float *)malloc(num * sizeof(float));
  memcpy(tmp_vsrc, vsrc, num * sizeof(float));
  int i = 0;
  for (i = 0; i < num; i ++) {
    if ( i <= 1) {
      vsrc[i] = 0.0f;
      continue;
    }

    if ( (num - 1) == i || (num - 2) == i) {
      vsrc[i] = 0.0f;
      continue;
    }

    vsrc[i] = -1. / 12 * tmp_vsrc[i - 2] + 4. / 3 * tmp_vsrc[i - 1] -
              2.5 * tmp_vsrc[i] + 4. / 3 * tmp_vsrc[i + 1] - 1. / 12 * tmp_vsrc[i + 2];
  }

  free(tmp_vsrc);
}

void transVsrc(std::vector<float> &vsrc, int nt, int ng) {
  std::vector<float> trans(nt * ng);
  matrix_transpose(&vsrc[0], &trans[0], ng, nt);
  for (int ig = 0; ig < ng; ig++) {
    second_order_virtual_source_forth_accuracy(&trans[ig * nt], nt);
  }

  matrix_transpose(&trans[0], &vsrc[0], nt, ng);
}

void cross_correlation(float *src_wave, float *vsrc_wave, float *image, int model_size, float scale) {
  for (int i = 0; i < model_size; i ++) {
    image[i] -= src_wave[i] * vsrc_wave[i] * scale;
  }

}

void calgradient(const Damp4t10d &fmMethod,
    const std::vector<float> &encSrc,
    const std::vector<float> &vsrc,
    std::vector<float> &g0,
    int nt, float dt,
		int shot_id)
{
  int nx = fmMethod.getnx();
  int nz = fmMethod.getnz();
  int ns = fmMethod.getns();
  int ng = fmMethod.getng();
  const ShotPosition &allGeoPos = fmMethod.getAllGeoPos();
  const ShotPosition &allSrcPos = fmMethod.getAllSrcPos();

  std::vector<float> bndr = fmMethod.initBndryVector(nt);
  std::vector<float> sp0(nz * nx, 0);
  std::vector<float> sp1(nz * nx, 0);
  std::vector<float> gp0(nz * nx, 0);
  std::vector<float> gp1(nz * nx, 0);


  ShotPosition curSrcPos = allSrcPos.clipRange(shot_id, shot_id);

  for(int it=0; it<nt; it++) {
    fmMethod.addSource(&sp1[0], &encSrc[it], curSrcPos);
    //printf("it = %d, forward 1\n", it);
    fmMethod.stepForward(&sp0[0], &sp1[0]);
    //printf("it = %d, forward 2\n", it);
    std::swap(sp1, sp0);
    //fmMethod.writeBndry(&bndr[0], &sp0[0], it); -test
		const int check_step = 5;
    if ((it > 0) && (it != (nt - 1)) && !(it % check_step)) {
      char check_file_name1[64];
      char check_file_name2[64];
      sprintf(check_file_name1, "./check_time_%d_1.su", it);
      sprintf(check_file_name2, "./check_time_%d_2.su", it);
			FILE *f1 = fopen(check_file_name1, "wb");
			FILE *f2 = fopen(check_file_name2, "wb");
			fwrite(&sp0[0], sizeof(float), nx * nz, f1);
			fwrite(&sp1[0], sizeof(float), nx * nz, f2);
			fclose(f1);
			fclose(f2);
    }
  }
	char check_file_name1[64] = "./check_time_last_1.su";
	char check_file_name2[64] = "./check_time_last_2.su";
	FILE *f1 = fopen(check_file_name1, "wb");
	FILE *f2 = fopen(check_file_name2, "wb");
	fwrite(&sp0[0], sizeof(float), nx * nz, f1);
	fwrite(&sp1[0], sizeof(float), nx * nz, f2);
	fclose(f1);
	fclose(f2);

  for(int it = nt - 1; it >= 0 ; it--) {
    //fmMethod.readBndry(&bndr[0], &sp0[0], it);	-test
		const int check_step = 5;
		if(it == nt - 1)
		{
			char check_file_name1[64] = "./check_time_last_1.su";
			char check_file_name2[64] = "./check_time_last_2.su";
			FILE *f1 = fopen(check_file_name1, "rb");
			FILE *f2 = fopen(check_file_name2, "rb");
			fread(&sp1[0], sizeof(float), nx * nz, f1);
			fread(&sp0[0], sizeof(float), nx * nz, f2);
			fclose(f1);
			fclose(f2);
		}
		else if ((check_step > 0) && !(it % check_step) && (it != 0)) {
			char check_file_name1[64];
			char check_file_name2[64];
			sprintf(check_file_name1, "./check_time_%d_1.su", it);
			sprintf(check_file_name2, "./check_time_%d_2.su", it);
			FILE *f1 = fopen(check_file_name1, "rb");
			FILE *f2 = fopen(check_file_name2, "rb");
			fread(&sp1[0], sizeof(float), nx * nz, f1);
			fread(&sp0[0], sizeof(float), nx * nz, f2);
			fclose(f1);
			fclose(f2);
		}

    //std::swap(sp0, sp1); -test
    fmMethod.stepBackward(&sp0[0], &sp1[0]);
    //fmMethod.subEncodedSource(&sp0[0], &encSrc[it]);
    std::swap(sp0, sp1);	//-test
    fmMethod.subSource(&sp0[0], &encSrc[it], curSrcPos);

    /**
     * forward propagate receviers
     */
    fmMethod.addSource(&gp1[0], &vsrc[it * ng], allGeoPos);
    //printf("it = %d, receiver 1\n", it);
    fmMethod.stepForward(&gp0[0], &gp1[0]);
    //printf("it = %d, receiver 2\n", it);
    std::swap(gp1, gp0);

    if (dt * it > 0.4) {
      //printf("it = %d, cross 1\n", it);
      cross_correlation(&sp0[0], &gp0[0], &g0[0], g0.size(), 1.0);
      //printf("it = %d, cross 2\n", it);
    } else if (dt * it > 0.3) {
      //printf("it = %d, cross 3\n", it);
      cross_correlation(&sp0[0], &gp0[0], &g0[0], g0.size(), (dt * it - 0.3) / 0.1);
      //printf("it = %d, cross 4\n", it);
    } else {
      //printf("it = %d, cross 5\n");
      break;
    }
 }
}

} /// end of namespace


FwiFramework::FwiFramework(Damp4t10d &method, const FwiUpdateSteplenOp &updateSteplenOp,
    const FwiUpdateVelOp &_updateVelOp,
    const std::vector<float> &_wlt, const std::vector<float> &_dobs) :
    fmMethod(method), updateStenlelOp(updateSteplenOp), updateVelOp(_updateVelOp), wlt(_wlt), dobs(_dobs),
    essRandomCodes(ESS_SEED),
    ns(method.getns()), ng(method.getng()), nt(method.getnt()),
    nx(method.getnx()), nz(method.getnz()), dx(method.getdx()), dt(method.getdt()),
    updateobj(0), initobj(0)
{
  g0.resize(nx*nz, 0);
  updateDirection.resize(nx*nz, 0);
}

/*
void FwiFramework::epoch(int iter) {
  // create random codes
  //const std::vector<int> encodes = essRandomCodes.genPlus1Minus1(ns);
  std::vector<int> encodes = essRandomCodes.genPlus1Minus1(ns);
	for(int i = 0 ; i < encodes.size() ; i ++)
		encodes[i] = 0;
	encodes[1] = 1;

  std::stringstream ss;
  std::copy(encodes.begin(), encodes.end(), std::ostream_iterator<int>(ss, " "));
  DEBUG() << "code is: " << ss.str();

  Encoder encoder(encodes);
  std::vector<float> encsrc  = encoder.encodeSource(wlt);
	for(int i = 0 ; i < wlt.size() ; i ++)
		if(encsrc[i*ns + 1] - wlt[i] >= 1e-6 || wlt[i] - encsrc[i*ns + 1] >= 1e-6)
			printf("src rss %d = %lf\n", i, encsrc[i*ns + 1] - wlt[i]);
  std::vector<float> encobs = encoder.encodeObsData(dobs, nt, ng);
	for(int i = 0 ; i < encobs.size() ; i ++)
		if(encobs[i] - dobs[encobs.size() * 1 + i] >= 1e-6 || dobs[encobs.size() * 1+ i ] - encobs[i] >= 1e-6)
		{
			printf("obj rss %d = %lf\n", i, encobs[i] - dobs[encobs.size() * 1 + i]);
			exit(1);
		}

  std::vector<float> dcal(nt * ng, 0);
  fmMethod.EssForwardModeling(encsrc, dcal);
  fmMethod.fwiRemoveDirectArrival(&encobs[0],1);
  fmMethod.fwiRemoveDirectArrival(&dcal[0],1);

	if(iter == 0)
	{
		std::vector<float> trans(ng * nt, 0);
    matrix_transpose(&dobs[encobs.size() * 1 + 0], &trans[0], ng, nt);
		FILE *f = fopen("shot_merge.bin","wb");
		fwrite(&trans[0], sizeof(float), nt * ng, f);
		fclose(f);
	}

  std::vector<float> vsrc(nt * ng, 0);
  vectorMinus(encobs, dcal, vsrc);
  float obj1 = cal_objective(&vsrc[0], vsrc.size());
  initobj = iter == 0 ? obj1 : initobj;
  DEBUG() << format("obj: %e") % obj1;

  //printf("check a\n");
  transVsrc(vsrc, nt, ng);
  //printf("check b\n");

  std::vector<float> g1(nx * nz, 0);
  //printf("check c\n");
  calgradient(fmMethod, encsrc, vsrc, g1, nt, dt);
  //printf("check d\n");

  DEBUG() << format("grad %.20f") % sum(g1);
  //printf("check e\n");

  fmMethod.scaleGradient(&g1[0]);
  fmMethod.maskGradient(&g1[0]);

  updateGrad(&g0[0], &g1[0], &updateDirection[0], g0.size(), iter);

  updateStenlelOp.bindEncSrcObs(encsrc, encobs);
  float steplen;
  updateStenlelOp.calsteplen(updateDirection, obj1, iter, steplen, updateobj);

  Velocity &exvel = fmMethod.getVelocity();
  updateVelOp.update(exvel, exvel, updateDirection, steplen);

  fmMethod.refillBoundary(&exvel.dat[0]);
}
*/

void FwiFramework::epoch(int iter) {
	std::vector<float> g2(nx * nz, 0);
	std::vector<float> encsrc(wlt);
	std::vector<float> encobs(ng * nt, 0);
	float obj1 = 0.0f;
	for(int is = 0 ; is < ns ; is ++) {
		INFO() << format("calculate gradient, shot id: %d") % is;
		memcpy(&encobs[0], &dobs[is * ng * nt], sizeof(float) * ng * nt);

		INFO() << "sum encobs: " << std::accumulate(encobs.begin(), encobs.begin() + ng * nt, 0.0f);
		INFO() << encsrc[0] << " " << encsrc[132];
		INFO() << "sum encsrc: " << std::accumulate(encsrc.begin(), encsrc.begin() + nt, 0.0f);

		std::vector<float> dcal(nt * ng, 0);
		fmMethod.FwiForwardModeling(encsrc, dcal, is);

		std::vector<float> trans_dcal(ng * nt, 0.0f);
		matrix_transpose(&dcal[0], &trans_dcal[0], ng, nt);

		/*
		if(iter == 1 && is == 0)
		{
			FILE *f_dcal = fopen("dcal.bin", "wb");
			fwrite(&trans_dcal[0], sizeof(float), ng * nt, f_dcal);
			fclose(f_dcal);
		}
		*/

		INFO() << dcal[0];
		INFO() << "sum dcal: " << std::accumulate(dcal.begin(), dcal.begin() + ng * nt, 0.0f);

		fmMethod.fwiRemoveDirectArrival(&encobs[0], is);
		fmMethod.fwiRemoveDirectArrival(&dcal[0], is);

		INFO() << "sum encobs2: " << std::accumulate(encobs.begin(), encobs.begin() + ng * nt, 0.0f);
		INFO() << "sum dcal2: " << std::accumulate(dcal.begin(), dcal.begin() + ng * nt, 0.0f);

		std::vector<float> vsrc(nt * ng, 0);
		vectorMinus(encobs, dcal, vsrc);
		obj1 += cal_objective(&vsrc[0], vsrc.size());
		initobj = iter == 0 ? obj1 : initobj;
		DEBUG() << format("obj: %e") % obj1;

		transVsrc(vsrc, nt, ng);

		INFO() << "sum vsrc: " << std::accumulate(vsrc.begin(), vsrc.begin() + ng * nt, 0.0f);

		std::vector<float> g1(nx * nz, 0);
		calgradient(fmMethod, encsrc, vsrc, g1, nt, dt, is);

		DEBUG() << format("grad %.20f") % sum(g1);

		//fmMethod.scaleGradient(&g1[0]);
		fmMethod.maskGradient(&g1[0]);

		/*
		char filename[20];
		sprintf(filename, "gradient%02d.bin", is);
		FILE *f = fopen(filename,"wb");
		fwrite(&g1[0], sizeof(float), nx * nz, f);
		fclose(f);
		*/

		std::transform(g2.begin(), g2.end(), g1.begin(), g2.begin(), std::plus<float>());

		DEBUG() << format("global grad %.20f") % sum(g2);
	}

  updateGrad(&g0[0], &g2[0], &updateDirection[0], g0.size(), iter);

	float steplen;
	float obj_val1 = 0, obj_val2 = 0, obj_val3 = 0;

	updateStenlelOp.calsteplen(dobs, updateDirection, obj1, iter, steplen, updateobj);
	
	float alpha1 = updateStenlelOp.alpha1;
	float alpha2 = updateStenlelOp.alpha2;
	float alpha3 = updateStenlelOp.alpha3;
	float obj_val1_sum = updateStenlelOp.obj_val1_sum;
	float obj_val2_sum = updateStenlelOp.obj_val2_sum;
	float obj_val3_sum = updateStenlelOp.obj_val3_sum;
	float maxAlpha3 = updateStenlelOp.maxAlpha3;
	bool	toParabolic = updateStenlelOp.toParabolic;
	updateStenlelOp.parabola_fit(alpha1, alpha2, alpha3, obj_val1_sum, obj_val2_sum, obj_val3_sum, maxAlpha3, toParabolic, iter, steplen, updateobj);

  INFO() << format("In calculate_steplen(): iter %d  steplen (alpha4) = %e") % iter % steplen;

  Velocity &exvel = fmMethod.getVelocity();
	INFO() << format("sum vel %f") % sum(exvel.dat);
  updateVelOp.update(exvel, exvel, updateDirection, steplen);
	INFO() << format("sum vel2 %f") % sum(exvel.dat);

  fmMethod.refillBoundary(&exvel.dat[0]);
}


void FwiFramework::writeVel(sf_file file) const {
  fmMethod.sfWriteVel(fmMethod.getVelocity().dat, file);
}

float FwiFramework::getUpdateObj() const {
  return updateobj;
}

float FwiFramework::getInitObj() const {
  return initobj;
}
