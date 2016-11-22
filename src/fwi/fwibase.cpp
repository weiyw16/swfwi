/*
 * fwibase.cpp
 *
 *  Created on: Nov 22, 2016
 *      Author: cbw
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
#include "fwibase.h"

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
  for (int ig = 0; ig < ng; ig++) {
    second_order_virtual_source_forth_accuracy(&vsrc[ig * nt], nt);
  }
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
		int shot_id, int rank)
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
      sprintf(check_file_name1, "./rank_%d_check_time_%d_1.su", rank, it);
      sprintf(check_file_name2, "./rank_%d_check_time_%d_2.su", rank, it);
			FILE *f1 = fopen(check_file_name1, "wb");
			FILE *f2 = fopen(check_file_name2, "wb");
			fwrite(&sp0[0], sizeof(float), nx * nz, f1);
			fwrite(&sp1[0], sizeof(float), nx * nz, f2);
			fclose(f1);
			fclose(f2);
    }
  }
	char check_file_name1[64];
	char check_file_name2[64];
	sprintf(check_file_name1, "./rank_%d_check_time_last_1.su", rank);
	sprintf(check_file_name2, "./rank_%d_check_time_last_2.su", rank);
	FILE *f1 = fopen(check_file_name1, "wb");
	FILE *f2 = fopen(check_file_name2, "wb");
	fwrite(&sp0[0], sizeof(float), nx * nz, f1);
	fwrite(&sp1[0], sizeof(float), nx * nz, f2);
	fclose(f1);
	fclose(f2);

  std::vector<float> vsrc_trans(ng * nt, 0.0f);
  matrix_transpose(const_cast<float*>(&vsrc[0]), &vsrc_trans[0], nt, ng);

  for(int it = nt - 1; it >= 0 ; it--) {
    //fmMethod.readBndry(&bndr[0], &sp0[0], it);	-test
		const int check_step = 5;
		if(it == nt - 1)
		{
			char check_file_name1[64];
			char check_file_name2[64];
			sprintf(check_file_name1, "./rank_%d_check_time_last_1.su", rank);
			sprintf(check_file_name2, "./rank_%d_check_time_last_2.su", rank);
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
			sprintf(check_file_name1, "./rank_%d_check_time_%d_1.su", rank, it);
			sprintf(check_file_name2, "./rank_%d_check_time_%d_2.su", rank, it);
			FILE *f1 = fopen(check_file_name1, "rb");
			FILE *f2 = fopen(check_file_name2, "rb");
			fread(&sp1[0], sizeof(float), nx * nz, f1);
			fread(&sp0[0], sizeof(float), nx * nz, f2);
			fclose(f1);
			fclose(f2);
		}

    fmMethod.stepBackward(&sp0[0], &sp1[0]);
    std::swap(sp0, sp1);	//-test
    fmMethod.subSource(&sp0[0], &encSrc[it], curSrcPos);

    /**
     * forward propagate receviers
     */
    fmMethod.addSource(&gp1[0], &vsrc_trans[it * ng], allGeoPos);
    fmMethod.stepForward(&gp0[0], &gp1[0]);
    std::swap(gp1, gp0);

    if (dt * it > 0.4) {
      cross_correlation(&sp0[0], &gp0[0], &g0[0], g0.size(), 1.0);
    } else if (dt * it > 0.3) {
      cross_correlation(&sp0[0], &gp0[0], &g0[0], g0.size(), (dt * it - 0.3) / 0.1);
    } else {
      break;
    }
 }
}

} /// end of namespace


FwiBase::FwiBase(Damp4t10d &method, const std::vector<float> &_wlt, const std::vector<float> &_dobs) :
    fmMethod(method), wlt(_wlt), dobs(_dobs),
    ns(method.getns()), ng(method.getng()), nt(method.getnt()),
    nx(method.getnx()), nz(method.getnz()), dx(method.getdx()), dt(method.getdt()),
    updateobj(0), initobj(0)
{
  g0.resize(nx*nz, 0);
  updateDirection.resize(nx*nz, 0);
}

void FwiBase::writeVel(sf_file file) const {
	fmMethod.sfWriteVel(fmMethod.getVelocity().dat, file);
}

float FwiBase::getUpdateObj() const {
	return updateobj;
}

float FwiBase::getInitObj() const {
	return initobj;
}
