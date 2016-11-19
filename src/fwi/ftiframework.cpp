/*
 * ftiframework.cpp
 *
 *  Created on: Mar 10, 2016
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
#include "ftiframework.h"

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

void cross_correlation_born(float *src_wave, float *vsrc_wave, float *image, int nx, int nz, float scale, int H) {
	float t_src_wave, t_vsrc_wave;
	for(int h = -H ; h < H ; h ++) {
		for (int i = 0; i < nx ; i ++) {
			for (int j = 0; j < nz ; j ++) {
				t_src_wave = i + h < nx ? src_wave[(i + h) * nz + j] : src_wave[(nx - 1) * nz + j];
				t_src_wave = i + h >= 0 ? t_src_wave : src_wave[0 * nz + j];
				t_vsrc_wave = i - h < nx ? vsrc_wave[(i - h) * nz + j] : vsrc_wave[(nx - 1) * nz + j];
				t_vsrc_wave = i - h >= 0 ? t_vsrc_wave: vsrc_wave[0 * nz + j];
				image[(h + H) * nx * nz + i * nz + j] += t_src_wave * t_vsrc_wave * scale;
			}
		}
	}
}

void image_born(const Damp4t10d &fmMethod,
    const std::vector<float> &encSrc,
    const std::vector<float> &vsrc,
    std::vector<float> &g0,
    int nt, float dt,
		int shot_id, int rank, int H)
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
    fmMethod.stepForward(&sp0[0], &sp1[0]);
    std::swap(sp1, sp0);
    fmMethod.writeBndry(&bndr[0], &sp0[0], it); //-test
  }

  std::vector<float> vsrc_trans(ng * nt, 0.0f);
  matrix_transpose(const_cast<float*>(&vsrc[0]), &vsrc_trans[0], nt, ng);

  for(int it = nt - 1; it >= 0 ; it--) {
    fmMethod.readBndry(&bndr[0], &sp0[0], it);	//-test
    std::swap(sp0, sp1); //-test
    fmMethod.stepBackward(&sp0[0], &sp1[0]);
    //fmMethod.subEncodedSource(&sp0[0], &encSrc[it]);
    fmMethod.subSource(&sp0[0], &encSrc[it], curSrcPos);

    /**
     * forward propagate receviers
     */
    fmMethod.addSource(&gp1[0], &vsrc_trans[it * ng], allGeoPos);
    fmMethod.stepForward(&gp0[0], &gp1[0]);
    std::swap(gp1, gp0);

    cross_correlation_born(&sp0[0], &gp0[0], &g0[0], nx, nz, 1.0, H);
 }
}

void calgradient(const Damp4t10d &fmMethod,
    const std::vector<float> &encSrc,
    const std::vector<float> &vsrc,
    std::vector<float> &I,
    std::vector<float> &gd,
    int nt, float dt,
		int shot_id, int rank, int H)
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

	std::vector<float> ps(nt * nx * nz, 0);
	std::vector<float> pg(nt * nx * nz, 0);

  ShotPosition curSrcPos = allSrcPos.clipRange(shot_id, shot_id);

  for(int it=0; it<nt; it++) {
    fmMethod.addSource(&sp1[0], &encSrc[it], curSrcPos);
    fmMethod.stepForward(&sp0[0], &sp1[0]);
    std::swap(sp1, sp0);
		for(int ix = 0 ; ix < nx ; ix ++)
			for(int iz = 0 ; iz < nz ; iz ++)
				ps[it * nx * nz + ix * nz + iz] = sp1[ix * nz + iz] - sp0[ix * nz + iz];
  }

  std::vector<float> vsrc_trans(ng * nt, 0.0f);
  matrix_transpose(const_cast<float*>(&vsrc[0]), &vsrc_trans[0], nt, ng);

  for(int it = nt - 1; it >= 0 ; it--) {
    //fmMethod.readBndry(&bndr[0], &sp0[0], it);	-test
    //std::swap(sp0, sp1); -test
    fmMethod.stepBackward(&sp0[0], &sp1[0]);
    //fmMethod.subEncodedSource(&sp0[0], &encSrc[it]);
    std::swap(sp0, sp1);	//-test
    fmMethod.subSource(&sp0[0], &encSrc[it], curSrcPos);

    /**
     * forward propagate receviers
     */
    fmMethod.addSource(&gp1[0], &vsrc_trans[it * ng], allGeoPos);
    fmMethod.stepForward(&gp0[0], &gp1[0]);
    std::swap(gp1, gp0);
		for(int ix = 0 ; ix < nx ; ix ++)
			for(int iz = 0 ; iz < nz ; iz ++)
				pg[it * nx * nz + ix * nz + iz] = gp1[ix * nz + iz] - gp0[ix * nz + iz];
	}

	sp0.assign(nx * nz, 0);
	sp1.assign(nx * nz, 0);
	gp0.assign(nx * nz, 0);
	gp1.assign(nx * nz, 0);

	const Velocity &exvel = fmMethod.getVelocity();

	for(int it=0; it<nt; it++) {
		for(int h = 0 ; h < H ; h ++)
			for(int ix = 0 ; ix < nx ; ix ++) 
				for(int iz = 0 ; iz < nz ; iz ++) 
					sp1[ix * nz + iz] += ps[it * nx * nz + (ix + 2 * h) * nz + iz] * h * h * I[h * nx * nz + (ix + h) * nz + iz];
		fmMethod.stepForward(&sp0[0], &sp1[0]);
		std::swap(sp1, sp0);
		for(int ix = 0 ; ix < nx ; ix ++) 
			for(int iz = 0 ; iz < nz ; iz ++) 
				gd[ix * nz + iz] += 2 * sp0[ix * nz + iz] * pg[it * nx * nz + ix * nz + iz] * exvel.dat[ix * nz + iz];
	}

	for(int it = nt - 1; it >= 0 ; it--) {
		for(int h = 0 ; h < H ; h ++)
			for(int ix = 0 ; ix < nx ; ix ++)
				for(int iz = 0 ; iz < nz ; iz ++)
					gp1[ix * nz + iz] += pg[(nt - (it - 1)) * nx * nz + (ix - 2 * h) * nz + iz] * h * h * I[h * nx * nz + (ix - h) * nz + iz];
    fmMethod.stepForward(&gp0[0], &gp1[0]);
    std::swap(gp1, gp0);
		for(int ix = 0 ; ix < nx ; ix ++) 
			for(int iz = 0 ; iz < nz ; iz ++) 
				gd[ix * nz + iz] += 2 * gp0[ix * nz + iz] * ps[it * nx * nz + ix * nz + iz] * exvel.dat[ix * nz + iz];
	}


}

} /// end of namespace


FtiFramework::FtiFramework(Damp4t10d &method, const FwiUpdateSteplenOp &updateSteplenOp,
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

void FtiFramework::epoch(int iter) {
	std::vector<float> encsrc(wlt);
	std::vector<float> encobs(ng * nt, 0);
	int rank, np, k, ntask, shot_begin, shot_end;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	k = std::ceil(ns * 1.0 / np);
	ntask = std::min(k, ns - rank*k);
	shot_begin = rank * k;
	shot_end = shot_begin + ntask;
	float local_obj1 = 0.0f, obj1 = 0.0f;
	int H = 10;
	std::vector<float> g1(2 * H * nx * nz, 0);
	std::vector<float> g2(2 * H * nx * nz, 0);

	sf_file sf_g2;
	if(rank == 0 && iter == 0)
	{
		sf_g2 = sf_output("g2.rsf");
		sf_putint(sf_g2, "n1", nz);
		sf_putint(sf_g2, "n2", nx);
		sf_putint(sf_g2, "n3", 2 * H);
	}

	for(int is = shot_begin ; is < shot_end ; is ++) {
		std::vector<float> encobs_trans(nt * ng, 0.0f);
		INFO() << format("calculate gradient, shot id: %d") % is;
		memcpy(&encobs_trans[0], &dobs[is * ng * nt], sizeof(float) * ng * nt);

		matrix_transpose(&encobs_trans[0], &encobs[0], ng, nt);

		INFO() << "sum encobs: " << std::accumulate(encobs.begin(), encobs.end(), 0.0f);
		INFO() << encsrc[0] << " " << encsrc[132];
		INFO() << "sum encsrc: " << std::accumulate(encsrc.begin(), encsrc.begin() + nt, 0.0f);

		std::vector<float> dcal(nt * ng, 0);
		std::vector<float> dcal_trans(ng * nt, 0.0f);
		fmMethod.FwiForwardModeling(encsrc, dcal_trans, is);
		matrix_transpose(&dcal_trans[0], &dcal[0], ng, nt);

		INFO() << dcal[0];
		INFO() << "sum dcal: " << std::accumulate(dcal.begin(), dcal.end(), 0.0f);

		fmMethod.fwiRemoveDirectArrival(&encobs[0], is);
		fmMethod.fwiRemoveDirectArrival(&dcal[0], is);

		INFO() << "sum encobs2: " << std::accumulate(encobs.begin(), encobs.end(), 0.0f);
		INFO() << "sum dcal2: " << std::accumulate(dcal.begin(), dcal.end(), 0.0f);

		std::vector<float> vsrc(nt * ng, 0);
		vectorMinus(encobs, dcal, vsrc);
		local_obj1 += cal_objective(&vsrc[0], vsrc.size());
		initobj = iter == 0 ? local_obj1 : initobj;
		//DEBUG() << format("obj: %e") % obj1;
		INFO() << "obj: " << local_obj1 << "\n";

		//transVsrc(vsrc, nt, ng);

		std::copy(encobs.begin(), encobs.end(), vsrc.begin());	//-test
		//std::copy(dcal.begin(), dcal.end(), vsrc.begin());	//-test

		INFO() << "sum vsrc: " << std::accumulate(vsrc.begin(), vsrc.end(), 0.0f);

		g1.assign(nx * nz, 0.0f);
		//std::vector<float> g1(nx * nz, 0);
		image_born(fmMethod, encsrc, vsrc, g1, nt, dt, is, rank, H);

		DEBUG() << format("grad %.20f") % sum(g1);

		//fmMethod.scaleGradient(&g1[0]);
		fmMethod.maskGradient(&g1[0]);

		std::transform(g2.begin(), g2.end(), g1.begin(), g2.begin(), std::plus<float>());

		DEBUG() << format("global grad %.20f") % sum(g2);
	}

	g1.assign(2 * H * nx * nz, 0.0f);
	MPI_Allreduce(&g2[0], &g1[0], g2.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&local_obj1, &obj1, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	if(rank == 0)
	{
		DEBUG() << format("****** global grad %.20f") % sum(g1);
		DEBUG() << format("****** sum obj: %.20f") % obj1;
	}

	if(rank == 0 && iter == 0)
	{
		sf_floatwrite(&g1[0], 2 * H * nx * nz, sf_g2);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	exit(1);

	std::vector<float> gd(nx * nz, 0);
	std::vector<float> grad(nx * nz, 0);
	for(int is = shot_begin ; is < shot_end ; is ++) {
		INFO() << format("Calculating gradient %d:") % is;
		std::vector<float> encobs_trans(nt * ng, 0.0f);
		std::vector<float> vsrc(ng * nt, 0);
		memcpy(&encobs_trans[0], &dobs[is * ng * nt], sizeof(float) * ng * nt);
		matrix_transpose(&encobs_trans[0], &encobs[0], ng, nt);
		std::copy(encobs.begin(), encobs.end(), vsrc.begin());	//-test
		calgradient(fmMethod, encsrc, vsrc, g1, gd, nt, dt, is, rank, H);
	}
	MPI_Allreduce(&gd[0], &grad[0], gd.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	if(rank == 0 && iter == 0) {
		 sf_file sf_g2 = sf_output("grad.rsf");
		 sf_putint(sf_g2, "n1", nz);
		 sf_putint(sf_g2, "n2", nx);
		 sf_floatwrite(&grad[0], nx * nz, sf_g2);
	}

	exit(1);

	updateGrad(&g0[0], &g1[0], &updateDirection[0], g0.size(), iter);

	/*
		 if(rank == 0 && iter == 0)
		 {
		 sf_file sf_ud = sf_output("updateDirection.rsf");
		 sf_putint(sf_ud, "n1", nz);
		 sf_putint(sf_ud, "n2", nx);
		 sf_floatwrite(&updateDirection[0], nx * nz, sf_ud);
		 }
		 MPI_Barrier(MPI_COMM_WORLD);
		 exit(1);
		 */

	float steplen;
	float obj_val1 = 0, obj_val2 = 0, obj_val3 = 0;

	updateStenlelOp.calsteplen(dobs, updateDirection, obj1, iter, steplen, updateobj, rank, shot_begin, shot_end);


	float alpha1 = updateStenlelOp.alpha1;
	float alpha2 = updateStenlelOp.alpha2;
	float alpha3 = updateStenlelOp.alpha3;
	float obj_val1_sum = updateStenlelOp.obj_val1_sum;
	float obj_val2_sum = updateStenlelOp.obj_val2_sum;
	float obj_val3_sum = updateStenlelOp.obj_val3_sum;
	float maxAlpha3 = updateStenlelOp.maxAlpha3;
	bool	toParabolic = updateStenlelOp.toParabolic;
	updateStenlelOp.parabola_fit(alpha1, alpha2, alpha3, obj_val1_sum, obj_val2_sum, obj_val3_sum, maxAlpha3, toParabolic, iter, steplen, updateobj);

	if(rank == 0)
	{
		INFO() << format("In calculate_steplen(): iter %d  steplen (alpha4) = %e") % iter % steplen;

		INFO() << format("steplen = %.20f") % steplen;
	}

	Velocity &exvel = fmMethod.getVelocity();

	/*
		 if(iter == 1)
		 {
		 sf_file sf_exvel = sf_output("exvel_before.rsf");
		 sf_putint(sf_exvel, "n1", nz);
		 sf_putint(sf_exvel, "n2", nx);
		 sf_floatwrite(&exvel.dat[0], nx * nz, sf_exvel);
		 exit(1);
		 }
		 */

	if(rank == 0)
		INFO() << format("sum vel %f") % sum(exvel.dat);

	updateVelOp.update(exvel, exvel, updateDirection, steplen);

	if(rank == 0)
		INFO() << format("sum vel2 %f") % sum(exvel.dat);

	/*
	if(rank == 0)
	{
		char f_name[64];
		sprintf(f_name, "exvel_after%02d.rsf", iter);
		sf_file sf_exvel2 = sf_output(f_name);
		sf_putint(sf_exvel2, "n1", nz);
		sf_putint(sf_exvel2, "n2", nx);
		sf_floatwrite(&exvel.dat[0], nx * nz, sf_exvel2);
		if(iter == 3)
		exit(1);
	}
	*/

	//fmMethod.refillBoundary(&exvel.dat[0]);

}


void FtiFramework::writeVel(sf_file file) const {
	fmMethod.sfWriteVel(fmMethod.getVelocity().dat, file);
}

float FtiFramework::getUpdateObj() const {
	return updateobj;
}

float FtiFramework::getInitObj() const {
	return initobj;
}
