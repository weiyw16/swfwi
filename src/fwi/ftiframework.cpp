/*
 * ftiframework.cpp
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
#include "ftiframework.h"

FtiFramework::FtiFramework(Damp4t10d &method, const FwiUpdateSteplenOp &updateSteplenOp,
    const FwiUpdateVelOp &_updateVelOp,
    const std::vector<float> &_wlt, const std::vector<float> &_dobs) :
    FwiFramework(method, updateSteplenOp, _updateVelOp, _wlt, _dobs)
{
}

void FtiFramework::epoch(int iter) {
	std::vector<float> encobs(ng * nt, 0);
	int rank, np, k, ntask, shot_begin, shot_end;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	k = std::ceil(ns * 1.0 / np);
	ntask = std::min(k, ns - rank*k);
	shot_begin = rank * k;
	shot_end = shot_begin + ntask;
	float local_obj1 = 0.0f, obj1 = 0.0f;
	int H = 60;
	std::vector<float> img(2 * H * nx * nz, 0);
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
		INFO() << format("calculate image, shot id: %d") % is;
		memcpy(&encobs_trans[0], &dobs[is * ng * nt], sizeof(float) * ng * nt);
		matrix_transpose(&encobs_trans[0], &encobs[0], ng, nt);
		fmMethod.fwiRemoveDirectArrival(&encobs[0], is);
		img.assign(2 * H * nx * nz, 0.0f);
		image_born(fmMethod, wlt, encobs, img, nt, dt, is, rank, H);
		DEBUG() << ("sum grad: ") << std::accumulate(&img[H * nx * nz], &img[(H + 1) * nx * nz], 0.0f);
		//fmMethod.bornScaleGradient(&img[0]);
		fmMethod.bornMaskGradient(&img[0], H);
		std::transform(g2.begin(), g2.end(), img.begin(), g2.begin(), std::plus<float>());
		DEBUG() << ("global grad: ") << std::accumulate(&g2[H * nx * nz], &g2[(H + 1) * nx * nz], 0.0f);
	}

	img.assign(2 * H * nx * nz, 0.0f);
	MPI_Allreduce(&g2[0], &img[0], g2.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&local_obj1, &obj1, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	if(rank == 0)
	{
		DEBUG() << ("****** global grad: ") << std::accumulate(&img[H * nx * nz], &img[(H + 1) * nx * nz], 0.0f);
		DEBUG() << format("****** sum obj: %.20f") % obj1;
	}

	if(rank == 0 && iter == 0)
	{
		sf_floatwrite(&img[0], 2 * H * nx * nz, sf_g2);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	exit(1);

	if(rank == 0 && iter == 0) {
		std::vector<float> img0(2 * H * nx * nz, 0.0f);
		sf_file sf_img0 = sf_input("g2.rsf");
		sf_floatread(&img0[0], 2 * H * nx * nz, sf_img0);

		sf_file sf_img1 = sf_output("g2_backup.rsf");
		sf_putint(sf_img1, "n1", nz);
		sf_putint(sf_img1, "n2", nx);
		sf_putint(sf_img1, "n3", 2 * H);
		sf_floatwrite(&img0[0], 2 * H * nx * nz, sf_img1);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	exit(1);

	std::vector<float> gd(nx * nz, 0);
	std::vector<float> grad(nx * nz, 0);
	for(int is = shot_begin ; is < shot_end ; is ++) {
		INFO() << format("************Calculating gradient %d:") % is;
		std::vector<float> encobs_trans(nt * ng, 0.0f);
		std::vector<float> vsrc(ng * nt, 0);
		memcpy(&encobs_trans[0], &dobs[is * ng * nt], sizeof(float) * ng * nt);
		matrix_transpose(&encobs_trans[0], &encobs[0], ng, nt);	//removeDirectArrival?
		calgradient(fmMethod, wlt, encobs, img, gd, nt, dt, is, rank, H);
	}
	MPI_Allreduce(&gd[0], &grad[0], gd.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	if(rank == 0 && iter == 0) {
		 sf_file sf_g2 = sf_output("grad.rsf");
		 sf_putint(sf_g2, "n1", nz);
		 sf_putint(sf_g2, "n2", nx);
		 sf_floatwrite(&grad[0], nx * nz, sf_g2);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	exit(1);

	updateGrad(&g0[0], &img[0], &updateDirection[0], g0.size(), iter);

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

void FtiFramework::calgradient(const Damp4t10d &fmMethod,
    const std::vector<float> &wlt,
    const std::vector<float> &vsrc,
    std::vector<float> &img,
    std::vector<float> &gd,
    int nt, float dt,
		int shot_id, int rank, int H)
{
  int nx = fmMethod.getnx();
  int nz = fmMethod.getnz();
  int ns = fmMethod.getns();
  int ng = fmMethod.getng();
	std::vector<float> gd0(nx * nz, 0.0f);
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
    fmMethod.addSource(&sp1[0], &wlt[it], curSrcPos);
    fmMethod.stepForward(&sp0[0], &sp1[0]);
    std::swap(sp1, sp0);
		for(int ix = 0 ; ix < nx ; ix ++)	//lack of the last timestep?
			for(int iz = 0 ; iz < nz ; iz ++)
				ps[it * nx * nz + ix * nz + iz] = sp1[ix * nz + iz] - sp0[ix * nz + iz];
  }

  std::vector<float> vsrc_trans(ng * nt, 0.0f);
  matrix_transpose(const_cast<float*>(&vsrc[0]), &vsrc_trans[0], nt, ng);

  for(int it = nt - 1; it >= 0 ; it--) {
		/*
    //fmMethod.readBndry(&bndr[0], &sp0[0], it);	-test
    //std::swap(sp0, sp1); -test
    fmMethod.stepBackward(&sp0[0], &sp1[0]);
    //fmMethod.subEncodedSource(&sp0[0], &wlt[it]);
    std::swap(sp0, sp1);	//-test
    fmMethod.subSource(&sp0[0], &wlt[it], curSrcPos);
		*/

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

	const Velocity &exvel = fmMethod.getVelocity();

	sp0.assign(nx * nz, 0);
	sp1.assign(nx * nz, 0);

	float ps_t = 0, pg_t = 0, img_t = 0;
	for(int it=0; it<nt; it++) {
		for(int h = -H ; h < H ; h ++) {
			int ind = h + H;
			for(int ix = 0 ; ix < nx ; ix ++) 
				for(int iz = 0 ; iz < nz ; iz ++) { 
				ps_t = ix + 2 * h >= 0 ? ps[it * nx * nz + (ix + 2 * h) * nz + iz] : 0;
				ps_t = ix + 2 * h < nx ? ps_t : 0;
				img_t = ix + h >= 0 ? img[ind * nx * nz + (ix + h) * nz + iz] : 0; 
				img_t = ix + h < nx ? img_t : 0; 
				sp1[ix * nz + iz] += ps_t * h * h * img_t;
			}
		}
		fmMethod.stepForward(&sp0[0], &sp1[0]);
		std::swap(sp1, sp0);
		for(int ix = 0 ; ix < nx ; ix ++) 
			for(int iz = 0 ; iz < nz ; iz ++) 
				gd0[ix * nz + iz] += 2 * sp0[ix * nz + iz] * pg[it * nx * nz + ix * nz + iz] * exvel.dat[ix * nz + iz];
	}

	gp0.assign(nx * nz, 0);
	gp1.assign(nx * nz, 0);

	for(int it = nt - 1; it >= 0 ; it--) {
		for(int h = -H ; h < H ; h ++) {
			int ind = h + H;
			for(int ix = 0 ; ix < nx ; ix ++) 
					for(int iz = 0 ; iz < nz ; iz ++) {
					pg_t = ix - 2 * h >= 0 ? pg[it * nx * nz + (ix - 2 * h) * nz + iz] : 0;
					pg_t = ix - 2 * h < nx ? pg_t : 0;
					img_t = ix - h >= 0 ? img[ind * nx * nz + (ix - h) * nz + iz] : 0;
					img_t = ix - h < nx ? img_t : 0;
					gp1[ix * nz + iz] +=  pg_t * h * h * img_t;
			}
		}
    fmMethod.stepForward(&gp0[0], &gp1[0]);
    std::swap(gp1, gp0);
		for(int ix = 0 ; ix < nx ; ix ++) 
			for(int iz = 0 ; iz < nz ; iz ++) 
				gd0[ix * nz + iz] += 2 * gp0[ix * nz + iz] * ps[it * nx * nz + ix * nz + iz] * exvel.dat[ix * nz + iz];
	}
	std::transform(gd.begin(), gd.end(), gd0.begin(), gd.begin(), std::plus<float>());
}

void FtiFramework::image_born(const Damp4t10d &fmMethod,
    const std::vector<float> &wlt,
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
    fmMethod.addSource(&sp1[0], &wlt[it], curSrcPos);
    fmMethod.stepForward(&sp0[0], &sp1[0]);
    std::swap(sp1, sp0);
    //fmMethod.writeBndry(&bndr[0], &sp0[0], it); //-test
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
    //fmMethod.readBndry(&bndr[0], &sp0[0], it);	//-test
    //std::swap(sp0, sp1); //-test
    fmMethod.stepBackward(&sp0[0], &sp1[0]);
    //fmMethod.subEncodedSource(&sp0[0], &wlt[it]);
    std::swap(sp0, sp1); //-test
    fmMethod.subSource(&sp0[0], &wlt[it], curSrcPos);

    /**
     * forward propagate receviers
     */
    fmMethod.addSource(&gp1[0], &vsrc_trans[it * ng], allGeoPos);
    fmMethod.stepForward(&gp0[0], &gp1[0]);
    std::swap(gp1, gp0);

    if (dt * it > 0.4) {
      //printf("it = %d, cross 1\n", it);
      cross_correlation(&sp0[0], &gp0[0], &g0[0], nx, nz, 1.0, H);
      //printf("it = %d, cross 2\n", it);
    } else if (dt * it > 0.3) {
      //printf("it = %d, cross 3\n", it);
      cross_correlation(&sp0[0], &gp0[0], &g0[0], nx, nz, (dt * it - 0.3) / 0.1, H);
      //printf("it = %d, cross 4\n", it);
    } else {
      //printf("it = %d, cross 5\n");
      break;
    }
    //cross_correlation_born(&sp0[0], &gp0[0], &g0[0], nx, nz, 1.0, H);
 }
}

void FtiFramework::cross_correlation(float *src_wave, float *vsrc_wave, float *image, int nx, int nz, float scale, int H) {
	float t_src_wave, t_vsrc_wave;
#pragma omp parallel for private(t_src_wave, t_vsrc_wave)
	for(int h = -H ; h < H ; h ++) {
		int ind = h + H;
		for (int i = 0; i < nx ; i ++) {
			for (int j = 0; j < nz ; j ++) {
				t_src_wave = i + h < nx ? src_wave[(i + h) * nz + j] : 0;
				t_src_wave = i + h >= 0 ? t_src_wave : 0;
				t_vsrc_wave = i - h < nx ? vsrc_wave[(i - h) * nz + j] : 0;
				t_vsrc_wave = i - h >= 0 ? t_vsrc_wave: 0;
				image[ind * nx * nz + i * nz + j] += t_src_wave * t_vsrc_wave * scale;
			}
		}
	}
}

