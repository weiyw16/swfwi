/*
 * docproduct.cpp
 * for making doc product experiments for fwi framework
 *
 *  Created on: Sep 25, 2017
 *      Author: weiyw
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

FwiFrameDoc::FwiFrameDoc(ForwardModeling &method, const FwiUpdateSteplenOp &updateSteplenOp,
    const FwiUpdateVelOp &_updateVelOp,
    const std::vector<float> &_wlt, const std::vector<float> &_dobs) :
    FwiBase(method, _wlt, _dobs), updateStenlelOp(updateSteplenOp), updateVelOp(_updateVelOp)
{
}

void FwiFrameDoc::epoch(int iter) {

	std::vector<float> g1(nx * nz, 0);
	std::vector<float> g2(nx * nz, 0);
	std::vector<float> encobs(ng * nt, 0);
	
	/* source parallel in mpi */
	int rank, np, k, ntask, shot_begin, shot_end;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	k = std::ceil(ns * 1.0 / np);
	ntask = std::min(k, ns - rank*k);
	shot_begin = rank * k;
	shot_end = shot_begin + ntask;

	float local_obj1 = 0.0f, obj1 = 0.0f;

	for(int is = shot_begin ; is < shot_end ; is ++) {
		std::vector<float> encobs_trans(nt * ng, 0.0f);
		INFO() << format("calculate gradient, shot id: %d") % is;
		memcpy(&encobs_trans[0], &dobs[is * ng * nt], sizeof(float) * ng * nt);/* copy dobs --> encobs_trans */

		matrix_transpose(&encobs_trans[0], &encobs[0], ng, nt);

		/*
			 if(iter == 1)
			 {
			 sf_file sf_encobs = sf_output("encobs.rsf");
			 sf_putint(sf_encobs, "n1", nt);
			 sf_putint(sf_encobs, "n2", ng);
			 sf_floatwrite(&encobs[0], nt * ng, sf_encobs);
			 }
			 */

		/*
			 sf_file sf_wlt = sf_output("wlt.rsf");
			 sf_putint(sf_wlt, "n1", nt);
			 sf_floatwrite(&wlt[0], nt, sf_wlt);
			 */

		INFO() << "sum encobs: " << std::accumulate(encobs.begin(), encobs.end(), 0.0f);
		INFO() << wlt[0] << " " << wlt[132];
		INFO() << "sum wlt: " << std::accumulate(wlt.begin(), wlt.begin() + nt, 0.0f);// why should add them all ? 

		std::vector<float> dcal(nt * ng, 0);
		std::vector<float> dcal_trans(ng * nt, 0.0f);
		fmMethod.FwiForwardModeling(wlt, dcal_trans, is);// where get fmMethod ??
		matrix_transpose(&dcal_trans[0], &dcal[0], ng, nt);//common.cpp


		/*
			 if(iter == 0)
			 {
			 char fg2[64];
			 sprintf(fg2, "dcal_%02d.rsf", is);
			 sf_file sf_dcal = sf_output(fg2);
			 sf_putint(sf_dcal, "n1", nt);
			 sf_putint(sf_dcal, "n2", ng);
			 sf_floatwrite(&dcal[0], nt * ng, sf_dcal);
			 }
			 */

		/*
			 if(iter == 1 && is == 0)
			 {
			 FILE *f_dcal = fopen("dcal.bin", "wb");
			 fwrite(&trans_dcal[0], sizeof(float), ng * nt, f_dcal);
			 fclose(f_dcal);
			 }
			 */

		INFO() << dcal[0];
		INFO() << "sum dcal: " << std::accumulate(dcal.begin(), dcal.end(), 0.0f);

		fmMethod.fwiRemoveDirectArrival(&encobs[0], is);
		fmMethod.fwiRemoveDirectArrival(&dcal[0], is);

		/*
			 sf_file sf_encobs = sf_output("encobs2.rsf");
			 sf_putint(sf_encobs, "n1", nt);
			 sf_putint(sf_encobs, "n2", ng);
			 sf_floatwrite(&encobs[0], nt * ng, sf_encobs);
			 exit(1);
			 */

		/*
			 if(iter == 1)
			 {
			 sf_file sf_dcal = sf_output("dcal2.rsf");
			 sf_putint(sf_dcal, "n1", nt);
			 sf_putint(sf_dcal, "n2", ng);
			 sf_floatwrite(&dcal[0], nt * ng, sf_dcal);
			 exit(1);
			 }
			 */

		INFO() << "sum encobs2: " << std::accumulate(encobs.begin(), encobs.end(), 0.0f);
		INFO() << "sum dcal2: " << std::accumulate(dcal.begin(), dcal.end(), 0.0f);

		std::vector<float> vsrc(nt * ng, 0);
		vectorMinus(encobs, dcal, vsrc);// common.h
		local_obj1 += cal_objective(&vsrc[0], vsrc.size());
		initobj = iter == 0 ? local_obj1 : initobj;
		//DEBUG() << format("obj: %e") % obj1;
		INFO() << "obj: " << local_obj1 << "\n";

		transVsrc(vsrc, nt, ng);

		INFO() << "sum vsrc: " << std::accumulate(vsrc.begin(), vsrc.end(), 0.0f);

		g1.assign(nx * nz, 0.0f);
		//std::vector<float> g1(nx * nz, 0);
		calgradient(fmMethod, wlt, vsrc, g1, nt, dt, is, rank);

		/*
			 sf_file sf_vsrc= sf_output("vsrc.rsf");
			 sf_putint(sf_vsrc, "n1", nt);
			 sf_putint(sf_vsrc, "n2", ng);
			 sf_floatwrite(&vsrc[0], nt * ng, sf_vsrc);
			 exit(1);
			 */

		DEBUG() << format("grad %.20f") % sum(g1);

		//fmMethod.scaleGradient(&g1[0]);
		fmMethod.maskGradient(&g1[0]);

		/*
			 char fg1[64];
			 sprintf(fg1, "g1_%02d.rsf", is);
			 sf_file sf_g1 = sf_output(fg1);
			 sf_putint(sf_g1, "n1", nz);
			 sf_putint(sf_g1, "n2", nx);
			 sf_floatwrite(&g1[0], nx * nz, sf_g1);
			 */

		/*
			 char filename[20];
			 sprintf(filename, "gradient%02d.bin", is);
			 FILE *f = fopen(filename,"wb");
			 fwrite(&g1[0], sizeof(float), nx * nz, f);
			 fclose(f);
			 */

		std::transform(g2.begin(), g2.end(), g1.begin(), g2.begin(), std::plus<float>());//g2 + g1 ??

		/*
			 sf_file sf_g2 = sf_output("g2.rsf");
			 sf_putint(sf_g2, "n1", nz);
			 sf_putint(sf_g2, "n2", nx);
			 sf_floatwrite(&g2[0], nx * nz, sf_g2);
			 exit(1);
			 */

		DEBUG() << format("global grad %.20f") % sum(g2);
	}

	g1.assign(nx * nz, 0.0f);
	MPI_Allreduce(&g2[0], &g1[0], g2.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&local_obj1, &obj1, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	/* int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
	                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) */
	if(rank == 0)
	{
		DEBUG() << format("****** global grad %.20f") % sum(g1);
		DEBUG() << format("****** sum obj: %.20f") % obj1;
	}


	/*
		 if(rank == 0 && iter == 0)
		 {
		 sf_file sf_g2 = sf_output("g2.rsf");
		 sf_putint(sf_g2, "n1", nz);
		 sf_putint(sf_g2, "n2", nx);
		 sf_floatwrite(&g1[0], nx * nz, sf_g2);
		 }
		 */

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

void FwiFrameDoc::calgradient(const ForwardModeling &fmMethod,
    const std::vector<float> &wlt,
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
    fmMethod.addSource(&sp1[0], &wlt[it], curSrcPos);
    //printf("it = %d, forward 1\n", it);
    fmMethod.stepForward(sp0,sp1);
    //printf("it = %d, forward 2\n", it);
    std::swap(sp1, sp0);
    fmMethod.writeBndry(&bndr[0], &sp0[0], it); //-test
		/*
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
		*/
  }
	/*
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
	*/

  std::vector<float> vsrc_trans(ng * nt, 0.0f);
  matrix_transpose(const_cast<float*>(&vsrc[0]), &vsrc_trans[0], nt, ng);

  for(int it = nt - 1; it >= 0 ; it--) {
    fmMethod.readBndry(&bndr[0], &sp0[0], it);	//-test
		/*
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
		*/

    std::swap(sp0, sp1); //-test
    fmMethod.stepBackward(&sp0[0], &sp1[0]);
    //fmMethod.subEncodedSource(&sp0[0], &wlt[it]);
    //std::swap(sp0, sp1);	//-test
    fmMethod.subSource(&sp0[0], &wlt[it], curSrcPos);

    /**
     * forward propagate receviers
     */
    fmMethod.addSource(&gp1[0], &vsrc_trans[it * ng], allGeoPos);
    //printf("it = %d, receiver 1\n", it);
    fmMethod.stepForward(gp0,gp1);
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

