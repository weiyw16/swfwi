/*
 * damp4t10d.h
 *
 *  Created on: Feb 29, 2016
 *      Author: rice
 */

#ifndef SRC_FM2D_DAMP4T10D_H_
#define SRC_FM2D_DAMP4T10D_H_

extern "C" {
#include <rsf.h>
#include "fdutil.h"
}
#include <boost/function.hpp>
#include "velocity.h"
#include "shot-position.h"

class Damp4t10d {
public:
  Damp4t10d(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float dt, float dx, float fm, int nb, int nt);

  Velocity expandDomain(const Velocity &vel);
  Velocity expandDomain_notrans(const Velocity &vel);


	void addBornwv(float *fullwv_t0, float *fullwv_t1, float *fullwv_t2, const float *exvel_m, float dt, int it, float *rp1) const;
  void stepForward(float *p0, float *p1) const;
  void stepBackward(float *p0, float *p1) const;
  void bindVelocity(const Velocity &_vel);
  void bindRealVelocity(const Velocity &_vel);
  void addSource(float *p, const float *source, const ShotPosition &pos) const;
  void addSource(float *p, const float *source, int is) const;
  void subSource(float *p, const float *source, const ShotPosition &pos) const;
  void addEncodedSource(float *p, const float *encsrc) const;
  void recordSeis(float *seis_it, const float *p) const;
  void bornMaskGradient(float *grad, int H) const;
  void bornScaleGradient(float *grad, int H) const;
  void maskGradient(float *grad) const;
  void scaleGradient(float *grad) const;
  void refillBoundary(float *vel) const;
  void sfWriteVel(const std::vector<float> &exvel, sf_file file) const;

  void fwiRemoveDirectArrival(float* data, int shot_id) const;
  void removeDirectArrival(float* data) const;
  void subEncodedSource(float *p, const float *source) const;
  void refillVelStencilBndry();

  std::vector<float> initBndryVector(int nt) const;
  void writeBndry(float* _bndr, const float* p, int it) const;
  void readBndry(const float* _bndr, float* p, int it) const;

  void FwiForwardModeling(const std::vector<float> &encsrc, std::vector<float> &dcal, int shot_id) const;
  void EssForwardModeling(const std::vector<float> &encsrc, std::vector<float> &dcal) const;
	void BornForwardModeling(const std::vector<float>& exvel, const std::vector<float>& encSrc, std::vector<float>& dcal, int shot_id) const;

public:
  const Velocity &getVelocity() const;
  Velocity &getVelocity();
	const std::vector<float> getVelocityDiff() const;
  const ShotPosition &getAllSrcPos() const;
  const ShotPosition &getAllGeoPos() const;
  int getns() const;
  int getng() const;
  float getdt() const;
  float getdx() const;
  int getnt() const;
  int getnx() const;
  int getnz() const;

private:
  void manipSource(float *p, const float *source, const ShotPosition &pos, boost::function2<float, float, float> op) const;
  void recordSeis(float *seis_it, const float *p, const ShotPosition &geoPos) const;
  void removeDirectArrival(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float* data, int nt, float t_width) const;

public:
	void GetXBoundaryMPos(int xPos, int zPos, int *xMPos, int *zMPos, int nx, int nz) const;
	void GetZBoundaryMPos(int xPos, int zPos, int *xMPos, int *zMPos, int nx, int nz) const;
	void initCPML(int nx, int nz);
	void applyCPML(float *uLa, float *u, float *uNe, const float *vel, int nx, int nz);
	void initFdUtil(sf_file &vinit, Velocity *v, int nb, float dx, float dt);

private:
  const static int EXFDBNDRYLEN = 6;

private:
  const Velocity *vel;
  const Velocity *vel_real;
  const ShotPosition *allSrcPos;
  const ShotPosition *allGeoPos;
  float dt;
  float dx;
  float fm;
  int bx0, bxn;
  int bz0, bzn;
  int nt;
  mutable int bndrSize;
  mutable int bndrWidth;


private:
  std::vector<float> bndr;

	std::vector<float> psiX, psiXLa, phiX, phiXLa, EtaX, EtaXLa, psi2X, psi2XLa, phi2X, phi2XLa, u020BXLa;
	std::vector<float> psiZ, psiZLa, phiZ, phiZLa, EtaZ, EtaZLa, psi2Z, psi2ZLa, phi2Z, phi2ZLa, u002BZLa;
	std::vector<float> ux, uz, uxLa, uzLa;
	int psixlen, psizlen;

	struct fdm2 *fd;
	struct spon *sp;
	struct abc2 *abc;

};

#endif /* SRC_FM2D_DAMP4T10D_H_ */
