/*
 * spongabc4d.cpp
 *
 *  Created on: Feb 28, 2016
 *      Author: rice
 */

#include <ostream>
#include <cmath>
#include "spongabc4d.h"
#include "sum.h"
#include "logger.h"

static void expand(Velocity &exvel, const Velocity &v0, int nb) {
  int nx = v0.nx;
  int nz = v0.nz;
  int nxpad = nx + 2 * nb;
  int nzpad = nz + nb;
  const std::vector<float> &a = v0.dat;
  std::vector<float> &b = exvel.dat;

  /// internal

  DEBUG() << format("000: %.20f") % sum(b);
  DEBUG() << format("000: sum a %.20f") % sum(a);
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      b[(nb + ix) * nzpad + iz] = a[ix * nz + iz];
//      DEBUG() << format("ix %d, iz %d, nb %d, aidx %d, vv[%d] %.20f") %
//          ix % iz % nb % (ix * nz + iz) % ((nb + ix) * nzpad + iz) % b[(nb + ix) * nzpad + iz];
    }
  }

  DEBUG() << format("111: %.20f") % sum(b);

  /// boundary
  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < nb; iz++) {
      b[ix * nzpad + (nzpad - iz - 1)] = b[ix * nzpad + (nzpad - nb - 1)];  /* bottom*/
    }
  }

  for (int ix = 0; ix < nb; ix++) {
    for (int iz = 0; iz < nzpad; iz++) {
      b[ix * nzpad + iz] = b[nb * nzpad + iz];/* left */
      b[(nxpad - ix - 1) * nzpad + iz] = b[(nxpad - nb - 1) * nzpad + iz]; /* right */
    }
  }
}

static void numericTrans(Velocity &v0, float dt) {
  int nx = v0.nx;
  int nz = v0.nz;
  std::vector<float> &vv = v0.dat;

  for(int ix=0;ix<nx;ix++){
    for(int iz=0;iz<nz;iz++){
      float tmp=vv[ix * nz + iz]*dt;
      vv[ix * nz + iz]=tmp*tmp;/* vv=vv^2*dt^2 */
    }
  }

}

SpongAbc4d::SpongAbc4d(float _dt, float _dx, float _dz, int _nb)
  : IModeling(), bndr(_nb), dt(_dt), dx(_dx), dz(_dz), nb(_nb)
{
  initCoeff();
  initbndr();
}

void SpongAbc4d::stepForward(float *p0, float *p1) const {
  int nx = vel->nx;
  int nz = vel->nz;
  const std::vector<float> &vv = vel->dat;

  for (int ix = 2; ix < nx - 2; ix++) {
    for (int iz = 2; iz < nz - 2; iz++) {
      float tmp = c0 * p1[ix * nz + iz] +
                  c11 * (p1[ix * nz + (iz - 1)] + p1[ix * nz + (iz + 1)]) +
                  c12 * (p1[ix * nz + (iz - 2)] + p1[ix * nz + (iz + 2)]) +
                  c21 * (p1[(ix - 1) * nz + iz] + p1[(ix + 1) * nz + iz]) +
                  c22 * (p1[(ix - 2) * nz + iz] + p1[(ix + 2) * nz + iz]);

      p0[ix * nz + iz] = 2 * p1[ix * nz + iz] - p0[ix * nz + iz] + vv[ix * nz + iz] * tmp;
    }
  }

  applySponge(p0);
  applySponge(p1);

}

void SpongAbc4d::initCoeff() {
  /*< initialize 4-th order FD coefficients >*/
  float tmp = 1.0/(dz*dz);
  c11 = 4.0*tmp/3.0;
  c12= -tmp/12.0;
  tmp = 1.0/(dx*dx);
  c21 = 4.0*tmp/3.0;
  c22= -tmp/12.0;
  c0=-2.0*(c11+c12+c21+c22);
}

void SpongAbc4d::applySponge(float* p) const {
  int nz = vel->nz;
  int nx = vel->nx;

  for(int ib=0; ib<nb; ib++) {
    float w = bndr[ib];

    int ibz = nz-ib-1;
    for(int ix=0; ix<nx; ix++) {
      p[ix * nz + ibz] *= w; /* bottom sponge */
    }

    int ibx = nx-ib-1;
    for(int iz=0; iz<nz; iz++) {
      p[ib  * nz + iz] *= w; /*   left sponge */
      p[ibx * nz + iz] *= w; /*  right sponge */
    }
  }
}

void SpongAbc4d::initbndr() {
  for(int ib=0;ib<nb;ib++){
    float tmp=0.015*(nb-ib);
    bndr[ib]=std::exp(-tmp*tmp);
  }
}



Velocity SpongAbc4d::transformVelocityForModeling(const Velocity& v0) const {
  int nxpad = v0.nx + 2 * nb;
  int nzpad = v0.nz + nb; // free surface
  Velocity exvel(nxpad, nzpad);

  expand(exvel, v0, nb);
  DEBUG() << format("before trans sum exvel: %.20f") % sum(exvel.dat);

  numericTrans(exvel, dt);

  return exvel;
}


void SpongAbc4d::addSource(float* p, const float* source, int ns,
    const int* sxz, int snz)
{
  int nzpad = snz + nb;

  for (int is = 0; is < ns; is++) {
    int sx = sxz[is] / snz + nb;
    int sz = sxz[is] % snz;
    p[sx * nzpad + sz] += source[is];
  }
}

void SpongAbc4d::addSource(float* p, const float* source, const ShotPosition& pos) {
  int nzpad = vel->nz;

  for (int is = 0; is < pos.ns; is++) {
    int sx = pos.getx(is) + nb;
    int sz = pos.getz(is);
    p[sx * nzpad + sz] += source[is];
//    DEBUG() << format("sx %d, sz %d, source[%d] %f") % sx % sz % is % source[is];
  }
}

void SpongAbc4d::recordSeis(float* seis_it, const float* p,
    const ShotPosition& geoPos) {

  int ng = geoPos.ns;
  int nzpad = vel->nz;

//  DEBUG() << format("ng %d") % ng;
//  float sum = 0;
  for (int ig = 0; ig < ng; ig++) {
    int gx = geoPos.getx(ig) + nb;
    int gz = geoPos.getz(ig);
    int idx = gx * nzpad + gz;
    seis_it[ig] = p[idx];
//    DEBUG() << format("ig %d, idx %d, v %.20f") % ig % idx % seis_it[ig];
  }

}

std::vector<float> SpongAbc4d::initBndryVector(int nt) const {
}