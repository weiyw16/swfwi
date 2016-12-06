/*
 * damp4t10d.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: rice
 */

#include <cmath>
#include <functional>
#include "damp4t10d.h"
#include "logger.h"
#include "sum.h"
#include "sfutil.h"
#include "common.h"

extern "C" {
#include <rsf.h>
#include "fdutil.h"
#include "fd4t10s-damp-zjh.h"
#include "fd4t10s-zjh.h"
#include "fd4t10s-sponge.h"
}

#define SIGN(a,b) (b >= 0 ? fabs(a) : -fabs(a))
//#define FREE

static void initbndr(std::vector<float> &bndr, int nb) {
  for(int ib=0;ib<nb;ib++){
    float tmp=0.0087*(nb-ib-1);
    bndr[ib]=std::exp(-tmp*tmp);
  }
}

static void applySponge(float* p, const float *vel, const float *bndr, int nx, int nz, int nb, float dt, float dx) {
  for(int ib=0; ib<nb; ib++) {
    float w = bndr[ib];

    int ibz = nz-ib-1;
    for(int ix=0; ix<nx; ix++) {
#ifdef FREE
#else
      p[ix * nz + ib] *= w;
#endif
      p[ix * nz + ibz] *= w;
    }

    int ibx = nx-ib-1;
    for(int iz=0; iz<nz; iz++) {
      p[ib  * nz + iz] *= w;
      p[ibx * nz + iz] *= w;
    }
  }
}

void Damp4t10d::initFdUtil(sf_file &vinit, Velocity *v, int nb, float dx, float dt) {
	const int ompchunk = 8;
	fd = fdutil_init(false, false, sf_iaxa(vinit, 1), sf_iaxa(vinit, 2), bx0, ompchunk); 
	sp = sponge_make(bx0);
	float *v2 = (float *)malloc(sizeof(float) * v->nx * v->nz);
	float **vv;
	for(int i = 0 ; i < v->nx * v->nz ; i ++)
		v2[i] = 1.0 / v->dat[i] * dx * dx / dt / dt;
	vv = f1dto2d(v2, v->nx, v->nz);
	abc = abcone2d_make(EXFDBNDRYLEN, dt, vv, false, fd);
}

void Damp4t10d::GetXBoundaryMPos(int xPos, int zPos, int *xMPos, int *zMPos, int nx, int nz) const 
{
	int d = EXFDBNDRYLEN;
	int nb = bx0 - d;
	int cnx = nx - bx0 - bxn;
	if(xPos - d < nb)
		*xMPos = xPos - d;
	else
		*xMPos = xPos - d - cnx;
	*zMPos = zPos - d;
}

void Damp4t10d::GetZBoundaryMPos(int xPos, int zPos, int *xMPos, int *zMPos, int nx, int nz) const
{
	int d = EXFDBNDRYLEN;
	int nb = bx0 - d;
	int cnz = nz - bz0 - bzn;
	*xMPos = xPos - d;
	if(zPos - d < nb)
		*zMPos = zPos - d;
	else
		*zMPos = zPos - d - cnz;
}

void Damp4t10d::initCPML(int nx, int nz) {
	int d = EXFDBNDRYLEN;
	int nb = bx0 - d;

	psiX.resize(2 * nb * (nz - d * 2));
	psiXLa.resize(2 * nb * (nz - d * 2));
	phiX.resize(2 * nb * (nz - d * 2));
	phiXLa.resize(2 * nb * (nz - d * 2));
	EtaX.resize(2 * nb * (nz - d * 2));
	EtaXLa.resize(2 * nb * (nz - d * 2));
	psi2X.resize(2 * nb * (nz - d * 2));
	psi2XLa.resize(2 * nb * (nz - d * 2));
	phi2X.resize(2 * nb * (nz - d * 2));
	phi2XLa.resize(2 * nb * (nz - d * 2));
	u020BXLa.resize(2 * nb * (nz - d * 2));

	psiX.assign(2 * nb * (nz - d * 2),0);
	psiXLa.assign(2 * nb * (nz - d * 2),0);
	phiX.assign(2 * nb * (nz - d * 2),0);
	phiXLa.assign(2 * nb * (nz - d * 2),0);
	EtaX.assign(2 * nb * (nz - d * 2),0);
	EtaXLa.assign(2 * nb * (nz - d * 2),0);
	psi2X.assign(2 * nb * (nz - d * 2),0);
	psi2XLa.assign(2 * nb * (nz - d * 2),0);
	phi2X.assign(2 * nb * (nz - d * 2),0);
	phi2XLa.assign(2 * nb * (nz - d * 2),0);
	u020BXLa.assign(2 * nb * (nz - d * 2),0);

	psiZ.resize(2 * nb * (nx - d * 2));
	psiZLa.resize(2 * nb * (nx - d * 2));
	phiZ.resize(2 * nb * (nx - d * 2));
	phiZLa.resize(2 * nb * (nx - d * 2));
	EtaZ.resize(2 * nb * (nx - d * 2));
	EtaZLa.resize(2 * nb * (nx - d * 2));
	psi2Z.resize(2 * nb * (nx - d * 2));
	psi2ZLa.resize(2 * nb * (nx - d * 2));
	phi2Z.resize(2 * nb * (nx - d * 2));
	phi2ZLa.resize(2 * nb * (nx - d * 2));
	u002BZLa.resize(2 * nb * (nx - d * 2));

	psiZ.assign(2 * nb * (nx - d * 2),0);
	psiZLa.assign(2 * nb * (nx - d * 2),0);
	phiZ.assign(2 * nb * (nx - d * 2),0);
	phiZLa.assign(2 * nb * (nx - d * 2),0);
	EtaZ.assign(2 * nb * (nx - d * 2),0);
	EtaZLa.assign(2 * nb * (nx - d * 2),0);
	psi2Z.assign(2 * nb * (nx - d * 2),0);
	psi2ZLa.assign(2 * nb * (nx - d * 2),0);
	phi2Z.assign(2 * nb * (nx - d * 2),0);
	phi2ZLa.assign(2 * nb * (nx - d * 2),0);
	u002BZLa.assign(2 * nb * (nx - d * 2),0);

	ux.resize(nx * nz);
	uxLa.resize(nx * nz);
	uz.resize(nx * nz);
	uzLa.resize(nx * nz);

	ux.assign(nx * nz,0);
	uxLa.assign(nx * nz,0);
	uz.assign(nx * nz,0);
	uzLa.assign(nx * nz,0);

	psixlen = nz - d * 2;
	psizlen = 2 * nb;
}

void Damp4t10d::applyCPML(float *uLa, float *u, float *uNe, const float *vel, int nx, int nz) {
	int d = EXFDBNDRYLEN;
	int nBMPosX, nBMPosZ;
	float lB, dDlB, dD2lB, alphaDlB, alphaD2lB;
	float DlB0; //Larger ValuedxB leads to stronger attenation
	float alphaDlB0; //Larger ValueAlphaxB leads to faster phase variation
	int cnx = nx - bx0 - bxn;
	int cnz = nz - bz0 - bzn;
	float blengthx = cnx * dx;
	float blengthz = cnz * dx;
	float u020 = 0.0f, u002 = 0.0f, aB, bB;
#ifdef USE_OPENMP
	#pragma omp parallel for private(nBMPosX, nBMPosZ, lB, dDlB, dD2lB, alphaDlB, alphaD2lB, DlB0, alphaDlB0, u020, u002, aB, bB)
#endif
	for(int ix = d ; ix < nx - d ; ix ++) { 
		for(int iz = d ; iz < nz - d ; iz ++) {
			ux[ix * nz + iz] = (2. / 3. * (u[(ix + 1) * nz + iz] - u[(ix - 1) * nz + iz]) - 1. / 12. * (u[(ix + 2) * nz + iz] - u[(ix - 2) * nz + iz])) / dx;
			uz[ix * nz + iz] = (2. / 3. * (u[ix * nz + iz + 1] - u[ix * nz + iz - 1]) - 1. / 12. * (u[ix * nz + iz + 2] - u[ix * nz + iz - 2])) / dx;

			if (ix < bx0 || ix >= nx - bxn || iz < bz0 || iz >= nz - bzn)
			{
				u020 = 1.0 / 12.0 * (-30 * u[ix * nz + iz] + 16 *(u[(ix - 1) * nz + iz]+u[(ix + 1) * nz + iz]) - (u[(ix - 2) * nz + iz]+u[(ix + 2) * nz + iz])) / dx / dx;	
				u002 = 1.0 / 12.0 * (-30 * u[ix * nz + iz] + 16 *(u[ix * nz + iz - 1]+u[ix * nz + iz + 1]) - (u[ix * nz + iz - 2]+u[ix * nz + iz + 2])) / dx / dx;

				//X boundaries
				if(ix < bx0 || ix >= nx - bxn)
				{
					DlB0 = 16000; //Larger ValuedxB leads to stronger attenation
					alphaDlB0 = 0; //Larger ValueAlphaxB leads to faster phase variation
					//Generate boundary coordinates
					GetXBoundaryMPos(ix, iz, &nBMPosX, &nBMPosZ, nx, nz);
					//Compute the convolutions
					if(ix < bx0)
						lB = (ix - bx0) * dx;
					else
						lB = (ix - (nx - bxn - 1)) * dx;
					dDlB = DlB0 * (lB / blengthx) * (lB / blengthx);
					dD2lB = DlB0 * 2 * lB / (blengthx * blengthx);
					//alphaDlB = alphaDlB0 * fabs(lB) / blengthx;
					//alphaD2lB = SIGN(alphaDlB0 / blengthx, lB);
					alphaDlB = alphaDlB0 * (1 - fabs(lB) / blengthx);
					alphaD2lB = -SIGN(alphaDlB0 / blengthx, lB);
					
					bB = exp(-(dDlB + alphaDlB) * dt);
					aB = (1 - bB) / (dDlB + alphaDlB);
					psi2X[nBMPosX * psixlen + nBMPosZ] = bB * psi2XLa[nBMPosX * psixlen + nBMPosZ] \
						+ aB * 0.5 * (u020 + u020BXLa[nBMPosX * psixlen + nBMPosZ]);
					phi2X[nBMPosX * psixlen + nBMPosZ] = bB * phi2XLa[nBMPosX * psixlen + nBMPosZ] \
						+ aB * 0.5 * (psi2X[nBMPosX * psixlen + nBMPosZ] + psi2XLa[nBMPosX * psixlen + nBMPosZ]);
					psiX[nBMPosX * psixlen + nBMPosZ] = bB * psiXLa[nBMPosX * psixlen + nBMPosZ] \
						+ aB * 0.5 * (ux[ix * nz + iz] + uxLa[ix * nz + iz]);
					phiX[nBMPosX * psixlen + nBMPosZ] = bB * phiXLa[nBMPosX * psixlen + nBMPosZ] \
						+ aB * 0.5 * (psiX[nBMPosX * psixlen + nBMPosZ] + psiXLa[nBMPosX * psixlen + nBMPosZ]);
					EtaX[nBMPosX * psixlen + nBMPosZ] = bB * EtaXLa[nBMPosX * psixlen + nBMPosZ] \
						+ aB * 0.5 * (phiX[nBMPosX * psixlen + nBMPosZ] + phiXLa[nBMPosX * psixlen + nBMPosZ]);

					//Update the former status
					u020BXLa[nBMPosX * psixlen + nBMPosZ] = u020;
					psi2XLa[nBMPosX * psixlen + nBMPosZ] = psi2X[nBMPosX * psixlen + nBMPosZ];
					phi2XLa[nBMPosX * psixlen + nBMPosZ] = phi2X[nBMPosX * psixlen + nBMPosZ];
					psiXLa[nBMPosX * psixlen + nBMPosZ] = psiX[nBMPosX * psixlen + nBMPosZ];
					phiXLa[nBMPosX * psixlen + nBMPosZ] = phiX[nBMPosX * psixlen + nBMPosZ];
					EtaXLa[nBMPosX * psixlen + nBMPosZ] = EtaX[nBMPosX * psixlen + nBMPosZ];

					u020 = u020\
								 - 2 * dDlB * psi2X[nBMPosX * psixlen + nBMPosZ] + (dDlB * dDlB) * phi2X[nBMPosX * psixlen + nBMPosZ]\
								 - dD2lB * psiX[nBMPosX * psixlen + nBMPosZ]\
								 + dDlB * (2 * dD2lB + alphaD2lB) * phiX[nBMPosX * psixlen + nBMPosZ]\
								 - (dDlB * dDlB) * (dD2lB+alphaD2lB) * EtaX[nBMPosX * psixlen + nBMPosZ];
				}

				//Z boundaries
				if(iz < bz0 || iz >= nz - bzn)
				{
					if(iz < bz0)
						DlB0 = 1000; //Larger ValuedxB leads to stronger attenation
					else
						DlB0 = 4000;
					alphaDlB0 = 0; //Larger ValueAlphaxB leads to faster phase variation
					//Generate boundary coordinates
					GetZBoundaryMPos(ix, iz, &nBMPosX, &nBMPosZ, nx , nz);
					//Compute the convolutions
					if(iz < bz0)
						lB = (iz - bz0) * dx;
					else
						lB = (iz - (nz - bzn - 1)) * dx;
					dDlB = DlB0 * (lB / blengthz) * (lB / blengthz);
					dD2lB = DlB0 * 2 * lB / (blengthz * blengthz);
					alphaDlB = alphaDlB0 * fabs(lB) / blengthz;
					alphaD2lB = SIGN(alphaDlB0 / blengthz, lB);
					
					bB = exp(-(dDlB + alphaDlB) * dt);
					aB = (1 - bB) / (dDlB + alphaDlB);
					psi2Z[nBMPosX * psizlen + nBMPosZ] = bB * psi2ZLa[nBMPosX * psizlen + nBMPosZ] \
						+ aB * 0.5 * (u002 + u002BZLa[nBMPosX * psizlen + nBMPosZ]);
					phi2Z[nBMPosX * psizlen + nBMPosZ] = bB * phi2ZLa[nBMPosX * psizlen + nBMPosZ] \
						+ aB * 0.5 * (psi2Z[nBMPosX * psizlen + nBMPosZ] + psi2ZLa[nBMPosX * psizlen + nBMPosZ]);
					psiZ[nBMPosX * psizlen + nBMPosZ] = bB * psiZLa[nBMPosX * psizlen + nBMPosZ] \
						+ aB * 0.5 * (uz[ix * nz + iz] + uzLa[ix * nz + iz]);
					phiZ[nBMPosX * psizlen + nBMPosZ] = bB * phiZLa[nBMPosX * psizlen + nBMPosZ] \
						+ aB * 0.5 * (psiZ[nBMPosX * psizlen + nBMPosZ] + psiZLa[nBMPosX * psizlen + nBMPosZ]);
					EtaZ[nBMPosX * psizlen + nBMPosZ] = bB * EtaZLa[nBMPosX * psizlen + nBMPosZ] \
						+ aB * 0.5 * (phiZ[nBMPosX * psizlen + nBMPosZ] + phiZLa[nBMPosX * psizlen + nBMPosZ]);

					//Update the former status
					u002BZLa[nBMPosX * psizlen + nBMPosZ] = u002;
					psi2ZLa[nBMPosX * psizlen + nBMPosZ] = psi2Z[nBMPosX * psizlen + nBMPosZ];
					phi2ZLa[nBMPosX * psizlen + nBMPosZ] = phi2Z[nBMPosX * psizlen + nBMPosZ];
					psiZLa[nBMPosX * psizlen + nBMPosZ] = psiZ[nBMPosX * psizlen + nBMPosZ];
					phiZLa[nBMPosX * psizlen + nBMPosZ] = phiZ[nBMPosX * psizlen + nBMPosZ];
					EtaZLa[nBMPosX * psizlen + nBMPosZ] = EtaZ[nBMPosX * psizlen + nBMPosZ];
					//Update u002
					u002 = u002 \
						- 2 * dDlB * psi2Z[nBMPosX * psizlen + nBMPosZ] + (dDlB * dDlB) * phi2Z[nBMPosX * psizlen + nBMPosZ] \
						- dD2lB * psiZ[nBMPosX * psizlen + nBMPosZ] \
						+ dDlB * (2 * dD2lB + alphaD2lB) * phiZ[nBMPosX * psizlen + nBMPosZ] \
						- (dDlB * dDlB) * (dD2lB + alphaD2lB) * EtaZ[nBMPosX * psizlen + nBMPosZ];
				}
				uNe[ix * nz + iz] = 2 * u[ix * nz + iz] - uLa[ix * nz + iz] + (1.0f / vel[ix * nz + iz]) * (u002 + u020) * dx * dx;
			}
		}
	}
	std::swap(ux, uxLa);
	std::swap(uz, uzLa);
}

static void fillForStencil(Velocity &exvel, int halo) {
  int nxpad = exvel.nx;
  int nzpad = exvel.nz;
  int nx = nxpad - 2 * halo;
  int nz = nzpad - 2 * halo;

  std::vector<float> &vel_e = exvel.dat;

  //expand z direction first
  for (int ix = halo; ix < nx + halo; ix++) {
    for (int iz = 0; iz < halo; iz++) {
      vel_e[ix * nzpad + iz] = vel_e[ix * nzpad + halo];               // top
    }
    for (int iz = nz + halo; iz < nzpad; iz ++) {
      vel_e[ix * nzpad + iz] = vel_e[ix * nzpad + (nz + halo - 1)]; // bottom
    }
  }

  //Then x direction
  for (int iz = 0; iz < nzpad; iz++) {
    for (int ix = 0; ix < halo; ix++) {
      vel_e[ix * nzpad + iz] = vel_e[halo * nzpad + iz];               // left
    }
    for (int ix = nx + halo; ix < nxpad; ix++) {
      vel_e[ix * nzpad + iz] = vel_e[(nx + halo - 1) * nzpad + iz]; // right
    }
  }
}

static void expandForStencil(Velocity &exvel, const Velocity &v0, int halo) {
  int nx = v0.nx;
  int nz = v0.nz;
  int nzpad = nz + 2 * halo;

  const std::vector<float> &vel = v0.dat;
  std::vector<float> &vel_e = exvel.dat;

  //copy the vel into vel_e
  for (int ix = halo; ix < nx + halo; ix++) {
    for (int iz = halo; iz < nz + halo; iz++) {
      vel_e[ix * nzpad + iz] = vel[(ix - halo) * nz +  (iz - halo)];
    }
  }

  fillForStencil(exvel, halo);

}

static void expandBndry(Velocity &exvel, const Velocity &v0, int nb) {
  int nx = v0.nx;
  int nz = v0.nz;
  int nxpad = nx + 2 * nb;
#ifdef FREE
  int nzpad = nz + nb;
#else
  int nzpad = nz + 2 * nb;
#endif
  const std::vector<float> &a = v0.dat;
  std::vector<float> &b = exvel.dat;

#ifdef FREE
  /// internal
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      b[(nb + ix) * nzpad + iz] = a[ix * nz + iz];
    }
  }

  /// boundary, free surface
  for (int ix = 0; ix < nb; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      b[ix * nzpad + iz] = a[iz];                              /* left */
      b[(nb + nx + ix) * nzpad + iz] = a[(nx - 1) * nz + iz];  /* right */
    }
  }

  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < nb; iz++) {
      b[ix * nzpad + (nz + iz)] = b[ix * nzpad + (nz - 1)];  /* bottom*/
    }
  }
#else
  /// internal
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      b[(nb + ix) * nzpad + (nb + iz)] = a[ix * nz + iz];
    }
  }

  /// boundary
  for (int ix = 0; ix < nb; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      b[ix * nzpad + nb + iz] = a[iz];                              /* left */
      b[(nb + nx + ix) * nzpad + nb + iz] = a[(nx - 1) * nz + iz];  /* right */
    }
  }

  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < nb; iz++) {
      b[ix * nzpad + iz] = b[ix * nzpad + nb];												 /* top*/
      b[ix * nzpad + (nb + nz + iz)] = b[ix * nzpad + (nb + nz - 1)];  /* bottom*/
    }
  }
#endif
}

static void transvel(std::vector<float> &vel, float dx, float dt) {
  for (size_t i = 0; i < vel.size(); i ++) {
    //vel[i] = (dx * dx) / (vel[i] * vel[i] * dt * dt);
    vel[i] = (dx * dx) / (dt * dt * vel[i] * vel[i]);
  }
}

static void recoverVel(std::vector<float> &vel, float dx, float dt) {
  for (size_t i = 0; i < vel.size(); i ++) {
    vel[i] = std::sqrt(dx*dx / (dt*dt*vel[i]));
  }
}

Velocity Damp4t10d::expandDomain_notrans(const Velocity& _vel) {
  // expand for boundary, free surface
  int nb = bx0 - EXFDBNDRYLEN;

	sf_file sf_v1 = sf_output("v0_before.rsf");
	sf_putint(sf_v1, "n1", _vel.nz);
	sf_putint(sf_v1, "n2", _vel.nx);
	sf_floatwrite(const_cast<float*>(&_vel.dat[0]), _vel.nz * _vel.nx, sf_v1);

#ifdef FREE
  Velocity exvelForBndry(_vel.nx + 2 * nb, _vel.nz + nb);
  expandBndry(exvelForBndry, _vel, nb);
#else
  Velocity exvelForBndry(_vel.nx + 2 * nb, _vel.nz + 2 * nb);
  expandBndry(exvelForBndry, _vel, nb);
#endif 

	sf_file sf_v2 = sf_output("v0_after.rsf");
	sf_putint(sf_v2, "n1", _vel.nz + 2 * nb);
	sf_putint(sf_v2, "n2", _vel.nx + 2 * nb);
	sf_floatwrite(const_cast<float*>(&exvelForBndry.dat[0]), (_vel.nz + 2 * nb) * (_vel.nx + 2 * nb), sf_v2);

	/*
	sf_file sf_v1 = sf_output("v0_bndry.rsf");
	sf_putint(sf_v1, "n1", _vel.nz + nb);
	sf_putint(sf_v1, "n2", _vel.nx + 2 * nb);
	sf_floatwrite(const_cast<float*>(&exvelForBndry.dat[0]), (_vel.nz + nb) * (_vel.nx + 2 * nb), sf_v1);
	exit(1);
	*/

  //transvel(exvelForBndry.dat, dx, dt);

	/*
	sf_file sf_v1 = sf_output("v0_before.rsf");
	sf_putint(sf_v1, "n1", _vel.nz + nb);
	sf_putint(sf_v1, "n2", _vel.nx + 2 * nb);
	sf_floatwrite(const_cast<float*>(&exvelForBndry.dat[0]), (_vel.nz + nb) * (_vel.nx + 2 * nb), sf_v1);
	exit(1);
	*/

  // expand for stencil
  Velocity ret(exvelForBndry.nx+2*EXFDBNDRYLEN, exvelForBndry.nz+2*EXFDBNDRYLEN);
  expandForStencil(ret, exvelForBndry, EXFDBNDRYLEN);

	/*
	sf_file sf_v1 = sf_output("v0_after.rsf");
	sf_putint(sf_v1, "n1", ret.nz);
	sf_putint(sf_v1, "n2", ret.nx);
	sf_floatwrite(const_cast<float*>(&ret.dat[0]), ret.nx * ret.nz, sf_v1);
	exit(1);
	*/

  return ret;
}

Velocity Damp4t10d::expandDomain(const Velocity& _vel) {
  // expand for boundary, free surface
  int nb = bx0 - EXFDBNDRYLEN;

	/*
	sf_file sf_v1 = sf_output("v0_before.rsf");
	sf_putint(sf_v1, "n1", _vel.nz);
	sf_putint(sf_v1, "n2", _vel.nx);
	sf_floatwrite(const_cast<float*>(&_vel.dat[0]), _vel.nz * _vel.nx, sf_v1);
	*/

#ifdef FREE
  Velocity exvelForBndry(_vel.nx + 2 * nb, _vel.nz + nb);
  expandBndry(exvelForBndry, _vel, nb);
#else
  Velocity exvelForBndry(_vel.nx + 2 * nb, _vel.nz + 2 * nb);
  expandBndry(exvelForBndry, _vel, nb);
#endif 

	/*
	sf_file sf_v2 = sf_output("v0_after.rsf");
	sf_putint(sf_v2, "n1", _vel.nz + 2 * nb);
	sf_putint(sf_v2, "n2", _vel.nx + 2 * nb);
	sf_floatwrite(const_cast<float*>(&exvelForBndry.dat[0]), (_vel.nz + 2 * nb) * (_vel.nx + 2 * nb), sf_v2);
	exit(1);
	*/

	/*
	sf_file sf_v1 = sf_output("v0_bndry.rsf");
	sf_putint(sf_v1, "n1", _vel.nz + nb);
	sf_putint(sf_v1, "n2", _vel.nx + 2 * nb);
	sf_floatwrite(const_cast<float*>(&exvelForBndry.dat[0]), (_vel.nz + nb) * (_vel.nx + 2 * nb), sf_v1);
	exit(1);
	*/

  transvel(exvelForBndry.dat, dx, dt);

	/*
	sf_file sf_v1 = sf_output("v0_before.rsf");
	sf_putint(sf_v1, "n1", _vel.nz + nb);
	sf_putint(sf_v1, "n2", _vel.nx + 2 * nb);
	sf_floatwrite(const_cast<float*>(&exvelForBndry.dat[0]), (_vel.nz + nb) * (_vel.nx + 2 * nb), sf_v1);
	exit(1);
	*/

  // expand for stencil
  Velocity ret(exvelForBndry.nx+2*EXFDBNDRYLEN, exvelForBndry.nz+2*EXFDBNDRYLEN);
  expandForStencil(ret, exvelForBndry, EXFDBNDRYLEN);

	/*
	sf_file sf_v1 = sf_output("v0_after.rsf");
	sf_putint(sf_v1, "n1", ret.nz);
	sf_putint(sf_v1, "n2", ret.nx);
	sf_floatwrite(const_cast<float*>(&ret.dat[0]), ret.nx * ret.nz, sf_v1);
	exit(1);
	*/

  return ret;
}

void Damp4t10d::addBornwv(float *fullwv_t0, float *fullwv_t1, float *fullwv_t2, const float *exvel_m, float dt, int it, float *rp1) const {
	int nx = vel->nx;
	int nz = vel->nz;

	if(it == 0) {
		for(int i = bx0 ; i < nx - bxn ; i ++)
			for(int j = bz0 ; j < nz - bzn ; j ++)
				rp1[i * nz + j] += 2 * (fullwv_t2[i * nz + j] - fullwv_t1[i * nz + j]) / vel->dat[i * nz + j] * exvel_m[i * nz + j] / dt;
	}
	else if(it == nt - 1) {
		for(int i = bx0 ; i < nx - bxn ; i ++)
			for(int j = bz0 ; j < nz - bzn ; j ++)
				rp1[i * nz + j] += 2 * (fullwv_t1[i * nz + j] - fullwv_t0[i * nz + j]) / vel->dat[i * nz + j] * exvel_m[i * nz + j] / dt;
	}
	else {
		for(int i = bx0 ; i < nx - bxn ; i ++)
			for(int j = bz0 ; j < nz - bzn ; j ++)
				rp1[i * nz + j] -= 2 * (fullwv_t2[i * nz + j] - 2 * fullwv_t1[i * nz + j] + fullwv_t0[i * nz + j]) / vel->dat[i * nz + j] * exvel_m[i * nz + j];
	}
}

void Damp4t10d::stepForward(float* p0, float* p1) const {
  static std::vector<float> u2(vel->nx * vel->nz, 0);

  static std::vector<float> pre(vel->nx * vel->nz, 0);
	memcpy(&pre[0], p0, sizeof(float) * vel->nx * vel->nz);
  //fd4t10s_damp_zjh_2d_vtrans(p0, p1, &vel->dat[0], &u2[0], vel->nx, vel->nz, bx0);
  fd4t10s_sponge_2d_vtrans(p0, p1, &vel->dat[0], &u2[0], vel->nx, vel->nz, bx0);
	//float **pp0 = f1dto2d(p0, vel->nx, vel->nz);
	//float **pp1 = f1dto2d(p1, vel->nx, vel->nz);
	//abcone2d_apply(pp0, pp1, EXFDBNDRYLEN, abc, fd);
	//sponge2d_apply(pp0, sp, fd);
	//sponge2d_apply(pp1, sp, fd);
	
	const_cast<Damp4t10d*>(this)->applyCPML(&pre[0], p1, p0, &vel->dat[0], vel->nx, vel->nz);
	//applySponge(&p0[0], &bndr[0], vel->nx, vel->nz, bx0);
	//applySponge(&p1[0], &bndr[0], vel->nx, vel->nz, bx0);
	//applySponge(&p0[0], &vel->dat[0], &bndr[0], vel->nx, vel->nz, bx0, dt, dx);
	//applySponge(&p1[0], &vel->dat[0], &bndr[0], vel->nx, vel->nz, bx0, dt, dx);

}

void Damp4t10d::bindVelocity(const Velocity& _vel) {
  this->vel = &_vel;
}

void Damp4t10d::bindRealVelocity(const Velocity& _vel) {
  this->vel_real = &_vel;
}

void Damp4t10d::recordSeis(float* seis_it, const float* p,
    const ShotPosition& geoPos) const {

  int ng = geoPos.ns;
  int nzpad = vel->nz;

//  DEBUG() << format("ng %d") % ng;
//  float sum = 0;
  for (int ig = 0; ig < ng; ig++) {
    int gx = geoPos.getx(ig) + bx0;
    int gz = geoPos.getz(ig) + bz0;	
    int idx = gx * nzpad + gz;
    seis_it[ig] = p[idx];
//    DEBUG() << format("ig %d, idx %d, v %.20f") % ig % idx % seis_it[ig];
  }
//  DEBUG() << format("sum %.20f") % sum;

}



const Velocity& Damp4t10d::getVelocity() const {
  return *vel;
}

void Damp4t10d::stepBackward(float* p0, float* p1) const {
  static std::vector<float> u2(vel->nx * vel->nz, 0);
  fd4t10s_zjh_2d_vtrans(p0, p1, &vel->dat[0], &u2[0], vel->nx, vel->nz);
}

void Damp4t10d::addSource(float* p, const float* source,
    const ShotPosition& pos) const
{
  manipSource(p, source, pos, std::plus<float>());
}

void Damp4t10d::subSource(float* p, const float* source,
    const ShotPosition& pos) const {
  manipSource(p, source, pos, std::minus<float>());
}

void Damp4t10d::manipSource(float* p, const float* source,
    const ShotPosition& pos, boost::function2<float, float, float> op) const {
  int nzpad = vel->nz;

  for (int is = 0; is < pos.ns; is++) {
    int sx = pos.getx(is) + bx0;
    int sz = pos.getz(is) + bz0; 
    p[sx * nzpad + sz] = op(p[sx * nzpad + sz], source[is]);
//    DEBUG() << format("sx %d, sz %d, source[%d] %f") % sx % sz % is % source[is];
  }
}

void Damp4t10d::bornMaskGradient(float* grad, int H) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;
#ifdef USE_OPENMP
	#pragma omp parallel for
#endif
	for (int h = -H ; h < H ; h ++) {
		int ind = h + H;
		for (int ix = 0; ix < nxpad; ix++) {
			for (int iz = 0; iz < nzpad; iz++) {
				if (ix < bx0 || iz < bz0 || ix >= nxpad - bxn || iz >= nzpad - bzn) {
					grad[ind * nxpad * nzpad + ix  * nzpad + iz] = 0.f;
				}
			}
		}
	}
}

void Damp4t10d::maskGradient(float* grad) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;
  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < nzpad; iz++) {
      if (ix < bx0 || iz < bz0 || ix >= nxpad - bxn || iz >= nzpad - bzn) {
        grad[ix  * nzpad + iz] = 0.f;
      }
    }
  }
}

void Damp4t10d::refillBoundary(float* gradient) const {
  int nzpad = vel->nz;
  int nxpad = vel->nx;

  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < bz0; iz++) {
      gradient[ix * nzpad + iz] = gradient[ix * nzpad + bz0];           // top
    }
    for (int iz = nzpad - bzn; iz < nzpad; iz++) {
      gradient[ix * nzpad + iz] = gradient[ix * nzpad + nzpad - bzn - 1];  // bottom
    }
  }

  for (int ix = 0; ix < bx0; ix++) {                              // left
    for (int iz = 0; iz < nzpad; iz++) {
      gradient[ix * nzpad + iz] = gradient[bx0 * nzpad + iz];
    }
  }

  for (int ix = nxpad - bxn; ix < nxpad; ix++) {                        // right
    for (int iz = 0; iz < nzpad; iz++) {
      gradient[ix * nzpad + iz] = gradient[(nxpad - bxn - 1) * nzpad + iz];
    }
  }
}

void Damp4t10d::sfWriteVel(const std::vector<float> &exvel, sf_file file) const {
  //assert(exvel.size() == vel->dat.size());
  int nzpad = vel->nz;
  int nxpad = vel->nx;
  int nz = nzpad - bz0 - bzn;

  std::vector<float> vv = exvel;
  recoverVel(vv, dx, dt);
  for (int ix = bx0; ix < nxpad - bxn; ix++) {
    sf_floatwrite(&vv[ix * nzpad + EXFDBNDRYLEN], nz, file);
  }
}

void Damp4t10d::refillVelStencilBndry() {
  Velocity &exvel = getVelocity();
  fillForStencil(exvel, EXFDBNDRYLEN);
}

void Damp4t10d::FwiForwardModeling(const std::vector<float>& encSrc,
    std::vector<float>& dcal, int shot_id) const {
  int nx = getnx();
  int nz = getnz();
  int ns = getns();
  int ng = getng();

  std::vector<float> p0(nz * nx, 0);
  std::vector<float> p1(nz * nx, 0);
  ShotPosition curSrcPos = allSrcPos->clipRange(shot_id, shot_id);

  /*
	sf_file sf_v0 = sf_input("/home/cbw/fwijob/fwi_test2/v0.rsf");
	sf_floatread(const_cast<float*>(&vel->dat[0]), nz * nx, sf_v0);
  */

  /*
	sf_file sf_v1 = sf_output("v1.rsf");
	sf_putint(sf_v1, "n1", nz);
	sf_putint(sf_v1, "n2", nx);
	sf_floatwrite(const_cast<float*>(&vel->dat[0]), nz * nx, sf_v1);
  exit(1);
  */

  /*
	sf_file sf_p0 = sf_input("/home/rice/cbw/pfwi/job/p0.rsf");
	sf_floatread(const_cast<float*>(&p0[0]), nz * nx, sf_p0);
  
	sf_file sf_p1 = sf_input("/home/rice/cbw/pfwi/job/p1.rsf");
	sf_floatread(const_cast<float*>(&p1[0]), nz * nx, sf_p1);
  */

  for(int it=0; it<nt; it++) {
    addSource(&p1[0], &encSrc[it], curSrcPos);

    /*
    sf_file sf_p0 = sf_output("pp0.rsf");
    sf_putint(sf_p0, "n1", nz);
    sf_putint(sf_p0, "n2", nx);
    sf_floatwrite(const_cast<float*>(&p0[0]), nz * nx, sf_p0);

    sf_file sf_p1 = sf_output("pp1.rsf");
    sf_putint(sf_p1, "n1", nz);
    sf_putint(sf_p1, "n2", nx);
    sf_floatwrite(const_cast<float*>(&p1[0]), nz * nx, sf_p1);
    exit(1);
    */

    stepForward(&p0[0], &p1[0]);

    /*
    sf_file sf_p0 = sf_output("pp0.rsf");
    sf_putint(sf_p0, "n1", nz);
    sf_putint(sf_p0, "n2", nx);
    sf_floatwrite(const_cast<float*>(&p0[0]), nz * nx, sf_p0);

    sf_file sf_p1 = sf_output("pp1.rsf");
    sf_putint(sf_p1, "n1", nz);
    sf_putint(sf_p1, "n2", nx);
    sf_floatwrite(const_cast<float*>(&p1[0]), nz * nx, sf_p1);
    exit(1);
    */

    std::swap(p1, p0);

    /*
    sf_file sf_p0 = sf_output("pp0.rsf");
    sf_putint(sf_p0, "n1", nz);
    sf_putint(sf_p0, "n2", nx);
    sf_floatwrite(const_cast<float*>(&p0[0]), nz * nx, sf_p0);

    sf_file sf_p1 = sf_output("pp1.rsf");
    sf_putint(sf_p1, "n1", nz);
    sf_putint(sf_p1, "n2", nx);
    sf_floatwrite(const_cast<float*>(&p1[0]), nz * nx, sf_p1);
    exit(1);
    */

    recordSeis(&dcal[it*ng], &p0[0]);

    /*
    sf_file sf_p0 = sf_output("pp0.rsf");
    sf_putint(sf_p0, "n1", nz);
    sf_putint(sf_p0, "n2", nx);
    sf_floatwrite(const_cast<float*>(&p0[0]), nz * nx, sf_p0);

    sf_file sf_p1 = sf_output("pp1.rsf");
    sf_putint(sf_p1, "n1", nz);
    sf_putint(sf_p1, "n2", nx);
    sf_floatwrite(const_cast<float*>(&p1[0]), nz * nx, sf_p1);
    exit(1);
    */
  }

  /*
  sf_file sf_dcal0 = sf_output("dcal0.rsf");
  sf_putint(sf_dcal0, "n1", ng);
  sf_putint(sf_dcal0, "n2", nt);
  sf_floatwrite(const_cast<float*>(&dcal[0]), ng * nt, sf_dcal0);
  exit(1);
  */

}

void Damp4t10d::BornForwardModeling(const std::vector<float> &exvel_m, const std::vector<float>& encSrc,
    std::vector<float>& dcal, int shot_id) const {
  int nx = getnx();
  int nz = getnz();
  int ns = getns();
  int ng = getng();

  std::vector<float> fullwv(3 * nz * nx, 0);
  std::vector<float> p0(nz * nx, 0);
  std::vector<float> p1(nz * nx, 0);
  std::vector<float> rp0(nz * nx, 0);
  std::vector<float> rp1(nz * nx, 0);
	float *fullwv_t0, *fullwv_t1, *fullwv_t2, *fullwv_t;	
	fullwv_t0 = &fullwv[0];
	fullwv_t1 = &fullwv[nz * nx];
	fullwv_t2 = &fullwv[2 * nz * nx];

  ShotPosition curSrcPos = allSrcPos->clipRange(shot_id, shot_id);
	int it = 0;
	for(int it0 = 0 ; it0 < nt + 1 ; it0 ++) {
		addSource(&p1[0], &encSrc[it0], curSrcPos);
		stepForward(&p0[0], &p1[0]);
		std::swap(p1, p0);
		swap3(fullwv_t0, fullwv_t1, fullwv_t2);
		std::copy(p0.begin(), p0.end(), fullwv_t2);

		it = it0 - 1;
		if(it < 0) 
			continue;
		addBornwv(fullwv_t0, fullwv_t1, fullwv_t2, &exvel_m[0], dt, it, &rp1[0]);
		//fmMethod.addSource(&p1[0], &wlt[it], curSrcPos);
		stepForward(&rp0[0], &rp1[0]);
		std::swap(rp1, rp0);
		recordSeis(&dcal[it*ng], &rp0[0]);
	}
}

void Damp4t10d::EssForwardModeling(const std::vector<float>& encSrc,
    std::vector<float>& dcal) const {
  int nx = getnx();
  int nz = getnz();
  int ns = getns();
  int ng = getng();

  std::vector<float> p0(nz * nx, 0);
  std::vector<float> p1(nz * nx, 0);

  for(int it=0; it<nt; it++) {
    addEncodedSource(&p1[0], &encSrc[it * ns]);
    stepForward(&p0[0], &p1[0]);
    std::swap(p1, p0);
    recordSeis(&dcal[it*ng], &p0[0]);
  }
}

void Damp4t10d::bornScaleGradient(float* grad, int H) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;
  int nz = nzpad - bz0 - bzn;

#ifdef USE_OPENMP
	#pragma omp parallel for
#endif
	for (int h = -H ; h < H ; h ++) {
		int ind = h + H;
		for (int ix = bx0; ix < nxpad - bxn; ix++) {
			for (int iz = 1; iz < nz; iz++) {
				grad[ind * nxpad * nzpad + ix*nzpad + iz+bz0] *= std::sqrt(static_cast<float>(iz));
			}
		}
	}
}

void Damp4t10d::scaleGradient(float* grad) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;
  int nz = nzpad - bz0 - bzn;

  for (int ix = bx0; ix < nxpad - bxn; ix++) {
    for (int iz = 1; iz < nz; iz++) {
      grad[ix*nzpad + iz+bz0] *= std::sqrt(static_cast<float>(iz));
    }
  }
}

void Damp4t10d::removeDirectArrival(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float* data, int nt, float t_width) const {
  int half_len = t_width / dt;
  int sx = allSrcPos.getx(0) + bx0;
  int sz = allSrcPos.getz(0) + bz0;
  int gz = allGeoPos.getz(0) + bz0; // better to assume all receivers are located at the same depth

  float vel_average = 0.0;
  int gmin = (sz < gz) ? sz : gz;
  int gmax = (sz > gz) ? sz : gz;

//  printf("dt %f, half_len %d, sx %d, selav %d, gelav %d\n", dt, half_len, sx, sz, gz);
//  printf("gmin %d, gmax %d\n", gmin, gmax);

  const std::vector<float> &vv = this->vel->dat;
  int nx = this->vel->nx;
  int nz = this->vel->nz;
  for (int i = 0; i < nx; i ++) {
    for (int k = gmin; k <= gmax; k ++) {
      vel_average += vv[i * nz + k];
    }
  }
  vel_average /= nx * (gmax - gmin + 1);

  //printf("vel_average: %.20f\n", vel_average);
  //exit(1);

  int ng = allGeoPos.ns;

  for (int itr = 0; itr < ng; itr ++) {
    int gx = allGeoPos.getx(itr) + bx0;
    int gz = allGeoPos.getz(itr) + bz0;

    float dist = (gx-sx)*(gx-sx) + (gz-sz)*(gz-sz);
    int t = (int)sqrt(dist * vel_average);
    int start = t;
    int end = ((t + 2 * half_len) > nt) ? nt : (t + 2 * half_len);

    for (int j = start; j < end; j ++) {
      data[itr * nt + j] = 0.f;
    }
  }

}

Damp4t10d::Damp4t10d(const ShotPosition& _allSrcPos, const ShotPosition& _allGeoPos,
    float _dt, float _dx, float _fm, int _nb, int _nt) :
      vel(NULL), allSrcPos(&_allSrcPos), allGeoPos(&_allGeoPos),
      dt(_dt), dx(_dx), fm(_fm),  nt(_nt)
{
#ifdef FREE
  bz0 = EXFDBNDRYLEN;
#else
  bz0 = _nb + EXFDBNDRYLEN;
#endif
  bx0 = bxn = bzn = _nb + EXFDBNDRYLEN;
  bndr.resize(bx0);
  initbndr(bndr, bndr.size());
}

void Damp4t10d::addEncodedSource(float* p, const float* encsrc) const {
  this->addSource(p, encsrc, *this->allSrcPos);
}

void Damp4t10d::subEncodedSource(float* p, const float* source) const {
  this->subSource(p, source, *this->allSrcPos);
}

void Damp4t10d::recordSeis(float* seis_it, const float* p) const {
  this->recordSeis(seis_it, p, *this->allGeoPos);
}

void Damp4t10d::fwiRemoveDirectArrival(float* data, int shot_id) const {
  float t_width = 1.5 / fm;
  ShotPosition curSrcPos = allSrcPos->clipRange(shot_id, shot_id);
  this->removeDirectArrival(curSrcPos, *this->allGeoPos, data, nt, t_width);
}

void Damp4t10d::removeDirectArrival(float* data) const {
  float t_width = 1.5 / fm;
  this->removeDirectArrival(*this->allSrcPos, *this->allGeoPos, data, nt, t_width);
}

void Damp4t10d::addSource(float* p, const float* source, int is) const {

}

int Damp4t10d::getns() const {
  return allSrcPos->ns;
}

int Damp4t10d::getng() const {
  return allGeoPos->ns;
}

float Damp4t10d::getdt() const {
  return dt;
}

float Damp4t10d::getdx() const {
  return dx;
}

int Damp4t10d::getnt() const {
  return nt;
}

const ShotPosition& Damp4t10d::getAllSrcPos() const {
  return *allSrcPos;
}

const ShotPosition& Damp4t10d::getAllGeoPos() const {
  return *allGeoPos;
}

int Damp4t10d::getnx() const {
  assert(vel != NULL);
  return vel->nx;
}

int Damp4t10d::getnz() const {
  assert(vel != NULL);
  return vel->nz;
}

int Damp4t10d::getbx0() const {
  assert(vel != NULL);
  return bx0;
}

Velocity& Damp4t10d::getVelocity() {
  return *const_cast<Velocity *>(vel);
}

const std::vector<float> Damp4t10d::getVelocityDiff() const
{
	std::vector<float> vel_m(vel->nx * vel->nz, 0.0f);
	vectorMinus(vel_real->dat, vel->dat, vel_m);
	return vel_m;
}

std::vector<float> Damp4t10d::initBndryVector(int nt) const {
  if (vel == NULL) {
    ERROR() << __PRETTY_FUNCTION__ << ": you should bind velocity first";
    exit(1);
  }
  int nxpad = vel->nx;
  int nzpad = vel->nz;
//  int nx = nxpad - 2*nb;
//  int nz = nzpad - nb;

//  int nx = nxpad - bx0 - bxn;
//  int nz = nzpad - bz0 - bzn;
//
//  bndrLen = FDLEN + 1;
//  bndrSize = bndrLen * (
//             nx +   /* bottom */
//             2 * nz /* left + right */
//             );

  bndrWidth = 6;
  int nx = nxpad - (bx0 - bndrWidth + bxn - bndrWidth);
  int nz = nzpad - bz0 - bzn;

  bndrSize = bndrWidth * (
      nx +   /* bottom */
      2 * nz /* left + right */
  );

  return std::vector<float>(nt*bndrSize, 0);
}

void Damp4t10d::writeBndry(float* _bndr, const float* p, int it) const {
  /**
     * say the FDLEN = 2, then the boundary we should save is mark by (*)
     * we omit the upper layer
     *
     *    **+-------------+**
     *    **-             -**
     *    **-             -**
     *    **-             -**
     *    **-             -**
     *    **-             -**
     *    **+-------------+**
     *    *******************
     *    *******************
     *
     */
    int nxpad = vel->nx;
    int nzpad = vel->nz;

    int nx = nxpad - (bx0 - bndrWidth + bxn - bndrWidth);
    int nz = nzpad - bz0 - bzn;

    float *bndr = &_bndr[it * bndrSize];

    for (int ix = 0; ix < nx; ix++) {
      for(int iz = 0; iz < bndrWidth; iz++) {
        bndr[iz + bndrWidth*ix] = p[(ix+bx0-bndrWidth)*nzpad + (nzpad - bzn + iz)]; // bottom
      }
    }

    for (int iz = 0; iz < nz; iz++) {
      for(int ix=0; ix < bndrWidth; ix++) {
        bndr[bndrWidth*nx+iz+nz*ix]         = p[(bx0-bndrWidth + ix)*nzpad + (bz0 + iz)];   // left
        bndr[bndrWidth*nx+iz+nz*(ix+bndrWidth)] = p[(nxpad - bxn + ix)*nzpad + (bz0 + iz)];  // right
      }
    }
}

void Damp4t10d::readBndry(const float* _bndr, float* p, int it) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;

  int nx = nxpad - (bx0 - bndrWidth + bxn - bndrWidth);
  int nz = nzpad - bz0 - bzn;
  const float *bndr = &_bndr[it * bndrSize];

  for (int ix = 0; ix < nx; ix++) {
    for(int iz = 0; iz < bndrWidth; iz++) {
      p[(ix+bx0-bndrWidth)*nzpad + (nzpad - bzn + iz)] = bndr[iz + bndrWidth*ix]; // bottom
    }
  }

  for (int iz = 0; iz < nz; iz++) {
    for(int ix=0; ix < bndrWidth; ix++) {
      p[(bx0-bndrWidth + ix)*nzpad + (bz0 + iz)] = bndr[bndrWidth*nx+iz+nz*ix];   // left
      p[(nxpad - bxn + ix)*nzpad + (bz0 + iz)] = bndr[bndrWidth*nx+iz+nz*(ix+bndrWidth)];  // right
    }
  }
}
