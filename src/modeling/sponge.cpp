#include "sponge.h"

void Sponge::initbndr(int nb) {
	bndr.resize(nb);
	bndr.assign(nb, 0.0f);
	float sponge_coef = .003737;
  for(int ib=0;ib<nb;ib++){
    float tmp=sponge_coef*(nb-ib-1);
    bndr[ib]=std::exp(-tmp*tmp);
  }
}

void Sponge::applySponge(float* p, const float *vel, int nx, int nz, int nb, float dt, float dx, int freeSurface) {
  for(int ib=0; ib<nb; ib++) {
    float w = bndr[ib];

    int ibz = nz-ib-1;
    for(int ix=0; ix<nx; ix++) {
			if(!freeSurface) {
				p[ix * nz + ib] *= w;
			}
      p[ix * nz + ibz] *= w;
    }

    int ibx = nx-ib-1;
    for(int iz=0; iz<nz; iz++) {
      p[ib  * nz + iz] *= w;
      p[ibx * nz + iz] *= w;
    }
  }
}

