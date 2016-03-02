/*
 * spongabc4d.h
 *
 *  Created on: Feb 28, 2016
 *      Author: rice
 */

#ifndef SRC_FM2D_SPONGABC4D_H_
#define SRC_FM2D_SPONGABC4D_H_

#include "velocity.h"
#include "shot-position.h"
#include "i-modeling.h"

class SpongAbc4d : public IModeling {
public:
  SpongAbc4d(float _dt, float _dx, float _dz, int _nb);

public: // overrice
  Velocity transformVelocityForModeling(const Velocity &v0) const;

public:
  void stepForward(float *p0, float *p1) const;
  void addSource(float *p, const float *source, int ns, const int *sxz, int snz);
  void addSource(float *p, const float *source, const ShotPosition &pos);
  void recordSeis(float *seis_it, const float *p, const ShotPosition &geoPos);
  std::vector<float> initBndryVector(int nt) const;

private:
  void applySponge(float *p) const;
  void initCoeff();
  void initbndr();

private:
  std::vector<float> bndr;
  float dt, dx, dz;
  int nb;
  float c0, c11, c12, c21, c22;
};

#endif /* SRC_FM2D_SPONGABC4D_H_ */