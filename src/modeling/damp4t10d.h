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
}
#include <boost/function.hpp>
#include "velocity.h"
#include "shot-position.h"

class Damp4t10d {
public:
  Damp4t10d(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float dt, float dx, int nb);

  Velocity expandDomain(const Velocity &vel);


  void stepForward(float *p0, float *p1) const;
  void stepBackward(float *p0, float *p1) const;
  void bindVelocity(const Velocity &_vel);
  void addSource(float *p, const float *source, const ShotPosition &pos) const;
  void addSource(float *p, const float *source, int is) const;
  void addEncodedSource(float *p, const float *encsrc) const;
  void recordSeis(float *seis_it, const float *p) const;
  const Velocity &getVelocity() const;
  void maskGradient(float *grad) const;
  void refillBoundary(float *vel) const;
  void sfWriteVel(sf_file file) const;

  void removeDirectArrival(float* data, int nt, float t_width) const;
  void subEncodedSource(float *p, const float *source) const;

public:
  int getTotalSrc() const;
  int getTotalGeo() const;

private:
  void manipSource(float *p, const float *source, const ShotPosition &pos, boost::function2<float, float, float> op) const;
  void recordSeis(float *seis_it, const float *p, const ShotPosition &geoPos) const;
  void subSource(float *p, const float *source, const ShotPosition &pos) const;
  void removeDirectArrival(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float* data, int nt, float t_width) const;

private:
  const static int FDLEN = 5;

private:
  const Velocity *vel;
  const ShotPosition *allSrcPos;
  const ShotPosition *allGeoPos;
  float dt;
  float dx;
  int nb;
};

#endif /* SRC_FM2D_DAMP4T10D_H_ */
