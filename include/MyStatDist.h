/**
 * Implementation of MCMC algorithm for Return Dynamics in C++
 *
 * Copyright (C) 2019 Minsung Jang < msjang at iastate dot edu >
 *
 */

#ifndef MYSTATDIST_H
#define MYSTATDIST_H

#include <cmath>
#include <vector>

class MyStatDist{
public:
    MyStatDist();
    virtual ~MyStatDist();
    virtual double inv_cdf(const double& quantile) const = 0;
};

class StandardNormalDistribution : public MyStatDist {
 public:
  StandardNormalDistribution();
  virtual ~StandardNormalDistribution();
  // Inverse cumulative distribution function
  virtual double inv_cdf(const double& quantile) const;
};

#endif
