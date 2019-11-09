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

  // Inverse cumulative distribution function (aka the probit function)
  virtual double inv_cdf(const double& quantile) const;
};

#endif
