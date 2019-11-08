/**
 * Implementation of MCMC algorithm for Return Dynamics in C++
 *
 * Copyright (C) 2019 Minsung Jang < msjang at iastate dot edu >
 *
 */

#ifndef RETURN_DYANMICS_MCMC_H
#define RETURN_DYNAMICS_MCMC_H

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>
#include <string>
#include <assert.h>
// Eigen
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

struct ParamSet{
    double E_mu;
    double beta;
    double sig_y;
    double sig_mu;
    double rho;
};

struct HyperParamSet{
    double e;
    double E;
    double f;
    double F;
    double alpha;
    double beta;
    double a;
    double b;
};

class ReturnDynamicsMCMC{
private:
    long T;
public:
    // Input sample
    vector<VectorXd> Y_sample;
    vector<VectorXd> mu_sample;

    // Simulation
    ParamSet* true_param;

    // Posterior
    double E_mu;
    double beta;
    double sig_y;
    double sig_mu;
    double rho;

    // Prior
    double prior_E;
    double prior_e;
    double prior_F;
    double prior_f;
    double prior_alpha;
    double prior_beta;
    double prior_a;
    double prior_b;

    // Contructor & Destructor
    ReturnDynamicsMCMC(){
        true_param=nullptr;
    }
    ~ReturnDynamicsMCMC(){}

    // Member Functions
    // ... Posterior: pre-computation
    VectorXd compute_tp1(VectorXd& vec);
    VectorXd compute_tp0(VectorXd& vec);
    VectorXd compute_disp_vec(VectorXd& vec1, VectorXd& vec2, double coeff=1.0);
    VectorXd compute_disp_const(VectorXd& vec, double coeff);
    // ... Posterior
    double posterior_E_mu(default_random_engine& generator, VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG=false);
    // ... MCMC
    void getIntervalSize(VectorXd& mu_vec, VectorXd& Y_vec);
    void getInitialGuess(ParamSet& init_guess);
    void getHyperParam(HyperParamSet& hyper_param);
    void initSimulation(int& n_path, int& length,ParamSet& simul_thetha, double mu_0=0.0, bool OUTPUT_FILE=true);
    void runMCMC(int& n_path, long& n_iter, ParamSet& init_guess, HyperParamSet& hyper_param, bool LOG_FILE=true);
};

#endif
