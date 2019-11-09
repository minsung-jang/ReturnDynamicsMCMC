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
// Custom
#include "MyStatDist.h"

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
    // AUX
    bool _CHECK_STREAM;
    default_random_engine main_generator;

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
        _CHECK_STREAM=true;
        long seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator(seed);
        main_generator = generator;
    }
    ~ReturnDynamicsMCMC(){}

    // Member Functions
    // ... Posterior: pre-computation
    double stdNormalCDF(double z);
    VectorXd compute_tp1(VectorXd& vec);
    VectorXd compute_tp0(VectorXd& vec);
    VectorXd compute_disp_vec(VectorXd& vec1, VectorXd& vec2, double coeff=1.0);
    VectorXd compute_disp_const(VectorXd& vec, double coeff);
    double compute_K(VectorXd& mu_vec, VectorXd& Y_vec);
    // ... Posterior
    double posterior_E_mu(default_random_engine& generator, VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG=false);
    double posterior_beta(default_random_engine& generator, VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG=false);
    double posterior_sig_y(default_random_engine& generator,
                           VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG=false);
    vector<double> posterior_phi_omega(default_random_engine& generator1,
                                       VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG=false);
    VectorXd posterior_mu(VectorXd& mu_vec, VectorXd& Y_vec);
    // ... MCMC
    void getIntervalSize(VectorXd& mu_vec, VectorXd& Y_vec);
    void getInitialGuess(ParamSet& init_guess);
    void getHyperParam(HyperParamSet& hyper_param);
    void loadCSV(string file_name);
    void initSimulation(int& n_path, int& length,ParamSet& simul_thetha, bool OUTPUT_FILE=true);
    void runMCMC(int& n_path, long& n_iter, ParamSet& init_guess, HyperParamSet& hyper_param, bool LOG_FILE=true);
};

#endif
