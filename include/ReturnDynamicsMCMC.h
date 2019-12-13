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
    string input_file_name;
    bool REFRESH_SIMUL;
    ParamSet simul_theta;
    ParamSet init_guess;
    HyperParamSet hyper_param;
    int n_simul_path;
    int n_simul_length;
    long n_iter;
    double burn_in_ratio;
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

    // TEST
    vector<double> E_mu_tmp;
    vector<double> beta_tmp;
    vector<double> sig_y_tmp;
    vector<double> sig_mu_tmp;
    vector<double> rho_tmp;
    VectorXd mu_tmp;

    // Contructor & Destructor
    ReturnDynamicsMCMC(){
        true_param=nullptr;
    }
    ~ReturnDynamicsMCMC(){
        delete true_param;
    }

    // Member Functions
    void getInitialSetup(bool DEBUG=false);
    void displaySimulationSetup();
    void displayCommonSetup();
    void getInitialGuess();
    void getHyperParam();

    // ... Posterior: pre-computation
    double stdNormalCDF(double z);
    VectorXd compute_tp1(VectorXd& vec);
    VectorXd compute_tp0(VectorXd& vec);
    VectorXd compute_disp_vec(VectorXd& vec1, VectorXd& vec2, double coeff=1.0);
    VectorXd compute_disp_const(VectorXd& vec, double coeff);
    double compute_K(VectorXd& mu_vec, VectorXd& Y_vec);

    // ... Posterior
    double posterior_E_mu(mt19937_64& generator, VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG=false);
    double posterior_beta(mt19937_64& generator, VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG=false);
    double posterior_sig_y(mt19937_64& generator,
                           VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG=false);
    vector<double> posterior_phi_omega(mt19937_64& generator1,
                                       VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG=false);
    void posterior_mu(mt19937_64& generator, VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG=false);

    // ... MCMC
    void getIntervalSize(VectorXd& mu_vec, VectorXd& Y_vec);
    void loadInputData(string& file_name);
    void initSimulation();
    void runMCMC(bool DEBUG=false);

    // ... TEST
    //void loadTestParam(string& file_name);
    //void loadTestmu(string& file_name);

};

#endif
