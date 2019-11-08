/**
 * Implementation of MCMC algorithm for Return Dynamics in C++
 *
 * Copyright (C) 2019 Minsung Jang < msjang at iastate dot edu >
 *
 */

#include "ReturnDynamicsMCMC.h"

int main()
{
    clock_t time_begin=clock();

    ReturnDynamicsMCMC* rdMCMC = new ReturnDynamicsMCMC();

    ParamSet simul_theta;
    ParamSet init_guess;
    HyperParamSet hyper_param;
    int n_simul_path = 1;
    int n_simul_length = 5000;
    long n_iter=2000;
    bool REFRESH_SIMUL=true;

    // Simulation
    if (REFRESH_SIMUL){
        // Set simulation parameter
        simul_theta.E_mu=0.1;
        simul_theta.beta=0.9;
        simul_theta.rho=-0.5;
        simul_theta.sig_y=0.15;
        simul_theta.sig_mu=0.1;
        //
        rdMCMC->initSimulation(n_simul_path, n_simul_length, simul_theta);
    }else {
        // Load csv data
    }

    // MCMC
    // Initial guess for MCMC
    init_guess.E_mu =0.05;
    init_guess.beta=0.5;
    init_guess.rho=-0.8;
    init_guess.sig_y=1.0;
    init_guess.sig_mu=1.0;

    // Set Hyperparameters
    hyper_param.e=0.0;
    hyper_param.E=1.0;
    hyper_param.f=0.0;
    hyper_param.F=1.0;
    hyper_param.alpha=100;
    hyper_param.beta=2;
    hyper_param.a=2;
    hyper_param.b=200;

    try {
        rdMCMC->runMCMC(n_simul_path, n_iter, init_guess, hyper_param);
    } catch (const exception& e) {
        cout << e.what() << endl;
    }

    clock_t time_end=clock();

    double elapsed_time = double(time_end - time_begin) / CLOCKS_PER_SEC;

    cout << "Elapsed Time: " << elapsed_time <<" [sec]" << endl;

    return 0;
}
