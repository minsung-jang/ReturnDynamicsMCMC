/**
 * Implementation of MCMC algorithm for Return Dynamics in C++
 *
 * Copyright (C) 2019 Minsung Jang < msjang at iastate dot edu >
 *
 */

#include "ReturnDynamicsMCMC.h"

void ReturnDynamicsMCMC::getIntervalSize(VectorXd& mu_vec, VectorXd& Y_vec)
{
    assert(Y_vec.size()==mu_vec.size());
    T = Y_vec.size()-1; // t=0,1,2,...,T, total data size: T+1
}

void ReturnDynamicsMCMC::getInitialGuess(ParamSet& init_guess)
{
    E_mu = init_guess.E_mu;
    beta = init_guess.beta;
    rho = init_guess.rho;
    sig_y = init_guess.sig_y;
    sig_mu = init_guess.sig_mu;
}

void ReturnDynamicsMCMC::getHyperParam(HyperParamSet& hyper_param)
{
    prior_e=hyper_param.e;
    prior_E=hyper_param.E;
    prior_f=hyper_param.f;
    prior_F=hyper_param.F;
    prior_alpha=hyper_param.alpha;
    prior_beta=hyper_param.beta;
    prior_a=hyper_param.a;
    prior_b=hyper_param.b;
}

void ReturnDynamicsMCMC::initSimulation(int& n_path, int& length,
                                        ParamSet& simul_theta, double mu_0,
                                        bool OUTPUT_FILE){

    cout << "========== Initialize a simulation ==========" << endl;
    cout << "... # Sample Path : " << n_path << endl;
    cout << "... # Length      : " << length << endl;
    cout << "... Parameters: E_mu = " << simul_theta.E_mu << ", " <<
            "beta = " << simul_theta.beta << ", " <<
            "sig_y = " << simul_theta.sig_y << ", " <<
            "sig_mu = " << simul_theta.sig_mu << ", " <<
            "rho = " << simul_theta.rho << endl;

    true_param = new ParamSet;
    true_param = &simul_theta;

    // Random Generator
    for (int k=0; k < n_path; k++){

        VectorXd Y(length+1);
        VectorXd mu(length+1);
        mu[0]=mu_0;
        Y[0]=mu_0;

        // Output file per path
        string file_name = "../simulation/MCMC_simulation_[" + to_string(k+1) + "].csv";
        ofstream outputFile(file_name);
        outputFile << "Parameters: E_mu = " << simul_theta.E_mu << ", " <<
                      "beta = " << simul_theta.beta << ", " <<
                      "sig_y = " << simul_theta.sig_y << ", " <<
                      "sig_mu = " << simul_theta.sig_mu << ", " <<
                      "rho = " << simul_theta.rho << endl;
        outputFile << "Y , mu" << endl;

        // Generate two indepedent standard normal
        default_random_engine generator1;
        default_random_engine generator2;
        normal_distribution<double> std_normal (0.0, 1.0);

        // Generate Y and mu
        for (int i=1; i <= length; i++){

            // Generate epsilon_Y and epsilon_mu using two independent stanard normals
            double z=std_normal(generator1);
            double eps_mu=std_normal(generator2);
            double eps_y=sqrt(1-pow(simul_theta.rho,2))*z+simul_theta.rho*eps_mu;

            // Conversion to Y and mu based on the model
            Y[i] = mu[i-1] + simul_theta.sig_y * eps_mu;
            mu[i] = simul_theta.E_mu + simul_theta.beta*(mu[i-1]-simul_theta.E_mu) + simul_theta.sig_y * eps_y;

            // Write to the output file
            outputFile << Y[i] << ", " << mu[i] << endl;
        }
        outputFile.close();

        // Save samples to memory
        Y_sample.push_back(Y);
        mu_sample.push_back(mu);
    }
}

void ReturnDynamicsMCMC::runMCMC(int& n_path, long& n_iter,
                                 ParamSet& init_guess, HyperParamSet& hyper_param, bool LOG_FILE){

    cout << "========== Run MCMC ==========" << endl;

    getInitialGuess(init_guess);
    getHyperParam(hyper_param);

    // Write output file
    ofstream outputFile("../MCMC_final_output.csv");
    outputFile << "Sample No. , E_mu, beta, sig_y, sig_mu, rho" << endl;
    if (true_param!=nullptr){
        cout << "... Parameters: E_mu = " << true_param->E_mu << ", " <<
                "beta = " << true_param->beta << ", " <<
                "sig_y = " << true_param->sig_y << ", " <<
                "sig_mu = " << true_param->sig_mu << ", " <<
                "rho = " << true_param->rho << endl;
        outputFile << " True value: ," << true_param->E_mu << ", " <<
                   true_param->beta <<  ", " <<
                   true_param->sig_y << ", " <<
                   true_param->sig_mu << ", " <<
                   true_param->rho << endl;
    }else {
        outputFile << "There is no true parameter...,N/A,N/A,N/A,N/A,N/A" << endl;
    }

    // Per sample
    if (Y_sample.empty()){
        outputFile.close();
        throw runtime_error("Input Y is empty");
    }

    for (unsigned int s=0; s < uint(n_path); s++){

        VectorXd Y_input = Y_sample[s];
        VectorXd mu_input = mu_sample[s];

        // Determine T
        getIntervalSize(mu_input,Y_input);

        // Set the parameters with the initial guess


        // Log file
        string logFileName = "../logfile/MCMC_log_[" + to_string(s+1) + "].csv";
        ofstream logFile(logFileName);
        logFile << "E_mu, beta, sig_y, sig_mu, rho" << endl;
        logFile << E_mu << ", " <<
                   beta << ", " <<
                   sig_y << ", " <<
                   sig_mu << ", " <<
                   rho << endl;

        long seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator_E_mu(seed);

        for (long k=0; k < n_iter; k++){
            E_mu = posterior_E_mu(generator_E_mu, mu_input, Y_input,false);
            if(isinf(E_mu))
                throw runtime_error("E_mu diverges to infinity");
            if(isnan(E_mu))
                throw runtime_error("E_mu is not real");

            logFile << E_mu << ", " <<
                       beta << ", " <<
                       sig_y << ", " <<
                       sig_mu << ", " <<
                       rho << endl;
        }
        logFile.close();

        // Write final output file
        outputFile << s+1 << " , " <<
                      E_mu << ", " <<
                      beta << ", " <<
                      sig_y << ", " <<
                      sig_mu << ", " <<
                      rho << ", " <<endl;
        cout << "... Sample path No." << s+1 << " is completed" << endl;
    }
    outputFile.close();
}
