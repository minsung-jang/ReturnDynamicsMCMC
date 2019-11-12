/**
 * Implementation of MCMC algorithm for Return Dynamics in C++
 *
 * Copyright (C) 2019 Minsung Jang < msjang at iastate dot edu >
 *
 */

#include "ReturnDynamicsMCMC.h"

void ReturnDynamicsMCMC::getIntervalSize(VectorXd& mu_vec, VectorXd& Y_vec)
{
    if(Y_vec.size()!=mu_vec.size()){
        cout << "Y length : " << Y_vec.size() << endl;
        cout << "mu length : " << mu_vec.size() << endl;
        throw runtime_error("The # of Y is inconsistent with that of mu");
    }

    T = Y_vec.size()-1; // t=0,1,2,...,T, total data size: T+1
}

void ReturnDynamicsMCMC::displaySimulationSetup()
{

    cout << "... # Path  : " << n_simul_path << endl;
    cout << "... # Length : " << n_simul_length << endl;
    cout << "... Parameters: " << endl;
    cout << "... --- E_mu = " << simul_theta.E_mu << endl;
    cout << "... --- beta = " << simul_theta.beta << endl;
    cout << "... --- sigma_y = " << simul_theta.sig_y << endl;
    cout << "... --- sigma_mu = " << simul_theta.sig_mu << endl;
    cout << "... --- rho = " << simul_theta.rho << endl;
}

void ReturnDynamicsMCMC::displayCommonSetup()
{
    cout << "... Total iteration : " << n_iter << endl;
    cout << "... Burn-in period ratio : " << burn_in_ratio << endl;
    cout << "... Initial guess: " << endl;
    cout << "... --- E_mu = " << init_guess.E_mu << endl;
    cout << "... --- beta = " << init_guess.beta << endl;
    cout << "... --- sigma_y = " << init_guess.sig_y << endl;
    cout << "... --- sigma_mu = " << init_guess.sig_mu << endl;
    cout << "... --- rho = " << init_guess.rho << endl;
    cout << "... Hyperparameters: " << endl;
    cout << "... --- e = " << hyper_param.e << endl;
    cout << "... --- E = " << hyper_param.E << endl;
    cout << "... --- f = " << hyper_param.f << endl;
    cout << "... --- F = " << hyper_param.F << endl;
    cout << "... --- alpha* = " << hyper_param.alpha << endl;
    cout << "... --- beta* = " << hyper_param.beta << endl;
    cout << "... --- a = " << hyper_param.a << endl;
    cout << "... --- b = " << hyper_param.b << endl;
}

void ReturnDynamicsMCMC::getInitialGuess()
{
    E_mu = init_guess.E_mu;
    beta = init_guess.beta;
    rho = init_guess.rho;
    sig_y = init_guess.sig_y;
    sig_mu = init_guess.sig_mu;
}

void ReturnDynamicsMCMC::getHyperParam()
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

void ReturnDynamicsMCMC::getInitialSetup(bool DEBUG)
{
    ifstream file("../setup.txt");
    if(!file.is_open())
        throw runtime_error("there is no setup file (setup.txt).");

    string str;
    string token;
    getline(file, str);
    getline(file, str); // skip the first, second line
    int line_count=2;
    while ( getline(file, str) ){
        line_count++;
        auto pos = str.find(':');
        if(pos != string::npos)
            str.erase(0, pos + 1);
        str.erase(remove(str.begin(),str.end(),' '),str.end());
        istringstream iss(str);
        while ( getline(iss, token) )
        {
            // Simulation
            // ... whether or not perform simulation
            if (line_count==3){
                string string_refesh_simul;
                string_refesh_simul=token;
                if (string_refesh_simul=="yes" || string_refesh_simul=="Yes" ||
                        string_refesh_simul=="y" || string_refesh_simul=="Y"){
                    REFRESH_SIMUL=true;
                    cout << "[ Simulation Mode ]" << endl;
                }

            }

            // Load input data
            if (line_count==4 && REFRESH_SIMUL==false){
                input_file_name=token;
                string input_file_name_full = "../" + input_file_name;
                ifstream input_file(input_file_name_full);

                if(!input_file.is_open())
                    throw runtime_error("there is no input data with such name, '"+input_file_name+"'");

                loadInputData(input_file_name_full);
                cout << "[ Use an observation ]" << endl;
                cout << "... Loaded input data: " << input_file_name_full << endl;
            }

            // ... Simulation Path #
            if (line_count==8)
            {
                n_simul_path = stoi(token);
            }

            // ... Simulation Length
            if (line_count==9)
                n_simul_length = stoi(token);

            // ... Parameters
            if (line_count==11)
                simul_theta.E_mu = stod(token);
            if (line_count==12)
                simul_theta.beta = stod(token);
            if (line_count==13)
                simul_theta.sig_y = stod(token);
            if (line_count==14)
                simul_theta.sig_mu = stod(token);
            if (line_count==15)
                simul_theta.rho = stod(token);

            // MCMC implementaion
            // ... Iteration #
            if (line_count==19)
                n_iter = stoi(token);
            if (line_count==20)
                burn_in_ratio = stod(token);
            if (line_count==22)
                init_guess.E_mu = stod(token);
            if (line_count==23)
                init_guess.beta = stod(token);
            if (line_count==24)
                init_guess.sig_y = stod(token);
            if (line_count==25)
                init_guess.sig_mu = stod(token);
            if (line_count==26)
                init_guess.rho = stod(token);

            // ... Hyperparameters
            if (line_count==28)
                hyper_param.e = stod(token);
            if (line_count==29)
                hyper_param.E = stod(token);
            if (line_count==30)
                hyper_param.f = stod(token);
            if (line_count==31)
                hyper_param.F = stod(token);
            if (line_count==32)
                hyper_param.alpha = stod(token);
            if (line_count==33)
                hyper_param.beta = stod(token);
            if (line_count==34)
                hyper_param.a = stod(token);
            if (line_count==35)
                hyper_param.b = stod(token);
        }
    }

    if (REFRESH_SIMUL){
        initSimulation();

        if (DEBUG)
            displaySimulationSetup();
    }

    if (DEBUG)
        displayCommonSetup();

    getInitialGuess();
    getHyperParam();

}

void ReturnDynamicsMCMC::loadInputData(string& file_name)
{
    vector<double> Y_tmp;

    ifstream file(file_name);

    Y_sample.clear();
    string str;
    getline(file, str); // skip the first line
    while (getline(file, str)){
        istringstream iss(str);
        vector<string> tokens;
        string token;
        while ( getline(iss, token,',') ){
            tokens.push_back(token);
        }
        Y_tmp.push_back(stod(tokens[0]));
        tokens.clear();
    }

    VectorXd out(Y_tmp.size());
    for (unsigned int i=0; i < Y_tmp.size(); i++)
    {
        out[i]=Y_tmp[i];
    }
    Y_sample.push_back(out);

}

void ReturnDynamicsMCMC::initSimulation(){

    true_param = new ParamSet;
    true_param = &simul_theta;

    // Random Generator
    for (int k=0; k < n_simul_path; k++){

        VectorXd Y(n_simul_length+1);
        VectorXd mu(n_simul_length+1);
        mu[0]=simul_theta.E_mu;
        Y[0]=simul_theta.E_mu;

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

        // t=0
        outputFile << Y[0] << ", " << mu[0] << endl;

        // Generate Y and mu for t=1,...,T
        for (int i=1; i <= n_simul_length; i++){

            // Generate epsilon_Y and epsilon_mu using two independent stanard normals
            double z=std_normal(generator1);
            double eps_y=std_normal(generator2);
            double eps_mu=sqrt(1-simul_theta.rho*simul_theta.rho)*z+simul_theta.rho*eps_y;

            // Conversion to Y and mu based on the model
            Y[i] = mu[i-1] + simul_theta.sig_y * eps_y;
            mu[i] = simul_theta.E_mu + simul_theta.beta*(mu[i-1]-simul_theta.E_mu) + simul_theta.sig_mu * eps_mu;

            // Write to the output file
            outputFile << Y[i] << ", " << mu[i] << endl;
        }
        outputFile.close();

        // Save samples to memory
        Y_sample.push_back(Y);
        mu_sample.push_back(mu);
    }
}

void ReturnDynamicsMCMC::runMCMC(bool DEBUG){

    clock_t time_begin=clock();

    cout << "[ ================= Run MCMC ================= ]" << endl;

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
        outputFile << "True parameters are unkown...,N/A,N/A,N/A,N/A,N/A" << endl;
    }

    // Check the input
    if (Y_sample.empty()){
        outputFile.close();
        throw runtime_error("Input Y is empty");
    }
    ulong n_path = Y_sample.size();

    // Per sample
    for (unsigned long s=0; s < n_path; s++){
        cout << ">> Sample path No." << s+1 << " has started" << endl;

        // Clock
        clock_t time_begin_per_sample=clock();

        // Input
        VectorXd Y_input = Y_sample[s];
        VectorXd mu_input = Y_input;

        Y_current = Y_input;
        mu_current = mu_input;

        // Determine T
        getIntervalSize(mu_input,Y_input);

        // Log file
        string logFileName = "../logfile/MCMC_log_[" + to_string(s+1) + "].csv";
        ofstream logFile(logFileName);
        logFile << "E_mu, beta, sig_y, sig_mu, rho" << endl;
        logFile << E_mu << ", " <<
                   beta << ", " <<
                   sig_y << ", " <<
                   sig_mu << ", " <<
                   rho << endl;

        random_device seed_gen;
        default_random_engine main_generator(seed_gen());

        for (long k=0; k < n_iter; k++){

            if (DEBUG)
                cout << endl << ">>> [ " << (k+1) << " ]-th iteration "<< endl << endl;

            // [1] E_mu
            E_mu = posterior_E_mu(main_generator, mu_input, Y_input,DEBUG);
            if(isinf(E_mu))
                throw runtime_error("E_mu diverges to infinity");
            if(isnan(E_mu))
                throw runtime_error("E_mu is not real");
            logFile << E_mu << ", " ;

            // [2] beta
            beta = posterior_beta(main_generator, mu_input, Y_input,DEBUG);
            if(isinf(beta))
                throw runtime_error("beta diverges to infinity");
            if(isnan(beta))
                throw runtime_error("beta is not real");
            logFile << beta << ", " ;

            // [3] sig_y
            sig_y = posterior_sig_y(main_generator, mu_input, Y_input,DEBUG);
            if(isinf(sig_y))
                throw runtime_error("sig_y diverges to infinity");
            if(isnan(sig_y))
                throw runtime_error("sig_y is not real");
            logFile << sig_y << ", ";

            //[4]  sig_mu and rho
            vector<double> phi_omega = posterior_phi_omega(main_generator,
                                                           mu_input, Y_input,DEBUG);
            sig_mu = phi_omega[0];
            if(isinf(sig_mu))
                throw runtime_error("sig_mu diverges to infinity");
            if(isnan(sig_mu))
                throw runtime_error("sig_mu is not real");
            rho = phi_omega[1];
            if(isinf(rho))
                throw runtime_error("rho diverges to infinity");
            if(isnan(rho))
                throw runtime_error("rho is not real");
            logFile << sig_mu << ", " << rho << endl;

            // [5] mu
            mu_input = posterior_mu(main_generator, mu_input, Y_input,DEBUG);
        }
        clock_t time_end_per_sample = clock();

        // Elapsed time
        double elapsed_time_per_sample = double(time_end_per_sample - time_begin_per_sample) / CLOCKS_PER_SEC;
        cout << "... Elapsed Time : " << elapsed_time_per_sample <<" [sec] for " <<
                n_iter << " iterations" << endl;
        logFile << endl << "Elapsed Time : " << elapsed_time_per_sample <<" [sec] for " <<
                n_iter << " iterations" << endl;

        logFile.close();

        // Write a final output file
//        outputFile << s+1 << " , " <<
//                      E_mu << ", " <<
//                      beta << ", " <<
//                      sig_y << ", " <<
//                      sig_mu << ", " <<
//                      rho << ", " <<endl;
        cout << ">> Sample path No." << s+1 << " has been completed" << endl;
    }

    clock_t time_end=clock();

    double elapsed_time = double(time_end - time_begin) / CLOCKS_PER_SEC;
    cout << "Elapsed Time in total: " << elapsed_time <<" [sec]" << endl;
//    outputFile << endl << "Elapsed Time in total: " << elapsed_time <<" [sec]" << endl;

//    outputFile.close();
}
