/**
 * Implementation of MCMC algorithm for Return Dynamics in C++
 *
 * Copyright (C) 2019 Minsung Jang < msjang at iastate dot edu >
 *
 */

#include "ReturnDynamicsMCMC.h"

int main()
{
    ReturnDynamicsMCMC* rdMCMC = new ReturnDynamicsMCMC();

    try {

        //Legacy commands
//        string ref_param_file = "../ref_result.csv";
//        string ref_mu_file ="../ref_mu.csv";
//        rdMCMC->loadTestParam(ref_param_file);
//        rdMCMC->loadTestmu(ref_mu_file);

        rdMCMC->getInitialSetup();
        rdMCMC->runMCMC();

    } catch (const exception& e) {
        cout << "ERROR OCCURRED: " << e.what() << endl;
    }

    return 0;
}
