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

        rdMCMC->getInitialSetup();
        rdMCMC->runMCMC();

    } catch (const exception& e) {
        cout << "ERROR OCCURRED: " << e.what() << endl;
    }

    return 0;
}
