#include "ReturnDynamicsMCMC.h"

VectorXd ReturnDynamicsMCMC::compute_tp1(VectorXd& vec)
{
    VectorXd rtn = vec;
    rtn = rtn.tail(T);

    return rtn;
}

VectorXd ReturnDynamicsMCMC::compute_tp0(VectorXd& vec)
{
    VectorXd rtn = vec;
    rtn = rtn.head(T);

    return rtn;
}

VectorXd ReturnDynamicsMCMC::compute_disp_vec(VectorXd& vec1, VectorXd& vec2, double coeff)
{
    VectorXd rtn=compute_tp1(vec1)-coeff*compute_tp0(vec2);

    return rtn;
}

VectorXd ReturnDynamicsMCMC::compute_disp_const(VectorXd& vec, double coeff)
{
    long N = vec.size();
    VectorXd rtn(N);
    rtn.setOnes();
    rtn = vec-coeff*rtn;

    return rtn;
}

double ReturnDynamicsMCMC::posterior_E_mu(default_random_engine& generator, VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG)
{
    if (DEBUG)
        cout << "[Posterior] E_mu" << endl;
    if (DEBUG)
        cout << "... input [ E_mu=" << E_mu << ", beta = " << beta << ", rho = " <<
                rho << ", sig_y = " << sig_y << ", sig_mu = " << sig_mu << " ] " <<endl;

    // Vectors
    VectorXd C_tp1=compute_disp_vec(Y_vec, mu_vec);

    VectorXd D_tp1=compute_disp_vec(mu_vec, mu_vec, beta);

    // Costants
    double B = (1.0 - beta) / sig_mu;
    double B_sq =B*B;
    double psi = 1.0 / (1.0 - rho*rho);
    double E_sq = prior_E*prior_E;
    double e_sq = prior_e*prior_e;

    // W,S
    double W =  1.0 / E_sq + psi*B_sq*double(T);
    double S = e_sq / E_sq - (rho*psi*B)*C_tp1.sum() / sig_y +(psi*B)*D_tp1.sum() / sig_mu;

    // Random Generator
    if (DEBUG)
        cout << "... dist: ~ N( " << S/W << ", " << 1.0/W << " )" << endl;

    normal_distribution<double> dist (S/W, sqrt(1.0/W));
    double rtn = dist(generator);
    if (DEBUG)
        cout << "... sampled : " << rtn << endl;

    return rtn;
}
