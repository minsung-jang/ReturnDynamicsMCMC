#include "ReturnDynamicsMCMC.h"

double ReturnDynamicsMCMC::stdNormalCDF(double z){
    return 0.5*(1+erf(z/sqrt(2)));
}

VectorXd ReturnDynamicsMCMC::compute_tp1(VectorXd& vec)
{
    VectorXd rtn = vec.tail(T);

    return rtn;
}

VectorXd ReturnDynamicsMCMC::compute_tp0(VectorXd& vec)
{
    VectorXd rtn = vec.head(T);

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

double ReturnDynamicsMCMC::compute_K(VectorXd& mu_vec, VectorXd& Y_vec)
{
    // Vectors
    VectorXd C_tp1=compute_disp_vec(Y_vec, mu_vec);
    C_tp1 = C_tp1 / sig_y;

    VectorXd mu_t=compute_tp0(mu_vec);
    VectorXd Delta_t =compute_disp_const(mu_t, E_mu);
    Delta_t = Delta_t / sig_mu;

    VectorXd mu_tp1=compute_tp1(mu_vec);
    VectorXd Delta_tp1 = compute_disp_const(mu_tp1, E_mu);
    Delta_tp1 = Delta_tp1 / sig_mu;

    VectorXd delta= Delta_tp1 - beta * Delta_t;

    double psi = 1.0 / (1.0 - rho*rho);

    double rtn = rho*psi* C_tp1.dot(delta);

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
    C_tp1 = C_tp1 / sig_y;

    VectorXd D_tp1=compute_disp_vec(mu_vec, mu_vec, beta);
    D_tp1 = D_tp1 / sig_mu;

    // Costants
    double B = (1.0 - beta) / sig_mu;
    double B_sq =B*B;
    double psi = 1.0 / (1.0 - rho*rho);
    double E_sq = prior_E*prior_E;

    // W,S
    double W = 1.0 / E_sq + psi*B_sq*double(T);
    double S = prior_e / E_sq - rho * psi * B *C_tp1.sum() + psi * B * D_tp1.sum();

    // Random Generator
    if (DEBUG)
        cout << "... dist: ~ N( " << S/W << ", " << sqrt(1.0/W) << "^2 )" << endl;

    normal_distribution<double> dist (S/W, sqrt(1.0/W));
    double rtn=-1.0;

    while(rtn<0){
        rtn=dist(generator);
    }

    if (DEBUG)
        cout << "... sampled : " << rtn << endl;

    if (DEBUG)
        cout << endl;

    return rtn;
}

double ReturnDynamicsMCMC::posterior_beta(default_random_engine& generator, VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG){

    if (DEBUG)
        cout << "[Posterior] beta" << endl;
    if (DEBUG)
        cout << "... input [ E_mu=" << E_mu << ", beta = " << beta << ", rho = " <<
                rho << ", sig_y = " << sig_y << ", sig_mu = " << sig_mu << " ] " <<endl;

    // Vectors
    VectorXd C_tp1=compute_disp_vec(Y_vec, mu_vec);
    C_tp1 = C_tp1 / sig_y;

    VectorXd ones(T);
    ones.setOnes();
    VectorXd Delta_t = mu_current.head(T)-E_mu*ones;
    Delta_t = Delta_t / sig_mu;
    VectorXd Delta_tp1 = mu_current.tail(T)-E_mu*ones;
    Delta_tp1 = Delta_tp1 / sig_mu;

    // Constants
    double psi = 1.0 / (1.0 - rho*rho);
    double F_sq = prior_F*prior_F;

    // W,S
    double W = 1.0 / F_sq + psi*Delta_t.dot(Delta_t);
    double S = prior_f / F_sq - rho*psi*C_tp1.dot(Delta_t) + psi*Delta_tp1.dot(Delta_t);

    double n_mean=S/W;
    double n_var=1.0/W;

    if (DEBUG)
        cout << "... dist: ~ N( " << n_mean << ", " << sqrt(n_var) << "^2 )" << endl;

    double z_upper = (1.0 - n_mean) / sqrt(n_var);
    double z_lower = (-1.0 - n_mean) / sqrt(n_var); //!!!
    if (DEBUG)
        cout << "... z_lower: " << z_lower << ", z_upper: " << z_upper << endl;

    // Random Generator for beta
    double rtn = 1.0;

    uniform_real_distribution<double> dist (0, 1.0);
    double u = dist(generator);
    // Avoid the cases where
    // ... stdNormalCDF(z_upper)=stdNormalCDF(z_lower)=0 --> 0-0=0
    // ... stdNormalCDF(z_upper)=stdNormalCDF(z_lower)=1 --> 1-1=0
    if ( (stdNormalCDF(z_upper) > 0.0) && (stdNormalCDF(z_lower) < 1.0) ){
        double prob = (stdNormalCDF(z_upper)-stdNormalCDF(z_lower))*u+stdNormalCDF(z_lower);
        StandardNormalDistribution std_normal_dist;
        rtn=std_normal_dist.inv_cdf(prob);
    }else{
        rtn=(z_upper - z_lower)*u+z_lower;
        if (DEBUG)
            cout << "!!! WARNING: beta was picked with an extreme value... apply an approximation" << endl;
    }
    // Transformation "rtn=z" to sigma * z + mu
    rtn = n_mean + sqrt(n_var) * rtn;

    if (fabs(rtn) >= 1.0)
        throw runtime_error("ABS(beta) >= 1");

    if (DEBUG)
        cout << "... sampled : " << rtn << endl;

    if (DEBUG)
        cout << endl;

    return rtn;
}

double ReturnDynamicsMCMC::posterior_sig_y(default_random_engine& generator,
                                           VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG)
{
    if (DEBUG)
        cout << "[Posterior] sigma_y" << endl;
    if (DEBUG)
        cout << "... input [ E_mu=" << E_mu << ", beta = " << beta << ", rho = " <<
                rho << ", sig_y = " << sig_y << ", sig_mu = " << sig_mu << " ] " <<endl;

    double K = compute_K(mu_vec, Y_vec);
    if (DEBUG)
        cout << "... K :" << K << endl;

    // Vectors
    VectorXd C_tp1=compute_disp_vec(Y_vec, mu_vec);

    // Constants
    double psi = 1.0 / (1.0 - rho*rho);
    double A = 0.5*T + prior_alpha;
    double B = 1.0 / prior_beta + 0.5*psi*C_tp1.dot(C_tp1);
    B = 1.0 / B;

    // Random generation for sig_y*
    double sig_y_new = 0.0;
    while (sig_y_new<=0.0){ // avoid selecting 0 or negative --> InvGamma != 1 / 0 and > 0
        gamma_distribution<double> dist_sig_y_sq(A,B);
        sig_y_new = dist_sig_y_sq(generator);
    }
    sig_y_new = 1.0 / sqrt(sig_y_new);

    if (DEBUG)
        cout << "... New: " << sig_y_new << ", Current: " << sig_y << endl;

    // Compute L(sigma) and apply the criterion that determines updating sig_y* or not.
    uniform_real_distribution<double> dist_unif (0.0, 1.0);

    double L_current = exp (K / sig_y);
    double L_new = exp (K / sig_y_new);
    double u=dist_unif(generator);

    // Otherwise,
    //double log_u=log(dist_unif(generator));
    //double sig_RHS = K*(sig_y - sig_y_new) / sig_y / sig_y_new;

    if (L_current*u < L_new)
    {
        if (DEBUG){
            cout << "... Update sig_y = " << sig_y_new << endl;
            cout << endl;

        }
        return sig_y_new;
    }
    else{
        if (DEBUG){
            cout << "... Keep sig_y = " << sig_y << endl;
            cout << endl;
        }
        return sig_y;
    }
}

vector<double> ReturnDynamicsMCMC::posterior_phi_omega(default_random_engine& generator,
                                                       VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG)
{
    if (DEBUG)
        cout << "[Posterior] (sigma_mu, rho)" << endl;
    if (DEBUG)
        cout << "... input [ E_mu=" << E_mu << ", beta = " << beta << ", rho = " <<
                rho << ", sig_y = " << sig_y << ", sig_mu = " << sig_mu << " ] " <<endl;

    // Variable Conversion
    double phi = sig_mu*rho;
    double omega = sig_mu*sig_mu*(1.0 - rho*rho);
    if (DEBUG)
        cout << "... reparametrization...input [ phi=" << phi << ", omega=" << omega << " ]" << endl;

    // Vectors
//    VectorXd C_tp1=compute_disp_vec(Y_vec, mu_vec);
    VectorXd C_tp1 = Y_current.tail(T) - mu_current.head(T);
    VectorXd gamma= C_tp1 / sig_y;

//    VectorXd mu_t=compute_tp0(mu_vec);
//    VectorXd Delta_t = compute_disp_const(mu_t, E_mu);
//    VectorXd mu_tp1=compute_tp1(mu_vec);
//    VectorXd Delta_tp1 = compute_disp_const(mu_tp1, E_mu);
    VectorXd ones(T);
    ones.setOnes();
    VectorXd Delta_t = mu_current.head(T)-E_mu*ones;
    VectorXd Delta_tp1 = mu_current.tail(T)-E_mu*ones;
    VectorXd delta = Delta_tp1 - beta * Delta_t;

    // S, W
    double S = delta.dot(gamma);
    double W = gamma.dot(gamma) + 2.0;
    double A = 0.5*T+0.5+prior_a;
    double B = 0.5*delta.dot(delta)+1.0/prior_b-0.5*S*S/W;
    B= 1.0 / B;

    // Random generation for omega
    double omega_new = 0.0;
    while (omega_new<=0.0){ // omega > 0
        gamma_distribution<double> dist_omega (A, B);
        omega_new = dist_omega(generator);
    }
    omega_new = 1.0 / omega_new;

    // Normal posterior for phi|omega
    double n_mean=S/W;
    double n_var=omega_new/W;
    if (DEBUG)
        cout << "... dist: ~ N( " << n_mean << ", " << sqrt(n_var) << "^2 )" << endl;

    // Random generation for phi|omega
    default_random_engine phi_generator;
    double phi_new = 0.0;
    while (phi_new==0.0){ // phi !=0
        normal_distribution<double> dist_phi(n_mean, sqrt(n_var));
        phi_new = dist_phi(phi_generator);
    }
    if (DEBUG)
        cout << "... reparametrization...update [ phi=" << phi_new << ", omega=" << omega_new << " ]" << endl;

    // Reparameterization
    double sig_mu_new = sqrt( omega_new + phi_new * phi_new);
    double rho_new = phi_new /sqrt( omega_new + phi_new * phi_new);

    if (DEBUG)
        cout << "... update [ sig=" << sig_mu_new << ", rho=" << rho_new << " ]" << endl;

    sig_mu = sig_mu_new;
    rho = rho_new;

    vector<double> rtn;
    rtn.clear();
    rtn.push_back(sig_mu);
    rtn.push_back(rho);

    if (DEBUG)
        cout << endl;

    return rtn;
}

VectorXd ReturnDynamicsMCMC::posterior_mu(default_random_engine& generator, VectorXd& mu_vec, VectorXd& Y_vec, bool DEBUG)
{
    assert(T >=2);

    VectorXd check_mu=mu_vec-mu_current;
    VectorXd check_y=Y_vec-Y_current;

    if (DEBUG)
        cout << "[Posterior] mu" << endl;
    if (DEBUG)
        cout << "... input [ E_mu=" << E_mu << ", beta = " << beta << ", rho = " <<
                rho << ", sig_y = " << sig_y << ", sig_mu = " << sig_mu << " ] " <<endl;

    VectorXd rtn(T+1);

    // Parameter
    double sig_y_sq=sig_y*sig_y;
    double sig_mu_sq=sig_mu*sig_mu;
    double sig_y_mu=sig_y*sig_mu;
    double beta_sq =beta*beta;

    // mu_0
    double y_1= Y_current[1];
    double mu_1 = mu_current[1];
    double W_0 = 1.0 / sig_y_sq - 2.0 * rho * beta / sig_y_mu + beta_sq/sig_mu_sq;
    double S_0 = y_1/sig_y_sq +
            (-rho*beta*y_1 - rho*mu_1 + rho*E_mu - rho*beta*E_mu) / sig_y_mu +
            (beta_sq*E_mu + beta*mu_1 - E_mu*beta ) / sig_mu_sq ;

    normal_distribution<double> dist_0(S_0/W_0, sqrt((1.0-rho*rho)/W_0));
    double mu_0_gen = dist_0(generator);
    rtn[0] = mu_0_gen;
    if (DEBUG){
        cout << "... mu_0 ~ N( " << S_0/W_0 << ", " << (1.0-rho*rho)/W_0 << " )" << endl;
        cout << "... mu_0 = " << mu_0_gen << " <= "<< mu_current[0] << endl;
    }

    // mu_T
    double y_T= Y_current[T];
    double mu_T1 = mu_current[T-1];
    double W_T = 1.0/sig_mu_sq;
    double S_T = (rho * y_T- rho * mu_T1) / sig_y_mu +
            (E_mu + beta *mu_T1 - beta * E_mu) / sig_mu_sq;

    normal_distribution<double> dist_T(S_T/W_T, sqrt((1.0-rho*rho)/W_T));
    double mu_T_gen = dist_T(generator);
    rtn[T] = mu_T_gen;
    if (DEBUG){
        cout << "... mu_T ~ N( " << S_T/W_T << ", " << (1.0-rho*rho)/W_T << " )" << endl;
        cout << "... mu_T = " << mu_T_gen << " <= "<< mu_current[T] << endl;
    }

    // mu_t
    double W_t=1.0 / sig_y_sq - 2.0 * rho * beta / sig_y_mu + (beta_sq + 1.0) / sig_mu_sq;
    for (int t=1; t<T; t++)
    {
        double y_tp1 = Y_current[t+1];
        double y_t = Y_current[t];
        double mu_tp1 = mu_current[t+1];
        double mu_tm1 = mu_current[t-1];
        double S_t = y_tp1 / sig_y_sq +
                (-rho* beta * y_tp1- rho * mu_tp1 + rho * E_mu - rho * beta * E_mu) / sig_y_mu +
                (beta_sq * E_mu + mu_tp1 * beta - beta * E_mu) / sig_mu_sq +
                (rho * y_t - rho * mu_tm1) / sig_y_mu +
                (E_mu + beta * mu_tm1 - beta * E_mu)/ sig_mu_sq;
        normal_distribution<double> dist_t(S_t/W_t, sqrt((1.0-rho*rho)/W_t));
        double mu_t_gen = dist_t(generator);
        rtn[t]=mu_t_gen;
//        if (DEBUG){
//            cout << "... mu_"<< t <<"~ N( " << S_t/W_t << ", " << (1.0-rho*rho)/W_t << " )" << endl;
//            cout << "... mu_"<< t <<"= " << mu_t_gen << " <= "<< mu_vec[t] << endl;
//        }
    }

    mu_current=rtn;

    if (DEBUG)
        cout << "... Mean of {mu_t}_t={0,...,T} : " << rtn.mean() << endl;

    if (DEBUG)
        cout << endl;

    return rtn;
}

