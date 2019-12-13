# ReturnDynamicsMCMC
### : Implementing MCMC algorithms to find parameters in Return Dynamics

## [ How to run ]
- **"setup.txt"** in the main directory is the setting file for how to implement MCMC 
- It contains the options that choose either performing simulation or using observed data,
  MCMC initial guess for paramters and hyperparameters
- Using the shell script below, you may implement MCMC and have figures that summarize final results

```
:~$ ./run.sh
```

## [ Summary for the parameters from MCMC implementations ]

![table](./plot_script/posterior_output.png)

## [ Comparison between the true and estimated mu ]
### The last mu

![mu_last](./plot_script/mu_last.png)

### The last 100 iteration-averaged mu

![mu_mean](./plot_script/mu_mean.png)

## [ Chainplot for each parameter in the total iterations ]

![E_mu](./plot_script/E_mu_total.png)

![beta](./plot_script/beta_total.png)

![sigma_y](./plot_script/sigma_y_total.png)

![sigma_mu](./plot_script/sigma_mu_total.png)

![rho](./plot_script/rho_total.png)

## [ Chainplot for each parameter in early stage of iterations ]

![E_mu_f](./plot_script/E_mu_first.png)

![beta_f](./plot_script/beta_first.png)

![sigma_y_f](./plot_script/sigma_y_first.png)

![sigma_mu_f](./plot_script/sigma_mu_first.png)

![rho_f](./plot_script/rho_first.png)


## [ Prerequisites for this code in Ubuntu ]
### Core library needed for C++
- **Eigen3**
```
:~$ sudo apt-get install -y apt-get install libeigen3-dev
```
### Graphic tool in Python script
- **Python3**
```
sudo apt-get install -y python3 python3-pip python3-dev python3-env
```
- **Numpy & Scipy**
```
sudo pip install numpy scipy
```
- **Matplotlib**
```
sudo pip install matplotlib
```
