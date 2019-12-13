#!/usr/bin/env python
# coding: utf-8

# # Plot the MCMC results
# 
# ## Import libraries

# In[1]:


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ## Table

# In[2]:


posterior_file_dir = '../posterior_output.csv'
setup_file_dir = '../setup.txt'


# ## Chain plots

# In[3]:


def plotChainplot(data, string_post, burnin):
    
    file_tag=""
    if(burnin):
        file_tag="_burnin"
    else:
        file_tag="_first"
        
    string_pre = "Chainplots of [ "
    
    n=len(data)    
    E_mu=np.array(data["E_mu"]).astype(float)
    beta=np.array(data["beta"]).astype(float)
    sigma_y=np.array(data["sigma_y"]).astype(float)
    sigma_mu=np.array(data["sigma_mu"]).astype(float)
    rho=np.array(data["rho"]).astype(float)
    t=np.arange(0,n,1)

    fig1=plt.figure(1, figsize=[12,4], dpi=240)
    plt.plot(t, E_mu)
    plt.plot(0,E_mu_init,"r*", markersize=12)
    plt.title(string_pre+'E_mu'+string_post)
    plt.show()
    fig1.savefig("E_mu"+file_tag+".png", dpi=fig1.dpi)

    fig2=plt.figure(2, figsize=[12,4], dpi=240)
    plt.plot(t, beta)
    plt.plot(0,beta_init,"r*", markersize=12)
    plt.title(string_pre+'beta'+string_post)
    plt.show()
    fig2.savefig("beta"+file_tag+".png", dpi=fig2.dpi)

    fig3=plt.figure(3, figsize=[12,4], dpi=240)
    plt.plot(t, sigma_y)
    plt.plot(0,sigma_y_init,"r*", markersize=12)
    plt.title(string_pre+'sigma_y'+string_post)
    plt.show()
    fig3.savefig("sigma_y"+file_tag+".png", dpi=fig3.dpi)

    fig4=plt.figure(4, figsize=[12,4], dpi=240)
    plt.plot(t, sigma_mu)
    plt.plot(0,sigma_mu_init,"r*", markersize=12)
    plt.title(string_pre+'sigma_mu'+string_post)
    plt.show()
    fig4.savefig("sigma_mu"+file_tag+".png", dpi=fig4.dpi)

    fig5=plt.figure(5, figsize=[12,4], dpi=240)
    plt.plot(t, rho)
    plt.plot(0,rho_init,"r*", markersize=12)
    plt.title(string_pre+'rho'+string_post)
    plt.show()
    fig5.savefig("rho"+file_tag+".png", dpi=fig5.dpi)


## Log file path
log_file_dir='../log/MCMC_log_simul_1.csv'
first_period=300
burnin_period=10000

## Read CSV
params=["E_mu", "beta", "sigma_y", "sigma_mu", "rho"]
data = pd.read_csv(log_file_dir,skiprows=1, delimiter=",", names=params)
data.drop(data.tail(1).index,inplace=True)

n=len(data)
E_mu_init=data.astype(float).at[0,"E_mu"]
beta_init=data.astype(float).at[0,"beta"]
sigma_y_init=data.astype(float).at[0,"sigma_y"]
sigma_mu_init=data.astype(float).at[0,"sigma_mu"]
rho_init=data.astype(float).at[0,"rho"]

## First period
data_first=data[1:first_period+1]
string_first = " ] for the first " + str(first_period)+ " iterations"
plotChainplot(data_first, string_first, False)

## After burn-in period
data_burnin=data.tail(n-burnin_period)
string_burnin = " ] after the burn-in period of " + str(burnin_period)
plotChainplot(data_burnin,  string_burnin, True)


# ## Latent Variable

# In[4]:


simul_file_dir ='../log/simul_path_1.csv'
latent_file_dir = '../log/MCMC_log_latent_1.csv'

simul_column=["Y", "mu"]
latent_column=["last","mean"]
data_simul = pd.read_csv(simul_file_dir,skiprows=2, delimiter=",", names=simul_column)
data_latent = pd.read_csv(latent_file_dir,skiprows=1,names=latent_column)

mu_true=np.array(data_simul["mu"]).astype(float)
mu_last=np.array(data_latent["last"]).astype(float)
mu_mean=np.array(data_latent["mean"]).astype(float)
t_mu=np.arange(0,len(mu_true),1)

fig_mu1=plt.figure(1, figsize=[10,6], dpi=240)
plt.plot(t_mu,mu_true,"tomato",linestyle="--",label="true mu")
plt.plot(t_mu,mu_last,"blue",label="last mu")
plt.title("The estimated mu from the last iteration vs. the true mu")
plt.legend()
plt.show()
fig_mu1.savefig("mu_last.png", dpi=fig_mu1.dpi)

fig_mu2=plt.figure(2, figsize=[10,6], dpi=240)
plt.plot(t_mu,mu_true,"tomato",linestyle="--",label="true mu")
plt.plot(t_mu,mu_mean,"blue",label="averaged mu")
plt.title("The averaged mu from the last 100 iterations vs. the true mu")
plt.legend()
plt.show()
fig_mu2.savefig("mu_mean.png", dpi=fig_mu2.dpi)

