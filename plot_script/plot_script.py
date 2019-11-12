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


# ## Plotting

# In[2]:


## Log file path
log_file_dir='../log/MCMC_log_obs_1.csv'

## Read CSV
params=["E_mu", "beta", "sigma_y", "sigma_mu", "rho"]
data = pd.read_csv(log_file_dir,skiprows=1, delimiter=",", names=params)
data.drop(data.tail(1).index,inplace=True)

n=len(data)

E_mu=np.array(data["E_mu"]).astype(float)
beta=np.array(data["beta"]).astype(float)
sigma_y=np.array(data["sigma_y"]).astype(float)
sigma_mu=np.array(data["sigma_mu"]).astype(float)
rho=np.array(data["rho"]).astype(float)
t=np.arange(0,n,1)

fig1=plt.figure(1, figsize=[12,4], dpi=240)
plt.plot(t, E_mu)
plt.title('E_mu')
#plt.show()
fig1.savefig("E_mu.png", dpi=fig1.dpi)

fig2=plt.figure(2, figsize=[12,4], dpi=240)
plt.plot(t, beta)
plt.title('beta')
#plt.show()
fig2.savefig("beta.png", dpi=fig2.dpi)

fig3=plt.figure(3, figsize=[12,4], dpi=240)
plt.plot(t, sigma_y)
plt.title('sigma_y')
#plt.show()
fig3.savefig("sigma_y.png", dpi=fig3.dpi)

fig4=plt.figure(4, figsize=[12,4], dpi=240)
plt.plot(t, sigma_mu)
plt.title('sigma_mu')
#plt.show()
fig4.savefig("sigma_mu.png", dpi=fig4.dpi)

fig5=plt.figure(5, figsize=[12,4], dpi=240)
plt.plot(t, rho)
plt.title('rho')
#plt.show()
fig5.savefig("rho.png", dpi=fig5.dpi)

