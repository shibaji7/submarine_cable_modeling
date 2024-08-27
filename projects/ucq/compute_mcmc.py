######################################################
# Run MCMC for SCUBAS model for sampling and estimating
######################################################

import numpy as np

def physical_model(params, x):
    # Example: A simple linear model
    a, b = params
    return a * x + b

def likelihood(params, x, y_obs, y_err):
    y_model = physical_model(params, x)
    return -0.5 * np.sum(((y_obs - y_model) / y_err)**2)


import pymc3 as pm

with pm.Model() as model:
    a = pm.Normal('a', mu=0, sigma=1)
    b = pm.Normal('b', mu=0, sigma=1)


with model:
    # Set up the sampler
    step = pm.Metropolis()

    # Run the MCMC sampler
    trace = pm.sample(draws=1000, tune=500, step=step)