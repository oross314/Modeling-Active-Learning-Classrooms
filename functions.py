import numpy as np
import matplotlib.pyplot
import pandas as pd
from scipy.spatial import ConvexHull


def effect_size(mean1, mean2, sd1, sd2):
    """Calculate Cohen's d effect size between two groups."""
    pooled_sd = np.sqrt((sd1**2 + sd2**2) / 2)
    return (mean2 - mean1) / pooled_sd


powers = np.array([ 2, 3,  4,])
def add_cross_terms(data):
    ## this functions adds polynomial and cross terms to the data
        
    pts, features = data.shape
    #
    cross_factors = (data[:, :features].reshape(pts, 1, features) *data[:, :features].reshape(pts, features, 1))#.reshape(pts, features * features)
    #keep only upper triangle

    cross_factors = cross_factors[:, np.triu_indices(features, k=1)[0], np.triu_indices(features, k=1)[1]]
    cross_21 = (data.reshape(pts, 1, features)**2 * data.reshape(pts, features, 1))
    cross_12 = (data.reshape(pts, 1, features) * data.reshape(pts, features, 1)**2)


    cross_factors = np.hstack((cross_factors, cross_21[:, np.triu_indices(features, k=1)[0], np.triu_indices(features, k=1)[1]]))
    cross_factors = np.hstack((cross_factors, cross_12[:, np.triu_indices(features, k=1)[0], np.triu_indices(features, k=1)[1]]))



    
    copy_data = data.copy()
    for i in powers:
        copy_data = np.hstack((copy_data, data**i))
    data = copy_data


    data = np.hstack((data, cross_factors))
    data = np.hstack((data, data[:, :1]*data[:, 1:2]* data[:, 2:3]))  # interaction term for first three features
    
    data = np.hstack((data, data[:, :1]*data[:, 1:2]* data[:, 3:4]))  # interaction term for first three features
    data = np.hstack((data, data[:, :1]*data[:, 2:3]* data[:, 3:4]))  # interaction term for first three features
    data = np.hstack((data, data[:, 1:2] * data[:, 2:3] * data[:, 3:4]))  # interaction term for second third and fourth features

    data = np.hstack((data, np.ones((data.shape[0], 1))))  # add bias term
    return data

def aic_correction(aic, data):
    #small sample size correction for AIC

    n = data.shape[0]  # number of observations
    k = data.shape[1]  # number of parameters
    return aic + (2 * k * (k + 1)) / (n - k - 1)
