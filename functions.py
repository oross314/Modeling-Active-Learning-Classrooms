import numpy as np
import matplotlib.pyplot
import pandas as pd
from sklearn.model_selection import train_test_split
from statsmodels.regression.linear_model import  WLS



def effect_size(mean1, mean2, sd1, sd2):
    """Calculate Cohen's d effect size between two groups."""
    pooled_sd = np.sqrt((sd1**2 + sd2**2) / 2)
    return (mean2 - mean1) / pooled_sd



def add_cross_terms(data):
    ## this functions adds polynomial and cross terms to the data
        
    pts, features = data.shape
    #
    cross_factors = (data[:, :features].reshape(pts, 1, features) *data[:, :features].reshape(pts, features, 1))#.reshape(pts, features * features)
    #keep only upper triangle

    #second order terms
    cross_factors = cross_factors[:, np.triu_indices(features, k=1)[0], np.triu_indices(features, k=1)[1]]
    cross_factors = np.hstack((cross_factors, data**2))

    #third order terms
    cross_21 = (data.reshape(pts, 1, features)**2 * data.reshape(pts, features, 1))
    cross_12 = (data.reshape(pts, 1, features) * data.reshape(pts, features, 1)**2)


    cross_factors = np.hstack((cross_factors, cross_21[:, np.triu_indices(features, k=1)[0], np.triu_indices(features, k=1)[1]]))
    cross_factors = np.hstack((cross_factors, cross_12[:, np.triu_indices(features, k=1)[0], np.triu_indices(features, k=1)[1]]))


    data = np.hstack((data, cross_factors))

    

    data = np.hstack((data, data[:, :1]*data[:, 1:2]* data[:, 2:3]))  # interaction term for first three features
    if features >= 4:
        data = np.hstack((data, data[:, :1]*data[:, 1:2]* data[:, 3:4]))  # interaction term for first second and fourth features
        data = np.hstack((data, data[:, :1]*data[:, 2:3]* data[:, 3:4]))  # interaction term for first third and fourth features
        data = np.hstack((data, data[:, 1:2] * data[:, 2:3] * data[:, 3:4]))  # interaction term for second third and fourth features
    
    if features >= 5:
        data = np.hstack((data, data[:, :1] * data[:, 1:2] * data[:, 4:5]))  # interaction term for first second and fifth features
        data = np.hstack((data, data[:, :1] * data[:, 2:3] * data[:, 4:5]))  # interaction term for first third and fifth features
        data = np.hstack((data, data[:, :1] * data[:, 3:4] * data[:, 4:5]))  # interaction term for first fourth and fifth features
        data = np.hstack((data, data[:, 1:2] * data[:,  2:3] * data[:, 4:5]))  # interaction term for second third and fifth features
        data = np.hstack((data, data[:, 1:2] * data[:, 3:4] * data[:, 4:5]))  # interaction term for second fourth and fifth features
        data = np.hstack((data, data[:, 2:3] * data[:, 3:4] * data[:, 4:5]))  # interaction term for third fourth and fifth features
    data = np.hstack((data, data[:, :4]**3))  #third order terms

    
    #fourth order terms
    data = np.hstack((data, data[:, :4]**4))

    data = np.hstack((data, np.ones((data.shape[0], 1))))  # add bias term
    return data


def aic_correction(aic, data):
    #small sample size correction for AIC

    n = data.shape[0]  # number of observations
    k = data.shape[1]  # number of parameters
    return aic + (2 * k * (k + 1)) / (n - k - 1)


def feature_selection(data, target, its=10000, prt = True, num_choose = 10, weights = False):
    if np.all(weights == False):
        weights = np.ones(data.shape[0])

    regr = WLS(target, data, weights=weights).fit()
    true_base_AIC = regr.bic
    AICs = np.zeros(its)
    keepers = np.zeros((its, data.shape[1]), dtype=bool)
    keepers[:, -1] = 1  #always keep the bias term

    for j in range(its):
        base_AIC = true_base_AIC * 1
        
        #randomly shuffle order of variables, skipp bias
        shuffle = np.append(np.array([data.shape[1] - 1]), np.random.permutation(data.shape[1] - 1))
        pruned_data = data[:, shuffle]

        i = 1
        keeping = []
        while i < pruned_data.shape[1]:

            #remove one feature
            cut_data = np.delete(pruned_data, i, axis=1)

            regr = WLS(target, cut_data, weights=weights).fit()
            new_AIC = regr.bic

            #if AIC doesn't increase by more than 2, remove the variable and update the base AIC
            ################################################3
            if new_AIC   < base_AIC  : #changed tolerance to 2
            #####################################################3
                pruned_data = np.delete(pruned_data, i, axis=1)
                base_AIC = new_AIC
            #otherwise keep the variable and move to the next one
            else:
                keepers[j, shuffle[i + (data.shape[1] - pruned_data.shape[1])]] = 1
                keeping.append(i + (data.shape[1] - pruned_data.shape[1]))
                i += 1
            
        AICs[j] = base_AIC

        if not j % 1000 and prt:
            print(f'Iteration {j + 1} finished with AIC {base_AIC.round(2)}')
    good_fits = np.argsort(AICs)[:num_choose]
    return keepers, AICs, good_fits

def bootstrap_LOO(data, target, keepers, good_fits, bootstrap_n, leave_out = False, weights = 1.0):

    if leave_out is False:
        #if a specific set isn't specificied to do LOOs on, do all of them
        leave_out = np.arange(data.shape[0])
    #if leave out is int make list
    if isinstance(leave_out, int):
        leave_out = [leave_out,]
    LOO_predictions = np.zeros((len(good_fits), bootstrap_n, len(leave_out)))

    # Store LOO weights for each fit and each data point
    LOO_weights = np.zeros((len(good_fits), bootstrap_n,  data.shape[1], len(leave_out),))
    #perform leave-one-out cross-validation for each of the best fits
    for i in range(len(leave_out)):
        for k in range(len(good_fits)):
            leave_out_data = np.delete(data[:, keepers[good_fits][k]], leave_out[i], axis=0)
            leave_out_truth = np.delete(target, leave_out[i], axis=0)
            if not np.all(weights == 1):
                leave_out_weights = np.delete(weights, leave_out[i], axis=0)
            else: leave_out_weights = np.ones(leave_out_data.shape[0])
            for l in range(bootstrap_n):
                #select .85 of the remaining data for training
                cut_data, _, cut_truth, _, cut_weights, _ = train_test_split(leave_out_data, leave_out_truth, leave_out_weights, test_size=.15)
                regr = WLS(cut_truth, cut_data, weights=cut_weights).fit()
                params = regr.params
                LOO_predictions[k, l, i] = data[leave_out[i], keepers[good_fits][k]] @ params
                # Store weights for this LOO fit
                LOO_weights[k, l, keepers[good_fits][k], i] = params


    return np.squeeze(LOO_predictions), np.squeeze(LOO_weights)

def feature_boot_LOO(data, target, ES_variance, bootstrap_n = 10, n_fits = 10, ):
#create empty arrays to store weights and predictions for each bootstrap sample, each best fit, and each data point
    feature_selection_its = int(data.shape[0]*5)
    LOO_predictions = np.zeros(( n_fits, bootstrap_n, data.shape[0]))
    LOO_weights = np.zeros(( n_fits, bootstrap_n, data.shape[0], data.shape[1]))
    good_fit_AICs = np.zeros(( n_fits, data.shape[0] ))

    #perform leave-one-out cross-validation for each of the best fits
    for k in range(data.shape[0]): 
        print(f'LOO iteration {k} of {data.shape[0]}')
        cut_data = np.delete(data, k, axis=0)
        cut_truth = np.delete(target, k, axis=0)
        cut_weights = np.delete(1/ES_variance, k, axis=0)
        keepers, AICs, good_fits = feature_selection(cut_data, cut_truth, its=feature_selection_its, prt = False, num_choose = n_fits, weights = cut_weights)
        #print(AICs[good_fits])

        good_fit_AICs[:, k] = AICs[good_fits]
        LOO_predictions[:, :, k], LOO_weights[:, :, k, :] = bootstrap_LOO(data, target, keepers, good_fits, 
                                                                        bootstrap_n=bootstrap_n, leave_out = k, weights = 1/ES_variance   )
    return LOO_predictions, LOO_weights, good_fit_AICs