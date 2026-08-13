import numpy as np
import matplotlib.pyplot
import pandas as pd
from sklearn.model_selection import train_test_split
from statsmodels.regression.linear_model import OLS



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
    
    #######################
    if features >= 5:
        data = np.hstack((data, data[:, :1] * data[:, 1:2] * data[:, 4:5]))  # interaction term for first second and fifth features
        data = np.hstack((data, data[:, :1] * data[:, 2:3] * data[:, 4:5]))  # interaction term for first third and fifth features
        data = np.hstack((data, data[:, :1] * data[:, 3:4] * data[:, 4:5]))  # interaction term for first fourth and fifth features
        data = np.hstack((data, data[:, 1:2] * data[:,  2:3] * data[:, 4:5]))  # interaction term for second third and fifth features
        data = np.hstack((data, data[:, 1:2] * data[:, 3:4] * data[:, 4:5]))  # interaction term for second fourth and fifth features
        data = np.hstack((data, data[:, 2:3] * data[:, 3:4] * data[:, 4:5]))  # interaction term for third fourth and fifth features
        ########################
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


def feature_selection(data, target, its=10000, n =1, prt = True, num_choose = 5):
    
    #Establish baseline AIC with all features
    AIC = 0
    for i in range(10):
        X_train, X_test, y_train, y_test = train_test_split(data, target, test_size= int(data.shape[0] * 0.05), )

        regr = OLS(y_train, X_train).fit()
        AIC += aic_correction(regr.aic, X_train)
    true_base_AIC = AIC / 10

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
            #initialize AIC
            AIC = 0
            #run regression n times with n different subsets of the training data
            for k in range(n):
                
                cut_data = np.delete(pruned_data, i, axis=1)

                delete_set = np.random.choice(cut_data.shape[0], size=int(cut_data.shape[0] * 0.05), replace=False)
                cut_data = np.delete(cut_data, delete_set, axis=0)
                cut_truth = np.delete(target, delete_set, axis=0)
                
                regr = OLS(cut_truth, cut_data).fit()
                AIC += aic_correction(regr.aic, cut_data)
            #take average AIC over n runs
            av_AIC = AIC / n
            #if AIC improves by removing variable, permanently remove it
            if av_AIC - 2 < base_AIC  : #changed tolerance to 2

                pruned_data = np.delete(pruned_data, i, axis=1)
                base_AIC = av_AIC
            #otherwise keep the variable and move to the next one
            else:
                keepers[j, shuffle[i + (data.shape[1] - pruned_data.shape[1])]] = 1
                keeping.append(i + (data.shape[1] - pruned_data.shape[1]))
                i += 1
            
        AICs[j] = base_AIC
        if not j % 1000 and prt:
            print(f'Iteration {j + 1} finished with AIC {base_AIC.round(2)}')
    good_fits = np.argsort(AICs)[:num_choose]
    #np.save('weights/keepers.npy', keepers)
    #np.save('weights/AICs.npy', AICs)
    return keepers, AICs, good_fits

def bootstrap_LOO(data, target, keepers, good_fits, bootstrap_n, leave_out = False):

    if leave_out is False:
        leave_out = np.arange(data.shape[0])
    LOO_predictions = np.zeros((len(good_fits), bootstrap_n, data.shape[0]))

    # Store LOO weights for each fit and each data point
    LOO_weights = np.zeros((len(good_fits), bootstrap_n, data.shape[0], data.shape[1]))
    #perform leave-one-out cross-validation for each of the best fits
    for k in range(len(good_fits)):
        for l in range(bootstrap_n):
            for i in leave_out:
                cut_data = np.delete(data[:, keepers[good_fits][k]], i, axis=0)
                cut_truth = np.delete(target, i, axis=0)
                #select .85 of the remaining data for training
                cut_data, _, cut_truth, _ = train_test_split(cut_data, cut_truth, test_size=.1)
                regr = OLS(cut_truth, cut_data).fit()
                params = regr.params
                LOO_predictions[k, l, i] = data[i, keepers[good_fits][k]] @ params
                # Store weights for this LOO fit
                LOO_weights[k, l, i, keepers[good_fits][k]] = params

    #average over feature sets
    LOO_predictions = np.mean(LOO_predictions, axis=0)

    #find mean and standard deviation of predictions for each data point across all LOO fits
    LOO_predictions_mean = np.mean(LOO_predictions, axis=0)
    LOO_predictions_std = np.std(LOO_predictions, axis=0)
    return LOO_predictions_mean, LOO_predictions_std, LOO_weights