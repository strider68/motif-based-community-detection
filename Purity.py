from sklearn import metrics
import numpy as np

#Compute purity
def Purity(y_true,y_pred):
    contigency_matrix=metrics.cluster.contingency_matrix(y_true,y_pred)
    return np.sum(np.amax(contigency_matrix,axis=0))/np.sum(contigency_matrix)
