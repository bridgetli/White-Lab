"""
Created on Tue Oct  6 16:12:36 2020
@author: Bridget_PC
Practice using scrikit-learn for PCA
"""

import numpy as np
from sklearn.decomposition import PCA

X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
pca = PCA(n_components=2)
pca.fit(X)
print(pca.explained_variance_ratio_)
print(pca.singular_values_)
