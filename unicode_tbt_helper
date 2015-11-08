import pylab as pl
import numpy as np

uni_matrix = np.loadtxt('tabelaeta.txt')

A_x = 1e-5
omega_x = np.unique(uni_matrix[:,8])

index = np.argwhere((uni_matrix[:,8] == omega_x[0]) & (uni_matrix[:,10] == A_x))


index = index.T
index = index[0]

x=uni_matrix[index,:]


np.savetxt('omega_x.txt', x)



