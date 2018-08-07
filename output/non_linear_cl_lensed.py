import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/michalychforever/Dropbox/Docs/science/class_public-master/output/non_linear_cl_lensed.dat', '/Users/michalychforever/Dropbox/Docs/science/class_public-master/output/linear_cl_lensed.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['non_linear_cl_lensed', 'linear_cl_lensed']

fig, ax = plt.subplots()
y_axis = [u'phiphi', u'Ephi']
tex_names = ['phiphi', 'Ephi']
x_axis = 'l'
ax.set_xlabel('$\ell$', fontsize=16)
plt.show()