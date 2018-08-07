import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/michalychforever/Dropbox/Docs/science/class_public-master/output/linear_cl_lensed.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['linear_cl_lensed']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'phiphi', u'Ephi']
tex_names = ['phiphi', 'Ephi']
x_axis = 'l'
ylim = []
xlim = []
ax.semilogx(curve[:, 0], curve[:, 5])
ax.semilogx(curve[:, 0], curve[:, 7])

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('$\ell$', fontsize=16)
plt.show()