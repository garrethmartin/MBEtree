import MBE
import numpy as np
import matplotlib.pyplot as plt

XX = np.load('example_data.npy')
X = XX[:,0:3]
w = XX[:,3]

rho = MBE.modified_breiman_density(X, weights=w)

fig, ax = plt.subplots(figsize=(8,7));
h = ax.scatter(X[:,0], X[:,1], s=4, c=np.log10(rho), cmap='jet');
ax.set_xticks([])
ax.set_yticks([])
ax.axis('equal')
plt.colorbar(h, label=r'$\rm log_{10}\ n$')
plt.savefig('density_estimation.png', bbox_inches='tight')