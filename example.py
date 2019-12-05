import MBEtree
import numpy as np
import matplotlib.pyplot as plt

XX = np.load('example_data.npy')
X = XX[:,0:2]
w = XX[:,3]

rho = MBEtree.modified_breiman_density(X, weights=w)

fig, ax = plt.subplots(figsize=(8,7));
h = ax.scatter(X[:,0], X[:,1], c=np.log10(rho), s=np.log10(w/np.min(w)), cmap='jet');
ax.set_xticks([])
ax.set_yticks([])
ax.axis('equal')
plt.colorbar(h, label=r'$\rm log_{10}\ n$')
plt.savefig('density_estimation.png', bbox_inches='tight')
