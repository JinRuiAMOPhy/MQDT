import numpy as np
import matplotlib.pyplot as plt

dfde = np.loadtxt("dfde2-plot.out", dtype = float, unpack = False)

tau = np.loadtxt('tau2-plot.out', dtype = float, unpack = False)

theo = np.loadtxt('theo-plot.out', dtype = float, unpack = False)

nu = np.loadtxt('nux_nuy.out', dtype = float, unpack = False)

fig, axs = plt.subplots(2,1)

axs[0].plot(nu[:,0], nu[:,1], color = 'gray', alpha = 0.5)
for i in range(1,6+1) :
 line = axs[0].plot(tau[:,0], tau[:,i], label = '{}'.format(i))
 color = axs[0].get_lines()[-1].get_c() 
 axs[0].scatter(theo[:,0], theo[:,i], c= color)
 axs[1].plot(tau[:,0], tau[:,i], c= color)
axs[0].legend(loc = 'upper right', ncol = 2,
              title ='Channel index')
axs[0].set_ylabel(r'$-\nu_{IP1}$ (mod 1)')
#axs[1].set_ylabel(r'$\frac{df}{dE}$')
axs[1].set_ylabel(r'Oscillator strength $f_n$')
axs[1].set_xlabel(r'$\nu_{IP2}$')
fig.savefig('result.pdf', format = 'pdf')
