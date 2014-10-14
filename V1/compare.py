#!/usr/bin/python

'''
comparative study between outputs
'''

import numpy as np
import matplotlib.pyplot as plt

# output files
out_alter_05 = np.genfromtxt('0.05_output_alternative.txt', skip_header = 1, usecols = (1, 2, 3))
out_alter_1 = np.genfromtxt('0.1_output_alternative.txt', skip_header = 1, usecols = (1, 2, 3))
out_orig_05 = np.genfromtxt('0.05_output_origin.txt', skip_header = 1, usecols = (1, 2, 3, 4))
out_orig_1 = np.genfromtxt('0.1_output_origin.txt', skip_header = 1, usecols = (1, 2, 3, 4))
out_gd = np.genfromtxt('output_geweke.txt', skip_header = 1, usecols = (1, 2, 3))

n_alter_05 = out_alter_05[0,0]
n_alter_1 = out_alter_1[0,0]
n_orig_05 = out_orig_05[0,0]
n_orig_1 = out_orig_1[0,0]
n_gd = out_gd[0,0]

## ess calculation
ess_alter_05 = out_alter_05[:,1:3]
ess_alter_1 = out_alter_1[:,1:3]
ess_orig_05 = out_orig_05[:,1:4]
ess_orig_1 = out_orig_1[:,1:4]
ess_gd = out_gd[:,1:3]

np.percentile(ess_alter_05[:,1], 50)
np.percentile(ess_alter_1[:,1], 50)
np.percentile(ess_orig_05[:,2], 50)
np.percentile(ess_orig_1[:,2], 50)
np.percentile(ess_gd[:,1], 50)

## ratio of mcse ##
mcse_alter_05 = np.genfromtxt('0.05_mcse_alternative.txt', skip_header = 1, usecols = 1)
mcse_alter_1 = np.genfromtxt('0.1_mcse_alternative.txt', skip_header = 1, usecols = 1)
mcse_orig_05 = np.genfromtxt('0.05_mcse_origin.txt', skip_header = 1, usecols = 1)
mcse_orig_1 = np.genfromtxt('0.1_mcse_origin.txt', skip_header = 1, usecols = 1)
mcse_gd = np.genfromtxt('mcse_geweke.txt', skip_header = 1, usecols = 1)

ratio_bm_1 = (mcse_alter_05*np.sqrt(n_alter_05))/(mcse_orig_05*np.sqrt(n_orig_05))
ratio_bm_2 = (mcse_alter_1*np.sqrt(n_alter_1))/(mcse_orig_1*np.sqrt(n_orig_1))

plt.subplots_adjust(hspace=.4)
plt.subplot(212)
plt.hist(ratio_bm_1, normed = 0, facecolor='red', alpha = 0.9, label = '$\epsilon = 0.05$')
plt.hist(ratio_bm_2, normed = 0, facecolor='blue', alpha=0.5, label = '$\epsilon = 0.1$')
plt.title('Ratios of variance estimators $\hat{\sigma}_{i} (n)$\'s')
plt.xlabel('ratio')
plt.ylabel('Frequency')
plt.legend()

## ratio mean ##
mean_alter_05 = np.genfromtxt('0.05_mean_alternative.txt', skip_header = 1, usecols = 1)
mean_alter_1 = np.genfromtxt('0.1_mean_alternative.txt', skip_header = 1, usecols = 1)
mean_orig_05 = np.genfromtxt('0.05_mean_origin.txt', skip_header = 1, usecols = 1)
mean_orig_1 = np.genfromtxt('0.1_mean_origin.txt', skip_header = 1, usecols = 1)

ratio_mean_1 = sorted(mean_alter_05/mean_orig_05)
ratio_mean_2 = sorted(mean_alter_1/mean_orig_1)

plt.subplot(211)
plt.hist(ratio_mean_1[int(0.05*len(ratio_mean_1)):int(0.95*len(ratio_mean_1))], normed = 0, facecolor='red', alpha = 0.9, label = '$\epsilon = 0.05$')
plt.hist(ratio_mean_2[int(0.05*len(ratio_mean_2)):int(0.95*len(ratio_mean_2))], normed = 0, facecolor='blue', alpha=0.5, label = '$\epsilon = 0.1$')
plt.title('Ratios of posterior mean estimators $Z_{i} (n)$\'s')
plt.xlabel('ratio')
plt.ylabel('Frequency')
plt.legend()

### ratio std #
sd_alter_05 = np.genfromtxt('0.05_sd_alternative.txt', skip_header = 1, usecols = 1)
sd_alter_1 = np.genfromtxt('0.1_sd_alternative.txt', skip_header = 1, usecols = 1)
sd_orig_05 = np.genfromtxt('0.05_sd_origin.txt', skip_header = 1, usecols = 1)
sd_orig_1 = np.genfromtxt('0.1_sd_origin.txt', skip_header = 1, usecols = 1)
sd_gd = np.genfromtxt('sd_geweke.txt', skip_header = 1, usecols = 1)
# 
# ratio_sd_1 = sorted(sd_alter_05/sd_orig_05)
# ratio_sd_2 = sorted(sd_alter_1/sd_orig_1)
# 
# plt.hist(ratio_sd_1, normed = 0, facecolor='red', alpha = 0.9, label = '$\epsilon = 0.05$')
# plt.hist(ratio_sd_2, normed = 0, facecolor='blue', alpha=0.5, label = '$\epsilon = 0.1$')
# plt.title('Ratios of posterior std estimator $\lambda_{n}$')
# plt.xlabel('ratio')
# plt.ylabel('Frequency')
# plt.legend()

## ratio w/sd ##
ratio_wsd_05= 2*1.96*mcse_alter_05/sd_orig_05
ratio_wsd_1= 2*1.96*mcse_alter_1/sd_orig_05
ratio_wsd_gd= 2*1.96*mcse_gd/sd_orig_05

plt.hist(ratio_wsd_gd, normed = 0, facecolor='red', alpha = 0.5, label = 'Geweke Diagnostic')
plt.hist(ratio_wsd_05, normed = 0, facecolor='green', alpha = 0.5, label = 'T($\epsilon$ = 0.05)')
plt.hist(ratio_wsd_1, normed = 0, facecolor='blue', alpha = 0.5, label = 'T($\epsilon$ = 0.1)')
plt.title('Comparison between GD and $T(\epsilon)$')
plt.xlabel('$w_{i} (n, \delta)/\hat{\lambda}_i (n)$')
plt.ylabel('Frequency')
plt.legend()