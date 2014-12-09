import numpy
from numpy import array
import dadi
import models


####DATA IMPORT
dd = dadi.Misc.make_data_dict('hunsteiniSnpFile.dadi')

####AFS generation
fs = dadi.Spectrum.from_data_dict(dd, ['1','2'],[19,17])
ns = fs.sample_sizes

####SET UP THE GRID
pts_l = [60,70,80]

####GET THE MODEL
func = models.modelFour

####MODEL PARAMETERS
# s, nu1, nu2, T, m12
params = array([0.5, 0.4, 0.4, 0.05, 1.2, 1.2])
upper_bound = [1, 100, 100, 10, 20, 20]
lower_bound = [1e-2, 1e-2, 0, 0, 0, 0]

####EXTRAPOLATING FUNCTION
func_ex = dadi.Numerics.make_extrap_log_func(func)

####MODEL AFS
model = func_ex(params, ns, pts_l)

####LIKELYHOOD OF DATA (given model)
ll_model = dadi.Inference.ll_multinom(model,fs)
print 'Model log-likelihood:', ll_model

####OPTIMAL THETA
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

####PARAMETER PERTURBATION
p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound)
popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        verbose = len(params),
        maxiter = 3)

print 'Optimized parameters', repr(popt)
model = func_ex(popt, ns, pts_l)
ll_opt = dadi.Inference.ll_multinom(model,fs)
print 'Optimized log-likelihood:', ll_opt

####PLOT A MODEL => DATA COMPARISON
import pylab
pylab.figure()
dadi.Plotting.plot_2d_comp_multinom(model, fs, vmin = 1, resid_range = 3,
        pop_ids = ('New Ireland','New Hanover'))

# ####MS STUFF
# mscore = ourDemographicModel.prior_onegrow_mig_mscore(params)
# # I find that it's most efficient to simulate with theta=1 and then scale up.
# mscommand = dadi.Misc.ms_command(1., ns, mscore, int(1e6))
# ## We use Python's os module to call this command from within the script.
# ## If you have ms installed, uncomment these lines to see the results.
# import os
# os.system('%s > test.msout' % mscommand)
# msdata = dadi.Spectrum.from_ms_file('test.msout')
# pylab.figure()
# dadi.Plotting.plot_2d_comp_multinom(model, theta*msdata, vmin=1,
#                                     pop_ids=('YRI','CEU'))
# pylab.show()








