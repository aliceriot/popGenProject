# Analyzing the hunsteini data using dadi
# need to already have run the 'makeSnpData' script!

import numpy
from numpy import array
import dadi


#in demographic_models we define the demographic functions we need
import ourDemographicModel


####DATA import
dd = dadi.Misc.make_data_dict('hunsteinAsnpFile.txt')


#which of these is right? dunno
#fs = dadi.Spectrum.from_data_dict(dd, ['1','2','3'],[20,20,20])

# I think we have 34 for '2' and 38 for '1'
newfs = dadi.Spectrum.from_data_dict(dd, ['1','2'],[38,34])
#newfs = dadi.Spectrum.from_data_dict(dd, ['1','2'],[19,17])
ns = newfs.sample_sizes



####SETTING STUFF UP
# I don't really get what these do...
pts_l = [40,50,60]

####EXAMPLE CODE w/ simple model
# The Demographics1D and Demographics2D modules contain a few simple models,
# mostly as examples. We could use one of those.
func = dadi.Demographics2D.split_mig
# ll for this model: -1136.61
params = array([1.792, 0.426, 0.309, 1.551])
upper_bound = [100, 100, 3, 20]


####REAL CODE w/ model from demographic models
func = ourDemographicModel.prior_onegrow_mig

# PARAMETERS! ok, so there are:
# recall that this model is for a model with growth, split, bottleneck
# in pop 2, recovery and then migration

#     nu1F: The ancestral population size after growth. (Its initial size is
#           defined to be 1.)
#     nu2B: The bottleneck size for pop2
#     nu2F: The final size for pop2
#     m: The scaled migration rate
#     Tp: The scaled time between ancestral population growth and the split.
#     T: The time between the split and present

#     n1,n2: Size of fs to generate.
#     pts: Number of points to use in grid for evaluation.
    
params = array([1.881, 0.0710, 1.845, 0.911, 0.355, 0.111])

#these help the process go more quickly? or something
upper_bound = [100, 100, 100, 100, 3, 3]

lower_bound = [1e-2, 1e-2, 1e-2, 0, 0, 0]

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)

# Calculate the model AFS, with the extrapolated model function from
# above
model = func_ex(params, ns, pts_l)

# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, newfs)
print 'Model log-likelihood:', ll_model
# The optimal value of theta given the model.
theta = dadi.Inference.optimal_sfs_scaling(model, newfs)

# Perturb our parameter array before optimization. This does so by taking each
# parameter a up to a factor of two up or down.
p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound)
# Do the optimization. By default we assume that theta is a free parameter,
# since it's trivial to find given the other parameters. If you want to fix
# theta, add a multinom=False to the call.
# (This is commented out by default, since it takes several minutes.)
# The maxiter argument restricts how long the optimizer will run. For production
# runs, you may want to set this value higher, to encourage better convergence.
popt = dadi.Inference.optimize_log(p0, newfs, func_ex, pts_l, 
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(params),
                                   maxiter=3)
print 'Optimized parameters', repr(popt)
model = func_ex(popt, ns, pts_l)
ll_opt = dadi.Inference.ll_multinom(model, newfs)
print 'Optimized log-likelihood:', ll_opt

# Plot a comparison of the resulting fs with the data.


####PLOTTING####

#here we generate a plot comparing the model AFS to the data AFS
import pylab
pylab.figure()
dadi.Plotting.plot_2d_comp_multinom(model, newfs, vmin=1, resid_range=3,
                                    pop_ids =('1','2'))
# This ensures that the figure pops up. It may be unecessary if you are using
# ipython.
#pylab.show()
pylab.savefig('YRI_CEU.png', dpi=50)

# Let's generate some data using ms, if you have it installed.
mscore = demographic_models.prior_onegrow_mig_mscore(params)
# I find that it's most efficient to simulate with theta=1 and then scale up.
mscommand = dadi.Misc.ms_command(1., ns, mscore, int(1e6))
## We use Python's os module to call this command from within the script.
## If you have ms installed, uncomment these lines to see the results.
#import os
#os.system('%s > test.msout' % mscommand)
#msdata = dadi.Spectrum.from_ms_file('test.msout')
#pylab.figure()
#dadi.Plotting.plot_2d_comp_multinom(model, theta*msdata, vmin=1,
#                                    pop_ids=('YRI','CEU'))
#pylab.show()

# Below here we compare uncertainty estimates from folded and unfolded spectra.
# Estimates are done using the hessian (Fischer Information Matrix).
# Due to linkage in the data, these are underestimates. It is still
# informative, however, to compare the two methods.

## These are the optimal parameters when the spectrum is folded. They can be
## found simply by passing fold=True to the above call to optimize_log.
#pfold =  array([1.907,  0.073,  1.830,  0.899,  0.425,  0.113])
#
## The interface to hessian computation is designed for general functions, so we
## need to define the specific functions of interest here. These functions
## calculate -ll given the logs of the parameters. (Because we work in log
## parameters, the uncertainties we estimate will be *relative* parameter
## uncertainties.)
#from dadi.Inference import ll_multinom
#func = lambda lp: -ll_multinom(func_ex(numpy.exp(lp), ns, pts_l), data)
#foldfunc = lambda lp: -ll_multinom(func_ex(numpy.exp(lp), ns, pts_l).fold(), 
#                                   data.fold()) 
#
## Calculate the two hessians
#h = dadi.Hessian.hessian(func, numpy.log(params), 0.05)
#hfold = dadi.Hessian.hessian(foldfunc, numpy.log(pfold), 0.05)
#
## Now we calculate the *relative* parameter uncertainties.
#uncerts = numpy.sqrt(numpy.diag(numpy.linalg.inv(h)))
#uncerts_folded = numpy.sqrt(numpy.diag(numpy.linalg.inv(hf)))
#
## The increase in uncertainty is not too bad. Tp increasing by 50% is the only
## substantial one.
#print uncerts_folded/uncerts - 1
