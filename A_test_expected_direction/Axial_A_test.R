#FOR A 'CLEAN' RUN, RESTART Rstudio
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 01 24
#     MODIFIED:	James Foster              DATE: 2024 01 24
#
#  DESCRIPTION: Generates axial data and finds and plots the assumed
#               means and mean vectors of a symmetrical axial distribution,
#               for both the maximum likelihood mean and a user-defined expected mean.
#               
#       INPUTS: 
#               
#      OUTPUTS: 
#
#	   CHANGES: - 
#	   
#   REFERENCES: Batschelet E (1981).
#               Multimodal Samples, Chap 1.6, p. 21
#               Chapter 1: Measures of Location
#               In: Circular Statistics in Biology
#               Academic Press (London)
#               
#               Fitak RR & Johnsen S. (2017)
#               Bringing the analysis of animal orientation data full circle: 
#               model-based approaches with maximum likelihood. 
#               J Exp Biol. 220(Pt 21):3878-3882.
#               doi: 10.1242/jeb.167056.
# 
#               Jammalamadaka, S. Rao and SenGupta, A. (2001).
#               Topics in Circular Statistics,
#               Section 2.2.4, Circular Distributions.
#               World Scientific Press, Singapore.
#
#    EXAMPLES:  Run section by section
#
# 

# seed random number generator for reproducibility ------------------------
set.seed(20171101)#publication date of Fitak & Johnsen 2017


# Input Variables ----------------------------------------------------------

##  User input -----------------------------------------------------------
expected_mean = 45#° a model-selection method will be used to check test this axis
csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"
angle_name = "angles" #The title of the column with angles; NO SPACES PLEASE
angle_unit = "degrees" # "degrees" or "radians"
point_col = 'darkblue' # colour the datapoints are plotted in
chosen_lines_col = adjustcolor('red', alpha.f = 0.5) # colour the mean chosen by model selection is plotted in
other_lines_col = adjustcolor('gray20', alpha.f = 0.5) # colour the other mean is plotted in
errbar_dist = 0.15 # distance between errorbars and the outside of the circle (in "lines of text")

#Check the operating system and assign a logical flag (T or F)
sys_win = Sys.info()[['sysname']] == 'Windows'
#On Windows computers, use user profile instead of home directory
if(sys_win){
  #get rid of all the backslashes
  ltp = gsub('\\\\', '/', Sys.getenv('USERPROFILE'))#Windows makes this difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp = Sys.getenv('HOME')#Easier on Mac
}

# Useful functions --------------------------------------------------------
## Load package ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    require(circular)#package for handling cirular data
    require(CircStats)#package for circular hypothesis tests
    # require(CircMLE)#package for fitting circular distributions
  }
)

# Generate samples from circular distributions ----------------------------

# To demonstrate the method, we can generate a random bimodal distribution.

#set n
n_angles = 20
#choose a random mean direction
suppressWarnings( {mu1 = rcircularuniform(1) })
#choose a mean vector length
mvl1 = 0.6 #quite oriented
#for reference the standard deviation should be nearly 60°
sqrt(-2*log(mvl1)) * 180/pi
# N.B. for this sample size p = 0.05 for a mean vector length of:
sqrt(-log(0.05)/n_angles)
# [1] 0.3870228
#derive the concentration parameter "kappa"
kappa1 = A1inv(mvl1)

## Generate a unimodal sample --------------------------------------------
angles_sim1 = suppressWarnings(
  rvonmises(n = n_angles,
            mu = mu1,
            kappa = kappa1,
            control.circular = list(units = angle_unit)
  )
)
# print(angles_sim1)
##   Circular Data: 
##    Type = angles 
##    Units = degrees 
##    Template = none 
##    Modulo = asis 
##    Zero = 0 
##    Rotation = counter 
##    [1] 288.27219 358.83880  50.34669  43.94500 319.60549 115.80873 354.46077
##    [8] 165.68954  86.77692  51.42809 122.55517 301.85103  83.04510 352.74657
##    [15]  58.70790  70.50992 106.40028  25.57120 356.10757  53.97658

## Generate complementary unimodal sample --------------------------------

angles_sim2 = suppressWarnings(
  rvonmises(n = n_angles,
            mu = mu1+pi, #the mean is 180° from mu1 in angles_sim1
            kappa = kappa1,
            control.circular = list(units = angle_unit)
  )
)

#N.B. Both samples individually are non-uniform according to a Rayleigh test

# rayleigh.test(angles_sim1)
##   Test Statistic:  0.5422 
##   P-value:  0.002 

# rayleigh.test(angles_sim2)
##   Test Statistic:  0.7164 
##   P-value:  0 

# Combine the two samples -------------------------------------------------


#the two samples now make up a single bimodal sample
angles_sim = c(angles_sim1, angles_sim2)
#we know that each of these distibutions is centred on a mean
#and that those two means are separated by 180°
#so we have an axial distribution by combining them

#using Batchelet's method of doubling angles does not identify a non-uniform distribution
# angles_doubled = angles_sim*2 #doubled angles map angles on either side of the circle onto one-another
# rayleigh.test(angles_doubled) #this method often works, but not always
#    Test Statistic:  0.2285 
#    P-value:  0.1238 



# Fit an axial von Mises distribution -------------------------------------

#It is more practical to fit the axial distribution using maximum likelihood methods.
#In this case I'll assume: 
#a "symmetric" distribution in which equal parts of the dataset
#are centred around each mean ()
#which is "axial", so the two means are separated by 180°.


## Set up modelling functions --------------------------------------------

#To find the _maximum_ likelihood mean and concentration parameters
#we can _minimise_ _negative_ log-likelihood 
#for a symmetrical axial von Mises distribution
#see also CircMLE::M3A by Fitak & Johnsen


### log-likelihood for the input parameters given the observed ang --------

SymAxLL = function(angles, #observed angles
                   params, #mean (mu) and concentration (kappa)
                   cc = list(type = 'angles',
                             unit = 'degrees',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = 'clock')
)
{
  suppressWarnings(
    {
  sum(log(#N.B. dmixedvonmises does not produce log scaled density because the two distributions need to be added
    dmixedvonmises(x = as.circular(x = angles, #angles are in optional format cc
                                   control.circular = cc),
                   mu1 = circular::rad(params[1]),#mean was in degrees, not in radians
                   mu2 = circular::rad(params[1])+pi, #mean 2 is exactly 180° away
                   kappa1 = params[2], #concentration between 0 and infinity
                   kappa2 = params[2], #concentration between 0 and infinity
                   prop = 0.5) #expect exactly half of the data centred on each mean (symmetric)
  ))
    }
  )
}


### negative log likelihood, and set reasonable parameter ranges ---------

LL_SymAxVM = function(params,  #proposed mean (mu) and concentration (kappa)
                      angles,  #observed angles
                      cc = list(type = 'angles',
                                unit = 'degrees',
                                modulo = '2pi',
                                zero = pi/2,
                                rotation = 'clock'),
                      mu_expected = NULL, #expected mean (if specified, always use this to calculate likelihood)
                      kappa_expected = NULL #expected concentration (if specified, always use this to calculate likelihood)
                      )
{ #set boundaries that indicate bad parameter values
  if((params[1] < 0 | params[1] > 360-1e-16 | #only look for means between 0° and 360°
      params[2] <=  0 | params[2] > 2e2)) #concentrations equivalent to sd of between 4° and infinite size 
  {ll = 1e9}else{#return large number if parameters out of bounds
  if(!is.null(mu_expected)){params[1] = mu_expected} #if expected mean is specified, overwrite parameter
  if(!is.null(kappa_expected)){params[2] = kappa_expected}#if expected kappa is specified, overwrite parameter
  ll = -SymAxLL(angles = angles, #return the negative log likelihood (optimiser will minimise this function)
                params = params) #for these observations and parameters
  }
  return(ll)
}

# #validation (not run)
# # two means 180° apart generate the same negative log likelihood
# LL_SymAxVM(angles = angles_sim,
#            params = c(circular::deg(mu1),
#                       kappa1),
#            mu_expected = expected_mean)
# LL_SymAxVM(angles = angles_sim,
#            params = c(circular::deg(mu1)+180,
#                       kappa1)
# )
# # a sequence of angles around a circle
# seq_deg = circular::deg(mu1) + seq(from = 0, to  = 360, by = 5)
# #paired with the true concentration kappa in a list
# seq_param = lapply(X = seq_deg, FUN = c, kappa1)
# # calculate the negative log likelihood for all parameter pairs
# ll_deg = sapply(X = seq_param, 
#               FUN = LL_SymAxVM,
#               angles = angles_sim)
# #a sinusoidal sequence that dips twice at two means 180° apart
# plot(x = seq_deg, 
#      ll_deg, 
#      type = 'l')
# #the true mean is at this red line
# abline(v = circular::deg(mu1 + c(0,pi)), lty = 2, col = 2)
# #the mean that would be recovered is at the blue line, nearby
# abline(v = seq_deg[which.min(ll_deg)] + c(0,180), 
#        lty = 1,
#        col = 4)


### Set up an optimisation function that will search for the maxi --------
OptFun = function(angles,
                  prm = c(mu = runif(n = 1, min = 0, max = 359), 
                          kappa = exp(runif(n = 1, min = -1, max = 1))
                          ),
                  iter = 5e3,
                  mu_expected = NULL,
                  kappa_expected = NULL)
{
  stats::optim(prm, 
               fn = LL_SymAxVM, 
               method = "BFGS",
               control = list(maxit = iter), 
               hessian = T,
               angles = angles,
               mu_expected = mu_expected,
               kappa_expected = kappa_expected)
}

#make a looping version
LoopOpt = function(i, ...)
            {OptFun(...)}

## Find the maximum likelihood parameter estimates -----------------------

#to sample efficiently around the circle,
#we will run the optimisation multiple times
#in different "chains", saving the parameters
#that have the highest likelihood across all chains
nchains = 8
#takes less than 1 second


### Find the maximum likelihood axial mean -------------------------------

#loop through all chains and save the output
# system.time(
#   {
mle_optim = sapply(X = 1:nchains,
             FUN = LoopOpt,
             angles = angles_sim)

#   }
# )

#each entry in mle_optim now contains the results of one chain
# print(mle_optim[1])
#   [[1]]
#   mu        kappa 
#   1.652114e+02 9.284753e-04 
# print(mle_optim['value',1])
#   $value
#   [1] 73.51508

#find the entry with the highest likelihood (minimum negative log-likelihood)
params_mle = mle_optim[, #extract the column
                       which.min(mle_optim['value',]) #with the lowest value
                       ]$par #and find the parameters

#recalculate the log-likelihood
ll_sym_ax = SymAxLL(angles = angles_sim, #for the observed angles
                    params = params_mle) #and the fitted parameters


### Find the maximum likelihood kappa at the expected mean --------

#loop through all chains and save the output
expected_optim = sapply(X = 1:nchains,
                        FUN = LoopOpt,
                        angles = angles_sim, #the observed angles
                        mu_expected = expected_mean) #the mean is expected and not estimated

#find the entry with the highest likelihood (minimum negative log-likelihood)
params_expect = c(mu = expected_mean, #mean is not estimated
                  kappa = expected_optim[,#extract the column
                                         which.min(expected_optim['value',]) #lowest negative log-likelihood
                                         ]$par['kappa'] #just the kappa parameter
                  )
#recalculate the log-likelihood
ll_expect = SymAxLL(angles = angles_sim,#for the observed angles
                    params = params_expect)#and the fitted parameters


### Find the likelihood of a uniform distribution ------------------------

# Just to confirm that the data are not uniformly distributed
ll_uniform = sum(log(dcircularuniform(x = angles_sim)) #return log probability
)


# Model comparison --------------------------------------------------------


## convert likelihood to deviance ----------------------------------------

#more likely = lower deviance
dev_symax_mle = -ll_sym_ax*2
dev_symax_expect = -ll_expect*2
dev_uniform = -ll_uniform*2 #deviance is -loglikelihood x 2


## Test for uniformity --------------------------------------------------

#if the uniform distribution has lower deviance
if(dev_uniform < dev_symax_mle) #test for significantly lower deviance
{
#Does the estimated mean significantly reduce deviance?
pLR_unif = pchisq(q = dev_symax_mle - dev_uniform, # likelihood ratio chi-squared
             df = 2, # two parameter difference (uniform distribution has no parameters)
             lower.tail = TRUE)#!(ll_sample_mean > ll_uniform)) # p(larger deviance): null hypothesis, expected mean is true mean
print(data.frame(chi_squared = round(dev_symax_mle - dev_uniform, 3),
                 d.f. = 2,
                 p = round(pLR_unif, 4),
                 h2 = if(pLR_unif <0.05)
                 {'the sample is uniformly distributed rather than axially distributed'}else
                 {'the sample may be uniformly distributed or saxially distributed'})
)
}

## Test for difference from expected mean --------------------------------


#compare expected mean with the 2nd most likely distribution
if(ll_sym_ax > ll_uniform)
{
  #calculate the likelihood ratio chi-squared statistic
  lr_stat = dev_symax_expect - dev_symax_mle # should be positive, sample mean always more likely
}else
{ # N.B. This is unnecessary, in this case all distributions have the same likelihood
  #if uniform is more likely, compare expected against uniform
  lr_stat = dev_uniform -  dev_symax_expect# should be positive
}
#Does the estimated mean significantly reduce deviance?
pLR = pchisq(q = lr_stat, # likelihood ratio chi-squared
             df = 1, # one parameter difference (the mean is only estimated if not expected)
             lower.tail = !(ll_sym_ax > ll_uniform)) # p(larger deviance): null hypothesis, expected mean is true mean
print(data.frame(chi_squared = round(lr_stat, 3),
                 d.f. = 1,
                 p = round(pLR, 4),
                 h2 = if(pLR <0.05)
                 {'the sample axial mean differs significantly from the expected mean'}else
                 {'the sample axial mean _does not_ differ significantly from the expected mean'})
)


# Plot the bimodal sample -------------------------------------------------

par(mar = c(0,0,0,0))# make space for plotting outside of the circle
#plot on a circular plot
plot.circular(x = circular(x = angles_sim, 
                           type = 'angles',
                           unit = 'degrees',
                           modulo = '2pi',
                           zero = pi/2,
                           rotation = 'clock'
                ),
                stack = TRUE,
                bins = 360/5,
                sep = 0.05,
                col = point_col,
                xlim = c(-1,1)*1.2,# make space for plotting outside of the circle
                ylim = c(-1,1)*1.2# make space for plotting outside of the circle
              )

#add the test results at the top
mtext(text = paste('chi2 =', round(lr_stat, 3), 
            ', d.f. = 1,',
             'p(expected axis) =', round(pLR, 4)),
      line = -1.5)

#plot the expected mean axis
arrows.circular(x = circular(x = c(expected_mean, expected_mean),
                             type = 'angles',
                             unit = 'degrees',
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = 'clock'
                            ),
                y = c(-1,1),
                length = 0,
                lty = 3, 
                col = 8
)

#plot the maximum likelihood mean direction and mean vector
arrows.circular(x = circular(x = params_mle[1],
                              type = 'angles',
                              unit = 'degrees',
                              modulo = '2pi',
                              zero = pi/2,
                              rotation = 'clock'
                              ) + c(0,pi),
                y = rep(x = A1(params_mle[2]),#convert kappa to mean vector
                        times = 2),
                lwd = 2, 
                col = if(pLR <0.05){chosen_lines_col}else
                                    {other_lines_col},
                length = 0.1
                )

#plot the expected direction and mean vector along that axis
arrows.circular(x = circular(x = c(expected_mean, expected_mean),
                             type = 'angles',
                             unit = 'degrees',
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = 'clock'
                            ),
                y = c(-1,1)* A1(params_expect[2]),
                length = 0.1,
                lwd = 2, 
                col = if(pLR <0.05){other_lines_col}else
                                   {chosen_lines_col}
)


