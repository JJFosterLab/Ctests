#FOR A 'CLEAN' RUN, RESTART Rstudio
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 01 24
#     MODIFIED:	James Foster              DATE: 2025 07 16
#
#  DESCRIPTION: Generates axial data and finds and plots the assumed
#               means and mean vectors of a symmetrical axial distribution,
#               for both the maximum likelihood mean and a user-defined expected mean.
#               
#       INPUTS: 
#               
#      OUTPUTS: 
#
#	   CHANGES: - Added plotting functions from GUV dances
#	            - Better handling of circular units in dmixedvonmises
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
#TODO   ---------------------------------------------
#TODO   
#- Test with simulated data +
#- Fix plotting +
#- Perform test    +
#- Read in data   
#- Save results 
#- Comment in detail

# seed random number generator for reproducibility ------------------------
set.seed(20171101)#publication date of Fitak & Johnsen 2017


# Input Variables ----------------------------------------------------------

##  User input -----------------------------------------------------------
expected_mean = 45#° a model-selection method will be used to check test this axis
csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"
angle_name = "angles" #The title of the column with angles; NO SPACES PLEASE
angle_unit = "degrees" # "degrees" or "radians"
angle_rot = "clock" # "clock" or "counter"
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


# Plot MLE parameters -----------------------------------------------------
#convert angles to signed angles in (-180, 180)
Mod360.180 = function(x)
{#use atan2 to convert any angle to the range (-180,180)
  deg(
    atan2(y = sin(rad(x)),
          x = cos(rad(x))
    )
  )
}


#generic mean angle simulator
MeanRvm = function(n, #representative sample size
                   mu = circular(0), #mean (defaults to 0rad)
                   kappa, #kappa required
                   au = 'degrees', #units
                   ar = 'clock') #rotation direction
{
  mean.circular(rvonmises(n = n, 
                          mu = circular(mu, units = au, rotation = ar), 
                          kappa = kappa,
                          control.circular = list(units = au, rotation = ar)))
}


#Simulate confidence intervals for a unimodal or bimodal distribution
#fitted to a vector of "angles"
CI_vM = function(angles, #vector of angles fitted (used for sample size)
                 m1, #primary mean
                 k1, #primary concentration
                 m2 = NA, #secondary mean (ignored if NULL or NA)
                 k2 = NA, #secondary kappa
                 w1 = 1, #weighting of primary mean
                 n = 1e4, #number of simulations
                 au = 'degrees', 
                 ar = 'clock',
                 calc_q = TRUE,
                 alternative = 'one.sided', #two.sided less conservative
                 interval = 0.95, #confidence interval to calculate
                 speedup_parallel = TRUE
)
{
  if(speedup_parallel) #3x faster
  {
    cl = parallel::makePSOCKcluster(parallel::detectCores()-1)
    parallel::clusterExport(cl = cl, 
                            varlist = c('mean.circular',
                                        'circular',
                                        'rvonmises'),
                            envir = .GlobalEnv
    )
    parallel::clusterExport(cl = cl, 
                            varlist = c('MeanRvm',
                                        'angles',
                                        'm1',
                                        'k1',
                                        'm2',
                                        'k2',
                                        'w1',
                                        'n',
                                        'au',
                                        'ar'),
                            envir = environment()
    )
    #simulate primary mean
    m1_est = 
      parallel::parSapply(cl = cl,
                          X = 1:n,
                          FUN = function(i)
                          {
                            eval.parent(
                              {
                                MeanRvm(n = round(length(angles)*w1), #estimate number of observations at primary mean
                                        mu = m1, 
                                        kappa = k1,
                                        au = au,
                                        ar = ar)
                              }
                            )
                          },
                          simplify = 'array' #return an array of simulated angles
      )
    if(!is.na(m2)) #if there is a valid secondary mean
    {
      m2_est = 
        parallel::parSapply(cl = cl,
                            X = 1:n,
                            FUN = function(i)
                            {
                              eval.parent(
                                {
                                  MeanRvm(n = round(length(angles)*(1-w1)), #estimate number of observations at secondary mean
                                          mu = m2, 
                                          kappa = k2,
                                          au = au,
                                          ar = ar)
                                }
                              )
                            },
                            simplify = 'array' #return an array of simulated angles
        )
    }
    parallel::stopCluster(cl)
  }else
  { #if not using parallel, use the slower version via replicate()
    m1_est = replicate(n = n, 
                       MeanRvm(n = round(length(angles)*w1), 
                               mu = m1, 
                               kappa = k1,
                               au = au,
                               ar = ar)
    )
    if(!is.na(m2))
    {
      m2_est = replicate(n = n, 
                         MeanRvm(n = round(length(angles)*(1-w1)), 
                                 mu = m2, 
                                 kappa = k2,
                                 au = au,
                                 ar = ar)
      )
    }
  }
  return(
    if(calc_q) #calculate quantiles only if requested
    {
      #either two-sided, symmetrical around mean change
      #or one-sided, from zero change towards mean change
      probs1 = switch(alternative,
                      two.sided = sort(c(c(0,1)+c(1,-1)*(1-interval)/2, 0.5)),
                      one.sided = sort(c(c(0,1)+
                                           (if(Mod360.180(m1)>0) #N.B. quantile.circular counts anticlockwise
                                           {c(1,0)}else
                                           {c(0,-1)}
                                           )*(1-interval), 0.5)),
                      sort(c(c(0,1)+ #default to one-sided
                               (if(Mod360.180(m1)>0)
                               {c(1,0)}else
                               {c(0,-1)}
                               )*(1-interval), 0.5))
      )
      if(is.na(m2))
      {
        Mod360.180(
          quantile.circular(x = circular(x = m1_est,
                                         units = au,
                                         rotation = ar),
                            probs = probs1)
        )
      }else
      {
        probs2 = switch(alternative,
                        two.sided = sort(c(c(0,1)+c(1,-1)*(1-interval)/2, 0.5)),
                        one.sided = sort(c(c(0,1)+
                                             (if(Mod360.180(m2)>0)
                                             {c(1,0)}else
                                             {c(0,-1)}
                                             )*(1-interval), 0.5)),
                        sort(c(c(0,1)+ #default to one-sided
                                 (if(Mod360.180(m2)<0)
                                 {c(1,0)}else
                                 {c(0,-1)}
                                 )*(1-interval), 0.5))
        )
        list(m1 = Mod360.180(
          quantile.circular(x = circular(x = m1_est,
                                         units = au,
                                         rotation = ar),
                            probs = probs1)
        ),
        m2 = Mod360.180(
          quantile.circular(x = circular(x = m2_est,
                                         units = au,
                                         rotation = ar),
                            probs = probs2)
        )
        )
      }
    }else
    { #if quantiles not requested, return the simulations (mainly for troubleshooting)
      if(is.na(m2))
      {
        m1_est = 
          sapply(X = m1_est,#warning, check units and direction
                 FUN = Mod360.180)
        # m1_est = 90 
      }else
      {
        list(
          m1_est = 
            sapply(X = m1_est,#warning, check units and direction
                   FUN = Mod360.180),     
          m2_est = 
            sapply(X = m2_est,#warning, check units and direction
                   FUN = Mod360.180),
        )
      }
    }
  )
}

PlotCI_vM = function(ci_vec,
                     col = 'salmon',
                     lwd = 2,
                     radius = 0.95,
                     ...)#passed to lines()
{
  angle_seq1 = 
    c(
      seq(from = ci_vec[1], #lower
          to = ci_vec[1] +
            Mod360.180(ci_vec[2]-ci_vec[1]), #median
          length.out =1e2/2),
      seq(from = ci_vec[2], #median
          to = ci_vec[2] +
            Mod360.180(ci_vec[3]-ci_vec[2]) , #upper
          length.out =1e2/2)
    )
  lines(x = radius*sin( rad(angle_seq1) ),
        y = radius*cos( rad(angle_seq1) ),
        col = col,
        lwd = 2,
        lend = 'butt',
        ...
  )
  if(!is.na(ci_vec[4]))
  {
    angle_seq2 = 
      c(
        seq(from = ci_vec[1+3],
            to = ci_vec[1+3] +
              Mod360.180(ci_vec[2+3]-ci_vec[1+3]),
            length.out =1e2/2),
        seq(from = ci_vec[2+3],
            to = ci_vec[2+3] +
              Mod360.180(ci_vec[3+3]-ci_vec[2+3]) ,
            length.out =1e2/2)
      )
    lines(x = radius*sin( rad(angle_seq2) ),
          y = radius*cos( rad(angle_seq2) ),
          col = col,
          lwd = 2,
          lend = 'butt',
          ...
    )
  }
}




PlotMV_circMLE = function(mod_par,
                          au = 'degrees',
                          ar = 'clock',
                          col1 = 'salmon',
                          col2 = col1,
                          ...) #passed to arrows.circular()
{
  m1 = with(mod_par, 
            {
              circular(mu1,
                       unit = au,
                       rotation = ar,
                       modulo = '2pi',
                       zero = pi/2
              )
            }
  )
  with(mod_par,
       {
         arrows.circular(x = m1,
                         y = A1(kappa1),
                         col = col1,
                         lwd = 5*weight1,
                         length = 0.1,
                         ...
         )
       }
  )
  if(!is.na(mod_par$mu2))
  {
    m2 = with(mod_par, 
              {
                circular(mu2,
                         unit = au,
                         rotation = ar,
                         modulo = '2pi',
                         zero = pi/2
                )
              }
    )
    with(mod_par,
         {
           arrows.circular(x = m2,
                           y = A1(kappa2),
                           col = col2,
                           lwd = 5*(1-weight1),
                           length = 0.1,
                           ...
           )
         }
    )
  }
  
}



    # # Generate samples from circular distributions  (not used)----------------------------
    # 
    # # To demonstrate the method, we can generate a random bimodal distribution.
    # 
    # #set n
    # n_angles = 20
    # #choose a random mean direction
    # suppressWarnings( {mu1 = rcircularuniform(1) })
    # # mu1 = circular(pi/2 +0.2, template ='none')
    # #choose a mean vector length
    # mvl1 = 0.6 #quite oriented
    # #for reference the standard deviation should be nearly 60°
    # sqrt(-2*log(mvl1)) * 180/pi
    # # N.B. for this sample size p = 0.05 for a mean vector length of:
    # sqrt(-log(0.05)/n_angles)
    # # [1] 0.3870228
    # #derive the concentration parameter "kappa"
    # kappa1 = A1inv(mvl1)

    # ## Generate a unimodal sample --------------------------------------------
    # adata1 = suppressWarnings(
    #   rvonmises(n = n_angles,
    #             mu = mu1,
    #             kappa = kappa1,
    #             control.circular = list(units = angle_unit)
    #   )
    # )
    # # print(adata1)
    # ##   Circular Data: 
    # ##    Type = angles 
    # ##    Units = degrees 
    # ##    Template = none 
    # ##    Modulo = asis 
    # ##    Zero = 0 
    # ##    Rotation = counter 
    # ##    [1] 288.27219 358.83880  50.34669  43.94500 319.60549 115.80873 354.46077
    # ##    [8] 165.68954  86.77692  51.42809 122.55517 301.85103  83.04510 352.74657
    # ##    [15]  58.70790  70.50992 106.40028  25.57120 356.10757  53.97658
    # 
    # ## Generate complementary unimodal sample --------------------------------
    # 
    # adata2 = suppressWarnings(
    #   rvonmises(n = n_angles,
    #             mu = mu1+pi, #the mean is 180° from mu1 in adata1
    #             kappa = kappa1,
    #             control.circular = list(units = angle_unit)
    #   )
    # )
    # 
    # #N.B. Both samples individually are non-uniform according to a Rayleigh test
    # 
    # # rayleigh.test(adata1)
    # ##   Test Statistic:  0.5422 
    # ##   P-value:  0.002 
    # 
    # # rayleigh.test(adata2)
    # ##   Test Statistic:  0.7164 
    # ##   P-value:  0 
    # 
    # # Combine the two samples -------------------------------------------------
    # 
    # 
    # #the two samples now make up a single bimodal sample
    # adata = c(adata1, adata2)
    # #we know that each of these distibutions is centred on a mean
    # #and that those two means are separated by 180°
    # #so we have an axial distribution by combining them
    # 
    # #using Batchelet's method of doubling angles does not identify a non-uniform distribution
    # # angles_doubled = adata*2 #doubled angles map angles on either side of the circle onto one-another
    # # rayleigh.test(angles_doubled) #this method often works, but not always
    # #    Test Statistic:  0.2285 
    # #    P-value:  0.1238 

# Select files ---------------------------------------------------------

# set path to files
if(sys_win){#choose.files is only available on Windows
  message('\n\nPlease select the ".csv" file\n\n')
  Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
  path_file  <- choose.files(
    default = file.path(ltp,'Documents', "*.csv"),#For some reason this is not possible in the "root" user
    caption = 'Please select the ".csv" file'
  )
}else{
  message('\n\nPlease select the ".csv" file\n\n')
  Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
  path_file <- file.choose(new=F)
}
#show the user the path they have selected
if(is.null(path_file))
{stop('No file selected.')}else
{print(path_file)}

#Check for Byte Order Marks, which can make a mess
if(grepl(x = readLines(path_file,
                       n = 1,
                       warn = F),#check the first line of the file, where "angle" should be written
         pattern = "ï|ÿ|þ")#common BOM renderings, are there any others?
)
{utf8BOM = T}else{utf8BOM = F}

# Read in file ------------------------------------------------------------
adata = read.table(file = path_file,#read from user-selected file
                   header = T,#read the file header to use for variable names
                   sep = csv_sep,#values are separated by the user-specified character
                   fileEncoding = ifelse(test = utf8BOM, #If the file contains Byte Order Markers
                                         yes = "UTF-8-BOM",#read in using the appropriate format
                                         no = "")#, #if not, R can guess
                   #other parameters can be added here for troubleshooting
)

View(adata)#show the user the data that was

angles_data = circular(x = adata,
                       units = angle_unit,
                       rotation = angle_rot)
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
                   cc = list(unit = 'degrees',
                             rotation = 'clock')
)
{
  rangles = circular(rad(angles),
                     template = 'none')
  sum(log(#N.B. dmixedvonmises does not produce log scaled density because the two distributions need to be added
    dmixedvonmises(x = rangles,
                   # x = as.circular(x = angles, #angles are in optional format cc
                   #              control.circular = cc),
                   mu1 = circular(x = rad(params[1]), #for some reason requires radians
                                  template = 'none'),#mean was in degrees, not in radians
                   mu2 = circular(x = rad(params[1]+180), #for some reason requires radians
                                  template = 'none'), #mean 2 is exactly 180° away
                   kappa1 = params[2], #concentration between 0 and infinity
                   kappa2 = params[2], #concentration between 0 and infinity
                   prop = 0.5) #expect exactly half of the data centred on each mean (symmetric)
  ))
}


### negative log likelihood, and set reasonable parameter ranges ---------

LL_SymAxVM = function(params,  #proposed mean (mu) and concentration (kappa)
                      angles,  #observed angles
                      cc = list(type = 'angles',
                                unit = 'degrees',
                                modulo = '2pi',
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
                params = params,
                cc = cc) #for these observations and parameters
  }
  return(ll)
}

# #validation (not run)
# # two means 180° apart generate the same negative log likelihood
# LL_SymAxVM(angles = adata,
#            params = c(circular::deg(mu1),
#                       kappa1),
#            mu_expected = expected_mean)
# LL_SymAxVM(angles = adata,
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
#               angles = adata)
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
                  kappa_expected = NULL,
                  cc = list( unit = 'degrees',
                            rotation = 'clock'))
{
  stats::optim(prm, 
               fn = LL_SymAxVM, 
               method = "BFGS",
               control = list(maxit = iter), 
               hessian = T,
               angles = angles,
               mu_expected = mu_expected,
               kappa_expected = kappa_expected,
               cc = cc)
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
             angles = angles_data,
             cc = list(units = angle_unit,
                       rotation = angle_rot,
                       type = angle_name,
                       zero = 0, 
                       template = 'none', 
                       modulo = 'asis'))

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
ll_sym_ax = SymAxLL(angles = angles_data, #for the observed angles
                    params = params_mle) #and the fitted parameters


### Find the maximum likelihood kappa at the expected mean --------

#loop through all chains and save the output
expected_optim = sapply(X = 1:nchains,
                        FUN = LoopOpt,
                        angles = angles_data, #the observed angles
                        mu_expected = expected_mean) #the mean is expected and not estimated

#find the entry with the highest likelihood (minimum negative log-likelihood)
params_expect = c(mu = expected_mean, #mean is not estimated
                  kappa = expected_optim[,#extract the column
                                         which.min(expected_optim['value',]) #lowest negative log-likelihood
                                         ]$par['kappa'] #just the kappa parameter
                  )
#recalculate the log-likelihood
ll_expect = SymAxLL(angles = angles_data,#for the observed angles
                    params = params_expect)#and the fitted parameters


### Find the likelihood of a uniform distribution ------------------------

# Just to confirm that the data are not uniformly distributed
ll_uniform = sum(log(dcircularuniform(x = angles_data)) #return log probability
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
                 {'the sample may be uniformly distributed or axially distributed'})
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

circMLE_params = data.frame(mu1 = params_mle[1],
                            kappa1 = params_mle[2],
                            mu2 = 180+params_mle[1],
                            kappa2 = params_mle[2],
                            weight1 = 0.5
                            )

par(pty = 's')
par(mar = c(0,0,0,0))# make space for plotting outside of the circle
#plot on a circular plot
plot.circular(x = circular(angles_data,
                           units = angle_unit,
                           rotation = angle_rot,
                           zero = pi/2),
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
PlotMV_circMLE(mod_par = circMLE_params,
               au = angle_unit,
               col1 = if(pLR >0.05){other_lines_col}else
                                   {chosen_lines_col}
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


# Add confidence interval -------------------------------------------------


#plot the symmetric axial distribution CI around the mean
mle_ci = with(circMLE_params,
     {
      CI_vM(angles = angles_data,
            m1 = mu1,
            k1 = kappa1,
            m2 = mu2, 
            k2 = kappa2,
            w1 = 0.5,
            au = angle_unit,
            alternative = 'two.sided'
            )
     }
)

PlotCI_vM(ci_vec = unlist(mle_ci),
          col = if(pLR >0.05){other_lines_col}else
                            {chosen_lines_col},
          radius = 1+0.15)
