#FOR A 'CLEAN' RUN, RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 03 22
#     MODIFIED:	James Foster              DATE: 2024 04 02
#
#  DESCRIPTION: Loads a text file with paired data and performs an 
#               "inverse" V test for directedness towards an expected mean angle. 
#               Could be called a paired "A" test.
#               
#       INPUTS: A ".csv" table with a column of angles ("angle").
#               User should specify test details (line 50).
#               
#      OUTPUTS: Results table (.csv).
#
#	   CHANGES: - Combined 3 model selection procedure
#	            - Test on paired differences
#	            - Single script for paired and unpaired
#
#   REFERENCES: Batschelet E (1981).
#               The Rayleigh test, Chap 4.2, p. 54
#               Chapter 4: Tests for Randomness and Goodness-of-fit
#               In: Circular Statistics in Biology
#               Academic Press (London)
#
#               Woolf, B. (1957). 
#               THE LOG LIKELIHOOD RATIO TEST (THE G‐TEST).
#               Annals of Human Genetics 21, 397–409.
#
#    EXAMPLES:  Fill out user input (lines 50-55), then press ctrl+shift+s to run
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Read in data   +
#- Perform test   + 
#- Test with simulated data +
#- Save results +
#- Accurate simulation +
#- Appropriate comparison strategy +
#- Neaten up model comparison +
#- A test on differences +
#- V-test on grand mean + 
#- Simplify to test of differences +
#- Rotation direction option (needed for plots) +
#- Troubleshoot plotting direction  +
#- Tie-breaking approach
#- Comment in detail
#- Optimisation method for pairs (ML too biased?)

# Useful functions --------------------------------------------------------

# . Load packages ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    require(circular)#package for handling circular data
    require(CircStats)#package for circular hypothesis tests
  }
)


# . General functions -----------------------------------------------------

Mod360.180 = function(x)
{#use atan2 to convert any angle to the range (-180,180)
  deg(
    atan2(y = sin(rad(x)),
          x = cos(rad(x))
    )
  )
}



# Input Variables ----------------------------------------------------------

#  .  User input -----------------------------------------------------------
paired_data = TRUE # Are the data in the two columns paired (each from the same animal or group)?
csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"
angle_name = "angles" #The title of the column with angles; NO SPACES PLEASE
angle_unit = "degrees" # "degrees" or "radians"
angle_rot = "clock" # "clock" or "counter"

#Check the operating system and assign a logical flag (T or F)
sys_win <- Sys.info()[['sysname']] == 'Windows'
#User profile instead of home directory
if(sys_win){
  #get rid of all the backslashes
  ltp <- gsub('\\\\', '/', Sys.getenv('USERPROFILE'))#Why does windows have to make this so difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp <- Sys.getenv('HOME')#Life was easier on Mac
}

  # # Simulate data (not used) ------------------------------------------------
  # n_angles = 44
  # # minimum discriminable angle appears to be approx 35°
  # mu_offset = rad(35)
  # kappa_both = A1inv(0.7) #concentration around each trial mean
  # logkappa_var = 1.0 #scale of random variation in concentration (log units)
  # if(paired_data)
  # {
  # kappa_indiv = A1inv(0.98) #concentration across individuals (pairs)
  # #mean angle in trail 1 for each individual (pair)
  # mu1_sim = rvonmises(n = n_angles,
  #                       mu = rcircularuniform(1),#random angle
  #                       kappa = kappa_indiv#the wider the distribution of individual biases, the greater the influence of pairing
  #                       )
  # #simulate the full dataset
  # sim = data.frame(
  #                  angle_1 = round(c(suppressWarnings( #rvonmises converts to circular and warns
  #                    sapply(X = mu1_sim,
  #                           FUN = rvonmises,
  #                           n = 1,
  #                            kappa = exp(log(kappa_both) + rnorm(n = 1, sd = logkappa_var))
  #                    )
  #                  ))*180/pi),#convert to angles and round to nearest degree
  #                  angle_2 = round(c(suppressWarnings( #rvonmises converts to circular and warns
  #                              sapply(X =mu1_sim + mu_offset,# true difference,
  #                                     FUN = rvonmises,
  #                                     n = 1,
  #                                     kappa = exp(log(kappa_both) + rnorm(n = 1, sd = logkappa_var))
  #                              )
  #                  ))*180/pi) #convert to angles and round to nearest degree
  #                  )
  # }else
  # {
  # n_angles2 = ceiling(0.75*n_angles)
  # mu1_sim = rcircularuniform(n = 1,control.circular = list(units = angle_unit))
  # sim = data.frame(
  #                  angle_1 = round(c(suppressWarnings( #rvonmises converts to circular and warns
  #                    rvonmises(n = n_angles,
  #                              mu = circular(x = mu1_sim, units = angle_unit),
  #                            kappa = exp(log(kappa_both) + rnorm(n = 1, sd = logkappa_var))
  #                    )
  #                  ))),#convert to angles and round to nearest degree
  #                  angle_2 = round(c(suppressWarnings( #rvonmises converts to circular and warns
  #                    rvonmises(n = n_angles2,
  #                              mu = circular(x = mu1_sim+deg(mu_offset), units = angle_unit),
  #                              kappa = exp(log(kappa_both) + rnorm(n = 1, sd = logkappa_var))
  #                    )
  #                  ),
  #                  circular(x = rep(x = NA, times = n_angles - n_angles2),
  #                          units = angle_unit) #convert to angles and round to nearest degree
  #                  ) )
  #                 )
  # }
  # #save somewhere the user likely keeps data
  # write.table(x = sim,
  #             file = file.path(ltp,'Documents', "simulated_angles.csv"),
  #             sep = csv_sep,
  #             row.names = FALSE
  #             )

# . Select files ---------------------------------------------------------

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
dt_dim = dim(adata)
if(dt_dim[2] != 2)
{stop('did not find two samples, please check that data are in two columns')}

ref_angle = with(adata, 
                 {
                   mean.circular(x = circular(x = c(angle_1,angle_2),
                                                          type = 'angles',
                                                          unit = angle_unit,
                                                          template = 'geographics',
                                                          modulo = '2pi',
                                                          zero = pi/2,
                                                          rotation = angle_rot),
                                 na.rm = TRUE
                               )
                 } 
                 )

# Plot the data -----------------------------------------------------------

par(mar =rep(0,4))
plot.circular(x = circular(x = adata$angle_1, 
                           type = 'angles',
                           unit = angle_unit,
                           template = 'geographics',
                           modulo = '2pi',
                           zero = pi/2,
                           rotation = angle_rot
),
stack = TRUE,
bins = 360/5,
sep = 0.5/dt_dim[1],
col = 'cyan4'
)
par(new = T)
plot.circular(x = circular(x = adata$angle_2, 
                           type = 'angles',
                           unit = 'degrees',
                           template = 'geographics',
                           modulo = '2pi',
                           zero = pi/2,
                           rotation = angle_rot
),
stack = TRUE,
bins = 360/5,
sep = -0.5/dt_dim[1],
col = 'darkblue',
shrink = 1.05,
axes = F
)
arrows.circular(x = circular(x = ref_angle, 
                             type = 'angles',
                             unit = angle_unit,
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot
),
y = 1,
lwd = 5, 
col = rgb(0,0,0,0.1),
length = 0
)
arrows.circular(x = circular(x = ref_angle, 
                             type = 'angles',
                             unit = angle_unit,
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot
                            ),
                y = mean(x = 
                  cos(x = with(adata, c(angle_1, angle_2) - ref_angle)* pi/180), 
                  na.rm = T)
                )

lines.circular(x = circular(x = 
                    seq(from = pi,#or: bearing - pi/2
                       to = -pi, #or: bearing + pi/2
                       length.out = 1e3)
                    ),
               y = -1+rep( x = 1.644854/(sqrt(2 * length(adata$angle_2))), times = 1e3),
               lty = 2
)

# Perform V test ----------------------------------------------------------
#convert angles to circular class
cangs = with(adata, circular(x = c(angle_1 , 
                                    angle_2), 
                             unit = angle_unit,
                             rotation = angle_rot) )
emean = circular(ref_angle, unit = angle_unit, rotation = angle_rot)

rayv_test = rayleigh.test( x = cangs,
                           mu = emean
)
print(rayv_test)


# Fit maximum likelihood models -----------------------------------------------
# The v-test has identified a non-uniform distribution with
# a large component in the expected direction.


# . Fit models ------------------------------------------------------------



if(!paired_data)
{
# . . Unpaired data -------------------------------------------------------
  #for unpaired data
  #find the maximum likelihood von Mises distribution for the full dataset
  #with fewer degrees of freedom log likelihood is lower than dividing data into trials
  ml_grand_mean = with(adata,
             mle.vonmises(x = circular(c(angle_1,
                                         angle_2),
                                         units = angle_unit,
                                       rotation = angle_rot
                                       ),
                                  bias = TRUE)# correct bias
              )
  
  #find the maximum likelihood von Mises distribution for each column
  ml_sample_mean = with(adata,
             apply(X =  cbind(circular(angle_1,units = angle_unit, rotation = angle_rot),
                              circular(angle_2,units = angle_unit, rotation = angle_rot)),
                   MARGIN = 2,
                   FUN = function(angs)
                     {
                     mle.vonmises(x = circular(angs,units = angle_unit, rotation = angle_rot),
                                  bias = TRUE)
                     }
                   ) # correct bias
              )
  #find the maximum likelihood von Mises distribution for each column
  ml_sample_kappa = with(adata,
             apply(X =  cbind(circular(angle_1,units = angle_unit, rotation = angle_rot),
                              circular(angle_2,units = angle_unit, rotation = angle_rot)),
                   MARGIN = 2,
                   FUN = function(angs)
                     {
                     mle.vonmises(x = circular(angs,units = angle_unit, rotation = angle_rot),
                                  mu = circular(ml_grand_mean$mu,
                                                units = angle_unit,
                                                rotation = angle_rot),
                                  bias = TRUE)
                     }
                   ) # correct bias
              )
}else
{
# . . Paired data ---------------------------------------------------------
  #if data are paired, calculate the differences between the pairs
  #pairwise angular differences, positive indicates a turn to the right between trials
  pair_diffs = with(adata, 
                     {
                       apply(X = cbind(circular(angle_1,units = angle_unit, rotation = angle_rot), # subtract 1st angle
                                        circular(angle_2,units = angle_unit, rotation = angle_rot)), #from the 2nd angle
                                          MARGIN = 1, 
                                          FUN = function(angs)
                                          {
                                            diff(x = circular(x = as.numeric(angs),
                                                              units = angle_unit,
                                                              rotation = angle_rot) )
                                            
                                          }
                             )
                     }
                    )
  #fit the maximum likelihood distribution for differences between pairs
  # if angles shift in one direction between trials,
  # this distribution should have significantly higher likelihood
  ml_diff = mle.vonmises(x = circular(x = unlist(pair_diffs),
                                      units = angle_unit,
                                      rotation = angle_rot),
                        bias = TRUE)
  
  #fit the maximum likelihood distribution for differences centred on zero
  # if there is no consistent shift between trials,
  # this distribution should have similar likelihood,
  # with one less free parameter (expected mean of zero) 
  ml_same = mle.vonmises(x = circular(x = unlist(pair_diffs),
                                      units = angle_unit,
                                      rotation = angle_rot),
                         mu = circular(x = 0, units = angle_unit, rotation = angle_rot),
                         bias = TRUE)
}

# . Calculate log likelihoods ---------------------------------------------

if(!paired_data)
{
# . . Unpaired data ---------------------------------------------------------
  # What is the likelihood of orientation towards the grand mean?
  ll_grand_mean = with(ml_grand_mean, #using the maximum likelihood von Mises parameters
                        sum( # add together
                          dvonmises(x = circular(with(adata,c(angle_1,angle_2)),
                                                 units = angle_unit,
                                                 rotation = angle_rot), # probability density for each observed angle
                                    mu = mu, # ML estimated mean
                                    kappa = kappa, # ML estimated concentration
                                    log = TRUE), # on a log scale (i.e. add instead of multiplying)
                        na.rm = TRUE
                          )
                      ) 
  # What is the likelihood of orientation towards each trial's mean?
  ll_sample_mean = with(ml_sample_mean[[1]], #using the maximum likelihood von Mises parameters
                        sum( # add together
                          dvonmises(x = circular(adata$angle_1,
                                                 units = angle_unit,
                                                 rotation = angle_rot), # probability density for each observed angle
                                    mu = mu, # ML estimated mean
                                    kappa = kappa, # ML estimated concentration
                                    log = TRUE), # on a log scale (i.e. add instead of multiplying)
                          na.rm = TRUE
                          )
                      ) + 
                    with(ml_sample_mean[[2]], #using the maximum likelihood von Mises parameters
                        sum( # add together
                          dvonmises(x = circular(adata$angle_2,
                                                 units = angle_unit,
                                                 rotation = angle_rot), # probability density for each observed angle
                                    mu = mu, # ML estimated mean
                                    kappa = kappa, # ML estimated concentration
                                    log = TRUE), # on a log scale (i.e. add instead of multiplying)
                          na.rm = TRUE
                          )
                    )
  # What is the likelihood of orientation towards the grand mean, but with different concentrations?
  ll_sample_kappa = with(ml_sample_kappa[[1]], #using the maximum likelihood von Mises parameters
                        sum( # add together
                          dvonmises(x = circular(adata$angle_1,
                                                 units = angle_unit,
                                                 rotation = angle_rot), # probability density for each observed angle
                                    mu = mu, # ML estimated mean
                                    kappa = kappa, # ML estimated concentration
                                    log = TRUE), # on a log scale (i.e. add instead of multiplying)
                          na.rm = TRUE
                          )
                      ) + 
                    with(ml_sample_kappa[[2]], #using the maximum likelihood von Mises parameters
                        sum( # add together
                          dvonmises(x = circular(adata$angle_2,
                                                 units = angle_unit,
                                                 rotation = angle_rot), # probability density for each observed angle
                                    mu = mu, # ML estimated mean
                                    kappa = kappa, # ML estimated concentration
                                    log = TRUE), # on a log scale (i.e. add instead of multiplying)
                          na.rm = TRUE
                          )
                        ) 
  
  # Just to confirm the V-test, what is the likelihood of a uniform distribution?
  ll_uniform = with(adata,
                    sum(dvonmises(x = circular(adata$angle_2,
                                               units = angle_unit,
                                               rotation = angle_rot), # probability density for each observed angle
                                  mu = circular(0,template = 'none'), # ML estimated mean
                                  kappa = 0, # ML estimated concentration
                                  log = TRUE), #return log probability
                        na.rm = TRUE
                        )
  )
}else
{
# . . Paired data ---------------------------------------------------------
  #if data are paired, use the distribution of differences rather than the samples
  
  # calculate the likelihood of pair differences centred on the ML difference
  ll_diff = with(ml_diff, #using the maximum likelihood von Mises parameters
                       sum( # add together
                         dvonmises(x = circular(x = pair_diffs,
                                                units = angle_unit,
                                                rotation = angle_rot), # probability density for each observed angle
                                   mu = mu, # ML estimated mean
                                   kappa = kappa, # ML estimated concentration
                                   log = TRUE) # on a log scale (i.e. add instead of multiplying)
                       )
  )
  
  # calculate the likelihood of pair differences centred on zero
  ll_same = with(ml_same, #using the maximum likelihood von Mises parameters
                       sum( # add together
                         dvonmises(x = circular(x = pair_diffs,
                                                units = angle_unit,
                                                rotation = angle_rot), # probability density for each observed angle
                                   mu = mu, # ML estimated mean
                                   kappa = kappa, # ML estimated concentration
                                   log = TRUE) # on a log scale (i.e. add instead of multiplying)
                       )
  ) 
}

# Add ML models to the figure ---------------------------------------------


# Add the two distributions to the figure
aa = circular(x = seq(from = -180,#or: bearing - pi/2
                      to = 180, #or: bearing + pi/2
                      length.out = 1e3),
              unit = angle_unit,
              rotation = angle_rot)

if(!paired_data)
{
  # . . Unpaired data -------------------------------------------------------
  #probability density grand mean
  with(ml_grand_mean,
       lines.circular(x = 90-aa,
                      y = dvonmises(x = aa,
                                    mu = mu,
                                    kappa = kappa,
                                    log = FALSE) - 1,
                      lty = 3,
                      lwd = 2,
                      col = 'green4')
  )
  #probability density trial 1
  with(ml_sample_mean[[1]],
       lines.circular(x = 90-aa,
                      y = dvonmises(x = aa,
                                    mu = mu,
                                    kappa = kappa,
                                    log = FALSE) - 1,
                      lty = 3,
                      lwd = 2,
                      col = 'cyan3')
  )
  #probability density trial 2
  with(ml_sample_mean[[2]],
       lines.circular(x = 90-aa,
                      y = dvonmises(x = aa,
                                    mu = mu,
                                    kappa = kappa,
                                    log = FALSE) - 1,
                      lty = 3,
                      lwd = 2,
                      col = 'blue3')
  )

legend(x = 'bottomright',
       inset = c(0.01,0.02),
       cex = 0.4,
       bg = gray(level = 1,alpha = 0),
       bty = 'n',
       legend = c('vector to grand mean',
                  'p(>V) < 0.05',
                  'prob. dens.: grand mean',
                  'prob. dens.: sample mean 1',
                  'prob. den.: sample mean 2'),
       col = c('black',
               'black',
               'green4',
               'cyan4',
               'blue3'),
       lty = c(1,
               2,
               3,
               3,
               3)
)
}else
{


# . . Plot paired differences -------------------------------------------------


#N.B. unless the two means are identical, the likelihood of the expected mean
#is always lower than the maximum likelihood mean (it has _maximum_ likelihood).
#Our question is, it significantly more likely than the distribution around the 
#expected mean, for which we need to fit one fewer parameter (we already know the mean).
  #bootstrap CI
  # pair_diffs_CI = mle.vonmises.bootstrap.ci(x = circular(pair_diffs,units = angle_unit),
  #                                           bias = TRUE,
  #                                           reps = 1e4)
#simulate CI
  #generic mean angle calculator
MeanRvm = function(n, mu = circular(0), kappa, au = 'degrees', ar = 'clock')
{
  mean.circular(rvonmises(n = n, 
                          mu = circular(mu, units = au, rotation = ar), 
                          kappa = kappa,
                          control.circular = list(units = au, rotation = ar)))
}
#simulate von Mises distributions with the ML parameters and calculate the mean
sim_means = with(ml_diff, 
                 replicate(n = 1e4, 
                           MeanRvm(n = length(pair_diffs), 
                                   mu = mu, 
                                   kappa = kappa,
                                   au = angle_unit,
                                   ar = angle_rot)
                 )
)
#find the 95% confidence interval
sim_diffs_CI = Mod360.180(quantile.circular(x = circular(sim_means, units = angle_unit, rotation = angle_rot),
                                 probs = c(0,1)+0.5*c(1,-1)*0.05) )
# pair_diffs_CI = mle.vonmises.bootstrap.ci(x = circular(pair_diffs,units = angle_unit),
#                                           bias = TRUE,
#                                           reps = 1e4)

#open plot
par(mar =rep(0,4))
plot.circular(x = circular(x = pair_diffs, 
                           type = 'angles',
                           unit = angle_unit,
                           modulo = '2pi',
                           zero = pi/2,
                           rotation = angle_rot
),
stack = TRUE,
bins = 360/5,
sep = 0.5/dt_dim[1],
col = 'orange'
)
arrows.circular(x = circular(x = 0, 
                             type = 'angles',
                             unit = 'degrees',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot
),
y = 1,
lwd = 5, 
col = rgb(0,0,0,0.1),
length = 0
)
mtext(text = 'change in heading\ntrial 2 - trial 1',
      side = 3, 
      line = -17)
#plot mean vector for zero mean difference
with(ml_same,
     {
arrows.circular(x = circular(x = mu, 
                             type = 'angles',
                             unit = 'degrees',
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot
                            ),
                            y = A1(kappa),
                col = 'black',
                lwd = 2,
                length = 0.1
)
     }
)
#plot the maximum likelihood mean vector for differences
with(ml_diff,
     {
arrows.circular(x = circular(x = mu, 
                             type = 'angles',
                             unit = angle_unit,
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot
                            ),
                            y = A1(kappa),
                col = 'orange3',
                lwd = 2,
                length = 0.1
)
     }
)
#plot the simulated confidence interval for the ML mean
lines.circular(x = circular(x = 
                              seq(from = max(sim_diffs_CI),
                                   to = min(sim_diffs_CI),
                                   length.out = 1e3),
                            zero = pi/2,
                            units = angle_unit, 
                            rotation = angle_rot),  
               y = rep(x = 1-1.05,times = 1e3),
               col = 'orange3',
               lwd = 5,
               lend = 'butt')

legend(x = 'bottom',
       inset = c(0.2,0.2),
       cex = 0.7,
       bg = gray(level = 1,alpha = 0),
       bty = 'n',
       legend = c('vector to ML mean',
                  '95% CI for ML',
                  'vector for zero difference'),
       col = c('orange3',
               'orange3',
               'black'),
       lwd = c(2,
               5,
               2)
)


      
# with(pair_diffs_CI,
#      {
#       lines.circular(x = circular(x = 90-seq(from = max(mu.ci),
#                                    to = min(mu.ci),
#                                    length.out = 1e3),
#                         units = angle_unit),  
#                      y = rep(x = 1-0.9,times = 1e3),
#                      col = 'orange', 
#                      lwd = 1)
#        arrows.circular(x = circular(x = 90 -rep(x = ml_diff$mu,times = 2),
#                                     units = angle_unit,
#                                     rotation = angle_rot),
#                        y = A1(kappa.ci),
#                        code = 3,
#                        angle = 90,
#                        length = 0.1,
#                        col = 'orange'
#                        )
#      }
# )

}     
# Perform likelihood ratio test -------------------------------------------

#The likelihood ratio test compares models by estimating the likelihood gained
#for each additional parameter.
#https://en.wikipedia.org/wiki/Likelihood-ratio_test
# In our case, the ML von Mises requires two parameters (mean and concentration), 
# while the von Mises around the expected mean requires one (concentration).
# How much does the extra parameter increase the likelihood?

#convert likelihood to deviance (more likely = lower deviance)
if(!paired_data)
{
#on the original data
  dev_grand_mean = -ll_grand_mean*2 #deviance is -loglikelihood x 2
  dev_sample_mean = -ll_sample_mean*2 #deviance is -loglikelihood x 2
  dev_sample_kappa = -ll_sample_kappa*2 #deviance is -loglikelihood x 2
  dev_uniform = -ll_uniform*2 #deviance is -loglikelihood x 2
}else
{
#on the distribution of pairs
  dev_diff = -ll_diff*2 #deviance is -loglikelihood x 2
  dev_same = -ll_same*2 #deviance is -loglikelihood x 2
}
#compile model details
mod_details = if(!paired_data)
              {
              data.frame(
                modnm = c('grand mean', 'trial mean', 'trial kappa', 'uniform'),
                ll = c(ll_grand_mean, ll_sample_mean, ll_sample_kappa, ll_uniform),
                deviance = c(dev_grand_mean, dev_sample_mean, dev_sample_kappa, dev_uniform),
                rnk = rank(c(dev_grand_mean, dev_sample_mean, dev_sample_kappa, dev_uniform)),#paired differences use fewer observations, don't include in ranking 
                df = c(2, 4, 3, 0)
                         )
              }else
              {
              data.frame(
                modnm = c('pairs diff', 'pairs same'),
                ll = c(ll_diff, ll_same),
                deviance = c(dev_diff, dev_same),
                rnk = rank(c(dev_diff, dev_same)),#paired differences use fewer observations, don't include in ranking 
                df = c(2, 1)
                         )
              }
#set up tests for three hypotheses
lr_tests = if(!paired_data)
            {
            c('uniformity', # data are uniformly or non-uniformly distributed
             'trials_same_mean') # trials share the same mean
            }else
            {'pairs_diff_zero'} #mean of pairs is zero

#function to prepare and perform the LR test
LR_calc = function(tst, mdt)
{
  #collect the deviances for h0 and h1 and calculate the difference in degrees of freedom
  lr_res = 
  with(mdt,
       {
  switch(EXPR = tst,
#test for uniformity (i.e. more comprehensive Rayleigh test)
         uniformity = data.frame(
                       dev0 = deviance[modnm == 'uniform'],#null hypothesis: uniform distribution
                       dev1 = deviance[rnk == 1], #lowest rank is most likely model, any non-uniform distribution
                       d.f. = df[rnk == 1] #uniform has 0 degrees of freedom, test degrees of freedom are best model - 0
                                 ),
#test for unpaired grand mean
         trials_same_mean = data.frame(
                       dev0 = deviance[rnk == 2], #null hypothesis: trials don't differ
                       dev1 = deviance[modnm == 'trial mean'], #within trial obs. share a mean, expect lower deviance with more params
                       d.f. = df[modnm == 'trial mean'] -
                               df[rnk == 2] #2nd best fitting model
                                 ),
#test for nonzero differences of pairs
        pairs_diff_zero = data.frame(
                       dev0 = deviance[modnm == 'pairs same'], #null hypothesis: trials don't differ
                       dev1 = deviance[modnm == 'pairs diff'], #within trial obs. share a mean, expect lower deviance with more params
                       d.f. = df[modnm == 'pairs diff'] -
                               df[modnm == 'pairs same'] #grand mean has half the number of params
                                 ),

         )
       }
  )
  #calculate the change in deviance (chi-squared distributed)
  lr_res = within(lr_res,
         {
           chi_squared = abs(unlist(dev0) - unlist(dev1))
         }
         )
  #calculate the p value
  lr_res = within(lr_res,
         {
           p = pchisq(q = unlist(chi_squared),
                      df = unlist(d.f.),
                      lower.tail = FALSE)
         }
        )
  #adjust for multiple comparisons (3 in this case)
   lr_res = within(lr_res,
                        {
                          p_adjusted = p.adjust(p = p,
                                                method = 'BH',
                                                n = length(p))#could this be flexible?
                        }
                   )
   return(lr_res)
         
}

#Interpret the hypothesis tests depending on size and sign of difference in deviance
H1label = function(tst, d0, d1, pa)
{
  switch(EXPR = tst,
         uniformity = if(d0 > d1 & pa <0.05) # data may be oriented, not significantly oriented, or significantly disoriented
         {'data are significantly oriented'}else
         {'data _are not_ significantly oriented'},
         trials_same_mean = if(d0 > d1 & pa <0.05) # trials significantly differ in mean, not significantly differ in mean, or share a significant mean
         {'trial means differ significantly'}else
         {'trial means _do not_ differ significantly'},
         pairs_diff_zero = if(d0 > d1 & pa <0.05) # trials significantly differ in mean, not significantly differ in mean, or share a significant mean
         {'paired trials differ significantly'}else
         {'paired trials _do not_ differ significantly'}
  )
}

#Perform collect likelihood ratio tests for each hypothesis
all_results = data.frame(tests = lr_tests,
                 t(
                   sapply(X =lr_tests,
                             FUN = LR_calc,
                             mdt = mod_details)
                       )
                      )

#add labels to interpret each test according to p value and sign of comparison
all_results = within(all_results,
                     {
                     result = mapply(tst = tests,
                                     d0 = unlist(dev0),
                                     d1 = unlist(dev1),
                                     pa = p_adjusted,
                                     FUN = H1label
                                     )
                      }
                     )
#add calculated means
Exmu = function(ml){round(ml$mu, 3)}
all_results = within(all_results,
                     {
                     directions = if(!paired_data)
                                    {
                                   c(paste0(Exmu(ml_grand_mean),'°'),
                                    paste0(sapply(ml_sample_mean,Exmu),'°, ', collapse = '') )
                                     }else
                                     {
                                   paste0(Exmu(ml_diff),'°, ', collapse = '')
                                     }
                     }
                     )
#print the results for the user                     
with(all_results,
     {
            print(
              data.frame(chi_squared = round(unlist(chi_squared), 3),
                        d.f = unlist(d.f.),
                         p = round(unlist(p_adjusted), 4),
                         result = result)
            )
      }
)
message(
      '\nestimated difference\n',
      if(!paired_data)
      {
      round(diff(circular(x = sapply(ml_sample_mean,Exmu), 
                          units = angle_unit,
                          rotation = angle_rot)),3) 
      }else
      {
       Exmu(ml_diff)
      }
      )

# Save result -------------------------------------------------------------
#save result
write.table(x = apply(X = all_results,
                      MARGIN = 1:2,
                      FUN = unlist),
            file = file.path(dirname(path_file),#same place as original with
                             paste(sub(pattern='.csv',#similar name
                                       replacement = '',
                                       x = basename(path_file)
                             ),
                             if(!paired_data)
                             {'2-sample A test result.csv'}else
                             {'Paired A test result.csv'}
                             )
            ),
            row.names = FALSE,#rows do not need names
            sep = csv_sep #Use same separator as original csv file
)
