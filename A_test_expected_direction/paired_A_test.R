#FOR A 'CLEAN' RUN, RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 03 22
#     MODIFIED:	James Foster              DATE: 2024 03 25
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
#
#   REFERENCES: Batschelet E (1981).
#               The Rayleigh test, Chap 4.2, p. 54
#               Chapter 4: Tests for Randomness and Goodness-of-fit
#               In: Circular Statistics in Biology
#               Academic Press (London)
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
#- Optimisation method for pairs (ML too biased?)
#- Tie-breaking approach
#- Comment in detail

# Useful functions --------------------------------------------------------
# . Load package ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    require(circular)#package for handling circular data
    require(CircStats)#package for circular hypothesis tests
  }
)




# Input Variables ----------------------------------------------------------

#  .  User input -----------------------------------------------------------
csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"
angle_name = "angles" #The title of the column with angles; NO SPACES PLEASE
angle_unit = "degrees" # "degrees" or "radians"
expected_mean_angle = 'angle_1' # 2nd set should be same as 1st mean

#Check the operating system and assign a logical flag (T or F)
sys_win <- Sys.info()[['sysname']] == 'Windows'
#On computers set up by JMU Würzburg, use user profile instead of home directory
if(sys_win){
  #get rid of all the backslashes
  ltp <- gsub('\\\\', '/', Sys.getenv('USERPROFILE'))#Why does windows have to make this so difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp <- Sys.getenv('HOME')#Life was easier on Mac
}

  # # Simulate data (not used) ------------------------------------------------
  # n_angles = 44
  # mu_offset = circular(x = rad(65) )
  # # minimum discriminable angle appears to be approx 35°
  # kappa_both = A1inv(0.7) #concentration around each trial mean
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
  #                            kappa = kappa_both
  #                    )
  #                  ))*180/pi),#convert to angles and round to nearest degree
  #                  angle_2 = round(c(suppressWarnings( #rvonmises converts to circular and warns
  #                              sapply(X =mu1_sim + mu_offset,# true difference,
  #                                     FUN = rvonmises,
  #                                     n = 1,
  #                                     kappa = kappa_both
  #                              )
  #                  ))*180/pi) #convert to angles and round to nearest degree
  #                  )
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
{stop('pairs not found, please check data are in two columns')}

ref_angle = switch(expected_mean_angle,
                   `angle_1` = mean.circular(x = circular(x = adata$angle_1, 
                                                          type = 'angles',
                                                          unit = 'degrees',
                                                          template = 'geographics',
                                                          modulo = '2pi',
                                                          zero = pi/2,
                                                          rotation = 'clock')
                                             ),
                   `angle_2` = mean.circular(x = circular(x = adata$angle_2, 
                                                          type = 'angles',
                                                          unit = 'degrees',
                                                          template = 'geographics',
                                                          modulo = '2pi',
                                                          zero = pi/2,
                                                          rotation = 'clock')
                                             ),
                   mean.circular(x = circular(x = adata$angle_1, 
                                                          type = 'angles',
                                                          unit = 'degrees',
                                                          template = 'geographics',
                                                          modulo = '2pi',
                                                          zero = pi/2,
                                                          rotation = 'clock')
                                             )
                   
                   )

# Plot the data -----------------------------------------------------------

par(mar =rep(0,4))
plot.circular(x = circular(x = adata$angle_1, 
                           type = 'angles',
                           unit = 'degrees',
                           template = 'geographics',
                           modulo = '2pi',
                           zero = pi/2,
                           rotation = 'clock'
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
                           rotation = 'clock'
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
                             unit = 'degrees',
                             #template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = 'clock'
),
y = 1,
lwd = 5, 
col = rgb(0,0,0,0.1),
length = 0
)
arrows.circular(x = circular(x = ref_angle, 
                             type = 'angles',
                             unit = 'degrees',
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = 'clock'
                            ),
                y = mean(x = 
                  cos(x = 
                    switch(EXPR = expected_mean_angle, 
                          `angle_1` = adata$angle_2 - ref_angle,
                           `angle_2` = adata$angle_1 - ref_angle, 
                               adata$angle_2 - ref_angle)* pi/180) , 
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
cangs = with(adata, circular(x = switch(EXPR = expected_mean_angle, 
                                    `angle_1` = angle_2,
                                    `angle_2` = angle_1 , 
                                    angle_2), 
                             unit = angle_unit) )
emean = circular(ref_angle, unit = angle_unit)

rayv_test = rayleigh.test( x = cangs,
                           mu = emean
)
print(rayv_test)


# Calculate log-likelihoods -----------------------------------------------
# The v-test has identified a non-uniform distribution with
# a large component in the expected direction.

#find the maximum likelihood von Mises distribution for the full dataset

#with fewer degrees of freedom log likelihood is lower than dividing data into trials
ml_up = with(adata,
             mle.vonmises(x = circular(c(angle_1,
                                         angle_2),
                                         units = angle_unit
                                       ),
                                  bias = TRUE)# correct bias
              )
#find the maximum likelihood von Mises distribution for each trial

#with fewer degrees of freedom log likelihood is lower than dividing data into pairs
ml_vm = with(adata,
             apply(X =  cbind(circular(angle_1,units = angle_unit),
                              circular(angle_2,units = angle_unit)),
                   MARGIN = 2,
                   FUN = function(angs)
                     {
                     mle.vonmises(x = circular(angs,units = angle_unit),
                                  bias = TRUE)
                     }
                   ) # correct bias
              )

#find the maximum likelihood von Mises distribution around the mean of both heading
# this is the expectation that each pair of observations are drawn from the same distrution
# N.B. that is even if each pair comes from a different distribution
ml_ex = with(adata, 
             apply(X =  cbind(circular(angle_1,units = angle_unit),
                              circular(angle_2,units = angle_unit)), 
                   MARGIN = 1, 
                   FUN = function(angs)
                   {
                     mle.vonmises(x = if(identical(angs[1],angs[2]))
                       {circular(x = as.numeric(angs),
                                 units = angle_unit)+ #tie break by separating angles 
                                 circular(x = c(-1,1)*0.5, # by 1 degree
                                          units = angle_unit) #TODO better tie break,
                        }else{circular(x = angs,
                                       units = angle_unit)
                        },
                        bias = TRUE)
                   }
             ) # correct bias
)

# What is the likelihood of orientation towards the grand mean?
ll_grand_mean = with(ml_up, #using the maximum likelihood von Mises parameters
                      sum( # add together
                        dvonmises(x = circular(with(adata,c(angle_1,angle_2)),
                                               units = angle_unit), # probability density for each observed angle
                                  mu = mu, # ML estimated mean
                                  kappa = kappa, # ML estimated concentration
                                  log = TRUE) # on a log scale (i.e. add instead of multiplying)
                      )
                    ) 
# What is the likelihood of orientation towards each trial's mean?
ll_sample_mean = with(ml_vm[[1]], #using the maximum likelihood von Mises parameters
                      sum( # add together
                        dvonmises(x = circular(adata$angle_1,
                                               units = angle_unit), # probability density for each observed angle
                                  mu = mu, # ML estimated mean
                                  kappa = kappa, # ML estimated concentration
                                  log = TRUE) # on a log scale (i.e. add instead of multiplying)
                      )
                    ) + 
  with(ml_vm[[2]], #using the maximum likelihood von Mises parameters
                      sum( # add together
                        dvonmises(x = circular(adata$angle_2,
                                               units = angle_unit), # probability density for each observed angle
                                  mu = mu, # ML estimated mean
                                  kappa = kappa, # ML estimated concentration
                                  log = TRUE) # on a log scale (i.e. add instead of multiplying)
                      )
) 

# What is the likelihood of orientation towards each pair's expected mean?
ll_expect_mean = sum(
                  mapply(ml = ml_ex,
                        a1 = adata$angle_1,
                        a2 = adata$angle_2,
                        FUN = function(ml, a1, a2)
                            {
                            sum(
                              ifelse(test = a1 == a2,
                                     yes = dvonmises(x = circular(x = c(a1+1, a2-1), #TODO better tie break
                                                                  units = angle_unit),
                                                     mu = ml$mu, #expected mean
                                                     kappa = ml$kappa, # ML estimated concentration
                                                     log = TRUE), #return log probability
                                     no = dvonmises(x = circular(x = c(a1, a2),
                                                     units = angle_unit),
                                        mu = ml$mu, #expected mean
                                        kappa = ml$kappa, # ML estimated concentration
                                        log = TRUE) #return log probability
                                    )
                              )
                            }
                      )
                    )
  

# Just to confirm the V-test, what is the likelihood of a uniform distribution?
ll_uniform = with(adata,
                  sum(log(dcircularuniform(x = circular(x = c(angle_1, angle_2),
                                                   units = angle_unit))) #return log probability
                  )
)

# Add ML models to the figure ---------------------------------------------


# Add the two distributions to the figure
# aa = circular(x = seq(from = -ref_angle*pi/180+pi,#or: bearing - pi/2
#                       to = -ref_angle*pi/180-pi, #or: bearing + pi/2
#                       length.out = 1e3))
aa = circular(x = seq(from = -180,#or: bearing - pi/2
                      to = 180, #or: bearing + pi/2
                      length.out = 1e3),
              unit = angle_unit)

#probability density across all pairs
pair_dens = sapply(X = ml_ex, 
                   FUN = function(params)
                     {
                     with(params,
                          {
                     dvonmises(x = aa,
                               mu = mu,
                               kappa = kappa,
                               log = FALSE)
                          }
                     )
                   }
                   )
#take maximum (representative for pointwise estimates of each pair)
max_pair_dens = apply(X = pair_dens, 
                       MARGIN = 1,
                       FUN = max
                       )
#probability density grand mean
with(ml_up,
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
with(ml_vm[[1]],
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
with(ml_vm[[2]],
     lines.circular(x = 90-aa,
                    y = dvonmises(x = aa,
                                  mu = mu,
                                  kappa = kappa,
                                  log = FALSE) - 1,
                    lty = 3,
                    lwd = 2,
                    col = 'blue3')
)
#maximum probability density for each pair
lines.circular(x = 90-aa,
               y = max_pair_dens -1, #give a rough impression of probability density corresponding to each point
               lty = 3,
               lwd = 2,
               col = 'pink')
#mean vector for each pair
for(i in 1:length(ml_ex))
{
  with(ml_ex[[i]],
       arrows.circular(
                       x = circular(mu,
                                    type = 'angles',
                                    unit = 'degrees',
                                    template = 'geographics',
                                    modulo = '2pi',
                                    zero = pi/2,
                                    rotation = 'clock'),
                      y = A1(kappa),
                      length = 0,
                      lty = 1,
                      col = adjustcolor('pink', alpha.f = 0.9) )
  )
}

legend(x = 'bottomright',
       inset = c(0.01,0.02),
       cex = 0.4,
       bg = gray(level = 1,alpha = 0),
       bty = 'n',
       legend = c('vector to expected mean',
                  'p(>V) < 0.05',
                  'prob. dens.: grand mean',
                  'prob. dens.: sample mean 1',
                  'prob. den.: sample mean 2',
                  'max. prob. dens.: pairs same mean',
                  'mean vector: pairs same mean'),
       col = c('black',
               'black',
               'green4',
               'cyan4',
               'blue3',
               'pink',
               'pink'),
       lty = c(1,
               2,
               3,
               3,
               3,
               3,
               1)
)


#N.B. unless the two means are identical, the likelihood of the expected mean
#is always lower than the maximum likelihood mean (it has _maximum_ likelihood).
#Our question is, it significantly more likely than the distribution around the 
#expected mean, for which we need to fit one fewer parameter (we already know the mean).


# Perform likelihood ratio test -------------------------------------------
#The likelihood ratio test compares models by estimating the likelihood gained
#for each additional parameter.
#https://en.wikipedia.org/wiki/Likelihood-ratio_test
# In our case, the ML von Mises requires two parameters (mean and concentration), 
# while the von Mises around the expected mean requires one (concentration).
# How much does the extra parameter increase the likelihood?

#convert likelihood to deviance (more likely = lower deviance)
dev_grand_mean = -ll_grand_mean*2 #deviance is -loglikelihood x 2
dev_sample_mean = -ll_sample_mean*2 #deviance is -loglikelihood x 2
dev_expect_mean = -ll_expect_mean*2 #deviance is -loglikelihood x 2
dev_uniform = -ll_uniform*2 #deviance is -loglikelihood x 2
#compile model details
mod_details = data.frame(
                modnm = c('grand mean', 'trial mean', 'pair mean', 'uniform'),
                ll = c(ll_grand_mean, ll_sample_mean, ll_expect_mean, ll_uniform),
                deviance = c(dev_grand_mean, dev_sample_mean, dev_expect_mean, dev_uniform),
                rnk = rank(c(dev_grand_mean, dev_sample_mean, dev_expect_mean, dev_uniform)),
                df = c(2, 4, length(adata$angle_1)*2, 0)
                         )
#set up tests for three hypotheses
lr_tests = c('uniformity', # data are uniformly or non-uniformly distributed
             'pairs_same_mean', # pairs share the same mean (independently of trials)
             'trials_same_mean') # trials share the same mean
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
#test for pairing (pairs have same mean)
         pairs_same_mean = data.frame(
                       dev0 = deviance[modnm == 'trial mean'], #null hypothesis: pairs differ across trials
                       dev1 = deviance[modnm == 'pair mean'], #pairs share a mean, expect lower deviance with more params
                       d.f. = df[modnm == 'pair mean'] -
                               df[modnm == 'trial mean'] #trial comparison has only two means, fewer params
                                 ),
#test for unpaired grand mean
         trials_same_mean = data.frame(
                       dev0 = deviance[modnm == 'grand mean'], #null hypothesis: trials don't differ
                       dev1 = deviance[modnm == 'trial mean'], #within trial obs. share a mean, expect lower deviance with more params
                       d.f. = df[modnm == 'trial mean'] -
                               df[modnm == 'grand mean'] #grand mean has half the number of params
                                 )
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
                                                n = 3)#could this be flexible?
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
         pairs_same_mean = if(d0 > d1 & pa <0.05) # pairs may share a mean, share no significant mean, or significantly differ in mean
         {'pairs are significantly oriented in the same direction'}else
         {'pairs _are not_ oriented in the same direction'},
         trials_same_mean = if(d0 > d1 & pa <0.05) # trials significantly differ in mean, not significantly differ in mean, or share a significant mean
         {'trial means differ significantly'}else
         {'trial means _do not_ differ significantly'}
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
                     directions = c(paste0(round(ml_up$mu,3),'°'),
                                    paste0(sapply(ml_ex,Exmu),'°, ', collapse = ''),
                                    paste0(sapply(ml_vm,Exmu),'°, ', collapse = '')
                                    )
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
message(c('trial means\n', 
      paste0(sapply(ml_vm,Exmu),'°, ', collapse = ''), 
      '\nestimated difference\n', 
      diff(sapply(ml_vm,Exmu))
      )
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
                             'Paired A test result.csv'#ending in mean test
                             )
            ),
            row.names = FALSE,#rows do not need names
            sep = csv_sep #Use same separator as original csv file
)


  # # Random Effects model WIP ------------------------------------------------
  # 
  # REcirc = function(a1, 
  #                   a2, 
  #                   au = 'degrees',
  #                   trc = 0
  # )
  # {
  #   Get_mu = function(m){m$mu}
  #   Get_kappa = function(m){m$kappa}
  #   grand_mod = mle.vonmises(circular(x = c(a1, a2), units = au))  
  #   trial_mod = apply(X =  cbind(circular(a1,units = au),
  #                                circular(a2,units = au)),
  #                     MARGIN = 2,
  #                     FUN = function(angs)
  #                     {
  #                       mle.vonmises(x = circular(angs,units = au),
  #                                    bias = TRUE)
  #                     }
  #   )
  #   pair_mod = apply(X =  cbind(circular(a1,units = au),
  #                               circular(a2,units = au)), 
  #                    MARGIN = 1, 
  #                    FUN = function(angs)
  #                    {
  #                      mle.vonmises(x = if(identical(angs[1],angs[2]))
  #                      {circular(x = as.numeric(angs),
  #                                units = au)+ #tie break by separating angles 
  #                          circular(x = c(-1,1)*0.5, # by 1 degree
  #                                   units = au) #TODO better tie break,
  #                      }else{circular(x = angs,
  #                                     units = au)
  #                      },
  #                      bias = TRUE)
  #                    }
  #   )
  #   indiv_mod = mle.vonmises(x = circular(x = 
  #                                           sapply(X = pair_mod, 
  #                                                  FUN = Get_mu),
  #                                         units = au)
  #   )
  #   grand_mu = Get_mu(grand_mod)
  #   grand_logkappa = log(Get_kappa(grand_mod))
  #   trial_mu = sapply(X = trial_mod, FUN = Get_mu) - grand_mu
  #   trial_logkappa = log(sapply(X = trial_mod, FUN = Get_kappa)) - grand_logkappa
  #   pair_mu = sapply(X = pair_mod, FUN = Get_mu) - grand_mu
  #   pair_logkappa = log(sapply(X = pair_mod, FUN = Get_kappa)) - grand_logkappa
  #   indiv_mu = Get_mu(indiv_mod)
  #   indiv_logkappa = log(Get_kappa(indiv_mod))
  #   
  #   start_par = c(as.numeric(grand_mu) %% 360,#1 
  #                 grand_logkappa,#2
  #                 as.numeric(trial_mu) %% 360, #34
  #                 trial_logkappa,#56 
  #                 as.numeric(pair_mu) %% 360, #6+i
  #                 pair_logkappa,#6+N+i 
  #                 indiv_logkappa #6+2N + 1
  #   )
  #   start_par[is.infinite(start_par)] = -10
  #   LLopt = function(a1, a2, au = au, prms)
  #   {
  #     mm = prms[c(1,
  #                 3,
  #                 4,
  #                 6+1:length(a1))]
  #     if(any(mm > 360-1e-16) |
  #        any(mm < 0))
  #     {
  #       neg_ll = 1e9
  #     }else
  #     {
  #       ll = sum(
  #         sapply(X = 1:length(a1),
  #                au = au,
  #                prms = prms,
  #                FUN =  function(i, prms, au)
  #                {
  #                  dvonmises(x = circular(a1[i], units = au), 
  #                            mu = circular(prms[1] + # grand mean 
  #                                            prms[3] + #trial 1 mean
  #                                            prms[6+i], units = au), #pair i mean
  #                            kappa = exp(prms[2] + # grand kappa 
  #                                          prms[5] + #trial 1 kappa
  #                                          prms[6+length(a1)+i]), #pair i kappa
  #                            log = TRUE
  #                  )
  #                  + dvonmises(x = circular(a2[i], units = au), 
  #                              mu = circular(prms[2] + # grand mean 
  #                                              prms[4] + #trial 2 mean
  #                                              prms[6+i], units = au), #pair i mean
  #                              kappa = exp(prms[2] + # grand kappa 
  #                                            prms[6] + #trial 1 kappa
  #                                            prms[6+length(a1)+i]), #pair i kappa
  #                              log = TRUE
  #                  )
  #                  
  #                }
  #         )
  #       )+ 
  #         sum(
  #           dvonmises(x = circular(prms[6+1:length(a1)], 
  #                                  units = au),
  #                     mu = circular(0),
  #                     kappa = exp(prms[6+2*length(a1)+1]),
  #                     log = TRUE
  #           )#indiv_kappa
  #         )
  #       neg_ll = -ll
  #       # if(is.infinite(neg_ll)){neg_ll = 1e9}
  #     }
  #     return(neg_ll)
  #   }
  #   # LLopt(a1, a2, au = au, prms = start_par)
  #   opt = optim(par = start_par,
  #               fn = LLopt,
  #               a1 = a1,
  #               a2 = a2,
  #               au = au,
  #               control = list(trace = trc)
  #   )
  #   return(opt)
  # }
  # 
  # opt_remod = REcirc(a1 = adata$angle_1, a2 = adata$angle_2)
  # ll_remod = with(adata,
  #                 {
  #                   sum(
  #                     sapply(X = 1:length(angle_1),
  #                            au = angle_unit,
  #                            prms = opt_remod$par,
  #                            FUN =  function(i, prms, au)
  #                            {
  #                              dvonmises(x = circular(angle_1[i], units = au), 
  #                                        mu = circular(prms[1] + # grand mean 
  #                                                        prms[3] + #trial 1 mean
  #                                                        prms[6+i], units = au), #pair i mean
  #                                        kappa = exp(prms[2] + # grand kappa 
  #                                                      prms[5] + #trial 1 kappa
  #                                                      prms[6+length(angle_1)+i]), #pair i kappa
  #                                        log = TRUE
  #                              )
  #                              + dvonmises(x = circular(angle_2[i], units = au), 
  #                                          mu = circular(prms[2] + # grand mean 
  #                                                          prms[4] + #trial 2 mean
  #                                                          prms[6+i], units = au), #pair i mean
  #                                          kappa = exp(prms[2] + # grand kappa 
  #                                                        prms[6] + #trial 1 kappa
  #                                                        prms[6+length(angle_1)+i]), #pair i kappa
  #                                          log = TRUE
  #                              )
  #                              
  #                            }
  #                     )
  #                   )
  #                 }
  # )

