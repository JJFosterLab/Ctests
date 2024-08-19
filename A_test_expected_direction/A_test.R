#FOR A 'CLEAN' RUN, RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2023 06 09
#     MODIFIED:	James Foster              DATE: 2023 06 09
#
#  DESCRIPTION: Loads a text file and performs an "inverse" V test for directedness 
#               towards an expected mean angle. Could be called an "A" test.
#               
#       INPUTS: A ".csv" table with a column of angles ("angle").
#               User should specify test details (line 50).
#               
#      OUTPUTS: Results table (.csv).
#
#	   CHANGES: - 
#
#   REFERENCES: Batschelet E (1981).
#               The Rayleigh test, Chap 4.2, p. 54
#               Chapter 4: Tests for Randomness and Goodness-of-fit
#               In: Circular Statistics in Biology
#               Academic Press (London)
#
#               Aneshansley DJ & Larkin TS (1981). 
#               V-test is not a statistical test of “homeward direction.” 
#                Nature 293, 239–239.
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
expected_mean_angle = 0 # only check for northwards orientation

#Check the operating system and assign a logical flag (T or F)
sys_win <- Sys.info()[['sysname']] == 'Windows'
#On computers set up by JMU Würzburg, use user profile instead of home directory
if(sys_win){
  #get rid of all the backslashes
  ltp <- gsub('\\\\', '/', Sys.getenv('USERPROFILE'))#Why does windows have to make this so difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp <- Sys.getenv('HOME')#Life was easier on Mac
}

# Simulate data (not used) ------------------------------------------------
# n_angles = 10
# angles_sim = suppressWarnings(
#               rvonmises(n = n_angles,
#                                 mu = rcircularuniform(1),
#                                 kappa = A1inv(0.6)
#                       )
#                     )
# sim = data.frame(
#                  angle = round(c(angles_sim)*180/pi)
#                  )
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
plot.circular(x = circular(x = adata$angle, 
                           type = 'angles',
                           unit = 'degrees',
                           template = 'geographics',
                           modulo = '2pi',
                           zero = pi/2,
                           rotation = 'clock'
),
stack = TRUE,
bins = 360/5,
sep = 0.5/dim(adata)[1], # better spacing?
)
arrows.circular(x = circular(x = expected_mean_angle, 
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
arrows.circular(x = circular(x = expected_mean_angle, 
                             type = 'angles',
                             unit = 'degrees',
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = 'clock'
),
y = mean(
  cos((adata$angle - expected_mean_angle) * pi/180), 
  na.rm = T)
)
suppressWarnings(
  {
    lines.circular(x = seq(from = -expected_mean_angle*pi/180+pi/6 +pi/2,#or: bearing - pi/2
                           to = -expected_mean_angle*pi/180-pi/6 +pi/2, #or: bearing + pi/2
                           length.out = 1e3),
                   y = -1+rep( x = 1.644854/(sqrt(2 * length(adata$angle))), times = 1e3),
                   lty = 2
    )
  }
)

# Perform V test ----------------------------------------------------------
#convert angles to circular class
cangs = with(adata, circular(angle, unit = angle_unit) )
emean = circular(expected_mean_angle, unit = angle_unit)

rayv_test = rayleigh.test( x = cangs,
                           mu = emean
                 )
print(rayv_test)


# Calculate log-likelihoods -----------------------------------------------
# The v-test has identified a non-uniform distribution with
# a large component in the expected direction.

#find the maximum likelihood von Mises distribution for this data
ml_vm = mle.vonmises(x = cangs,
                     bias = TRUE) # correct bias

#find the maximum likelihood von Mises distribution around the expected mean for this data
ml_ex = mle.vonmises(x = cangs,
                          mu = emean )

# What is the likelihood of orientation towards the sample mean?
ll_sample_mean = with(ml_vm, #using the maximum likelihood von Mises parameters
                      sum( # add together
                        dvonmises(x = cangs, # probability density for each observed angle
                                    mu = mu, # ML estimated mean
                                    kappa = kappa, # ML estimated concentration
                                    log = TRUE) # on a log scale (i.e. add instead of multiplying)
                      )
)

# What is the likelihood of orientation towards the expected mean?
ll_expect_mean = with(ml_ex, 
                      sum(dvonmises(x = cangs, # probability density for each observed angle
                                    mu = mu, #expected mean
                                    kappa = kappa, # ML estimated concentration
                                    log = TRUE) #return log probability
                      )
)

# Just to confirm the V-test, what is the likelihood of a uniform distribution?
ll_uniform = sum(dvonmises(x = cangs, # probability density for each observed angle
                                mu = circular(0,template = 'none'), #expected mean
                                kappa = 0, # ML estimated concentration
                                log = TRUE) #return log probability
                      )

# Add the two distributions to the figure
aa = circular(x = seq(from = -expected_mean_angle*pi/180+pi,#or: bearing - pi/2
         to = -expected_mean_angle*pi/180-pi, #or: bearing + pi/2
         length.out = 1e3))
  
suppressWarnings(
  {
    with(ml_vm,
    lines.circular(x = aa +pi/2,
                   y = dvonmises(x = -aa,
                                 mu = mu,
                                 kappa = kappa,
                                 log = FALSE) - 1,
                   lty = 3,
                   col = 'cyan3')
    )
    with(ml_ex,
    lines.circular(x = aa+pi/2,
                   y = dvonmises(x = -aa,
                                 mu = mu,
                                 kappa = kappa,
                                 log = FALSE) - 1,
                   lty = 1,
                   col = 'pink')
    )
  }
)
legend(x = 'bottomright',
       inset = -0.01,
       cex = 0.5,
       bg = gray(level = 1,alpha = 0),
       bty = 'n',
       legend = c('vector to expected mean',
                  'p(>V) < 0.05',
                  'probability density: sample mean',
                  'probability density: expected mean'),
       col = c('black',
               'black',
               'cyan4',
               'pink'),
       lty = c(1,
               2,
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
dev_sample_mean = -ll_sample_mean*2 #deviance is -loglikelihood x 2
dev_expect_mean = -ll_expect_mean*2 #deviance is -loglikelihood x 2
dev_uniform = -ll_uniform*2 #deviance is -loglikelihood x 2

#compare expected mean with the 2nd most likely distribution
if(ll_sample_mean > ll_uniform)
{
#calculate the likelihood ratio chi-squared statistic
lr_stat = dev_expect_mean - dev_sample_mean # should be positive, sample mean always more likely
}else
{ # N.B. This is unnecessary, in this case all distributions have the same likelihood
  #if uniform is more likely, compare expected against uniform
  lr_stat = dev_uniform -  dev_expect_mean# should be positive
}
#Does the estimated mean significantly reduce deviance?
pLR = pchisq(q = lr_stat, # likelihood ratio chi-squared
       df = 1, # one parameter difference
       lower.tail = !(ll_sample_mean > ll_uniform)) # p(larger deviance): null hypothesis, expected mean is true mean
print(data.frame(chi_squared = round(lr_stat, 3),
        d.f. = 1,
        p = round(pLR, 4),
        h2 = if(pLR <0.05)
          {'the sample mean differs significantly from the expected mean'}else
          {'the sample mean _does not_ differ significantly from the expected mean'})
      )

# Save result -------------------------------------------------------------
co_rayv_test = capture.output(rayv_test)
#save result
write.table(x = with(rayv_test,#save relevant info in a table
                     cbind(data = path_file,
                           test = paste('(V test)' , co_rayv_test[2]),
                           h0 = 'the distribution of angles is uniform',
                           h1 = 'angles are from a unimodal von Mises distribution centred on the expected mean angle',
                           h2 = 'angles are from a unimodal von Mises distribution _not_ centred on the sample mean angle',
                           R_statistic = statistic,
                           expected_mean = mu,
                           p_h0 = p.value,
                           LogLik_h1 = ll_expect_mean,
                           LogLik_h2 = ll_sample_mean,
                           chi_sq = lr_stat,
                           df = 1,
                           p_h2 = pLR
                     )
),
file = file.path(dirname(path_file),#same place as original with
                 paste(sub(pattern='.csv',#similar name
                           replacement = '',
                           x = basename(path_file)
                 ),
                 'A test result.csv'#ending in mean test
                 )
),
row.names = FALSE,#rows do not need names
sep = csv_sep #Use same separator as original csv file

)

