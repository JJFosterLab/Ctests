#FOR A 'CLEAN' RUN, RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2023 08 09
#     MODIFIED:	James Foster              DATE: 2023 08 09
#
#  DESCRIPTION: Generates axial data and finds and plots the assumed
#               means and mean vectors of a symmetrical axial distribution.
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

#  .  User input -----------------------------------------------------------
csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"
angle_name = "angles" #The title of the column with angles; NO SPACES PLEASE
angle_unit = "degrees" # "degrees" or "radians"
point_col = 'darkblue' # colour the datapoints are plotted in
lines_col = adjustcolor('gray20', alpha.f = 0.5) # colour the original parameters are plotted in
est_col = 'magenta3' # colour the fits are plotted in
errbar_dist = 0.15 # distance between errorbars and the outside of the circle (in "lines of text")

#Check the operating system and assign a logical flag (T or F)
sys_win = Sys.info()[['sysname']] == 'Windows'
#On Windows computers, use user profile instead of home directory
if(sys_win){
  #get rid of all the backslashes
  ltp = gsub('\\\\', '/', Sys.getenv('USERPROFILE'))#Why does windows have to make this so difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp = Sys.getenv('HOME')#Life was easier on Mac
}

# Useful functions --------------------------------------------------------
# . Load package ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    require(circular)#package for handling cirular data
    require(CircStats)#package for circular hypothesis tests
    require(CircMLE)#package for fitting circular distributions
  }
)

# Generate samples from circular distributions ----------------------------


# . Generate a unimodal sample --------------------------------------------

#set n
n_angles = 15
#choose a random mean direction
suppressWarnings( {mu1 = rcircularuniform(1) })
#choose a mean vector length
mvl1 = 0.72
# N.B. for this sample size p = 0.05 for a mean vector length of:
sqrt(-log(0.05)/n_angles)
# [1] 0.3870228

kappa1 = A1inv(mvl1)

angles_sim1 = suppressWarnings(
              rvonmises(n = n_angles,
                                mu = mu1,
                                kappa = kappa1
                      )
                    )

par(mar = c(0,0,0,0)) # make space for plotting outside of the circle

plot.circular(x = circular(x = deg(angles_sim1), 
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
              ylim = c(-1,1)*1.2,# make space for plotting outside of the circle
)
arrows.circular(x = mu1,
                y = mvl1,
                lwd = 5, 
                col = lines_col,
                length = 0,
)
arrows.circular(x = mean(circular(x = deg(angles_sim1),
                             type = 'angles',
                             unit = 'degrees',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = 'clock'
                            )
                         ),
                y = rho.circular(angles_sim1),
                lwd = 2, 
                col = est_col,
                length = 0.1,
)
lines.circular(x = seq(from = -pi, 
                       to = pi, 
                       length.out = 1e3),
               y = rep(x = sqrt(-log(0.05)/n_angles),
                       times = 1e3)-1,
               lty = 3)

rayleigh.test(angles_sim1)


# . Generate complementary unimodal sample --------------------------------
angles_sim2 = suppressWarnings(
  rvonmises(n = n_angles,
            mu = mu1+pi,
            kappa = kappa1
  )
)


plot.circular(x = circular(x = deg(angles_sim2), 
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
ylim = c(-1,1)*1.2,# make space for plotting outside of the circle
)
arrows.circular(x = mu1+pi,
                y = mvl1,
                lwd = 5, 
                col = lines_col,
                length = 0,
)
arrows.circular(x = mean(circular(x = deg(angles_sim2),
                                  type = 'angles',
                                  unit = 'degrees',
                                  modulo = '2pi',
                                  zero = pi/2,
                                  rotation = 'clock'
)
),
y = rho.circular(angles_sim2),
lwd = 2, 
col = est_col,
length = 0.1,
)
lines.circular(x = seq(from = -pi, 
                       to = pi, 
                       length.out = 1e3),
               y = rep(x = sqrt(-log(0.05)/n_angles),
                       times = 1e3)-1,
               lty = 3)

rayleigh.test(angles_sim2)


# Combine the two samples -------------------------------------------------



angles_sim = c(angles_sim1, angles_sim2)



plot.circular(x = circular(x = deg(angles_sim), 
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
ylim = c(-1,1)*1.2,# make space for plotting outside of the circle
)
arrows.circular(x = mu1+c(0, pi),
                y = c(mvl1,mvl1),
                lwd = 5, 
                col = lines_col,
                length = 0,
)
arrows.circular(x = mean(circular(x = deg(angles_sim*2),
                                  type = 'angles',
                                  unit = 'degrees',
                                  modulo = '2pi',
                                  zero = pi/2,
                                  rotation = 'clock'
                                  )
                          )*2 + 
                  circular(x = c(0,180),
                           type = 'angles',
                           unit = 'degrees',
                           modulo = '2pi',
                           zero = pi/2,
                           rotation = 'clock'
                  ),
                          y = rep(x = rho.circular(angles_sim*2),
                                  times = 2),
                          lwd = 2, 
                          col = est_col,
                          length = 0.1,
                )
lines.circular(x = seq(from = -pi, 
                       to = pi, 
                       length.out = 1e3),
               y = rep(x = sqrt(-log(0.05)/(n_angles*2) ),
                       times = 1e3)-1,
               lty = 3)

rayleigh.test(angles_sim*2)


# Find Maximum Likelihood Model -------------------------------------------

cmle = circ_mle(angles_sim)
# plot_circMLE(angles_sim, cmle)
est_kappa = cmle$results$k1[1]
est_mean = cmle$results$q1[1]
est_prop = cmle$results$lamda[1]
A1(cmle$results$k1[1])



plot.circular(x = circular(x = deg(angles_sim), 
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
ylim = c(-1,1)*1.2,# make space for plotting outside of the circle
)
arrows.circular(x = mu1+c(0, pi),
                y = c(mvl1,mvl1),
                lwd = 5, 
                col = lines_col,
                length = 0,
)
arrows.circular(x = circular(x = deg(est_mean),
                             type = 'angles',
                             unit = 'degrees',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = 'clock'
) + 
  circular(x = c(0,180),
           type = 'angles',
           unit = 'degrees',
           modulo = '2pi',
           zero = pi/2,
           rotation = 'clock'
  ),
y = rep(x = A1(est_kappa),
        times = 2),
lwd = 2, 
col = est_col,
length = 0.1,
)
lines.circular(x = seq(from = -pi, 
                       to = pi, 
                       length.out = 1e3),
               y = rep(x = sqrt(-log(0.05)/(n_angles*2) ),
                       times = 1e3)-1,
               lty = 3)
# What if the different means account for different proportions? ----------



# WIP!!!! -----------------------------------------------------------------



