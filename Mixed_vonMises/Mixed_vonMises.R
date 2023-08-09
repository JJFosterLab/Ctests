#FOR A 'CLEAN' RUN, RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2023 08 09
#     MODIFIED:	James Foster              DATE: 2023 08 09
#
#  DESCRIPTION: Generates axial data and finds and plots the assumed
#               means and mean vectors of a symmetrical axial distribution.
#               
#       INPUTS: A ".csv" table with a column of angles ("angle").
#               User should specify test details (line 50).
#               
#      OUTPUTS: Results table (.csv).
#
#	   CHANGES: - Suppressed package loading messages (upset users)
#             - Plot critical mean in expected direction
#             - Handle BOM in input
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
lines_col = adjustcolor('gray20', alpha.f = 0.5) # colour the mean and sd are plotted in
errbar_dist = 0.15 # distance between errorbars and the outside of the circle (in "lines of text")

#Check the operating system and assign a logical flag (T or F)
sys_win = Sys.info()[['sysname']] == 'Windows'
#On computers set up by JMU Würzburg, use user profile instead of home directory
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
  }
)

# Generate samples from circular distributions ----------------------------


# . Generate a unimodal sample --------------------------------------------

#set n
n_angles = 20
suppressWarnings( {mu1 = rcircularuniform(1) })
kappa1 = A1inv(0.5)
angles_sim1 = suppressWarnings(
              rvonmises(n = n_angles,
                                mu = mu1,
                                kappa = kappa1
                      )
                    )
plot(angles_sim1,
     template = 'geographics')
rho.circular(angles_sim1)


# . Generate complementary unimodal sample --------------------------------
angles_sim2 = suppressWarnings(
  rvonmises(n = n_angles,
            mu = mu1+pi,
            kappa = kappa1
  )
)
rho.circular(angles_sim2)

angles_sim = c(angles_sim1, angles_sim2)

plot(angles_sim,
     template = 'geographics')
rho.circular(angles_sim*2)
# angles_sim = suppressWarnings(
#               rmixedvonmises(n = n_angles,
#                                 mu1 = mu1,
#                                 kappa1 = A1inv(0.9),
#                                 mu2 = mu1 + pi,
#                                 kappa2 = A1inv(0.9),
#                                 prop = 0.5
#                       )
                    # )
# sim = data.frame(
#                  angle = round(c(angles_sim)*180/pi)
#                  )
# write.table(x = sim,
#             file = file.path(ltp,'Documents', "simulated_axial.csv"),
#             sep = csv_sep,
#             row.names = FALSE
#             )

# WIP!!!! -----------------------------------------------------------------



