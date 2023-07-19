#FOR A 'CLEAN' RUN, RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2023 07 18
#     MODIFIED:	James Foster              DATE: 2023 07 19
#
#  DESCRIPTION: Loads a text file with axial data and finds and plots the assumed
#               means and standard deviations of a symmetrical axial distribution
#               as calculated by doubling the data. N.B. If the two means are not 
#               equally concentrated, equal proportions of the distribution, or 
#               separated by exactly 180°, CircMLE may be more appropriate.
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
#    EXAMPLES:  Fill out user input (lines 50-55), then press ctrl+shift+s to run
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Read in data   +
#- Test with simulated data +
#- Maximum likelihood version
#- Bootstrapped confidence intervals

# Useful functions --------------------------------------------------------
# . Load package ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    require(circular)#package for handling cirular data
    require(CircStats)#package for circular hypothesis tests
  }
)



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

# Simulate data (not used) ------------------------------------------------
# n_angles = 20
# suppressWarnings( {mu1 = rcircularuniform(1) })
# angles_sim = suppressWarnings(
#               rmixedvonmises(n = n_angles,
#                                 mu1 = mu1,
#                                 kappa1 = A1inv(0.9),
#                                 mu2 = mu1 + pi,
#                                 kappa2 = A1inv(0.9),
#                                 prop = 0.5
#                       )
#                     )
# sim = data.frame(
#                  angle = round(c(angles_sim)*180/pi)
#                  )
# write.table(x = sim,
#             file = file.path(ltp,'Documents', "simulated_axial.csv"),
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



# Plot the data -----------------------------------------------------------
par(mar = c(0,0,0,0)) # make space for plotting outside of the circle

plot.circular(x = circular(x = adata$angle, 
                           type = 'angles',
                           unit = 'degrees',
                           # template = 'geographics',
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


# Calculate summary statistics for doubled angles -------------------------


doubled_mean = mean.circular(x = 2* circular(x = adata$angle, 
                                             type = 'angles',
                                             unit = 'degrees',
                                             template = 'geographics',
                                             modulo = '2pi',
                                             zero = pi/2,
                                             rotation = 'clock'
)
)

doubled_sd = circular::deg(
                            sd.circular(x = 2* circular(x = adata$angle, 
                                                       type = 'angles',
                                                       unit = 'degrees',
                                                       template = 'geographics',
                                                       modulo = '2pi',
                                                       zero = pi/2,
                                                       rotation = 'clock'
                                                        )
                                        )
)

# Plot the summary statistics ---------------------------------------------

arrows.circular(x = circular(x = doubled_mean/2, 
                             type = 'angles',
                             unit = 'degrees',
                             #template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = 'clock'
),
y = 1,
lwd = 5, 
col = lines_col,
length = 0
)
arrows.circular(x = circular(x = doubled_mean/2 +180, 
                             type = 'angles',
                             unit = 'degrees',
                             #template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = 'clock'
),
y = 1,
lwd = 5, 
col = lines_col,
length = 0
)


# Choose the correct plotting direction for errorbars ---------------------
ci = c(-1,1)*doubled_sd/2 + doubled_mean/2

if(sign(max(ci) - doubled_mean/2) == -1 | 
   sign(doubled_mean/2 - min(ci)) == -1 )
  {#check if mean is right or left of max CI
  #otherwise, plot in the reverse direction
  sq = seq( from = as.numeric(min(ci)),
            to = -as.numeric(360 - max(ci)), 
            length.out = 1e3)
  }else{ #if it is "little" plot from min to max
              sq = seq( from = as.numeric(min(ci)),
                        to = as.numeric(max(ci)), 
                        length.out = 1e3)
  }

#draw a line in coordinates of rcos(x) rsin(y)
#any other line arguments (i.e. col or lwd) input as "..."
lines.circular(x = circular(x = sq,
                            type = 'angles',
                            unit = 'degrees',
                            template = 'geographics',
                            modulo = '2pi',
                            zero = pi/2,
                            rotation = 'clock'
                        ), 
               y = rep(x = errbar_dist,
                   times = 1e3), 
               col = lines_col,
               lwd = 3)

lines.circular(x = circular(x = sq+180, #2nd mean
                            type = 'angles',
                            unit = 'degrees',
                            template = 'geographics',
                            modulo = '2pi',
                            zero = pi/2,
                            rotation = 'clock'
                        ), 
               y = rep(x = errbar_dist,
                   times = 1e3), 
               col = lines_col,
               lwd = 3)

