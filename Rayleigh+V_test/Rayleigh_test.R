#FOR A 'CLEAN' RUN, RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2021 07 19
#     MODIFIED:	James Foster              DATE: 2021 07 19
#
#  DESCRIPTION: Loads a text file and performs a Rayleigh test for uniformity.
#               
#       INPUTS: A ".csv" table with a column of angles ("angle").
#               User should specify test details (line 50).
#               
#      OUTPUTS: Results table (.csv).
#
#	   CHANGES: - 
#             - 
#             - 
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
#- Perform test     
#- Test with simulated data 
#- Save results

# Useful functions --------------------------------------------------------
# . Load package ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
require(circular)#package for handling cirular data
require(CircStats)#package for circular hypothesis tests







# Input Variables ----------------------------------------------------------

#  .  User input -----------------------------------------------------------
csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"
angle_name = "angles" #The title of the column with angles; NO SPACES PLEASE
angle_unit = "degrees" # "degrees" or "radians"

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


# Read in file ------------------------------------------------------------
adata = read.table(file = path_file,#read from user-selected file
                   header = T,#read the file header to use for variable names
                   sep = csv_sep#,#values are separated by the user-specified character
                   #other parameters can be added here for troubleshooting
)

View(adata)#show the user the data that was
plot.circular(x = circular(x = adata$angle, 
                           type = 'angles',
                           unit = 'degrees',
                           template = 'geographics',
                           modulo = '2pi',
                           zero = 'pi',
                           rotation = 'clock'
),
stack = TRUE,
bins = 360/5,
sep = 0.05
)

# Perform Rayleigh test ---------------------------------------------------
ray_test = with(adata,
                rayleigh.test( x = circular(angle, unit = angle_unit)
                )
)
print(ray_test)

# Save result -------------------------------------------------------------
co_ray_test = capture.output(ray_test)
#save result
write.table(x = with(ray_test,#save relevant info in a table
                     cbind(data = path_file,
                           test = co_ray_test[2],
                           h0 = 'the distribution of angles is uniform',
                           h1 = 'angles are from a unimodal von Mises distribution',
                           R_statistic = statistic,
                           # mean_angle = mu,
                           p = p.value
                     )
),
file = file.path(dirname(path_file),#same place as original with
                 paste(sub(pattern='.csv',#similar name
                           replacement = '',
                           x = basename(path_file)
                 ),
                 'Rayleigh test result.csv'#ending in mean test
                 )
),
row.names = FALSE,#rows do not need names
sep = csv_sep #Use same separator as original csv file

)

