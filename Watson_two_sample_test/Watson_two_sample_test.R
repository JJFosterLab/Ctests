#FOR A 'CLEAN' RUN, RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2021 07 19
#     MODIFIED:	James Foster              DATE: 2021 07 20
#
#  DESCRIPTION: Loads a text file and performs Watson's two-sample test
#               for equality of distributions.
#               
#       INPUTS: A ".csv" table with columns for experiment group ("stimulus") and
#               angle ("angle").
#               User should specify test details (line 50).
#               
#      OUTPUTS: Results table (.csv).
#
#	   CHANGES: - Suppressed pacakge loading messages (upsetting users)
#             - 
#             - 
#
#   REFERENCES: Jammalamadaka, S. Rao and SenGupta, A. (2001).
#               Topics in Circular Statistics, Section 7.5, 
#               World Scientific Press (Singapore)
#
#    EXAMPLES:  Fill out user input (lines 50-56), then press ctrl+shift+s to run
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
suppressMessages(#these are disturbing users unnecessarily
  {
  require(circular)#package for handling cirular data
  require(CircStats)#package for circular hypothesis tests
  }
)




# Input Variables ----------------------------------------------------------

#  .  User input -----------------------------------------------------------
csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"
group_factor = "stimulus" #The title of the column; NO SPACES PLEASE
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
    # n_phases = 2#must be exactly two groups
    # n_angles = 10
    # angles_sim = sapply(sort(rep(0:(n_phases - 1)*pi/2, n_angles)),
    #                     function(mm)
    #                     {suppressWarnings(rvonmises(n = 1,mu = mm, kappa = 3))}
    #                     )
    # sim = data.frame(stimulus = sort(rep(1:n_phases, n_angles)),
    #                  angle = round(angles_sim*180/pi)
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
# #TODO, make this more flexible
# mdata$ID<- mdata[,individual_factor]#column contains ID in some form
# mdata$correct_choices<- mdata[,correct_name]#column contains total positive per animal
# mdata$incorrect_choices<- mdata[,incorrect_name]#column contain

# Perform Watson's two-sample test -----------------------------------------
w2_test = with(adata,
                watson.two.test( x = circular(subset(angle, stimulus %in% unique(stimulus)[1]), 
                                              units = angle_unit
                                              ),
                                 y = circular(subset(angle, stimulus %in% unique(stimulus)[2]), 
                                              units = angle_unit
                                 ),
                )
)
print(w2_test)

# Save result -------------------------------------------------------------
co_w2_test <- capture.output(w2_test)
#save result
write.table(x = with(w2_test,#save relevant info in a table
                     cbind(data = path_file,
                           test = co_w2_test[2],
                           h0 = 'all angles come from the same distribution',
                           h1 = 'difference between groups',
                           U2_statistic = statistic,
                           n1 = nx,
                           n2 = ny,
                           p = co_w2_test[5]
                     )
),
file = file.path(dirname(path_file),#same place as original with
                 paste(sub(pattern='.csv',#similar name
                           replacement = '',
                           x = basename(path_file)
                 ),
                 'Watson two-sample test result.csv'#ending in mean test
                 )
),
row.names = FALSE,#rows do not need names
sep = csv_sep #Use same separator as original csv file

)

