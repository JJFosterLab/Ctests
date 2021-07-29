#FOR A 'CLEAN' RUN, RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2021 07 29
#     MODIFIED:	James Foster              DATE: 2021 07 29
#
#  DESCRIPTION: Loads a text file and performs a Prentice generalised
#               rank-sum test for equal rankings within blocks.
#               This can be used, for example, to compare groups of
#               mean vector lengths calculated for the same animal under
#               different conditions.
#               
#       INPUTS: A ".csv" table with a column of animal IDs, or another group ID,
#               ("ID"), a column of conditions ("condition") and a column of
#               mean vector lengths ("rho"), or other non-circular variable,
#               or angles ("angle") if a mean vector needs to be calculated.
#               User should specify test details (line 70).
#               
#      OUTPUTS: Results table (.csv).
#
#	   CHANGES: - 
#             - 
#
#   REFERENCES: Prentice M. J. (1979).
#               On the Problem of m Incomplete Rankings
#               Biometrika, 66(1), 167-170
#               
#               stats.stackexchange.com/questions/129348/is-there-a-non-parametric-repeated-measures-test-for-replicated-block-data
#               
#               Used in:
#               Foster JJ, Kirwan JD, el Jundi B, Smolka J, Khaldy L,
#               Baird E, Byrne MJ, Nilsson D-E, Johnsen S, Dacke M (2019). 
#               Orienting to polarized light at night –
#               matching lunar skylight to performance in a nocturnal beetle.
#               Journal of Experimental Biology 222, jeb188532
#               
#               Dacke, M., Bell, A.T.A., Foster, J.J., Baird, E.J.,
#               Strube-Bloss, M.F., Byrne, M.J., el Jundi, B., (2019). 
#               Multimodal cue integration in the dung beetle compass. 
#               PNAS 201904308
#               
#               Foster JJ, Tocco C, Smolka J, Khaldy L, Baird E, Byrne M, 
#               Nilsson D-E & Dacke M (2021)
#               Light pollution forces a change in dung beetle orientation behavior,
#               Current Biology
#               https://doi.org/10.1016/j.cub.2021.06.038
#
#    EXAMPLES:  Fill out user input (lines 70-75), then press ctrl+shift+s to run
# 
#TODO   ---------------------------------------------
#TODO   
#- Read in data   +
#- Perform test   +
#- Test with simulated data +
#- Save results +
#- Version with angles  +
#- Individual circular plots  

# Useful functions --------------------------------------------------------
# . Load package ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
  require(circular)#package for handling circular data
  require(muStat)#package for Prentice test
  }
)
# Input Variables ----------------------------------------------------------
#  .  User input -----------------------------------------------------------
csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"
calculate_rho = FALSE #TRUE or FALSE: are data angles from which mean vector should be calculated?
alternative = 'two.sided'# 'two.sided', 'greater' or 'less': expected difference of condition 1 versus others 
if(calculate_rho)
{angle_name = "angles" #The title of the column with angles; NO SPACES PLEASE
angle_unit = "degrees"} # "degrees" or "radians"

#Check the operating system and assign a logical flag (TRUE or FALSE)
sys_win <- Sys.info()[['sysname']] == 'Windows'
#On computers set up by JMU Würzburg, use user profile instead of home directory
if(sys_win){
  #get rid of all the backslashes
  ltp <- gsub('\\\\', '/', Sys.getenv('USERPROFILE'))#Why does windows have to make this so difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp <- Sys.getenv('HOME')#Life was easier on Mac
}

# Simulate data (not used) ------------------------------------------------
    # n_animals = 20
    # n_angles = 10
    # n_conditions = 3
    # angles_sim = sapply(1:(n_conditions*n_animals),
    #                     function(an)
    #                     {suppressWarnings(
    #                       rvonmises(n = n_angles,
    #                                 mu = rcircularuniform(1),
    #                                 kappa = 1+ an%%n_conditions
    #                                 # kappa = ifelse(test = an < n_animals,
    #                                 #                yes = 2.0,
    #                                 #                no = 1.0
    #                                 #                )
    #                                 )
    #                       )}
    #                     )
    # sim = data.frame(ID = sort(rep(1:n_animals, n_angles*n_conditions)),
    #                  condition = rep(sort(rep(1:n_conditions,n_angles)), n_animals),
    #                  angle = round(c(angles_sim)*180/pi)
    #                  )
    # if(!calculate_rho)
    # {
    #   sim = aggregate(angle~ID*condition,
    #                   data = sim,
    #                   FUN = function(aa)
    #                     {rho = rho.circular(
    #                       circular(aa, unit = 'degrees')
    #                     )}
    #                   )
    #   sim = within(sim,{rho = angle; rm(angle)})
    # }
    # write.table(x = sim,
    #             file = file.path(ltp,'Documents', "simulated_rho.csv"),
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

View(adata)#show the user the data that was read in
if(calculate_rho)
{
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
                sep = 0.05,
                shrink = sqrt(length(adata$angle))/15
                )
  adata0 = adata#save original data
  adata = aggregate(angle~ID*condition,#for each combination of ID & condition
                  data = sim,
                  FUN = function(aa)
                  {rho = rho.circular(#calculate mean vector
                    circular(aa, unit = angle_unit)#using the user-specified format
                  )}
  )
  adata = within(adata,{rho = angle; rm(angle)})
  
}
stripchart(round(rho,2)~condition,
           data = adata,
           ylab = 'mean vector length',
           xlab = 'condition',
           method = 'stack',
           vertical = TRUE,
           col = 2:(length(unique(adata$condition))+1),
           pch = 19,
           xlim = range(adata$condition)+c(-0.5,0.5),
           ylim = c(0,1)
)
abline(h = c(0,1))
invisible(
lapply(unique(adata$ID),
       function(id)
         {
          with(subset(adata, ID %in% id), 
                lines(x = match(condition, unique(adata$condition)),
                      y = rho,
                      col = rgb(0,0,0,0.3),
                      lwd = 2
                      )
          )
         }
       )
)

# Perform Prentice test --------------------------------------------------
prent_test = with(adata,
                  #check if mean vectors differed in rank across conditions
                  prentice.test(y = rho,
                                groups = condition,
                                blocks = ID,
                                alternative = alternative
                               )
)
print(prent_test)


# Calculate effect size ---------------------------------------------------
#Reviewer 2 for Dacke et al., 2019 was concerned about the effect size
# For a chi-squared statistic, we can calculate Cramér's V
# this is a correlation-coefficient-like value between 0 and 1
crV = with(adata,
           {
             round( 
               sqrt((prent_test$statistic/length(rho))/
                     min(	length(unique(condition))-1,	
                          length(unique(ID))-1	
                          )
                    ),
               3 # round to 3 decimal places
               )
           }
          )
names(crV) <- "Cramér's V"
print(crV)
## See rcompanion version for reference
    # Instalload('rcompanion')
    # cramerV(wr$condition, wr$rho)
    # V={\sqrt  {{\frac  {\chi ^{2}/n}{\min(k-1,r-1)}}}} 
    # \chi ^{2} is derived from Pearson's chi-squared test
    # cramerV#look at how this is coded in rcompanion
    # N = length(x)
    # Chi.sq = suppressWarnings(chisq.test(x, y, correct = FALSE, 
    # ...)$statistic)
    # Phi = Chi.sq/N
    # Row = length(unique(x))
    # C = length(unique(y))
    # CV = sqrt(Phi/min(Row - 1, C - 1))
    #looks the same, but how are rows and columns handled in the Prentice test?
    #looking at Prentice (1979) this seems correct, resulting Chisq has 1 d.f., which is nGroups - 1.

# Save result -------------------------------------------------------------
co_prent_test = capture.output(prent_test)
#save result
write.table(x = with(prent_test,#save relevant info in a table
                     cbind(data = path_file,
                      test = co_prent_test[2],
                      h0 = 'mean vector lengths are the same',
                      h1 = switch(alternative,
                                  two.sided = 'mean vector lengths differ by condition',
                                  less = 'mean vectors are shorter for condition 1',
                                  greater = 'mean vectors are shorter for condition 1'
                                  ),
                      condition_1 = unique(adata$condition)[1],
                      chi2_statistic = statistic,
                      df = parameters,
                      p = p.value,
                      cramer_v = crV
                           )
                    ),
file = file.path(dirname(path_file),#same place as original with
                 paste(sub(pattern='.csv',#similar name
                           replacement = '',
                           x = basename(path_file)
                 ),
                 'Prentice test result.csv'#ending in mean test
                 )
),
row.names = FALSE,#rows do not need names
sep = csv_sep #Use same separator as original csv file

)

