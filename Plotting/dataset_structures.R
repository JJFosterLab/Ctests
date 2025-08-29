#FOR A 'CLEAN' RUN, RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2025 08 28
#     MODIFIED:	James Foster              DATE: 2025 08 29
#
#  DESCRIPTION: Plot example structures of circular dataset
#               
#       INPUTS: 
#               
#      OUTPUTS: Figures
#
#	   CHANGES: - 
#
#   REFERENCES: Batschelet E (1981).
#               The Rayleigh test, Chap 4.2, p. 54
#               Chapter 4: Tests for Randomness and Goodness-of-fit
#               In: Circular Statistics in Biology
#               Academic Press (London)
#
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Example homeward direction plot  +
#- Example concentration change plot  +
#- Example individual headings plot correlated  +
#- Example individual headings plot uncorrelated  +
#- Example individual headings and concentration  +
#- Example of individual change in heading

set.seed(0120810506)#ISBN Batschelet, 1981


# Set up plot colours
col_kappa = "#1E78B5"
col_rho = '#F08024'
col_sd = '#E74A29'
col_sd2 = '#E57461'
col_pdf = adjustcolor(col = '#21A885',
                      alpha.f = 0.7)
col_obs = '#3E1F51'

# Set up functions --------------------------------------------------------


require(circular)
#function for estimating sd from kappa
MardiaSD = function(k,
                    method = 'analytical',
                    n = 1e4)
{
  switch(method,
         simulation = 
           circular::sd.circular(
             circular::rvonmises(n = n, 
                                 mu = circular::circular(pi, # this is directly in the middle of (0,2pi)
                                                         template = 'none'),
                                 kappa = k)
           ),
         analytical = sqrt( -2* log( circular::A1(k))),
         sqrt( -2* log( circular::A1(k)))
  )
}

#estimate normal SD from kappa
NormalSD = function(k,
                    n = 1e4)
{
  sd(
    as.numeric(
      circular::rvonmises(n = n, 
                          mu = circular::circular(pi, # this is directly in the middle of (0,2pi)
                                                  template = 'none'),
                          kappa = k)
    )
  )
}

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
#Simulate confidence intervals for a unimodal or bimodal distribution
#fitted to a vector of "angles"
CI_vM = function(angles, #vector of angles fitted (used for sample size)
                 m1, #primary mean
                 k1, #primary concentration
                 m2 = NA, #secondary mean (ignored if NULL or NA)
                 k2 = NA, #secondary kappa
                 w1 = 1, #weighting of primary mean
                 force_mu = FALSE, #force median at true mu?
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
        if(force_mu)
        {
          Mod360.180(
            quantile( Mod360.180(as.numeric(m1_est) - m1),
                      probs = probs1) + m1
          )
        }else
        {
          Mod360.180(
            quantile.circular(x = circular(x = m1_est,
                                           units = au,
                                           rotation = ar),
                              probs = probs1)
          )
        }
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
        list(m1 = if(force_mu)
        {
          Mod360.180(
            quantile( Mod360.180(as.numeric(m1_est) - m1),
                      probs = probs1) + m1
          )
        }else
        {
          Mod360.180(
            quantile.circular(x = circular(x = m1_est,
                                           units = au,
                                           rotation = ar),
                              probs = probs1)
          )
        },
        m2 = if(force_mu)
        {
          Mod360.180(
            quantile( Mod360.180(as.numeric(m2_est) - m2),
                      probs = probs2) + m2
          )
        }else
        {
          Mod360.180(
            quantile.circular(x = circular(x = m2_est,
                                           units = au,
                                           rotation = ar),
                              probs = probs2)
          )
        }
        )
      }
    }else
    { #if quantiles not requested, return the simulations (mainly for troubleshooting)
      if(is.na(m2))
      {
        m1_est = 
          sapply(X = m1_est, FUN = Mod360.180)
      }else
      {
        list(
          m1_est = 
            sapply(X = m1_est, FUN = Mod360.180),     
          m2_est = 
            sapply(X = m2_est, FUN = Mod360.180),
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
  ci_vec = as.numeric(ci_vec)#remove circular formatting!
  #changed on 20250815, plotting issues near median
  angle_seq1.1 = 
    seq(from = ci_vec[1], #lower
        to = ci_vec[1] +
          Mod360.180(ci_vec[2]-ci_vec[1]), #median
        length.out =1e2/2)
  angle_seq1.2 = 
    seq(from = ci_vec[2], #median
        to = ci_vec[2] +
          Mod360.180(ci_vec[3]-ci_vec[2]) , #upper
        length.out =1e2/2)
  lines(x = radius*sin( rad(angle_seq1.1) ),
        y = radius*cos( rad(angle_seq1.1) ),
        col = col,
        lwd = lwd,
        lend = 'butt',
        ...
  )
  lines(x = radius*sin( rad(angle_seq1.2) ),
        y = radius*cos( rad(angle_seq1.2) ),
        col = col,
        lwd = lwd,
        lend = 'butt',
        ...
  )
  if(!is.na(ci_vec[4]))
  {
    #changed on 20250815
    angle_seq2.1 = 
      seq(from = ci_vec[1+3],
          to = ci_vec[1+3] +
            Mod360.180(ci_vec[2+3]-ci_vec[1+3]),
          length.out =1e2/2)
    
    angle_seq2.2 = 
      seq(from = ci_vec[2+3],
          to = ci_vec[2+3] +
            Mod360.180(ci_vec[3+3]-ci_vec[2+3]) ,
          length.out =1e2/2)
    lines(x = radius*sin( rad(angle_seq2.1) ),
          y = radius*cos( rad(angle_seq2.1) ),
          col = col,
          lwd = lwd,
          lend = 'butt',
          ....)
    lines(x = radius*sin( rad(angle_seq2.1) ),
          y = radius*cos( rad(angle_seq2.1) ),
          col = col,
          lwd = lwd,
          lend = 'butt',
          ...)
  }
}


DescriptCplot = function(k,
                         m = 0,#in degrees
                         ndata = 20,
                         lw = 3,
                         pcol = col_obs,
                         mvcol = col_rho,
                         sdcol = col_sd,
                         cicol = col_rho,
                         denscol = col_pdf,
                         bins = 360/5-1,
                         stack = TRUE,
                         sep = 0.05,
                         shrink = 1.25,
                         refline = NA,
                         save_sample = FALSE,
                         ...#passed to points.circular
)
{
  #convert mu to circular class
  cm = circular(x = m,
                units = 'degrees',
                rotation = 'clock',
                zero = pi/2)
  #kappa = 0 v. random, make repeatable
  set.seed(20250815)
  #generate dataset
  cd = rvonmises(mu = cm,
                 kappa = k,
                 n = ndata)
  #Generate a vector of angles
  ra = seq(from = -pi, to = pi, by = 0.01)
  #plot dataset
  plot.circular(x = circular(NA,
                             units = 'degrees',
                             rotation = 'clock',
                             zero = pi/2),
                stack = stack,
                bins = bins,
                shrink = shrink,
                sep = sep,
                col = pcol,
                axes = FALSE
  )
  #add reference line
  if(!is.na(refline))
  {
    arrows.circular(x = circular(refline,
                                 units = 'degrees',
                                 rotation = 'clock',
                                 zero = pi/2),
                    col = 'gray',
                    length = 0,
                    lend = 'butt',
                    lwd = lw)
  }
  #add probability density
  lines.circular(x = ra,
                 y = dvonmises(x = circular(ra,
                                            units = 'radians',
                                            rotation = 'clock',
                                            zero = pi),
                               mu = cm,
                               kappa = k)/shrink,
                 col = denscol,
                 lwd = lw)
  points.circular(cd,
                  stack = stack,
                  bins = bins,
                  shrink = shrink,
                  sep = sep,
                  col = pcol,
                  ...
  )
  #add the mean vector
  if(k>0)
  {
    arrows.circular(x = cm,
                    y = A1(k),
                    col = mvcol,
                    lwd = lw,
                    length = 0.1/shrink#,shrink = 1/shrink
    )
  }else
  {
    points(x = 0,
           y = 0, 
           col = mvcol,
           pch = 19,
           lwd = lw)
  }
  #add SD as lines
  #estimate Mardia SD
  msd = MardiaSD(k = k)
  arrows.circular(x = cm+circular(deg(msd),
                                  units = 'degrees',
                                  rotation = 'clock',
                                  zero = pi/2),
                  y = 1+5*sep*shrink,
                  col = sdcol,
                  lwd = lw,
                  lty = 3,
                  length = 0)
  arrows.circular(x = cm-circular(deg(msd),
                                  units = 'degrees',
                                  rotation = 'clock',
                                  zero = pi/2),
                  y = 1+5*sep*shrink,
                  col = sdcol,
                  lwd = lw,
                  lty = 3,
                  length = 0)
  # lines.circular(x = circular(rep(0-deg(msd), 2),
  #                             units = 'degrees',
  #                             rotation = 'clock',
  #                             zero = pi/2),
  #                 y = c(2,3)*sep*shrink,
  #                col = sdcol,
  #                lwd = lw)
  #add CI of the mean
  #calculate vector of estimates
  civ = CI_vM(angles = cd,
              m1 = cm,
              k1 = k,
              alternative = 'two.sided',
              force_mu = if(k == 0 ){TRUE}else{FALSE})
  #plot quantiles of estimates
  PlotCI_vM(ci_vec = civ,
            col = cicol, lwd = lw,
            radius = 1+5*sep*shrink)
  #add a title
  mtext(text = paste0('κ = ', k),
        side = 1,
        line = -3)
  if(save_sample)
  {return(cd)}
}

PCfun = function(angles,
                 col,
                 shrink = 1.5,
                 title = '',
                 side = 1)
{
  ca = circular(x = angles,
                units = 'degrees',
                rotation = 'clock')
  plot.circular(x = ca,
                col = col,
                stack = TRUE,
                sep = 0.1,
                bins = 355/5,
                units = 'degrees',
                rotation = 'clock',
                zero = pi/2,
                shrink = shrink)
  mtext(text = title,
        side = side,
        line = -2)
  lines(x = c(0,0),
        y = c(-1,1),
        col = 'gray')
  arrows.circular(x = mean.circular(ca),
                  y = rho.circular(ca),
                  zero = pi/2,
                  rotation = 'clock',
                  col = col,
                  length =0.1)
}


#invert the softplus link
#https://en.wikipedia.org/wiki/Softplus
#we are using this as our _inverse_ link function for kappa,
#maps almost 1:1 but keeps values >0 for low estimates
softplus = function(x)
{
  log(exp(x)+1) 
}
#this would return our kappa estimates back to the original scale
inv_softplus = function(x)
{
  log(exp(x)-1) 
}


# Simulate divergence from home direction -------------------------------------------
set.seed(0120810506)#ISBN Batschelet, 1981

par(pty = 's')
par(mar = c(0,0,0,0))

DescriptCplot(m = -15,
              k = 10,
              refline = 0, 
              ndata = 20,
              sdcol = NA,
              denscol = NA)

# Simulate reduced concentration -------------------------------------------
set.seed(0120810506)#ISBN Batschelet, 1981

par(pty = 's')
par(mar = c(0,0,0,0))
par(mfrow = c(1,2))

DescriptCplot(m = 0,
              k = 3,
              ndata = 20,
              sdcol = NA,
              denscol = NA)
lines(x = c(0,0),
      y = c(0,1),
      col = adjustcolor(col = 'gray', 
                        alpha.f = 150/255),
      lwd = 3,
      lend = 'butt')
DescriptCplot(m = 0,
              k = 0.5,
              ndata = 20,
              sdcol = NA,
              denscol = NA)
lines(x = c(0,0),
      y = c(0,1),
      col = adjustcolor(col = 'gray', 
                        alpha.f = 150/255),
      lwd = 3,
      lend = 'butt')


# Simulate datasets with low interindiv correlation ---------------------------------------
set.seed(0120810506)#ISBN Batschelet, 1981
kappa_mu = 0.1
kappa_id = 5.0

ndata = 10 # moderate sample size
# list of circular datasets
c0 = circular(x = 0,
              units = 'degrees',
              rotation = 'clock',
              zero = pi/2)
dt2 = rvonmises(n = 10,
                mu = c0,
                kappa = kappa_mu)

par(pty = 's')
par(mar = c(0,0,0,0))
par(mfrow = c(3,4))
dt_id2 = lapply(X = dt2, #rev(kd), # should this be largest to smallest?
       FUN = DescriptCplot,
       save_sample = TRUE,
       k = kappa_id,
       ndata = 20,
       refline = 0,
       sdcol = NA,
       denscol = NA)
#Add the population of biases
DescriptCplot(k = kappa_mu,
              ndata = 10,
              refline = 0,
              sdcol = NA,
              denscol = NA,
              pcol = NA,
              cicol = col_sd,
              mvcol = col_sd
              )
points.circular(dt2,
                bins = 360/5-1,
                stack = TRUE,
                sep = 0.05,
                shrink = 1.25,
                col = col_rho
)

dt_comb2 = do.call(what = c,
                  args = dt_id2)
PCfun(angles = dt_comb2,
      col = 'gray25',
      shrink = 3.0)
mle_comb2 = mle.vonmises(x = dt_comb2,bias = TRUE)
ci_comb2 = with(mle_comb2,
               CI_vM(angles = dt_comb2,
                     m1 = mu,
                     k1 = kappa,
                     alternative = 'two.sided')
              )
with(mle_comb2,
     {
arrows.circular(x = circular(mu,
                             units = 'degrees',
                             rotation = 'clock',
                             zero = pi/2),
                y = A1(kappa),
                lwd = 3,
                col = col_pdf,
                length = 0.1
                )
     }
)
PlotCI_vM(ci_vec = ci_comb2,
          col = col_pdf)
rayleigh.test(dt_comb2)
rayleigh.test(dt_comb2[1+0:9 * 20])

# Simulate dataset with €variable individual parametes ---------------------------------------
set.seed(0120810506)#ISBN Batschelet, 1981
kappa_mu_var = 1.0
kappa_var_mean = 1.5
kappa_var_sd = 2.0
kappa_id_var = rnorm(n = ndata,
                     mean = kappa_var_mean,
                     sd = kappa_var_sd)
#rectified
kappa_id_var[kappa_id_var<0] = 0

# list of circular datasets
dt_var = rvonmises(n = ndata,
                mu = c0,
                kappa = kappa_mu_var)

par(pty = 's')
par(mar = c(0,0,0,0))
par(mfrow = c(3,4))
dt_id_var = mapply(m = dt_var, 
                   k = round(kappa_id_var,2), 
       FUN = DescriptCplot,
       save_sample = TRUE,
       ndata = 20,
       refline = 0,
       sdcol = NA,
       denscol = NA,
       SIMPLIFY = FALSE)
#Add the population of biases
DescriptCplot(k = kappa_mu_var,
              ndata = ndata,
              refline = 0,
              sdcol = NA,
              denscol = NA,
              pcol = NA,
              cicol = col_sd,
              mvcol = col_sd
              )
points.circular(dt_var,
                bins = 360/5-1,
                stack = TRUE,
                sep = 0.05,
                shrink = 1.25,
                col = col_rho
)
mean_kappa_id_var = softplus(mean(inv_softplus(kappa_id_var)))
#Add decription of the average individual
DescriptCplot(k = kappa_var_mean,
              ndata = ndata,
              refline = 0,
              sdcol = NA,
              denscol = NA,
              pcol = NA,
              cicol = col_rho,
              mvcol = col_rho
              )
kappa_id_var_ci = kappa_var_mean + 
                    kappa_var_sd * 
                      qnorm(c(0,1) + c(1,-1)*0.05/2)
#rectify
kappa_id_var_ci[kappa_id_var_ci<0] = 0


arrows(x0 = sin(c0),
       x1 = sin(c0),
       y0 = A1(kappa_id_var_ci[1]),
       y1 = A1(kappa_id_var_ci[2]),
       lwd = 7,
       col = adjustcolor(col = col_sd2,
                        alpha.f = 100/255),
       length = 0.05,
       angle = 90,
       code = 3,
       lend = 'butt'
       )
mtext(text = paste0('(',paste(signif(kappa_id_var_ci, 2), collapse = ' '), ')'),
      side = 1,
      line = -1)
#v-test is not significant (or Rayleigh test)
rayleigh.test(dt_var, mu = c0)
rt_lst = lapply(X = lapply(dt_id_var, circular, units = 'degrees'),
                 FUN = rayleigh.test)
rt_lst_print = do.call(what = rbind, 
                       args = rt_lst)
print(rt_lst_print)
