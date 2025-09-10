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


# Simulate continuous sequences -------------------------------------------

#Generate a vector of angles
ra = seq(from = -pi, to = pi, by = 0.01)
#Generate a vector of plausible kappas
kk = seq(from = 0, to  = 100, by = 0.01)
#Convert to an estimate of mean vector length
rr = A1(kk)
#Convert to estimate of circular SD
# 5x faster with parallel
sm = parallel::parSapply(X = kk, FUN = MardiaSD, n = 1e4,
                         cl = parallel::makePSOCKcluster(parallel::detectCores() - 1)
                         )
#Convert to estimate of normal SD
sn = parallel::parSapply(
                         X = kk, FUN = NormalSD, n = 1e4,
                         cl = parallel::makePSOCKcluster(parallel::detectCores() - 1)
                        )

#estimate difference in degrees
sn_sm_absdiff = deg(abs(sn-sm))
#largest normal sd >1° from correct
sn_best_idx = max(which(sn_sm_absdiff>1))
sn_ok_idx = max(which(sn_sm_absdiff>5))

#for low k SD is unstable, generate more
kk_small = seq(from = 0.01, to = 0.05, length.out = 1e3)
sm_small = parallel::parSapply(X = kk_small, FUN = MardiaSD, method = 'simulation', n = 1e4,
                         cl = parallel::makePSOCKcluster(parallel::detectCores() - 1)
)
sn_small = parallel::parSapply(X = kk_small, FUN = NormalSD, n = 1e4,
                         cl = parallel::makePSOCKcluster(parallel::detectCores() - 1)
)

# Calculate discrete benchmarks -------------------------------------------

#Discrete sequence of benchmarks 
kd = c(0, 0.1,  c(0.3, 0.5, 1.0) %*% t(10^c(0:2)) )
#mean vector lengths
rd = round(A1(kd), 3)
#Mardia SD
ds = round(deg(sapply(X = kd, MardiaSD, n = 1e4)))


# Simulate datasets from parameters ---------------------------------------
ndata = 20 # moderate sample size
# list of circular datasets
cdata = lapply(X = kd, #rev(kd), # should this be largest to smallest?
               FUN = rvonmises,
               mu = circular(x = 0,
                             units = 'degrees',
                             rotation = 'clock',
                             zero = pi/2),
               n = ndata)

# Plots of value sequence -------------------------------------------------
par(mar = c(1,1,0,1)*4.0)
#set up plot
plot(x = NULL,
     ylim = c(0,1),
     xlim = c(0.01,100),
     type = 'l',
     lwd = 5,
     xlab = 'κ',
     ylab = 'mean vector length',
     log = 'x',
     axes = F)
axis(side = 2, col = col_rho, lwd = 2)
axis(side = 1,
     at = c(0, c(1,5) %*% t(10^c(-2:2)) ),
     labels = sapply(X = c(0, c(1,5) %*% t(10^c(-2:2)) ),
                 FUN = paste0),
     lwd = 2
    )
axis(side = 4,
     at = 0:12/12,
     labels = paste0(0:12*15, '°'),
     lwd = 2,
     col = col_sd
    )
mtext(text = 'circular SD',
      side = 4,
      line = 3)
abline(h = c(0,1),
       v = 0,
       lwd = 2)
#add data
#rho
lines(x = kk,
      y = rr,
      col = col_rho,
      lwd = 7)
#sd needs some smoothing for small values
  # #some smoothing for large kappa
  # lines(x = kk[kk>max(kk_small)],
  #       y = sm[kk>max(kk_small)]/pi,
  #       col = col_sd,
  #       lwd = 7)
#some smoothing for large kappa
lines(x = kk[sm<pi],
      y = sm[sm<pi]/pi,
      col = col_sd,
      lwd = 7)
# loess for small kappa
lines(x = kk_small,
      y = predict(
                  loess(formula = sm_small/pi ~ kk_small,
                        span = 1.0)
                  ),
      col = col_sd2,
      lwd = 5,
      lty = 3,
      lend = 'butt')


# Comparison of Mardia SD and Gaussian SD -----
par(mar = c(1,1,0,0)*4.0)
plot(x = NULL,
     xlab = 'normal SD',
     ylab = 'circular SD',
     xlim = c(0,180),
     ylim = c(0,180),
     col = col_sd,
     type = 'l',
     lwd = 7,
     axes = FALSE)
axis(side = 1, 
     at = 0:12*15,
     labels = paste0(0:12*15, '°'),
     col = 'black', 
     lwd = 2)
axis(side = 2, 
     at = 0:12*15,
     labels = paste0(0:12*15, '°'),
     col = col_sd, 
     lwd = 2)
abline(h = c(0,180),
       v = c(0,180),
       lwd = 2)

# lines(x = deg(sn),
#       y = deg(sm),
#       col = col_sd,
#       lwd = 7)

lines(x = deg(sn[sm<pi]),
      y = deg(
        predict(
          loess(formula = sm[sm<pi]~sn[sm<pi],
                span = 0.01)
        )
      ),
      col = col_sd,
      lwd = 7)

# loess for small kappa
lines(x = deg(sn_small[sm_small<pi]),
      y = deg(
        predict(
        loess(formula = sm_small[sm_small<pi] ~ sn_small[sm_small<pi],
              span = 1.0)
        )
      ),
      col = col_sd2,
      lwd = 5,
      lty = 3,
      lend = 'butt')

abline(a = 0,
       b = 1,
       col = 'gray',
       lty = 1,
       lwd = 3)
abline(a = c(-1)*5,
       b = c(1),
       col = 'gray',
       lty = c(2),
       lwd = 3)
abline(a = c(1)*5,
       b = c(1),
       col = 'gray',
       lty = c(2),
       lwd = 3)


segments(x0 = deg( sn[c(sn_best_idx, sn_ok_idx)] ),
         y0 = rep(0, length(c(sn_best_idx, sn_ok_idx))),
         y1 = deg( sm[c(sn_best_idx, sn_ok_idx)] ),
         col = col_sd2,
         lty = 3)
segments(x0 = rep(0, length(length(c(sn_best_idx, sn_ok_idx)))),
         x1 = deg( sn[c(sn_best_idx, sn_ok_idx)] ),
         y0 = deg( sm[c(sn_best_idx, sn_ok_idx)] ),
         col = col_sd2,
         lty = 3)
text(x = deg( sn[c(sn_best_idx, sn_ok_idx)] )+2,
     y = deg( sm[c(sn_best_idx, sn_ok_idx)] )-10,
     labels = paste0( round(deg(sn[c(sn_best_idx, sn_ok_idx)]) ), '°\n',
                      'bias ', c(1,5), '°'),
     col = col_sd2,
     adj = c(0,1),
     cex = 0.7
)

#out of curiousity, how do rho and sd scale?
plot(x = rr,
     y = deg(sm),
     pch = 19,
     col = gray(0, 0.3),
     axes = FALSE,
     xlab = 'mean vector length',
     ylab = 'circular SD',
     type = 'l',
     lwd = 7
     )
abline(v = c(0,1), h = c(0,180))
axis(side = 2, 
     at = 0:12*15,
     labels = paste0(0:12*15, '°'),
     col = col_sd, 
     lwd = 2)
axis(1, col = col_rho)


# Plot benchmarks -----
#for rho
par(mar = c(1,1,0,0)*4.0)
#set up plot
plot(x = NULL,
     ylim = c(0,1),
     xlim = range(kd[kd>0]) + c(0,30),
     type = 'l',
     lwd = 5,
     xlab = 'κ',
     ylab = 'mean vector length',
     log = 'x',
     axes = F)
axis(side = 2, col = col_rho, lwd = 2)
axis(side = 1,
     at = kd,
     labels = sapply(X = kd,
                     FUN = paste0),
     lwd = 2
)
abline(h = c(0,1),
       v = 0,
       lwd = 2)
#rho
lines(x = kk,
      y = rr,
      col = col_rho,
      lwd = 7)
segments(x0 = kd,
         y0 = rep(0, length(kd)),
         y1 = rd,
         col = col_rho,
         lty = 3)
segments(x0 = rep(0.007, length(kd)),
         x1 = kd,
         y0 = rd,
         col = col_rho,
         lty = 3)
text(x = exp(log(kd)+0.01),
     y = rd-0.01,
     labels = paste0(round(rd, 3)),
     col = col_rho,
     adj = c(0,1),
     cex = 0.7
)

#for sd
par(mar = c(1,0,0,1)*4.0)
#set up plot
plot(x = NULL,
     ylim = c(0,180),
     xlim = range(kd[kd>0]) + c(0,10),
     type = 'l',
     lwd = 5,
     xlab = 'κ',
     ylab = '',
     log = 'x',
     axes = F)
axis(side = 4, 
     at = 0:12*15,
     labels = paste0(0:12*15, '°'),
     col = col_sd, 
     lwd = 2)
mtext(text = 'circular SD',
      side = 4,
      line = 3)
axis(side = 1,
     at = kd,
     labels = sapply(X = kd,
                     FUN = paste0),
     lwd = 2
)
abline(h = c(0,180),
       v = 0,
       lwd = 2)

#sd needs some smoothing for small values
#no smoothing for large kappa
lines(x = kk[kk>max(kk_small)],
      y = deg( sm[kk>max(kk_small)] ),
      col = col_sd,
      lwd = 7)

# loess for small kappa
lines(x = kk_small,
      y = predict(
        loess(formula = sm_small/pi ~ kk_small,
              span = 1.0)
      ),
      col = col_sd2,
      lwd = 5,
      lty = 3,
      lend = 'butt')

segments(x0 = kd,
         y0 = rep(0, length(kd)),
         y1 = ds,
         col = col_sd,
         lty = 3)
segments(x0 = rep(130, length(kd)),
         x1 = kd,
         y0 = ds,
         col = col_sd,
         lty = 3)
text(x = kd,
     y = ds+7,
     labels = paste0(ds, '°'),
     col = col_sd,
     adj = c(0,1),
     cex = 0.7
)


# Plot panels -------------------------------------------------------------

# lapply(X = cdata,
#        FUN = plot.circular,
#        stack = TRUE,
#        bins = 360/5 -1,
#        col = 'darkblue',
#        sep = 0.05,
#        shrink = 1.25)
# 
# lines.circular(x = ra,
#                y = dvonmises(x = circular(ra,
#                                           template = 'none'),
#                              mu = circular(0,
#                                            units = 'degrees',
#                                            rotation = 'clock',
#                                            zero = pi/2),
#                              kappa = kd[length(kd)]),
#                col = 'red')

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

}

par(pty = 's')
par(mfrow = c(3,4),
    mar = c(1,0,0,0)*1)


lapply(X = kd,
       FUN = DescriptCplot,
       shrink = 2.4
       )
plot(x =NULL,
     axes = FALSE,
     xlab = '',
     ylab = '',
     main = '',
     xlim = c(0,0),
     ylim = c(0,0))
legend(x = 'left',
       legend = c('probability density',
                  'observed sample (n = 20)',
                  'mean vector',
                  '95% CI of the mean',
                  'mean ± SD'),
       col = c(col_pdf,
               col_obs,
               col_rho,
               col_rho,
               col_sd
               ),
       lty = c(1,
               NA,
               1,
               1,
               3),
       pch = c(NA,
               19,
               NA,
               NA,
               NA
               ),
       lwd = 5,
       cex = 1.0,
       bty = 'n'
       )


# Effect of kappa on PD at 0° ---------------------------------------------

P0by5deg = function(k, interval = 10)
{
  diff(
    pvonmises(q = circular(c(-1, 1)*0.5*interval, #the interval around the mean
                           units = 'degrees'),
              mu = circular(0, 
                            units = 'degrees'),
              kappa = k
    )
    )
}

intrvl = 10
# pd0 = sapply(X = kd,
#              FUN = dvonmises,
#              x = circular(0, template = 'none'),
#              mu = circular(0, template = 'none'))

p0b5 = sapply(X = kd,
             FUN = P0by5deg,
             interval = 10)
pp = sapply(X = kk,
             FUN = P0by5deg,
            interval = 10)


par(mfrow = c(1,1),
    mar = c(4,4,0,0))
plot(x = kk,
     y = pp,
     ylab = paste0('p(0°±',intrvl/2,'°)'),
     xlab = 'κ',
     ylim = c(0,1)*1.0,
     xlim = c(0.01,100),
     col = col_pdf,
     log = 'x',
     type = 'l',
     lwd = 5,
     axes = F)

segments(x0 = kd,
         y0 = rep(0, length(kd)),
         y1 = p0b5,
         col = col_pdf,
         lty = 3)
segments(x0 = rep(0.007, length(kd)),
         x1 = kd,
         y0 = p0b5,
         col = col_pdf,
         lty = 3)
text(x = exp(log(kd)+0.01),
     y = p0b5+0.01,
     labels = paste0(round(p0b5, 3)),
     col = col_pdf,
     adj = c(1,0),
     cex = 0.7
)

axis(side = 2, col = col_pdf, lwd = 2)
axis(side = 1,
     at = c(0, c(1,5) %*% t(10^c(-2:2)) ),
     labels = sapply(X = c(0, c(1,5) %*% t(10^c(-2:2)) ),
                     FUN = paste0),
     lwd = 2
)
abline(h = c(0,1),
       v = 0,
       lwd = 2)
abline(h = intrvl/360,
       lty = 2, 
       col = 'gray',
       lwd = 2)
