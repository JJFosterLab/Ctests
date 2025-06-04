require(circular)
#sample size
nn  =  10
#iterations for estimation
iter  =  1e6

#a function that generates angles from a circular uniform (high measurement precision)
#then uses them to calculate mean vector length, without any bias correction
Rhogen = function(n)
{
  ang = rcircularuniform(n = n)
  rho = rho.circular(x = ang)
  return(rho)
}

#replicate the call for iter iterations
#takes around 10s for 1e5 with "replicate()"
# with parallel::parSapply takes 20s for 1e6
# with parallel::parSapply takes 240s for 1e7
system.time(
  {
    #set up a parallel cluster using all but 1 core
cl = parallel::makePSOCKcluster(parallel::detectCores()-1)
#export all variables used to each cluster
parallel::clusterExport(cl = cl, 
                        varlist = c('Rhogen',
                                    'rcircularuniform',
                                    'rho.circular',
                                    'nn',
                                    'iter'
                                    ),
                        envir = environment()
                        )
#run Rhogen for the same sample size iter times
rps = parallel::parSapply(cl = cl,
                          X = rep(x = nn, 
                                  times = iter),
                          FUN = Rhogen)
  }
)

parallel::stopCluster(cl)
#find the bottom 5%
print(
  quantile(x = rps,
           p = 1-0.05)
)
#with 1e5 iterations appears accurate to ~3 decimal places
#with 1e6 iterations appears accurate to ~4 decimal places

hist(x = rps,
     breaks = 1e2,
     probability = TRUE,
     border = NA,
     xlim = c(0,1),
     main = paste0('H0 mean vector length, n = ', 
                   nn, 
                   ', ', iter, ' estimates'),
     xlab = 'mean vector length',
     ylab = 'probability density')



# Estimate beta distribution ----------------------------------------------

Beta_LP = function(x, par)
{
  alpha = par[1]
  beta = par[2]
  dd = dbeta(x = x,
             shape1 = alpha,
             shape2 = beta,
             log = TRUE)
  dd = ifelse(is.finite(dd),
              dd,
              -1e16)
  ll = sum(
            x = dd,
            na.rm = TRUE
          )
  return(-ll)
}

bp = optim(par = c(1.0,
              1.0),
      fn = Beta_LP,
      x = rps,
      method = "L-BFGS-B",
      lower = c(0,0),
      upper = c(1e3, 1e3)
      )

xx = seq(from = 0, to = 1, length.out = 1e3)

lines(x = xx, 
      y = dbeta(x = xx, 
                shape1 = bp$par[1],
                shape2 = bp$par[2]
      ),
      col = 'cyan4'
)


abline(v = c(quantile(x = rps,
                      p = 1-0.05),
             sqrt(-log(0.05)/nn),
             qbeta(p = 1-0.05,
                   shape1 = bp$par[1],
                   shape2 = bp$par[2])
              ),
       col = c(1,2, 'cyan4'),
       lty = c(1,1,3),
       lwd = c(1,1,3)
)
legend(x = 'topright',
       legend = c('p(>rho) = 0.05',
                  '√(-ln(0.05))/n',
                  'p(>beta) = 0.05'
                  ),
       col = c(1,2, 'cyan4'),
       lty = c(1,1, 3),
       lwd = c(3,3, 3),
       cex = 0.7
)
