########################################################################
#                                                                      #
#       Instructions                                                   #
#                                                                      #
########################################################################
# COPYRIGHT
#
# (c) Copyright Lawrence Joseph, 1994 - 1997.
#
# tt1, tt2, and tt3 are programs written by Lawrence Joseph,
# at the Division of Clinical Epidemiology, Department of Medicine,
# Montreal General Hospital.  These programs are an implementation of the
# manuscripts Bayesian Estimation of Disease Prevalence and the parameters
# of Diagnostic Tests in the Absence of a Gold Standard, by L. Joseph, T.
# Gyorkos, and L. Coupal, American Journal of Epidemiology,
# 1995;141:263-72.
#
# You are free to use these programs, for non-commercial purposes only,
# under two conditions:
#
# (1) This note is not to be removed;
# (2) Publications using tt1, tt2, or tt3 results should reference the
#     manuscript mentioned above.
#----------------------------------------------------------------------
# The manuscript mentioned above should be read carefully prior to using
# the program.
#----------------------------------------------------------------------
# This file contains the program tt1.gibbs and all required subroutines
# to calculate Bayesian posterior distributions via Gibbs sampling for
# the prevalence of a disease, and  sensitivity, specificity, positive
# and negative predictive values of a test for that disease, in the
# absence of a gold standard.  The program is written in S-PLUS version
# 3.1.  In order to run the program, type from the Splus prompt:
#
# tt1.gibbs(tp, tn, astart, cstart, sensstart, specstart, prevstart,
# alphaprev, betaprev, alphasens, betasens, alphaspec, betaspec, size,
#  throwaway, skip).
#
# The output is summarized using quantiles for each variable, although
# by changing line number 229, the full gibbs output is available.
#
# The parameters are defined as follows:
#
# tp == the total number of positive outcomes in the population tested by
# the (non-gold standard) test
#
# tn == the total number of negative outcomes in the population tested by
# the (non-gold standard) test
#
# astart == starting value for a, the unobserved number of true positives
#
# cstart == starting value for c, the unobserved number of false negatives
#
# sensstart == starting value for the sensitivity of the test
#
# specstart == starting value for the specificity of the test
#
# prevstart == starting value for the prevalencein the population
#
# alphaprev == first coefficient of the Beta prior distribution for the
# prevalence
#
# betaprev == second coefficient of the Beta prior distribution for the
# prevalence
#
# alphasens == first coefficient of the Beta prior distribution for the
# sensitivity
#
# betasens == second coefficient of the Beta prior distribution for the
# sensitivity
#
# alphaspec == first coefficient of the Beta prior distribution for the
# specificity
#
# betaspec == second coefficient of the Beta prior distribution for the
# specificity
#
# size == total number of Gibbs iterations
#
# throwaway == number of Gibbs iterations used for assessing convergence
#
# skip == step size for Gibbs iterates.  skip==1 means use all iterations,
# skip == 2 means use every second iterate, etc.
#
#---------------------------------------------------------------
#
# The general setup is:
#
#
#                   Truth (Unknown)
#                  |  +  |  -  |
#           T ---------------------------
#           e  +   |  a  |  b  | (a+b)
#           s ---------------------------
#           t  -   |  c  |  d  | (c+d)
#             ---------------------------
#                  |(a+c)|(b+d)|   N
#
# Note that only (a+b)=tp and (c+d)=tn are observed in the experiment.
#
#
#######################################################################
# After running the program, the output variables are:
#
# prev == the posterior prevalence
# sens == the posterior sensitivity
# spec == the posterior specificity
# ppv== the posterior positive predictive value
# npv == the posterior negative value
# lrp == the posterior likelihood ratio of a positive test
# lrn == the posterior likelihood ratio of a negative test
# last.values == the vector of last values of the Gibbs sampler,
# in case more iterations are needed.
#
#  A "q" in front of each variable indicates that the quantiles
#  are given.  For example, qprev indicates the posterior quantiles
#  for the prevalence.
#
######################################################################
#  Please let me know if you experience any problems while using     #
#  these functions:  Lawrence Joseph, joseph@binky.ri.mgh.mcgill.ca  #
######################################################################
#
tt1.a <-
  function(tp, sens, spec, prev)
  {
    if(tp == 0) {
      return(0)
    }
    p1 <- prev * sens
    p2 <- (1 - prev) * (1 - spec)
    p <- p1/(p1 + p2)
    nexta <- rbinom(1, tp, p)
    if(p == 0) {
      return(0)
    }
    if(p == 1) {
      return(tp)
    }
    return(nexta)
  }
tt1.c <-
  function(tn, sens, spec, prev)
  {
    if(tn == 0) {
      return(0)
    }
    p1 <- prev * (1 - sens)
    p2 <- (1 - prev) * spec
    p <- p1/(p1 + p2)
    nextc <- rbinom(1, tn, p)
    if(p == 0) {
      return(0)
    }
    if(p == 1) {
      return(tn)
    }
    return(nextc)
  }
tt1.sens <-
  function(a, c, alphasens, betasens)
  {
    return(rbeta(1, a + alphasens, c + betasens))
  }
tt1.spec <-
  function(tp, tn, a, c, alphaspec, betaspec)
  {
    return(rbeta(1, tn - c + alphaspec, tp - a + betaspec))
  }
tt1.prev<-function(tp, tn, sens, spec, a, c, alphaprev, betaprev)
{
  nextprev <- rbeta(1, a + c + alphaprev, tp - a + tn - c + betaprev)
  return(nextprev)
}
tt1.gibbs <-
  function(tp, tn, astart, cstart, sensstart, specstart, prevstart, alphaprev,
           betaprev, alphasens, betasens, alphaspec, betaspec, size, throwaway,
           skip)
  {
    throw <- throwaway + 1
    a.samp <- rep(-1, size)
    c.samp <- rep(-1, size)
    prev.samp <- rep(-1, size)
    sens.samp <- rep(-1, size)
    spec.samp <- rep(-1, size)
    ppv.samp <- rep(-1, size)
    npv.samp <- rep(-1, size)
    prev.samp[1] <- prevstart
    a.samp[1] <- astart
    c.samp[1] <- cstart
    sens.samp[1] <- sensstart
    spec.samp[1] <- specstart
    ppv.samp[1] <- a.samp[1]/tp
    npv.samp[1] <- (tn - c.samp[1])/tn
    for(i in 2:size) {
      a.samp[i] <- tt1.a(tp, sens.samp[i - 1], spec.samp[i - 1],
                         prev.samp[i - 1])
      c.samp[i] <- tt1.c(tn, sens.samp[i - 1], spec.samp[i - 1],
                         prev.samp[i - 1])
      sens.samp[i] <- tt1.sens(a.samp[i], c.samp[i], alphasens,
                               betasens)
      spec.samp[i] <- tt1.spec(tp, tn, a.samp[i], c.samp[i],
                               alphaspec, betaspec)
      prev.samp[i] <- tt1.prev(tp, tn, sens.samp[i], spec.samp[i],
                               a.samp[i], c.samp[i], alphaprev, betaprev)
    }
    ppv.samp<-sens.samp*prev.samp/(sens.samp*prev.samp + (1-spec.samp)*(1-prev.samp))
    npv.samp<-spec.samp*(1-prev.samp)/(spec.samp*(1-prev.samp) + prev.samp*(1-sens.samp))
    lrp.samp <- sens.samp/(1 - spec.samp)
    lrn.samp <- (1 - sens.samp)/(spec.samp)
    
    #
    #   Add lines below for automatic graphics of Gibbs sampler output
    #
    #        openlook()
    #        par(mfrow = c(3, 3))
    #        plot(1:throwaway, prev.samp[1:throwaway], type = "l")
    #        plot(throw:size, prev.samp[throw:size], type = "l")
    #        hist(a.samp[seq(throw, size, by = skip)])
    #        hist(c.samp[seq(throw, size, by = skip)])
    #        hist(sens.samp[seq(throw, size, by = skip)])
    #        hist(spec.samp[seq(throw, size, by = skip)])
    #        hist(ppv.samp[seq(throw, size, by = skip)])
    #        hist(npv.samp[seq(throw, size, by = skip)])
    #        hist(prev.samp[seq(throw, size, by = skip)])
    #
    # Add line below to return Gibbs output vectors, rather than summaries.
    #
    #
    #       return(a.samp, c.samp, sens.samp, spec.samp, prev.samp)
    #
    #
    qprev <- quantile(prev.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                                0.25, 0.5, 0.75, 0.95, 0.975))
    qsens <- quantile(sens.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                                0.25, 0.5, 0.75, 0.95, 0.975))
    qspec <- quantile(spec.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                                0.25, 0.5, 0.75, 0.95, 0.975))
    qppv <- quantile(ppv.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                              0.25, 0.5, 0.75, 0.95, 0.975))
    qnpv <- quantile(npv.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                              0.25, 0.5, 0.75, 0.95, 0.975))
    qlrp <- quantile(lrp.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                              0.25, 0.5, 0.75, 0.95, 0.975))
    qlrn <- quantile(lrn.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                              0.25, 0.5, 0.75, 0.95, 0.975))
    last.values <- c(a.samp[size], c.samp[size], sens.samp[size], spec.samp[size], prev.samp[size])
    list(qprev=qprev, qsens=qsens, qspec=qspec, qppv=qppv, qnpv=qnpv, qlrp=qlrp, qlrn=qlrn, last.values=last.values)
  }
mu.to.beta <-
  function(mu, sd)
  {
    var <- sd^2
    alpha <-  - (mu * (var + mu^2 - mu))/var
    beta <- ((mu - 1) * (var + mu^2 - mu))/var
    return(alpha, beta)
  }
beta.to.mu <-
  function(alpha, beta)
  {
    mean <- alpha/(alpha + beta)
    sd <- sqrt((alpha * beta)/((alpha + beta)^2 * (alpha + beta + 1)))
    return(mean, sd)
  }