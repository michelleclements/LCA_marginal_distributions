########################################################################
#                                                                      #
#       Instructions                                                   #
#                                                                      #
########################################################################
# COPYRIGHT
#
# (c) Copyright Lawrence Joseph, 1994 - 1997
#
# tt1, tt2, and tt3 are programs written by Lawrence Joseph,
# at the Division of Clinical Epidemiology, Department of Medicine,
# Montreal General Hospital.  These programs are an implementation of the
# manuscript Bayesian Estimation of Disease Prevalence and the Parameters
# of Diagnostic Tests in the Absence of a Gold Standard, by L. Joseph,
# T.  Gyorkos, and L. Coupal, American Journal of Epidemiology,
# 1995;141:263-72.
#
# You are free to use these programs, for non-commercial purposes only,
# under two conditions:
#
# (1) This note is not to be removed;
# (2) Publications using tt1, tt2, or tt3 results should reference the
#     manuscript mentioned above.
# ----------------------------------------------------------------------
# The manuscript mentioned above should be read carefully prior to using
# the program.
#----------------------------------------------------------------------
#
# This file contains the program tt2.gibbs and all required subroutines
# to calculate Bayesian posterior distributions via Gibbs sampling for
# the prevalence of a disease, and  sensitivity, specificity, positive
# and negative predictive values for each of two tests for that disease,
# in the absence of a gold standard.  The program is writtem in S-PLUS
# version 3.1.  In order to run the program, type from the
# Splus prompt:
#
# tt2.out<-tt2.gibbs(u, v, w, x, y1start, y2start,y3start, y4start,
# sens1start, spec1start, sens2start, spec2start, prevstart, alphaprev,
# betaprev, alphasens1, betasens1, alphaspec1, betaspec1, alphasens2,
# betasens2, alphaspec2, betaspec2, size).
#
# We can combine different runs (continued from the previous or not) by
# typing
#
# tt2.out3<-tt2.comb(tt2.out1, tt2.out2, overlap)
#
# The output is summarized (using quantiles for each variable) by typing
#
# tt2.sum(tt2.out3, throwaway, skip)
#
#  Likelihood ratios are then separately summarized by typing
#
# lr.sum(tt2.out3, throwaway, skip)
#
# The parameters are defined as follows:
#
# u, v, w, and x == see the table of data below
#
# y(1234)start == starting value for the unobserved number of true
#                 positives in u, v, w, or x, respectively.
#
# sens(12)start == starting value for the sensitivity of the test
#
# spec(12)start == starting value for the specificity of the test
#
# prevstart == starting value for the prevalence in the population
#
# alphaprev == first coefficient of the Beta prior distribution for the
#              prevalence
#
# betaprev == second coefficient of the Beta prior distribution for the
#            prevalence
#
# alphasens(12) == first coefficient of the Beta prior distribution for
#                  the sensitivity
#
# betasens(12) == second coefficient of the Beta prior distribution for
#                 the sensitivity
#
# alphaspec(12) == first coefficient of the Beta prior distribution for
#                  the specificity
#
# betaspec(12) == second coefficient of the Beta prior distribution for
#                 the specificity
#
# size == total number of Gibbs iterations
#
# throwaway == number of Gibbs iterations used for assessing convergence
#
# skip == step size for Gibbs iterates.  skip == 1 means use all
#         iterations, skip == 2 means use every second iterate, etc.
#
# overlap == Use overlap == 1 if runs are not continuous to each other,
#            or use overlap == 2 if they are consecutive.  Can also be
#            used for throwing away from the second of two independent
#            runs.
#--------------------------------------------------------------
# The general setup is that we observe the entirety of the following
# table:
#
#
#                   Test 2
#                  |  +  |  -  |
#           T ---------------------------
#           e  +   |  u  |  v  | (u+v)
#           s ---------------------------
#           t  -   |  w  |  x  | (w+x)
#           1 ---------------------------
#                  |(u+w)|(v+x)|   N
#
# However, there are latent data y1, y2, y3, and y4 which represent the
# unobserved number of true positives in each of the above cells, u, v,
# w, and x respectively.
#
#######################################################################
# After running the program, the output variables are:
#
# prev == the posterior prevalence
# sens1 == the posterior sensitivity from the first test
# spec1 == the posterior specificity from the first test
# ppv1 == the posterior positive predictive value from the first test
# npv1 == the posterior negative value from the first test
# lrp1 == the posterior likelihood ratio of a positive test
#          from the first test
# lrn1 == the posterior likelihood ratio of a negative test
#          from the first test
# sens2 == the posterior sensitivity from the second test
# spec2 == the posterior specificity from the second test
# ppv2 == the posterior positive predictive value from the second test
# npv2 == the posterior negative value from the second test
# lrp2 == the posterior likelihood ratio of a positive test
#          from the second test
# lrn2 == the posterior likelihood ratio of a negative test
#          from the second test
# lrpp == the posterior likelihood ratio when both tests are
#          positive
# lrpn == the posterior likelihood ratio when the first test is
#          positive and the second test is negative
# lrnp == the posterior likelihood ratio when the first test is
#          negative and the second test is positive
# lrnn == the posterior likelihood ratio when both tests are
#          negative
# last.values == the vector of last values of the Gibbs sampler,
# in case more iterations are needed.
#
#  A "q" in front of each variable indicates that the quantiles
#  are given.  For example, qprev indicates the posterior quantiles
#  for the prevalence.  This is seen in the output of the tt2.sum
#  and lr.sum commands.
#
######################################################################
#  Please let me know if you experience any problems while using     #
#  these functions:  Lawrence Joseph, joseph@binky.ri.mgh.mcgill.ca  #
######################################################################
#
tt2.y1 <-
  function(u, prev, sens1, spec1, sens2, spec2)
  {       if(u == 0) {
    return(0)
  }
    p1 <- prev * sens1 * sens2
    p2 <- (1 - prev) * (1 - spec1) * (1 - spec2)
    p <- p1/(p1 + p2)
    nexty1 <- rbinom(1, u, p)
    if(p == 0) {
      return(0)
    }
    if(p == 1) {
      return(u)
    }
    return(nexty1)
  }
tt2.y2 <-
  function(v, prev, sens1, spec1, sens2, spec2)
  {       if(v == 0) {
    return(0)
  }
    
    p1 <- prev * sens1 * (1 - sens2)
    p2 <- (1 - prev) * (1 - spec1) * spec2
    p <- p1/(p1 + p2)
    nexty2 <- rbinom(1, v, p)
    if(p == 0) {
      return(0)
    }
    if(p == 1) {
      return(v)
    }
    return(nexty2)
  }
tt2.y3 <-
  function(w, prev, sens1, spec1, sens2, spec2)
  {       if(w == 0) {
    return(0)
  }
    
    p1 <- prev * (1 - sens1) * sens2
    p2 <- (1 - prev) * spec1 * (1 - spec2)
    p <- p1/(p1 + p2)
    nexty3 <- rbinom(1, w, p)
    if(p == 0) {
      return(0)
    }
    if(p == 1) {
      return(w)
    }
    return(nexty3)
  }
tt2.y4 <-
  function(x, prev, sens1, spec1, sens2, spec2)
  {       if(x == 0) {
    return(0)
  }
    
    p1 <- prev * (1 - sens1) * (1 - sens2)
    p2 <- (1 - prev) * spec1 * spec2
    p <- p1/(p1 + p2)
    nexty4 <- rbinom(1, x, p)
    if(p == 0) {
      return(0)
    }
    if(p == 1) {
      return(x)
    }
    return(nexty4)
  }
tt2.prev <-
  function(u, v, w, x, y1, y2, y3, y4, alphaprev, betaprev)
  {
    nextprev <- rbeta(1, y1 + y2 + y3 + y4 + alphaprev, u + v + w + x - (y1 +
                                                                           y2 + y3 + y4) + betaprev)
    return(nextprev)
  }
tt2.sens1 <-
  function(y1, y2, y3, y4, alphasens1, betasens1)
  {
    nextsens1 <- rbeta(1, y1 + y2 + alphasens1, y3 + y4 + betasens1)
    return(nextsens1)
  }
tt2.sens2 <-
  function(y1, y2, y3, y4, alphasens2, betasens2)
  {
    nextsens2 <- rbeta(1, y1 + y3 + alphasens2, y2 + y4 + betasens2)
    return(nextsens2)
  }
tt2.spec1 <-
  function(u, v, w, x, y1, y2, y3, y4, alphaspec1, betaspec1)
  {
    nextspec1 <- rbeta(1, w + x - (y3 + y4) + alphaspec1, u + v - (y1 + y2) +
                         betaspec1)
    return(nextspec1)
  }
tt2.spec2 <-
  function(u, v, w, x, y1, y2, y3, y4, alphaspec2, betaspec2)
  {
    nextspec2 <- rbeta(1, v + x - (y2 + y4) + alphaspec2, u + w - (y1 + y3) +
                         betaspec2)
    return(nextspec2)
  }
tt2.gibbs <-
  function(u, v, w, x, y1start, y2start, y3start, y4start, sens1start, spec1start,
           sens2start, spec2start, prevstart, alphaprev, betaprev, alphasens1,
           betasens1, alphaspec1, betaspec1, alphasens2, betasens2, alphaspec2,
           betaspec2, size)
  {
    y1.samp <- rep(-1, size)
    y2.samp <- rep(-1, size)
    y3.samp <- rep(-1, size)
    y4.samp <- rep(-1, size)
    prev.samp <- rep(-1, size)
    sens1.samp <- rep(-1, size)
    spec1.samp <- rep(-1, size)
    ppv1.samp <- rep(-1, size)
    npv1.samp <- rep(-1, size)
    sens2.samp <- rep(-1, size)
    spec2.samp <- rep(-1, size)
    ppv2.samp <- rep(-1, size)
    npv2.samp <- rep(-1, size)
    prev.samp[1] <- prevstart
    y1.samp[1] <- y1start
    y2.samp[1] <- y2start
    y3.samp[1] <- y3start
    y4.samp[1] <- y4start
    sens1.samp[1] <- sens1start
    spec1.samp[1] <- spec1start
    sens2.samp[1] <- sens2start
    spec2.samp[1] <- spec2start
    for(i in 2:size) {
      y1.samp[i] <- tt2.y1(u, prev.samp[i - 1], sens1.samp[i - 1],
                           spec1.samp[i - 1], sens2.samp[i - 1], spec2.samp[i - 1]
      )
      y2.samp[i] <- tt2.y2(v, prev.samp[i - 1], sens1.samp[i - 1],
                           spec1.samp[i - 1], sens2.samp[i - 1], spec2.samp[i - 1]
      )
      y3.samp[i] <- tt2.y3(w, prev.samp[i - 1], sens1.samp[i - 1],
                           spec1.samp[i - 1], sens2.samp[i - 1], spec2.samp[i - 1]
      )
      y4.samp[i] <- tt2.y4(x, prev.samp[i - 1], sens1.samp[i - 1],
                           spec1.samp[i - 1], sens2.samp[i - 1], spec2.samp[i - 1]
      )
      sens1.samp[i] <- tt2.sens1(y1.samp[i], y2.samp[i], y3.samp[i],
                                 y4.samp[i], alphasens1, betasens1)
      sens2.samp[i] <- tt2.sens2(y1.samp[i], y2.samp[i], y3.samp[i],
                                 y4.samp[i], alphasens2, betasens2)
      spec1.samp[i] <- tt2.spec1(u, v, w, x, y1.samp[i], y2.samp[i],
                                 y3.samp[i], y4.samp[i], alphaspec1, betaspec1)
      spec2.samp[i] <- tt2.spec2(u, v, w, x, y1.samp[i], y2.samp[i],
                                 y3.samp[i], y4.samp[i], alphaspec2, betaspec2)
      prev.samp[i] <- tt2.prev(u, v, w, x, y1.samp[i], y2.samp[i],
                               y3.samp[i], y4.samp[i], alphaprev, betaprev)
    }
    
    ppv1.samp <- sens1.samp*prev.samp/(sens1.samp*prev.samp + (1-prev.samp)*(1-spec1.samp))
    npv1.samp <- spec1.samp*(1-prev.samp)/(spec1.samp*(1-prev.samp) + prev.samp*(1-sens1.samp))
    ppv2.samp <- sens2.samp*prev.samp/(sens2.samp*prev.samp + (1-prev.samp)*(1-spec2.samp))
    npv2.samp <- spec2.samp*(1-prev.samp)/(spec2.samp*(1-prev.samp) + prev.samp*(1-sens2.samp))
    last.values <- c(y1.samp[size], y2.samp[size], y3.samp[size], y4.samp[
      size], sens1.samp[size], spec1.samp[size], sens2.samp[size],
      spec2.samp[size], prev.samp[size])
    list(y1.samp=y1.samp, y2.samp=y2.samp, y3.samp=y3.samp, y4.samp=y4.samp, sens1.samp=sens1.samp, spec1.samp=spec1.samp,
         ppv1.samp=ppv1.samp, npv1.samp=npv1.samp, sens2.samp=sens2.samp, spec2.samp=spec2.samp, ppv2.samp=ppv2.samp,
         npv2.samp=npv2.samp, prev.samp=prev.samp, last.values=last.values)
  }
tt2.comb <-
  function(a.out, b.out, overlap)
  {
    size <- length(b.out$prev.samp)
    y1.samp <- c(a.out$y1.samp, b.out$y1.samp[overlap:size])
    y2.samp <- c(a.out$y2.samp, b.out$y2.samp[overlap:size])
    y3.samp <- c(a.out$y3.samp, b.out$y3.samp[overlap:size])
    y4.samp <- c(a.out$y4.samp, b.out$y4.samp[overlap:size])
    sens1.samp <- c(a.out$sens1.samp, b.out$sens1.samp[overlap:size])
    spec1.samp <- c(a.out$spec1.samp, b.out$spec1.samp[overlap:size])
    ppv1.samp <- c(a.out$ppv1.samp, b.out$ppv1.samp[overlap:size])
    npv1.samp <- c(a.out$npv1.samp, b.out$npv1.samp[overlap:size])
    sens2.samp <- c(a.out$sens2.samp, b.out$sens2.samp[overlap:size])
    spec2.samp <- c(a.out$spec2.samp, b.out$spec2.samp[overlap:size])
    ppv2.samp <- c(a.out$ppv2.samp, b.out$ppv2.samp[overlap:size])
    npv2.samp <- c(a.out$npv2.samp, b.out$npv2.samp[overlap:size])
    prev.samp <- c(a.out$prev.samp, b.out$prev.samp[overlap:size])
    return(y1.samp, y2.samp, y3.samp, y4.samp, sens1.samp, spec1.samp,
           ppv1.samp, npv1.samp, sens2.samp, spec2.samp, ppv2.samp,
           npv2.samp, prev.samp)
  }
tt2.sum <-
  function(tt2.out, throwaway, skip)
  {
    throw <- throwaway + 1
    size <- length(tt2.out$prev.samp)
    #
    #   Add following lines for graphical output
    #
    #        openlook()
    #        par(mfrow = c(4, 4))
    #        plot(1:throwaway, tt2.out$prev.samp[1:throwaway], type = "l")
    #        plot(throw:size, tt2.out$prev.samp[throw:size], type = "l")
    #        hist(tt2.out$y1.samp[seq(throw, size, by = skip)])
    #        hist(tt2.out$y2.samp[seq(throw, size, by = skip)])
    #        hist(tt2.out$y3.samp[seq(throw, size, by = skip)])
    #        hist(tt2.out$y4.samp[seq(throw, size, by = skip)])
    #        hist(tt2.out$sens1.samp[seq(throw, size, by = skip)])
    #        hist(tt2.out$spec1.samp[seq(throw, size, by = skip)])
    #        hist(tt2.out$ppv1.samp[seq(throw, size, by = skip)])
    #        hist(tt2.out$npv1.samp[seq(throw, size, by = skip)])
    #        hist(tt2.out$sens2.samp[seq(throw, size, by = skip)])
    #        hist(tt2.out$spec2.samp[seq(throw, size, by = skip)])
    #        hist(tt2.out$ppv2.samp[seq(throw, size, by = skip)])
    #        hist(tt2.out$npv2.samp[seq(throw, size, by = skip)])
    #        hist(tt2.out$prev.samp[seq(throw, size, by = skip)])
    qprev <- quantile(tt2.out$prev.samp[seq(throw, size, by = skip)], c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
    qsens1 <- quantile(tt2.out$sens1.samp[seq(throw, size, by = skip)], c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
    qspec1 <- quantile(tt2.out$spec1.samp[seq(throw, size, by = skip)], c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
    qppv1 <- quantile(tt2.out$ppv1.samp[seq(throw, size, by = skip)], c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
    qnpv1 <- quantile(tt2.out$npv1.samp[seq(throw, size, by = skip)], c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
    qsens2 <- quantile(tt2.out$sens2.samp[seq(throw, size, by = skip)], c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
    qspec2 <- quantile(tt2.out$spec2.samp[seq(throw, size, by = skip)], c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
    qppv2 <- quantile(tt2.out$ppv2.samp[seq(throw, size, by = skip)], c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
    qnpv2 <- quantile(tt2.out$npv2.samp[seq(throw, size, by = skip)], c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
    list(size=size, qprev=qprev, qsens1=qsens1, qspec1=qspec2, qppv1=qppv1, qnpv1=qnpv1, qsens2=qsens2,
         qspec2=qspec2,qppv2=qppv2, qnpv2=qnpv2)
  }
lr.sum <-
  function(tt2.out, throwaway, skip)
  {
    throw <- throwaway + 1
    size <- length(tt2.out$prev.samp)
    lrp1.samp <- tt2.out$sens1.samp/(1 - tt2.out$spec1.samp)
    lrn1.samp <- (1 - tt2.out$sens1.samp)/(tt2.out$spec1.samp)
    lrp2.samp <- tt2.out$sens2.samp/(1 - tt2.out$spec2.samp)
    lrn2.samp <- (1 - tt2.out$sens2.samp)/(tt2.out$spec2.samp)
    lrpp.samp <- lrp1.samp * lrp2.samp
    lrpn.samp <- lrp1.samp * lrn2.samp
    lrnp.samp <- lrn1.samp * lrp2.samp
    lrnn.samp <- lrn1.samp * lrn2.samp
    lrp1.mean <- mean(lrp1.samp)
    lrn1.mean <- mean(lrn1.samp)
    lrp2.mean <- mean(lrp2.samp)
    lrn2.mean <- mean(lrn2.samp)
    lrpp.mean <- mean(lrpp.samp)
    lrpn.mean <- mean(lrpn.samp)
    lrnp.mean <- mean(lrnp.samp)
    lrnn.mean <- mean(lrnn.samp)
    qlrp1 <- quantile(lrp1.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                                0.25, 0.5, 0.75, 0.95, 0.975))
    qlrn1 <- quantile(lrn1.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                                0.25, 0.5, 0.75, 0.95, 0.975))
    qlrp2 <- quantile(lrp2.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                                0.25, 0.5, 0.75, 0.95, 0.975))
    qlrn2 <- quantile(lrn2.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                                0.25, 0.5, 0.75, 0.95, 0.975))
    qlrpp <- quantile(lrpp.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                                0.25, 0.5, 0.75, 0.95, 0.975))
    qlrpn <- quantile(lrpn.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                                0.25, 0.5, 0.75, 0.95, 0.975))
    qlrnp <- quantile(lrnp.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                                0.25, 0.5, 0.75, 0.95, 0.975))
    qlrnn <- quantile(lrnn.samp[seq(throw, size, by = skip)], c(0.025, 0.05,
                                                                0.25, 0.5, 0.75, 0.95, 0.975))
    #
    #
    #      Add the following lines for graphics (histograms)
    #
    #
    #
    #
    #        par(mfrow=c(4,2))
    #        hist(lrp1.samp)
    #        hist(lrn1.samp)
    #        hist(lrp2.samp)
    #        hist(lrn2.samp)
    #        hist(lrpp.samp)
    #        hist(lrpn.samp)
    #        hist(lrnp.samp)
    #        hist(lrnn.samp)
    #
    #
    #
    #
    return(size, qlrp1, lrp1.mean, qlrn1, lrn1.mean, qlrp2, lrp2.mean,
           qlrn2,  lrn2.mean, qlrpp,  lrpp.mean, qlrpn, lrpn.mean, qlrnp, lrnp.mean,
           qlrnn, lrnn.mean)
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
