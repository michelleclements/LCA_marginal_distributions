The following file contains an example of running the program tt2.gibbs.
Other programs are run in a similar fashion.

To obtain the output for the AJE paper, I ran the program listed below
in R.  The program was run 4 times using 4 different starting
values, to check for convergence, although only one run is shown here.
Sensitivity runs for the prior distribution may also be important,
depending on your application.  Answers will not be identical from run
to run, due to the Monte Carlo variance.

R program (for 2 by 2 table, other situations are similar):
  
  First run the command

> gibbs.sampler.out<-tt2.gibbs(38,87,2,35,5,5,2,30,0.5,0.5,0.5,0.5,
                               0.5,1,1,21.96,5.49,4.1,1.76,4.44,13.31,71.25,3.75,20500)

to perform the Gibbs sampler.

As described in the tt2.doc file, the first 4 parameters are the data,
38,87,2,35 from the 2 by 2 table (see AJE paper for source of this
                                  data).  The next 4 parameters, 5,5,2,30, are arbitrary starting values
for the number of true positive subjects in each cell.  Other starting
values can also be used, provided they are in the correct range. For
example, all are numbers between 0 and the cell size, eg, 38 for the
first entry, etc.  The next 5 parameters, 0.5,0.5,0.5,0.5,0.5 are
starting values for the two sets of sensitivities and specificities,
and the prevalence.  Again, these are arbitrary values, but must be
within the range 0 to 1.  Although they are theoretically arbitrary,
these starting values may affect rates of convergence.  Therefore in
practice, you may wish to use your "best guess" values for each
parameter. The next 10 parameters,
1,1,21.96,5.49,4.1,1.76,4.44,13.31,71.25,3.75, are the 10 parameters of
the 5 Beta prior distributions for the five parameters of interest, the
prevalence of disease, and the sensitivities and specificities of each
test.  For example, the beta prior parameters for the prevalence here
was a Beta(1,1) distribution, ie, a uniform prior.
Similarly, the next two parameters indicate a Beta(21.96,5.49) prior
for the sensitivity of test 1.  Please refer to my AJE paper for the
source of these prior distributions, and for references on how to
derive prior densities for your own particular situation.  The built-in
subroutines beta.to.mu and mu.to.beta may be helpful in this regard.
For example:
  
  > beta.to.mu(21.96,5.49)
$mean:
  [1] 0.8

$sd:
  [1] 0.07499268

>  mu.to.beta(0.8,0.07499268)
$alpha:
  [1] 21.96

$beta:
  [1] 5.49

This allows conversion between beta parameters and their means and
standard deviations, and vice versa, as indicated by the example.

The last parameter, 20500, indicates the number of iterations of the
Gibbs sampler to run.  Depending on your computer setup, this number
may or may not be feasible in one run.  I would suggest starting with
small values, say 100, and gradually increasing the number until you
find the limits of your computer setup.

The line above creates output which can be summarized by next running the
command:
  
  > tt2.sum(gibbs.sampler.out,500,1)

The first parameter is the file named in the first command above, which
stores the vectors.  The second command is the number of iterations
before convergence.  Inference is then based on all of the output after
the first 500 iterations, since the skip parameter is 1.  The output
you should obtain is:
  
  $size:
  [1] 20500

$qprev:
  2.5%      5.0%     25.0%     50.0%     75.0%    95.0%     97.5%
0.5260121 0.5750092 0.7097868 0.7751134 0.830716 0.9004895 0.9230528

$qsens1:
  2.5%      5.0%     25.0%     50.0%     75.0%     95.0%     97.5%
0.7900652 0.8067705 0.8568488 0.8874314 0.9144861 0.9458535 0.9544222

$qspec1:
  2.5%      5.0%     25.0%     50.0%     75.0%     95.0%     97.5%
0.3731874 0.4141996 0.5772343 0.7041374 0.8210112 0.9343055 0.9570939

$qppv1:
  2.5%      5.0%    25.0%     50.0%     75.0%    95.0%    97.5%
0.6285226 0.6900756 0.8487309 0.9189233 0.9603374 0.9887083 0.9931422

$qnpv1:
  2.5%      5.0%     25.0%     50.0%     75.0%     95.0%     97.5%
0.273309 0.3407579 0.5270505 0.6340492 0.7269634 0.8297648 0.8571152

$qsens2:
  2.5%      5.0%     25.0%     50.0%     75.0%     95.0%     97.5%
0.2226895 0.2347676 0.2741186 0.3048363 0.3390376 0.399007 0.4248368

$qspec2:
  2.5%      5.0%     25.0%     50.0%     75.0%     95.0%    97.5%
0.9071047 0.9173849 0.9453828 0.9605303 0.9726511 0.9854504 0.9886054

$qppv2:
  2.5%      5.0%    25.0%     50.0%     75.0%     95.0%     97.5%
0.8735122 0.8978902 0.9447286 0.9647505 0.9785255 0.9902885 0.9929387

$qnpv2:
  2.5%      5.0%     25.0%     50.0%     75.0%     95.0%     97.5%
0.09755736 0.1269252 0.2156531 0.2844137 0.3663822 0.529226 0.5874898

These are all interpreted as quantiles for the various parameters.
For example, the output line

$qprev:
  2.5%      5.0%     25.0%     50.0%     75.0%    95.0%     97.5%
0.5260121 0.5750092 0.7097868 0.7751134 0.830716 0.9004895 0.9230528

says that the 2.5% quantile of the prevalence is estimated to be
0.5260121, and the uper 97.5% quantile is estimated to be 0.9230528.
Therefore a 95% credible interval is ( 0.5260121, 0.9230528 ). The
median (50 percentile) is 0.7751134, and the inter-quartile range is
( 0.7097868, 0.830716 ).


The output is given in the order:
  
  $size:    the number of gibbs iterations used for inference
$qprev:   prevalence of the condition
$qsens1:  sensitivity of test 1
$qspec1:  specificity of test 1
$qppv1:   positive predictive value of test 1
$qnpv1:   negative predictive value for test 1
$qsens2:  sensitivity of test 2
$qspec2:  specificity of test 2
$qppv2:   positive predictive value for test 2
$qnpv2:   negative predictive value for test 2
