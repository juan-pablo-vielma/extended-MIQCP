extended-MIQCP
==============

This is the code for the paper [Extended Formulations in Mixed Integer Conic Quadratic Programming](http://www.optimization-online.org/DB_HTML/2015/01/4736.html) by [Juan Pablo Vielma](http://web.mit.edu/jvielma/www/), [Iain Dunning](http://iaindunning.com/), [Joey Huchette](http://www.mit.edu/~huchette/) and [Miles Lubin](http://www.mit.edu/~mlubin/). Bellow are instructions on how to re-run the computational experiments and how to generate the tables and graphics in the paper. 

## Required Software

## Downloading the Code

## Running the Experiments

## Generating Tables

Code for generating the tables can be found in the **tables** folder. To generate all tables simply call the code with the name of the results file as the only argument. For instance, to generate the tables from the papers results call

``` $ julia createtables.jl ../results.csv ```

To generate a custom table use as arguments: the result file name, instance size, instance class (Mark, Short or Robust) and the list of solvers. The code then generates two files named **test_time.tex** and **test_quality.tex**, which contain the tables. For instance calling 

```$ julia createtables.jl ../results.csv 30 Mark CplexQCP GurobiQCP ```

generates the tables for the Classical instances for n=30 and for CPLEX's and Gurobi's QCP solvers.

## Generating Box Plots

The code for generating the box plots can be found in the **boxplot** folder. The file **makeboxplots.nb** contains [Mathematica](http://www.wolfram.com/mathematica/) code to generate the graphs. 

## Generating Performance Profiles




