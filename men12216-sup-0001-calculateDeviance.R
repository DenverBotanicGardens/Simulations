# calculateDeviance.R; Patrick G. Meirmans (2013) Non-convergence in Bayesian estimation of migration rates, Molecular Ecology Resources; Supplementary Material

# This script will calculate the Bayesian Deviance from the output of BayesAss (Wilson & Rannala 2003)
# For more information on the use of the deviance, see Faubet et al. (2007)

# To use this script, run BayesAss version 3 with the -t flag to produce a trace-file
# Then set the working directory of R to the folder with the output from BayesAss

# Change this value to the actual burnin used for your MCMC run
burnin = 1000000
# Change this value to the actual sampling interval used for your MCMC run
sampling.interval = 2000


# Read the data from the trace file
trace=read.table("BA3trace.txt", header=TRUE)

# Plotting the likelihoods is always a good idea
plot(trace$State,trace$LogProb, xlab="State", ylab="LogProb", type="l", lwd=2, col="firebrick3", font.lab=2)
# Draw a vertical line to indicate the end of the burnin
abline(v=burnin, col="grey70", lty=2)

# Calculate the deviance
range = (trace$State > burnin & trace$State %% sampling.interval == 0)
D = -2*mean(trace$LogProb[range])

# Print the result to the console
print(D)