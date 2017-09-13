# Simulate population genetics in a metapopulation of 115 demes, corresponding to Mediterranean MPAs
# In this file you will find all what you need to perform the calculations presented in the manuscript
# This file is divided into three parts
# 1) Creating the initial dataset
# 2) Performing the simulation
# 3) Analysing the results


# Clear the workspace
rm(list=ls())

# Load the MetaPopGen package
library(MetaPopGen)

#########################################################################################################################
#																				#
# Part 1: Creating the initial dataset	#
#																				#
#########################################################################################################################


	# General parameters
		n <- 115	# Number of demes
		l <- 2		# Number of alleles
		m <- 3		# Nmber of genotypes
		z <- 5		# Number of age-classes
		tmax <- 200	# Number of years of simulation


	# Initializes variables
	# This initialises the variables using the correct sizes for each dimension
		sigma_M <- array(NA,dim=c(m,n,z,tmax))
		sigma_F <- array(NA,dim=c(m,n,z,tmax))
		phi_M <- array(0,dim=c(m,n,z,tmax))
		phi_F <- array(0,dim=c(m,n,z,tmax))


	# This assigns values to the variables
		sigma_M[,,,] <- 0.8
		sigma_F[,,,] <- 0.8
		phi_M[,,3:5,] <- 1000000
		phi_F[,,3:5,] <- 200000


	# Load dispersal matrix
	# This file contains a (115*115) matrix of dispersal probabilities
	# issued from the biophysical simulations of Andrello et al. (2013, Plos One)
		load("Dispersal.RData")
								
	# Here the dispersal matrix is repeated for each year of simulation
	# i.e. we use the same dispersal matrix for each year
		delta <- array(dispersal,dim=c(n,n,tmax)) 
								


	# Initial population structure
	# Genotypes number are assigned at random, by drawing the number of individuals of each genotype from a uniform distribution between [0,300]
		N1_M <- array(NA,dim=c(m,n,z))		# Initial number of male individuals
		N1_F <- array(NA,dim=c(m,n,z))		# Initial number of female individuals
		for (deme in 1 : n){
			for (age in 1 : z) {
				for (genotype in 1 : m) {
					N1_M[genotype,deme,age] <- round(runif(1)*300)
					N1_F[genotype,deme,age] <- round(runif(1)*300)
				}
			}
		}


	# Mutation matrix
		mu <- array(1e-6,dim=c(l,l))
		mu[1,1] <- 1 - mu[2,1]
		mu[2,2] <- 1 - mu[1,2]


	# Maximum recruitment
		kappa0 <- array(1000,c(n,tmax))


	# Save dataset
		save(sigma_M,sigma_F,N1_M,N1_F,phi_M,phi_F,mu,delta,kappa0,file="Dataset.Mediterranean.MPAs.RData")


############################################### End of Part 1 ##########################################################


#########################################################################################################################
#																				#
# Part 2: Performing the simulation			#
#																				#
#########################################################################################################################

	# Clear the workspace
		rm(list=ls())

	# Load dataset (created in Part 1)
		load("Dataset.Mediterranean.MPAs.RData")

	# Simulation
	# We run 30 replicates 
	# The simulation results will be saved on disk, there will be one directory for each replicate, named with date-time
		nrepl <- 30

	# Create the directory for the results and set it as working directory
		dir.create(paste(getwd(),"/Simulation",sep=""))
		setwd(paste(getwd(),"/Simulation",sep=""))
	
	# Run the simulations
	for (i in 1 : nrepl) {
	  
	  cat("Replicate",i,"\n") # Just to prnt on screen the advancement of the simulation
	  
	  sim.metapopgen.dioecious(
		input.type="array",
		N1_M = N1_M, N1_F = N1_F,
		sigma_M = sigma_M, sigma_F = sigma_F,
		phi_M = phi_M, phi_F = phi_F,
		mu = mu, delta = delta, kappa0 = kappa0,
		save.res = T)
	}

############################################### End of Part 2 ##########################################################




#########################################################################################################################
#  																			#
# Part 3: Analysing the results					#
#																				#
#########################################################################################################################

	# Be sure that the working directory is the one where the results were saved; if you are running Part 3 just after Part 2, you are automatically in the right directory

	# Clear the workspace
		rm(list=ls())

	# General parameters
		n     <- 115	# Number of demes
		nrepl <- 30   # Number of replicates
		tmax <- 200   # Number of years of simulation

	# Part 3.1 : ####################################################
	# Calculate global FST for each replicate at each time step
	#################################################################

		# The female population is considered

		# Get the list of simulation directories (one directory per replicate, named after the date and time of simulation)
			list.dir <- dir()

		# Create the array for fst
			fst <- array(NA,dim=c(nrepl,tmax))
		
		# Calculate Fst
			for (i in 1 : nrepl) {
			  
			  dir.name <- list.dir[i]
			  
			  # Show the advancement of the calculation
			  print(dir.name)
			  flush.console()
			  
			  for (t in 1 : tmax){
				fst[i,t] <- fst.global.dioecious(dir.name,"F",t)
			  }
			}

		# Save fst array
			save(fst,file="Global.fst.in.time.RData")


		# Plot fst

			#tiff(filename = "Figure 4.tiff",
			#     width = 160, height = 160, units = "mm",
			#     compression = "lzw", res = 300)

			plot(fst[1,],ylim=c(0,0.02),type="l",xlab="Time",ylab="Global FST")
			for(i in 2: 30){
			  matplot(fst[i,],add=T,type="l")
			}
			matplot(colMeans(fst),add=T,type="l",col="red",lwd=2)
			#dev.off()


		# Calculate the final mean fst
			colMeans(fst)[tmax]

		# Calculate the trend in mean FST with time
			(colMeans(fst)[tmax] - colMeans(fst)[20]) / (tmax-20)

		# Calculate the variance among replicates and its increase with time
			plot(apply(fst,2,var),xlab="Time",ylab="Variance in FST among replicates")
			summary(lm(apply(fst,2,var)~seq(1,tmax)))



	# Part 3.2 : ####################################################
	# Calculate pairwise FST at the end of the simulation (t = 200)
	#################################################################


		# Create a big array that will contain the 30 pairwise fst matrices (one matrix per replicate)
		# Each matrix is a 115*115 matrix
			fst.matrices <- array(NA,dim=c(n,n,nrepl))
		
		# Calculate pairwise fst
			for (k in 1 : nrepl) {
			  
			  dir.name <- list.dir[k]
			  
			  # Show the advancement of the calculation
			  print(dir.name)
			  flush.console()
			  
			  # Create a pairwise fst matrix
			  fst.matrix <- array(NA,dim=c(n,n))
			  
			  for (i in 1 : n) {
				for (j in 1 : n) {
				  
				  # This calculates pairwise fst between the selected demes for age class one at the end of the simulation, between females
				  fst.matrix[i,j] <- fst.pairwise.dioecious(dir.name,i,j,1,1,tmax,tmax,"F","F")    
				  
				}
			  }
			  fst.matrices[,,k] <- fst.matrix
			}


		# Calculate the mean pairwise fst matrix
			fst.matrix <- array(NA,dim=c(n,n))
			for (i in 1 : n){
			  for (j in 1 : n){
				fst.matrix[i,j] <- mean(fst.matrices[i,j,])
			  }
			}

		# Plot the mean pairwise fst matrix
			#tiff(filename = "Figure 5.tiff",
			#     width = 160, height = 160, units = "mm",
			#     compression = "lzw", res = 300)
			library(fields)
			colors <- tim.colors(100)
			levelplot(fst.matrix,col.regions=colors,xlab="Deme",ylab="Deme",main="FST")
			#dev.off()

		# Statistics on pairwise fst
			quantile(fst.matrix,na.rm=T)



	# Part 3.3 : ####################################################
	# Test relationship between pairwise FST and a) spatial distance and b) probability of dispersal
	#################################################################

		setwd(paste(getwd(),"..",sep="/"))
		library(ape)


		# Load pairwise distances between demes
			load("Spatial.distance.RData")

		# Remove NaN from the fst matrix and the corresponding demes in the distance matrix
			fst.no.NaN      <- fst.matrix[-c(6,75,78),-c(6,75,78)] 
			sp.dist.no.NaN  <- sp.dist.matrix[-c(6,75,78),-c(6,75,78)]

		# Test hypothesis: fst increases with spatial distance using Mantel test
			mantel.test(sp.dist.no.NaN, fst.no.NaN,nperm=9999, alternative="greater",graph=T)


		# Load dispersal probabilities
			load("Dispersal.RData")

		# Remove the demes corresponding to NaN in the dispersal matrix
			dispersal.no.NaN  <- dispersal[-c(6,75,78),-c(6,75,78)]

		# Test hypothesis: fst decreases with increasing dispersal probability using Mantel test
			mantel.test(dispersal.no.NaN, fst.no.NaN,nperm=9999, alternative="less",graph=T)


############################################### End of Part 3 ##########################################################
