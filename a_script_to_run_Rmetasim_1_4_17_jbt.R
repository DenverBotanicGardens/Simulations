# Building an rmetasim-running script to explore Phacelia populations of Northern Colorado

# Jessie Berta-Thompson & Michelle DePrenger-Levin
# Wed Jan 4, 2016

# Source materials (occasionally quoted)
# rmetasim code and documentation, vignettes, Metasim publication Strand 2002, Kernelpop publication and documentation, Strand and Niehaus 2007


# This script contains a workflow for running rmetasim for
# 3 populations of Phacelia with light migration between them, 
# with rough informed guesses for parameters to be explored in more detail later.
# but this script is designed as a tool so that parameters can be easily changed for subsequent runs


# To start with a clean R slate, clear environment:

rm(list =ls())

# load rmetasim (if installed)
library(rmetasim)


# Create an rmetasim landscape
# Parameters for rmetasim simulations are stored in a landscape object (an organized list of lists and variables)
# We will fill in 6 pieces of the landscape in order (earlier ones influence later ones)

# "The six main sections of the landscape:
# intparam  - integer parameters describing landscape
# switchparam -  boolean switches turning features on and off
# floatparam - floating point paramters
# demography - lists describing the demography of the simulation (complex set of parameters)
# loci lists - describing the genetic components of the simulation
# individuals - a matrix that contains all individuals, their demographic states, ages and ids. It also contains the genetic information.
# The last three are values can change through the course of the simulation. Especially the last two"


#start by making an empty landscape object
rland <- landscape.new.empty()

# for readability & repeat use, using long form variable names to set parameters, then passing them to rmetasim functions.
# comment structure: # explain parameter choice for this run, default values in [], "vignette/publication parameter definitions"


# Landscape Building Part 1. Fill in parameters with integer values. 
# The function landscape.new.intparam(rland, ...) sets the values of a set of parameters that are coded with integers

habitats <- 3             # 3 Phacelia populations                                                                 [1]     "the number of different habitats or subpopulations within the landscape"                                             
stages <- 3               # 3 life stages: seedlings, still juvenile year-old rosettes, flowering                  [1]     "the number of stages in the life cycle of the organism; also the size of the S, R, and M matricies for local demography
locusnum <- 12            # 12 microsatellite loci used in genotyping study                                        [0]     "the number of different loci currently implemented for the simulation"
numepochs <- 0            # one migration scenario over time                                                       [0]     "the number of different migration scenerios to occur during the simulation"
currentgen <- 0           # start at generation 0                                                                  [0]     "the current generation the simulation has reached"                              
currentepoch <- 0         # start at epoch 0                                                                       [0]     "the current epoch the simulation has reached"                                   
totalgens <- 100000       # run for at most 100,000 generations; in practice this is set later: .simulate(numit)   [1000]  "the total number of generations to simulate"                           
numdemos <- 0             # only one within-population demography (life cycle) for all populations                 [0]     "the number of within population demographies"     
maxlandsize <- 1e+06      # overall landscape limit set to be bigger than populations expected to go               [2e+05] "the maximum number of individuals that can exist in the simulation"

#run the integer-populating function
rland <- landscape.new.intparam(rland, h=habitats, s=stages, cg=currentgen, ce=currentepoch, totgen=totalgens, maxland=maxlandsize)

# a few parameters don't have assigned variables to feed function but belong to this part of the landscape - easiest to fill in individually after.
rland$intparam$locusnum <- locusnum
rland$intparam$numepochs <- numepochs
rland$intparam$numdemos <- numdemos

#if desired, inspect results # print(rland$intparam)
cat(sprintf("The integer parameters have been set as follows:\n Number of habitats (populations) - %s\n Number of life stages - %s\n Number of loci - %s\n Number of epochs -  %s\n Current generation - %s\n Current epoch - %s\n Total generations - %s\n Number of within-population demographies - %s\n Maximum number of individuals in the simulation - %s\n ", 
            rland$intparam$habitats, rland$intparam$stages, rland$intparam$locusnum, rland$intparam$numepochs, rland$intparam$currentgen, rland$intparam$currentepoch, rland$intparam$totalgens, rland$intparam$numdemos, rland$intparam$maxlandsize))



# Landscape Building Part 2: Set on/off values for parameters that are switches

# function landscape.new.switchparam(rland, ...) sets boolean values for parameters that are switches. 4 values of 1 or 0, in order, set parameters.

randepoch <- 0       # only one epoch - no need for random transitions                                   [0] "1=choose the next epoch randomly using their individual probabilities, 0=choose the next epoch according to their ordering"     
randdemo <- 0        # only one local demography across the board - no need for random assignments       [0] "1=assign demographies at random, 0=assign demographies in order"                                                                
multp <- 1           # multiple paternity on, females can be pollinated by more than one male plant      [1] "1=multiple paternity, 0=entire families from a single mating"                                                                   
densdepdemo <- 0     # a single demography will apply regardless of density                              [0] "a flag to turn density-dependent population regulation on and off"                                                        

#run the switch parameter-populating function
rland <- landscape.new.switchparam(rland, randepoch, randdemo, multp, densdepdemo)

#if desired, inspect results #print(rland$switchparam)
cat(sprintf("the switches are set as follows:\n random epoch - %s\n random demography - %s\n multiple paternity - %s\n density dependent demographies - %s", rland$switchparam$randepoch, rland$switchparam$randdemo,rland$switchparam$multp, rland$switchparam$densdepdemo))




# Landscape Building Part 3: fill in parameters with float values

# there's only one parameter with a simple float value - the selfing rate: to set, rland <- landscape.new.floatparam(rland, float)
# We don't know selfing rate; Phacelia attracts pollinators (looking for outside pollen) 
#but doesn't have morphology to absolutely prevent selfing (more research required + run various  to assess impact)

selfing <- 0.5      # guess with some selfing, some outside pollen      [0] "This value is the selfing rate of the species (between 0 and 1 inclusive) assuming a mixed mating model where selfing occurs S proportion of the time and random mating occurs 1-S proportion of the time"

#execute command to set selfing rate
rland <- landscape.new.floatparam(rland, selfing)

# if desired inspect results #print(rland$floatparam)
print(sprintf("the selfing rate is now set to %s", rland$floatparam$selfing))


## Landscape Building Part 4: (a) local demographics (life stages) and (b) between-population demographics (migration)

# (a) The local demography - within-population life cycle dynamics
# Local demography is composed of survival, reproduction and male contribution matrices, which are combined to advance the population each generation
# There can be a different one for each subpopulation (chosen at random from the set or assigned in order)
# Local demography S, R, M matrices must be square n x n, where n is the number of stages in each population  

# Explanation of demography choice for a slightly complicated biennial
# stage 1 = seedlings, new rosettes
# stage 2 = juvenile rosettes created from last year's seedlings and waiting until the next spring to flower (an extra non reproductive year)
# stage 3 = flowering plants that reproduce then die, created from last years seedlings and year old juvenile rosettes.
# no stage is allowed to linger - individuals either move to another stage or die every year.
# values set to give a stable population that is 90% seedlings, 1% juveniles, 9% flowering plants

#Survival matrix 
#       from   
#       1  2  3
#to   1
#     2        
#     3        

localS <- matrix(c(            0,     0,   0,
                    0.0111111111,     0,   0,
                    0.0916666667,  0.75,   0) , byrow = TRUE, nrow = 3)
# reading down each column (transitions go from column to row)
# No seedlings remain seedlings, ~1% become juveniles, ~9% become flowering plants next year [~90% die]
# No juveniles become seedlings, no juveniles remain juveniles, 75% of juveniles survive to become flowering plants next year [25% die]
# No mature flowering plants become seedlings or juveniles and none survive another year [100% die]


#Reproduction - Maternal production of offspring
localR <- matrix(c(  0,  0,  10,
                     0,  0,  0,
                     0,  0,  0) , byrow = TRUE, nrow = 3)
print(localR)
#this means each flowering plant contributes 10 seedlings to the next generation, on average (drawn from poisson about 10)


# Male contribution
# Probability that pollen comes from a particular stage class.
# all pollen comes from flowering plants
# ****** Question ******** 
# should this value be 1 or 1-contents of migration matrices? how are these matrices combined? 100% percent from flowering plants, but e.g. 90% from within population, 10% from other populations.
localM <- matrix(c(  0,  0,  0,
                     0,  0,  0,
                     0,  0,  1) , byrow = TRUE, nrow = 3)

# set these matrices into the landscape
rland <- landscape.new.local.demo(rland,localS,localR,localM)

# repeating this command with more matrices puts more demographies into the localdem list (e.g. one for each population)

# see what this made
print(rland$demography$localdem)

# To enact density dependent local demography 
# 1) switch flag to on in switchparam
# 2) add local demographies in pairs, with flag k=0, k=1
# The demography for each population entered as above represents the low-density case (k=0 flag default) 
# You can enter a second demography to represent a different demography at carrying capacity (and these pairs can be different for each population)
# For a given population size, the rates will be linearly interpolated between the matrices. 
# To enter the carrying-capacity demographies, add flag k = 1 to local dem command: rland <- landscape.new.local.demo(rland,Sk,Rk,Mk, k=1)

# if desired, take a look at these, in their own list:
# print(rland$demography$localdemK)

# To change local demographies over time, I think you would run simulation for a time, change the landscape object at the demography matrices, then run it some more.
# to assign different local demographies to different populations, while adding demographies, keep in order of habitats as entered (1, 2, 3)
# to assign local demographies probabilistically, use randdemo switch and make several demographies (confusing, method for doing this comes later, localprob in new epoch.)
#***************Question***************** how do you assign different probabilities for local demographies if one value stands for epoch? could that be a vector? would make more sense to put probabilities with local demographies, not epochs.



# (b) The interpopulation movement of individuals 
# Migration matrices, each set of which is called an epoch 
# In one epoch, you'd have one kind of migration; mountain range appears, pollinator goes extinct, new kind of migration in the next epoch.
# Another set of S, R and M matrices describing movement between populations, with life-stage resolution (one set of three for each epoch)
# Each matrix is a square h x s on a side (# of habitats * number of stages) - a matrix of life stage matrices, for each habitat population
# For S and R matrices, values represent migration rates: fraction of FROM/source population that gets added to TO/target population (extrapolating from how local demography matrix multiplications work)
# For M matrices, values represent likelihood of male parentage of the TO/target population coming from the FROM/source population



#Survival: No juvenile or mature plants travel anywhere (unless a movement-prone clonal reproduction path exists?)
#        from  
#                  pop 1       pop 2       pop 3
#               juv  mat    juv  mat    juv  mat
# to   pop1 juv   0    0      0    0      0    0
#           mat   0    0      0    0      0    0 

#      pop2 juv   0    0      0    0      0    0
#           mat   0    0      0    0      0    0 

#      pop3 juv   0    0      0    0      0    0
#           mat   0    0      0    0      0    0 

#Reproduction: seeds might travel a little between populations (mature members of one population contribute juveniles to another population)
#        from  
#                  pop 1       pop 2       pop 3
#               juv  mat    juv  mat    juv  mat
# to   pop1 juv   0    0      0 0.01      0 0.01
#           mat   0    0      0    0      0    0 

#      pop2 juv   0 0.01      0    0      0 0.01
#           mat   0    0      0    0      0    0 

#      pop3 juv   0 0.01      0 0.01      0    0
#           mat   0    0      0    0      0    0 

# Male contribution: pollen can travel between populations (here 10% probability that pollen for a given seed will come from a particular population)
#        from  
#                  pop 1       pop 2       pop 3
#               juv  mat    juv  mat    juv  mat
# to   pop1 juv   0    0      0    0      0    0
#           mat   0    0      0  0.1      0  0.1

#      pop2 juv   0    0      0    0      0    0
#           mat   0  0.1      0    0      0  0.1 

#      pop3 juv   0    0      0    0      0    0
#           mat   0  0.1      0  0.1      0    0 

# These are going to be sparse (you could have a lot of migration, but probably only between a few lifestage & SRM forms)

# For these migration scenarios, sub matrices on the matrix-of-matrices diagonal (migration from a population to itself) are all zeros in examples,
# perhaps because shifts within populations are covered by the local demographies.


# Construct these matrices (m for migration)
# no plants travel (survival matrices all zeroes)
Sm <- matrix(c(rep(0,81)), nrow=9)

# values here based on replacement rates (gene flow without affecting size of population, each population gives and takes the same) for populations with sizes 5000,5000,18000
# Reproductive rates (number of matures * rate = number of next years seedlings) contributed by other populations
# ********Question********* Are these added to the within-population demographies? 
# or both added to and taken from? does this mean extra reproduction for migration?
# if extra, subtract from 10, giving each their own local demographies.
#Rm <- matrix(c(   0,    0,    0,     0,   0, 0.05,     0,  0, 0.0125,
#                  0,    0,    0,     0,   0,    0,     0,  0,      0,
#                  0,    0,    0,     0,   0,    0,     0,  0,      0,
#                  
#                  0,    0, 0.05,     0,   0,    0,     0,  0, 0.0125,
#                  0,    0,    0,     0,   0,    0,     0,  0,      0,
#                  0,    0,    0,     0,   0,    0,     0,  0,      0,
#                  
#                  0,    0, 0.05,     0,   0, 0.05,     0,  0,      0,
#                  0,    0,    0,     0,   0,    0,     0,  0,      0,
#                  0,    0,    0,     0,   0,    0,     0,  0,      0,
#                  ), byrow=TRUE, nrow=6)

#forget about seeds for now.
Rm <- matrix(c(rep(0,81)), nrow=9)

#5% probability that pollen comes from FROM population, when pollinating female parents in TO population
#********Question******** does this happen after local demography or should these rates be subtracted from 1 in lifestage within-population demography?
#********Question******** are diagonals supposed to be zero? how is this combined with local matrices.
Mm <- matrix(c(   0,  0,    0,     0,  0,    0,    0,  0,    0,
                  0,  0,    0,     0,  0,    0,    0,  0,    0,
                  0,  0,    0,     0,  0, 0.05,    0,  0, 0.05,
                  
                  0,  0,    0,     0,  0,    0,    0,  0,    0,
                  0,  0,    0,     0,  0,    0,    0,  0,    0,
                  0,  0, 0.05,     0,  0,    0,    0,  0, 0.05,
                  
                  0,  0,    0,     0,  0,    0,    0,  0,    0,
                  0,  0,    0,     0,  0,    0,    0,  0,    0,
                  0,  0, 0.05,     0,  0, 0.05,    0,  0,    0), byrow=TRUE, nrow=9)

#take a look
print(Sm)
print(Rm)
print(Mm)

#After the matrices, the epoch set-up command also accepts other arguments about demography

RndChooseProb <- 1            # random flag not turned on, unimportant, leave default           [1]                 "the probability of randomly choosing this epoch if rland$switchparam$randepoch==1"
StartGen <- 0                 # this epoch will staert at generation 0                          [0]                 "the generation this epoch begins if rland$switchparam$randepoch==0"
Extinct <- c(0,0,0)           # no extinction in model (could happen if populations shrink)     [c(0,0,0)]          "a vector of values giving the extinction probability for each subpopulation"
Carry <- c(1e+5,1e+5,1e+5)    # carrying capacities set to higher than expected so not limiting [c(1000,1000,1000)] "a vector of values denoting the maximum individuals for each subpopulation"
Localprob <- 1                # the probability of choosing the one local demography is 1       [1]                 "the probability of choosing this local demography if rland$switchparam$randdemo==1"

#populate an epoch with parameters

rland <- landscape.new.epoch(rland, S=Sm, R=Rm, M=Mm, epochprob = RndChooseProb, startgen = StartGen, extinct = Extinct, carry = Carry, localprob = Localprob)

#repeat command to have more epochs

#take a look
rland$demography$epochs


#5. Loci
#create lists of loci and alleles
#run default rland <- landscape.new.locus(rland)
#locus properties fed into landscape.new.locus
locustype <- 1          # strict stepwise model for mutation allele = integer, +-1 each mutation  [0? 251?] 0 = infinite alleles (alleles represented by integers, every mutation creates a new allele, no back mutation), 1 = strict step wise (alleles represented by integers, increase or decrease by 1 every mutation - like differences in length for simple sequence repeats - used for microsatellites with repeat regions, co-dominant, high mutation rates), type 2 = DNA base substitution (no mention of model yet, probably JC)
ploidy <- 2             # diploid                                                                 [1]       "This can be 1 for haploid and 2 for diploid"
mutationrate <- 0.01  # per-allele mutation rate 0.0001                                         [0]       "the per-allele mutation rate (not per site) range [0 1] (average # of mutations per allele per generation)"
transmission <- 0       # biparental inheritance                                                  [1]       "The mode of inheritance. This can be 0 for biparental or 1 for maternal."

# Each locus also has an object called alleles, which is a list of alleles for the locus, each of which is a list with more information     
# alleles defaults list of 2 alleles, aindex 1 and 2, birth=0, prop=0.5, state 1 and 2 respectively

# allele properties at this locus are also fed into landscape.new.locus to populate initial allele frequencies
# allelesize <- 50  # no allele size used for strict stepwise model [50] for type==2 DNA base substitution/variable length sequence state, enter custom DNA length here (one for locus, not per allele, so not quite what we want for microsatellite-like evolution)      

numalleles <- 100  # allow there to be 100 possible states of the microsatellite at this locus [2] number of alleles at this locus
frequencylist <- c(rep(0,100)) 
frequencylist[46] <- 0.125
frequencylist[47]<- 0.125
frequencylist[48]<- 0.125
frequencylist[49]<- 0.125
frequencylist[50] <- 0.125
frequencylist[51] <- 0.125
frequencylist[52] <- 0.125
frequencylist[53] <- 0.125
# only populate 8 alleles to start (then mutation free to explore the rest) [c(0.5,0.5)] Allele frequencies at locus (Same length as number of alleles)
# might be fun to make this more random.
print(frequencylist)

# statelist <- # the starting alleles themselves - no need to enter, they are generated by program (although you could) for type 0 or 1 this fills in integers [1:numalleles]; for type 2, generates DNA of length allelesize sampling ACTG from uniform distribution


# landscape.new.locus This function populates the loci parameters and makes an the allele list of length numalleles so that each allele gets the following properties:
# aindex index number of allele
# birth generation the allele first appeared
# prop - frequency of this allele (from frequency list)
# state - an integer for loci of type 0 or 1 (infinite allele or strict stepwise mutation), a string representing DNA sequence for loci type 2 (ACTG)
# possible entries: rland <- landscape.new.locus(rland, type=locustype, ploidy=ploidy, mutationrate=mutationrate, transmission=transmission, numalleles=numalleles, allelesize=allelesize,frequencies=frequencylist, states=statelist)


#If for some reason we wanted to get into allele construction, code to do so (from landscape.new.locus) would look like this
#rland$loci[[locusnum]]$alleles <- makealleles(type, numalleles, allelesize, frequencies, states)

#we'll probably want several loci every time.
numloci = 12

#implement locus builder function in a loop to make these 12 identical loci to start
for (l in c(1:numloci)){
  print(paste("building locus #",l))
  rland <- landscape.new.locus(rland, type=locustype, ploidy=ploidy, mutationrate=mutationrate, transmission=transmission, numalleles=numalleles, frequencies=frequencylist)
}


#see what's in there
print(rland$loci)

#6. Individuals
#a large matrix  with one individual per row
#structure of matrix - each individual gets these in their row
#stage notused birthdate ID matID patID loc1a1 loc1a2 loc2a1 loc2a2 loc3a1
#The first six columns are always present and the last five shown here are a product of the choice of number of loci. 
#The first column contains the individuals subpopulation and lifecycle stage (subpopulation = floor(x/rland\$intparam\$stages), lifecycle stage = x mod rland\$intparam\$stages). 
#The second column is currently unused, always 0. 
#The third column contains the generation in which the individual was born or created. 
#The next three contain numerical ids for the individual, its mother and its father. 
#After the first six columns the indivduals genetic code begins. 
#The loci are shown in order with 2 columns for diploid loci and 1 column for haploid loci. 
#The value of these columns represent the allele index of the allele the individual carries.

#To create individuals, enter a vector of population sizes
pop1seedling <- 900
pop1juv      <- 10
pop1mat      <- 90

pop2seedling <- 900
pop2juv      <- 10
pop2mat      <- 90

pop3seedling <- 900
pop3juv      <- 10
pop3mat      <- 90


initial_pop_structure = c(pop1seedling, pop1juv, pop1mat, pop2seedling, pop2juv, pop2mat, pop3seedling, pop3juv, pop3mat)

#these 9 categores become known as subpopulations 0-8 in dataset

rland <- landscape.new.individuals(rland,initial_pop_structure)

print(rland$individuals)
unique(rland$individuals[,1])
table(rland$individuals[,1])


# structure of individuals table (heart of data)
# "The first six columns are always present and the last n shown here are a product of the choice of number of loci. 
# The first column contains the individuals subpopulation and lifecycle stage:
# (subpopulation = floor(x/rland\$intparam\$stages), lifecycle stage = x mod rland\$intparam\$stages). 
# The second column is currently unused, always 0. 
# The third column contains the generation in which the individual was born or created. 
# The next three contain numerical ids for the individual, its mother and its father. 
# After the first six columns the indivduals genetic code begins. The loci are shown in order with 2 columns for diploid loci and 1 column for haploid loci. The value of these columns represent the allele index of the allele the individual carries."


#the allelic composition of individuals is chosen based on the allele frequencies (which are single per-locus values for all populations to start)
#****************Question**************is there a way to start with population structure? (different allele frequencies for different populations)

#simulate! 

#retval <- matrix(0,ncol=2,nrow=101) #store the results 1 = generation 2 = number of total individuals
#print(retval)
#retval[1,] <- c(0,nrow(rland$individuals)) #before any evolution occurs. Both populations from same source

for (gen in 1:100)
{
  print("generation")
  print(gen)
  rland <- landscape.extinct(rland)
  print("After extinct")
  print("Remaining subpopulations")
  print(unique(rland$individuals[,1]))
  print("Total individuals")
  print(nrow(rland$individuals))
  print("Number of individuals in each subpopulation")
  print(table(rland$individuals[,1]))
  print("")
  
  rland <- landscape.survive(rland)
  print("After survive")
  print("Remaining subpopulations")
  print(unique(rland$individuals[,1]))
  print("Total individuals")
  print(nrow(rland$individuals))
  print("Number of individuals in each subpopulation")
  print(table(rland$individuals[,1]))
  print("")
  
  
  
  rland <- landscape.reproduce(rland)
  print("After reproduce")
  print("Remaining subpopulations")
  print(unique(rland$individuals[,1]))
  print("Total individuals")
  print(nrow(rland$individuals))
  print("Number of individuals in each subpopulation")
  print(table(rland$individuals[,1]))
  print("")
  
  
  rland <- landscape.carry(rland)
  print("After carry")
  print("Remaining subpopulations")
  print(unique(rland$individuals[,1]))
  print("Total individuals")
  print(nrow(rland$individuals))
  print("Number of individuals in each subpopulation")
  print(table(rland$individuals[,1]))
  print("")
  
  
  rland <- landscape.advance(rland)
  print("After advance")
  print("Remaining subpopulations")
  print(unique(rland$individuals[,1]))
  print("Total individuals")
  print(nrow(rland$individuals))
  print("Number of individuals in each subpopulation")
  print(table(rland$individuals[,1]))
  print("")
  
  print("allele counts")
  print(table(rland$individuals[,7]))
  print(table(rland$individuals[,12]))
  print("")
}







table(rland$individuals[,7])


for (i in 2:30)#how many decades
{
  rland <- landscape.simulate(rland,1)
#  retval[i,] <- c(rland$intparam$currentgen, nrow(rland$individuals)) 
  print(sprintf("finished simulating decade %s of 100, generation %s", i-1,rland$intparam$currentgen ))
  print(unique(rland$individuals[,1]))
  print(nrow(rland$individuals))
  print(table(rland$individuals[,1]))
}






















#print(retval)


rland <- landscape.simulate(rland,1)
unique(rland$individuals[,1])

retval2 <- matrix(0,ncol=2,nrow=101) #store the results 1 = generation 2 = number of total individuals
#print(retval)
retval2[1,] <- c(0,nrow(rland$individuals)) #before any evolution occurs. Both populations from same source
for (i in 2:101)
{
  rland <- landscape.simulate(rland,10)
  retval2[i,] <- c(rland$intparam$currentgen, nrow(rland$individuals)) 
  print(sprintf("finished simulating decade %s of 100, generation %s", i-1,rland$intparam$currentgen ))
}
print(retval2)





#numyears = 1000
#rland <- landscape.simulate(rland,numyears)
# once these are time sampled, it would be nice to show compute time (Rtools for that?)

landscape.amova(rland)

rland.df <- rland$individuals
library(adegenet)

rland.df <- as.data.frame(rland.df)

pops <- seq(0,4,2)

alleles <- seq(7,29,2)
rland.2genind <- sapply(alleles, function(x){
  paste(rland.df[,x], rland.df[,x+1], sep = "-")
})


rland.df$pop[rland.df$V1 == 0 | rland.df$V1 == 1 ] <- "pop1"
rland.df$pop[rland.df$V1 == 2 | rland.df$V1 == 3 ] <- "pop2"
rland.df$pop[rland.df$V1 == 4 | rland.df$V1 == 5 ] <- "pop3"

# following doesn't work yet!
rland.df$pop <- sapply(pops, function(x){
  pop <- if(rland.df$V1 == x | rland.df$V1 == x+1){
    paste("pop",x,sep="")
  }
  pop
})

str(rland.df)

rland.df2genind <- data.frame(rland.df[,c(3,31)], rland.2genind)

library(strataG)

head(rland.df2genind)

rland.genind <- df2genind(rland.df2genind[,-(1:2)], pop = rland.df2genind$pop, type = "codom", sep="-")
class(rland.genind)
rland.genind

rland.gtypes <- genind2gtypes(rland.genind)

overallTest(rland.gtypes)
pairwiseTest(rland.gtypes)


"""
instead of simulate, run steps separately and save individuals and loci matrix at each step (or subset of steps)
rland <- landscape.new.example()
for (gen in 1:10)
{
  rland <- landscape.extinct(rland)
  rland <- landscape.reproduce(rland)
  rland <- landscape.survive(rland)
  rland <- landscape.carry(rland)
  rland <- landscape.advance(rland)
}
landscape.amova(rland)

another way is to simulate, save, then simulate some more.
rland <- create.land()
retval <- matrix(0,ncol=4,nrow=21) #store the results col1 = gen, cols2-3 PhiST
retval[1,] <- c(0,landscape.amova(rland)) #before any evolution occurs. Both populations from same source
for (i in 2:21)
    {
        rland <- landscape.simulate(rland,5)
        retval[i,] <- c(rland$intparam$currentgen,landscape.amova(rland)) 
    }
"""
