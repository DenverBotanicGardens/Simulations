
#Building an rmetasim script to explore initial parameter choices for Phacelia populations

# Jessie Berta-Thompson & Michelle DePrenger-Levin
# Wed Dec 14, 2016

# Source materials (occasionally quoted)
# rmetasim code and documentation, vignettes, Metasim publication Strand 2002, Kernelpop publication and documentation, Strand and Niehaus 2007



# This script contains a workflow for running rmetasim for
# 3 populations of Phacelia with light migration between them, with rough informed guesses for parameters to be explored in more detail later.


# To start with a clean R slate, clear environment:
rm(list =ls())

# load rmetasim (if installed)
library(rmetasim)


# Create an rmetasim landscape
# Parameters for rmetasim simulations are stored in a landscape object (an organized list of lists and variables)
# We will fill in 6 pieces in order (earlier ones influence later ones)

# "The six main sections of the landscape:
# intparam  - integer parameters describing landscape
# switchparam -  boolean switches turning features on and off
# floatparam - floating point paramters
# demography - lists describing the demography of the simulation (complex set of parameters)
# loci lists - describing the genetic components of the simulation
# individuals - a matrix that contains all individuals, their demographic states, ages and ids. It also contains the genetic information.
# The last three are values can change through the course of the simulation. Especially the last two"


#start by making an empty landscape
rland <- landscape.new.empty()

# for readability & repeat use, using long form variable names to set parameters, then pass to rmetasim functions.
# comment structure: explain parameter choice for this run, defaults in [], "vignette/publication parameter definitions", additional notes in ()


# Landscape Building Part 1. Fill in parameters with integer values. 
# The function landscape.new.intparam(rland, ...) sets the values of a set of parameters that are coded with integers

habitats <- 3             # 3 Phacelia populations                                                                 [1]     "the number of different habitats or subpopulations within the landscape"                                             
stages <- 3               # 3 life stages: seedlings, still juvenile year-old rosettes, flowering                  [1]     "the number of stages in the life cycle of the organism; also the size of the S, R, and M matricies for local demography
locusnum <- 12            # 12 microsatellite loci used in genotyping study                                        [0]     "the number of different loci currently implemented for the simulation"
numepochs <- 0            # one migration scenario over time                                                       [0]     "the number of different migration scenerios to occur during the simulation"
currentgen <- 0           # start at generation 0                                                                  [0]     "the current generation the simulation has reached"                              
currentepoch <- 0         # start at epoch 0                                                                       [0]     "the current epoch the simulation has reached"                                   
totalgens <- 100000       # run for at most 100,000 generations; in practice this is set later: .simulate(numit)   [1000]  "the total number of generations to simulate"                           
numdemos <- 0             # only one within-population demography for all                                          [0]     "the number of within population demographies"     
maxlandsize <- 1e+06      # overall landscape limit set to be bigger than populations expected to go               [2e+05] "the maximum number of individuals that can exist in the simulation"

#run the integer-populating function
rland <- landscape.new.intparam(rland, h=habitats, s=stages, cg=currentgen, ce=currentepoch, totgen=totalgens, maxland=maxlandsize)

# a few parameters don't have assigned variables to feed function but belong to this part of the landscape - easiest to fill in individually after.
rland$intparam$locusnum <- locusnum
rland$intparam$numepochs <- numepochs
rland$intparam$numdemos <- numdemos

#if desired, inspect results
print(rland$intparam)




# Landscape Building Part 2: Set on/off values for parameters that are switches

# function landscape.new.switchparam(rland, ...) sets boolean values for parameters that are switches
# 4 values of 1 or 0, in order, set parameters

randepoch <- 0       # only one epoch - no need for random transitions                                   [0] "1=choose the next epoch randomly using their individual probabilities, 0=choose the next epoch according to their ordering"     
randdemo <- 0        # only one local demography across the board                                        [0] "1=assign demographies at random, 0=assign demographies in order"                                                                
multp <- 1           # multiple paternity on, females can be pollinated by more than one male plant      [1] "1=multiple paternity, 0=entire families from a single mating"                                                                   
densdepdemo <- 0     # a single demography will apply regardless of density                              [0] "a flag to turn density-dependent population regulation on and off"                                                        

#run the switch parameter-populating function
rland <- landscape.new.switchparam(rland, randepoch, randdemo, multp, densdepdemo)

#if desired, inspect results
print(rland$switchparam)




# Landscape Building Part 3: fill in parameters with float values

# there's only one parameter with a float value built in - the selfing rate: to set, rland <- landscape.new.floatparam(rland, float)

selfing <- 0.5 # We don't know selfing rate; attracts pollinators but doesn't have morphology to prevent selfing (more research required + run various  to assess impact)     [0] "This value is the selfing rate of the species (between 0 and 1 inclusive) assuming a mixed mating model where selfing occurs S proportion of the time and random mating occurs 1-S proportion of the time"




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

localS <- matrix(c(            0,     0,   0
                    0.0111111111,     0,   0
                    0.0916666667,  0.75,   0) , byrow = TRUE, nrow = 3)

# No seedlings remain seedlings, 1% become juveniles, 9% become flowering plants [90% die]
# No juveniles become seedlings, no juveniles remain juveniles, 75% of juveniles survive to become flowering plants [25% die]
# No mature flowering plants become seedlings or juveniles and none survive another year [100% die]


#Reproduction - Maternal production of offspring
localR <- matrix(c(  0,  0, 10
                     0,  0,  0
                     0,  0,  0) , byrow = TRUE, nrow = 3)

#this means each flowering plant contributes 10 seedlings to the next generation, on average (drawn from poisson about 10)


# Male contribution
# Probability that pollen comes from a particular stage class.
localR <- matrix(c(  0,  0,  0
                     0,  0,  0
                     0,  0,  1) , byrow = TRUE, nrow = 3)

# set these matrices into the landscape
rland <- landscape.new.local.demo(rland,localS,localR,localM)

# repeating this command with more matrices puts more demographies into the localdem list (e.g. one for each population)

# see what this made
print(rland$demography$localdem)

# Density Dependent local demography 1) switch flag turned on in switchparam and 2) local demographies added in pairs, with flag k=0, k=1
# The demography for each population entered above represents the low-density case (k=0 flag default)
# You can enter a second demography to represent demography at carrying capacity (and these pairs can be different for each population)
# For a given population size, the rates will be linearly interpolated between the matrices. 
# To enter the carrying-capacity demographies, add flag k = 1 to local dem command: rland <- landscape.new.local.demo(rland,Sk,Rk,Mk, k=1)

# if desired, take a look at these, in their own list:
# print(rland$demography$localdemK)

# To change local demographies over time, I think you would run simulation for a time, change the landscape object, then run it some more.
# to assign local demographies to populations in order, keep in order
# to assign local demographies probabilistically, *******





# (b) The interpopulation movement of individuals 
# Migration matrices, each set of which is called an epoch (in one epoch, you'd have one kind of migration; mountain range appears, pollinator goes extinct, new kind of migration in the next epoch) 
#Another set of S, R and M matrices describing movement between populations, with life-stage resolution.
# Each matrix is a square h x s on a side (# of habitats * number of stages) - a matrix of life stage matrices, for each habitat population

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

# Male contribution: pollen can travel between populations (maybe matures contributing to gametes which contribute to reproduction?)
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
# For Survival matrix this means individuals travelling, for Reproduction offpsring travelling, for Male contribution pollen travelling.
# Need to figure out how to make this accurately reflect travelling pollen (above just a guess based on form of examples) - same issue as above with local demography matrices - I'm not totally grasping M matrices.
# For these migration scenarios, sub matrices on the matrix-of-matrices diagonal (migration from a population to itself) are all zeros in given examples- shifts within populations are covered by the other set of matrices
# These values represent migration rates: fraction of FROM/source population that gets added to TO/target population
# likelihood of male parentage from population 1, 2 or 3

#Construct these matrices (m for migration)

Sm <- matrix(c(rep(0,36)), nrow=6)

Rm <- matrix(c(   0,   0,   0,0.01,   0,0.01,
                  0,   0,   0,   0,   0,   0,
                  0,0.01,   0,   0,   0,0.01,
                  0,   0,   0,   0,   0,   0,
                  0,0.01,   0,0.01,   0,   0,
                  0,   0,   0,   0,   0,   0), byrow=TRUE, nrow=6)

Mm <- matrix(c(   0,   0,   0, 0.1,   0, 0.1,
                  0,   0,   0,   0,   0,   0,
                  0, 0.1,   0,   0,   0, 0.1,
                  0,   0,   0,   0,   0,   0,
                  0, 0.1,   0, 0.1,   0,   0,
                  0,   0,   0,   0,   0,   0), byrow=TRUE, nrow=6)

#take a look
print(Sm)
print(Rm)
print(Mm)

#After the matrices, the epoch set-up command also accepts other arguments about demography

RndChooseProb <- 1                   # [1]                 the probability of randomly choosing this epoch if rland$switchparam$randepoch==1
StartGen <- 0                   # [0]                 the generation this epoch begins if rland$switchparam$randepoch==0 
Extinct <- c(0,0,0)            # [c(0,0,0)]          a vector of values giving the extinction probability for each subpopulation
Carry <- c(1e+5,1e+5,1e+5)   # [c(1000,1000,1000)] a vector of values denoting the maximum individuals for each subpopulation
Localprob <- 1                   # [1]                 the probability of choosing this local demography if rland$switchparam$randdemo==1

#populate an epoch with parameters

rland <- landscape.new.epoch(rland, S=Sm, R=Rm, M=Mm, epochprob = RndChooseProb, startgen = StartGen, extinct = Extinct, carry = Carry, localprob = Localprob)

#repeat command to have more epochs

#take a look
rland$demography$epochs


#5. Loci
#create lists of loci and alleles
#run default rland <- landscape.new.locus(rland)
#locus properties fed into landscape.new.locus
locustype <- 1       # [0? 251?] 0 = infinite alleles (alleles represented by integers, every mutation creates a new allele, no back mutation), 1 = strict step wise (alleles represented by integers, increase or decrease by 1 every mutation - like differences in length for simple sequence repeats - used for microsatellites with repeat regions, co-dominant, high mutation rates), type 2 = DNA base substitution (no mention of model yet, probably JC)
ploidy <- 2       # [1] This can be 1 for haploid and 2 for diploid
mutationrate <- 0.0001  # [0] the per-allele mutation rate (not per site) range [0 1] (average # of mutations per allele per generation)
transmission <- 0       # [1] The mode of inheritance. This can be 0 for biparental or 1 for maternal.

#Each locus also has an object called alleles, which is a list of alleles for the locus, each of which is a list with more information     
#alleles defaults list of 2 alleles, aindex 1 and 2, birth=0, prop=0.5, state 1 and 2 respectively]

#allele properties at this locus are also fed into landscape.new.locus to populate initial allele frequencies
#allelesize <- 50  # [50] for type==2 DNA base substitution/variable length sequence state, enter custom DNA length here (one for locus, not per allele, so not quite what we want for microsatellite-like evolution)      
numalleles <- 8 # [2] number of alleles at this locus
frequencylist <- c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125) # [c(0.5,0.5)] Allele frequencies at locus (Same length as number of alleles)
#statelist <- # the starting alleles themselves - no need to enter, they are generated by program (although you could) for type 0 or 1 this fills in integers [1:numalleles]; for type 2, generates DNA of length allelesize sampling ACTG from uniform distribution

#all possible entries are here: rland <- landscape.new.locus(rland, type=locustype, ploidy=ploidy, mutationrate=mutationrate, transmission=transmission, numalleles=numalleles, allelesize=allelesize,frequencies=frequencylist, states=statelist)

#This function populates the loci parameters and makes an the allele list of length numalleles so that each allele gets the following properties:
# aindex index number of allele
# birth generation the allele first appeared
# prop - frequency of this allele (from frequency list)
# state - an integer for loci of type 0 or 1 (infinite allele or strict stepwise mutation), a string representing DNA sequence for loci type 2 (ACTG)


#If for some reason we wanted to get into allele construction, code to do so (from landscape.new.locus) would look like this
#rland$loci[[locusnum]]$alleles <- makealleles(type, numalleles, allelesize, frequencies, states)

#we'll probably want several loci every time.
numloci = 12

#implement locus builder function in a loop to make these 12 identical loci
for (l in c(1:numloci)){
  print(paste("building locus #",l))
  rland <- landscape.new.locus(rland, type=locustype, ploidy=ploidy, mutationrate=mutationrate, transmission=transmission, numalleles=numalleles, frequencies=frequencylist)
}


#see what's in there
#print(rland$loci)

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

pop1juv <- 3000
pop1mat <- 2

pop2juv <- 3
pop2mat <- 4

pop3juv <- 5
pop3mat <- 6


initial_pop_structure = c(pop1juv, pop1mat, pop2juv, pop2mat, pop3juv, pop3mat)

#these 6 categores become known as subpopulations 0-5 in dataset (I'm guessing on order/structuring here). 

rland <- landscape.new.individuals(rland,initial_pop_structure)


print(rland$individuals)

#simulate! 
numyears = 10
rland <- landscape.simulate(rland,numyears)
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
