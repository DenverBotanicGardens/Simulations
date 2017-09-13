
#Warming up with rmetasim

#Jessie Berta-Thompson 
#Wed Dec 7, 2016

#workflow primarily from vignette, with some slight variations to fit project-relevant information in:
#https://github.com/stranda/rmetasim/blob/master/vignettes/CreatingLandscapes.Rmd

#ideally this script will set up a workflow for running rmetasim in the future.


#if Global Environment auto-loads materials from skelesim, conflicts can occur with 2 versions of rmetasim; clear environment before loading rmetasim with:
rm(list =ls())

#load rmetasim (if installed)
library(rmetasim)

#Parameters for Rmetasim simulations are stored in a landscape object (organized list) with 6 pieces that we'll fill in order (order matters)

#start by making an empty landscape
rland <- landscape.new.empty()



#1. Fill in parameters with integer values. 

#There is a function to fill these in - landscape.new.intparam(); a few parameters don't have assigned variables to feed function, had trouble changing them through ordered entry - easiest to fill in individually after.
#for readability & repeat use, I'm using long form variable names (same as within landscape) to feed into landscape.new.intparam()
#landscape.new.intparam defaults included in comments in []; pure default run just  landscape name; rland <- landscape.new.intparam(rland)
#comments start with vignette definitions, then in () my notes

    habitats <- 3         # [1]     the number of different habitats or subpopulations within the landscape        (how many populations to model - starting with 3 sets )
      stages <- 2         # [1]     the number of life stages in the transition matrices                           (often: juvenile, adult)
    locusnum <- 12        # [0]     the number of different loci currently implemented for the simulation          (e.g. number of microsatellite markers; 12 for phacelia project)
   numepochs <- 0         # [0]     the number of different migration scenarios to occur during the simulation     (could have one regime/set of migration matrices for a time, then another)
  currentgen <- 0         # [0]     the current generation the simulation has reached                              (should start at 0)
currentepoch <- 0         # [0]     the current epoch the simulation has reached                                   (should start at 0)
   totalgens <- 1000      # [1000]  the total number of generations to simulate                                    (important to run long enough that population properties reach equilibrium (if they can) and burnin past effects of arbitrary initial conditions) 
    numdemos <- 0         # [0]     the number of within population demographies                                   (probably one set of life cycle demography matrices will work for all of our populations over time)
 maxlandsize <- 1e+04     # [2e+05] the maximum number of individuals that can exist in the simulation             (carrying capacity across all the populations)

#run the integer-populating function (those values + the name of the lanscape)
rland <- landscape.new.intparam(rland, h=habitats, s=stages, cg=currentgen, ce=currentepoch, totgen=totalgens, maxland=maxlandsize)

#fill in the stragglers as needed
rland$intparam$locusnum <- locusnum
rland$intparam$numepochs <- numepochs
rland$intparam$numdemos <- numdemos

#make sure it looks ok
print(rland$intparam)


#2. Fill in parameters that are switches

# with function landscape.new.switchparam(rland, ...)
#rland <- landscape.new.switchparam(rland) fills in defaults; 4 values of 1 or 0, in order, following rland, set parameters

#set up variables with values, again, longform for readability, happen to match internal variable names; in comments # [default] vignette definition (my notes)

  randepoch <- 0     # [0] 1=choose the next epoch randomly using their individual probabilities, 0=choose the next epoch according to their ordering     (we'll probably only have one migration scheme over time)
   randdemo <- 0     # [0] 1=assign demographies at random, 0=assign demographies in order                                                                (we'll probably only have one demography - not important)
      multp <- 1     # [1] 1=multiple paternity, 0=entire families from a single mating                                                                   (for pollen-based mating, multiple appropriate)
densdepdemo <- 0     # [0] a flag to turn density-dependent population regulation on and off                                                              (within-populaton demography regulated by population density;to use this, later on you'll define 2 demographies, one for when the population is 0 (low) and one for when it is at carrying capacity - no doubt there's some probability distribution to pick between them for intermediate values 


#run the switch parameter-populating command
rland <- landscape.new.switchparam(rland, randepoch, randdemo, multp, densdepdemo)

#make sure it looks ok
print(rland$switchparam)


#3. Fill in parameters with float values
#rland <- landscape.new.floatparam(rland, float) allows you to set value of selfing rate (there's only one parameter built in here), rland <- landscape.new.floatparam(rland) sets up defaults

selfing <- 0.5 # [0] This value is the selfing rate of the species (between 0 and 1 inclusive) assuming a mixed mating model where selfing occurs S proportion of the time and random mating occurs 1-S proportion of the time (need to do some more lit review - known for some phacelia, not necessarily ours)



#4. Fill in demography information (in a few parts)

# A. The local demography - within-population life cycle dynamics
#survival, reproduction and male contribution matrices; can have one for each subpopulation (could get more complicated, e.g. time)
# these matrices must be square with the number of stages in each population to a side (they describe transitions between and within life stages)
#imagining stage 1 = juvenile, stage 2 = adult, these matrices would be something like:

#Survival
#        from   
#        juv  mat
#to  juv 0.4    0
#    mat 0.5    0

#this means
#40% of juveniles remain juveniles (each year)
#50% move up to become mature (values chose because biennial takes 2 years to mature, so ~50% of juveniles will mature each year depending on if it's their year or not, but not all will survive in a given year (this parameter set puts total survival at 90% - unaccounted for 20% die))
#0% of matures become juvenile (they do reproduce, but that's in the next matrix)
#0% of matures survive (biennial dies after reproducing)


#Reproduction
#        from   
#        juv  mat
#to  juv   0    5
#    mat   0    0

#this means each mature individual contributes 5 juveniles to the next generation, on average (in fact drawn from poisson about 5)
#no juveniles make juveniles, no matures make matures, no juveniles make matures (they do become mature, but that's in the other matrix)


#Male contribution
#        from   
#        juv  mat
#to  juv   0    0
#    mat   0    1

# I don't know what this means or how it gets combined with others (need to figure it out), but examples (vignette + published) seem to put a one in the lower right
# l in variable names for life cycle
# now implement these ideas
Sl <- matrix(c(0.4, 0, 0.5, 0), byrow = TRUE, nrow = 2)
Rl <- matrix(c(0, 5, 0, 0), byrow = TRUE, nrow = 2)
Ml <- matrix(c(0, 0, 0, 1), byrow = TRUE, nrow = 2)

#take a look at the matrices
print(Sl)
print(Rl)
print(Ml)

# set these matrices into the landscape
rland <- landscape.new.local.demo(rland,Sl,Rl,Ml)

# repeating this command can put more demographies in (building up a list of local demographies for each population, each item of which is a list of these three kinds of matrices)

# see what this made
print(rland$demography$localdem)

#Program can make local demographies density-dependent (if switch flag turned on and extra local demographies added)
#The demography for each population entered above represents the low-density case; You can enter a second demography for each population to represent demography at carrying capacity.
#For a given population size, the rates will be linearly interpolated between the matrices. 
#to set these carrying-capacity demographies, just add flag k = 1 to above command (defaults to 0)

#e.g.
#rland <- landscape.new.local.demo(rland,Sk,Rk,Mk, k=1l)

#keep adding more in order if different ones required for different populations
#see what you got with
#print(rland$demography$localdemK)


# B. The interpopulation movement of individuals 
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
        Carry <- c(1000,1000,1000)   # [c(1000,1000,1000)] a vector of values denoting the maximum individuals for each subpopulation
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

pop1juv <- 1
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

