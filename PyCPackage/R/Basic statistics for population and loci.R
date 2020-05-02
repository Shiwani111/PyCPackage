# Pyrus calleryana (callery pear) trees were imported from China to the US to combat fireblight disease in European pear in 1917. Those callery pear trees were hugely used as rootstock for several edible pears. Now, these trees are invasive species in most of the eastern and southern states of the US. For our study, P. calleryana samples were collected from China, Japan, Korea and US. Initially, we had 100 samples but at the end, we were limited with only 57 samples as other samples were cutoff due to several reasons. Our study aims at studying the genetic diversity of those trees and evaluate the relatedness of US species with the Asian ones. Microsatellite loci are used for the study. 
library(poppr)
library(hierfstat)
library(mmod)
library(ape)
library(magrittr)
library(adegenet)
library(MASS)
library(geiger)
library(pegas)
library(treemap)
library(car)
library(agricolae)
library(ggplot2)
library(foreign)
library(plyr)
library(caTools)
library(reshape2)
library(ecodist)
library(polysat)
library(radiator)
library(diveRsity)
library(doParallel)
library(PopGenReport)
library(devtools)
library(roxygen2)
library(pegas)
library(lattice)
library(ade4)
library(vegan)

#Loading Genalex datafile
#' Loading Genalex datafile
#' 
#' Generates file with details about individuals, loci and allele sizes
#' @export
tks <- read.genalex("C:/Users/ssapkot1/Desktop/GenalexAsian.csv", ploidy = 2, geo = FALSE, region = FALSE, genclone = FALSE)
tks

# Doing Clone correction
#' Doing Clone correction
#' 
#' Checks if there are any clones and if there are any clones then corrects them
#' @export
missingno(tks, type = "loci", cutoff = 0.99, quiet = FALSE,
          freq = FALSE)
setPop(tks) <- ~Pop 
tks<-clonecorrect(tks, strata = ~Pop, combine = FALSE, keep = 1)
tks

# Displaying the class of datafile
#' Displaying the class of datafile
#' 
#' @export
class(tks)

# Displaying information about loci and populations
#' Displaying information about loci and populations
#' 
#' Displays the table with the details of genetic diversity parameters of loci and populations
#'  @export
info_table(tks, plot = TRUE)

# Basic statistics
#' Displaying basic statistics
#' 
#' Summarizes the file with genetic diversity parameters
#' @export
summary(tks)

# Getting per locus information
#' Displaying per locus information
#' 
#' Summarizes the locus information with its genetic diversity details
#' @export
BS <- basic.stats(tks)
BS

# Testing Hardy-Weinberg Equillibrium
#' Testing for Hardy-Weinberg Equillibrium
#' 
#' Looks if the loci follow HWE with alpha 0.05
#' @export
hw.test(tks, B = 1000) #Rough test. Pr.exact = 0 -> Locus not in HWE.
(tks.HWE <- hw.test(tks, B = 1000)) # performs 1000 permuatations
(tks.HWE.pop <- seppop(tks) %>% lapply(hw.test, B = 1000)) #by population
(tks.HWE.mat <- sapply(tks.HWE.pop, "[", i = TRUE, j = 3)) # Take the third column with all rows
alpha  <- 0.05
tksHWEmat <- tks.HWE.mat
tksHWEmat[tksHWEmat > alpha] <- 1
levelplot(t(tksHWEmat), scales=list(x=list(rot=43)))

#Looking how uninformatic loci are
#' Showing how uninformatic the used loci are
#' 
#' Displays how uninformatic loci are with cutoff value
#' @export
nLoc(tks)
i_tks <- informloci(tks)

#Multi Locus Genotype accumulation curve
#' Generating multi locus genotype accumulation curve
#' 
#' Checks how much power the used loci have to identify unique MLGs
#' @export
gac <- genotype_curve(tks, sample = 1000, quiet = TRUE)
gac
P.tab<-mlg.table(tks)

#Linkage disequilibrium (LD)
#The index of association
#' Displaying Index of Association
#' 
#' Checks how the populations are associated using IoA
#' @export
tks2 <- popsub(tks)
set.seed(1044)
ia(tks2, hist = TRUE, sample = 999)

#Pairwise LD over all loci
#' Displaying Linkage Disequillibrium over all loci
#' 
#' Summarizes if the used loci are linked or not
#' @export
tkspair <- tks %>% clonecorrect(strata = ~Pop) %>% pair.ia
plot(tkspair, label=TRUE, low = "orange", high = "purple")

#Minimum spanning network
#' Displaying minimum spanning network
#' 
#' Summarizes how the populations used are spanned in a network
#' @export
tks_sub <- popsub(tks, blacklist = character(0))
tks_dist <- diss.dist(tks_sub, percent = FALSE, mat = FALSE)
min_span_net <- poppr.msn(tks_sub, tks_dist, showplot = FALSE, include.ties = TRUE)
set.seed(69)
plot_poppr_msn(tks, min_span_net, inds = "ALL", mlg = FALSE, gadj = 6, nodescale = 10,
               palette = seasun, cutoff = NULL, quantiles = FALSE, beforecut = TRUE, pop.leg = TRUE,
               size.leg = FALSE, scale.leg = TRUE, layfun = igraph::layout_nicely)

#AMOVA
#' Displaying AMOVA 
#' 
#' Summarizes if there exists some genetic structure in populations or not
#' @export
tks <- as.genclone(tks)
tks
table (strata(tks, ~Pop))
tksamova <- poppr.amova(tks, ~Pop)
tksamova

#Significance testing
#' Displaying significance test for AMOVA
#' 
#' Significance test is done for AMOVA using 999 number of repeats
#' @export
set.seed(1999)
tkssignif <- randtest(tksamova, nrepet = 999)
plot(tkssignif)
tkssignif

#Discriminant analysis of principal components (DAPC)
#' Displaying DAPC to categorize the populations in different clusters 
#' 
#' Gives information about how the populations are related and clustered across two principal components (PCs)
#' @export
l.col <- seasun(40)
pop(tks) <- tks$pop
dapc.tks <- dapc(tks, var.contrib = TRUE, scale = FALSE, n.pca = 19, n.da = nPop(tks) - 1)
scatter(dapc.tks, col=l.col, cell = 0, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)

#DAPC cross-valiation
set.seed(1044)
pramx <- xvalDapc(tab(tks, NA.method = "mean"), pop(tks))
system.time(pramx <- xvalDapc(tab(tks, NA.method = "mean"), pop(tks),
                              n.pca = 5:25, n.rep = 10,
                              parallel = "multicore", ncpus = 4L))
names(pramx)
l.col <- seasun(40)
pramx[-1]
max(pramx$`Mean Successful Assignment by Number of PCs of PCA`)
contrib <- loadingplot(dapc.tks$var.contr, axis = 2, thres = 0.055, lab.jitter = 1)
contrib
scatter(dapc.tks, cex = 2, label.inds = FALSE, col = l.col, legend = TRUE, clabel = TRUE, posi.leg = "bottomleft", scree.pca = FALSE, scree.da = TRUE, posi.pca = "topleft" , posi.da = "bottomright", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1, cstar=0)
scatter(dapc.tks, cex = 2, label.inds = NULL, col = l.col, legend = FALSE, clabel = TRUE, posi.leg = "topleft", scree.pca = FALSE, scree.da = TRUE, posi.pca = "topright" , posi.da = "bottomright", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1, cstar=0)