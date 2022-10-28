### NS-WL Linkages MS
# KE O'Reilly et al.
# Code to run MixSIAR Bayesian mixing models

# Mixing models were run separately at each of the 7 sites
# and run separtely for the following Groups:
# NS Large Perch, WL Large Perch, NS Small Perch, WL Small Perch

#Load packages
library(tidyverse)
library(FSA)
library(ggrepel)
library(gridExtra)
library(ggpubr)
library(ggsignif)
library(MixSIAR)
library(mcmc)
library(sp)
library(splancs)

############Burns Harbor###################
### South Region
## Only have yellow perch from NS habitat (n=9)
## CHECK MIXSIAR MODEL FOR BH (following Smith et al. 2013)
sources <- read.table("BH_MixSIAR_Sources_test.csv",header=T,sep=",") #always put 13C(x) before 15N(y)
mixture <- read.table("BH_MixSIAR_Consumers_test.csv",header=T,sep=",")
TEF <-  read.table("BH_MixSIAR_TDF_test.csv", header=T,sep=",")
its <- 1500  #specify number of iterations
min_C <- -50  #specify dimensions and resolution for the mixing region figure
max_C <- -20    #choose values outside the 95% mixing region
min_N <- -2  
max_N <- 10
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C)
N_g <- seq(min_N,max_N,by=step_N)
mgrid <- function(a,b) {
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3))) 
p <- array(0, c(its,(nrow(mixture))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) { 
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TEF),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  
    f[j,1] <- rnorm(1, mean=TEF[j,1], sd=TEF[j,2])  
    f[j,2] <- rnorm(1, mean=TEF[j,3], sd=TEF[j,4])  
  }
  V <- v+f
  hull <- chull(V)
  hull_a <- append(hull,hull[1])
  P <- point.in.polygon(mixture[,1], mixture[,2], V[hull_a,1], V[hull_a,2])
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])
  m$y_f <- m$y[res:1,]
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2])
  m_r_s <- matrix(m_r,nrow=res,byrow=F)
  m_r_s[m_r_s > 1] <- 1
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)
  Par_values[i,] <- vals
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\n")) 
}

##FIGURE 1: variance of polygon area during simulation
Iterations <- Par_values[,ncol(Par_values)-1]
Variance <- Par_values[,ncol(Par_values)]
plot(Iterations, Variance, type="n")
lines(Iterations, Variance, lty=1, lwd=1.5, col="blue")
##FIGURE 2: proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1
Probabilities <- colSums(p)/its
print(Probabilities)
windows() 
barplot(Probabilities, xlab="Consumer", ylab="Probability consumer is within mixing polygon", 
        ylim=c(0,1), names.arg=seq(1,nrow(mixture),by=1))
##FIGURE 3: mixing region, consumers, average enriched source signatures
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TEF
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(mixture, pch=19, cex=1.3)
dev.copy2pdf(file="Mix_Region.pdf")
windows()  #create colour bar for figure 3
cust_color <- colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))
z <- matrix(1:100, nrow=1)
x <- 1
y <- seq(0,1,len=100)
image(x,y,z,col=colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))(100), 
      xaxt="n", xlab="", ylab="", useRaster=TRUE, bty="n", las=1)


#MixSIAR Model run for: BH NS LARGE PERCH
# Consumer Group: 1 = NS Large, n=3 individuals
mix.filename <- ("BH_MixSIAR_Consumers_LARGE_ONLY.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Sources: WL and NS prey fish, WL and NS invertebrates, WL seston
source.filename <- ("BH_MixSIAR_Sources_3.3.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load the discrimination data
# Used McCutchan's et al. (2003) TDF values
discr.filename <- ("BH_MixSIAR_TDF_3.3.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_BH_YP_LARGE_3.4.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
bh_ns_large_3.4 <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Examine output
output_JAGS(bh_ns_large_3.4, mix, source, output_options)
# NOTE: There is no global/overall estimated diet b/c we fit Group as a fixed effect instead of a random

#MixSIAR for: BH NS SMALL PERCH
# Consumer Group 1 = NS Small fish ("BH_MixSIAR_Consumers_SMALL_ONLY.csv"), n=6 individuals
mix.filename <- ("BH_MixSIAR_Consumers_SMALL_ONLY.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Sources are: WL and NS prey fish, wL and NS invertebrates, WL seston
source.filename <- ("BH_MixSIAR_Sources_3.3.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load the discrimination/TDF data
# Used McCutchan's et al. (2003) TDF values
discr.filename <- ("BH_MixSIAR_TDF_3.3.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_BH_YP_SMALL_3.4.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
bh_ns_small_3.4 <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(bh_ns_small_3.4, mix, source, output_options)
# NOTE: There is no global/overall estimated diet b/c we fit Group as a fixed effect instead of a random


################Calumet#####################
### South Region
## 2 habitats (NS, WL), n=34
## CHECK MIXSIAR MODEL FOR CA (this was repeated for all groups)
sources <- read.table("CA_MixSIAR_NS_Sources_test_2.25.csv",header=T,sep=",") #always put 13C(x) before 15N(y)
mixture <- read.table("CA_MixSIAR_NS_Consumers_test.csv",header=T,sep=",")
TEF <-  read.table("CA_MixSIAR_NS_TDF_test_2.25.csv", header=T,sep=",")
its <- 1500  #specify number of iterations
min_C <- -50  #specify dimensions and resolution for the mixing region figure
max_C <- -20    #choose values outside the 95% mixing region
min_N <- -2  
max_N <- 10
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C)
N_g <- seq(min_N,max_N,by=step_N)
mgrid <- function(a,b) {
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3))) 
p <- array(0, c(its,(nrow(mixture))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) { 
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TEF),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  
    f[j,1] <- rnorm(1, mean=TEF[j,1], sd=TEF[j,2])  
    f[j,2] <- rnorm(1, mean=TEF[j,3], sd=TEF[j,4])  
  }
  V <- v+f
  hull <- chull(V)
  hull_a <- append(hull,hull[1])
  P <- point.in.polygon(mixture[,1], mixture[,2], V[hull_a,1], V[hull_a,2])
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])
  m$y_f <- m$y[res:1,]
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2])
  m_r_s <- matrix(m_r,nrow=res,byrow=F)
  m_r_s[m_r_s > 1] <- 1
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)
  Par_values[i,] <- vals
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\n")) 
}

##FIGURE 1: variance of polygon area during simulation
Iterations <- Par_values[,ncol(Par_values)-1]
Variance <- Par_values[,ncol(Par_values)]
plot(Iterations, Variance, type="n")
lines(Iterations, Variance, lty=1, lwd=1.5, col="blue")
##FIGURE 2: proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1
Probabilities <- colSums(p)/its
print(Probabilities)
windows() 
barplot(Probabilities, xlab="Consumer", ylab="Probability consumer is within mixing polygon", 
        ylim=c(0,1), names.arg=seq(1,nrow(mixture),by=1))
##FIGURE 3: mixing region, consumers, average enriched source signatures
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TEF
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(mixture, pch=19, cex=1.3)
dev.copy2pdf(file="Mix_Region.pdf")
windows()  #create colour bar for figure 3
cust_color <- colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))
z <- matrix(1:100, nrow=1)
x <- 1
y <- seq(0,1,len=100)
image(x,y,z,col=colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))(100), 
      xaxt="n", xlab="", ylab="", useRaster=TRUE, bty="n", las=1)
#MixSIAR for: CA NS Large Perch
# Consumer Group: 1 = NS fish; n=4
mix.filename <- ("CA_MixSIAR_NS_Consumers.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Sources are: NS Alewives, NS Invertebrates, NS Seston, NS Shiners, WL Phragmites Invertebrates, WL SAV Invertebrates
source.filename <- ("CA_MixSIAR_NS_Sources_2.26.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load the discrimination/TDF data
# Used McCutchan's et al. (2003) TDF values
discr.filename <- ("CA_MixSIAR_NS_TDF_2.26.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_CA_NS_2.26.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
jags.1 <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(jags.1, mix, source, output_options)

# NOTE: No global/overall estimated diet b/c we fit Group as a fixed effect instead of a random

#MixSIAR for: CA WL Small Perch
# Consumer Group: 1 = WL Small fish (<150 mm TL)
mix.filename <- ("CA_MixSIAR_WL_SMALL_Consumers.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Sources are: WL, NS, prey fish and WL, NS, invertebrates, WL, NS Seston
source.filename <- ("CA_MixSIAR_WL_SMALL_Sources.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load the discrimination/TDF data
# Used McCutchan's et al. (2003) TDF values
discr.filename <- ("CA_MixSIAR_WL_SMALL_TDF.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_CA_WL_SMALL.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
jags.1 <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(jags.1, mix, source, output_options)

# NOTE: There is no global/overall estimated diet b/c we fit Group as a fixed effect instead of a random

# MixSIAR for: CA WL Large Perch
# Consumer Group: 1 = WL Large fish (>150 mm TL)
mix.filename <- ("CA_MixSIAR_WL_LARGE_Consumers.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Sources are: WL, NS, prey fish and WL, NS, invertebrates, WL, NS Seston
source.filename <- ("CA_MixSIAR_WL_LARGE_Sources.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load the discrimination/TDF data
# Used McCutchan's et al. (2003) TDF values
discr.filename <- ("CA_MixSIAR_WL_LARGE_TDF.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_CA_WL_LARGE.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
jags.1 <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(jags.1, mix, source, output_options)
# NOTE: There is no global/overall estimated diet b/c we fit Group as a fixed effect instead of a random

###############Cedar River####################
## West Region
# Only have individuals from WL habitat
# CHECK MIXSIAR MODEL FOR CR
sources <- read.table("CR_MixSIAR_SOURCES_test.csv",header=T,sep=",") #always put 13C(x) before 15N(y)
mixture <- read.table("CR_MixSIAR_Consumers_test.csv",header=T,sep=",")
TEF <-  read.table("CR_MixSIAR_TDF_test.csv", header=T,sep=",")
its <- 1500  #specify number of iterations
min_C <- -50  #specify dimensions and resolution for the mixing region figure
max_C <- -20    #choose values outside the 95% mixing region
min_N <- -2  
max_N <- 10
res <- 250 #resolution of the mixing region figure; reducing this improves performance

## RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C)
N_g <- seq(min_N,max_N,by=step_N)
mgrid <- function(a,b) {
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3))) 
p <- array(0, c(its,(nrow(mixture))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) { 
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TEF),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  
    f[j,1] <- rnorm(1, mean=TEF[j,1], sd=TEF[j,2])  
    f[j,2] <- rnorm(1, mean=TEF[j,3], sd=TEF[j,4])  
  }
  V <- v+f
  hull <- chull(V)
  hull_a <- append(hull,hull[1])
  P <- point.in.polygon(mixture[,1], mixture[,2], V[hull_a,1], V[hull_a,2])
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])
  m$y_f <- m$y[res:1,]
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2])
  m_r_s <- matrix(m_r,nrow=res,byrow=F)
  m_r_s[m_r_s > 1] <- 1
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)
  Par_values[i,] <- vals
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\n")) 
}

##FIGURE 1: variance of polygon area during simulation
Iterations <- Par_values[,ncol(Par_values)-1]
Variance <- Par_values[,ncol(Par_values)]
plot(Iterations, Variance, type="n")
lines(Iterations, Variance, lty=1, lwd=1.5, col="blue")
##FIGURE 2: proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1
Probabilities <- colSums(p)/its
print(Probabilities)
windows() 
barplot(Probabilities, xlab="Consumer", ylab="Probability consumer is within mixing polygon", 
        ylim=c(0,1), names.arg=seq(1,nrow(mixture),by=1))
##FIGURE 3: mixing region, consumers, average enriched source signatures
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TEF
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(mixture, pch=19, cex=1.3)
dev.copy2pdf(file="Mix_Region.pdf")
windows()  #create colour bar for figure 3
cust_color <- colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))
z <- matrix(1:100, nrow=1)
x <- 1
y <- seq(0,1,len=100)
image(x,y,z,col=colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))(100), 
      xaxt="n", xlab="", ylab="", useRaster=TRUE, bty="n", las=1)

#MixSIAR for CR WL LARGE Only
# Consumer Group: 1 = Large WL Fish ("CR_MixSIAR_WL_LARGE_Consumers.csv"), n=3
mix.filename <- ("CR_MixSIAR_WL_LARGE_Consumers.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Four Sources: WL Inverts, WL Prey Fish, NS Alewife, NS Inverts
source.filename <- ("CR_MixSIAR_WL_SOURCES.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load discrimination/TDF data
discr.filename <- ("CR_MixSIAR_WL_TDF.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=FALSE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_CR_WL_3.10_large.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
cr_wl_3.10_large <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(cr_wl_3.10_large, mix, source, output_options)

# NOTE: There is no global/overall estimated diet b/c we fit Group as a fixed effect instead of a random


############Peshtigo###################
## West Region
# Only have individuals from WL habitat
## CHECK PE MIXSIAR MODEL
sources <- read.table("PE_MixSIAR_SOURCES_test.csv",header=T,sep=",") #always put 13C(x) before 15N(y)
mixture <- read.table("PE_MixSIAR_Consumers_test.csv",header=T,sep=",")
TEF <-  read.table("PE_MixSIAR_TDF_test.csv", header=T,sep=",")
its <- 1500  #specify number of iterations
min_C <- -50  #specify dimensions and resolution for the mixing region figure
max_C <- -20    #choose values outside the 95% mixing region
min_N <- -2  
max_N <- 10
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##RUN simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C)
N_g <- seq(min_N,max_N,by=step_N)
mgrid <- function(a,b) {
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3))) 
p <- array(0, c(its,(nrow(mixture))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) { 
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TEF),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  
    f[j,1] <- rnorm(1, mean=TEF[j,1], sd=TEF[j,2])  
    f[j,2] <- rnorm(1, mean=TEF[j,3], sd=TEF[j,4])  
  }
  V <- v+f
  hull <- chull(V)
  hull_a <- append(hull,hull[1])
  P <- point.in.polygon(mixture[,1], mixture[,2], V[hull_a,1], V[hull_a,2])
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])
  m$y_f <- m$y[res:1,]
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2])
  m_r_s <- matrix(m_r,nrow=res,byrow=F)
  m_r_s[m_r_s > 1] <- 1
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)
  Par_values[i,] <- vals
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\n")) 
}

##FIGURE 1: variance of polygon area during simulation
Iterations <- Par_values[,ncol(Par_values)-1]
Variance <- Par_values[,ncol(Par_values)]
plot(Iterations, Variance, type="n")
lines(Iterations, Variance, lty=1, lwd=1.5, col="blue")
##FIGURE 2: proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1
Probabilities <- colSums(p)/its
print(Probabilities)
windows() 
barplot(Probabilities, xlab="Consumer", ylab="Probability consumer is within mixing polygon", 
        ylim=c(0,1), names.arg=seq(1,nrow(mixture),by=1))
##FIGURE 3: mixing region, consumers, average enriched source signatures
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TEF
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(mixture, pch=19, cex=1.3)
dev.copy2pdf(file="Mix_Region.pdf")
windows()  #create colour bar for figure 3
cust_color <- colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))
z <- matrix(1:100, nrow=1)
x <- 1
y <- seq(0,1,len=100)
image(x,y,z,col=colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))(100), 
      xaxt="n", xlab="", ylab="", useRaster=TRUE, bty="n", las=1)

#MixSIAR for: PE WL LARGE Perch
# Consumer Group: 1 = Large WL Fish ("CR_MixSIAR_WL_LARGE_Consumers.csv"), n=3
mix.filename <- ("PE_MixSIAR_WL_LARGE_Consumers.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Four Sources: WL Inverts, WL YOY Perch, NS Alewife, WL Pumpkinseed
source.filename <- ("PE_MixSIAR_WL_SOURCES.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load the discrimination/TDF data
discr.filename <- ("PE_MixSIAR_WL_TDF.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=FALSE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_PE_WL_3.10_large.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
pe_wl_3.10_large <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(pe_wl_3.10_large, mix, source, output_options)

# NOTE: There is no global/overall estimated diet b/c we fit Group as a fixed effect instead of a random


############Little Sturgeon###############
## West Region
# 2 habitats (NS, WL); n=18
## CHECK LS MIXSIAR MODEL (this was repeated for all groups)
sources <- read.table("LS_MixSIAR_NS_Sources_test.csv",header=T,sep=",") #always put 13C(x) before 15N(y)
mixture <- read.table("LS_MixSIAR_NS_Consumers_test.csv",header=T,sep=",")
TEF <-  read.table("LS_MixSIAR_NS_TDF_test.csv", header=T,sep=",")
its <- 1500  #specify number of iterations
min_C <- -50  #specify dimensions and resolution for the mixing region figure
max_C <- -20    #choose values outside the 95% mixing region
min_N <- -2  
max_N <- 10
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##RUN simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C)
N_g <- seq(min_N,max_N,by=step_N)
mgrid <- function(a,b) {
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3))) 
p <- array(0, c(its,(nrow(mixture))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) { 
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TEF),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  
    f[j,1] <- rnorm(1, mean=TEF[j,1], sd=TEF[j,2])  
    f[j,2] <- rnorm(1, mean=TEF[j,3], sd=TEF[j,4])  
  }
  V <- v+f
  hull <- chull(V)
  hull_a <- append(hull,hull[1])
  P <- point.in.polygon(mixture[,1], mixture[,2], V[hull_a,1], V[hull_a,2])
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])
  m$y_f <- m$y[res:1,]
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2])
  m_r_s <- matrix(m_r,nrow=res,byrow=F)
  m_r_s[m_r_s > 1] <- 1
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)
  Par_values[i,] <- vals
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\n")) 
}

##FIGURE 1: variance of polygon area during simulation
Iterations <- Par_values[,ncol(Par_values)-1]
Variance <- Par_values[,ncol(Par_values)]
plot(Iterations, Variance, type="n")
lines(Iterations, Variance, lty=1, lwd=1.5, col="blue")
##FIGURE 2: proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1
Probabilities <- colSums(p)/its
print(Probabilities)
windows() 
barplot(Probabilities, xlab="Consumer", ylab="Probability consumer is within mixing polygon", 
        ylim=c(0,1), names.arg=seq(1,nrow(mixture),by=1))
##FIGURE 3: mixing region, consumers, average enriched source signatures
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TEF
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(mixture, pch=19, cex=1.3)
dev.copy2pdf(file="Mix_Region.pdf")
windows()  #create colour bar for figure 3
cust_color <- colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))
z <- matrix(1:100, nrow=1)
x <- 1
y <- seq(0,1,len=100)
image(x,y,z,col=colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))(100), 
      xaxt="n", xlab="", ylab="", useRaster=TRUE, bty="n", las=1)

# MixSIAR Model for: NS Small Perch
# Consumer Group: 1 = Small  ("LS_MixSIAR_NS_SMALL_Consumers.csv"), n=3 individuals
mix.filename <- ("LS_MixSIAR_NS_SMALL_Consumers.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Sources are: WL, NS prey fish and WL, NS, invertebrates
source.filename <- ("LS_MixSIAR_NS_Sources_3.2.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load discrimination/TDF data
# Used McCutchan's et al. (2003) TDF values
discr.filename <- ("LS_MixSIAR_NS_TDF_3.2.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_LS_NS_Small_3.3.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
mod_3.3_small_LSNS <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(mod_3.3_small_LSNS, mix, source, output_options)
# NOTE: There is no global/overall estimated diet b/c we fit Group as a fixed effect instead of a random

# MixSIAR model for: WL Small Perch
# Consumer Group: 1 = WL Small  ("LS_MixSIAR_WL_SMALL_Consumers.csv"), n=4 individuals
mix.filename <- ("LS_MixSIAR_WL_SMALL_Consumers.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Sources are: WL, NS prey fish and WL, NS, invertebrates
source.filename <- ("LS_MixSIAR_NS_Sources_3.2.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load discrimination/TDF data
# Used McCutchan's et al. (2003) TDF values
discr.filename <- ("LS_MixSIAR_NS_TDF_3.2.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_LS_WL_Small_3.3.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
mod_3.3_small_LSWL <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(mod_3.3_small_LSWL, mix, source, output_options)
# NOTE: There is no global/overall estimated diet b/c we fit Group as a fixed effect instead of a random

# MixSIAR Model: WL Large Perch
# Consumer Group: 1 = WL Large  ("LS_MixSIAR_WL_LARGE_Consumers.csv"), n=3 individuals
mix.filename <- ("LS_MixSIAR_WL_LARGE_Consumers.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Sources are: WL, NS prey fish and WL, NS, invertebrates
source.filename <- ("LS_MixSIAR_NS_Sources_3.2.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load discrimination/TDF data
# Used McCutchan's et al. (2003) TDF values
discr.filename <- ("LS_MixSIAR_NS_TDF_3.2.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_LS_WL_Large_3.2.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
mod_3.2_large_LSWL <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(mod_3.2_large_LSWL, mix, source, output_options)
# NOTE: There is no global/overall estimated diet b/c we fit Group as a fixed effect instead of a random

## MixSIAR model for: NS Large Perch
# Consumer Group: 1 = Large  ("LS_MixSIAR_NS_LARGE_Consumers.csv"), n=3 individuals
mix.filename <- ("LS_MixSIAR_NS_LARGE_Consumers.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Sources are: WL, NS prey fish and WL, NS, invertebrates
source.filename <- ("LS_MixSIAR_NS_Sources_3.2.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load discrimination/TDF data
# Used McCutchan's et al. (2003) TDF values
discr.filename <- ("LS_MixSIAR_NS_TDF_3.2.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_LS_NS_Large_3.2.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
mod_3.2_large_LSNS <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(mod_3.2_large_LSNS, mix, source, output_options)
# NOTE: There is no global/overall estimated diet b/c we fit Group as a fixed effect instead of a random

#############Muskegon####################
## East Region
## CHECK MIXSIAR MODEL FOR MU (run for each group)
sources <- read.table("MU_MixSIAR_NS_Sources_test.csv",header=T,sep=",") #always put 13C(x) before 15N(y)
mixture <- read.table("MU_MixSIAR_NS_Consumers_test.csv",header=T,sep=",")
TEF <-  read.table("MU_MixSIAR_NS_TDF_test.csv", header=T,sep=",")
its <- 1500  #specify number of iterations
min_C <- -50  #specify dimensions and resolution for the mixing region figure
max_C <- -20    #choose values outside the 95% mixing region
min_N <- -2  
max_N <- 10
res <- 250 #resolution of mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C)
N_g <- seq(min_N,max_N,by=step_N)
mgrid <- function(a,b) {
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3))) 
p <- array(0, c(its,(nrow(mixture))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) { 
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TEF),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  
    f[j,1] <- rnorm(1, mean=TEF[j,1], sd=TEF[j,2])  
    f[j,2] <- rnorm(1, mean=TEF[j,3], sd=TEF[j,4])  
  }
  V <- v+f
  hull <- chull(V)
  hull_a <- append(hull,hull[1])
  P <- point.in.polygon(mixture[,1], mixture[,2], V[hull_a,1], V[hull_a,2])
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])
  m$y_f <- m$y[res:1,]
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2])
  m_r_s <- matrix(m_r,nrow=res,byrow=F)
  m_r_s[m_r_s > 1] <- 1
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)
  Par_values[i,] <- vals
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\n")) 
}

##FIGURE 1: variance of polygon area during simulation
Iterations <- Par_values[,ncol(Par_values)-1]
Variance <- Par_values[,ncol(Par_values)]
plot(Iterations, Variance, type="n")
lines(Iterations, Variance, lty=1, lwd=1.5, col="blue")
##FIGURE 2: proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1
Probabilities <- colSums(p)/its
print(Probabilities)
windows() 
barplot(Probabilities, xlab="Consumer", ylab="Probability consumer is within mixing polygon", 
        ylim=c(0,1), names.arg=seq(1,nrow(mixture),by=1))
##FIGURE 3: mixing region, consumers, average enriched source signatures
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TEF
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(mixture, pch=19, cex=1.3)
dev.copy2pdf(file="Mix_Region.pdf")
windows()  #create colour bar for figure 3
cust_color <- colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))
z <- matrix(1:100, nrow=1)
x <- 1
y <- seq(0,1,len=100)
image(x,y,z,col=colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))(100), 
      xaxt="n", xlab="", ylab="", useRaster=TRUE, bty="n", las=1)

#MixSIAR for: MU NS Large Perch
# Consumer Group: 1 = Large NS Fish ("MU_MixSIAR_NS_Consumers.csv"), removed 2 perch w/ low probabiltiies
mix.filename <- ("MU_MixSIAR_NS_Consumers_v2.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Five Sources: WL Seston, DRM Seston, NS Shiner, WL Invert, DRM Shiner
source.filename <- ("MU_MixSIAR_NS_SOURCES.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load the discrimination/TDF data
discr.filename <- ("MU_MixSIAR_NS_TDF.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=FALSE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_MU_NS_3.9_v2.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
mu_ns_3.9_v2 <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(mu_ns_3.9_v2, mix, source, output_options)
# NOTE: There is no global/overall estimated diet b/c we fit Group as a fixed effect instead of a random

#MixSIAR for: MU WL SMALL Perch
# Consumer Group: 1 = Small WL Fish ("MU_MixSIAR_WL_SMALL_Consumers.csv"), removed 2 perch w/ low probabiltiies
mix.filename <- ("MU_MixSIAR_WL_SMALL_Consumers_v2.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Sources: WL Seston, DRM Seston, NS Shiner, DRM Shiner
source.filename <- ("MU_MixSIAR_WL_SOURCES_SMALL_PERCH_3.8_4source.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load discrimination/TDF data
discr.filename <- ("MU_MixSIAR_WL_TDF_SMALL_PERCH_3.8_4source.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=FALSE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_MU_WL_SMALL_3.8_4source.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
mu_wl_small_3.8_4source <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(mu_wl_small_3.8_4source, mix, source, output_options)
# NOTE: There is no global/overall estimated diet b/c we fit Group as a fixed effect instead of a random

# MixSIAR for: MU WL LARGE Perch
# Consumer Group: 1 = LARGE WL Fish ("MU_MixSIAR_WL_LARGE_Consumers.csv"), n=13 individuals
mix.filename <- ("MU_MixSIAR_WL_LARGE_Consumers.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Sources: WL Small Perch, DRM Seston, NS Shiner, NS Invert
source.filename <- ("MU_MixSIAR_WL_SOURCES_3.7.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load discrimination/TDF data
discr.filename <- ("MU_MixSIAR_WL_TDF_3.7.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=FALSE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_MU_WL_LARGE_3.7.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
mu_wl_large_3.7 <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(mu_wl_large_3.7, mix, source, output_options)
# NOTE: There is no global/overall estimated diet b/c we fit Group as a fixed effect instead of a random

##############Pentwater###################
## East Region
# 1 habitat (WL); n=10
## CHECK MIXSIAR MODEL FOR PW
sources <- read.table("PW_MixSIAR_WL_Sources_test.csv",header=T,sep=",") #always put 13C(x) before 15N(y)
mixture <- read.table("PW_MixSIAR_WL_Consumers_test.csv",header=T,sep=",")
TEF <-  read.table("PW_MixSIAR_TDF_test.csv", header=T,sep=",")
its <- 1500  #specify number of iterations
min_C <- -50  #specify dimensions and resolution for the mixing region figure
max_C <- -20    #choose values outside 95% mixing region
min_N <- -2  
max_N <- 10
res <- 250 #resolution of mixing region figure; reducing this improves performance

##RUN simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C)
N_g <- seq(min_N,max_N,by=step_N)
mgrid <- function(a,b) {
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3))) 
p <- array(0, c(its,(nrow(mixture))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) { 
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TEF),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  
    f[j,1] <- rnorm(1, mean=TEF[j,1], sd=TEF[j,2])  
    f[j,2] <- rnorm(1, mean=TEF[j,3], sd=TEF[j,4])  
  }
  V <- v+f
  hull <- chull(V)
  hull_a <- append(hull,hull[1])
  P <- point.in.polygon(mixture[,1], mixture[,2], V[hull_a,1], V[hull_a,2])
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])
  m$y_f <- m$y[res:1,]
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2])
  m_r_s <- matrix(m_r,nrow=res,byrow=F)
  m_r_s[m_r_s > 1] <- 1
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)
  Par_values[i,] <- vals
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\n")) 
}

##FIGURE 1: variance of polygon area during simulation
Iterations <- Par_values[,ncol(Par_values)-1]
Variance <- Par_values[,ncol(Par_values)]
plot(Iterations, Variance, type="n")
lines(Iterations, Variance, lty=1, lwd=1.5, col="blue")
##FIGURE 2: proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1
Probabilities <- colSums(p)/its
print(Probabilities)
windows() 
barplot(Probabilities, xlab="Consumer", ylab="Probability consumer is within mixing polygon", 
        ylim=c(0,1), names.arg=seq(1,nrow(mixture),by=1))
##FIGURE 3: mixing region, consumers, average enriched source signatures
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TEF
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(mixture, pch=19, cex=1.3)
dev.copy2pdf(file="Mix_Region.pdf")
windows()  #create colour bar for figure 3
cust_color <- colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))
z <- matrix(1:100, nrow=1)
x <- 1
y <- seq(0,1,len=100)
image(x,y,z,col=colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))(100), 
      xaxt="n", xlab="", ylab="", useRaster=TRUE, bty="n", las=1)

#MixSIAR for: PW WL SMALL Perch
# Consumer Groups: 1 = WL fish <150 mm ("PW_MixSIAR_WL_SMALL_Consumers_3.4.csv"), n=6 individuals
mix.filename <- ("PW_MixSIAR_WL_SMALL_Consumers_3.4.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Sources are: WL, NS prey fish and WL, NS, DRM invertebrates 
source.filename <- ("PW_MixSIAR_WL_Sources_3.4.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load discrimination/TDF data
discr.filename <- ("PW_MixSIAR_WL_TDF_3.4.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_PW_WL_SMALL_YP_3.5.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
pw_wl_small_3.5 <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(pw_wl_small_3.5, mix, source, output_options)
# NOTE: There is no global/overall estimated diet b/c we fit Group as a fixed effect instead of a random

#MixSIAR for: PW WL LARGE Perch
# Consumer Groups: 1 = WL fish >150 mm ("PW_MixSIAR_WL_LARGE_Consumers_3.5.csv"), n=4 individuals
mix.filename <- ("PW_MixSIAR_WL_LARGE_Consumers_3.5.csv")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="GROUP", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# Sources are: WL, NS prey fish and WL, NS, DRM invertebrates 
source.filename <- ("PW_MixSIAR_WL_Sources_3.5.csv")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

# Load discrimination/TDF data
discr.filename <- ("PW_MixSIAR_WL_TDF_3.5.csv")
discr <- load_discr_data(filename=discr.filename, mix)

### Plot data as an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

### Calculate convex hull area
# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s) as in Brett (2014).
# Discrimination SD is added to the source SD (see ?calc_area for details)
# If source data are by factor, computes area for each polygon (one for each of 3 regions in wolf ex)
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

### Plot prior
# RED = prior
# DARK GREY = "uninformative"/generalist (alpha = 1)
# LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source)

### Write JAGS model files
model_filename <- "MixSIAR_model_PW_WL_LARGE_YP_3.5.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run model
#Good idea to use run = "test" first to check if 1) the data are loaded correctly and 2) the model is specified correctly:
jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

# Once that works, try with MCMC value that will converge
# extreme = 3,000,000 chains
pw_wl_large_3.5 <- run_model(run="extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

#Set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Look at output
output_JAGS(pw_wl_large_3.5, mix, source, output_options)
# NOTE: No global/overall estimated diet b/c we fit Group as a fixed effect instead of a random

###########Comparing models#################
### LOO and WAIC are "methods for estimating pointwise out-of-sample prediction accuracy from a
# fitted Bayesian model using the log-likelihood evaluated at the posterior simulations of the parameter values". 
# See Vehtari, Gelman, & Gabry (2017)
# LOO and WAIC are preferred over AIC or DIC
# LOO is more robust than WAIC
#'loo' estimates standard errors for the difference in LOO/WAIC between two models
# We can calculate the relative support for each model using LOO/WAIC weights

#Compare NS model 1 vs. NS model 2
mod1<-jags.1
mod2<-jags.1

#Output: 
# Model: names of x (input list)
# LOOic / WAIC: LOO information criterion or WAIC
#se_LOOic / se_WAIC: standard error of LOOic / WAIC
#dLOOic / dWAIC: difference between each model and the model with lowest LOOic/WAIC --> Best model has dLOOic = 0.
# se_dLOOic / se_dWAIC: standard error of the difference between each model and the model with lowest LOOic/WAIC
# weight: relative support for each model, calculated as Akaike weights (p.75 Burnham & Anderson 2002). 
# Interpretation: "an estimate of the probability that the model will make the
# best predictions on new data, conditional on the set of models considered" (McElreath 2015).

x <- list(mod1, mod2)
names(x) <- c("With NS Fish","With NS Alewife + Shiner")
compare_models(x) #With NS Fish had more support, but models were very similar

