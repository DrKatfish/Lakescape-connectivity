### NS-WL Linkages MS
# KE O'Reilly et al.
# Code to recreate figures

#Load packages
library(tidyverse)
library(skimr)
library(RColorBrewer)
library(viridis)
library(FSA)
library(ggrepel)
library(gridExtra)
library(ggpubr)
library(ggsignif)
library(effsize)
library(MixSIAR)
library(mcmc)
library(vegan)
library(SIBER)
library(sp)
library(splancs)
library(patchwork)
library(ggeffects) #plots mixed models
library(ggbeeswarm)
library(emmeans)
library(cowplot)
library(raster)
library(maps)
library(effects)
library(maptools)
library(mapdata)
library(scales)

############Figure 1#############
## Figure 1. Map
##Lake Michigan wetland-nearshore sites
map("usa", col="grey", xlim=c(-95,-80), ylim=c(40,50), fill=TRUE)
lmich <- readShapePoly("hydro_p_LakeMichigan.shp")
#Plot shapefile
lmich_tidy <- tidy(lmich)
plot(lmich, add=TRUE, xlim=c(-90,-85), ylim=c(40,47), col=alpha("blue", 0.6), border=FALSE)

samps <- read.csv("Site_Lat_Lon.csv")   #data for sampling sites, contains column of "lat" and column of "lon" with GPS points in decimal degrees

ggplot(samps,
       aes(x = lon,
           y = lat))+
  borders("world",
          xlim = c(-100, -50),
          ylim = c(40,50)) + 
  geom_point()+
  theme_bw()

##Final Map
a<-ggplot(samps,
          aes(x = lon,
              y = lat))+
  borders("world",
          xlim = c(-100, -50),
          ylim = c(40,50)) +
  geom_polygon(aes(x=long, y=lat, group=group), data=lmich_tidy, fill=alpha("blue", 0.6))+
  borders("state", fill = "grey")+
  geom_point(aes(shape = Region,
                 fill = Region), size = 10)+
  geom_text_repel(label=samps$Site, size=8, point.padding = 1)+
  scale_fill_manual(values = c("red",
                               "yellow",
                               "green"),
                    name = "")+
  scale_shape_manual(values=c(21,22,24), name="")+
  ##Zoom to specific x and y limits
  coord_cartesian(ylim = c(40, 47),
                  xlim = c(-90,-83))+
  xlab("Longitude")+ylab("Latitude")+
  theme_bw()+
  theme(axis.title.x=element_text(size=20, face="bold"), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=16), 
        axis.text.y=element_text(size=16),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14, face="bold"))

a

############Figure 2#############
## Figure 2
# Fig. 2A. DIC Differences
dic<-read_csv("DIC_Sites_Compare.csv")
a<-ggplot(dic, aes(fct_reorder(Site, regnum), y=d13C, fill=Habitat))
DIC.fig<-a+geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=d13C-SD, ymax=d13C+SD), colour="black", width=.2,
                position=position_dodge(.9))+ 
  scale_fill_manual(values=c("orange", "steelblue", "forestgreen"))+
  xlab("Site")+
  ylab(expression(paste(delta^{13}, "C (\u2030)")))+theme_bw()+
  theme(axis.title.x=element_text(size=20, face="bold"), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=18), 
        axis.text.y=element_text(size=18),
        legend.position='none')
DIC.fig

# Figures 2B,C, and D: Prey Fish, Invert, Seston
df<-read_csv("Source_Compare_Isotope.csv")
# FACTORS:
# Three regions (Region): East, South, and West
# Seven sites (Site): See Table 1 for site codes
# Habitat = NS (nearshore), WL (wetland), DRM (drowned river-mouth lake)
# Type = Invertebrate, Prey Fish, or Seston
# Group = corresponds with MixSIAR group
# Taxon = Specific taxon name
# FFG = functional feeding group
df$Region<-as.factor(df$Region)
df$Site<-as.factor(df$Site)
df$Habitat<-as.factor(df$Habitat)
df$Group<-as.factor(df$Group)
df$Type<-as.factor(df$Type)
df$Taxon<-as.factor(df$Taxon)
df$FFG<-as.factor(df$FFG)

### Isotopic differences among prey fish groups
# Biplot with mean (+/- SD) isotope values for each group
prey.df<-filter(df, Type=="Prey Fish")

# Figure 2B: Plot mean d13C and d15N of fish taxa across Region
PREYFISHmeanisotope<-prey.df %>%
  group_by(Site, Habitat) %>%
  summarise(mean_d13C= mean(d13C) , sd_d13C = sd(d13C), mean_d15N = mean(d15N) , sd_d15N = sd(d15N))

# Y axis = mean d13C/d15N of prey fish (+/- SD)
#### SET JITTER
myjit <- ggproto("fixJitter", PositionDodge,
                 width = 0.3,
                 dodge.width = 0.1,
                 jit = NULL,
                 compute_panel =  function (self, data, params, scales) 
                 {
                   
                   #Generate Jitter if not yet
                   if(is.null(self$jit) ) {
                     self$jit <-jitter(rep(0, nrow(data)), amount=self$dodge.width)
                   }
                   
                   data <- ggproto_parent(PositionDodge, self)$compute_panel(data, params, scales)
                   
                   data$x <- data$x + self$jit
                   #For proper error extensions
                   if("xmin" %in% colnames(data)) data$xmin <- data$xmin + self$jit
                   if("xmax" %in% colnames(data)) data$xmax <- data$xmax + self$jit
                   data
                 } )

PREYFISHmeanisotope$SiteOrder <- factor(PREYFISHmeanisotope$Site,levels = c("BH", "CA", "MU", "PW", "CR", "LS", "PE"))

prey.fig<-ggplot(PREYFISHmeanisotope, aes(SiteOrder, mean_d13C))+geom_point(aes(color=Habitat, shape=Habitat),
                                                                            size=9, position=myjit, alpha=0.85)+
  geom_errorbar(aes(ymin=mean_d13C-sd_d13C, ymax=mean_d13C+sd_d13C, color=Habitat), width=0.5, position=myjit)+ 
  xlab("Site")+ ylab(expression(paste(delta^{13}, "C (\u2030)")))+theme_bw()+
  scale_color_manual(values=c("orange", "steelblue", "forestgreen"))+
  theme(axis.title.x=element_text(size=20, face="bold"), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=18, vjust=0.5), 
        axis.text.y=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"))
prey.fig

## d15N
c<-ggplot(tmeanisotope, aes(Region, mean_d15N))+geom_point(aes(color=Taxon, shape=Region),size=8, position=position_dodge(width=0.5))+
  facet_grid(.~Taxon)+
  geom_errorbar(aes(ymin=mean_d15N-sd_d15N, ymax=mean_d15N+sd_d15N), colour="black", position=position_dodge(width=0.5))+ 
  xlab("Region")+ ylab(expression(paste(delta^{15}, "N (\u2030)")))+theme_bw()+
  theme(axis.title.x=element_text(size=20, face="bold"), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=18, angle=45, vjust=0.5), 
        axis.text.y=element_text(size=18),
        strip.text.x = element_text(size = 14, colour = "black"),
        legend.position = "none")
ggarrange(b+ rremove("x.text")+rremove("xlab"), c, 
          ncol = 1, nrow = 2)

### Invertebrates
# Biplot with mean (+/- SD) isotope values for each group
invert.df<-filter(df, Type=="Invertebrate")

# Figure Plot mean d13C and d15N of Taxa across Region
INVERTmeanisotope<-invert.df %>%
  group_by(Site, Habitat) %>%
  summarise(mean_d13C= mean(d13C) , sd_d13C = sd(d13C), mean_d15N = mean(d15N) , sd_d15N = sd(d15N))

# Y axis = mean d13C/d15N of an invert taxa (+/- SD)
INVERTmeanisotope$SiteOrder <- factor(INVERTmeanisotope$Site,levels = c("BH", "CA", "MU", "PW", "CR", "LS", "PE"))

invert.fig<-ggplot(INVERTmeanisotope, aes(SiteOrder, mean_d13C))+geom_point(aes(color=Habitat, shape=Habitat),size=9, position=myjit, alpha=0.85)+
  geom_errorbar(aes(ymin=mean_d13C-sd_d13C, ymax=mean_d13C+sd_d13C, color=Habitat), width=0.5, position=myjit)+ 
  xlab("Site")+ ylab(expression(paste(delta^{13}, "C (\u2030)")))+theme_bw()+
  scale_color_manual(values=c("orange", "steelblue", "forestgreen"))+
  theme(axis.title.x=element_text(size=20, face="bold"), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=18, vjust=0.5), 
        axis.text.y=element_text(size=18),
        legend.position='none')
invert.fig

#d15N
c<-ggplot(tmeanisotope, aes(Region, mean_d15N))+geom_point(aes(color=Taxon, shape=Region),size=8)+
  facet_grid(.~Taxon)+
  geom_errorbar(aes(ymin=mean_d15N-sd_d15N, ymax=mean_d15N+sd_d15N), colour="black")+ 
  xlab("Region")+ ylab(expression(paste(delta^{15}, "N (\u2030)")))+theme_bw()+
  theme(axis.title.x=element_text(size=20, face="bold"), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=18, angle=45, vjust=0.5), 
        axis.text.y=element_text(size=18),
        strip.text.x = element_text(size = 14, colour = "black"),
        legend.position = "none")
ggarrange(b+ rremove("x.text")+rremove("xlab"), c, 
          ncol = 1, nrow = 2)

### Seston
# Biplot with mean (+/- SD) isotope values for each group
seston.df<-filter(df, Type=="Seston")

# Figure Plot mean d13C and d15N of Taxa across Region
SESTONmeanisotope<-seston.df %>%
  group_by(Site, Habitat) %>%
  summarise(mean_d13C= mean(d13C) , sd_d13C = sd(d13C), mean_d15N = mean(d15N) , sd_d15N = sd(d15N))

# Y axis = mean d13C/d15N of seston (+/- SD)
SESTONmeanisotope$SiteOrder <- factor(SESTONmeanisotope$Site,levels = c("BH", "CA", "MU", "PW", "CR", "LS", "PE"))

seston.fig<-ggplot(SESTONmeanisotope, aes(SiteOrder, mean_d13C))+geom_point(aes(color=Habitat, shape=Habitat),size=9,alpha=0.85, position=myjit)+
  geom_errorbar(aes(ymin=mean_d13C-sd_d13C, ymax=mean_d13C+sd_d13C, color=Habitat), width=0.5, position=myjit)+ 
  xlab("Site")+ ylab(expression(paste(delta^{13}, "C (\u2030)")))+theme_bw()+
  scale_color_manual(values=c("orange", "steelblue", "forestgreen"))+
  theme(axis.title.x=element_text(size=20, face="bold"), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=18, vjust=0.5), 
        axis.text.y=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"))
seston.fig

resource.fig<-(DIC.fig | prey.fig)/(invert.fig | seston.fig)
resource.fig

############Figure 3#############
#### Figure 3. Food webs
#Load prey resources and consumers (perch):
resource<-read_csv("MixSIAR_Resources_Means_for_Plotting_7.30.21.csv")
perch<-read_csv("MixSIAR_Consumers_for_plotting_7.30.21.csv")

#Set colors for habitats (NS=nearshore, WL = wetland, and DRM = drowned river-mouth lake)
habitat.colors<-c(NS="steelblue", WL="forestgreen", DRM="orange")

#Set shapes for groups (NS/WL Large Perch, NS/WL Small Perch)
shapes <- c(16,17,15,18)
names(shapes) <- c("NS Large", "WL Large", "WL Small", "NS Small")

## The following code creates d13C/d15N biplots for each site
# Squares (mean d13C/d15N +/- SD) = prey resources
# Individual points = perch grouped by habitat of collection and size (Large or Small)

## Little Sturgeon
# Resources:
LS_prey<-filter(resource, Site=="LS")
# Consumers: NS Large Perch, WL Large Perch, NS Small Perch, WL Small Perch
LS_perch<-filter(perch, Site=="LS")

jitter <- position_jitter(width = 0.1, height = 0.1)

lsfoodweb<-ggplot(data=LS_prey, aes(meand13C, meand15N))+geom_point(data=LS_prey, aes(fill=Habitat), shape=22, color="black", size=8)+
  geom_errorbar(data=LS_prey, aes(ymin=meand15N-sd15N, ymax=meand15N+sd15N), size=1)+
  geom_errorbarh(data=LS_prey, aes(xmin=meand13C-sd13C, xmax=meand13C+sd13C), size=1)+
  geom_point(data=LS_perch, aes(d13C, d15N, shape=Group, color=Habitat), size=9, alpha=0.6, position=jitter)+
  theme_classic()+xlab(expression(paste(delta^{13}, "C (\u2030)")))+ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  theme(axis.title.x=element_text(size=20, face="bold"), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=16), axis.text.y=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"))+
  scale_colour_manual(values=habitat.colors, guide=guide_legend(title="Habitat"))+
  scale_fill_manual(values=habitat.colors, guide=guide_legend(title="Habitat"))+
  scale_shape_manual(values=shapes)

lsfoodweb

## Peshtigo
# Resources:
PE_prey<-filter(resource, Site=="PE")
# Consumers: WL Large Perch
PE_perch<-filter(perch, Site=="PE")

jitter <- position_jitter(width = 0.1, height = 0.1)

pefoodweb<-ggplot(data=PE_prey, aes(meand13C, meand15N))+geom_point(data=PE_prey, aes(fill=Habitat), shape=22, color="black", size=8)+
  geom_errorbar(data=PE_prey, aes(ymin=meand15N-sd15N, ymax=meand15N+sd15N), size=1)+
  geom_errorbarh(data=PE_prey, aes(xmin=meand13C-sd13C, xmax=meand13C+sd13C), size=1)+
  geom_point(data=PE_perch, aes(d13C, d15N, shape=Group, color=Habitat), size=9, alpha=0.6, position=jitter)+
  theme_classic()+xlab(expression(paste(delta^{13}, "C (\u2030)")))+ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  theme(axis.title.x=element_text(size=20, face="bold"), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=16), axis.text.y=element_text(size=16),
        legend.position="none")+
  scale_colour_manual(values=habitat.colors, guide=guide_legend(title="Habitat"))+
  scale_fill_manual(values=habitat.colors, guide=guide_legend(title="Habitat"))+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_shape_manual(values=shapes)

pefoodweb

## Cedar River
# Resources:
CR_prey<-filter(resource, Site=="CR")
# Consumers: WL Large Perch
CR_perch<-filter(perch, Site=="CR")

jitter <- position_jitter(width = 0.1, height = 0.1)

crfoodweb<-ggplot(data=CR_prey, aes(meand13C, meand15N))+geom_point(data=CR_prey, aes(fill=Habitat), shape=22, color="black", size=8)+
  geom_errorbar(data=CR_prey, aes(ymin=meand15N-sd15N, ymax=meand15N+sd15N), size=1)+
  geom_errorbarh(data=CR_prey, aes(xmin=meand13C-sd13C, xmax=meand13C+sd13C), size=1)+
  geom_point(data=CR_perch, aes(d13C, d15N, shape=Group, color=Habitat), size=9, alpha=0.6, position=jitter)+
  theme_classic()+xlab(expression(paste(delta^{13}, "C (\u2030)")))+ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  theme(axis.title.x=element_text(size=20, face="bold"), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=16), axis.text.y=element_text(size=16),
        legend.position = "none")+
  scale_colour_manual(values=habitat.colors, guide=guide_legend(title="Habitat"))+
  scale_fill_manual(values=habitat.colors, guide=guide_legend(title="Habitat"))+
  scale_shape_manual(values=shapes)

crfoodweb

#Combine all WEST sites using patchwork
westsites<-crfoodweb+pefoodweb+lsfoodweb

## Burns Harbor
# Resources:
BH_prey<-filter(resource, Site=="BH")
# Consumers: NS Large Perch, NS Small Perch
BH_perch<-filter(perch, Site=="BH")

jitter <- position_jitter(width = 0.1, height = 0.1)

bhfoodweb<-ggplot(data=BH_prey, aes(meand13C, meand15N))+geom_point(data=BH_prey, aes(fill=Habitat), shape=22, color="black", size=8)+
  geom_errorbar(data=BH_prey, aes(ymin=meand15N-sd15N, ymax=meand15N+sd15N), size=1)+
  geom_errorbarh(data=BH_prey, aes(xmin=meand13C-sd13C, xmax=meand13C+sd13C), size=1)+
  geom_point(data=BH_perch, aes(d13C, d15N, shape=Group, color=Habitat), size=9, alpha=0.6, position=jitter)+
  theme_classic()+xlab(expression(paste(delta^{13}, "C (\u2030)")))+ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  theme(axis.title.x=element_text(size=20, face="bold"), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=16), axis.text.y=element_text(size=16),
        legend.position="none")+
  scale_colour_manual(values=habitat.colors, guide=guide_legend(title="Habitat"))+
  scale_fill_manual(values=habitat.colors, guide=guide_legend(title="Habitat"))+
  scale_shape_manual(values=shapes)

bhfoodweb

## Calumet
# Resources:
CA_prey<-filter(resource, Site=="CA")
# Consumers: NS Large Perch, WL Large Perch, WL Small Perch
CA_perch<-filter(perch, Site=="CA")

jitter <- position_jitter(width = 0.2, height = 0.1)

cafoodweb<-ggplot(data=CA_prey, aes(meand13C, meand15N))+geom_point(data=CA_prey, aes(fill=Habitat), shape=22, color="black", size=8)+
  geom_errorbar(data=CA_prey, aes(ymin=meand15N-sd15N, ymax=meand15N+sd15N), size=1)+
  geom_errorbarh(data=CA_prey, aes(xmin=meand13C-sd13C, xmax=meand13C+sd13C), size=1)+
  geom_point(data=CA_perch, aes(d13C, d15N, shape=Group, color=Habitat), size=9, alpha=0.6, position=jitter)+
  theme_classic()+xlab(expression(paste(delta^{13}, "C (\u2030)")))+ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  theme(axis.title.x=element_text(size=20, face="bold"), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=16), axis.text.y=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"))+
  scale_colour_manual(values=habitat.colors, guide=guide_legend(title="Habitat"))+
  scale_fill_manual(values=habitat.colors, guide=guide_legend(title="Habitat"))+
  scale_shape_manual(values=shapes)

cafoodweb

#Combine both SOUTH sites using Patchwork
southsites<-bhfoodweb+cafoodweb

## Pentwater
# Resources:
PW_prey<-filter(resource, Site=="PW")
# Consumers: WL Large Perch, WL Small Perch
PW_perch<-filter(perch, Site=="PW")

jitter <- position_jitter(width = 0.1, height = 0.1)

pwfoodweb<-ggplot(data=PW_prey, aes(meand13C, meand15N))+geom_point(data=PW_prey, aes(fill=Habitat), shape=22, color="black", size=8)+
  geom_errorbar(data=PW_prey, aes(ymin=meand15N-sd15N, ymax=meand15N+sd15N), size=1)+
  geom_errorbarh(data=PW_prey, aes(xmin=meand13C-sd13C, xmax=meand13C+sd13C), size=1)+
  geom_point(data=PW_perch, aes(d13C, d15N, shape=Group, color=Habitat), size=9, alpha=0.6, position=jitter)+
  theme_classic()+xlab(expression(paste(delta^{13}, "C (\u2030)")))+ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  theme(axis.title.x=element_text(size=20, face="bold"), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=16), axis.text.y=element_text(size=16),
        legend.position="none")+
  scale_colour_manual(values=habitat.colors, guide=guide_legend(title="Habitat"))+
  scale_fill_manual(values=habitat.colors, guide=guide_legend(title="Habitat"))+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_shape_manual(values=shapes)

pwfoodweb

## Muskegon
# Resources:
MU_prey<-filter(resource, Site=="MU")
# Consumers: WL Large Perch, WL Small Perch
MU_perch<-filter(perch, Site=="MU")

jitter <- position_jitter(width = 0.1, height = 0.1)

mufoodweb<-ggplot(data=MU_prey, aes(meand13C, meand15N))+geom_point(data=MU_prey, aes(fill=Habitat), shape=22, color="black", size=8)+
  geom_errorbar(data=MU_prey, aes(ymin=meand15N-sd15N, ymax=meand15N+sd15N), size=1)+
  geom_errorbarh(data=MU_prey, aes(xmin=meand13C-sd13C, xmax=meand13C+sd13C), size=1)+
  geom_point(data=MU_perch, aes(d13C, d15N, shape=Group, color=Habitat), size=9, alpha=0.6, position=jitter)+
  theme_classic()+xlab(expression(paste(delta^{13}, "C (\u2030)")))+ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  theme(axis.title.x=element_text(size=20, face="bold"), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=16), axis.text.y=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"))+
  guides(color=FALSE)+
  scale_colour_manual(values=habitat.colors, guide=guide_legend(title="Habitat"))+
  scale_fill_manual(values=habitat.colors, guide=guide_legend(title="Habitat"))+
  scale_shape_manual(values=shapes)

mufoodweb

#Combine both EAST sites using Patchwork
eastsites<-pwfoodweb+mufoodweb

############Figure 4##############
##Otolith microchemistry NMDS
# See Ns-WL Linkages Otolith Analysis.R #
## Multivariate analysis of edge chemistry
# Data are chemistry (Mg, Sr, Ba) of otolith edge (i.e., mean value of last 10 um) for adult yellow perch collected in 2014
# Sites are Muskegon, Calumet, and Little Sturgeon
#Load packages
library(pairwiseAdonis)
library(car)
library(ROCR)

#Load data
df<-read_csv("Otolith_Edge_Chemistry_Multivariate_Analysis.csv")

#Subset data
mu<-filter(df, Site=="MU") #Muskegon
ca<-filter(df, Site=="CA") #Calumet
ls<-filter(df, Site=="LS") #Little Sturgeon

# Non-Metric Multidimensional Scaling (NMDS)
df.nmds<-metaMDS(mu[,7:10], plot=TRUE)
plot(df.nmds, type = "t") #plotting in base R
str(df.nmds) # gives stress value for plot
stressplot(df.nmds) # gain stress plot for stress values

## Plot using ggplot2
data.scores <- as.data.frame(scores(df.nmds))
data.scores$site <- mu$ID
data.scores$grp <- mu$Habitat
head(data.scores)
species.scores <- as.data.frame(scores(df.nmds, "species"))  

#Uses the scores function from vegan to extract the species (i.e., elements) scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # creates a column of species, from the rownames of species.scores
head(species.scores)

## Create hulls
grp.nearshore <- data.scores[data.scores$grp == "NS", ][chull(data.scores[data.scores$grp == 
                                                                            "NS", c("NMDS1", "NMDS2")]), ]
grp.wetland <- data.scores[data.scores$grp == "WL", ][chull(data.scores[data.scores$grp == 
                                                                          "WL", c("NMDS1", "NMDS2")]), ] 

hull.data <- rbind(grp.nearshore, grp.wetland) #combine groups
hull.data

## Plot with hulls
mu.nmds.plot<-ggplot(data=data.scores, aes(x=NMDS1, y=NMDS2, shape=grp, colour=grp)) + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.30) + # add convex hulls
  geom_point(size=6) + # add point markers
  scale_colour_manual(values=c("NS" = "steelblue", "WL" = "forestgreen")) +
  scale_fill_manual(values=c("NS" = "steelblue", "WL" = "forestgreen")) +
  coord_equal() +
  theme_bw() + 
  guides(col=guide_legend("Habitat"), shape=guide_legend("Habitat"), fill=FALSE)+
  theme(axis.title.x = element_text(size=18, face="bold"), # remove x-axis labels
        axis.title.y = element_text(size=18, face="bold"),
        axis.text.x=element_text(size=16),
        axis.text.y = element_text(size=16),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.text = element_text(size=20),
        legend.title = element_text(size=24, face="bold"))

mu.nmds.plot

### PERMANOVA to test for differences
df.dist<-vegdist(mu[,7:10])
df.pcoa=cmdscale(df.dist)
plot(df.pcoa[,1],df.pcoa[,2])
df.perm<-adonis(df.dist~mu$Habitat) # p=0.145
df.perm
scatterplotMatrix(mu[7:10])

# Linear Discriminant Analysis with Jacknifed Prediction
# The code above performs an LDA, using listwise deletion of missing data. 
# CV=TRUE generates jacknifed (i.e., leave one out) predictions. 
fit <- lda(Habitat ~lnBa+Sr, data=mu, #removing Mg increases prediction accuracy
           na.action="na.omit", CV=TRUE)
fit # show results
# Assess the accuracy of the prediction
# percent correct for each category of Habitat
ct <- table(mu$Habitat, fit$class)
diag(prop.table(ct, 1)) # 85% accurate for NS, 67% accurate for WL
# total percent correct
sum(diag(prop.table(ct))) #79% accuracy discriminating between WL-NS fish at MU

# Overall Accuracy with Sr, Mg, ln-trans Ba and ln+1 -trans Mn= 79%
# NS accuracy = 85%, WL accuracy = 67%

#Overall accuracy with Ba and Sr = 84%
# NS accuracy = 85%, WL accuracy = 83%

### Split data
## 60% of the sample size
smp_size <- floor(0.6 * nrow(mu))

## set the seed to make partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(mu)), size = smp_size)

train <- mu[train_ind, ]
test <- mu[-train_ind, ]
#Train model using "known" data
train.fit <- lda(Habitat ~ Ba+Sr+Mg, data=train)
train.fit
#Predict
lda.pred <- predict(train.fit, test)
names(lda.pred)
table(lda.pred$class, test$Habitat)
mean(lda.pred$class==test$Habitat)
ldahist(lda.pred$x[,1], g= lda.pred$class)

predmodel.test.lda = predict(train.fit, newdata=test)
table(Predicted=predmodel.test.lda$class, Habitat=test$Habitat)

### CONSTRUCT ROC AUC PLOT:
# Get posteriors as a dataframe
lda.pred.posteriors <- as.data.frame(lda.pred$posterior)
# Evaluate model
pred <- prediction(lda.pred.posteriors[,2], test$Habitat)
roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
auc.train <- performance(pred, measure = "auc")
auc.train <- auc.train@y.values
# Plot
plot(roc.perf)
abline(a=0, b= 1)
text(x = .25, y = .65 ,paste("AUC = ", round(auc.train[[1]],3), sep = ""))

##Compare AUCs of MU Large Perch models:
## Mg, Sr, Ba: 0.857
# Mg only: 0.571
# Sr only: 0.429
# Ba only: 0.857
## Ba is driving separation between NS and WL large perch at Muskegon

## Calumet
# These data are the otolith edge (i.e., last 10 um) chemistry (mg, Sr, Ba) of adult yellow perch collected in 2014
ca<-read_csv("multivariate_otolith_edge_chem_CA_ONLY.csv") # Load data

#Convert Site, Habitat, and Region to factors
ca$Site<-as.factor(ca$Site)
ca$Habitat<-as.factor(ca$Habitat)
ca$Group<-as.factor(ca$Group)
ca$Size<-as.factor(ca$Size)

# Non-Metric Multidimensional Scaling (NMDS)
df.nmds<-metaMDS(ca[,7:10], plot=TRUE)
plot(df.nmds, type = "t") #plotting in base R
str(df.nmds) # gives stress value for plot
stressplot(df.nmds) # To gain the stress plot for stress values for your MDS

## To plot using ggplot2
data.scores <- as.data.frame(scores(df.nmds))
data.scores$site <- ca$ID
data.scores$grp <- ca$Habitat
head(data.scores)
species.scores <- as.data.frame(scores(df.nmds, "species"))  

#Uses the scores function from vegan to extract the species (i.e., elements) scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # creates a column of species, from the rownames of species.scores
head(species.scores)

### Plot
ggplot() + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=5) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site, colour=grp),size=5,vjust=2) +  # add the site labels
  scale_colour_manual(values=c("NS" = "steelblue", "WL"="forestgreen")) +
  coord_equal() +
  theme_bw()+
  theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
        axis.title.y = element_text(face='bold',size=16,vjust=1),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size=12,color='black'))

## Create hulls
grp.nearshore <- data.scores[data.scores$grp == "NS", ][chull(data.scores[data.scores$grp == 
                                                                            "NS", c("NMDS1", "NMDS2")]), ]
grp.wetland <- data.scores[data.scores$grp == "WL", ][chull(data.scores[data.scores$grp == 
                                                                          "WL", c("NMDS1", "NMDS2")]), ]

hull.data <- rbind(grp.nearshore, grp.wetland) #combine groups
hull.data

## Plot with hulls
ca.nmds.plot<-ggplot(data=data.scores, aes(x=NMDS1, y=NMDS2, shape=grp, colour=grp)) + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.30) + # add convex hulls
  geom_point(size=10) + # add point markers
  scale_colour_manual(values=c("NS" = "steelblue", "WL" = "forestgreen")) +
  scale_fill_manual(values=c("NS" = "steelblue", "WL" = "forestgreen")) +
  coord_fixed(ratio=3) +
  theme_bw() + 
  guides(col=guide_legend("Habitat"), shape=guide_legend("Habitat"), fill=FALSE)+
  theme(axis.title.x = element_text(size=18, face="bold"), # remove x-axis labels
        axis.title.y = element_text(size=18, face="bold"),
        axis.text.x=element_text(size=16),
        axis.text.y = element_text(size=16),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.position = "none")

ca.nmds.plot

### PERMANOVA
df.dist<-vegdist(ca[,7:10])
df.pcoa=cmdscale(df.dist)
plot(df.pcoa[,1],df.pcoa[,2])
#ANOSIM
df.perm<-adonis(df.dist~ca$Habitat)
df.perm #p=0.001

### Linear Discriminant Function Analysis (LDA)
scatterplotMatrix(ca[7:10])

# Linear Discriminant Analysis with Jacknifed Prediction
# Code above performs an LDA using listwise deletion of missing data 
# CV=TRUE generates jacknifed (i.e., leave one out) predictions 
ca.fit <- lda(Habitat ~ Sr+Mn+Mg+Ba, data=ca, 
              na.action="na.omit", CV=TRUE)
ca.fit # show results
# Assess the accuracy of the prediction
# percent correct for each category of Habitat
ct <- table(ca$Habitat, ca.fit$class)
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))

## LDA accuracy
#Single Elements
# Sr only: 100%
# Ba only: 93% (100% NS, 91% WL)
# Mg only: 87% (75% NS, 91% WL)
# Mn only: 87% (50% NS, 100% WL)

# LDA accuracy
#Element combos
# Ba+Mn: 87% (50% NS, 100% WL)
# Ba+Mg: 93% (75% NS, 100% WL)
# Mn+Mg: 100%

### Split data
## 60% of the sample size
smp_size <- floor(0.6 * nrow(ca))

## set the seed to make partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(ca)), size = smp_size)

train <- ca[train_ind, ]
test <- ca[-train_ind, ]
#Train model using "known" data
train.fit <- lda(Habitat ~ Ba+Sr+Mg+Mn, data=train)
train.fit
#Predict
lda.pred <- predict(train.fit, test)
names(lda.pred)
table(lda.pred$class, test$Habitat)
mean(lda.pred$class==test$Habitat)

### CONSTRUCT ROC AUC PLOT:
# Posteriors as a dataframe
lda.pred.posteriors <- as.data.frame(lda.pred$posterior)
# Evaluate  model
pred <- prediction(lda.pred.posteriors[,2], test$Habitat)
roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
auc.train <- performance(pred, measure = "auc")
auc.train <- auc.train@y.values
# Plot
plot(roc.perf)
abline(a=0, b= 1)
text(x = .25, y = .65 ,paste("AUC = ", round(auc.train[[1]],3), sep = ""))

###Compare AUCs of CA Large Perch models for Calumet
## Mg, Sr, Ba: 1
# Mg only: 1
# Sr only: 1
# Ba only: 1

## Little Sturgeon
#Convert Site, Habitat, and Region to factors
ls$Site<-as.factor(ls$Site)
ls$Habitat<-as.factor(ls$Habitat)
ls$Group<-as.factor(ls$Group)
ls$Size<-as.factor(ls$Size)

# Non-Metric Multidimensional Scaling (NMDS)
ls.nmds<-metaMDS(ls[,7:10], plot=TRUE)
plot(ls.nmds, type = "t") #plotting in base R
str(ls.nmds) # gives stress value for plot
stressplot(ls.nmds) # Gain the stress plot for stress values for your MDS

## Plot using ggplot2
data.scores <- as.data.frame(scores(ls.nmds))
data.scores$site <- ls$ID
data.scores$grp <- ls$Habitat
head(data.scores)
species.scores <- as.data.frame(scores(ls.nmds, "species"))  

#Uses the scores function from vegan to extract the species (i.e., elements) scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # creates a column of species, from the rownames of species.scores
head(species.scores)

## Create hulls
grp.nearshore <- data.scores[data.scores$grp == "NS", ][chull(data.scores[data.scores$grp == 
                                                                            "NS", c("NMDS1", "NMDS2")]), ]
grp.wetland <- data.scores[data.scores$grp == "WL", ][chull(data.scores[data.scores$grp == 
                                                                          "WL", c("NMDS1", "NMDS2")]), ]

hull.data <- rbind(grp.nearshore, grp.wetland) #combine groups
hull.data

## Plot with hulls
ls.nmds.plot<-ggplot(data=data.scores, aes(x=NMDS1, y=NMDS2, shape=grp, colour=grp)) + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.30) + # add convex hulls
  geom_point(size=6) + # add point markers
  scale_colour_manual(values=c("NS" = "steelblue", "WL" = "forestgreen")) +
  scale_fill_manual(values=c("NS" = "steelblue", "WL" = "forestgreen")) +
  coord_equal() +
  theme_bw() + 
  guides(col=guide_legend("Habitat"), shape=guide_legend("Habitat"), fill=FALSE)+
  theme(axis.title.x = element_text(size=18, face="bold"), # remove x-axis labels
        axis.title.y = element_text(size=18, face="bold"),
        axis.text.x=element_text(size=16),
        axis.text.y = element_text(size=16),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.text = element_text(size=20),
        legend.title = element_text(size=24, face="bold"))

ls.nmds.plot

### PERMANOVA
df.dist<-vegdist(ls[,7:10])
df.pcoa=cmdscale(df.dist)
plot(df.pcoa[,1],df.pcoa[,2])
df.perm<-adonis(df.dist~ls$Habitat) # p=0.7
df.perm 

#### Linear Discriminant Function Analysis (LDA)
scatterplotMatrix(ls[7:10])

# Linear Discriminant Analysis with Jacknifed Prediction
# Code above performs an LDA using listwise deletion of missing data 
# CV=TRUE generates jacknifed (i.e., leave one out) predictions 
ls.fit <- lda(Habitat ~ Ba+Mg, data=ls, 
              na.action="na.omit", CV=TRUE)
ls.fit # show results
# Assess the accuracy of the prediction
# % correct for each category of Habitat
ct <- table(ls$Habitat, ls.fit$class)
diag(prop.table(ct, 1))
# total % correct
sum(diag(prop.table(ct)))

## LDA accuracy
#Single Elements
# Sr only: 20% (33% NS, 0% WL)
# Ba only: 60% (66.7% NS, 50 % WL)
# Mg only: 100% (100% NS, 100% WL)
# lnMn only: 80% (100% NS, 50% WL)

# LDA accuracy
#Element combos
# Sr+Ba+lnMn+Mg: 20% (0% NS, 50% WL)
# Sr+Ba: 40% (33% NS, 50% WL)
# Sr+Mg: 80% (100% NS, 50% WL)
# Sr+lnMn:60% (66.7% NS, 50% WL)
# Ba+lnMg: 80% (100% NS, 50% WL)
# Ba+Mg:100% (100% NS, 100% WL)


############Supplemental Figures#################
###Visualize MixSIAR Results###
#Supplemental Figure 2 - posterior estimates of NS Large Perch Fish Contributions
lp<-read_csv("MixSIAR_NS_Large_Perch_Contributions_All.csv")
ls<-filter(lp, Site=="Little Sturgeon")
ca<-filter(lp, Site=="Calumet")
bh<-filter(lp, Site=="Burns Harbor")
mu<-filter(lp, Site=="Muskegon")

#Little Sturgeon Large NS
a<-ggplot(ls, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("steelblue", "forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")

#Muskegon Large NS
b<-ggplot(mu, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("steelblue", "forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")

#Burns Harbor Large NS
c<-ggplot(bh, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("steelblue", "forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")

#Calumet Large NS
d<-ggplot(ca, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("steelblue", "forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")

#Combine
ns.lg.fig<-a+b+c+d
ns.lg.fig

#Supplemental Figure 3 - pos. estimates of NS Small Perch Fish Contributions
smns<-read_csv("MixSIAR_NS_Small_Perch_Contributions_All.csv")
ls<-filter(smns, Site=="Little Sturgeon")
bh<-filter(smns, Site=="Burns Harbor")

#Little Sturgeon Small NS
a<-ggplot(ls, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("steelblue", "forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")

#Burns Harbor Small NS
b<-ggplot(bh, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("steelblue", "forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")
#Combine
ns.sm.fig<-a+b
ns.sm.fig

# Supplemental Figure 4 - estimates of WL Large Perch Fish Contribution to Diet
wllg<-read_csv("MixSIAR_WL_Large_Perch_Contributions_All.csv")

ls<-filter(wllg, Site=="Little Sturgeon")
cr<-filter(wllg, Site=="Cedar River")
pe<-filter(wllg, Site=="Peshtigo")
pw<-filter(wllg, Site=="Pentwater")
ca<-filter(wllg, Site=="Calumet")
mu<-filter(wllg, Site=="Muskegon")

#Little Sturgeon Large WL
a<-ggplot(ls, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("steelblue", "forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")
# Cedar River Large WL
b<-ggplot(cr, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("steelblue", "forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")
# Peshtigo Large WL
c<-ggplot(pe, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("steelblue", "forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")

#Calumet Large WL
d<-ggplot(ca, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")

#Muskegon Large WL
e<-ggplot(mu, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("orange", "steelblue", "forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")
#Pentwater Large WL
f<-ggplot(pw, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("orange", "forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")
#Combine
wl.lg.fig<-a+b+c+d+e+f
wl.lg.fig

# Supplemental Figure 5 - pos. estimates of WL Small Perch Fish Contribution to Diet
wlsm<-read_csv("MixSIAR_WL_Small_Perch_Contributions_All.csv")

ls<-filter(wlsm, Site=="Little Sturgeon")
pw<-filter(wlsm, Site=="Pentwater")
ca<-filter(wlsm, Site=="Calumet")
mu<-filter(wlsm, Site=="Muskegon")

#Little Sturgeon Small WL
a<-ggplot(ls, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("steelblue", "forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")

#Calumet Small WL
b<-ggplot(ca, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")

#Muskegon Small WL
c<-ggplot(mu, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("orange", "steelblue", "forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")
#Pentwater Small WL
d<-ggplot(pw, aes(x=Resource, ymin = Min, lower = Low, middle = Mid, upper = Top, ymax = Max, fill=Habitat)) +
  geom_boxplot(stat = "identity") + xlab("Resource")+ ylab("Contribution (%)")+ theme_bw()+
  scale_fill_manual(values=c("orange", "forestgreen"))+
  theme(axis.title.x = element_text(face='bold',size=20,hjust=0.5, vjust=-1),
        axis.title.y = element_text(face='bold',size=20,vjust=1),
        axis.text.x = element_text(size=17,color='black', angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=18,color='black'),
        legend.text=element_text(size=16),
        legend.position = "none")
#Combine
wl.sm.fig<-a+b+c+d
wl.sm.fig
