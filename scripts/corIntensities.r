##---------------------------------------------------------------------#
##
## Title: epicV2 and epicV2 correlations
##
## Purpose: Take samples that have been run on both array type and:
##
## 1. plot/correlate intensity values
## 2. plot distributions of raw and normalised beta values
##
## Author: Emma Walker
##
## Date Created: 28/02/2024
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#


#----------------------------------------------------------------------#
# Set up
#----------------------------------------------------------------------#

library(dplyr)
library(bigmelon)
library(ggplot2)
library(reshape2)
library(minfi)

#dataDir <- ""

#setwd(dataDir)
dataDir<-getwd()

v1gdsFile <- file.path(dataDir, "data/epicV1/2_gds/raw.gds")
v2gdsFile <- file.path(dataDir, "data/epicV2/2_gds/raw.gds")

v1norm <- file.path(dataDir, "data/epicV1/3_normalised/normalised.rdata")
v2norm <- file.path(dataDir, "data/epicV2/2_gds/rawNorm.gds")

# load v1 raw betas
v1gfile <- openfn.gds(v1gdsFile, readonly = FALSE)
v1betas <- betas(v1gfile)[,]
closefn.gds(v1gfile)

# load v1 norm betas
load(v1norm)
v1normBetas <- normbeta
rm(normbeta)


# load v2 raw betas
v2gfile <- openfn.gds(v2gdsFile, readonly = FALSE)
v2betas <- betas(v2gfile)[,]
closefn.gds(v2gfile)

# load v2 norm betas
v2gfileNorm <- openfn.gds(v2norm, readonly = FALSE)
v2normbetas <- index.gdsn(v2gfileNorm, "normbeta")[,]
closefn.gds(v2gfileNorm)


# load QC data

load("data/epicV1/2_gds/QCmetrics/QCmetrics.rdata")
v1 <- QCmetrics

load("data/epicV2/2_gds/QCmetrics/QCmetrics.rdata")
v2 <- QCmetrics


# load samplesheets to get IDs for matching

samp1 <- read.csv("data/epicV1/0_metadata/EPIC_v1_v2.csv", stringsAsFactors = F)
colnames(samp1)[1] <- "Individual_ID"
samp1$newID <- paste0(samp1$Individual_ID, "_", samp1$Cell.Fraction)
samp1$Basename <- paste0(samp1$Chip.ID, "_", samp1$Chip.Location)


samp2 <- read.csv("data/epicV2/0_metadata/v2.csv", stringsAsFactors = F)[,1:6]
colnames(samp2)[1] <- "Individual_ID"
samp2$newID <- paste0(samp2$Individual_ID, "_", samp2$Cell.Fraction)
samp2$Basename <- paste0(samp2$Chip.ID, "_", samp2$Chip.Location)

# remove samples not in both datasets
samp2 <- samp2[!samp2$newID == setdiff(samp2$newID, samp1$newID),]

# add qcmetrics data to samplesheets and match IDs
s1 <- left_join(v1, samp1 %>% dplyr::select(Basename, newID))
s2 <- left_join(v2, samp2 %>% dplyr::select(Basename, newID))

s1 <- s1[match(s2$newID, s1$newID),]


# load pass/fail status and add to samplesheets

passV1 <- read.csv("data/epicV1/2_gds/QCmetrics/passQCStatusStage3AllSamples.csv", stringsAsFactors = F)
passV2 <- read.csv("data/epicV2/2_gds/QCmetrics/passQCStatusStage3AllSamples.csv", stringsAsFactors = F)

s1 <- left_join(s1, passV1 %>% dplyr::select(Basename, passQCS3))
s2 <- left_join(s2, passV2 %>% dplyr::select(Basename, passQCS3))

# find samples failed in either
fail <- rep(FALSE, nrow(s1))
fail[which(s1$passQCS3 == F | s2$passQCS3 == F)] <- TRUE


#----------------------------------------------------------------------#
# Intensity values
#----------------------------------------------------------------------#

# all samples

plotdf <- as.data.frame(cbind(s1$M.median, s1$U.median, s2$M.median, s2$U.median))
colnames(plotdf) <- c("v1M.median", "v1U.median", "v2M.median", "v2U.median")
plotdf <- mutate_all(plotdf, function(x) as.numeric(as.character(x)))
plotdf <- cbind(s1$newID, plotdf, fail)



# passed
  
  # s1[passV1$Basename[!passV1$passQCS3],] - code to get passed samples

#----------------------------------------------------------------------#
# Intensity boxplot
#----------------------------------------------------------------------#

meltdf <- reshape2::melt(plotdf)

ggplot(meltdf, aes(x=variable, y=value)) + 
         geom_boxplot()+
  ylab("intensity")


t.test(plotdf$v1M.median, plotdf$v2M.median, paired=TRUE)
# p-value = 0.2131
t.test(plotdf$v1U.median, plotdf$v2U.median, paired=TRUE)
# p-value = 0.4057


#----------------------------------------------------------------------#
# Intensity correlation plots
#----------------------------------------------------------------------#

# colour by cell type
plotdf$cellType <- gsub(".*_", "", plotdf$`s1$newID`)

# M intens
ggplot(plotdf, aes(x=v1M.median, y=v2M.median, colour=cellType))+
  geom_point()

cor.test(plotdf$v1M.median, plotdf$v2M.median)
# p-value = 0.01515
# r = 0.3


# U intens
ggplot(plotdf, aes(x=v1U.median, y=v2U.median, colour=cellType))+
  geom_point()

cor.test(plotdf$v1U.median, plotdf$v2U.median)
# p-value = 0.0009215
# r = 0.4



# colour by pass/fail

# M intens
ggplot(plotdf, aes(x=v1M.median, y=v2M.median, colour=fail))+
  geom_point()

# U intens
ggplot(plotdf, aes(x=v1U.median, y=v2U.median, colour=fail))+
  geom_point()


cor.test(plotdf$v1M.median[plotdf$fail == FALSE],
         plotdf$v2M.median[plotdf$fail == FALSE])


#----------------------------------------------------------------------#
# raw betas density plots 
#----------------------------------------------------------------------#

v1betas <- v1betas[,s1$Basename]

densityPlot(v1betas, sampGroups=s2$Cell_Type, legend=FALSE)
densityPlot(v2betas, sampGroups=s2$Cell_Type, legend=FALSE)


#----------------------------------------------------------------------#
# normalised betas density plots 
#----------------------------------------------------------------------#

## leave for now and do violin plots instead

# cbind betas vals
# rbind sample sheets (basenames, IDs, tech type)
# make sure colnames=basenames
# plot with samp groups as tech type


#----------------------------------------------------------------------#
# violin plots normalised betas
#----------------------------------------------------------------------#

## v1

# subset to QC passed samples and wrangle for plotting
v1Pass <- v1[v1$Basename %in% colnames(v1normBetas),]

v1normPlotdf <- as.data.frame(t(v1normBetas))
v1normPlotdf$Basename <- rownames(v1normPlotdf)

v1normPlotdf$cellType <- v1Pass$Cell_Type
v1normMelt<- reshape2::melt(v1normPlotdf)

# plot
ggplot(v1normMelt, aes(x=cellType, y=value)) + 
  geom_violin()


## v2

# subset to QC passed samples and wrangle for plotting
v2Pass <- v2[v2$Basename %in% colnames(v2normbetas),]

v2normPlotdf <- as.data.frame(t(v2normbetas))
v2normPlotdf$Basename <- rownames(v2normPlotdf)

v2normPlotdf$cellType <- v2Pass$Cell_Type
v2normMelt<- reshape2::melt(v2normPlotdf)

#plot
ggplot(v2normMelt, aes(x=cellType, y=value)) + 
  geom_violin()

