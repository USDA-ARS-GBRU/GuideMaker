#!/usr/bin/env Rscript
options(warn=-1)
args = commandArgs(trailingOnly=TRUE)




suppressWarnings(suppressMessages(library(karyoploteR)))
suppressWarnings(suppressMessages(library(regioneR)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(base)))
suppressWarnings(suppressMessages(library(zoo)))
suppressWarnings(suppressMessages(library(rtracklayer)))
suppressWarnings(suppressMessages(library(GenomeInfoDb)))



#library(IRanges)
#library(regioneR)
##library(GenomicRanges)
#library(stats4)
#library(BiocGenerics)
#library(parallel)
library(karyoploteR)
library(regioneR)
library(GenomeInfoDbData)
library(zoo)
library(dplyr)

# source("http://bioconductor.org/biocLite.R")
# biocLite("devtools")
# biocLite("Bioconductor/GenomeInfoDb")




###########
data.points.colors <- c("#FFBD07", "#00A6ED",  "#FF3366", "#8EA604", "#C200FB")


#args[1] # args.outdir

filepath = paste0(args[1],"/","targets.csv")
filepath = paste0("test","/","targets.csv")

data <- read.csv(filepath, header = TRUE) %>%
  dplyr::select(Accession,Feature.start,Feature.end,Feature.strand,Guide.start,Guide.end,Guide.strand,PAM)
  
# number of chromosomes
chrom_len <- data %>%
  dplyr::select(Accession, Feature.end) %>%
  dplyr::group_by(Accession) %>%
  dplyr::summarise(max = max(Feature.end))


WINDOW_SIZE = max(chrom_len$max)/ 100

WINDOW_SIZE = round(max(chrom_len$max)/ 100) ### divide into 100 bins
# Save plot

filename = paste0(args[1],"/","GuideMakerPlot.pdf")

pdf(file=filename)

y <- Seqinfo(seqnames=chrom_len$Accession ,
             seqlengths=chrom_len$max)


pp <- getDefaultPlotParams(plot.type=2)

#Change the ideogramheight param to create thicker ideograms 
pp$ideogramheight <- 15


# https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotCoverage/PlotCoverage.html
suppressWarnings(suppressMessages(library(regioneR)))
suppressWarnings(warning("filterChromosomes"))
suppressWarnings(library(karyoploteR))

kp <- plotKaryotype(plot.type = 2, y, main="GuideMaker Plot", plot.params=pp, cex = 0.5)
kpAddBaseNumbers(kp, tick.dist = round(max(chrom_len$max)/5), add.units = TRUE)

#kpAddLabels(kp, labels = "Trait 1", srt=90, pos=3, cex=1.8, label.margin = 0.025)
#kpText(kp, chr=seqlevels(kp$genome), y=1, x=1,  r0=0.2, r1=0, col="#444444", label="PAM", cex=0.8, pos=2)

# density
df_density <- data.frame(chr=data$Accession, start=data$Guide.start, end =data$Guide.end,strand=data$Guide.strand)
df_density = makeGRangesFromDataFrame(df_density)  #
#kp <- kpPlotDensity(kp, df_density)
kp = kpPlotDensity(kp, data=df_density, r0=0, r1=0.3, window.size = WINDOW_SIZE, border="blue",col="orchid")
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.5, cex=0.5)
kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.5)

WINDOW_SIZE_kb = WINDOW_SIZE/1000
WINDOW_SIZE_label = paste0("Density at ",WINDOW_SIZE_kb, "KB")

kpAddLabels(kp, labels = WINDOW_SIZE_label, srt=90, pos=1, cex=0.5, label.margin = 0.01,side="right", data.panel=1, r0=0.0, r1=0.5)
# datapoints

df <- data.frame(chr=data$Accession, start=data$Guide.start, end =data$Guide.start,strand=data$Guide.strand)
data.points = makeGRangesFromDataFrame(df)  #
mcols(data.points) <- data.frame(y=rnorm(n = dim(df)[1], 0.5, sd = 0.1))
data$PAM = as.factor(data$PAM)
data["pam_color"] = data.points.colors[data$PAM]


point_size = 80000 / (max(chrom_len$max))
point_size = ifelse(point_size < 0.08, 0.08, point_size)
#Data points
#kpAxis(kp, ymin = 0, ymax = 1, r0=0, r1=0.8, numticks = 2, col="#666666", cex=0.5,labels="PAM",label.margin = 5)
kpPoints(kp, data=data.points, pch=16, cex=point_size, col=data$pam_color, r0=0.6, r1=1)
kpAddLabels(kp, labels = "PAM", srt=90, pos=3, cex=0.5, label.margin = 0.025, r0=0.6, r1=1)



###plot gene as region
df <- data.frame(chr=data$Accession, start=data$Feature.start, end=data$Feature.end, strand=data$Feature.strand)
df <- unique(df)
kpPlotRegions(kp, data=df, col="#FFEECC", border="#FFCCAA", r0=0.0, r1=0.2, data.panel=2)
kpAddLabels(kp, labels = "Features", srt=90, pos=3, cex=0.5, label.margin = 0.025, data.panel=2, r0=0.0, r1=0.2)

invisible(dev.off())
