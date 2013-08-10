##################################################
##  Script to summarize the SVs validated by PCR    
#  			
##  Aparicio Lab WSOP 2013-001 developed by			
##	 Dr Damian Yap , Research Associate				
##	 dyap@bccrc.ca  Version 1.0 (Aug 2013) 			
##################################################


# These commands must be specifed in order for this script to work
# source("http://www.bioconductor.org/biocLite.R"); source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); library('BSgenome.Hsapiens.UCSC.hg19')
# install.packages("XLConnect")
# source("http://www.bioconductor.org/biocLite.R"); biocLite("limma")


library(Biostrings)
library("IRanges")
library("GenomicRanges")
library(Rsamtools)
library('BSgenome.Hsapiens.UCSC.hg19')
library(limma)
library(foreign)
library(lattice)
library(XLConnect)
library(cluster) 
library(fpc)
library(gplots)
library(RColorBrewer)


#################################################
# Directory structure and file names

file="/share/lustre/backup/dyap/Projects/Tumour_Evol/Rearrangements/SV qPCR results.xls"
outdir="/share/lustre/backup/dyap/Projects/Tumour_Evol/Rearrangements"

######################################################

# workbook <- system.file(file, package = "XLConnect")
# Load workbook
wb <- loadWorkbook(file, create = FALSE)
# Query available worksheets
sheets <- getSheets(wb)

# This command reads all the sheets in a workbook into a data_list
sv <- readWorksheet(wb, sheets[1])

# This is the full set but lots of NAs for sample without full set

outdf <- data.frame( Pos = rep("", nrow(sv)),
		      ID = rep("", nrow(sv)),
		     Sample = rep("", nrow(sv)),
                     Type = rep("", nrow(sv)),
                     Tumor = rep("", nrow(sv)),
                     Normal = rep("", nrow(sv)),
                     Xeno1 = rep("", nrow(sv)),
                     Xeno2 = rep("", nrow(sv)),
                     Xeno3 = rep("", nrow(sv)),
                     Xeno4 = rep("", nrow(sv)),
                     Xeno5 = rep("", nrow(sv)),
                     Xeno6 = rep("", nrow(sv)),
                     stringsAsFactors = FALSE)
                     

for (i in seq(1:nrow(sv)))
	{

	if (is.na(sv$Col2)[i]) next

	outdf$ID[i] <- strsplit(sv$Col2[i], split="_")[[1]][1]
	outdf$Sample[i] <- strsplit(strsplit(sv$Col2[i], split="_")[[1]][2], split="-")[[1]][1]
	outdf$Type[i] <- strsplit(strsplit(sv$Col2[i], split="_")[[1]][2], split="-")[[1]][2]

	outdf$Pos[i] <- sv$Col1[i]



if (sv$T[i] == "-") outdf$Tumor[i] <- 0
if (sv$T[i] > 0) outdf$Tumor[i] <- 1
if (sv$T[i] == "n/a") outdf$Tumor[i] <-"N/A"

if (sv$N[i] == "-") outdf$Normal[i] <- 0
if (sv$N[i] > 0) outdf$Normal[i] <- 1
if (sv$N[i] == "n/a") outdf$Normal[i] <- "N/A"


if (sv$X1[i] == "-") outdf$Xeno1[i] <- 0
if (sv$X1[i] > 0) outdf$Xeno1[i] <- 1
if (sv$X1[i] == "n/a") outdf$Xeno1[i] <- "N/A"

if (sv$X2[i] == "-") outdf$Xeno2[i] <- 0
if (sv$X2[i] > 0) outdf$Xeno2[i] <- 1
if (sv$X2[i] == "n/a") outdf$Xeno2[i] <- "N/A"

if (sv$X3[i] == "-") outdf$Xeno3[i] <- 0
if (sv$X3[i] > 0) outdf$Xeno3[i] <- 1
if (sv$X3[i] == "n/a") outdf$Xeno3[i] <- "N/A"

if (sv$X4[i] == "-") outdf$Xeno4[i] <- 0
if (sv$X4[i] > 0) outdf$Xeno4[i] <- 1
if (sv$X4[i] == "n/a") outdf$Xeno4[i] <- "N/A"

if (sv$X5[i] == "-") outdf$Xeno5[i] <- 0
if (sv$X5[i] > 0) outdf$Xeno5[i] <- 1
if (sv$X5[i] == "n/a") outdf$Xeno5[i] <- "N/A"

if (sv$X6[i] == "-") outdf$Xeno6[i] <- 0
if (sv$X6[i] > 0) outdf$Xeno6[i] <- 1
if (sv$X6[i] == "n/a") outdf$Xeno6[i] <- "N/A"

	}

# Summary of data
table(outdf$Sample)
table(outdf$Type)


###################
# Assign colour to all types of rearrangements

ty<-names(table(outdf$Type))

rowcol<-brewer.pal(length(ty), "Accent")

Colassign<-data.frame(Type = rep("", length(ty)),
		      Col = rep("", length(ty)),
                     stringsAsFactors = FALSE)

for ( i in seq(length(ty)))
	{
	Colassign$Type[i] <- ty[i]
	Colassign$Col[i] <- rowcol[i]
	}	 

Colassign$Type[1]<-"No data"
# white =0 (no PCR product)
# black = any PCR product

hmcols<-colorRampPalette(c("white","black"))(2)


#############################
# CHANGE THIS FOR EACH SAMPLE
# subsetting data

sample_list=c("SA493","SA494","SA495","SA499","SA500")

for (sample in sample_list)
{
xaxislab=paste(sample, "Sample", sep=" ")
pdffile=paste(paste(outdir,sample,sep="/"),"SV.pdf",sep="_")

sub <- outdf[ which(outdf$Sample==sample), ]

# Summary of the data
# table(SA499$Type)

data <- data.matrix(sub[,5:12])


# Remove columns with all NAs
mat <- data[,colSums(is.na(data))<nrow(data)]

row.names(mat)<-make.unique(sub$ID)

rearr <- data.frame(Type = rep("", nrow(sub)),
		      Col = rep("", nrow(sub)),
                     stringsAsFactors = FALSE)

row.names(rearr)<-make.unique(sub$ID)

# Get the colors for different types of rearrangements
for (j in seq(nrow(sub)))
	{
	match <- sub$Type[j]
	rearr$Type[j]<- match
	rearr$Col[j] <- Colassign[Colassign[,1] %in% c(match),2]
	}

title="PCR validated Structural Variants"

pdf(pdffile, width=7, height=8)

heatmap.2(mat, main=title, xlab=sample, ylab="Rearrangement ID", scale="none", key = FALSE, Colv=FALSE
, cexCol=0.8, cexRow=0.6, col = hmcols, RowSideColors=rearr$Col, trace="none")

legend("top",legend=Colassign$Type, fill=Colassign$Col, border=TRUE, bty="o", y.intersp = 0.7, cex=0.7)

dev.off()

}
################
