# dont have any singletons because remove low read taxa (threshold=10), correct? may mess up diversity analysis
# cant do fisher
# can use prevelence filtering to reduce complexity and contam without effecting downstream analysis, https://www.frontiersin.org/articles/10.3389/fmicb.2020.607325/full 
# use different normalisation processes for beta vs DAo: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13115

### ======IMPORT AND FILTER================================ ----
### packages ----
library(plyr)
library(metagenomeSeq)
library(DESeq2)
library(phyloseq)
library(vegan)
library(microViz)
library(strex)
library(edgeR)

library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(microbiome)
library(decontam)

Packages <- c("tidyverse","readxl","cowplot","ggrepel","ggpubr", # generic
              "reshape2", # for melting in HM
              "DescTools", # for overlap with telo/centro
              "dplyr", # for kw test
              "GenomicRanges", "rtracklayer", # for hg conversions
              "msigdbr", # for hallmark gene lists
              "SciViews", "vegan", # for shannon calculations
              "pvclust", # clustering
              "fastDummies","betareg","car","frmselection", # building model
              "caret", #LOOCV
              "BSgenome", # for hg19 coordinate in chromomap
              "chromoMap", # visualising bin location
              "org.Hs.eg.db", # swtiching between gene IDS
              "limma","rrvgo", # GO
              "survival","survminer", # for survival
              "TCGAbiolinks", # pull tcga rna data
              "TCGAutils", # for TCGA UUID barcode translator
              "DESeq2","powerjoin", # for RNA expression
              "ape"
)

lapply(Packages, library, character.only = TRUE)

### functions ----
`permanova_pairwise` <- function(x,
                                 grp,
                                 permutations = 999,
                                 method = 'bray',
                                 padj = 'bonferroni', ...) {
  f     <- as.factor(grp)
  if (!all(table(f) > 1)) warning('factor has singletons! perhaps lump them?')
  co    <- combn(unique(as.character(f)),2)
  nco   <- NCOL(co)
  out   <- data.frame(matrix(NA, nrow=nco, ncol=5))
  dimnames(out)[[2]] <- c('pairs', 'SumOfSqs', 'F.Model', 'R2', 'pval')
  if (!inherits(x, 'dist')) {
    D <- vegan::vegdist(x, method=method)
  } else {
    D <- x
  }
  cat('Now performing', nco, 'pairwise comparisons. Percent progress:\n')
  for(j in 1:nco) {
    cat(round(j/nco*100,0),'...  ')
    ij  <- f %in% c(co[1,j],co[2,j])
    Dij <- as.dist(as.matrix(D)[ij,ij])
    fij <- data.frame(fij = f[ij])
    a   <- vegan::adonis2(Dij ~ fij, data=fij, permutations = permutations,
                          ...)
    out[j,1] <- paste(co[1,j], 'vs', co[2,j])
    out[j,2] <- a$SumOfSqs[1]
    out[j,3] <- a$F[1]
    out[j,4] <- a$R2[1]
    out[j,5] <- a$`Pr(>F)`[1]
  }
  cat('\n')
  out$p.adj <- p.adjust(out$pval, method=padj)
  attr(out, 'p.adjust.method') <- padj
  cat('\np-adjust method:', padj, '\n\n')
  return(out)
}
GeomSplitViolin <- ggplot2::ggproto("GeomSplitViolin", ggplot2::GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x'])
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- ggplot2::GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", ggplot2::GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

theme_custom <- function(){
  theme_classic() %+replace%  
    theme(
      plot.margin=margin(t=0.2,r=0.5,b=0.5,l=0.5,"cm"),
      panel.background = element_blank(),
      plot.title = element_text(size=28, colour='black',face='bold', hjust=0),
      axis.line = element_line(size = 0.5, colour = "black"),
      axis.title = element_text(size=28, colour='black'),
      axis.text = element_text(size=28, colour='black'),
      axis.ticks.length=unit(0.2, "cm")
    )
}

'%!in%' <- function(x,y)!('%in%'(x,y))

# Function to pull basic info from rawdata eg number of samples/bins/patients
# Requires: chr | start | stop | sample1CN | sample2CN
# Requires sample name in colname, in structure: patientID.sampleNumber
PullDataInfo <- function(rawdata) {
  
  # Dataframe identifying start and stop codon, and chromosome for each bin
  start.stop <- rawdata[,c(1:3)]
  start.stop$bin <- 1:nrow(start.stop)
  start.stop <- start.stop[, c(4, 1, 2, 3)]
  
  # Vector holding sample ID
  sampleIDs <- colnames(rawdata)[-c(1:3)] 
  
  # Vector holding patient identifiers, ie the string before the '.' in sample IDs
  patientIDs <- unique(sub("\\..*", "", colnames(rawdata)))[-c(1:3)] 
  
  # Number of samples
  noSamples <- length(sampleIDs)
  
  # Number of patients
  noPatients <- length(patientIDs)
  
  # list of number of samples per patient
  sampPerPatient <- list()
  for ( i in 1:noPatients ) {
    wd <- rawdata[, which(sub("\\..*", "", colnames(rawdata))==patientIDs[i])]
    sampPerPatient[[i]] <- ncol(wd)
  }
  
  # Number of bins
  noBins <- length(start.stop$bin)
  
  # Visualisation may require a vector to identify of chromosome ends and chromosome midpoints 
  chr.ends <- cumsum(table(start.stop$chr))
  list <- list()
  l <- 1
  for ( i in 1:length(chr.ends) ) {  
    if ( i == 1 ) { 
      list[[l]] <- chr.ends[i]/2 
      l <- l+1
    }
    else { 
      list[[l]] <- chr.ends[i-1] + ((chr.ends[i]-chr.ends[i-1])/2)
      l <- l+1
    }
  }
  chr.mid <- unlist(list)
  chr.ends <- data.frame(start=c(0,chr.ends[-22]), end=(chr.ends), col=c('G','W'))
  
  #average number of bins per patient
  binsPerPatient <- list()
  for (i in 1:length(patientIDs)) {
    wd <- data.frame(rawdata[,which(sub("\\..*", "", colnames(rawdata))==patientIDs[i])])
    binsPerPatient[[i]] <- round(mean(colSums(!is.na(wd))))
  }
  
  
  # Return data
  newData <- list('start.stop'=start.stop, 'sampleIDs'=sampleIDs, 'patientIDs'=patientIDs, 'noSamples'=noSamples, 'noPatients'=noPatients,
                  'sampPerPatient'=sampPerPatient, 'noBins'=noBins,'chr.mid'=chr.mid,  'chr.end'=chr.ends, 'binsPerPatient'=binsPerPatient)
  return(newData)
}

# Function to create 5 new data structures relating to clonality
# Requires: chr | start | stop | sample1CN | sample2CN
# dataInfo must be previously generated from rawdata using PullDataInfo
PullDataClonality <- function(rawdata, dataInfo) {
  
  # A dataframe with a bin per row and a patient per column, with values indicating clonality. 
  # 0=notCNA, 1=subclonalCNA, 2=clonalCNA
  clonal.data <- rawdata[,-c(1:3)]
  l <- 1
  clonal <- list()
  for ( k in 1:nrow(clonal.data) ) {
    for ( i in 1:length(dataInfo$patientIDs) ) {
      wd <- clonal.data[k, which(sub("\\..*", "", colnames(clonal.data))==dataInfo$patientIDs[i])]
      wd <- wd[, !is.na(wd)]
      
      if ( 1 %in% wd | 3 %in% wd ) { #if one of the samples has a mutation, proceed
        if ( length(unique(t(wd)))==1 ) { #if all the same then clonal
          clonal[[l]] <- 2
          l <- l + 1
        }
        else { #different = subclonal
          clonal[[l]] <- 1
          l <- l + 1
        }
      }
      else if ( length(wd)==0 ) { #all MR for that bit are NA
        clonal[[l]] <- NA
        l <- l + 1
      }
      else { #neither sample has a mutation
        clonal[[l]] <- 0 
        l <- l + 1
      }
    }
  }
  
  clonal.data <- data.frame(t(matrix(unlist(clonal), ncol=dataInfo$noBins)))
  colnames(clonal.data) <- dataInfo$patientIDs
  clonal.data[] <- lapply(clonal.data, factor, levels=unique(unlist(clonal.data)))
  
  
  # Dataframe detailing the counts of gains/losses and whether they are subclonal or clonal
  CNA.clo.counts <- data.frame(bin = dataInfo$start.stop$bin, chr1 = NA, chr2 = NA,
                               clonal.aneu = NA, subclonal.aneu = NA, gain = NA, loss = NA,
                               clonal.gain = NA, clonal.loss = NA, clonal.noCNA = NA, 
                               subclonal.gain = NA, subclonal.loss = NA)
  
  CNA.clo.counts$chr1 <- as.numeric(dataInfo$start.stop[match(CNA.clo.counts$bin, dataInfo$start.stop$bin), 2]) 
  CNA.clo.counts$chr2 <- as.numeric(dataInfo$start.stop[match(CNA.clo.counts$bin, dataInfo$start.stop$bin), 2]) 
  CNA.clo.counts$chr2 <- factor(CNA.clo.counts$chr2,levels = rev(seq(1:22))) 
  
  data <- rawdata[,-c(1:3)]
  for ( k in 1:nrow(data) ) { # for a bin
    clonal.all <- clonal.gain <- clonal.loss <- clonal.noCNA <- subclonal.all <- subclonal.gain <- subclonal.loss <- subclonal.noCNA <- 0
    
    for ( i in 1:dataInfo$noPatients ) { # for a patient
      
      # If its NA
      if ( is.na(clonal.data[k,i]) == TRUE ) {
        next
      }
      
      # If its clonal
      if ( clonal.data[k,i] == 2 ) {
        wd <- data[k, which(sub("\\..*", "", colnames(data))==dataInfo$patientIDs[i])]
        wd <- wd[, is.na(wd)!=TRUE]
        
        if ( 3 %in% wd ) { # if both are gains
          clonal.gain <- clonal.gain + 1
        }
        else if ( 1 %in% wd ) { # if both are losses
          clonal.loss <- clonal.loss + 1
        }
      }
      
      # if its subclonal
      else if ( clonal.data[k,i] == 1 ) {
        wd <- data[k, which(sub("\\..*", "", colnames(data))==dataInfo$patientIDs[i])]
        wd <- wd[, is.na(wd)!=TRUE]
        
        if ( 3 %in% wd ) { # if one is a gain
          subclonal.gain <- subclonal.gain + 1
        }
        if ( 1 %in% wd ) { # if one is a loss
          subclonal.loss <- subclonal.loss + 1
        }
      }
      
      # if its no CNA
      else if ( clonal.data[k,i] == 0 ) {
        clonal.noCNA <- clonal.noCNA + 1
      }
    }
    
    CNA.clo.counts$clonal.gain[k] <- clonal.gain
    CNA.clo.counts$clonal.loss[k] <- clonal.loss
    CNA.clo.counts$clonal.noCNA[k] <- clonal.noCNA
    CNA.clo.counts$subclonal.gain[k] <- subclonal.gain
    CNA.clo.counts$subclonal.loss[k] <- subclonal.loss
  }
  
  CNA.clo.counts$clonal.aneu <- CNA.clo.counts$clonal.gain + CNA.clo.counts$clonal.loss
  CNA.clo.counts$subclonal.aneu <- CNA.clo.counts$subclonal.gain + CNA.clo.counts$subclonal.loss
  CNA.clo.counts$gain <- CNA.clo.counts$clonal.gain + CNA.clo.counts$subclonal.gain
  CNA.clo.counts$loss <- CNA.clo.counts$clonal.loss + CNA.clo.counts$subclonal.loss
  
  # dataframe showing: bin | countGain/NoPatient | countLoss/NoPatient | 
  CloFreq <- cbind(bin=CNA.clo.counts$bin, gain=CNA.clo.counts$gain/dataInfo$noPatients, loss=CNA.clo.counts$loss/dataInfo$noPatients)
  
  # dataframe showing what percent of gain and loss are subclonal. On a patient basis: bin | gain | loss
  pcSubclonal <- data.frame(bin=1:dataInfo$noBins, gain=CNA.clo.counts$subclonal.gain / CNA.clo.counts$gain, loss=CNA.clo.counts$subclonal.loss / CNA.clo.counts$loss)
  
  # Count of noCNA, subclonalCNA, and clonalCNA by patient
  patientClo <- as.data.frame(t(sapply(clonal.data, table))) 
  patientClo <- patientClo[, c('0','1','2')]
  colnames(patientClo) <- c('noCNA','subclonal','clonal')
  patientClo$CNA <- patientClo$subclonal + patientClo$clonal
  patientClo$patient <- rownames(patientClo)
  
  # Return data
  newData <- list('clonal.data'=clonal.data, 'CNA.clo.counts'=CNA.clo.counts, 
                  'CloFreq'=CloFreq, 'pcSubclonal'=pcSubclonal, 'patientClo'=patientClo)
  return(newData)
}

# Function to calc pic score per patient across multiple samples
PIC <- function(data, sample.no, index) { 
  # based on table of unique patient nos
  # returns PIC per row (bin), calc across cols as given by index
  # PIC formula: 1- ((CN1/n)^2 + (CN2/n)^2 + (CN3/n)^2), where CN1 is no. of counts of copy number1
  PIC <- 1 - ((rowSums(data[,index]==1, na.rm = TRUE)/sample.no)^2 + 
                (rowSums(data[,index]==2, na.rm = TRUE)/sample.no)^2 + 
                (rowSums(data[,index]==3, na.rm = TRUE)/sample.no)^2)
  return(PIC)
}

# Function to create 4 new data structures relating to diversity
# Requires: chr | start | stop | sample1CN | sample2CN
# There must only be 2 multiregion samples per patient
# dataInfo must be previously generated from rawdata using PullDataInfo
PullDataDiversity <- function(rawdata, dataInfo) {
  
  # Cannot use averaged raw CN across the MR samples in case of variable number of samples per patient.
  # Will therefore can use a PIC.frac for each patient. 
  # This is a continuous, not binary, measure of clonality
  # Define the maximum pic score depending on the number of samples (up to 13)
  max.pics <- list()
  for ( d in 1:50 ) {
    if ( d %% 3 == 0 ) {
      n1 <- n2 <- n3 <- d/3
    }
    else if ( d %% 3 == 1 ) {
      n1 <- ((d-1)/3) + 1
      n2 <- ((d-1)/3)
      n3 <- ((d-1)/3)
    }
    else if ( d %% 3 == 2 ) {
      n1 <- ((d-2)/3) + 1
      n2 <- ((d-2)/3) + 1
      n3 <- ((d-2)/3)
    }
    max.pics[[d]] <- pic.score <- 1 - ( (n1/d)^2 + (n2/d)^2 + (n3/d)^2 )
  }
  
  # A dataframe with a bin per row and a patient per column, with values indicating pic score
  # And a dataframe showing pic.frac per patient
  pic <- list()
  pic.frac <- list()
  for ( i in 1:length(dataInfo$patientIDs)) {
    # Set working data as the cols holding the samples
    wd <- rawdata[, which(sub("\\..*", "", colnames(rawdata))==dataInfo$patientIDs[i])]
    
    # Record the number of samples
    upto <- ncol(wd)
    
    # Use PIC function on wd
    pic[[i]] <- PIC(wd, upto, c(1:upto))
    pic[[i]] <- na.omit(pic[[i]])
    
    # Define the max possible diversity given the number of sample, as maxPIC*number of bins
    # max.ITH <- max.pics[[upto]] * length(pic[[i]])
    max.ITH <- max.pics[[upto]] * dataInfo$binsPerPatient[[i]]
    
    # Store as a dataframe
    pic.frac[[i]] <- data.frame(pic.frac = sum(pic[[i]], na.rm = TRUE)/max.ITH)
  }
  
  pic.data <- as.data.frame(do.call('cbind',pic))
  colnames(pic.data) <- dataInfo$patientIDs
  pic.frac <- do.call('rbind',pic.frac)
  
  # Calulate average PIC per bin as a measure of average bin hetero across patients 
  ave.pic <- dataInfo$start.stop[,c(1:2)]
  ave.pic$avePic <- rowSums(pic.data)/(ncol(pic.data))
  
  # A dataframe of the propotion of genome gained/lost per sample, alongside ith
  pga <- data.frame(t(rawdata[,-c(1:3)]), check.names = FALSE)
  pga <- data.frame(prop.gain=apply(pga,1,function(x) sum(x == 3, na.rm = TRUE)/ncol(pga)),
                    prop.loss=apply(pga,1,function(x) sum(x == 1, na.rm = TRUE)/ncol(pga)))
  pga$prop.aneu <- pga$prop.gain + pga$prop.loss
  pga <- cbind(as.data.frame(lapply(pic.frac, rep, dataInfo$sampPerPatient)),
               pga)
  
  # Return data
  newData <- list('pic.data'=pic.data, 'pic.frac'=pic.frac, 'ave.pic'=ave.pic, 'pga'=pga)
  return(newData)
}

# Function to plot a split violin plot
GeomSplitViolin <- ggplot2::ggproto("GeomSplitViolin", ggplot2::GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x'])
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- ggplot2::GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", ggplot2::GeomPolygon$draw_panel(newdata, ...))
  }
})


geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

### choose confidence level in kraken2 ----

# import biom, convert to phylo, fix data, add metadata, only bacteria
list_of_files <- list.files(path = "~/Documents/Microbial/EPICCbraken.conf",
                            recursive = TRUE,
                            pattern = "\\.biom$",
                            full.names = TRUE)

conf.phylos <- list()
for (i in 1:length(list_of_files)) {
  x <- import_biom(list_of_files[i], parseFunction=parse_taxonomy_greengenes)
  x <- phyloseq_validate(x, remove_undetected = TRUE)
  conf.phylos[[i]] <- tax_fix(x)
  
  conf.phylos[[i]] <- conf.phylos[[i]] %>%
    ps_mutate(
      Id = substr(Id, 7, nchar(Id)-23),
      patient = str_before_nth(Id,"_",1),
      region = substr(Id,1,6),
      subregion = str_before_nth(Id,"_",2),
      sample = str_before_nth(Id,"_",3),
      type = substr(str_after_nth(Id,"_",2),1,1),
      group = ifelse( (substr(Id,6,6) %in% c("F","G","H")), "adenoma", 
                      ifelse( (substr(Id,6,6) %in% c("E")), "normal",
                              ifelse( (substr(Id,6,6) %in% c("A","B")), "carcinoma", 
                                      ifelse( (substr(Id,6,6) %in% c("C","D") & patient == "C516"), "adenoma",
                                              ifelse( (substr(Id,6,6) %in% c("C","D") & patient != "C516"), "carcinoma", 
                                                      ifelse( (substr(Id,6,6) %in% c("Z")), "blood", 
                                                              ifelse( (substr(Id,6,6) %in% c("W")), "normalInCancer", NA))))))),
      
    )
  
  conf.phylos[[i]] <- subset_taxa(conf.phylos[[i]], Kingdom=="Bacteria")
}

# quantify richness (n species) and evenness (shannon)
rich <- list()
even <- list()
for (i in 1:length(conf.phylos)) {
  print(i)
  rich[[i]] <- estimate_richness(conf.phylos[[i]], measures = c("Chao1"))[1]
  rich[[i]]$sample <- rownames(rich[[i]])
  colnames(rich[[i]]) <- c(str_before_first(str_after_last(list_of_files[i],"/"), "_"),"sample")
  rich[[i]] <- rich[[i]][c(2,1)]
  even[[i]] <- evenness(conf.phylos[[i]], 'pielou')
  even[[i]]$sample <- rownames(even[[i]])
  colnames(even[[i]]) <- c(str_before_first(str_after_last(list_of_files[i],"/"), "_"),"sample")
  even[[i]] <- even[[i]][c(2,1)]
}

rich <- reduce(rich, full_join)
rich <- gather(rich, key="key", value="value", -sample)
rich$key <- factor(rich$key, levels = c(seq(0,1,0.1)))
rich$sample <- substr(rich$sample, 7, 16)

even <- reduce(even, full_join)
even <- gather(even, key="key", value="value", -sample)
even$key <- factor(even$key, levels = c(seq(0,1,0.1)))

jpeg('tempfig.jpeg', width = (40), height = (20), units = "cm", res = 300)
ggplot(rich, aes(x=key, y=value, group=sample, colour=sample)) +
  geom_line(size=1, alpha=0.8) +
  geom_point(size=4, alpha=1) +
  scale_color_manual(values =  c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                       "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6",
                       "#6A3D9A", "#FFFF99", "#B15928","#FF9933")) +
  ylab("Chao1 richness") +
  xlab("Kraken2 confidence level") +
  theme_custom() +
  theme(legend.text = element_text(size=20), legend.title = element_blank(),
        legend.position = "bottom")
dev.off()

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(even, aes(x=key, y=value, group=sample, color=sample)) +
  geom_line(size=1, alpha=0.8) +
  geom_point(size=4, alpha=1) +
  scale_color_manual(values =  c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                                 "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6",
                                 "#6A3D9A", "#FFFF99", "#B15928","#FF9933")) +
  ylab("Pielou eveness") +
  xlab("Kraken2 confidence level") +
  theme_custom() +
  theme(legend.text = element_text(size=24), legend.title = element_text(size=24),
        legend.position = "none")
dev.off()


### data import and fix ----

    # import biom data
#data <- import_biom("~/Documents/Microbial/bracken_all_26.6.23.biom", parseFunction=parse_taxonomy_greengenes)
data <- import_biom("~/Documents/Microbial/bracken_all_conf0.1.biom", parseFunction=parse_taxonomy_greengenes)

    # checks that your sample and taxa names are consistent across the different slots of the phyloseq object
    # removes taxa that sum to zero across all samples (shouldnt be any)
data <- phyloseq_validate(data, remove_undetected = TRUE)

    # genus and fmaily rank are sometimes NA
tax_table(data)[30:54, 3:7]

    # fix the NAs (take name of higher rank)
data_fixed <- tax_fix(data)

    # check again
tax_table(data_fixed)[40:54, 4:7]

### metadata ----

    # add identifiers to metadata
data_fixed <- data_fixed %>%
  ps_mutate(
    Id = substr(Id, 7, nchar(Id)-23),
    patient = str_before_nth(Id,"_",1),
    region = substr(Id,1,6),
    subregion = str_before_nth(Id,"_",2),
    sample = str_before_nth(Id,"_",3),
    type = substr(str_after_nth(Id,"_",2),1,1),
    group = ifelse( (substr(Id,6,6) %in% c("F","G","H")), "adenoma", 
                   ifelse( (substr(Id,6,6) %in% c("E")), "normal",
                           ifelse( (substr(Id,6,6) %in% c("A","B")), "carcinoma", 
                                   ifelse( (substr(Id,6,6) %in% c("C","D") & patient == "C516"), "adenoma",
                                           ifelse( (substr(Id,6,6) %in% c("C","D") & patient != "C516"), "carcinoma", 
                                                   ifelse( (substr(Id,6,6) %in% c("Z")), "blood", 
                                                           ifelse( (substr(Id,6,6) %in% c("W")), "normalInCancer", NA))))))),
    
  )

### check gland duplicates ----
    
    # calculate largest library size
library_sizes <- data.frame(otu_table(data_fixed))
library_sizes <- data.frame(Id = sample_data(data_fixed)$Id,
                            librarySize = colSums(library_sizes))

    # add to sample data
sample_data(data_fixed)$librarySize <- colSums(data.frame(otu_table(data_fixed)))

keep <- list()
for (i in 1:length(unique(sample_data(data_fixed)$sample))) {
  
  samples <- sample_data(data_fixed)$Id[which(sample_data(data_fixed)$sample == unique(sample_data(data_fixed)$sample)[i])]
  wd <- library_sizes[which(library_sizes$Id == samples),]
  keep[[i]] <- wd$Id[which.max(wd$librarySize)]

}

data_fixed <- data_fixed %>%
  ps_mutate(
    duplicate = ifelse(Id %in% keep, FALSE, TRUE)
  )

data_noDups <- subset_samples(data_fixed, duplicate==FALSE)

### bacteria vs viruses (bacteria) ----

input_raw <- subset_taxa(data_noDups, Kingdom=="Bacteria")
#viruses <- subset_taxa(data_noDups, Kingdom=="Viruses")

### sample counts ----
x <- sample_data(input_raw)
length(unique(x$patient))
nrow(x)
mean(table(x$patient))
min(table(x$patient))
max(table(x$patient))

### remove blood ----

input_raw <- subset_samples(input_raw, group != "blood")

### library size filter (<1000) ----

# add library size to sample_data
jpeg('tempfig.jpeg', width = (40), height = (20), units = "cm", res = 300)
ggplot(sample_data(input_raw), aes(x=reorder(patient, log(librarySize), mean), y=log(librarySize), fill=patient)) +
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  geom_hline(yintercept = log(1000)) +
  geom_jitter(shape=21, size=4, alpha = 0.8, width=0.2, height=0) +
  scale_y_continuous(name="Library size", 
                     breaks = c(5,7.5,10,12.5,15),
                     labels = c(round(exp(5),0),round(exp(7.5),0),round(exp(10),0),
                                  round(exp(12.5),0),round(exp(15),0))) +
  xlab(NULL) +
  theme_custom() +
  theme(
    axis.text.x = element_text(size=24, angle=90, vjust=0.5),
    panel.grid.major.y = element_line(size=0.25),
    panel.grid.minor.y = element_line(size=0.25),
    legend.position = "none")
dev.off()

# filter samples with <1000
input_raw <- subset_samples(input_raw, librarySize > 1000 )


### taxonomic/prevalence filtering ----

# remove taxa with zero counts
input_raw <- phyloseq_validate(input_raw, remove_undetected = TRUE)

# read count per phyum
phylacount <- data.frame(table(tax_table(input_raw)[, "Phylum"], exclude = NULL))

# note phyla with one species
filterPhyla <- phylacount$Var1[phylacount$Freq<2]

# species prevalence across samples
prevdf <- apply(X = otu_table(input_raw),
                MARGIN = ifelse(taxa_are_rows(input_raw), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(input_raw),
                     tax_table(input_raw))

# compute average (and total) prevalences of the species in each phylum.
phyla.prevdf <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
colnames(phyla.prevdf) <- c("phyla","meanPrev","sumPrev")

# check the prevelence of phyla to be filtered out
phyla.prevdf[phyla.prevdf$phyla %in% filterPhyla,]

# relationship of prevalence and total read count for each feature to remove outliers
# Subset to the remaining phyla
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(input_raw, "Phylum"))

# each point is different species
# find the minimum prevalence criteria
plot <- ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(input_raw),color=Phylum)) +
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 5, alpha = 0.7) + # Include a guess for parameter
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum, nrow=3) + 
  theme(legend.position="none") +
  theme_custom() +
  theme(legend.position = "none", 
        panel.grid.major = element_line(size=0.25),
        panel.grid.minor = element_line(size=0.25),
        strip.text = element_text(size=24))

jpeg('tempfig.jpeg', width = (80), height = (40), units = "cm", res = 300)
plot
dev.off()

#jpeg('tempfig.jpeg', width = (3*37.795*15), height = (3*37.795*10))
#plot
#dev.off()

### subset out noise ----

# subset single-species phyla
input_raw <- subset_taxa(input_raw, !Phylum %in% filterPhyla)

# set prevalence threshold (will likely need adjusting)
prevalenceThreshold <- 0.01 * nsamples(input_raw)
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
input_raw <- prune_taxa(keepTaxa, input_raw)

### ======NORMALISATION==================================== ----
### library size stats ----

# average library size
mean(sample_data(input_raw)$librarySize)
min(sample_data(input_raw)$librarySize)
max(sample_data(input_raw)$librarySize)


### quick tsm on options ----

input_tsm <- transform_sample_counts(input_raw, function(ASV) (ASV/sum(ASV) *10000) )

### abundance value transformation (tsm) ----

# New methods include: UQ, CSS (metagennomeSeq), DESeq-VS, edgeR-TMM
# but these standardize the within-sample variance across samples (forcing each sample to have the same distribution of reads)
# which is important for DA, but problematic when comparing entire communities, because variance and evenness are tightly linked
# even=low variance
# so these methods suppress evennesss
# only TSM and rare guarentee number of reads equal across samples, important for bray curtis

# log trans +1 make sense for genes: reduce effect of dominant to see rare. important because of housekeeping
# reducing the importance of dominant OTUs and amplifying the importance of rare OTUs may be misleadingx of the differences among communities.

# in general use TSM (pref) or rare for community-level analysis. If some species as dominant in ALL samples, use CSS to focus on rarer ones
# but must interpret in this context: "uncommon members of the community differ after reducing the importance of the common members"

# also micrbiome analyses would normally violate DESeq and edgeR-TMM normalization assumptions of a constant abundance of a majority of species


# total sum normalization
# Proportions are criticized because they do not account for heteroskedasticity and are bad for DA testing
# alpha requires integer, so round
input_tsm <- transform_sample_counts(input_raw, function(ASV) (ASV/sum(ASV) *10000) )
input_tsm.round <- transform_sample_counts(input_raw, function(ASV) round(ASV/sum(ASV) *10000) )

# rarefying by randomly subsampling each sample to the lowest read depth of any sample
# discards potentially useful data
# bad for DA testing
set.seed(1782) 
input_rarefy <- rarefy_even_depth(input_raw, rngseed = TRUE, replace = FALSE) 

# deseq2: (best dispersion plot when duplicates removed)
# needs fittype=local
# dds1 is now the 0.05 prev filter: 14 ddint converge, fewer on line
# dds2 0.05, keep dups: 57 row, refitting for 2122 genes. fewer points on line
# dds3 0.01 keep dups: 52 row, refitting for 3193 genes. fewer points on line
# dds0 0.05 dup removed, <1000 libsize removed
# no dispersion, only needed for DA. otherwise use DESeq(dds0, fitType="local") + getVarianceStabilizedData(dds)
# alpha requires integer, so round
dds <- phyloseq_to_deseq2(input_raw , ~ subregion)
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))} # calculate geometric means prior to estimate size factors: necessary because lots of zeros
geoMeans <- apply(counts(dds), 1, gm_mean)
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
vst <- varianceStabilizingTransformation(dds, blind=F, fitType = "local") #variance normalizes the data (this includes a log2(x+1) transformation)
input_DESeqVS <- input_raw # create a new phyloseq object with the vst counts
input_DESeqVS.round <- input_raw
otu_table(input_DESeqVS) <- otu_table(as.matrix(assay(vst)), taxa_are_rows = TRUE)
otu_table(input_DESeqVS.round) <- otu_table(round(as.matrix(assay(vst))), taxa_are_rows = TRUE)
sizeFactors(dds)

# edgeR-TMM
# alpha requires integer, so round
wd <- as.data.frame(otu_table(input_raw))
tmm <- edgeR::calcNormFactors(wd, method="TMM") 

counts <- NULL
for(i in 1:length(tmm)){
  col.i <- wd[,i]
  col.i <- (col.i/(tmm[i]*sum(col.i)))*10000
  counts <- cbind(counts, col.i)
}
colnames(counts) <- colnames(wd)
rownames(counts) <- rownames(wd)
input_edgeRTMM <- input_raw # create a new phyloseq object with the tmm counts
input_edgeRTMM.round <- input_raw
otu_table(input_edgeRTMM) <- otu_table(counts, taxa_are_rows = TRUE)
otu_table(input_edgeRTMM.round) <- otu_table(round(counts), taxa_are_rows = TRUE)

# CSS (metagenomics) doesnt work
# alpha requires integer, so round
perc <- cumNormStatFast(obj = as.matrix(otu_table(input_raw))) # percentile for which to scale
css <- cumNormMat(as.matrix(otu_table(input_raw)), perc, 1000)
input_CSS <- input_raw
input_CSS.round <- input_raw
otu_table(input_CSS) <- otu_table(css, taxa_are_rows = TRUE)
otu_table(input_CSS.round) <- otu_table(round(css), taxa_are_rows = TRUE)

# UQ - doesnt work
input_raw <- as.data.frame(otu_table(input_raw))
uq <- edgeR::calcNormFactors(input_raw, method ="upperquartile")
counts <- NULL
for(i in 1:length(uq)){
  col.i <- input_raw[,i]
  col.i <- col.i/(uq[i]*sum(col.i))
  counts <- cbind(counts,col.i)}
colnames(counts) <- colnames(input_raw)
rownames(counts) <- rownames(input_raw)
input_UQ <- input_raw # create a new phyloseq object with the uq counts
otu_table(input_UQ) <- otu_table(counts, taxa_are_rows = TRUE)

### effect of normalization on library size ----
library_sizes <- data.frame(Id = sample_data(input)$Id,
                            Raw = colSums(data.frame(otu_table(input))),
                            TSM = colSums(data.frame(otu_table(input_tsm))),
                            Rarefy = colSums(data.frame(otu_table(input_rarefy))),
                            DESeq2VS = colSums(data.frame(otu_table(input_DESeqVS))),
                            #CSS = colSums(data.frame(otu_table(input_CSS))),
                            EdgeRTMM = colSums(data.frame(otu_table(input_edgeRTMM)))
)

library_sizes.long <- gather(library_sizes, key="key", value="value", -Id)
library_sizes.long$Id <- factor(library_sizes.long$Id,levels = library_sizes$Id[order(library_sizes$Raw)])
library_sizes.long$key <- factor(library_sizes.long$key, levels=c("Raw","TSM","Rarefy","DESeq2VS","EdgeRTMM","CSS"))

ggplot(library_sizes.long, aes(x=Id , y=log(abs(value)), group=key, color=key)) +
  geom_line() +
  theme_custom() +
  ylab("log library size") +
  xlab(NULL) +
  theme(axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8))

### effect of normalization on alpha diversity ----

# plot shannon diversity for raw (ground truth). must use rouded version where necessary
shannon <- data.frame(Id = sample_data(input)$Id,
                      #subregion = sample_data(bacteria_prev2Filter)$subregion,
                      Raw = unlist(estimate_richness(input, measures=c("Shannon"))),
                      TSM = unlist(estimate_richness(input_tsm.round, measures=c("Shannon"))),
                      Rarefy = unlist(estimate_richness(input_rarefy, measures=c("Shannon"))),
                      DESeq2VS = unlist(estimate_richness(input_DESeqVS.round, measures=c("Shannon"))),
                      EdgeRTMM = unlist(estimate_richness(input_edgeRTMM.round, measures=c("Shannon"))),
                      #CSS = unlist(estimate_richness(input_CSS.round, measures=c("Shannon")))
)

shannon.long <- gather(shannon, key="key", value="value", -Id)
shannon.long$Id <- factor(shannon.long$Id,levels = shannon$Id[order(shannon$Raw)])
shannon.long$key <- factor(shannon.long$key, levels=c("Raw","TSM","Rarefy","DESeq2VS","EdgeRTMM","CSS"))

ggplot(shannon.long, aes(x=Id , y=(value), group=key, color=key)) +
  geom_line(alpha=0.5, size=0.8) +
  #facet_grid(rows=vars(key)) +
  ylab("shannon") +
  xlab(NULL) +
  theme_custom() +
  theme(axis.text = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8))

# plot simpson diversity for raw (ground truth). must use rounded version where necessary
simpson <- data.frame(Id = sample_data(input)$Id,
                      #subregion = sample_data(bacteria_prev2Filter)$subregion,
                      Raw = unlist(estimate_richness(input, measures=c("Simpson"))),
                      TSM = unlist(estimate_richness(input_tsm.round, measures=c("Simpson"))),
                      Rarefy = unlist(estimate_richness(input_rarefy, measures=c("Simpson"))),
                      DESeq2VS = unlist(estimate_richness(bacteria_DESeqVS.round, measures=c("Simpson"))),
                      EdgeRTMM = unlist(estimate_richness(input_edgeRTMM.round, measures=c("Simpson"))),
                      #CSS = unlist(estimate_richness(input_CSS.round, measures=c("Simpson")))
)

simpson.long <- gather(simpson, key="key", value="value", -Id)
simpson.long$Id <- factor(simpson.long$Id,levels = simpson$Id[order(simpson$Raw)])
simpson.long$key <- factor(simpson.long$key, levels=c("Raw","TSM","Rarefy","DESeq2VS","EdgeRTMM","CSS"))

ggplot(simpson.long, aes(x=Id , y=(value), group=key, color=key)) +
  geom_line(alpha=0.5, size=1) +
  #facet_grid(rows=vars(key)) +
  theme_custom() +
  theme(axis.text = element_text(size=5, angle=90, vjust=0.5))


### effect of normalization on phyla abundances ----

# define function to calculate prevalence (fraction of samples) 
getPrevalent <- function(phyloseqOb, level, n) {
  phyloseqOb <- phyloseqOb %>%
    aggregate_top_taxa(level = level, top = n) %>%
    microbiome::transform(transform = "compositional")
  x <- apply(X = otu_table(phyloseqOb),
             MARGIN = ifelse(taxa_are_rows(phyloseqOb), yes = 1, no = 2),
             FUN = function(x){sum(x > 0)})
  x <- data.frame(Prevalence = x,
                  TotalAbundance = taxa_sums(phyloseqOb),
                  tax_table(phyloseqOb))
  
  x <- plyr::ddply(x, level, function(df1){cbind(sum(df1$Prevalence),sum(df1$TotalAbundance))})   # compute average (and total) prevalences of the species in each phylum.
  colnames(x) <- c(level,"Prevalence",'TotalAbundance')
  return(x)
}


# compare top 20 genera (prevalence and abundance) across norm methods
# DESeq2 should be ignored

topGen_prev2Filter <- getPrevalent(bacteria_prev2Filter, "Genus", 10)
topGen_tsm <- getPrevalent(bacteria_tsm, "Genus", 10)
topGen_rarefy <- getPrevalent(bacteria_rarefy, "Genus", 10)
topGen_DESeqVS <- getPrevalent(bacteria_DESeqVS, "Genus", 10)
topGen_edgeRTMM <- getPrevalent(bacteria_edgeRTMM, "Genus", 10)
topGen_CSS <- getPrevalent(bacteria_CSS, "Genus", 10)

abundance_comparison <- data.frame(none=topGen_prev2Filter$Genus[order(topGen_prev2Filter$TotalAbundance, decreasing = TRUE)],
                                   tsm=topGen_tsm$Genus[order(topGen_tsm$TotalAbundance, decreasing = TRUE)],
                                   rarefy=topGen_rarefy$Genus[order(topGen_rarefy$TotalAbundance, decreasing = TRUE)],
                                   DESeq2VS=topGen_DESeqVS$Genus[order(topGen_DESeqVS$TotalAbundance, decreasing = TRUE)],
                                   edgeRRMM=topGen_edgeRTMM$Genus[order(topGen_edgeRTMM$TotalAbundance, decreasing = TRUE)],
                                   CSS=topGen_CSS$Genus[order(topGen_CSS$TotalAbundance, decreasing = TRUE)])


# consider relationship between prevalence and relative abundance with each method
# rarefy thins out the data
# DESeq2 completely flattens
# tsm, CSS, and edgeR produce same relative abundance/prevalence at genus level - persue these options

x <- bacteria_tsm %>%
  aggregate_taxa(level = "Species") %>%
  microbiome::transform(transform = "compositional")
phylacount <- data.frame(table(tax_table(x)[, "Phylum"], exclude = NULL))

prevdf <- apply(X = otu_table(x),
                MARGIN = ifelse(taxa_are_rows(x), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prevdf,
                     RelativeAbundance = taxa_sums(x),
                     tax_table(x))

ggplot(prevdf, aes(RelativeAbundance, Prevalence / nsamples(x),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 1, alpha = 0.7) + # Include a guess for parameter
  scale_x_log10() +  xlab("Relative Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# tsm, CSS, and edgeR have same rel abundances
# would have same shannon but needed to round  
# but CSS and edgeR dont guarentee equal library size which will mess up BC

# tsm is best for stable library size, but has delftia in top 10 (as does rarefy)
# deseq2/edgeR as a no
# css top genera makes more sense, lib size is more stable than deseq2, but still not guarenteed

### ======CONTAMINATION ANALYSIS=========================== ----
### evidence contamination at phylum level----

x <- input_tsm %>%
  aggregate_top_taxa("Phylum", top=10) %>%
  microbiome::transform(transform = "compositional")

x <- data.frame(otu_table(x))
x$phylum <- rownames(x)
x <- gather(x[which(x$phylum!="Other"),], key="key", value="value", -phylum)

nrow(x[which(x$phylum=="Proteobacteria" & x$value>0.7),])/length(unique(x$key))

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(x=reorder(phylum,-value,mean,na.rm=TRUE) , y=(value))) +
  geom_boxplot(aes(fill=phylum), alpha=0.5) +
  geom_jitter(aes(fill=phylum), alpha=0.7, shape=21, width = 0.23) +
  xlab(NULL) +
  ylab("Relative abundance \n(phylum level)") +
  theme_custom() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=22),
        axis.text.y = element_text(size=22),
        axis.title = element_text(size=22),
        legend.position = "none")
dev.off()

### evidence contamination at genus level----

x <- input_tsm %>%
  aggregate_top_taxa("Genus", top=10) %>%
  microbiome::transform(transform = "compositional")

x <- data.frame(otu_table(x))
x$genus <- rownames(x)
x <- gather(x[which(x$genus!="Other"),], key="key", value="value", -genus)

nrow(x[which(x$genus=="Pseudomonas" & x$value>0.7),])/length(unique(x$key))

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(x=reorder(genus,-value,mean,na.rm=TRUE) , y=(value))) +
  geom_boxplot(aes(fill=genus), alpha=0.5) +
  geom_jitter(aes(fill=genus), alpha=0.7, shape=21, width = 0.23) +
  xlab(NULL) +
  ylab("Relative abundance \n(genus level)") +
  theme_custom() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=22),
        axis.text.y = element_text(size=22),
        axis.title = element_text(size=22),
        legend.position = "none")
dev.off()

jpeg('tempfig.jpeg', width = (40), height = (60), units = "cm", res = 300)
input_tsm %>%
  ps_arrange(desc(patient), desc(region)) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 15,
    label = "subregion",
    bar_outline_colour = NA,
    sample_order = "default",
    bar_width = 0.7,
    merge_other = TRUE,
    taxon_renamer = toupper
  ) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  facet_grid(rows=vars(patient), scales = "free",space="free", switch = "y") +
  theme(axis.text.x = element_text(size=24),
        axis.title.x = element_text(size=26),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(size=22, angle=0),
        strip.background = element_rect(colour="white", fill="white"),
        #panel.spacing = unit(.05, "lines"),
        legend.position = "top",
        legend.text = element_text(size=18),
        legend.title = element_text(size=22),
        panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
dev.off()

#jpeg('tempfig.jpeg', width = (3*37.795*30), height = (3*37.795*15))
#plot
#dev.off()


### define contamination genera ----

# removal of potential contam done here: https://www.sciencedirect.com/science/article/pii/S1567134820300411
# genera from here: https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-014-0087-z

# contam genera
contam <- c('Afipia', 'Aquabacterium', 'Asticcacaulis', 'Aurantimonas', 'Beijerinckia', 'Bosea', 'Bradyrhizobium', 'Brevundimonas', 
            'Caulobacter', 'Craurococcus', 'Devosia', 'Hoeflea', 'Mesorhizobium', 'Methylobacterium', 'Novosphingobium', 'Ochrobactrum', 
            'Paracoccus', 'Pedomicrobium', 'Phyllobacterium', 'Rhizobium', 'Roseomonas', 'Sphingobium', 'Sphingomonas', 'Sphingopyxis',
            'Acidovorax', 'Azoarcus', 'Azospira', 'Burkholderia', 'Comamonas', 'Cupriavidus', 'Curvibacter', 'Delftia', 'Duganella', 
            'Herbaspirillum', 'Janthinobacterium', 'Kingella', 'Leptothrix', 'Limnobacter', 'Massilia', 'Methylophilus', 
            'Methyloversatilis', 'Oxalobacter', 'Pelomonas', 'Polaromonas', 'Ralstonia', 'Schlegelella', 'Sulfuritalea', 'Undibacterium', 
            'Variovorax',
            'Acinetobacter', 'Enhydrobacter', 'Enterobacter', 'Escherichia', 'Nevskia', 'Pseudomonas', 'Pseudoxanthomonas', 'Psychrobacter', 
            'Stenotrophomonas', 'Xanthomonas',
            'Aeromicrobium', 'Arthrobacter', 'Beutenbergia', 'Brevibacterium', 'Corynebacterium', 'Curtobacterium', 'Dietzia', 'Geodermatophilus', 
            'Janibacter', 'Kocuria', 'Microbacterium', 'Micrococcus', 'Microlunatus', 'Patulibacter', 'Propionibacterium', 'Rhodococcus', 'Tsukamurella',
            'Abiotrophia', 'Bacillus', 'Brevibacillus', 'Brochothrix', 'Facklamia', 'Paenibacillus', 'Streptococcus',
            'Chryseobacterium', 'Dyadobacter', 'Flavobacterium', 'Hydrotalea', 'Niastella', 'Olivibacter', 'Pedobacter', 'Wautersiella',
            'Deinococcus')
contam <- c(contam, "Cutibacterium")

### UHGP and IGC bacteria ----

# pull list of UHGP bacteria
UHGP <- read_excel("~/Documents/Microbial/data/UHGP_2021_gutSpecies.xlsx", skip = 1)
UHGP$NCBIgenus <- str_before_first( str_after_first(UHGP$`Taxonomy lineage (NCBI)`, pattern = "g__"), pattern = ";s__")
UHGP$NCBIspecies <- str_after_first(UHGP$`Taxonomy lineage (NCBI)`, pattern = "s__")
UHGP$NCBIspeciesOnly <- str_after_first(UHGP$NCBIspecies, pattern = " ")
UHGP$GTDBgenus <- str_before_first( str_after_first(UHGP$`Taxonomy lineage (GTDB)`, pattern = "g__"), pattern = ";s__")
UHGP$GTDBspecies <- str_after_first(UHGP$`Taxonomy lineage (GTDB)`, pattern = "s__")
UHGP$GTDBspeciesOnly <- str_after_first(UHGP$GTDBspecies, pattern = " ")
UHGP <- UHGP[c(25:32)]

# IGC
IGC <- read_excel("~/Documents/Microbial/data/wang2014_IGC.xlsx", skip = 2)
output <- list()
for (i in 1:nrow(IGC)) { # get genus
  wd <- IGC$`Detailed taxon level`[i]
  indexG <- unlist(gregexpr("\\|g", wd)) + 1
  indexPipe <- unlist(gregexpr("\\|", wd))
  index <- sum(indexPipe<indexG)
  output[[i]] <- str_before_nth( str_after_nth(IGC$`Detailed taxon name`[i], '\\|', index), '\\|', 1)
}
IGC$Genus <- unlist(output)

output <- list()
for (i in 1:nrow(IGC)) { # get species
  wd <- IGC$`Detailed taxon level`[i]
  indexS <- unlist(gregexpr("\\|species", wd)) + 1
  indexPipe <- unlist(gregexpr("\\|", wd))
  index <- sum(indexPipe<indexS)
  if( sum(indexPipe>indexS) > 0 ) { #ie there are more pipes
    output[[i]] <- str_after_nth( str_before_nth( str_after_nth(IGC$`Detailed taxon name`[i], '\\|', index), '\\|', 1), ' ',1)
  } else {
    output[[i]] <- str_after_nth(IGC$`Detailed taxon name`[i], ' ',1)
  }
}

IGC$SpeciesOnly <- unlist(output)

# find matches
tax <- data.frame(tax_table(input_tsm))
tax$GS <- paste(tax$Genus, tax$Species, sep=" ")
tax$inUHGP <- ifelse(tax$GS %in% c(UHGP$GTDBspecies, UHGP$NCBIspecies), TRUE, FALSE)
tax$inIGC <- ifelse(tax$GS %in% c(IGC$Species), TRUE, FALSE)
sum(tax$inUHGP)
sum(tax$inIGC)

# make lists of gut genera
gutG <- unique(c(IGC$Genus, UHGP$GTDBgenus, UHGP$NCBIgenus))
gutnocontamG <- gutG[gutG %!in% contam]
nocontamG <- tax$Genus[tax$Genus %!in% contam]
gutS <- unique(c(IGC$Species, UHGP$GTDBspecies, UHGP$NCBIspecies))
gutnocontamS <- unique(c(IGC$Species[which(IGC$Genus %!in% contam)], 
                         UHGP$GTDBspecies[which(UHGP$GTDBgenus %!in% contam)], 
                         UHGP$NCBIspecies[which(UHGP$NCBIgenus %!in% contam)]))


### variance of contam vs noncontam genera ----

# get otu table at genus level and make compositional
x <- input_tsm %>%
  aggregate_taxa("Genus") %>%
  transform_sample_counts(function(ASV) (ASV/sum(ASV) *10000) )
otu <- data.frame(otu_table(x))

# for each patient, get variance per genus
output <- list()
for ( i in 1:length(unique(sample_data(x)$patient)) ) {
  wd <- data.frame(otu[, which(sample_data(x)$patient == unique(sample_data(x)$patient)[i])])
  output[[i]] <- data.frame(patient = unique(sample_data(x)$patient)[i], 
                            Variance = apply(wd, 1, var), 
                            St.dev = apply(wd, 1, sd), 
                            mean = apply(wd, 1, mean))
  output[[i]]$genus <- rownames(output[[i]])
  output[[i]]$COV <- output[[i]]$St.dev / output[[i]]$mean
  output[[i]]$group <- unlist(ifelse(output[[i]]$genus %in% nocontamG, "Unlikely contaminating \ngenera", "Potential contaminating \ngenera"))
}

output <- do.call(rbind, output)


#hm
jpeg('tempfig.jpeg', width = (40), height = (65), units = "cm", res = 300)
ggplot(output, aes(x=patient, y=genus, fill = COV)) +
  geom_tile() +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  facet_grid(group ~ patient, scales='free', space = 'free', switch = "y") +
  theme_custom() +
  theme(
    legend.position = "top",
    legend.text = element_text(size=24),
    legend.title = element_text(size=24),
    legend.key.width = unit(4, "cm"),
    panel.spacing.y = unit(1.2, "lines"),
    panel.spacing.x = unit(0.5, "lines"),
    panel.background = element_blank(),
    strip.text.y.left =  element_text(size=24),
    strip.text.x = element_blank(),
    axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=24),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size = 1)) +
  guides(fill=guide_colourbar(title="Coefficient of variance"))
dev.off()


# plot with each dot = one genus variance of a patient
jpeg('tempfig.jpeg', width = (40), height = (30), units = "cm", res = 300)
ggplot(output, aes(x=group , y=(St.dev/mean), fill=patient)) +
  geom_boxplot(alpha=0.6, width=0.7, position = position_dodge(0.9), outlier.shape = NA) +
  geom_jitter(shape=21, position = position_dodge(0.9)) +
  ylim(0,5) +
  ylab("Coefficient of variance") +
  stat_compare_means(method = "t.test", size=9, comparisons = list(c("Unlikely contaminating \ngenera","Potential contaminating \ngenera"))) +
  theme_custom() +
  theme(axis.title.x = element_blank(),
        legend.text = element_text(size=28), legend.title = element_text(size=28), 
        legend.position = "top")
dev.off()

output$stdevMean <- output$St.dev/output$mean
x <- aggregate(output$stdevMean, list(output$patient, output$group), FUN=mean, na.rm=TRUE)
colnames(x) <- c("patient","group","AverageStdevMean")

mean(x$AverageStdevMean[which(x$group=="Potential contaminating \ngenera")])
mean(x$AverageStdevMean[which(x$group=="Unlikely contaminating \ngenera")])


# plot with each dot = one patient per group
output$stdevMean <- output$St.dev/output$mean
x <- aggregate(output$stdevMean, list(output$patient, output$group), FUN=mean, na.rm=TRUE)
colnames(x) <- c("patient","group","AverageStdevMean")
mean(x$AverageStdevMean[which(x$group=="Potential contaminating \ngenera")])
mean(x$AverageStdevMean[which(x$group=="Unlikely contaminating \ngenera")])

jpeg('tempfig.jpeg', width = (30), height = (20), units = "cm", res = 300)
ggplot(x, aes(x=group , y=AverageStdevMean)) +
  geom_violin(aes(fill=group), alpha = 0.7) +
  geom_boxplot(alpha=0.8, width=0.5) +
  geom_line(aes(group=patient)) +
  geom_jitter(size = 3, shape=21, alpha=1, width = 0.1, height = 0, aes(fill=group)) +
  ylim(0,3) +
  #geom_text_repel(aes(label=patient)) +
  #ylab("Standard deviation / mean \n(averaged per patient)") +
  ylab("Average coefficient of \nvariation per patient") +
  scale_fill_manual(values = c("#9966CC","#FF9933")) +
  stat_compare_means(method = "t.test", paired = TRUE, comparisons = list(c("Potential contaminating \ngenera","Unlikely contaminating \ngenera")), size=9) +
  theme_custom() +
  theme(legend.position = "none", axis.title.x = element_blank())
dev.off()









# plot each patient separately             
ggplot(output, aes(x=group , y=(St.dev/mean), fill=group)) +
  geom_boxplot(alpha=0.8) +
  facet_grid(cols=vars(patient), switch = "both") +
  xlab(NULL) +
  #geom_text_repel(aes(label=genus)) +
  ylab("Standard deviation / mean") +
  scale_fill_manual(values = c("#9966CC","#FF9933")) +
  #geom_jitter() +
  stat_compare_means(method = "t.test", comparisons = list(c("Potential \nContaminents","Unlikely \nContaminents")),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) +
  theme_custom() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.margin=unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

ggplot(output, aes(x=group , y=(St.dev/mean), color=patient),) +
  geom_boxplot(aes(x=group , y=(St.dev/mean), fill=group, color=NULL), width=0.9, alpha=0.6) +
  geom_boxplot(alpha=0.6, width=0.7, position = position_dodge(0.85)) +
  #geom_text_repel(aes(label=genus)) +
  ylab("Standard deviation / mean") +
  scale_fill_manual(values = c("#9966CC","#FF9933")) +
  stat_compare_means(method = "t.test", comparisons = list(c("Potential \nContaminents","Unlikely \nContaminents"))) +
  theme_custom() 



# plot with each dot = one genus, average variance across patients
output$stdevMean <- output$St.dev/output$mean
x <- aggregate(output$stdevMean, list(output$genus), FUN=mean, na.rm=TRUE)
x$group <- output[match(x$Group.1, output$genus), 6]
colnames(x) <- c("genus","AverageStdevMean", "group")

ggplot(x, aes(x=group , y=AverageStdevMean)) +
  geom_violin(aes(fill=group)) +
  #geom_line(aes(group=patient)) +
  geom_boxplot(alpha=0.8) +
  geom_jitter(shape=21, aes(fill=group)) +
  geom_text_repel(aes(label=genus)) +
  ylab("Standard deviation / mean \n(averaged per patient)") +
  scale_fill_manual(values = c("#9966CC","#FF9933")) +
  stat_compare_means(method = "t.test",comparisons = list(c("Potential \nContaminent","Unlikely \nContaminent"))) +
  theme_custom() 


### variance of random vs random genera across patients ----

# get otu table at genus level and make compositional
x <- input_tsm %>%
  aggregate_taxa("Genus") %>%
  transform_sample_counts(function(ASV) (ASV/sum(ASV) *10000) )
otu <- data.frame(otu_table(x))

# for each patient, get variance per genus
output <- list()
for ( i in 1:length(unique(sample_data(x)$patient)) ) {
  wd <- data.frame(otu[, which(sample_data(x)$patient == unique(sample_data(x)$patient)[i])])
  output[[i]] <- data.frame(patient = unique(sample_data(x)$patient)[i], 
                            Variance = apply(wd, 1, var), 
                            St.dev = apply(wd, 1, sd), 
                            mean = apply(wd, 1, mean))
  output[[i]]$genus <- rownames(output[[i]])
  output[[i]]$COV <- output[[i]]$St.dev / output[[i]]$mean
  output[[i]]$group <- unlist(ifelse(output[[i]]$genus %in% nocontamG, "Unlikely contaminating \ngenera", "Potential contaminating \ngenera"))
}

output <- do.call(rbind, output)





# select random genus
t <- list()
for (i in 1:10000) {
  set.seed(i)
  random <- sample(output$genus, 93, replace = FALSE)
  x <- output
  x$group <- unlist(ifelse(x$genus %in% random, "Potential contaminating \ngenera", "Unlikely contaminating \ngenera"))
  x <- aggregate(x$COV, list(x$patient, x$group), FUN=mean, na.rm=TRUE)
  colnames(x) <- c("patient","group","AveCOV")
  
  ttest <- t.test(x$AveCOV[which(x$group == "Potential contaminating \ngenera")],
                  x$AveCOV[which(x$group == "Unlikely contaminating \ngenera")],
                  paired = TRUE)
  
  #if( (ttest$estimate[1] < ttest$estimate[2]) & (ttest$p.value<0.05) ) {
  if(  (ttest$p.value<0.05) ) {
    t[[i]] <- TRUE
  } else {
    t[[i]] <- FALSE
  }
}

sum(unlist(t))/10000

### variance of IGC gut vs not gut species (x) ----

#across patient
x <- input_tsm %>%
  aggregate_taxa("Species") 
otu <- data.frame(otu_table(x))
tax <- data.frame(tax_table(x))
rownames(otu) <- paste(tax$Genus, tax$Species, sep=" ")

output <- list()
for ( i in 1:length(unique(sample_data(x)$patient)) ) {
  wd <- data.frame(otu[, which(sample_data(x)$patient == unique(sample_data(x)$patient)[i])])
  output[[i]] <- data.frame(patient = unique(sample_data(x)$patient)[i], 
                            Variance = apply(wd, 1, var), 
                            St.dev = apply(wd, 1, sd), 
                            mean = apply(wd, 1, mean))
  output[[i]]$GS <- rownames(output[[i]])
  output[[i]]$group <- unlist(ifelse(output[[i]]$GS %in% IGC$Species, "In ICG", "Not in ICG"))
}

output <- do.call(rbind, output)

ggplot(output, aes(x=group , y=(St.dev/mean), fill=group)) +
  #geom_violin() +
  geom_boxplot(alpha=0.8) +
  facet_grid(cols=vars(patient), switch = "both") +
  xlab(NULL) +
  #geom_text_repel(aes(label=genus)) +
  ylab("Standard deviation / mean") +
  scale_fill_manual(values = c("#9966CC","#FF9933")) +
  #geom_jitter() +
  stat_compare_means(method = "t.test", comparisons = list(c("In ICG", "Not in ICG")),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) +
  theme_custom() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.margin=unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

ggplot(output, aes(x=group , y=(St.dev/mean))) +
  geom_violin(aes(fill=group)) +
  geom_boxplot(alpha=0.8) +
  geom_jitter(shape=21, aes(fill=group)) +
  #geom_text_repel(aes(label=genus)) +
  ylab("Standard deviation / mean") +
  scale_fill_manual(values = c("#9966CC","#FF9933")) +
  stat_compare_means(method = "t.test", comparisons = list(c("In ICG", "Not in ICG"))) +
  theme_custom() 




# by patient
x <- input_tsm %>%
  subset_samples(group=="carcinoma" & type=="G") %>%
  aggregate_taxa("Species")
otu <- data.frame(otu_table(x))
tax <- data.frame(tax_table(x))
rownames(otu) <- paste(tax$Genus, tax$Species, sep=" ")

output <- data.frame(Variance = apply(otu, 1, var), 
                     St.dev = apply(otu, 1, sd), 
                     mean = apply(otu, 1, mean))
output$GS <- rownames(output)
output$group <- unlist(ifelse(output$GS %in% IGC$Species, "In ICG", "Not in ICG"))

ggplot(output, aes(x=group , y=(St.dev/mean))) +
  geom_violin(aes(fill=group)) +
  geom_boxplot(alpha=0.8) +
  geom_jitter(shape=21, aes(fill=group)) +
  #geom_text_repel(aes(label=genus)) +
  ylab("Standard deviation / mean") +
  scale_fill_manual(values = c("#9966CC","#FF9933")) +
  stat_compare_means(method = "t.test", comparisons = list(c("In ICG", "Not in ICG"))) +
  theme_custom() 

### variance of UHGP gut vs not gut species (x) ----

# across patients
x <- input_tsm %>%
  aggregate_taxa("Species") %>%
  subset_samples(group=="carcinoma" & type=="G")
otu <- data.frame(otu_table(x))
tax <- data.frame(tax_table(x))
rownames(otu) <- paste(tax$Genus, tax$Species, sep=" ")

output <- list()
for ( i in 1:length(unique(sample_data(x)$patient)) ) {
  wd <- otu[, which(sample_data(x)$patient == unique(sample_data(x)$patient)[i])]
  output[[i]] <- data.frame(patient = unique(sample_data(x)$patient)[i], 
                            Variance = apply(wd, 1, var), 
                            St.dev = apply(wd, 1, sd), 
                            mean = apply(wd, 1, mean))
  output[[i]]$GS <- rownames(output[[i]])
  output[[i]]$group <- unlist(ifelse(output[[i]]$GS %in% c(UHGP$GTDBspecies, UHGP$NCBIspecies), "In UHGP", "Not in UHGP"))
}

output <- do.call(rbind, output)

ggplot(output, aes(x=group , y=(St.dev/mean), fill=group)) +
  #geom_violin() +
  geom_boxplot(alpha=0.8) +
  facet_grid(cols=vars(patient), switch = "both") +
  xlab(NULL) +
  #geom_text_repel(aes(label=genus)) +
  ylab("Standard deviation / mean") +
  scale_fill_manual(values = c("#9966CC","#FF9933")) +
  #geom_jitter() +
  stat_compare_means(method = "t.test", comparisons = list(c("In UHGP", "Not in UHGP")),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) +
  theme_custom() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.margin=unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

# by patient
x <- input_tsm %>%
  subset_samples(group=="carcinoma" & type=="G") %>%
  aggregate_taxa("Species")
otu <- data.frame(otu_table(x))
tax <- data.frame(tax_table(x))
rownames(otu) <- paste(tax$Genus, tax$Species, sep=" ")

output <- data.frame(Variance = apply(otu, 1, var), 
                     St.dev = apply(otu, 1, sd), 
                     mean = apply(otu, 1, mean))
output$GS <- rownames(output)
output$group <- unlist(ifelse(output$GS %in% c(UHGP$GTDBspecies, UHGP$NCBIspecies), "In UHGP", "Not in UHGP"))

ggplot(output, aes(x=group , y=(St.dev/mean))) +
  geom_violin(aes(fill=group)) +
  geom_boxplot(alpha=0.8) +
  geom_jitter(shape=21, aes(fill=group)) +
  #geom_text_repel(aes(label=genus)) +
  ylab("Standard deviation / mean") +
  scale_fill_manual(values = c("#9966CC","#FF9933")) +
  stat_compare_means(method = "t.test", comparisons = list(c("In UHGP", "Not in UHGP"))) +
  theme_custom() 

### set up other options for contamination removal ----

tax <- data.frame(tax_table(input_tsm))

# just exclude contam genera
keep <- rownames(tax[which(tax$Genus %in% nocontamG),])
bacteria_nocontamG <- prune_taxa(keep, input_tsm)

# genera in ICG/UHGP
keep <- rownames(tax[which(tax$Genus %in% gutG),])
bacteria_gutG <- prune_taxa(keep, input_tsm)

# genera in ICG/UHGP contam removed
keep <- rownames(tax[which(tax$Genus %in% gutnocontamG),])
bacteria_gutnocontamG <- prune_taxa(keep, input_tsm)

# species in ICG/UHGP
keep <- rownames(tax[which(paste(tax$Genus, tax$Species, sep=' ') %in% gutS),])
bacteria_gutS <- prune_taxa(keep, input_tsm)

### species in ICG/UHGP contam removed ----
tax <- data.frame(tax_table(input_tsm))
keep <- rownames(tax[which(paste(tax$Genus, tax$Species, sep=' ') %in% gutnocontamS),])
input_tsm_decon <- prune_taxa(keep, input_tsm)
input_raw_decon <- prune_taxa(keep, input_raw)


### ======EFFECT OF COMTAM REMOVAL=========================
### boxplot of phyla ----
x <- input_tsm_decon %>%
  aggregate_top_taxa("Phylum", top=10) %>%
  microbiome::transform(transform = "compositional")

x <- data.frame(otu_table(x))
x$phylum <- rownames(x)
x <- gather(x[which(x$phylum!="Other"),], key="key", value="value", -phylum)

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(x=reorder(phylum,-value,median,na.rm=TRUE) , y=(value))) +
  geom_boxplot(aes(fill=phylum), alpha=0.5) +
  geom_jitter(aes(fill=phylum), alpha=0.7, shape=21, width = 0.23) +
  xlab(NULL) +
  ylab("Relative abundance \n(phylum level)") +
  theme_custom() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=22),
        axis.text.y = element_text(size=22),
        axis.title = element_text(size=22),
        legend.position = "none")
dev.off()


x <- data.frame(aggregate(x$value, list(x$phylum), FUN=mean), 
                aggregate(x$value, list(x$phylum), FUN=min),
                aggregate(x$value, list(x$phylum), FUN=max))


### boxplot of genera ----
x <- input_tsm_decon %>%
  aggregate_top_taxa("Genus", top=10) %>%
  microbiome::transform(transform = "compositional")

x <- data.frame(otu_table(x))
x$genus <- rownames(x)
x <- gather(x[which(x$genus!="Other"),], key="key", value="value", -genus)

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(x=reorder(genus,-value,median,na.rm=TRUE) , y=(value))) +
  geom_boxplot(aes(fill=genus), alpha=0.5) +
  geom_jitter(aes(fill=genus), alpha=0.7, shape=21, width = 0.23) +
  xlab(NULL) +
  ylab("Relative abundance \n(genus level)") +
  theme_custom() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=22),
        axis.text.y = element_text(size=22),
        axis.title = element_text(size=22),
        legend.position = "none")
dev.off()


x <- data.frame(aggregate(x$value, list(x$genus), FUN=mean), 
                aggregate(x$value, list(x$genus), FUN=min),
                aggregate(x$value, list(x$genus), FUN=max))

### boxplot of species ----
x <- input_tsm_decon %>%
  aggregate_taxa("Species") %>%
  microbiome::transform(transform = "compositional")

tax <- data.frame(tax_table(x))
x <- data.frame(otu_table(x))
x$genus <- tax[match(rownames(x), tax$unique), 6]
x$species <- tax[match(rownames(x), tax$unique), 7]
keep <- names(head(sort(rowMeans(x[,-which(colnames(x) %in% c("genus","species"))]), decreasing = TRUE), 20))
x <- gather(x[which(x$species %in% keep),], key="key", value="value", -species, -genus)

ggplot(x, aes(x=reorder(species,-value,median,na.rm=TRUE) , y=(value),  color=species)) +
  #geom_violin() +
  geom_boxplot() +
  xlab(NULL) +
  ylab("Phylogenetic abundance \n(species level)") +
  #geom_jitter() +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        legend.position = "none")

### PCAs to visualize any separation of contamination proportion ----

# calculate proportion of genera are contam or notcontam&gut
y <- input_tsm %>%
  aggregate_taxa("Genus") %>%
  subset_samples(group=="carcinoma" & type=="G") 
  
x <- t(data.frame(otu_table(y)))
propNocontam <- rowSums(x[,which(colnames(x) %!in% contam)]) / rowSums(x)
propGutG <- rowSums(x[,which(colnames(x) %in% gutG)]) / rowSums(x)
propGutnocontamG <- rowSums(x[,which(colnames(x) %in% gutnocontamG)]) / rowSums(x)

# needs log to help spread
x <- y %>%
  transform_sample_counts(function(x) log(1 + x))

pca_input <- t(data.frame(otu_table(x)))
pca <- prcomp(pca_input)
percentVar <- round(100 * summary(pca)$importance[2,])
df <- cbind(sample_data(x), propNocontam, propGutG, propGutnocontamG,pca$x)

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(df, aes(x=PC1, y=PC2, color=1-propGutnocontamG)) + 
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("Proportion of genera \nlikely to be contamination") +
  scale_colour_gradient(
    low = "#FF9933",
    high = "#663399") +
  #coord_fixed() +
  theme_custom() +
  theme(legend.position = "right",
        panel.grid.major = element_line(size=0.25),
        panel.grid.minor = element_line(size=0.25),
        legend.text = element_text(size=28),  
        legend.title = element_text(size=28), legend.key.height = unit(2,"cm")
  ) +
  guides(color=guide_colourbar(title="",
                               title.position = "top", limits = c(0,1), breaks = c(0, 1)))
dev.off()


# what is PC1
x <- data.frame(rotation=pca$rotation[,1])
y <- data.frame(tax_table(input_tsm))
#x$genus <- y[match(rownames(x), y$Species), 6]
#x$species <- y[match(rownames(x), rownames(y)), 7]
x$contam <- ifelse(rownames(x) %in% contam, TRUE, FALSE)
x$genus <- rownames(x)
#colnames(x) <- c("rotation","genus","species","contam")
x <- slice_max(x, order_by = abs(x$rotation), n=200)

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(x=reorder(genus, rotation), y=rotation, fill=contam)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#FF9933","#663399")) +
  ylab("PC1 rotation") +
  xlab("Genera") +
  theme_custom() +
  coord_flip() +
  theme(axis.text.y = element_blank(),
        legend.position = "top", legend.text = element_text(size=24),
        legend.title = element_text(size=24)) +
  guides(fill=guide_legend(title="Contaminating genera",
                               title.position = "top", limits = c(0,1), breaks = c(0, 1)))
dev.off()


### plot detailed abundance of genera per sample before decontam----

# plot abundance (not relative) in top 10 genera in carcinoma glands per sample (duplicates removed)
# does work in species because some species assigned to >1 genus

plot <- input_tsm %>%
  #ps_filter(type =="G") %>%
  #ps_filter(patient %!in% c("C543","C549")) %>%
  #tax_fix(unknowns = c("africanus", "albus", "anatis", "aurantiaca", "avium", "baltica", "canadensis", "capsulatus", "colombiense", "denitrificans", "denticola", "desulfuricans", "donghaensis", "elongata", "faecium", "felis", "fermentans", "ferrooxidans", "gingivalis", "haemolyticus", "halophilus", "hominis", "hongkongensis", "hydrothermalis", "intermedia", "jejuni", "koreensis", "litoralis", "marina", "marinum", "marinus", "maris", "massiliensis", "michiganensis", "mobilis", "muelleri", "oceani", "parvum", "pneumoniae", "putrefaciens", "roseum", "ruber", "salexigens", "salmonicida", "sputigena", "succinogenes", "thermophila", "thermophilum", "thermophilus", "vaginalis", "versatilis", "vulgaris")) %>%
  ps_arrange(desc(patient), desc(region)) %>%
  comp_barplot(
    #tax_level = "Species", n_taxa = 10,
    tax_level = "Genus", n_taxa = 20,
    #tax_level = "Phylum", n_taxa = 5,
    #facet_by = "patient",
    label = "subregion",
    bar_outline_colour = NA,
    sample_order = "default",
    bar_width = 0.7,
    merge_other = TRUE,
    taxon_renamer = toupper
  ) + 
  scale_y_continuous(expand = c(0, 0)) +
  #coord_flip() +
  facet_grid(cols=vars(patient), scales = "free",space="free") +
  #facet_grid(facets = "patient", scales = "free",space="free" ) +
  theme(axis.text.y = element_text(size=26),
        axis.title.y = element_text(size=26),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=10),
        #strip.background = element_blank(),
        strip.text = element_text(size=26),
        #panel.spacing = unit(.05, "lines"),
        legend.position = "top",
        panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

jpeg('tempfig.jpeg', width = (3*37.795*30), height = (3*37.795*15))
plot
dev.off()


### plot detailed abundance of genera per sample after decontam----

# plot abundance (not relative) in top 10 genera in carcinoma glands per sample (duplicates removed)
# does work in species because some species assigned to >1 genus

plot <- input_tsm_decon %>%
  #ps_filter(type =="G") %>%
  #ps_filter(patient %!in% c("C543","C549")) %>%
  #tax_fix(unknowns = c("africanus", "albus", "anatis", "aurantiaca", "avium", "baltica", "canadensis", "capsulatus", "colombiense", "denitrificans", "denticola", "desulfuricans", "donghaensis", "elongata", "faecium", "felis", "fermentans", "ferrooxidans", "gingivalis", "haemolyticus", "halophilus", "hominis", "hongkongensis", "hydrothermalis", "intermedia", "jejuni", "koreensis", "litoralis", "marina", "marinum", "marinus", "maris", "massiliensis", "michiganensis", "mobilis", "muelleri", "oceani", "parvum", "pneumoniae", "putrefaciens", "roseum", "ruber", "salexigens", "salmonicida", "sputigena", "succinogenes", "thermophila", "thermophilum", "thermophilus", "vaginalis", "versatilis", "vulgaris")) %>%
  ps_arrange(desc(patient), desc(region)) %>%
  comp_barplot(
    #tax_level = "Species", n_taxa = 10,
    tax_level = "Genus", n_taxa = 20,
    #tax_level = "Phylum", n_taxa = 5,
    #facet_by = "patient",
    label = "subregion",
    bar_outline_colour = NA,
    sample_order = "default",
    bar_width = 0.7,
    merge_other = TRUE,
    taxon_renamer = toupper
  ) + 
  scale_y_continuous(expand = c(0, 0)) +
  #coord_flip() +
  facet_grid(cols=vars(patient), scales = "free",space="free") +
  #facet_grid(facets = "patient", scales = "free",space="free" ) +
  theme(axis.text.y = element_text(size=26),
        axis.title.y = element_text(size=26),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=10),
        #strip.background = element_blank(),
        strip.text = element_text(size=26),
        #panel.spacing = unit(.05, "lines"),
        legend.position = "top",
        panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

jpeg('tempfig.jpeg', width = (3*37.795*30), height = (3*37.795*15))
plot
dev.off()


### PCA to ID outlier samples ----

# needs log to help spread
x <- input_tsm_decon %>%
  aggregate_taxa("Genus") %>%
  subset_samples(group=="carcinoma" & type=="G") %>%
  transform_sample_counts(function(x) log(1 + x))

pca_input <- t(data.frame(otu_table(x)))
pca <- prcomp(pca_input)
percentVar <- round(100 * summary(pca)$importance[2,])
df <- cbind(sample_data(x),pca$x)

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(df, aes(x=PC1, y=PC2, color=patient)) + 
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_custom() +
  theme(legend.position = "none",
        panel.grid.major = element_line(size=0.25),
        panel.grid.minor = element_line(size=0.25),
        legend.text = element_text(size=28),  
        legend.title = element_blank()
  )
dev.off()

# not needed!!!!!!!
agg <- input_tsm_decon %>%
  aggregate_taxa("Phylum")
agg <- t(data.frame(otu_table(agg)))
df <- cbind(df, agg)
agg <- input_tsm_decon %>%
  aggregate_taxa("Genus")
agg <- t(data.frame(otu_table(agg)))
df <- cbind(df, agg)

# plot
ggplot(df, aes(x=PC1, y=PC2, color=patient)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #geom_text_repel(aes(label=subregion)) +
  #coord_fixed() +
  #ggtitle("Number of species") +
  #stat_ellipse() +
  theme_custom() +
  theme(legend.position = "right",
        panel.grid.major = element_line(size=0.25),
        panel.grid.minor = element_line(size=0.25),
        legend.text = element_text(size=28),  
        legend.title = element_blank()
  )

# split by PC and compare genera
investigate <- df$Id[which(df$PC1<=-15 & df$PC2<=5)]

input <- input %>%
  ps_mutate(
    investigate = ifelse(Id %in% investigate, TRUE, FALSE)
  )

plot <- input %>%
  ps_arrange(desc(patient), desc(region)) %>%
  #subset_samples(patient=="C550") %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 40,
    #tax_level = "Phylum", n_taxa = 5,
    #facet_by = "patient",
    label = "subregion",
    bar_outline_colour = NA,
    sample_order = "default",
    bar_width = 0.7,
    merge_other = TRUE,
    taxon_renamer = toupper
  ) + 
  scale_y_continuous(expand = c(0, 0)) +
  #coord_flip() +
  facet_grid(cols=vars(investigate), scales = "free",space="free") +
  #facet_grid(facets = "patient", scales = "free",space="free" ) +
  theme(axis.text.y = element_text(size=26),
        axis.title.y = element_text(size=26),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=10),
        #strip.background = element_blank(),
        strip.text = element_text(size=26),
        #panel.spacing = unit(.05, "lines"),
        legend.position = "top",
        panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))


### DA on normal v cancer ----

# requires raw (not tsm) values.
#dds_input <- input_raw
dds_input <- input_raw_decon
dds <- phyloseq_to_deseq2(dds_input, ~  patient + group)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))} # calculate geometric means prior to estimate size factors: necessary because lots of zeros
geoMeans <- apply(counts(dds), 1, gm_mean)
deseq2 <- estimateSizeFactors(dds, geoMeans = geoMeans)
deseq2 <- DESeq(deseq2, fitType="local")
plotDispEsts(deseq2)

res = results(deseq2, contrast = c("group","carcinoma","normal"))
res = res[order(res$padj, na.last=NA), ]
sigtab2 = res[(res$padj < 0.05), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(dds_input)[rownames(sigtab2), ], "matrix"))

# plot
sigtabgen = subset(sigtab2, !is.na(Species))

# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$name <- paste(substring(sigtabgen$Genus,1,1),sigtabgen$Species,sep='. ')

jpeg('tempfig.jpeg', width = (30), height = (20), units = "cm", res = 300)
ggplot(sigtabgen, aes(y=reorder(name,log2FoldChange) , x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  ylab(NULL) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  theme_custom() +
  theme(legend.text = element_text(size=24), legend.title = element_text(size=24),
        axis.text.y = element_text(size=24))
dev.off()

### ======DIVERSITY INTRO================================== ----
### number of species ----
x <- data.frame(tax_table(subset_samples(input_tsm_decon, group %in% c("carcinoma","normal"))))

### boxplot of species ----
x <- input_tsm_decon %>%
  aggregate_taxa("Species") %>%
  microbiome::transform(transform = "compositional")

tax <- data.frame(tax_table(x))
x <- data.frame(otu_table(x))
x$genus <- tax[match(rownames(x), tax$unique), 6]
x$species <- tax[match(rownames(x), tax$unique), 7]
keep <- names(head(sort(rowMeans(x[,-which(colnames(x) %in% c("genus","species"))]), decreasing = TRUE), 20))
x <- gather(x[which(x$species %in% keep),], key="key", value="value", -species, -genus)

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(x=reorder(species,-value,median,na.rm=TRUE) , y=(value))) +
  geom_boxplot(aes(fill=species), alpha=0.5) +
  geom_jitter(aes(fill=species), alpha=0.7, shape=21, width = 0.23) +
  xlab(NULL) +
  ylab("Relative abundance \n(species level)") +
  theme_custom() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=22),
        axis.text.y = element_text(size=22),
        axis.title = element_text(size=22),
        legend.position = "none")
dev.off()


jpeg('tempfig.jpeg', width = (40), height = (60), units = "cm", res = 300)
input_tsm_decon %>%
  ps_arrange(desc(patient), desc(region)) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 15,
    label = "subregion",
    bar_outline_colour = NA,
    sample_order = "default",
    bar_width = 0.7,
    merge_other = TRUE,
    taxon_renamer = toupper
  ) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  facet_grid(rows=vars(patient), scales = "free",space="free", switch = "y") +
  theme(axis.text.x = element_text(size=24),
        axis.title.x = element_text(size=26),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(size=22, angle=0),
        strip.background = element_rect(colour="white", fill="white"),
        #panel.spacing = unit(.05, "lines"),
        legend.position = "top",
        legend.text = element_text(size=18),
        legend.title = element_text(size=22),
        panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
dev.off()


### number of species per carcinoma gland ----
y <- data.frame(otu_table(input_tsm_decon))

x <- aggregate_taxa(input_tsm_decon, "Species")
x <- data.frame(otu_table(x))
x <- colSums(x>0)
x <- cbind(sample_data(input_tsm_decon), data.frame(nSpecies=x))

mean(x$nSpecies)
min(x$nSpecies)
max(x$nSpecies)

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(x=reorder(patient, nSpecies, median, na.rm=T), y=nSpecies, fill=patient)) +
  geom_boxplot(alpha=0.6, outlier.size = 0) +
  geom_jitter(shape=21, alpha=0.5, width = 0.2, size=3) +
  xlab(NULL) +
  geom_hline(yintercept = 5, linetype="dashed") +
  scale_y_continuous(name = "Number of species", limits = c(0,140), breaks = c(seq(0,140,20))) +
  theme_custom() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=90, vjust=0.5, size=20),
        panel.grid.major.y = element_line(size=0.5))
dev.off()

x <- aggregate_taxa(input_tsm_decon, "Species")
x <- data.frame(otu_table(x))
x <- colSums(x>0)
x <- cbind(sample_data(input_tsm_decon), data.frame(nSpecies=x), 
           data.frame(shannon=unlist(estimate_richness(input_tsm_decon, measures=c("Shannon")))),
           data.frame(simpson=unlist(estimate_richness(input_tsm_decon, measures=c("Simpson")))))

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(y=nSpecies, x=librarySize)) +
  geom_point(alpha=0.6, shape=21, fill="#000099", size=3) +
  geom_smooth(method='lm', formula= y~x, color="#000099") +
  geom_hline(yintercept = 5, linetype="dashed") +
  ylab("Number of species") +
  xlab("Library size") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "none")
dev.off()

### similar number species across G/B/L ----
nSpecies <- data.frame(Id = sample_data(input_tsm_decon)$Id,
                       patient = sample_data(input_tsm_decon)$patient,
                       type = sample_data(input_tsm_decon)$type,
                       group = sample_data(input_tsm_decon)$group,
                       nSpecies = colSums(data.frame(otu_table(aggregate_taxa(input_tsm_decon, "Species")))>0))

nSpecies$type <- factor(nSpecies$type, levels = c("G","B","L"))
nSpecies <- nSpecies[which(nSpecies$group %in% c("carcinoma","normal")),]
nSpecies$librarySize <- library_sizes[match(nSpecies$Id, library_sizes$Id), 2]
nSpecies$relnspecies <- nSpecies$nSpecies/nSpecies$librarySize

library(viridis)
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(nSpecies, aes(x=type , y=nSpecies,  fill=group)) +
  geom_violin(alpha=0.3, width=1, position = position_dodge(1)) +
  geom_boxplot(alpha=0.6, width=0.3, position = position_dodge(1)) +
  geom_jitter(alpha=0.5, shape=21, 
              position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1) ) +
  stat_compare_means(method = "wilcox.test", size=9, comparisons = list(c("B","L"), c("B","G"), c("L","G"))) +
  stat_compare_means(aes(group = group), method = "wilcox.test", label = "p.format",size=9, label.y = -15) +
  scale_fill_viridis(option = "D", discrete = TRUE) +
  scale_x_discrete(labels = c("Gland","Bulk","Leftover")) +
  xlab(NULL) +
  ylab("Number of species") +
  ylim(-20,220) +
  theme_custom() +
  theme(legend.position = "bottom", legend.text = element_text(size=28),
        legend.title = element_blank())
dev.off()

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(nSpecies, aes(x=type , y=log(relnspecies),  fill=group)) +
  geom_violin(alpha=0.3, width=1, position = position_dodge(1)) +
  geom_boxplot(alpha=0.6, width=0.3, position = position_dodge(1)) +
  geom_jitter(alpha=0.5, shape=21, 
              position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1) ) +
  scale_y_continuous(name="Number of species/Library size", limits = c(-12.5,-3.5),
                     breaks = c(-10,-7.5,-5),
                     labels = c(round(exp(-10),6),round(exp(-7.5),5),round(exp(-5),2))) +
  stat_compare_means(method = "wilcox.test", size=9, comparisons = list(c("B","L"), c("B","G"), c("L","G"))) +
  stat_compare_means(aes(group = group), method = "wilcox.test", label = "p.format",size=9, label.y = -12) +
  scale_fill_viridis(option = "D", discrete = TRUE) +
  scale_x_discrete(labels = c("Gland","Bulk","Leftover")) +
  xlab(NULL) +
  ylab("Number of species") +
  theme_custom() +
  theme(legend.position = "bottom", legend.text = element_text(size=28),
        legend.title = element_blank())
dev.off()

mean(nSpecies$nSpecies[which(nSpecies$type=="G" & nSpecies$group=="carcinoma")], na.rm=T)
mean(nSpecies$nSpecies[which(nSpecies$type=="G" & nSpecies$group=="normal")], na.rm=T)

mean(nSpecies$relnspecies[which(nSpecies$type=="G")], na.rm=T)
mean(nSpecies$relnspecies[which(nSpecies$type=="B")], na.rm=T)

mean(nSpecies$nSpecies[which(nSpecies$type=="B" & nSpecies$group=="carcinoma")], na.rm=T)
mean(nSpecies$nSpecies[which(nSpecies$type=="B" & nSpecies$group=="normal")], na.rm=T)

mean(nSpecies$relnspecies[which(nSpecies$type=="B")], na.rm=T)
mean(nSpecies$relnspecies[which(nSpecies$type=="G")], na.rm=T)

t.test(log(nSpecies$nSpecies[which(nSpecies$type=="B" & nSpecies$group=="carcinoma")]),
       log(nSpecies$nSpecies[which(nSpecies$type=="B" & nSpecies$group=="normal")]))


mean(nSpecies$nSpecies[which(nSpecies$type=="B")], na.rm=T)
mean(nSpecies$nSpecies[which(nSpecies$type=="G")], na.rm=T)
mean(nSpecies$nSpecies[which(nSpecies$type=="L")], na.rm=T)


### glands have biggest library across G/B/L ----

library("viridis") 
x <- sample_data(input_tsm_decon)
x$type <- factor(x$type, levels = c("G","B","L"))
x <- x[which(x$group %in% c("carcinoma","normal")),]

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(x=type , y=log(librarySize),  fill=group)) +
  geom_violin(alpha=0.3, width=1, position = position_dodge(1)) +
  geom_boxplot(alpha=0.6, width=0.3, position = position_dodge(1)) +
  geom_jitter(alpha=0.5, shape=21, 
              position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1) ) +
  scale_y_continuous(name="Library size", limits = c(7.5,17.5),
                     breaks = c(7.5,10,12.5,15,17.5),
                     labels = c(round(exp(7.5),0),round(exp(10),0),round(exp(12.5),0),
                                round(exp(15),0), round(exp(17.5),0))) +
  stat_compare_means(method = "wilcox.test", size=9, comparisons = list(c("B","L"), c("B","G"), c("L","G"))) +
  stat_compare_means(aes(group = group), method = "wilcox.test", label="p.format",size=9, label.y.npc = "bottom") +
  scale_fill_viridis(option = "D", discrete = TRUE) +
  scale_x_discrete(labels = c("Gland","Bulk","Leftover")) +
  xlab(NULL) +
  ylab("Library size") +
  theme_custom() +
  theme(legend.position = "bottom", legend.text = element_text(size=28),
        legend.title = element_blank())
dev.off()


mean((x$librarySize[which(x$type=="B")]), na.rm=T)
mean((x$librarySize[which(x$type=="G")]), na.rm=T)
mean(log(x$librarySize[which(x$type=="L")]), na.rm=T)

t.test(x$librarySize[which(x$type=="G")], x$librarySize[which(x$type=="B")])



### subset to remove samples with species <5 ----
x <- aggregate_taxa(input_tsm_decon, "Species")
x <- data.frame(otu_table(x))
x <- colSums(x>0)
x <- cbind(sample_data(input_tsm_decon), data.frame(nSpecies=x))

keep <- rownames(x)[which(x$nSpecies>5)]
keep <- substr(keep, 7, nchar(keep)-23)

input_tsm_decon_rem <- subset_samples(input_tsm_decon, Id %in% keep)
input_raw_decon_rem <- subset_samples(input_raw_decon, Id %in% keep)


### similar alpha diversity  across G/B/L ----

# dont need tsm norm for alpha diverstiy, but will use because has 2removed samples
shannon_type <- data.frame(Id = sample_data(input_tsm_decon)$Id,
                           patient = sample_data(input_tsm_decon)$patient,
                           type = sample_data(input_tsm_decon)$type,
                           group = sample_data(input_tsm_decon)$group,
                           shannon = unlist(estimate_richness(input_tsm_decon, measures=c("Shannon"))))

shannon_type$type <- factor(shannon_type$type, levels = c("G","B","L"))
shannon_type <- shannon_type[which(shannon_type$group %in% c("carcinoma","normal")),]

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(shannon_type, aes(x=type , y=shannon,  fill=group)) +
  geom_violin(alpha=0.3, width=1, position = position_dodge(1)) +
  geom_boxplot(alpha=0.6, width=0.3, position = position_dodge(1)) +
  geom_jitter(alpha=0.5, shape=21, 
              position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1) ) +
  stat_compare_means(method = "wilcox.test", size=9, comparisons = list(c("B","L"), c("B","G"), c("L","G"))) +
  stat_compare_means(aes(group = group), method = "wilcox.test", label="p.format",size=9, label.y.npc = "bottom") +
  scale_fill_viridis(option = "D", discrete = TRUE) +
  scale_x_discrete(labels = c("Gland","Bulk","Leftover")) +
  ylim(0,5) +
  xlab(NULL) +
  ylab("Shannon diversity") +
  theme_custom() +
  theme(legend.position = "bottom", legend.text = element_text(size=28),
        legend.title = element_blank())
dev.off()


mean((shannon_type$shannon[which(shannon_type$type=="B")]), na.rm=T)
mean((shannon_type$shannon[which(shannon_type$type=="G")]), na.rm=T)
mean((shannon_type$shannon[which(shannon_type$type=="L")]), na.rm=T)

t.test(shannon_type$shannon[which(shannon_type$type=="B" & shannon_type$group=="carcinoma")],
shannon_type$shannon[which(shannon_type$type=="B" & shannon_type$group=="normal")])

#simpsons
simpson_type <- data.frame(Id = sample_data(input_tsm_decon)$Id,
                           patient = sample_data(input_tsm_decon)$patient,
                           type = sample_data(input_tsm_decon)$type,
                           group = sample_data(input_tsm_decon)$group,
                           simpson = unlist(estimate_richness(input_tsm_decon, measures=c("Simpson"))))

simpson_type$type <- factor(simpson_type$type, levels = c("G","B","L"))
simpson_type <- simpson_type[which(simpson_type$group %in% c("carcinoma","normal")),]

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(simpson_type, aes(x=type , y=simpson,  fill=group)) +
  geom_violin(alpha=0.3, width=1, position = position_dodge(1)) +
  geom_boxplot(alpha=0.6, width=0.3, position = position_dodge(1)) +
  geom_jitter(alpha=0.5, shape=21, 
              position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1) ) +
  stat_compare_means(method = "wilcox.test", size=9, comparisons = list(c("B","L"), c("B","G"), c("L","G"))) +
  stat_compare_means(aes(group = group), method = "wilcox.test", label="p.format",size=9, label.y.npc = "bottom") +
  scale_fill_viridis(option = "D", discrete = TRUE) +
  scale_x_discrete(labels = c("Gland","Bulk","Leftover")) +
  xlab(NULL) +
  ylim(0,1.3) +
  ylab("Simpson diversity") +
  theme_custom() +
  theme(legend.position = "bottom", legend.text = element_text(size=28),
        legend.title = element_blank())
dev.off()

mean((simpson_type$simpson[which(simpson_type$type=="B" & simpson_type$group=="normal")]), na.rm=T)
mean((simpson_type$simpson[which(simpson_type$type=="B" & simpson_type$group=="carcinoma")]), na.rm=T)


mod <- lm(shannon ~ nSpecies + librarySize, data=x)
summary(mod)

#pielou
pielou_type <- data.frame(Id = sample_data(input_tsm_decon)$Id,
                           patient = sample_data(input_tsm_decon)$patient,
                           type = sample_data(input_tsm_decon)$type,
                           group = sample_data(input_tsm_decon)$group,
                          pielou = evenness(input_tsm_decon, 'pielou'))


pielou_type$type <- factor(pielou_type$type, levels = c("G","B","L"))
pielou_type <- pielou_type[which(pielou_type$group %in% c("carcinoma","normal")),]

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(pielou_type, aes(x=type , y=pielou,  fill=group)) +
  geom_violin(alpha=0.3, width=1, position = position_dodge(1)) +
  geom_boxplot(alpha=0.6, width=0.3, position = position_dodge(1)) +
  geom_jitter(alpha=0.5, shape=21, 
              position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1) ) +
  stat_compare_means(method = "wilcox.test", size=9, comparisons = list(c("B","L"), c("B","G"), c("L","G"))) +
  stat_compare_means(aes(group = group), method = "wilcox.test", label="p.format",size=9, label.y.npc = "bottom") +
  scale_fill_viridis(option = "D", discrete = TRUE) +
  scale_x_discrete(labels = c("Gland","Bulk","Leftover")) +
  xlab(NULL) +
  ylim(0,1.3) +
  ylab("Pielou evenness") +
  theme_custom() +
  theme(legend.position = "bottom", legend.text = element_text(size=28),
        legend.title = element_blank())
dev.off()
### correlate number species, library size, shannon ----

x <- aggregate_taxa(input_tsm_decon_rem, "Species")
x <- data.frame(otu_table(x))
x <- colSums(x>0)
x <- cbind(sample_data(input_tsm_decon_rem), data.frame(nSpecies=x), 
           data.frame(shannon=unlist(estimate_richness(input_tsm_decon_rem, measures=c("Shannon")))),
           data.frame(simpson=unlist(estimate_richness(input_tsm_decon_rem, measures=c("Simpson")))))

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(y=simpson, x=nSpecies)) +
  geom_point(alpha=0.6, shape=21, fill="#000099",size=3) +
  geom_smooth(method='lm', formula= y~x, color="#000099") +
  xlab("Number of species") +
  ylab("Simpson diversity") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'right', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "none", plot.margin = margin(8,40,1,1))
dev.off()

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(y=shannon, x=librarySize)) +
  geom_point(alpha=0.6, shape=21, fill="#000099",size=3) +
  geom_smooth(method='lm', formula= y~x, color="#000099") +
  xlab("Library size") +
  ylab("Shannon diversity") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "none", plot.margin = margin(1,40,1,1))
dev.off()

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(y=simpson, x=librarySize)) +
  geom_point(alpha=0.6, shape=21, fill="#000099", size=3) +
  geom_smooth(method='lm', formula= y~x, color="#000099") +
  xlab("Library size") +
  ylab("Simpsons diversity") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'right', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "none", plot.margin = margin(10,40,1,1))
dev.off()

ggplot(x, aes(y=simpson, x=shannon)) +
  geom_point(alpha=0.6, shape=21, fill="#000099") +
  geom_smooth(method='lm', formula= y~x, color="#000099") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "none")

### phyla/genera by region (bar chart) (x) ----

# phylum
input_tsm_decon_rem %>%
  aggregate_top_taxa(level = "Phylum",top=4) %>%
  microbiome::transform(transform = "compositional") %>%
  ps_arrange(desc(patient), desc(region)) %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 4,
    #facet_by = "patient",
    label = "sample",
    bar_outline_colour = NA,
    sample_order = "default",
    bar_width = 0.7,
    merge_other = TRUE,
    taxon_renamer = toupper
  ) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  facet_grid(facets = "patient", scales = "free",space="free" ) +
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

# genus
input_tsm_decon_rem %>%
  aggregate_top_taxa(level = "Genus",top=20) %>%
  microbiome::transform(transform = "compositional") %>%
  ps_arrange(desc(patient), desc(region)) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 20,
    #facet_by = "patient",
    label = "sample",
    bar_outline_colour = NA,
    sample_order = "default",
    bar_width = 0.7,
    merge_other = TRUE,
    taxon_renamer = toupper
  ) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  facet_grid(facets = "patient", scales = "free",space="free" ) +
  theme(axis.text = element_text(size=10),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))



prevdf <- apply(X = otu_table(input_tsm_decon_rem),
                MARGIN = ifelse(taxa_are_rows(input_tsm_decon_rem), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(input_tsm_decon_rem),
                     tax_table(input_tsm_decon_rem))
prevdf$avePersample <- prevdf$TotalAbundance/prevdf$Prevalence
prevdf <- prevdf[which(prevdf$Genus %in% c("Fusobacterium","Aquabacterium","Escherichia","Delftia","Enterbacter","Klebsiella","Pseudomonas")),]
prevdf1 <- subset(prevdf, Genus %in% get_taxa_unique(input_tsm_decon_rem, "Genus"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(input_tsm_decon_rem),color=Genus)) +
  geom_hline(yintercept = 0.1, alpha = 0.5, linetype = 2) +  geom_point(size = 1, alpha = 0.7) + # Include a guess for parameter
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")

### DA sample type ----
dds_input <- input_tsm_decon_rem_car %>%
  aggregate_taxa("Genus")

sample_data(dds_input)$MSI <- unlist(metadata_tsm.G_patient[match(sample_data(dds_input)$patient, metadata_tsm.G_patient$patient), 5])

dds_input <- subset_samples(dds_input, type!="L")

dds <- phyloseq_to_deseq2(dds_input, ~  type)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))} # calculate geometric means prior to estimate size factors: necessary because lots of zeros
geoMeans <- apply(counts(dds), 1, gm_mean)
deseq2 <- estimateSizeFactors(dds, geoMeans = geoMeans)
deseq2 <- DESeq(deseq2, fitType="local")
plotDispEsts(deseq2)

res = results(deseq2, contrast = c("type","G","B")) #B is base
res = res[order(res$padj, na.last=NA), ]
sigtab2 = res[(res$padj < 0.05), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(dds_input)[rownames(sigtab2), ], "matrix"))

# plot
sigtabgen = subset(sigtab2, !is.na(Species))

# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$name <- paste(substring(sigtabgen$Genus,1,1),sigtabgen$Species,sep='. ')


ggplot(sigtabgen, aes(y=reorder(name,log2FoldChange) , x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  ylab(NULL) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  theme_custom() +
  theme(legend.text = element_text(size=24), legend.title = element_text(size=24),
        axis.text.y = element_text(size=20))

### DA on MSS v MSI ----
dds_input <- input_tsm_decon_rem_car

sample_data(dds_input)$MSI <- unlist(metadata_tsm.G_patient[match(sample_data(dds_input)$patient, metadata_tsm.G_patient$patient), 5])

dds_input <- subset_samples(dds_input, !is.na(MSI))

dds <- phyloseq_to_deseq2(dds_input, ~  MSI)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))} # calculate geometric means prior to estimate size factors: necessary because lots of zeros
geoMeans <- apply(counts(dds), 1, gm_mean)
deseq2 <- estimateSizeFactors(dds, geoMeans = geoMeans)
deseq2 <- DESeq(deseq2, fitType="local")
plotDispEsts(deseq2)

res = results(deseq2, contrast = c("MSI","MSS","MSI")) #MSI is base
res = res[order(res$padj, na.last=NA), ]
sigtab2 = res[(res$padj < 0.05), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(dds_input)[rownames(sigtab2), ], "matrix"))

# plot
sigtabgen = subset(sigtab2, !is.na(Species))

# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$name <- paste(substring(sigtabgen$Genus,1,1),sigtabgen$Species,sep='. ')


ggplot(sigtabgen, aes(y=reorder(name,log2FoldChange) , x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  ylab(NULL) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  theme_custom() +
  theme(legend.text = element_text(size=24), legend.title = element_text(size=24),
        axis.text.y = element_text(size=20))




### ======METADATA FOR CARCINOMA GLANDS==================== ----
### subset carcinoma  ----
input_tsm_decon_rem_car <- subset_samples(input_tsm_decon_rem, group=="carcinoma")
input_raw_decon_rem_car <- subset_samples(input_raw_decon_rem, group=="carcinoma")

### subset carcinoma glands ----
input_tsm_decon_rem_carG <- subset_samples(input_tsm_decon_rem, group=="carcinoma" & type=="G")
input_raw_decon_rem_carG <- subset_samples(input_raw_decon_rem, group=="carcinoma" & type=="G")

### pull data per sample ----
bacteria_tsm.G <- input_tsm_decon_rem_car

# pull data on CIN per sample
dwgs.perpatient <- readRDS("~/Documents/SCAA/Data/dwgs.perpatient.rds")
dwgs.ploidyRecentre.perpatient <- readRDS("~/Documents/SCAA/Data/dwgs.ploidyRecentre.perpatient.rds")
ploidy.perpatient <- readRDS("~/Documents/SCAA/Data/ploidy.perpatient.rds")
ploidy.persample <- readRDS("~/Documents/SCAA/Data/ploidy.persample.rds")
totalSNVcount <- read.delim("~/Documents/SCAA/Data/total_number_of_mutations.txt")
scaaPP <- readRDS("~/Documents/SCAA/Data/scaaPP.rds")

metadata_tsm.G_sample <- data.frame(sample_data(bacteria_tsm.G))
metadata_tsm.G_sample$ploidy <- unlist(ploidy.persample[match(metadata_tsm.G_sample$sample, gsub("\\.", "_", ploidy.persample$sample)), 3])
metadata_tsm.G_sample$WGD <- ifelse(metadata_tsm.G_sample$ploidy>2.5, TRUE, FALSE)
metadata_tsm.G_sample$MSI <- ifelse(metadata_tsm.G_sample$patient %in% c("C552","C548","C518","C562","C516","C536"), "MSI", "MSS")
metadata_tsm.G_sample$totalSNVperpatient <- totalSNVcount[match(metadata_tsm.G_sample$patient, totalSNVcount$Patient), 2]

# calculate PGA per sample and average for each patient (no recentering)
persample <- list()
for (i in 1:length(dwgs.perpatient) ) {
  wd <- dwgs.perpatient[[i]]
  wd[,-c(1:3)] <- round(wd[,-c(1:3)])
  weights <- wd$stop-wd$start
  weights <- weights / sum(weights)
  
  if (ncol(wd)==4) {
    y <- weights * ifelse(wd[4]==2, 1, 0)
    persample[[i]] <- mean(y)
  } else {
    y <- apply(wd[,-c(1:3)], 2, function(x) {
      weights * ifelse(x!=2, 1, 0)
    })
    persample[[i]] <- colSums(y, na.rm = T)
  }
}

pga <- data.frame(pga = unlist(persample))
metadata_tsm.G_sample$pga <- pga[match(metadata_tsm.G_sample$sample, gsub("\\.", "_", rownames(pga))), 1]

# # calculate PGA per sample and average for each patient (with recentering)
persample <- list()
for (i in 1:length(dwgs.ploidyRecentre.perpatient) ) {
  wd <- dwgs.ploidyRecentre.perpatient[[i]]
  wd[,-c(1:3)] <- round(wd[,-c(1:3)])
  weights <- wd$stop-wd$start
  weights <- weights / sum(weights)
  
  if (ncol(wd)==4) {
    y <- weights * ifelse(wd[4]==2, 1, 0)
    persample[[i]] <- mean(y)
  } else {
    y <- apply(wd[,-c(1:3)], 2, function(x) {
      weights * ifelse(x!=2, 1, 0)
    })
    persample[[i]] <- colSums(y, na.rm = T)
  }
}

pga_recentre <- data.frame(pga_recentre = unlist(persample))
metadata_tsm.G_sample$pga_recentre <- pga_recentre[match(metadata_tsm.G_sample$sample, gsub("\\.", "_", rownames(pga_recentre))), 1]

# calculate number of segments pp
persample <- list()
for (i in 1:length(dwgs.perpatient) ) {
  wd <- dwgs.perpatient[[i]]
  x <- apply(wd[-c(1:3)],2, function(x) {length((rle(x))[[1]])})
  persample[[i]] <- setNames(x, colnames(wd)[-c(1:3)])
}

nSegments <- data.frame(nSegments = unlist(persample))
metadata_tsm.G_sample$nSegments <- nSegments[match(metadata_tsm.G_sample$sample, gsub("\\.", "_", rownames(nSegments))), 1]

# add in alpha diversity metric
metadata_tsm.G_sample$shannon = unlist(estimate_richness(bacteria_tsm.G, measures=c("Shannon")))
metadata_tsm.G_sample$simpson = unlist(estimate_richness(bacteria_tsm.G, measures=c("Simpson")))

# pull quanTIseq data per sample
quanTIseq_sample <- read.table("~/Documents/Microbial/data/EPICC_quanTIseq_cell_type_proportions.txt", header = TRUE)
qc <- read.table("~/Documents/Microbial/data/ListRNAPass.EPICC.txt")
quanTIseq_sample <- quanTIseq_sample[which(quanTIseq_sample$Sample %in% qc$V1),]
quanTIseq_sample$patient <- str_before_nth(quanTIseq_sample$Sample,"_",1)
quanTIseq_sample$region <- substr(str_before_nth(quanTIseq_sample$Sample,"_",2), 1, 6)
quanTIseq_sample$subregion <- str_before_nth(quanTIseq_sample$Sample,"_",2)
quanTIseq_sample <- quanTIseq_sample[c(13:15,1:12)]
colnames(quanTIseq_sample)[4] <- "sample"

# merge with immune
metadata_tsm.G_sample <- merge(metadata_tsm.G_sample, quanTIseq_sample, all = TRUE, by = c("patient","region","subregion","sample"))
metadata_tsm.G_sample$lymphocytes <- metadata_tsm.G_sample$B.cells + metadata_tsm.G_sample$NK.cells + metadata_tsm.G_sample$T.cells.CD4 +
  metadata_tsm.G_sample$T.cells.CD8 + metadata_tsm.G_sample$Tregs

# add in two columns to ID if immune/microbiome data available
metadata_tsm.G_sample$microbiomeData <- ifelse(!is.na(metadata_tsm.G_sample$librarySize), TRUE, FALSE)
metadata_tsm.G_sample$quanTIseqData <- ifelse(!is.na(metadata_tsm.G_sample$Other), TRUE, FALSE)

# add firmicutes/bacteroidetes ratio
FBratio <- microbiome::aggregate_top_taxa(bacteria_tsm.G, "Phylum", top = 19) %>%
  microbiome::transform(transform = "compositional")
FBratio <- cbind(sample_data(FBratio), data.frame(t(otu_table(FBratio))))
FBratio$FBratio <- FBratio$Firmicutes/FBratio$Bacteroidetes
metadata_tsm.G_sample$FBratio <- FBratio[match(metadata_tsm.G_sample$Id, FBratio$Id), 15]

### pull data per region ----

# collapse into per region
metadata_tsm.G_region <- data.frame(region = unique(metadata_tsm.G_sample$region),
                                    librarySize = unlist(aggregate(metadata_tsm.G_sample$librarySize, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2]),
                                    ploidy = unlist(aggregate(metadata_tsm.G_sample$ploidy, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2]))
metadata_tsm.G_region$patient <- str_before_nth(metadata_tsm.G_region$region, "_", 1)
metadata_tsm.G_region <- metadata_tsm.G_region[c(4,1:3)]
metadata_tsm.G_region$WGD <- ifelse(metadata_tsm.G_region$ploidy>2.5, TRUE, FALSE)
metadata_tsm.G_region$MSI <- ifelse(metadata_tsm.G_region$patient %in% c("C552","C548","C518","C562","C516","C536"), "MSI", "MSS")
metadata_tsm.G_region$totalSNVperpatient <- unlist(aggregate(metadata_tsm.G_sample$totalSNVperpatient, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_region$pga <- unlist(aggregate(metadata_tsm.G_sample$pga, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_region$nSegments <- unlist(aggregate(metadata_tsm.G_sample$nSegments, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2])

# add alpha
metadata_tsm.G_region$shannon <- unlist(aggregate(metadata_tsm.G_sample$shannon, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_region$simpson <- unlist(aggregate(metadata_tsm.G_sample$simpson, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2]) 

# add immune
metadata_tsm.G_region$Bcells <- unlist(aggregate(metadata_tsm.G_sample$B.cells, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_region$M1 <- unlist(aggregate(metadata_tsm.G_sample$Macrophages.M1, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_region$M2 <- unlist(aggregate(metadata_tsm.G_sample$Macrophages.M2, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_region$Mono <- unlist(aggregate(metadata_tsm.G_sample$Monocytes, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_region$Neutro <- unlist(aggregate(metadata_tsm.G_sample$Neutrophils, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_region$NK <- unlist(aggregate(metadata_tsm.G_sample$NK.cells, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_region$CD4 <- unlist(aggregate(metadata_tsm.G_sample$T.cells.CD4, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_region$CD8 <- unlist(aggregate(metadata_tsm.G_sample$T.cells.CD8, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_region$Treg <- unlist(aggregate(metadata_tsm.G_sample$Tregs, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_region$Dend <- unlist(aggregate(metadata_tsm.G_sample$Dendritic.cells, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2]) 
metadata_tsm.G_region$Other <- unlist(aggregate(metadata_tsm.G_sample$Other, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2]) 
metadata_tsm.G_region$lymphocytes <- unlist(aggregate(metadata_tsm.G_sample$lymphocytes, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2]) 
metadata_tsm.G_region$NLR <- metadata_tsm.G_region$Neutro / (metadata_tsm.G_region$Bcells+metadata_tsm.G_region$NK+
                                                               metadata_tsm.G_region$CD4+metadata_tsm.G_region$CD8+
                                                               metadata_tsm.G_region$Treg)

# add in two columns to ID if immune/microbiome data available
metadata_tsm.G_region$microbiomeData <- ifelse(!is.na(metadata_tsm.G_region$librarySize), TRUE, FALSE)
metadata_tsm.G_region$quanTIseqData <- ifelse(!is.na(metadata_tsm.G_region$Other), TRUE, FALSE)

# add average firmicutes/bacteroidetes ratio
metadata_tsm.G_region$FBratio <- unlist(aggregate(metadata_tsm.G_sample$FBratio, list(metadata_tsm.G_sample$region), FUN=mean, na.rm =TRUE)[2]) 

# grab all pairwise bray curtis measure within a patient
braycurtis_dist <- phyloseq::distance(bacteria_tsm.G, method="bray")
braycurtis_dist <- as.matrix(braycurtis_dist)

sub_dist <- list()
regions_all <- sample_data(bacteria_tsm.G)$region

for (i in 1:length(unique(regions_all)) ) { 
  region <- unique(regions_all)[i]
  index <- which(regions_all == region)
  sample_group <- sample_names(bacteria_tsm.G)[index]
  sub_dist[[region]] <- braycurtis_dist[sample_group, sample_group]
  sub_dist[[region]][!lower.tri(sub_dist[[region]])] <- NA
}

braygroups <- melt(sub_dist)
df.bray_region <- braygroups[complete.cases(braygroups), ]
df.bray_region$region <- factor(df.bray_region$L1, levels=names(sub_dist))

# make DISsimilarity
#df.bray_region$value <- 1-df.bray_region$value

# bray is between individual samples, calculate average distances within group
x <- aggregate(df.bray_region$value, list(df.bray_region$region), mean)
metadata_tsm.G_region$braycurtis <- x[match(metadata_tsm.G_region$region, x$Group.1), 2]


### pull data per patient ----

# collapse into per patient
metadata_tsm.G_patient <- data.frame(patient = unique(metadata_tsm.G_sample$patient),
                                     librarySize = unlist(aggregate(metadata_tsm.G_sample$librarySize, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2]),
                                     ploidy = unlist(aggregate(metadata_tsm.G_sample$ploidy, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2]))
metadata_tsm.G_patient$WGD <- ifelse(metadata_tsm.G_patient$ploidy>2.5, TRUE, FALSE)
metadata_tsm.G_patient$MSI <- ifelse(metadata_tsm.G_patient$patient %in% c("C552","C548","C518","C562","C516","C536"), "MSI", "MSS")
metadata_tsm.G_patient$totalSNVperpatient <- unlist(aggregate(metadata_tsm.G_sample$totalSNVperpatient, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_patient$pga <- unlist(aggregate(metadata_tsm.G_sample$pga, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_patient$nSegments <- unlist(aggregate(metadata_tsm.G_sample$nSegments, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])

# add alpha
metadata_tsm.G_patient$shannon <- unlist(aggregate(metadata_tsm.G_sample$shannon, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_patient$simpson <- unlist(aggregate(metadata_tsm.G_sample$simpson, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2]) 

# add immune
metadata_tsm.G_patient$Bcells <- unlist(aggregate(metadata_tsm.G_sample$B.cells, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_patient$M1 <- unlist(aggregate(metadata_tsm.G_sample$Macrophages.M1, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_patient$M2 <- unlist(aggregate(metadata_tsm.G_sample$Macrophages.M2, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_patient$Mono <- unlist(aggregate(metadata_tsm.G_sample$Monocytes, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_patient$Neutro <- unlist(aggregate(metadata_tsm.G_sample$Neutrophils, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_patient$NKP <- unlist(aggregate(metadata_tsm.G_sample$NK.cells, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_patient$CD4 <- unlist(aggregate(metadata_tsm.G_sample$T.cells.CD4, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_patient$CD8 <- unlist(aggregate(metadata_tsm.G_sample$T.cells.CD8, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_patient$Treg <- unlist(aggregate(metadata_tsm.G_sample$Tregs, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])
metadata_tsm.G_patient$Dend <- unlist(aggregate(metadata_tsm.G_sample$Dendritic.cells, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])  
metadata_tsm.G_patient$lymphocytes <- unlist(aggregate(metadata_tsm.G_sample$lymphocytes, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2])  
metadata_tsm.G_patient$Other <- unlist(aggregate(metadata_tsm.G_sample$Other, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2]) 

# add in two columns to ID if immune/microbiome data available
metadata_tsm.G_patient$microbiomeData <- ifelse(!is.na(metadata_tsm.G_patient$librarySize), TRUE, FALSE)
metadata_tsm.G_patient$quanTIseqData <- ifelse(!is.na(metadata_tsm.G_patient$Other), TRUE, FALSE)

# add average firmicutes/bacteroidetes ratio
metadata_tsm.G_patient$FBratio <- unlist(aggregate(metadata_tsm.G_sample$FBratio, list(metadata_tsm.G_sample$patient), FUN=mean, na.rm =TRUE)[2]) 

# grab all pairwise bray curtis measure within a patient
braycurtis_dist <- phyloseq::distance(bacteria_tsm.G, method="bray")
braycurtis_dist <- as.matrix(braycurtis_dist)

sub_dist <- list()
patient_all <- sample_data(bacteria_tsm.G)$patient

for (i in 1:length(unique(patient_all)) ) { 
  patient <- unique(patient_all)[i]
  index <- which(patient_all == patient)
  sample_group <- sample_names(bacteria_tsm.G)[index]
  sub_dist[[patient]] <- braycurtis_dist[sample_group, sample_group]
  sub_dist[[patient]][!lower.tri(sub_dist[[patient]])] <- NA
}

braygroups <- melt(sub_dist)
df.bray_patient <- braygroups[complete.cases(braygroups), ]
df.bray_patient$patient <- factor(df.bray_patient$L1, levels=names(sub_dist))

# make DISsimilarity
#df.bray_patient$value <- 1-df.bray_patient$value

# bray is between individual samples, calculate average distances within group
x <- aggregate(df.bray_patient$value, list(df.bray_patient$patient), mean)
metadata_tsm.G_patient$braycurtis <- x[match(metadata_tsm.G_patient$patient, x$Group.1), 2]

# calculate mean distance between patient and other patients (with sampling)
mean_dist <- list()
patient_all <- sample_data(bacteria_tsm.G)$patient

for (i in 1:length(unique(patient_all)) ) { 
  patient <- unique(patient_all)[i]
  index <- which(patient_all == patient)
  antiIndex <- which(patient_all != patient)
  sample_group <- sample_names(bacteria_tsm.G)[index]
  other_group <- sample_names(bacteria_tsm.G)[antiIndex]
  dist <- braycurtis_dist[sample_group, other_group]
  dist[!lower.tri(dist)] <- NA
  mean_dist[[i]] <- mean(dist, na.rm=TRUE)
}

x <- data.frame(patient=unique(patient_all), braycurtis.allOther=unlist(mean_dist))
metadata_tsm.G_patient$braycurtis.allOther <- x[match(metadata_tsm.G_patient$patient, x$patient), 2]


### ======ITH OF ALPHA DIVERSITY=========================== ----
### alpha diversity (per sample, by patient) ----

mean(metadata_tsm.G_sample$shannon[which(metadata_tsm.G_sample$microbiomeData==TRUE)])
min(metadata_tsm.G_sample$shannon[which(metadata_tsm.G_sample$microbiomeData==TRUE)])
max(metadata_tsm.G_sample$shannon[which(metadata_tsm.G_sample$microbiomeData==TRUE)])

#coefficient of varaicen per patient
output <- list()
for ( i in 1:length(unique(metadata_tsm.G_sample$patient)) ) {
  wd <- metadata_tsm.G_sample$shannon[which(metadata_tsm.G_sample$patient==unique(metadata_tsm.G_sample$patient)[i])]
  wd <- na.omit(wd)
  output[[i]] <- data.frame(patient = unique(metadata_tsm.G_sample$patient)[i], 
                            Variance = var(wd), 
                            St.dev = sd(wd), 
                            mean = mean(wd))
  output[[i]]$COV <- output[[i]]$St.dev / output[[i]]$mean
}

output <- do.call(rbind, output)
metadata_tsm.G_sample$COV <- output[match(metadata_tsm.G_sample$patient, output$patient),5]

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(metadata_tsm.G_sample[which(metadata_tsm.G_sample$microbiomeData==TRUE),], aes(x=reorder(patient, shannon, median, na.rm=T) , y=shannon,  fill=patient)) +
  #geom_bar(aes(y=COV), stat = "summary", fun = "mean", alpha=0.5) + 
  geom_boxplot(alpha=0.6, outlier.size = 0) +
  geom_jitter(shape=21, alpha=0.5, width = 0.2, size=3) +
  xlab(NULL) +
  #scale_y_continuous(name = "Shannon diversity",
   #                  sec.axis = sec_axis(~., name="Second Axis")) +
  theme_custom() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=90, vjust=0.5, size=20),
        panel.grid.major.y = element_line(size=0.5))
dev.off()


ggplot(output, aes(x=reorder(patient, COV, median, na.rm=T) , y=COV,  fill=patient)) +
  geom_bar(stat="identity")
  
mean(metadata_tsm.G_sample$simpson[which(metadata_tsm.G_sample$microbiomeData==TRUE)])
min(metadata_tsm.G_sample$simpson[which(metadata_tsm.G_sample$microbiomeData==TRUE)])
max(metadata_tsm.G_sample$simpson[which(metadata_tsm.G_sample$microbiomeData==TRUE)])

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(metadata_tsm.G_sample[which(metadata_tsm.G_sample$microbiomeData==TRUE),], aes(x=reorder(patient, simpson, median, na.rm=T) , y=simpson,  fill=patient)) +
  geom_boxplot(alpha=0.6, outlier.size = 0) +
  geom_jitter(shape=21, alpha=0.5, width = 0.2, size=3) +
  xlab(NULL) +
  scale_y_continuous(name = "Simpson's diversity") +
  theme_custom() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=90, vjust=0.5, size=20),
        panel.grid.major.y = element_line(size=0.5))
dev.off()




### alpha diversity by WGD ----

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(metadata_tsm.G_sample[!is.na(metadata_tsm.G_sample$WGD),], aes(x=WGD , y=simpson,  fill=WGD)) +
  geom_violin(alpha=0.5) +
  geom_boxplot(alpha=0.6, outlier.size = 0, width=0.5) +
  geom_jitter(size = 3, shape=21, alpha=0.8, width = 0.19, height = 0) +
  xlab(NULL) +
  ylim(0.4,1.01) +
  ylab("Simpson's diversity") +
  scale_fill_manual(values =  c("#F66666", "#009900")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("TRUE","FALSE")), size=9) +
  scale_x_discrete(labels = c("diploid","WGD")) +
  theme_custom() +
  theme(axis.text.x = element_text(), legend.position = "none")
dev.off()

t.test(metadata_tsm.G_sample$simpson[which(metadata_tsm.G_sample$WGD==TRUE)],
       metadata_tsm.G_sample$simpson[which(metadata_tsm.G_sample$WGD==FALSE)])

### alpha diversity by MSI ----

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(metadata_tsm.G_sample[!is.na(metadata_tsm.G_sample$MSI),], aes(x=MSI , y=shannon,  fill=MSI)) +
  geom_violin(alpha=0.5) +
  geom_boxplot(alpha=0.6, outlier.size = 0, width=0.5) +
  geom_jitter(size = 3, shape=21, alpha=0.8, width = 0.19, height = 0) +
  xlab(NULL) +
  #ylim(0,1.01) +
  ylab("Shannon diversity") +
  scale_fill_manual(values =  c("#339999", "#CC6600")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("MSI","MSS")), size=9) +
  theme_custom() +
  theme(axis.text.x = element_text(), legend.position = "none")
dev.off()

  
t.test(metadata_tsm.G_sample$shannon[which(metadata_tsm.G_sample$WGD==TRUE)],
       metadata_tsm.G_sample$shannon[which(metadata_tsm.G_sample$WGD==FALSE)])

t.test(metadata_tsm.G_sample$simpson[which(metadata_tsm.G_sample$MSI=="MSI")],
       metadata_tsm.G_sample$simpson[which(metadata_tsm.G_sample$MSI=="MSS")])




### alpha diversity range per patient ----

x <- data.frame(aggregate(metadata_tsm.G_sample$shannon, list(metadata_tsm.G_sample$patient), FUN=function(x) c(mn = mean(x, na.rm=TRUE), n = sd(x, na.rm=TRUE) )))
x <- do.call(data.frame, x)
colnames(x) <- c("patient","mean","stDev")
x$stdevMean <- x$stDev/x$mean
x$MSI <- metadata_tsm.G_patient[match(x$patient, metadata_tsm.G_patient$patient), 5]
ggplot(x, aes(x=reorder(patient,stdevMean,median) , y=stdevMean,  fill=MSI)) +
  geom_point(shape=21)


x <- data.frame(aggregate(metadata_tsm.G_sample$simpson, list(metadata_tsm.G_sample$patient), FUN=function(x) c(mn = mean(x, na.rm=TRUE), n = sd(x, na.rm=TRUE) )))
x <- do.call(data.frame, x)
colnames(x) <- c("patient","mean","stDev")
x$stdevMean <- x$stDev/x$mean
x$MSI <- metadata_tsm.G_patient[match(x$patient, metadata_tsm.G_patient$patient), 5]
ggplot(x, aes(x=reorder(patient,stdevMean,median) , y=stdevMean,  fill=MSI)) +
  geom_point(shape=21)

### ======ALPHA DIVERSITY AND CIN========================== ----
### association between alpha diversity and PGA (per sample)----
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(metadata_tsm.G_patient[!is.na(metadata_tsm.G_patient$shannon),], aes(x=(pga) , y=(shannon), color=MSI)) +
  geom_smooth(method='lm', formula= y~x, aes(color=MSI, fill=MSI)) +
  geom_point(shape=21, aes(fill=MSI), color="black", size=4, alpha=0.8) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  scale_color_manual(values =  c("#339999", "#CC6600")) +
  scale_fill_manual(values =  c("#339999", "#CC6600")) +
  ylab("Average Shannon diversity per patient") +
  xlab("Proportion of Genome Altered") +
  theme_custom() +
  theme(legend.text = element_text(size=28), legend.title = element_blank(), legend.position = "top")
dev.off()

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(metadata_tsm.G_patient[!is.na(metadata_tsm.G_patient$simpson),], aes(x=(pga) , y=(simpson), color=MSI)) +
  geom_smooth(method='lm', formula= y~x, aes(color=MSI, fill=MSI)) +
  geom_point(shape=21, aes(fill=MSI), color="black", size=4, alpha=0.8) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  scale_color_manual(values =  c("#339999", "#CC6600")) +
  scale_fill_manual(values =  c("#339999", "#CC6600")) +
  ylab("Average Simpson's diversity per patient") +
  xlab("Proportion of Genome Altered") +
  theme_custom() +
  theme(legend.text = element_text(size=28), legend.title = element_blank(), legend.position = "top")
dev.off()

ggplot(metadata_tsm.G_sample, aes(x=(pga_recentre) , y=(shannon), color=MSI)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

jpeg('tempfig.jpeg', width = (20), height = (170), units = "cm", res = 300)
ggplot(metadata_tsm.G_sample, aes(x=(pga) , y=(shannon), color=MSI)) +
  geom_smooth(method='lm', formula= y~x, aes(color=MSI, fill=MSI)) +
  geom_point(shape=21, aes(fill=MSI), color="black", size=4, alpha=0.8) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  scale_color_manual(values =  c("#339999", "#CC6600")) +
  scale_fill_manual(values =  c("#339999", "#CC6600")) +
  facet_grid(rows=vars(patient)) +
  ylab("Shannon diversity") +
  xlab("Proportion of Genome Altered") +
  theme_custom() +
  theme(legend.text = element_text(size=28), legend.title = element_blank(), legend.position = "top")
dev.off()

mod <- lm(shannon ~ patient + pga, data=metadata_tsm.G_sample[which(metadata_tsm.G_sample$MSI=="MSS"),])
summary(mod)

### association between alpha diversity and number of segments (per sample) (x)----
ggplot(metadata_tsm.G_patient, aes(x=(nSegments) , y=(shannon), color=MSI)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() 

ggplot(metadata_tsm.G_patient, aes(x=nSegments , y=simpson, color=MSI)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() 



### association between alpha diversity and ploidy (per sample/region)----
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(metadata_tsm.G_patient[!is.na(metadata_tsm.G_patient$ploidy),], aes(x=(ploidy) , y=(shannon), color=MSI)) +
  geom_point(shape=21, aes(fill=MSI), color="black", alpha=0.8, size=4) +
  geom_smooth(method='lm', formula= y~x, aes(color=MSI, fill=MSI)) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  scale_color_manual(values =  c("#339999", "#CC6600")) +
  scale_fill_manual(values =  c("#339999", "#CC6600")) +
  ylab("Average Shannon diversity per patient") +
  xlab("Ploidy") +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=28), legend.title = element_blank())
dev.off()

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(metadata_tsm.G_patient[!is.na(metadata_tsm.G_patient$ploidy),], aes(x=(ploidy) , y=(simpson), color=MSI)) +
  geom_point(shape=21, aes(fill=MSI), color="black", alpha=0.8, size=4) +
  geom_smooth(method='lm', formula= y~x, aes(color=MSI, fill=MSI)) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  scale_color_manual(values =  c("#339999", "#CC6600")) +
  scale_fill_manual(values =  c("#339999", "#CC6600")) +
  ylab("Average Simpson's diversity per patient") +
  xlab("Ploidy") +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=28), legend.title = element_blank())
dev.off()



ggplot(metadata_tsm.G_region, aes(x=(ploidy) , y=(shannon), color=MSI)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

### alpha and WGD (per sample)----
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(metadata_tsm.G_patient[which(metadata_tsm.G_patient$microbiomeData==TRUE & !is.na(metadata_tsm.G_patient$WGD)),], aes(x=WGD , y=(shannon),  fill=WGD)) +
  geom_violin(alpha=0.5) +
  geom_boxplot(alpha=0.6, outlier.size = 0, width=0.5) +
  geom_jitter(size = 3, shape=21, alpha=0.8, width = 0.19, height = 0) +
  xlab(NULL) +
  #ylim(0.4,1.01) +
  ylab("Average Shannon diversity per patient") +
  scale_fill_manual(values =  c("#F66666", "#009900")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("TRUE","FALSE")), size=9) +
  scale_x_discrete(labels = c("diploid","WGD")) +
  theme_custom() +
  theme(axis.text.x = element_text(), legend.position = "none")
dev.off()

t.test(metadata_tsm.G_sample$simpson[which(metadata_tsm.G_sample$microbiomeData==TRUE & metadata_tsm.G_sample$WGD==TRUE)],
       metadata_tsm.G_sample$simpson[which(metadata_tsm.G_sample$microbiomeData==TRUE & metadata_tsm.G_sample$WGD==FALSE)])

ggplot(metadata_tsm.G_sample[which(metadata_tsm.G_sample$microbiomeData==TRUE & !is.na(metadata_tsm.G_sample$WGD)),], aes(x=WGD , y=(simpson),  fill=WGD)) +
  geom_violin(alpha=0.3) +
  geom_boxplot(alpha=0.6) +
  geom_jitter(alpha=0.5, shape=21) +
  stat_compare_means(method = "t.test", comparisons = list(c("TRUE","FALSE"))) +
  scale_fill_manual(values =  c("#339999", "#CC6600")) +
  ylab("Simpon diversity") +
  xlab(NULL) +
  scale_x_discrete(labels=c("notWGD", "WGD")) +
  theme_custom() +
  theme(legend.position = "none")

t.test(metadata_tsm.G_sample$simpson[which(metadata_tsm.G_sample$microbiomeData==TRUE & metadata_tsm.G_sample$WGD==TRUE)],
       metadata_tsm.G_sample$simpson[which(metadata_tsm.G_sample$microbiomeData==TRUE & metadata_tsm.G_sample$WGD==FALSE)])

### association between number of mutations and average alpha (per patient)----

ggplot(metadata_tsm.G_patient[which(metadata_tsm.G_patient$microbiomeData==TRUE),], aes(x=log(totalSNVperpatient) , y=(shannon))) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5))

ggplot(metadata_tsm.G_patient[which(metadata_tsm.G_patient$microbiomeData==TRUE),], aes(x=log(totalSNVperpatient) , y=(simpson), color=MSI)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5))


### MSS and alpha diversity (per patient) ----

ggplot(metadata_tsm.G_sample[which(metadata_tsm.G_sample$microbiomeData==TRUE),], aes(x=MSI , y=shannon,  color=MSI)) +
  geom_violin() +
  geom_boxplot(outlier.shape = 0) +
  geom_jitter() +
  stat_compare_means(method = "t.test", comparisons = list(c("MSI","MSS"))) +
  theme_custom() +
  theme(legend.position = "none")

ggplot(metadata_tsm.G_sample[which(metadata_tsm.G_sample$microbiomeData==TRUE),], aes(x=MSI , y=simpson,  color=MSI)) +
  geom_violin() +
  geom_boxplot(outlier.shape = 0) +
  geom_jitter() +
  stat_compare_means(method = "t.test", comparisons = list(c("MSI","MSS"))) +
  theme_custom() +
  theme(legend.position = "none")


### ======BETA DISSIMILARITY=============================== ----
### bray PCoA by region (x) ----

# bray on species-level on NOT log-transformed data separates best
ordination <- ordinate(bacteria_tsm.G, method="PCoA", distance="bray")
plot_ordination(bacteria_tsm.G, ordination, color="region") + 
  stat_ellipse() +
  theme(aspect.ratio=1)

### bray PCoA per patient (x) ----
# bray on species-level NOT log-transformed data separates best
ordination <- ordinate(bacteria_tsm.G, method="PCoA", distance="bray")
plot_ordination(bacteria_tsm.G, ordination, color="patient") + 
  #stat_ellipse() +
  theme_custom() +
  theme(legend.position = "right",
        panel.grid.major = element_line(size=0.25),
        panel.grid.minor = element_line(size=0.25),
        legend.text = element_text(size=28),  
        legend.title = element_blank()
  )

theme(aspect.ratio=1)
jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
plot
dev.off()



### beta diversity (per region, by patient) ----

# plotting all the BC differences between sample per patient. one dot=one comparison. only car glands
jpeg('tempfig.jpeg', width = (35), height = (20), units = "cm", res = 300)
ggplot(metadata_tsm.G_region[which(!is.na(metadata_tsm.G_region$braycurtis)),], aes(x=reorder(patient, braycurtis, mean, na.rm=T) , y=braycurtis,  fill=patient)) +
  geom_boxplot(alpha=0.6, outlier.size = 0) +
  geom_jitter(shape=21, alpha=0.5, width = 0.2, size=3) +
  ylab("Pairwise Bray-Curtis") +
  xlab(NULL) +
  ylim(0,0.8) +
  theme_custom() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=90, vjust=0.5, size=20),
        panel.grid.major.y = element_line(size=0.5))
dev.off()

# by MSI
ggplot(metadata_tsm.G_region[which(!is.na(metadata_tsm.G_region$braycurtis)),], aes(x=MSI , y=braycurtis,  fill=MSI, shape = type)) +
  geom_boxplot(alpha=0.6, outlier.size = 0) +
  geom_jitter(shape=21, alpha=0.5, width = 0.2) +
  xlab(NULL) +
  ylab("Pairwise Bray-Curtis dissimilarity") +
  stat_compare_means(method = "t.test", comparisons = list(c("MSI","MSS"))) +
  scale_fill_viridis(option = "D", discrete = TRUE) +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5), legend.position = "none")

### ======BETA DISSIMILARITY AND CIN======================= ----
### association between beta dissimilarity and PGA (per region/patient)----
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(metadata_tsm.G_patient, aes(x=(pga) , y=(braycurtis), color=MSI)) +
  geom_smooth(method='lm', formula= y~x, ) +
  geom_point(shape=21, aes(fill=MSI), color="black", size=4, alpha=0.8) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  scale_color_manual(values =  c("#339999", "#CC6600")) +
  scale_fill_manual(values =  c("#339999", "#CC6600")) +
  ylab("Average Bray-Curtis diversity per patient") +
  xlab("Proportion of Genome Altered") +
  theme_custom() +
  theme(legend.text = element_text(size=28), legend.title = element_blank(), legend.position = "top")
dev.off()

### association between beta dissimilarity and number of segments in MSS (per region)----
ggplot(metadata_tsm.G_patient, aes(x=(nSegments) , y=(braycurtis), color=MSI)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5))

### positive change in beta dissimilarity with ploidy (per region)----
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(metadata_tsm.G_patient, aes(x=(ploidy) , y=(braycurtis), color=MSI)) +
  geom_smooth(method='lm', formula= y~x, aes(color=MSI, fill=MSI)) +
  geom_point(shape=21, aes(fill=MSI), color="black", size=4, alpha=0.8) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  scale_color_manual(values =  c("#339999", "#CC6600")) +
  scale_fill_manual(values =  c("#339999", "#CC6600")) +
  ylab("Average Bray-Curtis diversity per patient") +
  xlab("Ploidy") +
  theme_custom() +
  theme(legend.text = element_text(size=28), legend.title = element_blank(), legend.position = "top")
dev.off()

### beta dissimilarity with WGD (per region/patient)----
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(metadata_tsm.G_patient[which(!is.na(metadata_tsm.G_patient$WGD)),], aes(x=WGD , y=(braycurtis),  fill=WGD)) +
  geom_violin(alpha=0.5) +
  geom_boxplot(alpha=0.6, outlier.size = 0, width=0.5) +
  geom_jitter(size = 3, shape=21, alpha=0.8, width = 0.19, height = 0) +
  xlab(NULL) +
  #ylim(0.4,1.01) +
  ylab("Average Bray-Curtis diversity per patient") +
  scale_fill_manual(values =  c("#F66666", "#009900")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("TRUE","FALSE")), size=9) +
  scale_x_discrete(labels = c("diploid","WGD")) +
  theme_custom() +
  theme(axis.text.x = element_text(), legend.position = "none")
dev.off()

t.test(metadata_tsm.G_patient$braycurtis[which(metadata_tsm.G_patient$WGD==TRUE)],
       metadata_tsm.G_patient$braycurtis[which(metadata_tsm.G_patient$WGD==FALSE)])

### association between number of mutations and beta dissimilarity in MSS (per patient)----
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(metadata_tsm.G_patient[which(metadata_tsm.G_patient$microbiomeData==TRUE),], aes(x=(totalSNVperpatient) , y=(braycurtis), color=MSI)) +
  geom_smooth(method='lm', formula= y~x, aes(color=MSI, fill=MSI)) +
  geom_point(shape=21, aes(fill=MSI), color="black", size=4, alpha=0.8) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  scale_color_manual(values =  c("#339999", "#CC6600")) +
  scale_fill_manual(values =  c("#339999", "#CC6600")) +
  ylab("Average Bray-Curtis diversity per patient") +
  xlab("Total SNV count") +
  theme_custom() +
  theme(legend.text = element_text(size=28), legend.title = element_blank(), legend.position = "top")
dev.off()

### association between MSI and beta dissimilarity (per patient) ----

ggplot(metadata_tsm.G_patient[which(metadata_tsm.G_patient$microbiomeData==TRUE),], aes(x=MSI , y=(braycurtis),  color=MSI)) +
  geom_violin() +
  geom_boxplot(outlier.shape = 0) +
  geom_jitter() +
  stat_compare_means(method = "t.test", comparisons = list(c("MSI","MSS"))) +
  theme_custom() 

### no: pic and beta ----

epicc.raw <- readRDS("~/Documents/CNA/Github/Data/epicc.raw.rds")

epicc.info <- PullDataInfo(rawdata = epicc.raw) # hg19
epicc.clonality <- PullDataClonality(rawdata = epicc.raw, dataInfo = epicc.info)
epicc.diversity <- PullDataDiversity(rawdata = epicc.raw, dataInfo = epicc.info)
epicc.actualITH <- data.frame(patient = lapply(data.frame(patient=epicc.info$patientIDs), rep, epicc.info$sampPerPatient),
                              sample = epicc.info$sampleIDs,
                              actual = lapply(data.frame(actual=epicc.diversity$pic.frac$pic.frac), rep, epicc.info$sampPerPatient))

metadata_tsm.G_patient$picActual <- epicc.actualITH[match(metadata_tsm.G_patient$patient, epicc.actualITH$patient), 3]
metadata_tsm.G_patient$picPred <- epicc.diversity$pic.frac[match(metadata_tsm.G_patient$patient, epicc.diversity$pic.frac), 3]

ggplot(metadata_tsm.G_patient[which(metadata_tsm.G_patient$MSI=="MSI"),], aes(x=(picActual) , y=(braycurtis))) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom()

ggplot(metadata_tsm.G_patient, aes(x=(picActual) , y=(simpson))) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(axis.text.x = element_text( angle=90, vjust=0.5))

### ======IMMUNE ACTIVATION================================ ----
### immune composition (per sample, by patient) ----
quanTIseq_qcPass.long <- gather(quanTIseq_sample, key="key", value="value", -sample, -patient, -region, -subregion)

jpeg('tempfig.jpeg', width = (40), height = (40), units = "cm", res = 300)
ggplot(quanTIseq_qcPass.long, aes(x=patient, y=log(value),  color=patient)) +
  geom_violin() +
  facet_grid(rows = vars(key)) +
  geom_boxplot(outlier.shape = 0) +
  geom_jitter() +
  theme_custom() +
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5), legend.position = "none")
dev.off()

### association between alpha and immune (per region) ----

# remove sample that is 100% B cells
wd <- metadata_tsm.G_region[which(metadata_tsm.G_region$Bcells!=1),]

# shannon
x <- gather(wd[,which(colnames(wd) %in% c("region","shannon","Bcells","M1","M2","Mono",
                                          "Neutro","NK","CD4","CD8","Treg","Dend","MSI"))], key="key", value="value", -shannon, -region, -MSI)

jpeg('tempfig.jpeg', width = (72), height = (20), units = 'cm', res = 300)
ggplot(x, aes(x=(shannon) , y=(value), color=key)) +
  geom_point(size=3) +
  facet_grid(rows=vars(MSI)) +
  xlab("Ave. Shannon per region") +
  ylab("Immune cell proportion") +
  facet_grid(cols=vars(key), rows=vars(MSI)) +
  geom_smooth(method='lm', formula= y~x) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*",stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=6, ) +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        legend.position="none",
        strip.text = element_text(size=24))
dev.off()

# simpson
x <- gather(wd[,which(colnames(wd) %in% c("region","simpson","Bcells","M1","M2","Mono",
                                          "Neutro","NK","CD4","CD8","Treg","Dend","MSI"))], key="key", value="value", -simpson, -region, -MSI)

jpeg('tempfig.jpeg', width = (72), height = (20), units = 'cm', res = 300)
ggplot(x, aes(x=(simpson) , y=(value), color=key)) +
  geom_point(size=3) +
  facet_grid(rows=vars(MSI)) +
  xlab("Ave. Simpson's per region") +
  ylab("Immune cell proportion") +
  facet_grid(cols=vars(key), rows=vars(MSI)) +
  geom_smooth(method='lm', formula= y~x) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*",stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=5, ) +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        legend.position="none",
        strip.text = element_text(size=24))
dev.off()


### association between ploidy and immune (per region) ----

# remove sample that is 100% B cells
wd <- metadata_tsm.G_region[which(metadata_tsm.G_region$Bcells!=1),]
wd$totalImmune <- 1-wd$Other
#wd$Bcells <- wd$Bcells/wd$totalImmune
#wd$CD4 <- wd$CD4/wd$totalImmune
#wd$CD8 <- wd$CD8/wd$totalImmune
#wd$Dend <- wd$Dend/wd$totalImmune
#wd$M1 <- wd$M1/wd$totalImmune
#wd$M2 <- wd$M2/wd$totalImmune
#wd$Mono <- wd$Mono/wd$totalImmune
#wd$Neutro <- wd$Neutro/wd$totalImmune
#wd$NK <- wd$NK/wd$totalImmune
#wd$Treg <- wd$Treg/wd$totalImmune


x <- gather(wd[,which(colnames(wd) %in% c("region","ploidy","Bcells","M1","M2","Mono",
                                          "Neutro","NK","CD4","CD8","Treg","Dend","MSI","lymphocytes"))], key="key", value="value", -ploidy, -region, -MSI)
x$WGD <- ifelse(x$ploidy>2.5, TRUE, FALSE)
x <- na.omit(x)

ggplot(x, aes(x=(WGD) , y=(value), color=key)) +
  geom_jitter() +
  geom_boxplot() +
  xlab("Ave. ploidy per region") +
  facet_grid(cols=vars(key)) +
  #geom_smooth(method='lm', formula= y~x) +
  #ggpmisc::stat_poly_eq(formula = y ~ x, 
  #                      aes(label = paste(stat(adj.rr.label), "*\", \"*",stat(p.value.label), "*\"\"", sep = "")),
  #                      parse = TRUE, label.x = 'left', label.y = 0.8, size=2, ) +
  stat_compare_means(method = "t.test", size=5, comparisons = list(c("TRUE","FALSE"))) +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        legend.position="none",
        strip.text = element_text(size=24))

### association between beta dissimilarity and immune (per patient) (x) ----

# remove sample that is 100% B cells
wd <- metadata_tsm.G_patient[which(metadata_tsm.G_patient$Bcells!=1),]

x <- gather(wd[,which(colnames(wd) %in% c("patient","braycurtis","Bcells","M1","M2","Mono",
                                          "Neutro","NK","CD4","CD8","Treg","Dend","MSI"))], key="key", value="value", -braycurtis, -patient, -MSI)

jpeg('tempfig.jpeg', width = (72), height = (20), units = 'cm', res = 300)
ggplot(x, aes(x=(braycurtis) , y=(value), color=key)) +
  geom_point(size=3) +
  facet_grid(rows=vars(MSI)) +
  xlab("Bray curtis per patient") +
  ylab("Immune cell proportion") +
  facet_grid(cols=vars(key), rows=vars(MSI)) +
  geom_smooth(method='lm', formula= y~x) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*",stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=2, ) +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        legend.position="none",
        strip.text = element_text(size=24))
dev.off()
### ======SEPARATION OF REGIONS============================ ----
### PCoA ----

#create permanova p
permanovaP <- list()
patients <- unique(sample_data(input_tsm_decon_rem)$patient)
for (i in 1:length(patients)) {
  wd <- subset_samples(input_tsm_decon_rem_car, patient==patients[i])
  regionFreq <- data.frame(table(sample_data(wd)$region))
  keep <- regionFreq$Var1[which(regionFreq$Freq>1)]
  
  if (length(keep)<2) {
    permanovaP[[i]] <- NA
    next
  }
  
  wd <- subset_samples(wd, region %in% keep)
  wd <- microbiome::transform(wd, "compositional")
  
  tax <- data.frame(tax_table(wd))
  tax$name <- paste(substr(tax$Genus,1,1), tax$Species, sep = ". ")
  
  otu <- abundances(wd)
  meta <- meta(wd)
  #meta$region[which(meta$region %in% c("C562_B","C562_D"))] <- "C562_BD"
  meta$region <- factor(meta$region)
  
  set.seed(123)
  permanova <- adonis(t(otu) ~ region,
                      data = meta, permutations=99, method = "bray")
  
  permanovaP[[i]] <- print(as.data.frame(permanova$aov.tab)["region", "Pr(>F)"])
  
  #permutest(betadisper(vegdist(t(otu), method = "bray"), meta$region), pairwise = TRUE) # Null hypothesis = no difference in dispersion between groups.
}



#plot
metadata_tsm.G_sample$r <- str_after_first(metadata_tsm.G_sample$region,"_")
sample_data(input_tsm_decon_rem)$r <- str_after_first(sample_data(input_tsm_decon_rem)$region,"_")

mds_plot.perpatient <- list()
l <- 1
for(i in 1:length(unique(sample_data(input_tsm_decon_rem)$patient))) {
  
  wd <- subset_samples(input_tsm_decon_rem, patient==unique(sample_data(input_tsm_decon_rem)$patient)[i])
  
  if (ncol(otu_table(wd))==1) {
    next
  }
  
  ordination <- ordinate(wd, method="PCoA", distance="bray")
  mds_plot.perpatient[[l]] <- plot_ordination(wd, ordination) + 
    geom_point(size = 8, aes(color = as.factor(r), shape = type)) +
    scale_color_manual(name = 'region', 
                       values = c("A" = '#99CC00', "B" = '#9966CC', "C" = '#336666', "D" = '#FF6666', 
                                  "E" = "#00CCCC", "F" = "#FF9933", "G" = "#FFFF00", "H" = "#6600FF", "W" = "#993399"),
                       drop = FALSE) +
    scale_shape_manual(values = c("B" = 15, "G" = 16, "L" = 17),
                     drop = FALSE) +
    theme_custom() +
    labs(title = unique(sample_data(input_tsm_decon_rem)$patient)[i], 
         subtitle = paste("PERMANOVA p",permanovaP[[i]], sep = "=")) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size=24),
          legend.position = "right",
          panel.grid.major = element_line(size=0.25),
          panel.grid.minor = element_line(size=0.25),
          legend.text = element_text(size=28),  
          legend.title = element_blank() 
    ) 
  l <- l + 1
}


mds_plot.perpatient[[1]]
#create lgend
df <- data.frame(A=c(1,2,1),
                 B=c(2, 3,1),
                 C=c(1,1,1),
                 D=c(1,1,1),
                 E=c(1,1,1),
                 F=c(1,1,1),
                 G=c(1,1,1),
                 H=c(1,1,1),
                 W=c(1,1,1))
rownames(df) <- c("B","G","L")
df$this <- rownames(df)
df <- gather(df, key="key", value="value", -this)

jpeg('tempfig.jpeg', width = (40), height = (30), units = 'cm', res = 300)
ggplot(df, aes(x=key, y=this, color=key, shape=this)) +
  geom_point(size=5) +
  scale_color_manual(name = 'Region', 
                     values = c("A" = '#99CC00', "B" = '#9966CC', "C" = '#336666', "D" = '#FF6666', 
                                "E" = "#00CCCC", "F" = "#FF9933", "G" = "#FFFF00", "H" = "#6600FF", "W" = "#993399")) +
  scale_shape_manual(name = 'Sample type', values = c("B" = 15, "G" = 16, "L" = 17)) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_text(size=20))
dev.off()

#jpeg('tempfig.jpeg', width = (3*37.795*25), height = (3*37.795*19))
jpeg('tempfig.jpeg', width = (70), height = (70), units = 'cm', res = 300)
ggarrange(plotlist=mds_plot.perpatient, ncol = 5, nrow=6, common.legend=TRUE, legend='bottom')
dev.off()

### BC between carcinoma samples of same v different region ----

# grab all pairwise bray curtis measure within a patient
braycurtis_dist <- phyloseq::distance(input_tsm_decon_rem, method="bray")
braycurtis_dist <- as.matrix(braycurtis_dist)

sub_dist <- list()
patient_all <- sample_data(input_tsm_decon_rem)$patient

for (i in 1:length(unique(patient_all)) ) { 
  patient <- unique(patient_all)[i]
  index <- which(patient_all == patient)
  sample_group <- sample_names(input_tsm_decon_rem)[index]
  sub_dist[[patient]] <- braycurtis_dist[sample_group, sample_group]
  sub_dist[[patient]][!lower.tri(sub_dist[[patient]])] <- NA
}

braygroups <- melt(sub_dist)
df.bray_patient <- braygroups[complete.cases(braygroups), ]
df.bray_patient$patient <- factor(df.bray_patient$L1, levels=names(sub_dist))

# is pairwise c-c, n-n, c-n
df.bray_patient$Id1 <- substr(df.bray_patient$Var1, 7, nchar(as.character(df.bray_patient$Var1))-23)
df.bray_patient$region1 <- substr(df.bray_patient$Id1,1,6)
df.bray_patient$Id2 <- substr(df.bray_patient$Var2, 7, nchar(as.character(df.bray_patient$Var2))-23)
df.bray_patient$region2 <- substr(df.bray_patient$Id2,1,6)

# plot per comparison
df.bray_patient$comp <- ifelse(df.bray_patient$region1==df.bray_patient$region2, "sameRegion", "differentRegion")

mean(df.bray_patient$value[which(df.bray_patient$comp=="sameRegion")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="differentRegion")])

# add sample type
x <- sample_data(input_tsm_decon_rem)
df.bray_patient$type1 <- x[match(df.bray_patient$Id1, x$Id),7]
df.bray_patient$type2 <- x[match(df.bray_patient$Id2, x$Id),7]
df.bray_patient$cancerOnly <- ifelse(df.bray_patient$type1=="carcinoma" & df.bray_patient$type2=="carcinoma", TRUE, FALSE)


jpeg('tempfig.jpeg', width = (22), height = (20), units = 'cm', res = 300)
ggplot(df.bray_patient[which(df.bray_patient$cancerOnly==T),], aes(x=reorder(comp, value, median) , y=value,  fill=comp)) +
  geom_violin(alpha=0.3, width=0.8, position = position_dodge(0.8)) +
  geom_boxplot(alpha=0.6, width=0.3, position = position_dodge(0.8)) +
  geom_jitter(alpha=0.2, shape=21, size=3,
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2) ) +
  stat_compare_means(method = "wilcox.test", size=9, comparisons = list(c("sameRegion","differentRegion")), label = "p.format") +
  scale_fill_manual(values =  c("#CC6699", "#CC6600")) +
  scale_x_discrete(labels=c("Within-region \n(carcinoma only)","Between-regions \n(carcinoma only)")) +
  xlab(NULL) +
  ylim(0,1.1) +
  ylab("Pairwise Bray-Curtis dssimilarity") +
  theme_custom() +
  theme(legend.position = "none", legend.text = element_text(size=24),
        legend.title = element_blank())
dev.off()

mean(df.bray_patient$value[which(df.bray_patient$cancerOnly==T & df.bray_patient$comp=="differentRegion")])
mean(df.bray_patient$value[which(df.bray_patient$cancerOnly==T & df.bray_patient$comp=="sameRegion")])

### BC between carcinoma samples vs bc between carcinoma and normal ----

# grab all pairwise bray curtis measure within a patient
braycurtis_dist <- phyloseq::distance(input_tsm_decon_rem, method="bray")
braycurtis_dist <- as.matrix(braycurtis_dist)

sub_dist <- list()
patient_all <- sample_data(input_tsm_decon_rem)$patient

for (i in 1:length(unique(patient_all)) ) { 
  patient <- unique(patient_all)[i]
  index <- which(patient_all == patient)
  sample_group <- sample_names(input_tsm_decon_rem)[index]
  sub_dist[[patient]] <- braycurtis_dist[sample_group, sample_group]
  sub_dist[[patient]][!lower.tri(sub_dist[[patient]])] <- NA
}

braygroups <- melt(sub_dist)
df.bray_patient <- braygroups[complete.cases(braygroups), ]
df.bray_patient$patient <- factor(df.bray_patient$L1, levels=names(sub_dist))

# is pairwise c-c, n-n, c-n
df.bray_patient$Id1 <- substr(df.bray_patient$Var1, 7, nchar(as.character(df.bray_patient$Var1))-23)
df.bray_patient$group1 <- ifelse( (substr(df.bray_patient$Id1,6,6) %in% c("F","G","H")), "adenoma", 
                                  ifelse( (substr(df.bray_patient$Id1,6,6) %in% c("E")), "normal",
                                          ifelse( (substr(df.bray_patient$Id1,6,6) %in% c("A","B")), "carcinoma", 
                                                  ifelse( (substr(df.bray_patient$Id1,6,6) %in% c("C","D") & patient == "C516"), "adenoma",
                                                          ifelse( (substr(df.bray_patient$Id1,6,6) %in% c("C","D") & patient != "C516"), "carcinoma", 
                                                                  ifelse( (substr(df.bray_patient$Id1,6,6) %in% c("Z")), "blood", 
                                                                          ifelse( (substr(df.bray_patient$Id1,6,6) %in% c("W")), "normalInCancer", NA)))))))

df.bray_patient$Id2 <- substr(df.bray_patient$Var2, 7, nchar(as.character(df.bray_patient$Var2))-23)
df.bray_patient$group2 <- ifelse( (substr(df.bray_patient$Id2,6,6) %in% c("F","G","H")), "adenoma", 
                                  ifelse( (substr(df.bray_patient$Id2,6,6) %in% c("E")), "normal",
                                          ifelse( (substr(df.bray_patient$Id2,6,6) %in% c("A","B")), "carcinoma", 
                                                  ifelse( (substr(df.bray_patient$Id2,6,6) %in% c("C","D") & patient == "C516"), "adenoma",
                                                          ifelse( (substr(df.bray_patient$Id2,6,6) %in% c("C","D") & patient != "C516"), "carcinoma", 
                                                                  ifelse( (substr(df.bray_patient$Id2,6,6) %in% c("Z")), "blood", 
                                                                          ifelse( (substr(df.bray_patient$Id2,6,6) %in% c("W")), "normalInCancer", NA)))))))


# plot per comparison
df.bray_patient$comp <- paste(df.bray_patient$group1, df.bray_patient$group2, sep = " v ")
df.bray_patient <- df.bray_patient[which(df.bray_patient$comp %!in% c("normalInCancer v normal", "normalInCancer v carcinoma")),]

mean(df.bray_patient$value[which(df.bray_patient$comp=="normal v normal")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="adenoma v normal")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="normal v carcinoma")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="adenoma v normal")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="adenoma v adenoma")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="adenoma v carcinoma")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="normal v carcinoma")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="carcinoma v carcinoma")])

comparisons <- list(
  c("normal v normal","adenoma v adenoma"),
  c("normal v normal","adenoma v normal"),
  c("normal v normal","carcinoma v carcinoma"),
  c("normal v normal","normal v carcinoma"),
  c("normal v normal","adenoma v carcinoma"),
  
  c("adenoma v adenoma","adenoma v normal"),
  c("adenoma v adenoma","carcinoma v carcinoma"),
  c("adenoma v adenoma","normal v carcinoma"),
  c("adenoma v adenoma","adenoma v carcinoma"),
  
  c("adenoma v normal","carcinoma v carcinoma"),
  c("adenoma v normal","normal v carcinoma"),
  c("adenoma v normal","adenoma v carcinoma"),
  
  c("carcinoma v carcinoma","normal v carcinoma"),
  c("carcinoma v carcinoma","adenoma v carcinoma"),
  
  c("normal v carcinoma","adenoma v carcinoma")
)

jpeg('tempfig.jpeg', width = (40), height = (40), units = 'cm', res = 300)
ggplot(df.bray_patient, aes(x=reorder(comp, value, mean) , y=value,  fill=comp)) +
  geom_violin(alpha=0.3, width=1, position = position_dodge(1)) +
  geom_boxplot(alpha=0.6, width=0.3, position = position_dodge(1)) +
  geom_jitter(alpha=0.5, shape=21, size=3,
              position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1) ) +
  stat_compare_means(method = "t.test", size=9, comparisons = comparisons, label = "p.format", step.increase = 0.18) +
  scale_fill_viridis(option = "D", discrete = TRUE) +
  #scale_x_discrete(labels = c("Gland","Bulk","Leftover")) +
  xlab(NULL) +
  ylim(0,3.5) +
  ylab("Pairwise Bray-Curtis dssimilarity") +
  theme_custom() +
  theme(legend.position = "none", legend.text = element_text(size=28),
        legend.title = element_blank(),
        axis.text.x = element_text(angle=45, vjust=1, size=28, hjust=1))
dev.off()

### BC between carcinoma samples vs bc between G and B ----

# grab all pairwise bray curtis measure within a patient
braycurtis_dist <- phyloseq::distance(input_tsm_decon_rem, method="bray")
braycurtis_dist <- as.matrix(braycurtis_dist)

sub_dist <- list()
patient_all <- sample_data(input_tsm_decon_rem)$patient

for (i in 1:length(unique(patient_all)) ) { 
  patient <- unique(patient_all)[i]
  index <- which(patient_all == patient)
  sample_group <- sample_names(input_tsm_decon_rem)[index]
  sub_dist[[patient]] <- braycurtis_dist[sample_group, sample_group]
  sub_dist[[patient]][!lower.tri(sub_dist[[patient]])] <- NA
}

braygroups <- melt(sub_dist)
df.bray_patient <- braygroups[complete.cases(braygroups), ]
df.bray_patient$patient <- factor(df.bray_patient$L1, levels=names(sub_dist))

# is pairwise c-c, n-n, c-n
df.bray_patient$Id1 <- substr(df.bray_patient$Var1, 7, nchar(as.character(df.bray_patient$Var1))-23)
df.bray_patient$type1 <- substr(str_after_nth(df.bray_patient$Id1,"_",2),1,1)
df.bray_patient$group1 <- substr(str_after_nth(df.bray_patient$Id1,"_",1),1,1)
df.bray_patient$Id2 <- substr(df.bray_patient$Var2, 7, nchar(as.character(df.bray_patient$Var2))-23)
df.bray_patient$type2 <- substr(str_after_nth(df.bray_patient$Id2,"_",2),1,1)
df.bray_patient$group2 <- substr(str_after_nth(df.bray_patient$Id2,"_",1),1,1)

# plot per comparison
df.bray_patient$comp <- ifelse(df.bray_patient$group1==df.bray_patient$group2 & df.bray_patient$type1==df.bray_patient$type2,
                               "sameRsameT", ifelse(df.bray_patient$group1==df.bray_patient$group2 & df.bray_patient$type1!=df.bray_patient$type2,
                                                    "sameRdifferentT",ifelse(df.bray_patient$group1!=df.bray_patient$group2 & df.bray_patient$type1==df.bray_patient$type2,
                                                                             "diffRsameT","diffRdiffT")))


mean(df.bray_patient$value[which(df.bray_patient$comp=="normal v normal")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="adenoma v normal")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="normal v carcinoma")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="adenoma v normal")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="adenoma v adenoma")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="adenoma v carcinoma")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="normal v carcinoma")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="carcinoma v carcinoma")])

comparisons <- list(
  c("sameRsameT","sameRdifferentT"),
  c("sameRsameT","diffRsameT"),
  c("sameRsameT","diffRdiffT"),
  c("sameRdifferentT","diffRsameT"),c("sameRdifferentT","diffRdiffT"),
  c("diffRdiffT","diffRsameT")
)

jpeg('tempfig.jpeg', width = (40), height = (30), units = 'cm', res = 300)
ggplot(df.bray_patient, aes(x=reorder(comp, value, median) , y=value,  fill=comp)) +
  geom_violin(alpha=0.3, width=1, position = position_dodge(1)) +
  geom_boxplot(alpha=0.6, width=0.3, position = position_dodge(1)) +
  geom_jitter(alpha=0.5, shape=21, size=3,
              position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1) ) +
  stat_compare_means(method = "wilcox.test", size=9, comparisons = comparisons, label = "p.format", step.increase = 0.18) +
  scale_fill_viridis(option = "D", discrete = TRUE) +
  scale_x_discrete(labels = c("Same region,\nsame type","Same region,\ndifferent type",
  "Different region,\nsame type","Different region,\ndifferent type")) +
  xlab(NULL) +
  ylim(0,2) +
  ylab("Pairwise Bray-Curtis dssimilarity") +
  theme_custom() +
  theme(legend.position = "none", legend.text = element_text(size=28),
        legend.title = element_blank(),
        axis.text.x = element_text(angle=45, vjust=1, size=28, hjust=1))
dev.off()

### BC between WGDvWGD compared to diploid-WGD ----

# grab all pairwise bray curtis measure within a patient
braycurtis_dist <- phyloseq::distance(input_tsm_decon_rem_car, method="bray")
braycurtis_dist <- as.matrix(braycurtis_dist)

sub_dist <- list()
patient_all <- sample_data(input_tsm_decon_rem_car)$patient

for (i in 1:length(unique(patient_all)) ) { 
  patient <- unique(patient_all)[i]
  index <- which(patient_all == patient)
  sample_group <- sample_names(input_tsm_decon_rem_car)[index]
  sub_dist[[patient]] <- braycurtis_dist[sample_group, sample_group]
  sub_dist[[patient]][!lower.tri(sub_dist[[patient]])] <- NA
}

braygroups <- melt(sub_dist)
df.bray_patient <- braygroups[complete.cases(braygroups), ]
df.bray_patient$patient <- factor(df.bray_patient$L1, levels=names(sub_dist))

# is pairwise c-c, n-n, c-n
df.bray_patient$Id1 <- substr(df.bray_patient$Var1, 7, nchar(as.character(df.bray_patient$Var1))-23)
df.bray_patient$group1 <- ifelse(df.bray_patient$Id1 %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$WGD==TRUE)], TRUE, FALSE )

df.bray_patient$Id2 <- substr(df.bray_patient$Var2, 7, nchar(as.character(df.bray_patient$Var2))-23)
df.bray_patient$group2 <- ifelse(df.bray_patient$Id2 %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$WGD==TRUE)], TRUE, FALSE )

# plot per comparison
df.bray_patient$comp <- paste(df.bray_patient$group1, df.bray_patient$group2, sep = " v ")
df.bray_patient$comp <- ifelse(df.bray_patient$comp=="FALSE v FALSE", "diploid v diploid", 
                               ifelse(df.bray_patient$comp=="TRUE v TRUE", "WGD v WGD", "diploid v WGD"))

mean(df.bray_patient$value[which(df.bray_patient$comp=="diploid v diploid")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="diploid v WGD")])
mean(df.bray_patient$value[which(df.bray_patient$comp=="WGD v WGD")])

comparisons <- list(
  c("diploid v diploid","diploid v WGD"),
  c("diploid v diploid","WGD v WGD"),
  c("diploid v WGD","WGD v WGD"))

jpeg('tempfig.jpeg', width = (30), height = (20), units = 'cm', res = 300)
ggplot(df.bray_patient, aes(x=reorder(comp, value, mean) , y=value,  fill=comp)) +
  geom_violin(alpha=0.3, width=1, position = position_dodge(1)) +
  geom_boxplot(alpha=0.6, width=0.5, position = position_dodge(1)) +
  geom_jitter(alpha=0.5, shape=21, size=3,
              position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2) ) +
  stat_compare_means(method = "t.test", size=9, comparisons = comparisons, label = "p.format", step.increase = 0.18) +
  scale_fill_viridis(option = "D", discrete = TRUE) +
  #scale_x_discrete(labels = c("Gland","Bulk","Leftover")) +
  xlab(NULL) +
  ylim(0,1.5) +
  ylab("Pairwise Bray-Curtis dssimilarity") +
  theme_custom() +
  theme(legend.position = "none", legend.text = element_text(size=28),
        legend.title = element_blank(),
        axis.text.x = element_text(angle=45, vjust=1, size=28, hjust=1))
dev.off()


### DA on diploid v WGD ----
dds_input <- input_tsm_decon_rem_carG
  
sample_data(dds_input)$ploidy <- unlist(ploidy.persample[match(sample_data(dds_input)$sample, gsub("\\.", "_", ploidy.persample$sample)), 3])
sample_data(dds_input)$WGD <- ifelse(sample_data(dds_input)$ploidy>2.5, TRUE, FALSE)
sample_data(dds_input)$MSI <- unlist(metadata_tsm.G_patient[match(sample_data(dds_input)$patient, metadata_tsm.G_patient$patient), 5])

dds_input <- subset_samples(dds_input, !is.na(WGD))
dds_input <- subset_samples(dds_input, MSI=="MSS")

dds <- phyloseq_to_deseq2(dds_input, ~ patient + WGD)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))} # calculate geometric means prior to estimate size factors: necessary because lots of zeros
geoMeans <- apply(counts(dds), 1, gm_mean)
deseq2 <- estimateSizeFactors(dds, geoMeans = geoMeans)
deseq2 <- DESeq(deseq2, fitType="local")
plotDispEsts(deseq2)

resultsNames(deseq2)
#res = results(deseq2, contrast = c("MSI","MSS","MSI")) #MSI is base
res = results(deseq2)
res = res[order(res$padj, na.last=NA), ]
sigtab2 = res[(res$padj < 0.01), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(dds_input)[rownames(sigtab2), ], "matrix"))

# plot
sigtabgen = subset(sigtab2, !is.na(Species))

# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$name <- paste(substring(sigtabgen$Genus,1,1),sigtabgen$Species,sep='. ')


jpeg('tempfig.jpeg', width = (30), height = (20), units = "cm", res = 300)
ggplot(sigtabgen, aes(y=reorder(name,log2FoldChange) , x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  ylab(NULL) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  theme_custom() +
  theme(legend.text = element_text(size=24), legend.title = element_text(size=24),
        axis.text.y = element_text(size=24))
dev.off()


### beta diversity (per region, by patient) and by MSI/WGD ----

# grab all pairwise bray curtis measure within a patient
braycurtis_dist <- phyloseq::distance(input_tsm_decon_rem_car, method="bray")
braycurtis_dist <- as.matrix(braycurtis_dist)

sub_dist <- list()
patient_all <- sample_data(input_tsm_decon_rem_car)$patient

for (i in 1:length(unique(patient_all)) ) { 
  patient <- unique(patient_all)[i]
  index <- which(patient_all == patient)
  sample_group <- sample_names(input_tsm_decon_rem_car)[index]
  sub_dist[[patient]] <- braycurtis_dist[sample_group, sample_group]
  sub_dist[[patient]][!lower.tri(sub_dist[[patient]])] <- NA
}

braygroups <- melt(sub_dist)
df.bray_patient <- braygroups[complete.cases(braygroups), ]
df.bray_patient$patient <- factor(df.bray_patient$L1, levels=names(sub_dist))

df.bray_patient$WGD <- metadata_tsm.G_patient[match(df.bray_patient$patient, metadata_tsm.G_patient$patient), 4]
df.bray_patient$MSI <- metadata_tsm.G_patient[match(df.bray_patient$patient, metadata_tsm.G_patient$patient), 5]

x <- aggregate(df.bray_patient$value, list(df.bray_patient$patient), FUN=median) 
x$MSI <- metadata_tsm.G_patient[match(x$Group.1, metadata_tsm.G_patient$patient), 5]

# plotting all the BC differences between sample per patient. one dot=one comparison. only car glands
ggplot(df.bray_patient, aes(x=reorder(patient, value, mean, na.rm=T) , y=value,  fill=patient)) +
  geom_boxplot(alpha=0.6, outlier.size = 0) +
  geom_jitter(shape=21, alpha=0.5, width = 0.2) +
  xlab(NULL) +
  ylab("Pairwise Bray-Curtis dissimilarity") +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5), legend.position = "none")

# by MSI
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(x=MSI , y=x,  fill=MSI)) +
  geom_boxplot(alpha=0.6, outlier.size = 0) +
  geom_jitter(shape=21, alpha=0.5, width = 0.2, size=3) +
  xlab(NULL) +
  ylim(0, 1.1) +
  ylab("Pairwise Bray-Curtis dissimilarity") +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("MSI","MSS")), size=9) +
  scale_fill_manual(values =  c("#339999", "#CC6600")) +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5), legend.position = "none")
dev.off()

t.test(df.bray_patient$value[which(df.bray_patient$MSI=="MSI")],
       df.bray_patient$value[which(df.bray_patient$MSI=="MSS")])

# by WGD
ggplot(df.bray_patient, aes(x=WGD , y=value,  fill=WGD)) +
  geom_boxplot(alpha=0.6, outlier.size = 0) +
  geom_jitter(shape=21, alpha=0.5, width = 0.2) +
  xlab(NULL) +
  ylab("Pairwise Bray-Curtis dissimilarity") +
  stat_compare_means(method = "t.test", comparisons = list(c("TRUE","FALSE"))) +
  scale_fill_viridis(option = "D", discrete = TRUE) +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5), legend.position = "none")

t.test(df.bray_patient$value[which(df.bray_patient$WGD=="TRUE")],
       df.bray_patient$value[which(df.bray_patient$WGD=="FALSE")])


### ======DA BY CIN======================================== ----
### DE metadata ----

# split samples by pga recentre (gene dysbosis hence recentre) 
jpeg('tempfig.jpeg', width = (15), height = (20), units = 'cm', res = 300)
ggplot(metadata_tsm.G_sample, aes(x=pga)) + 
  geom_histogram(data=subset(metadata_tsm.G_sample, pga4group == "1"), color="black", fill = "#1B9E77", alpha = 0.5, binwidth = 0.01) +
  geom_histogram(data=subset(metadata_tsm.G_sample, pga4group == '2'), color="black", fill = "#E6AB02", alpha = 0.5, binwidth = 0.01) +
  geom_histogram(data=subset(metadata_tsm.G_sample, pga4group == "3"), color="black", fill = "#9966cc", alpha = 0.5, binwidth = 0.01) +
  geom_histogram(data=subset(metadata_tsm.G_sample, pga4group == '4'), color="black", fill = "#006699", alpha = 0.5, binwidth = 0.01) +
  theme_custom()
dev.off()



# make pga groups (metadata_tsm.G_sample on only carcinoma glands)
metadata_tsm.G_sample$pga2group <- ntile(metadata_tsm.G_sample$pga, n=2)
metadata_tsm.G_sample$pga3group <- ntile(metadata_tsm.G_sample$pga, n=3)
metadata_tsm.G_sample$pga4group <- ntile(metadata_tsm.G_sample$pga, n=4)

metadata_tsm.G_sample$pga_recentre2group <- ntile(metadata_tsm.G_sample$pga_recentre, n=2)
metadata_tsm.G_sample$pga_recentre3group <- ntile(metadata_tsm.G_sample$pga_recentre, n=3)
metadata_tsm.G_sample$pga_recentre4group <- ntile(metadata_tsm.G_sample$pga_recentre, n=4)

metadata_tsm.G_sample$pga3group_2 <- ifelse(metadata_tsm.G_sample$pga<0.2, 1, ifelse(metadata_tsm.G_sample$pga>0.8, 3, 2))
metadata_tsm.G_sample$pga2group_2 <- ifelse(metadata_tsm.G_sample$pga<0.5, 1, 2)

#add clinical data
clinical_data <- read_excel("~/Documents/Microbial/data/clinical_data.xlsx", skip = 0)
clinical_data$ID <- paste("C",clinical_data$ID,sep='')
clinical_data$location_recode <- ifelse(clinical_data$Location=="Ascending colon", "right", 
                                        ifelse(clinical_data$Location=="Rectosigmoid", "left",
                                               ifelse(clinical_data$Location=="Transverse colon", "Transverse",
                                                      ifelse(clinical_data$Location=="Distal Sigmoid", "left",
                                                             ifelse(clinical_data$Location=="Hepatic flexure", "right",
                                                                    ifelse(clinical_data$Location=="Distal sigmoid/upper rectum", "left",
                                                                           ifelse(clinical_data$Location=="Sigmoid", "left",
                                                                                  ifelse(clinical_data$Location=="Sigmoid, mesenteric aspect", "left",
                                                                                         ifelse(clinical_data$Location=="Rectum", "left",
                                                                                                ifelse(clinical_data$Location=="Caecum", "right",
                                                                                                       ifelse(clinical_data$Location=="Lower Rectum", "left",
                                                                                                              ifelse(clinical_data$Location=="Ascending", "right",
                                                                                                                     ifelse(clinical_data$Location=="Proximal rectum", "left",
                                                                                                                            ifelse(clinical_data$Location=="Distal sigmoid", "left",
                                                                                                                                   ifelse(clinical_data$Location=="sigmoid", "left",
                                                                                                                                          ifelse(clinical_data$Location=="Hepatic Flexure", "right",NA))))))))))))))))

metadata_tsm.G_sample$side <- unlist(clinical_data[match(metadata_tsm.G_sample$patient, clinical_data$ID),36])

# add to metadata
wd <- subset_samples(input_raw_decon_rem_car, Id %in% metadata_tsm.G_sample$Id[!is.na(metadata_tsm.G_sample$pga)])

wd <- wd %>%
  ps_mutate(
    pga4group = as.factor(ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga4group==1)],1,
                                 ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga4group==2)], 2,
                                        ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga4group==3)], 3,
                                               ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga4group==4)], 4, NA) ) ) )),
    pga3group = as.factor(ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga3group==1)],1,
                                 ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga3group==2)], 2,
                                        ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga3group==3)], 3, NA) ) )), 
    pga3group_2 = as.factor(ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga3group_2==1)],1,
                                   ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga3group_2==2)], 2,
                                          ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga3group==3)], 3, NA) ) )), 
    pga2group_2 = as.factor(ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga2group_2==1)],1,2 )), 
    pga2group = as.factor(ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga2group==1)],1,
                                 ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga2group==2)], 2, NA) ) ),
    pga_recentre4group = as.factor(ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga_recentre4group==1)],1,
                                          ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga_recentre4group==2)], 2,
                                                 ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga_recentre4group==3)], 3,
                                                        ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga_recentre4group==4)], 4, NA) ) ) )),
    pga_recentre3group = as.factor(ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga_recentre3group==1)],1,
                                          ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga_recentre3group==2)], 2,
                                                 ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga_recentre3group==3)], 3, NA) ) )), 
    pga_recentre2group = as.factor(ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga_recentre2group==1)],1,
                                          ifelse(Id %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$pga_recentre2group==2)], 2, NA) ) ),
    ploidy = metadata_tsm.G_sample[match(Id, metadata_tsm.G_sample$Id), 11],
    WGD = metadata_tsm.G_sample[match(Id, metadata_tsm.G_sample$Id), 12],
    MSI = metadata_tsm.G_sample[match(Id, metadata_tsm.G_sample$Id), 13],
    side = metadata_tsm.G_sample[match(Id, metadata_tsm.G_sample$Id), 40]
  )

### DE on high v low on pga per sample ----

# remove MSI
#input <- subset_samples(input, MSI=="MSS")

# run DE
dds <- phyloseq_to_deseq2(wd, ~ patient + pga2group_2)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))} # calculate geometric means prior to estimate size factors: necessary because lots of zeros
geoMeans <- apply(counts(dds), 1, gm_mean)
deseq2 <- estimateSizeFactors(dds, geoMeans = geoMeans)
deseq2 <- DESeq(deseq2, fitType="local")
#plotDispEsts(deseq2)

res = results(deseq2, contrast = c("pga2group","2","1"))
#res = results(deseq2, contrast = c("MSI","MSI","MSS"))
#res = results(deseq2, contrast = c("side","right","left"))
#res = results(deseq2)
res = res[order(res$padj, na.last=NA), ]
sigtab2 = res[(res$padj < 0.05), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(input_raw_decon_rem_car)[rownames(sigtab2), ], "matrix"))

# plot
sigtabgen = subset(sigtab2, !is.na(Species))

# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$name <- paste(substring(sigtabgen$Genus,1,1),sigtabgen$Species,sep='. ')

jpeg('tempfig.jpeg', width = (25), height = (20), units = 'cm', res = 300)
ggplot(sigtabgen, aes(y=reorder(name,log2FoldChange) , x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  ylab(NULL) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  theme_custom() +
  theme(legend.text = element_text(size=24), legend.title = element_text(size=24),)
dev.off()

#jpeg('tempfig.jpeg', width = (3*37.795*10), height = (3*37.795*10))
#plot
#dev.off()

### ======PHYLO DISTANCE AND BC============================ ----
### phylo dist corr BC ----

# add tree
#random_tree = rtree(ntaxa(input_tsm_decon), rooted=TRUE, tip.label=taxa_names(input_tsm_decon))
#plot(random_tree)
#input_tsm_decon = merge_phyloseq(input_tsm_decon, random_tree)

# choose distance
#braycurtis_dist <- phyloseq::distance(input_tsm_decon, method="unifrac", weighted=FALSE)
braycurtis_dist <- phyloseq::distance(input_tsm_decon_rem, method="bray")
braycurtis_dist <- as.matrix(braycurtis_dist)

# calculate pairwise distance
sub_dist <- list()
patient_all <- sample_data(input_tsm_decon_rem)$patient

for (i in 1:length(unique(patient_all)) ) { 
  patient <- unique(patient_all)[i]
  index <- which(patient_all == patient)
  sample_group <- sample_names(input_tsm_decon_rem)[index]
  sub_dist[[patient]] <- braycurtis_dist[sample_group, sample_group]
  sub_dist[[patient]][!lower.tri(sub_dist[[patient]])] <- NA
}

braygroups <- melt(sub_dist)
df.bray_patient <- braygroups[complete.cases(braygroups), ]
df.bray_patient$patient <- factor(df.bray_patient$L1, levels=names(sub_dist))

df.bray_patient$Var1 <- str_before_first(as.character(df.bray_patient$Var1), ".k")
df.bray_patient$Var2 <- str_before_first(as.character(df.bray_patient$Var2), ".k")
colnames(df.bray_patient)[3] <- "BC"

# make comparison name variable
df.bray_patient$comp <- paste(pmin(df.bray_patient$Var1,df.bray_patient$Var2), pmax(df.bray_patient$Var1,df.bray_patient$Var2), sep = ":")


# tree distances per patient
library(adephylo)
final_tree_set <- readRDS("~/Documents/Microbial/data/final_tree_set.rds")
phylodist <- list()
rootdist <- list()
for (i in 1:length(unique(df.bray_patient$patient))) {
  
  patient <- unique(df.bray_patient$patient)[i]
  index <- grep(patient, names(final_tree_set))
  keep <- unique(c(df.bray_patient$Var1, df.bray_patient$Var2))
  
  tree <- final_tree_set[[index]]
  n <- Ntip(tree)
  m <- Nnode(tree)
  
  x <- data.frame(distRoot(
    tree,
    tips = "all",
    method = c("patristic")
  ))
  y <- data.frame( x[which(rownames(x) %in% keep),], row.names = rownames(x)[which(rownames(x) %in% keep)])
  rootdist[[i]] <- y
  phylodist[[i]] <- dist.nodes(final_tree_set[[index]])
  row.names(phylodist[[i]]) <- colnames(phylodist[[i]]) <- c(tree$tip.label,c((n+1):(n+m)))
  y <- phylodist[[i]]
  phylodist[[i]] <- phylodist[[i]][which(rownames(phylodist[[i]]) %in% c(df.bray_patient$Var1)), which(colnames(phylodist[[i]]) %in% c(df.bray_patient$Var2))]
  
  if ( length(phylodist[[i]]) == 0) {
    next
  }
  
  #phylodist[[i]][!lower.tri(phylodist[[i]])] <- NA
  phylodist[[i]] <- melt(phylodist[[i]], measure.vars = var1)
  x <- phylodist[[i]]
  phylodist[[i]] <- phylodist[[i]][complete.cases(phylodist[[i]]), ]
  
}

phylodist <- do.call("rbind", phylodist)
colnames(phylodist)[3] <- "phylodist"
rootdist <- do.call("rbind", rootdist)
colnames(rootdist)[1] <- "rootdist"
phylodist$Var1 <- as.character(phylodist$Var1)
phylodist$Var2 <- as.character(phylodist$Var2)
phylodist$comp <- paste(pmin(phylodist$Var1,phylodist$Var2), pmax(phylodist$Var1,phylodist$Var2), sep = ":")

# combine
distComp <- merge(df.bray_patient, phylodist)
distComp <- distComp[complete.cases(distComp),]
distComp <- distComp[which(distComp$BC!=0),]
distComp <- distComp[!duplicated(distComp),]

# add in if comparison is between wgd v diploid
distComp$Var1 <- substr(distComp$Var1, 7, nchar(as.character(distComp$Var1)))
distComp$Var2 <- substr(distComp$Var2, 7, nchar(as.character(distComp$Var2)))
distComp$WGD1 <- ifelse(distComp$Var1 %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$WGD==TRUE)], TRUE, FALSE )
distComp$WGD2 <- ifelse(distComp$Var2 %in% metadata_tsm.G_sample$Id[which(metadata_tsm.G_sample$WGD==TRUE)], TRUE, FALSE )
distComp$compWGD <- paste(distComp$WGD1, distComp$WGD2, sep = " v ")
distComp$compWGD <- ifelse(distComp$compWGD=="FALSE v FALSE", "diploid v diploid", 
                               ifelse(distComp$compWGD=="TRUE v TRUE", "WGD v WGD", "diploid v WGD"))

#plot across patient, normalise phylodist by maximum dist per patient
maxPhylodist <- aggregate(distComp$phylodist, list(distComp$patient), FUN=max) 
maxBCdist <- aggregate(distComp$BC, list(distComp$patient), FUN=max) 
colnames(maxPhylodist) <- c("patient","maxPhylodist")
colnames(maxBCdist) <- c("patient","maxBCdist")
distComp <- merge(distComp, maxPhylodist)
distComp <- merge(distComp, maxBCdist)


Npairwise <- aggregate(distComp$phylodist, list(distComp$patient), FUN=length) 
colnames(Npairwise) <- c("patient","Npairwise")
distComp <- merge(distComp, Npairwise)

distComp$normPhylodist <- distComp$phylodist/distComp$maxPhylodist
distComp$normBCdist <- distComp$BC/distComp$maxBCdist

distComp$MSI <- ifelse(distComp$patient %in% metadata_tsm.G_sample$patient[which(metadata_tsm.G_sample$MSI=="MSS")], "MSS", "MSI")

distComp$type1 <- substr(str_after_nth(distComp$Var1,"_",2),1,1)
distComp$type2 <- substr(str_after_nth(distComp$Var2,"_",2),1,1)
distComp$Gonly <- ifelse(distComp$type1=="G" & distComp$type2=="G",TRUE, FALSE)

jpeg('tempfig.jpeg', width = (20), height = (20), units = 'cm', res = 300)
ggplot(distComp[which(distComp$Npairwise>3),], aes(x=(normPhylodist), y=(normBCdist), color=MSI)) +
  geom_smooth(method='lm', formula= y~x, aes(color=MSI, fill=MSI)) +
  geom_point(shape=21, aes(fill=MSI), color="black", size=4, alpha=0.8) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'right', label.y = 'bottom', size=8, ) +
  scale_color_manual(values =  c("#339999", "#CC6600")) +
  scale_fill_manual(values =  c("#339999", "#CC6600")) +
  ylab("Average Bray-Curtis diversity per patient") +
  xlab("Total SNV count") +
  theme_custom() +
  theme(legend.text = element_text(size=28), legend.title = element_blank(), legend.position = "top")
dev.off()


# for patient with normal, consider root dist
ggplot(distComp[which(distComp$Npairwise>3),], aes(x=(normPhylodist), y=(normBCdist))) +
  geom_point(alpha=0.6, shape=21, fill="#000099",size=3) +
  geom_smooth(method='lm', formula= y~x, color="#000099") +
  scale_y_continuous(name="Normalised pairwise Bray-Curtis", limits = c(0,1.1), breaks = seq(0,1,0.2)) +
  xlab("Normalised phylogenetic distance") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "none")

# plot per patient
fig <- list()
p <- list()
c <- list()
for (i in 1:length(unique(distComp$patient))) {
  patient <- unique(distComp$patient)[i]
  wd <- distComp[which(distComp$patient==patient),]
  if(nrow(wd)<0) {
    next
  }
  fig[[i]] <- ggplot(wd, aes(x=normPhylodist, y=BC)) +
    geom_point() +
    ggtitle(patient) +
    geom_smooth(method='lm', formula= y~x) +
    scale_x_continuous(name="Normalised phylogenetic Distance", limits = c(0,1), breaks = seq(0,1,0.5)) +
    scale_y_continuous(name="Pairwise Bray-Curtis", limits = c(0,1), breaks = seq(0,1,0.5)) +
    ggpmisc::stat_poly_eq(formula = y ~ x, 
                          aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                          parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
    theme_custom() +
    theme(axis.title = element_blank())
  if (nrow(distComp[which(distComp$patient==patient),])==1) {
    p[[i]] <-NA
    next
  } else {
    mod <- lm(BC~normPhylodist, distComp[which(distComp$patient==patient),])
    p[[i]] <- setNames(summary(mod)$coefficients[2,4], patient)
    c[[i]] <- setNames(summary(mod)$coefficients[2,1], patient)
  }
}

fig[[26]]

#jpeg('tempfig.jpeg', width = (3*37.795*20), height = (3*37.795*15))
jpeg('tempfig.jpeg', width = (50), height = (70), units = 'cm', res = 300)
annotate_figure(
ggarrange(fig[[1]],fig[[2]],fig[[4]],fig[[5]],
          fig[[6]],fig[[7]],fig[[8]],fig[[9]],fig[[10]],
          fig[[11]],fig[[12]],fig[[13]],fig[[14]],fig[[15]],
          fig[[16]],fig[[17]],fig[[19]],fig[[20]],
          fig[[22]],fig[[23]],fig[[24]],fig[[25]],
          fig[[26]], nrow=6, ncol=4, common.legend=TRUE, legend='bottom'), 
left = text_grob("Pairwise Bray-Curtis", size=28, rot = 90), bottom = text_grob("Normalised phylogenetic distance", size=28))
#ggarrange(plotlist=fig,ncol = 5, common.legend=TRUE, legend='bottom')
dev.off()

p <- data.frame(p=unlist(p))
p$patient <- rownames(p)
c <- data.frame(c=unlist(c))
c$patient <- rownames(c)
y <- merge(c,p)
y$adj <- p.adjust(y$p, "fdr")
y$sig <- ifelse(y$p<0.1, TRUE, FALSE)
y$adjsig <- ifelse(y$adj<0.05, TRUE, FALSE)

x <- data.frame(table(distComp$patient))
y$nComp <- x[match(y$patient, x$Var1), 2]

ggplot(y, aes(x=reorder(patient,nComp), y=nComp, fill=sig)) +
  geom_bar(stat="identity")

# variance in ploidy
x <- metadata_tsm.G_sample[,which(colnames(metadata_tsm.G_sample) %in% c("patient","ploidy"))]
x <- x[complete.cases(x),]
x <- aggregate(x$ploidy, list(x$patient), FUN=sd)
y$ploidyVar <- x[match(y$patient, x$Group.1),2]
y$ploidy <- metadata_tsm.G_patient[match(y$patient, metadata_tsm.G_patient$patient), 3]
ggplot(y, aes(x=p, y=(ploidy))) +
  geom_point()

### phylo dist and diff PGA ----
distComp$var1New <- str_before_last(distComp$Var1,"_")
distComp$var1New <- sub("_", ".", distComp$var1New)

distComp$var2New <- str_before_last(distComp$Var2,"_")
distComp$var2New <- sub("_", ".", distComp$var2New)

output <- list()
for (i in 1:nrow(distComp)) {
  output[[i]] <- abs(pga$pga[which(rownames(pga)==distComp$var1New[i])] - pga$pga[which(rownames(pga)==distComp$var2New[i])])
}

output[lengths(output) == 0] <- NA

distComp$diffPGA <- unlist(output)


#jpeg('tempfig.jpeg', width = (20), height = (20), units = 'cm', res = 300)
ggplot(distComp, aes(x=BC, y=(diffPGA))) +
  geom_point(alpha=0.6, shape=21, fill="#000099",size=3) +
  geom_smooth(method='lm', formula= y~x, color="#000099") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "none")
#dev.off()

### library size diff and BC ----

braycurtis_dist <- phyloseq::distance(input_tsm_decon, method="bray")
braycurtis_dist <- as.matrix(braycurtis_dist)

# calculate pairwise distance
sub_dist <- list()
patient_all <- sample_data(input_tsm_decon)$patient

for (i in 1:length(unique(patient_all)) ) { 
  patient <- unique(patient_all)[i]
  index <- which(patient_all == patient)
  sample_group <- sample_names(input_tsm_decon)[index]
  sub_dist[[patient]] <- braycurtis_dist[sample_group, sample_group]
  sub_dist[[patient]][!lower.tri(sub_dist[[patient]])] <- NA
}

braygroups <- melt(sub_dist)
df.bray_patient <- braygroups[complete.cases(braygroups), ]
df.bray_patient$patient <- factor(df.bray_patient$L1, levels=names(sub_dist))

df.bray_patient$Var1 <- str_before_first(as.character(df.bray_patient$Var1), ".k")
df.bray_patient$Var2 <- str_before_first(as.character(df.bray_patient$Var2), ".k")
colnames(df.bray_patient)[3] <- "BC"

# make comparison name variable
df.bray_patient$comp <- paste(pmin(df.bray_patient$Var1,df.bray_patient$Var2), pmax(df.bray_patient$Var1,df.bray_patient$Var2), sep = ":")

# get library size
output <- list()
for (i in 1:nrow(df.bray_patient)) {
  wd <- sample_data(input_tsm)
  x <- wd$librarySize[which(wd$Id == substr(df.bray_patient$Var1[i], 7, nchar(df.bray_patient$Var1[i])))]
  y <- wd$librarySize[which(wd$Id == substr(df.bray_patient$Var2[i], 7, nchar(df.bray_patient$Var2[i])))]
  output[[i]] <- abs(x-y)
}
df.bray_patient$librarySizeDiff <- unlist(output)

jpeg('tempfig.jpeg', width = (20), height = (20), units = 'cm', res = 300)
ggplot(df.bray_patient, aes(y=librarySizeDiff, x=BC)) +
  geom_point(alpha=0.6, shape=21, fill="#000099",size=3) +
  geom_smooth(method='lm', formula= y~x, color="#000099") +
  #scale_y_continuous(name="Normalised pairwise Bray-Curtis", limits = c(0,1.1), breaks = seq(0,1,0.2)) +
  xlab("Pairwise Bray-Curtis distance") +
  ylab("Pairwise difference in library size") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "none")
dev.off()

### phylotree and heatmap ----

unique(sample_data(input_tsm_decon_rem)$patient)

# phylotree
final_tree_set <- readRDS("~/Documents/Microbial/data/final_tree_set.rds")

i <- 30
patient <- names(final_tree_set)[i]
plot(final_tree_set[[i]])
i <- final_tree_set[[i]]
taxa_names(i) <- str_after_nth(taxa_names(i), "_", 2)
taxa_names(i) <- str_before_nth(taxa_names(i), "_", 2)
i <- drop.tip(i, length(i[["tip.label"]]))
#i <- drop.tip(i,1)

cols <- substr(taxa_names(i),1,1)
annoCol<-data.frame(Region=c(A="#cc0000", B="#6699cc", C="#669900", D="#7570B3",
                       E="#ffccff", F="#cc6600", G="#669999", H="#ffcc33", W="#0066cc"))
cols <- annoCol[match(cols, row.names(annoCol)),1]

jpeg('tempfig.jpeg', width = (10), height = (10), units = 'cm', res = 300)
plot(i, use.edge.length=F, no.margin = F,
       font=2, main = "Phylogenetic tree",
       tip.color=cols)
dev.off()


# BC hm
braycurtis_dist <- phyloseq::distance(subset_samples(input_tsm_decon_rem, patient =="C561") , method="bray")
braycurtis_dist <- as.matrix(braycurtis_dist)
colnames(braycurtis_dist) <- str_after_nth(colnames(braycurtis_dist), "_", 2)
colnames(braycurtis_dist) <- str_before_nth(colnames(braycurtis_dist), "_", 2)
rownames(braycurtis_dist) <- str_after_nth(rownames(braycurtis_dist), "_", 2)
rownames(braycurtis_dist) <- str_before_nth(rownames(braycurtis_dist), "_", 2)
braycurtis_dist[braycurtis_dist==0] <- NA

anno<-data.frame(row.names=rownames(braycurtis_dist), Region=substr(rownames(braycurtis_dist),1,1))
annoCol<-list(Region=c(A="#cc0000", B="#6699cc", C="#669900", D="#7570B3",
                      E="#ffccff", F="#cc6600", G="#669999", H="#ffcc33", W="#0066cc"))

#library("pheatmap")

jpeg('tempfig.jpeg', width = (11), height = (7), units = 'cm', res = 300)
pheatmap(braycurtis_dist, 
         #cutree_rows = 6, 
         cluster_rows = F, cluster_cols = F,
         treeheight_col = 0,
         na_col="white", colorRampPalette(rev(brewer.pal(6, "PuBuGn")))(5),
         annotation_colors = annoCol,
         annotation_row = anno, legend = T, annotation_legend = F,
         main = "Pairwise Bray-Curtis")
dev.off()


### all species hm per patient ----
wd <- subset_samples(input_tsm_decon_rem, patient=="C522")
wd <- wd %>%
  aggregate_taxa("Species") %>%
  microbiome::transform(transform = "compositional") %>%
  transform_sample_counts(function(x) log(1 + x))

wd <- data.frame(otu_table(wd))
colnames(wd) <- str_after_nth(colnames(wd), "_", 2)
colnames(wd) <- str_before_nth(colnames(wd), "_", 2)
wd <- wd[rowSums(wd)>0,]
wd <- t(wd)

anno<-data.frame(row.names=rownames(wd), Region=substr(rownames(wd),1,1))
annoCol<-list(Region=c(A="#cc0000", B="#6699cc", C="#669900", D="#7570B3",
                      E="#ffccff", F="#cc6600", G="#669999", H="#ffcc33", W="#0066cc"))

jpeg('tempfig.jpeg', width = (20), height = (7), units = 'cm', res = 300)
pheatmap::pheatmap(wd, show_colnames = F, 
                   cutree_rows = 5, cluster_rows = T, cluster_cols = F,
                   treeheight_col = 0,
                   border_color = "NA",
                   clustering_method = "ward.D2",
                   color = hcl.colors(50, "Purple-Blue", rev = T),
                   annotation_colors = annoCol,
                   annotation_row = anno,
                   main = "Species abundance")
dev.off()






# plot per patient
fig <- list()
p <- list()
c <- list()
for (i in 1:length(unique(distComp$patient))) {
  patient <- unique(distComp$patient)[i]
  fig[[i]] <- ggplot(distComp[which(distComp$patient==patient),], aes(x=phylodist, y=BC)) +
    geom_point() +
    #ggtitle(patient) +
    geom_smooth(method='lm', formula= y~x) +
    scale_x_continuous(name="Phylogenetic Distance", n.breaks = 4) +
    ylab("Pairwise Bray-Curtis") +
    ggpmisc::stat_poly_eq(formula = y ~ x, 
                          aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                          parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
    theme_custom() 
  if (nrow(distComp[which(distComp$patient==patient),])==1) {
    p[[i]] <-NA
    next
  } else {
    mod <- lm(BC~normPhylodist, distComp[which(distComp$patient==patient),])
    p[[i]] <- setNames(summary(mod)$coefficients[2,4], patient)
    c[[i]] <- setNames(summary(mod)$coefficients[2,1], patient)
  }
}

jpeg('tempfig.jpeg', width = (20), height = (7), units = 'cm', res = 300)
fig[[2]]
dev.off()
### variance of DA species per patient ----


dds_input <- subset_samples(input_tsm_decon_rem_car, patient=="C562")
sample_data(dds_input)$region <- str_after_first(sample_data(dds_input)$region, "_")

dds <- phyloseq_to_deseq2(dds_input, ~ region)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))} # calculate geometric means prior to estimate size factors: necessary because lots of zeros
geoMeans <- apply(counts(dds), 1, gm_mean)
deseq2 <- estimateSizeFactors(dds, geoMeans = geoMeans)
deseq2 <- DESeq(deseq2, fitType="local")
plotDispEsts(deseq2)

resultsNames(deseq2)
res = results(deseq2, contrast = list("region_D_vs_A")) #MSI is base
res = res[order(res$padj, na.last=NA), ]
sigtab2 = res[(res$padj < 0.05), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(dds_input)[rownames(sigtab2), ], "matrix"))

# plot
sigtabgen = subset(sigtab2, !is.na(Species))

# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$name <- paste(substring(sigtabgen$Genus,1,1),sigtabgen$Species,sep='. ')

#jpeg('tempfig.jpeg', width = (30), height = (20), units = "cm", res = 300)
ggplot(sigtabgen, aes(y=reorder(name,log2FoldChange) , x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  ylab(NULL) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  theme_custom() +
  theme(legend.text = element_text(size=24), legend.title = element_text(size=24),
        axis.text.y = element_text(size=24))
#dev.off()

### variance of contam vs noncontam genera ----

# get otu table at genus level and make compositional
x <- input_tsm_decon_rem_car %>%
  aggregate_taxa("Species") %>%
  transform_sample_counts(function(ASV) (ASV/sum(ASV) *100) )
otu <- data.frame(otu_table(x))

# for each patient, get variance per genus
output <- list()
for ( i in 1:length(unique(sample_data(x)$patient)) ) {
  wd <- data.frame(otu[, which(sample_data(x)$patient == unique(sample_data(x)$patient)[i])])
  output[[i]] <- data.frame(patient = unique(sample_data(x)$patient)[i], 
                            Variance = apply(wd, 1, var), 
                            St.dev = apply(wd, 1, sd), 
                            mean = apply(wd, 1, mean))
  output[[i]]$species <- rownames(output[[i]])
  output[[i]]$COV <- (output[[i]]$St.dev / output[[i]]$mean)
}

output <- do.call(rbind, output)
output <- na.omit(output)
tax <- data.frame(tax_table(input_tsm_decon_rem_car))
output$genus <- tax[match(output$species, tax$Species), 6]
output$name <- paste(substring(output$genus,1,1),output$species,sep='. ')
drop <- list()
for (i in 1:nrow(output)) {
    drop[[i]] <- grepl("_",output$name[i])
}
output$drop <- unlist(drop)
output <- output[which(output$drop==FALSE),]
x <- aggregate(output$mean, list(output$name), median)
keep1 <- x$Group.1[which(x$x>0.5)]
output_sub <- output[which(output$name %in% keep1),]
x <- aggregate(output_sub$COV, list(output_sub$name), mean)
x <- data.frame(table(output_sub$name))
keep2 <- x$Var1[which(x$Freq>10)]
output_sub <- output[which(output$name %in% keep2),]

output_sub$WGD <- metadata_tsm.G_patient[match(output_sub$patient, metadata_tsm.G_patient$patient),4]
output_sub$MSI <- metadata_tsm.G_patient[match(output_sub$patient, metadata_tsm.G_patient$patient),5]

jpeg('tempfig.jpeg', width = (30), height = (20), units = "cm", res = 300)
ggplot(data=output_sub, aes(x=reorder(name, COV, median, na.rm=T), y=COV)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(size=3, shape=21, aes(fill=patient), width = 0.1, height=0, alpha=0.8) +
  coord_flip() +
  theme_custom() +
  theme(axis.title.y = element_blank(),
        legend.text = element_text(size=28), legend.title = element_text(size=28))
dev.off()

ggplot(data=output_sub, aes(x=MSI, y=COV, fill=MSI)) +
  geom_boxplot(alpha=0.6, outlier.size = 0) +
  geom_jitter(shape=21, alpha=0.5, width = 0.3, size=3) +
  geom_text_repel(aes(label=name), max.overlaps = 12) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("MSI","MSS")), size=9) +
  scale_fill_manual(values =  c("#339999", "#CC6600")) +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5), legend.position = "none")


### ======================================================= ----
### DE ----
otu <- data.frame(otu_table(input_tsm_decon_rem_car))
tax <- data.frame(tax_table(input_tsm_decon_rem_car))
x <- otu[which(rownames(otu)=="851"),]
hist(as.numeric(x), breaks = 100)
highFN <- str_before_first(str_after_first(colnames(x)[which(x>0)], "_") , ".k")
lowFN <- str_before_first(str_after_first(colnames(x)[which(x==0)], "_") , ".k")

exp <- readRDS("~/Documents/CNA/Data/EPICC/EPICC_expression/All_EPICC_counts.rds")
keep <- readRDS("~/Documents/CNA/Data/EPICC/EPICC_expression/gene_clustering_and_id_conversion.rds")


meta <- metadata_tsm.G_sample[which(metadata_tsm.G_sample$microbiomeData==TRUE),]
meta <- meta[which(meta$patient=="C562"),]
meta <- meta[!is.na(meta$pga),]
exp_sub <- exp[colnames(exp) %in% meta$sample]
exp_sub <- cbind(exp[1], exp_sub)

meta <- meta[which(meta$sample %in% colnames(exp)),]
meta$pga4group <- factor(meta$pga4group)
meta$FNhigh <- ifelse(meta$Id %in% highFN, TRUE, FALSE)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = exp_sub, 
                                      colData = meta, 
                                      design = ~ patient + FNhigh, tidy = TRUE)

dde <- DESeq2::DESeq(dds)

#res <- DESeq2::results(dde, contrast=c("FNhigh",TRUE,FALSE))
res <- DESeq2::results(dde)
res <- res[order(res$padj),]
head(res)

jpeg('tempfig.jpeg', width = (20), height = (20), units = 'cm', res = 300)
par(mfrow=c(1,1))
par(mar=c(5,5,1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=NULL, xlim=c(-3,3), cex=1, cex.axis=2, cex.lab=2, cex.main=2))
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

res.sig <- as.data.frame(res)
res.sig <- res.sig[which(res.sig$padj<0.05),]

res.sig$ensemble <- rownames(res.sig)
res.sig$gene_id <- mapIds(org.Hs.eg.db, keys = res.sig$ensemble,
                          column = c('SYMBOL'), keytype = 'ENSEMBL')
res.sig$entrez <- mapIds(org.Hs.eg.db, keys = res.sig$ensemble,
                         column = c('ENTREZID'), keytype = 'ENSEMBL')

# GO of DE genes
Up <- res.sig[which(res.sig$log2FoldChange>0),]
Down <- res.sig[which(res.sig$log2FoldChange<0),]
GOanno <- limma::goana(list(Up=Up$entrez,Down=Down$entrez))

GOanno$p.Up.adj = p.adjust(GOanno$P.Up, method = "fdr")
GOanno$p.Down.adj = p.adjust(GOanno$P.Down, method = "fdr")

data <- GOanno[which(GOanno$p.Up.adj <= 0.05),]
data$ID <- rownames(data)
simMatrix_up <- rrvgo::calculateSimMatrix(data$ID, orgdb="org.Hs.eg.db", ont = c("BP"), method="Rel")
scores <- setNames(-log10(data$p.Up.adj), data$ID)
reducedTerms_up <- rrvgo::reduceSimMatrix(simMatrix_up, scores = scores, threshold=0.6, orgdb="org.Hs.eg.db")
parentGO_DEup <- data.frame(tapply(reducedTerms_up$score, reducedTerms_up$parentTerm, FUN=sum))

data <- GOanno[which(GOanno$p.Down.adj <= 0.05),]
data$ID <- rownames(data)
simMatrix_down <- rrvgo::calculateSimMatrix(data$ID, orgdb="org.Hs.eg.db", ont = c("BP"), method="Rel")
scores <- setNames(-log10(data$p.Down.adj), data$ID)
reducedTerms_down <- rrvgo::reduceSimMatrix(simMatrix_down, scores = scores, threshold=0.6, orgdb="org.Hs.eg.db")
parentGO_DEdown <- data.frame(tapply(reducedTerms_down$score, reducedTerms_down$parentTerm, FUN=sum))

library(clusterProfiler)
geneList <- res.sig$log2FoldChange
names(geneList) <- res.sig$gene_id
geneList <- sort(geneList, decreasing = TRUE)

gene <- names(geneList)
bitr_kegg(gene, fromType = "kegg", toType = "Module", organism = "hsa")

set.seed(123456)
kk2 <- gseKEGG(geneList     = geneList,
               organism     = "human",
               keyType       = "ncbi-geneid",
               nPerm = 1000,
               minGSSize    = 2,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

### permanova per patient ----

patients <- unique(sample_data(input_tsm_decon_rem_car)$patient)
wd <- subset_samples(input_tsm_decon_rem_car, patient==patients[1])
regionFreq <- data.frame(table(sample_data(wd)$region))
keep <- regionFreq$Var1[which(regionFreq$Freq>0)]

#if (length(keep)<1) {
#  next
#}

wd <- subset_samples(wd, region %in% keep)
wd <- microbiome::transform(wd, "compositional")

tax <- data.frame(tax_table(wd))
tax$name <- paste(substr(tax$Genus,1,1), tax$Species, sep = ". ")

otu <- abundances(wd)
meta <- meta(wd)
#meta$region[which(meta$region %in% c("C562_B","C562_D"))] <- "C562_BD"
meta$region <- factor(meta$region)

set.seed(123)
permanova <- adonis(t(otu) ~ region,
                    data = meta, permutations=99, method = "bray")

print(as.data.frame(permanova$aov.tab)["region", "Pr(>F)"])

permutest(betadisper(vegdist(t(otu), method = "bray"), meta$region), pairwise = TRUE) # Null hypothesis = no difference in dispersion between groups.







permanova_pairwise(t(otu), meta$region,
                   permutations = 999,
                   method = "bray",
                   padj = "bonferroni"
)

C(meta$region, contr.sum) 

coef <- coefficients(permanova)["region2",]
names(coef) <- tax[match(rownames(tax), names(coef)),8]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa", )



### RNA is plastic, does it associated with microbiome? ----
exp <- readRDS("~/Documents/CNA/Data/EPICC/EPICC_expression/All_EPICC_counts.rds")
exp <- exp[which(rowSums(exp[-1])>0),]

### CMS----
devtools::install_github("Lothelab/CMScaller")
library(Biobase)
library(CMScaller)

par(mfrow=c(1,2))
exp <- readRDS("~/Documents/CNA/Data/EPICC/EPICC_expression/All_EPICC_counts.rds")
rownames(exp) <- exp$GeneID
exp <- exp[-1]
CMS <- CMScaller(exp, RNAseq=TRUE, doPlot=TRUE, rowNames = "ensg")

head(res)
CMSdf <- data.frame(CMS)

x <- input_tsm_decon_rem_car %>%
  aggregate_taxa("Species") %>%
  transform_sample_counts(function(x) log(1 + x))

pca_input <- t(data.frame(otu_table(x)))
pca <- prcomp(pca_input)
percentVar <- round(100 * summary(pca)$importance[2,])
df <- cbind(sample_data(x),pca$x)
df$CMS <- CMSdf[match(df$sample, rownames(CMSdf)), 1]

ggplot(df[!is.na(df$CMS),], aes(x=PC1, y=PC2, color=CMS)) + 
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_custom() +
  theme(legend.position = "right",
        panel.grid.major = element_line(size=0.25),
        panel.grid.minor = element_line(size=0.25),
        legend.text = element_text(size=28),  
        legend.title = element_blank()
  )

x <- input_tsm_decon_rem_car
sample_data(x)$CMS <- CMSdf[match(df$sample, sample_data(x)$sample), 1]
x <- subset_samples(x, !is.na(CMS))
ordination <- ordinate(x, method="PCoA", distance="bray")
plot_ordination(x, ordination, color="CMS") + 
  theme(aspect.ratio=1, legend.position = "none") 

### CMS DA ----

x <- input_tsm_decon_rem_car %>%
  aggregate_taxa("Genus") 

sample_data(x)$CMS <- CMSdf[match(sample_data(x)$sample, rownames(CMSdf)), 1]
sample_data(x)$MSI <- unlist(metadata_tsm.G_patient[match(sample_data(x)$patient, metadata_tsm.G_patient$patient), 5])

sample_data(x)$CMS1 <- ifelse(sample_data(x)$CMS=="CMS1", TRUE, FALSE)
sample_data(x)$CMS2 <- ifelse(sample_data(x)$CMS=="CMS2", TRUE, FALSE)
sample_data(x)$CMS3 <- ifelse(sample_data(x)$CMS=="CMS3", TRUE, FALSE)
sample_data(x)$CMS4 <- ifelse(sample_data(x)$CMS=="CMS4", TRUE, FALSE)

x <- subset_samples(x, !is.na(CMS))
dds <- phyloseq_to_deseq2(x, ~   patient + CMS1)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))} # calculate geometric means prior to estimate size factors: necessary because lots of zeros
geoMeans <- apply(counts(dds), 1, gm_mean)
deseq2 <- estimateSizeFactors(dds, geoMeans = geoMeans)
deseq2 <- DESeq(deseq2, fitType="local")
plotDispEsts(deseq2)

resultsNames(deseq2)
#res = results(deseq2, contrast = c("CMS","CMS4","CMS3")) #MSI is base
res = results(deseq2)
res = res[order(res$padj, na.last=NA), ]
sigtab2 = res[(res$padj < 0.01), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(x)[rownames(sigtab2), ], "matrix"))

sigtabgen = subset(sigtab2, !is.na(Species))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$name <- paste(substring(sigtabgen$Genus,1,1),sigtabgen$Species,sep='. ')

#jpeg('tempfig.jpeg', width = (30), height = (20), units = "cm", res = 300)
ggplot(sigtabgen, aes(y=reorder(name,log2FoldChange) , x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  ylab(NULL) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  theme_custom() +
  theme(legend.text = element_text(size=24), legend.title = element_text(size=24),
        axis.text.y = element_text(size=24))
#dev.off()


### ======================================================= ----
### ======WIP============================================== ----
### phylo dist and HC ----

library("dendextend")
ps_rel_abund <- microbiome::aggregate_taxa(input_tsm_decon, "Species") %>%
  subset_samples(patient=="C542") %>%
  microbiome::transform(transform = "compositional")

ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abund))
ps_rel_otu <- t(ps_rel_otu)
meta <- data.frame(phyloseq::sample_data(ps_rel_abund))

bc_dist <- ps_rel_otu %>% 
  vegan::vegdist(method = "bray") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram() 

dendextend::labels(bc_dist) <- meta$subregion[order.dendrogram(bc_dist)]

colors_to_use <- as.numeric(as.factor(meta$subregion[order.dendrogram(bc_dist)]))
palette <- sample(c("#666600","#CC9933","#CC0000","#993333","#FF9933",
                    "#660066","#6600FF","#9933CC","#6666FF","#0066FF",
                    "#33CCFF","#003333","#00CC99","#003300","#339966",
                    "#CC99CC","#CCCC99","#00FF33","#FF3300","#330033", 
                    "#CC00FF","#330066","#CC9966","#6633CC","#6666CC",
                    "#00CCFF","#CC9933","#CC0000","#993333","#FF9933"), length(unique(colors_to_use)))

bc_dist <- dendextend::color_labels(bc_dist,col=palette[colors_to_use])

plot(bc_dist)

dend <- as.dendrogram(bc_dist)
labels = dend %>% labels
heights = dend %>% hang.dendrogram %>% get_leaves_attr("height")
data.frame(labels = labels,
           heights = heights)

# choose height that keeps within cluster <250000
bc_dist <- ps_rel_otu %>% 
  vegan::vegdist(method = "bray") %>%
  hclust(method = "ward.D2")

for (i in 1:length(bc_dist$height)) {
  h <- bc_dist$height[i]
  cut <- cutree(bc_dist, h=h)
  ps_rel_otu$cluster <- unlist(cut)
}


### hierarchical clustering and heatmap on top 15 genus per sample ----
library("dendextend")
ps_rel_abund <- microbiome::aggregate_taxa(input_t, "Species") %>%
  microbiome::transform(transform = "compositional")

ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abund))
ps_rel_otu <- t(ps_rel_otu)
meta <- data.frame(phyloseq::sample_data(ps_rel_abund))

bc_dist <- ps_rel_otu %>% 
  vegan::vegdist(method = "bray") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram() 

dendextend::labels(bc_dist) <- meta$patient[order.dendrogram(bc_dist)]

colors_to_use <- as.numeric(as.factor(meta$patient[order.dendrogram(bc_dist)]))
palette <- sample(c("#666600","#CC9933","#CC0000","#993333","#FF9933",
                    "#660066","#6600FF","#9933CC","#6666FF","#0066FF",
                    "#33CCFF","#003333","#00CC99","#003300","#339966",
                    "#CC99CC","#CCCC99","#00FF33","#FF3300","#330033", 
                    "#CC00FF","#330066","#CC9966","#6633CC","#6666CC",
                    "#00CCFF","#CC9933","#CC0000","#993333","#FF9933"), length(unique(colors_to_use)))

bc_dist <- dendextend::color_labels(bc_dist,col=palette[colors_to_use])

plot(bc_dist)


### consider the other distance measures ----

    # remove the two distance-methods that require a tree, and the generic custom method that requires user-defined distance arguments
dist_methods <- unlist(distanceMethodList)
dist_methods <- dist_methods[-(1:3)]
dist_methods = dist_methods[-which(dist_methods=="ANY")]

    # Loop through each distance method, save each plot to a list
phyloseqOb <- bacteria_tsm
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(phyloseqOb, method=i)
  # Calculate ordination
  iMDS  <- ordinate(phyloseqOb, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(phyloseqOb, iMDS, color="subregion")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=subregion))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics")

jpeg('tempfig.jpeg', width = (3*37.795*10), height = (3*37.795*10))
p
dev.off()


### heatmap of relative abundance of prevalent phyla ----

old <- prune_taxa(names(sort(taxa_sums(bacteria_prev.subset),TRUE)[1:10]), bacteria_prev.subset)
new <- prune_taxa(names(sort(taxa_sums(bacteria_vst),TRUE)[1:10]), bacteria_vst)

plot_heatmap(old, taxa.order = "Phylum", taxa.label = "Species", sample.label = "region", 
             sample.order = "patient")
plot_heatmap(new, taxa.order = "Phylum", taxa.label = "Species", sample.label = "region", 
             sample.order = "patient")


# plot abundance (not relative) in top 10 genera in carcinoma glands per sample (duplicates removed)
# does work in species because some species assigned to >1 genus
bacteria_vst %>%
  ps_arrange(desc(patient), desc(region)) %>%
  comp_barplot(
    tax_level = "Order", n_taxa = 41,
    #facet_by = "patient",
    label = "region",
    bar_outline_colour = NA,
    sample_order = "default",
    bar_width = 0.7,
    merge_other = TRUE,
    taxon_renamer = toupper
  ) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  facet_grid(facets = "patient", scales = "free",space="free" ) +
  theme(axis.text = element_text(size=5),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

comp_barplot(bacteria_vst,
  tax_level = "Genus", n_taxa = 20,
  #facet_by = "patient",
  label = "region",
  bar_outline_colour = NA,
  sample_order = "default",
  bar_width = 0.7,
  merge_other = TRUE,
  taxon_renamer = toupper
) 

    # and subset by genera/species


