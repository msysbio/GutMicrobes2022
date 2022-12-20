##################################

# R script to reproduce figures shown in van de Velde et al.:
# "Technical versus biological variability in a synthetic human gut community"

# Explanations on data see below function definitions. Note that CellScanner abundances predicted
# for experiment 2 are less reliable than for experiment 1, since no monoculture data were collected
# for the second experiment. CellScanner total counts are not scaled to cells/ml and thus not directly
# comparable to 16S total counts.

# Author: Karoline Faust

# Date: 31st August 2022

# Load required packages
# seqgroup can be installed from: https://hallucigenia-sparsa.github.io/seqgroup/
library(seqgroup)
# Optionally, networks can be sent directly to Cytoscape when installing RCy3
# BiocManager::install("RCy3")

# Set the working directory to the location of the result file or adjust the file path
base::load("Results.RData")

########### Functions ########

# Plot the time series of taxa across groups representing replicates.
# Error bars are drawn to represent standard deviation across replicates.
# abundances: a matrix where rows are taxa and columns samples
# groups: vector giving for each sample its group
# time: vector giving for each sample its time
# title: plot title
# ylab: y axis label (x axis label is time)
# colors: color map specifying a color per species
# returnData: return the reformatted data used for plotting
plotGroupRepTS<-function(abundances, groups=c(), time=c(), title="", ylab="", colors=NULL, returnData=FALSE){

  bioreplicates=unique(groups)
  time.points=unique(time)
  global.df=matrix(nrow=(length(time.points)*nrow(abundances)),ncol=4)
  colnames(global.df)=c("Taxon","Time","Mean","sd")

  row.counter=1
  # loop over taxa
  for(index in 1:nrow(abundances)){
    taxon=rownames(abundances)[index]
    # loop over time points
    for(time.point in time.points){
      tp.indices=which(time==time.point)
      abundances.tp=abundances[index,tp.indices] # get all values for the time point
      global.df[row.counter,1]=taxon
      global.df[row.counter,2]=time.point
      global.df[row.counter,3]=mean(abundances.tp,na.rm = TRUE)
      global.df[row.counter,4]=sd(abundances.tp, na.rm = TRUE)
      row.counter=row.counter+1
    } # end loop over time points
  }

  global.df=as.data.frame(global.df)
  global.df$Mean=as.numeric(global.df$Mean)
  global.df$sd=as.numeric(global.df$sd)
  global.df$Time=as.numeric(global.df$Time)
  if(!is.null(colors)){
    print("User-specified colors")
    ggplot(global.df, aes(x=Time, y=Mean, color=Taxon)) + geom_line() + theme_classic() + geom_pointrange(aes(ymin=Mean-sd, ymax=Mean+sd)) + ggtitle(title) + ylab(ylab) + xlab("Time in hours") +  scale_colour_manual(values=colors)
  }else{
    ggplot(global.df, aes(x=Time, y=Mean, color=Taxon)) + geom_line() + theme_classic() + geom_pointrange(aes(ymin=Mean-sd, ymax=Mean+sd)) + ggtitle(title) + ylab(ylab) + xlab("Time in hours")
  }
  if(returnData){
    return(global.df)
  }
}

# Plot the time series of a selected taxon across groups representing replicates.
# Error bars are drawn to represent standard deviation across replicates and significance of
# correlation with time is assessed.
# abundances: a matrix where rows are taxa and columns samples
# groups: vector giving for each sample its group
# time: vector giving for each sample its time
# species.index: index of taxon to plot
# title: plot title
# xlab: x axis label (if empty, set to: Time in hours)
# ylab: y axis label
# return: when true, plotting is suppressed and the data frame object used for plotting is returned
# showPval: when true, Pearson's r and p-value of correlation with time is shown
plotSingleRepTS<-function(abundances, groups=c(), time=c(), species.index=1, title="", xlab="", ylab="", return=FALSE, showPval=TRUE){
  bioreplicates=unique(groups)
  time.points=unique(time)
  global.df=matrix(nrow=length(time.points),ncol=3)
  colnames(global.df)=c("Time","Mean","sd")

  if(xlab==""){
    xlab="Time in hours"
  }

  row.counter=1
  # loop over time points
  for(time.point in time.points){
    tp.indices=which(time==time.point)
    abundances.tp=abundances[species.index,tp.indices] # get all values for the time point
    global.df[row.counter,1]=time.point
    global.df[row.counter,2]=mean(abundances.tp,na.rm = TRUE)
    global.df[row.counter,3]=sd(abundances.tp, na.rm = TRUE)
    row.counter=row.counter+1
  } # end loop over time points

  # prepare global data set
  # http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
  global.df=as.data.frame(global.df)
  global.df$Mean=as.numeric(global.df$Mean)
  global.df$sd=as.numeric(global.df$sd)
  global.df$Time=as.numeric(global.df$Time)

  out.cor=cor.test(global.df$Time,global.df$Mean)

  if(title!="" && showPval){
    title=paste(title,", ",sep="")
  }
  if(showPval){
    title=paste(title, "r=",round(out.cor$estimate,3),", p-value=",round(out.cor$p.value,3),sep="")
  }

  if(return){
    return(global.df)
  }else{
    ggplot(global.df, aes(x=Time, y=Mean)) + geom_line() + theme_classic() + geom_pointrange(aes(ymin=Mean-sd, ymax=Mean+sd)) + ggtitle(title) + ylab(ylab) + xlab(xlab)
  }
}

# ambr_taxa: relative taxon abundances, column names follow the template vX{separator}Y{separator}Z, where X is the vessel number, Y the time point in hours and Z the technical replicate
# selected indices, e.g. for a technical replicate
# counts: total counts reported in the same order as samples in ambr_taxa
# removeStr: a string to remove to handle other naming formats more flexibly
# startIndex: start index for the vessel, can be 1 or 7 (since there were 6 monocultures)
getAbsoluteCounts<-function(ambr_taxa, selected.indices=c(), counts=c(), separator="\\.", removeStr="", startIndex=7){
  if(length(selected.indices)==0){
    selected.indices=c(1:ncol(ambr_taxa))
  }
  # assign counts to ambr2 tech rep 1
  ambr_taxa_abs=ambr_taxa[,selected.indices]
  for(i in 1:ncol(ambr_taxa[,selected.indices])){
    colname=colnames(ambr_taxa[,selected.indices])[i]
    if(removeStr!=""){
      colname=gsub(removeStr,"",colname)
    }
    splitted=strsplit(colname,separator)
    vessel=splitted[[1]][1]
    if(removeStr==""){
      vessel=strsplit(vessel,"v")[[1]][2]
    }
    if(startIndex==7){
      vessel=paste("co",(as.numeric(vessel)-6),sep="")
    }else{
      vessel=paste("co",vessel,sep="")
    }
    tp=as.numeric(splitted[[1]][2])
    tp.indices=which(counts[,2]==tp)
    vessel.indices=which(counts[,1]==vessel)
    index=intersect(tp.indices,vessel.indices)
    count=counts[index,3]
    ambr_taxa_abs[,i]=ambr_taxa[,selected.indices[i]]*count
  }
  return(ambr_taxa_abs)
}

# Get back indices of samples selected for a given metadata item (vessel or time point).
# data: taxon abundances
# metadata: data frame containing specified metadata item
# metadataItem: the metadata item of interest (vessel or time point)
# selectedMetadataValues: the metadata values to keep for the selected metadata item
selectDataSubSet<-function(data, metadata, metadataItem="Vessel", selectedMetadataValues=c()){
 metadataValues=metadata[[metadataItem]]
 indices=c()
 for(i in 1:length(metadataValues)){
   if(metadataValues[i] %in% selectedMetadataValues){
     indices=c(indices,i)
   }
 }
 return(indices)
}

# Compute overall and technical variability as community dissimilarity
# ambr_taxa: 16S taxon abundances including technical replicates
# ambr_metadata: specification of time point, vessel and technical replicate per sample in ambr_taxa
# dis: dissimilarity/distance supported by vegdist
# reftechrep: technical replicate used to compute dissimilarity for biological replicates, can be r1, r2 or r3
# Method returns a list with entries tech.dissims and bio.dissims
getTechVersusVesselDis<-function(ambr_taxa, ambr_metadata, dis="bray", reftechrep="r1"){
  dissimMat=as.matrix(vegdist(t(ambr_taxa),method=dis))
  # expected: 3 per vessel and time point, so 3 x 6 x 6 = 108
  tech.dissims=c()
  # expected: 15 (6*5/2) per time point, so 15 * 6 = 90
  bio.dissims=c()
  # loop over time points
  for(tp in unique(ambr_metadata$Timepoint)){
    print(paste("Processing time",tp))
    tp.indices=which(ambr_metadata$Timepoint==tp)
    # loop over biological replicates
    for(vessel in unique(ambr_metadata$Vessel)){
      # grep all technical replicates of vessel and tp
      vessel.indices=which(ambr_metadata$Vessel==vessel)
      vessel.indices=intersect(tp.indices,vessel.indices)
      target.indices=c()
      for(vessel.index in vessel.indices){
        target.index=which(colnames(ambr_taxa)==ambr_metadata$X[vessel.index])
        target.indices=c(target.indices,target.index)
      }
      #print(colnames(ambr_taxa)[target.indices])
      submat=dissimMat[target.indices,target.indices]
      tech.dissims=c(tech.dissims,submat[lower.tri(submat)]) # collect dissimilarities of technical replicates per vessel and time point
    }
    # all biological replicates at current time point
    biorep.indices=intersect(tp.indices,which(ambr_metadata$Replicate==reftechrep))
    target.indices=c()
    for(biorep.index in biorep.indices){
      target.index=which(colnames(ambr_taxa)==ambr_metadata$X[biorep.index])
      target.indices=c(target.indices,target.index)
    }
    #print(colnames(ambr_taxa)[target.indices])
    submat=dissimMat[target.indices,target.indices]
    bio.dissims=c(bio.dissims,submat[lower.tri(submat)])
  }
  res=list(tech.dissims=tech.dissims,bio.dissims=bio.dissims)
  return(res)
}



######### Analysis prep ####################

# Note that the 6 community replicates are either annotated from 1 to 6 or from 7 to 12
# Time points are: 12, 25, 37, 49, 61 and 67h
ambr_taxa = results$seqData # Normalised 16S data for exp1, including the 3 technical replicates
ambr_metadata = results$seqMetadata # Specification of vessel, time point and technical replicate for exp1
metabolites = results$metabolites # Metabolite data (version 2)
heterogeneities = results$heterogeneities # Flow cytometry heterogeneities per vessel and time point for exp1
cs = results$csData # CellScanner taxon counts for exp1 (automated gating with CellScanner on 50K flow cytometry events, not yet converted into counts/mL)
counts = results$counts # Total counts for exp1 (obtained through manual gating)
# validation experiment, time points are: 0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180, 192, 204, 216, 228, 240 and 252 hours
ambr_taxa2_abs = results$seqData2 # 16S data for 3 control vessels (exp2) multiplied with flow cytometry counts obtained through manual gating
ambr_metadata2 = results$seqMetadata2 # Specification of vessel and time point for exp2
cs2 = results$csData2 # CellScanner taxon counts for exp2 (automated gating with CellScanner on 50K flow cytometry events, not yet converted into counts/mL)
cs2.metadata = results$csMetadata2 # CellScanner metadata for exp2

# Convert CS data to relative abundances
cs.proportions=seqgroup::normalize(cs)

# Convert CS2 data to relative abundances
cs2.proportions=seqgroup::normalize(cs2)

ambr2.totalcounts=colSums(ambr_taxa2_abs)
# Convert validation experiment (exp2) 16S data to relative abundances
ambr_taxa2=seqgroup::normalize(ambr_taxa2_abs)

# Select one technical replicate
tech.rep="1"
tech.replicates=c()
for(i in 1:ncol(ambr_taxa)){
  sample.name=colnames(ambr_taxa)[i]
  splitted=strsplit(sample.name,"\\.")
  tech.replicates=c(tech.replicates,splitted[[1]][3])
}
selected.indices=which(tech.replicates==tech.rep)

# Get absolute 16S data for selected technical replicate for exp1
ambr_taxa_abs=getAbsoluteCounts(ambr_taxa,selected.indices=selected.indices,counts=counts)

# time has to be numeric for sorting in ascending order per group (needed for time arrow plotting in seqPCoA)
time.raw = ambr_metadata$Timepoint
temp=strsplit(time.raw,"t")
time=matrix(nrow=length(time.raw),ncol=2,unlist(temp),byrow = TRUE)[,2]
sorted.indices=sortSamples(ambr_taxa,groups = ambr_metadata$Vessel, time = time)

# define a color scheme for time series plots
group.colors <- c("Blautia_hydrogenotrophica" = "#ED7D31", "Bacteroides_thetaiotaomicron" = "#A5A5A5", "Roseburia_intestinalis" ="#4472C3", "Collinsella_aerofaciens" = "#FFBF01")
group.colors.cs=group.colors
group.colors.cs[["Unknown"]]="black"

###### Time series plots (Figure 2) #########

xlab="Time in hours"
ylab="Relative abundance"

# Relative abundances with 16S for selected tech replicate, skipping Prevotella and Pseudomonas
title="16S"
seq.data = plotGroupRepTS(ambr_taxa[1:4, selected.indices],groups=ambr_metadata$Vessel[selected.indices],time=as.numeric(time[selected.indices]), ylab=ylab, title=title, colors = group.colors, returnData = TRUE)
p1=ggplot(seq.data, aes(x=Time, y=Mean, color=Taxon)) + geom_line() + theme_classic() + geom_pointrange(aes(ymin=Mean-sd, ymax=Mean+sd)) + ggtitle(title) + ylab(ylab) + xlab(xlab) +  scale_colour_manual(values=group.colors)
plot(p1)

# Absolute abundances with 16S for selected tech replicate, skipping Prevotella and Pseudomonas
title="16S counts"
seq.data.abs = plotGroupRepTS(ambr_taxa_abs[1:4, ],groups=ambr_metadata$Vessel[selected.indices],time=as.numeric(time[selected.indices]), ylab=ylab, title=title, colors = group.colors, returnData = TRUE)
p2=ggplot(seq.data.abs, aes(x=Time, y=Mean, color=Taxon)) + geom_line() + theme_classic() + geom_pointrange(aes(ymin=Mean-sd, ymax=Mean+sd)) + ggtitle(title) + ylab("Cells/ml") + xlab(xlab) +  scale_colour_manual(values=group.colors)
plot(p2)

# Relative abundances with CellScanner
title="CellScanner"
cs.data=plotGroupRepTS(cs.proportions,groups=ambr_metadata$Vessel,time=as.numeric(time[selected.indices]), returnData = TRUE)
p3=ggplot(cs.data, aes(x=Time, y=Mean, color=Taxon)) + geom_line() + theme_classic() + geom_pointrange(aes(ymin=Mean-sd, ymax=Mean+sd)) + ggtitle(title) + ylab(ylab) + xlab(xlab) +  scale_colour_manual(values=group.colors.cs)
plot(p3)

# Metabolite data
title="Metabolites"
met.data=plotGroupRepTS(t(metabolites),groups=ambr_metadata$Vessel,time=as.numeric(time[selected.indices]), returnData = TRUE)
colnames(met.data)=c("Metabolite","Time","Mean","sd")
p4=ggplot(met.data, aes(x=Time, y=Mean, color=Metabolite)) + geom_line() + theme_classic() + geom_pointrange(aes(ymin=Mean-sd, ymax=Mean+sd)) + ggtitle(title) + ylab("mg/L") + xlab(xlab)
plot(p4)

###### Ordination (Figure 4A, 4B and 4C) #########

arrowFactor=0.3
xmax=0.32
# create a factor that combines vessel and time point
tp.vessel.combi=c()
tp.vessel.combi.cs=c()
for(i in 1:length(ambr_metadata$Replicate)){
  item=paste(ambr_metadata$Vessel[i], ambr_metadata$Timepoint[i],sep="_")
  tp.vessel.combi=c(tp.vessel.combi,item)
  if(i <= nrow(cs)){
    tp.vessel.combi.cs=c(tp.vessel.combi.cs,paste("V",cs[i,2],cs[i,1],sep="_"))
  }
}

# sum metabolites according to location in ordination plot for a clearer plot
separate.metabolites=c(4,5)
formpyr=c(3,6)
sugars=c(1,2)
scfa=c(7,8,9,10)
grouped.metabolites = metabolites[,separate.metabolites]
grouped.metabolites[["SUGARS"]]=rowSums(metabolites[,sugars])
grouped.metabolites[["SCFAs"]]=rowSums(metabolites[,scfa])
grouped.metabolites[["FORMIC/PYRUVIC"]]=rowSums(metabolites[,formpyr])

# uncomment to abbreviate taxon names in ordination plot
#abbreviated.taxonnames=c("RI","BH","BT","CA","PC","Pseudomonas")
#rownames(ambr_taxa)=abbreviated.taxonnames

# ordination 16S selected tech rep
seqPCoA(ambr_taxa[,selected.indices],groups=ambr_metadata$Vessel[selected.indices], main=paste("PCoA 16S, tech replicate",tech.rep), metadata=metabolites, topTaxa=nrow(ambr_taxa), topMetadata = ncol(metabolites), metadataFactor = 0.25, labels = ambr_metadata$Timepoint[selected.indices], arrowFactor = arrowFactor, time=TRUE,  xlim=c(-0.41,xmax), ylim=c(-0.4,0.4), hideGroupLegend=TRUE)
# grouped metabolites, group legend visible, taxon arrows in red:
seqPCoA(ambr_taxa[,selected.indices],groups=ambr_metadata$Vessel[selected.indices], main=paste("PCoA 16S, tech replicate",tech.rep), metadata=grouped.metabolites, topTaxa=nrow(ambr_taxa), taxonColor="red", topMetadata = ncol(metabolites), metadataFactor = 0.25, labels = ambr_metadata$Timepoint[selected.indices], arrowFactor = arrowFactor, time=TRUE,  xlim=c(-0.41,xmax), ylim=c(-0.4,0.4), hideGroupLegend=FALSE)

# ordination across all technical replicates
seqPCoA(ambr_taxa[,sorted.indices],groups=ambr_metadata$Vessel[sorted.indices], main="PCoA 16S with 3 technical replicates", errorbars=tp.vessel.combi[sorted.indices], topTaxa=0, errorbarLabels = c(), arrowFactor = 0.2,ellipseConf = 0.75,xlim=c(-0.3,0.2),ylim=c(-0.4,0.4),sizes=as.numeric(time[sorted.indices]),
        size.legend = FALSE, hideGroupLegend=TRUE)

# ordination CS
seqPCoA(cs.proportions,groups=ambr_metadata$Vessel[selected.indices], time=TRUE, main="PCoA CS", metadata=metabolites, topTaxa=0, topMetadata = ncol(metabolites), metadataFactor = 0.25, labels = ambr_metadata$Timepoint[selected.indices], arrowFactor = arrowFactor, xlim=c(-0.35,xmax), ylim=c(-0.4,0.42), hideGroupLegend=TRUE)

# Supplement: ordination across the three technical replicates for the last four time points
lastfourtimepoints=c(7:18,25:36,43:54,61:72,79:90,97:108)
sizes=as.numeric(time[sorted.indices[lastfourtimepoints]])
seqPCoA(ambr_taxa[,sorted.indices[lastfourtimepoints]],groups=ambr_metadata$Vessel[sorted.indices[lastfourtimepoints]], main="PCoA CS with the four last time points", errorbars=tp.vessel.combi[sorted.indices[lastfourtimepoints]], topTaxa=0, errorbarLabels = c(),ellipseConf = 0.75,
        xlim=c(-0.35,0.4),ylim=c(-0.55,0.4),
        sizes=sizes,
        size.legend = FALSE, hideGroupLegend=TRUE)


######## Variability and heterogeneity (Figure 4D and Figure 6) ######################

res1=getTechVersusVesselDis(ambr_taxa = ambr_taxa,ambr_metadata = ambr_metadata)
res2=getTechVersusVesselDis(ambr_taxa = ambr_taxa,ambr_metadata = ambr_metadata,reftechrep = "r2")
res3=getTechVersusVesselDis(ambr_taxa = ambr_taxa,ambr_metadata = ambr_metadata, reftechrep = "r3")

# Make a violin plot for dissimilarities
na.num=length(res1$tech.dissims)-length(res1$bio.dissims)
# technical dissimilarities are identical across the 3 runs, so only 1 is considered
df=cbind(c(res1$bio.dissims,rep(NA,na.num)), c(res2$bio.dissims,rep(NA,na.num)), c(res3$bio.dissims,rep(NA,na.num)),res1$tech.dissims)
unit1="Var.T1"
unit2="Var.T2"
unit3="Var.T3"
unit4="Technical.Variability"
colnames(df)=c(unit1,unit2,unit3,unit4)
df.melted=melt(df)
combinations=list()
combinations[[paste(unit1,unit2,sep="")]]=c(unit1,unit2)
combinations[[paste(unit1,unit3,sep="")]]=c(unit1,unit3)
combinations[[paste(unit1,unit4,sep="")]]=c(unit1,unit4)
combinations[[paste(unit2,unit3,sep="")]]=c(unit2,unit3)
combinations[[paste(unit2,unit4,sep="")]]=c(unit2,unit4)
ggplot(df.melted, aes(Var2, value)) +geom_violin()+ ggpubr::stat_compare_means(comparisons=combinations, method="wilcox.test", p.adjust.method = "bh") + xlab("Relative abundance")+ylab("Bray Curtis dissimilarity")


# Make a violin plot for coefficient of variation for relative abundances (standard deviation divided by mean)
# 24 values for each species and time point, computed across the six vessels
reorder=c(19:24,1:6,7:12,13:18)  # reorder species and skip unknowns
cv.cs.values=cs.data[reorder,4]/cs.data[reorder,3] # standard deviation divided by mean
cv.16S.values=seq.data[,4]/seq.data[,3] # selected technical replicate
wilcox.test(cv.cs.values,cv.16S.values, paired=TRUE) # paired Wilcoxon p-value is highly significant
df=cbind(cv.16S.values,cv.cs.values)
unit1="16S"
unit2="CellScanner"
colnames(df)=c(unit1,unit2)
df.melted=melt(df)
combinations=list()
combinations[[paste(unit1,unit2,sep="")]]=c(unit1,unit2) # a single combination is tested, multiple-testing correction not required
ggplot(df.melted, aes(Var2, value)) +geom_violin()+ ggpubr::stat_compare_means(comparisons=combinations, method="wilcox.test", p.adjust.method = "none") + xlab("")+ylab("Coefficient of variation")

# same for metabolites
cv.met=met.data[,4]/met.data[,3]
wilcox.test(cv.met,cv.16S.values, alternative="less") # cannot be paired

# Figure 6 - manual gating only
counts.exp1=plotSingleRepTS(abundances=t(as.matrix(counts[,3])),groups=counts$Vessel, return=TRUE,time=counts$Timepoint,title="Cell count per mL experiment 1")
counts.exp2=plotSingleRepTS(abundances=t(as.matrix(ambr2.totalcounts)),groups=ambr_metadata2$Vessel, return=TRUE,time=ambr_metadata2$Time,title="Cell count per mL experiment 2")
counts.exp1[["Experiment"]]=rep("exp1",nrow(counts.exp1))
counts.exp2[["Experiment"]]=rep("exp2",nrow(counts.exp2))
counts.ambr=rbind(counts.exp1, counts.exp2)
ggplot(counts.ambr, aes(x=Time, y=Mean, color=Experiment)) + ggtitle("Total cell counts") + xlab("Time in hours") +  ylab("Cell count per mL") + geom_line() + theme_classic() + geom_pointrange(aes(ymin=Mean-sd, ymax=Mean+sd))

h.vec=as.vector(heterogeneities) # appended column-wise, so all vessels for 1 time point
plotSingleRepTS(abundances=t(as.matrix(h.vec)),groups=counts$Vessel,time=counts$Timepoint,title="Heterogeneity",ylab="Sum of mean range of events across channels", showPval = FALSE)

# Beta diversity (supplement)
# get vessels and technical replicate of samples
tech.reps=c()
vessels=c()
for(i in 1:ncol(ambr_taxa)){
  sample.index=which(ambr_metadata$X==colnames(ambr_taxa)[i])
  tech.reps=c(tech.reps,ambr_metadata$Replicate[sample.index])
  vessels=c(vessels,ambr_metadata$Vessel[sample.index])
}
compareGroups(abundances=ambr_taxa, groups = tech.reps, property="beta", pvalViz = TRUE)
compareGroups(abundances=ambr_taxa, groups = vessels, property="beta", pvalViz = TRUE)

########## Statistics ##########

total.counts.cs2=colSums(cs2)

# Table 1
# Agreement CS and 16S
ri.cor=cor.test(cs.proportions[4,],ambr_taxa[1,selected.indices]) # Roseburia
bh.cor=cor.test(cs.proportions[1,],ambr_taxa[2,selected.indices]) # Blautia
bt.cor=cor.test(cs.proportions[2,],ambr_taxa[3,selected.indices]) # Bacteroides
ca.cor=cor.test(cs.proportions[3,],ambr_taxa[4,selected.indices]) # Collinsella
print(paste("Corrected p-value threshold (Bonferroni): ",0.05/4)) # Bonferroni (division of threshold p-value by number of comparisons)
print(paste("Roseburia:", round(ri.cor$estimate,3), ri.cor$p.value))
print(paste("Blautia:", round(bh.cor$estimate,3), bh.cor$p.value))
print(paste("Bacteroides:",round(bt.cor$estimate,3), bt.cor$p.value))
print(paste("Collinsella:", round(ca.cor$estimate,3),ca.cor$p.value))

# Absolute counts vs time
# main experiment
cor.test(as.numeric(time[selected.indices]),colSums(ambr_taxa_abs))
# r: -0.76, p-val < 0.00001

# validation experiment
cor.test(ambr_metadata2$Time,ambr2.totalcounts)
# -0.02, p-val = 0.862

# BT absolute counts vs time for first experiment
t12=c(1,7,13,19,25,31)
t25=c(2,8,14,20,26,32)
# absolute abundances BT at 12h and 25h
bt.12h=ambr_taxa_abs[3,t12]
bt.25h=ambr_taxa_abs[3,t25]
wilcox.test(bt.12h,bt.25h)
# p-value = 0.002

# BT absolute counts vs time for second experiment
t36=c(4,26,48)
t60=c(6,28,50)
# absolute abundances BT at 36h and 60h
bt.36h=ambr_taxa2_abs[3,t36]
bt.60h=ambr_taxa2_abs[3,t60]
wilcox.test(bt.36h,bt.60h)
# p-value = 0.1 (lower number of vessels)


# total counts first experiment
counts.12h=colSums(ambr_taxa_abs)[t12]
counts.25h=colSums(ambr_taxa_abs)[t25]
wilcox.test(counts.12h,counts.25h)
# W = 36, p-value = 0.002165

# total counts second experiment
counts.36h=total.counts.cs2[t36]
counts.60h=total.counts.cs2[t60]
wilcox.test(counts.36h,counts.60h) # not significant
boxplot(counts.36h,counts.60h,names=c("36h","60h"))

# heterogeneities 25h vs 67h
wilcox.test(heterogeneities[,2],heterogeneities[,6])
# W = 7, p-value = 0.09307

########## Correlation networks (Supplementary Figure) ########

net.16S=barebonesCoNet(abundances=ambr_taxa_abs[1:4,], metadata=metabolites,methods=c("spearman"), min.occ = 0,init.edge.num="all", pval.cor =TRUE)
# createNetworkFromIgraph(net.16S) # for annotation in Cytoscape
plot(net.16S, main="16S absolute abundance Spearman network")
net.cs=barebonesCoNet(abundances=cs, metadata=metabolites,methods=c("spearman"),min.occ=0, init.edge.num = "all", pval.cor =TRUE)
plot(net.cs, main="CS absolute abbundance Spearman network")

####### Plots for 2nd experiment ###########

group.colors <- c("Blautia" = "#ED7D31", "Bacteroides" = "#A5A5A5", "Roseburia" ="#4472C3", "Faecalibacterium" = "#FFBF01")
group.colors.cs=group.colors
group.colors.cs[["Unknown"]]="black"

# plot absolute species abundances over time
title="16S counts exp2"
ylab="Cells/ml"
seq2.data.abs = plotGroupRepTS(ambr_taxa2_abs,groups=ambr_metadata2$Vessel,time=ambr_metadata2$Time, ylab=ylab, title=title, returnData = TRUE, colors=group.colors)
p2=ggplot(seq2.data.abs, aes(x=Time, y=Mean, color=Taxon)) + geom_line() + theme_classic() + geom_pointrange(aes(ymin=Mean-sd, ymax=Mean+sd)) + ggtitle(title) + ylab(ylab) + xlab(xlab) + scale_colour_manual(values=group.colors)
plot(p2)

# plot relative species abundances over time
title="16S exp2"
ylab="Relative abundance"
seq2.data = plotGroupRepTS(ambr_taxa2,groups=ambr_metadata2$Vessel,time=ambr_metadata2$Time, ylab=ylab, title=title, returnData = TRUE, colors=group.colors)
p2=ggplot(seq2.data, aes(x=Time, y=Mean, color=Taxon)) + geom_line() + theme_classic() + geom_pointrange(aes(ymin=Mean-sd, ymax=Mean+sd)) + ggtitle(title) + ylab(ylab) + xlab(xlab) + scale_colour_manual(values=group.colors)
plot(p2)

# plot relative abundances with CS
# WARN: since no monoculture data were collected for exp2, CS data are less reliable here
rownames(cs2.proportions)=c("Blautia","Bacteroides","Faecalibacterium","Roseburia","Unknown")
title="CS exp2"
ylab="Relative abundance"
cs2.data = plotGroupRepTS(cs2.proportions,groups=cs2.metadata$Vessel,time=cs2.metadata$Time, ylab=ylab, title=title, returnData = TRUE, colors=group.colors.cs)
p2=ggplot(cs2.data, aes(x=Time, y=Mean, color=Taxon)) + geom_line() + theme_classic() + geom_pointrange(aes(ymin=Mean-sd, ymax=Mean+sd)) + ggtitle(title) + ylab(ylab) + xlab(xlab) + scale_colour_manual(values=group.colors.cs)
plot(p2)

# Comparison of coefficient of variation per time point across biological replicates for 16S experiments 1 and 2 (with selected technical replicate for exp1)
cv.16S2.values=seq2.data[,4]/seq2.data[,3]
df2=cbind(cv.16S.values,cv.16S2.values)
unit1="16Sexp1"
unit2="16Sexp2"
colnames(df2)=c(unit1,unit2)
df2.melted=melt(df2)
combinations=list()
combinations[[paste(unit1,unit2,sep="")]]=c(unit1,unit2) # a single combination is tested, multiple-testing correction not required
ggplot(df2.melted, aes(Var2, value)) +geom_violin()+ ggpubr::stat_compare_means(comparisons=combinations, method="wilcox.test", p.adjust.method = "none") + xlab("Relative abundances")+ylab("Coefficient of variation across vessels per time point")

# ordination plot with time series (relative abundances) for exp2; time is already given in ascending order per group
# vessel trajectories are much closer to each other initially than in exp1
seqPCoA(ambr_taxa2,groups=ambr_metadata2$Vessel, topTaxa=nrow(ambr_taxa2), main="PCoA exp2 for 16S", labels = ambr_metadata2$Time, arrowFactor = arrowFactor, time=TRUE,  xlim=c(-0.41,xmax), ylim=c(-0.4,0.4), hideGroupLegend=TRUE)

# Interpretation for initial similarity: Blautia and then Bacteroides dominates first time points to a large extent

# ordination plot without first seven time points
selected.indices=selectDataSubSet(data=ambr_taxa2,metadata=ambr_metadata2,metadataItem = "Time", selectedMetadataValues = ambr_metadata2$Time[8:22])
seqPCoA(ambr_taxa2[,selected.indices],groups=ambr_metadata2$Vessel[selected.indices], main="PCoA exp2 for 16S without first 7 time points", topTaxa=nrow(ambr_taxa2), labels = ambr_metadata2$Time[selected.indices], arrowFactor = arrowFactor, time=TRUE,  xlim=c(-0.35,0.35), ylim=c(-0.4,0.4), hideGroupLegend=TRUE)

# beta diversity per vessel for exp2
compareGroups(abundances=ambr_taxa2, groups = ambr_metadata2$Vessel, property="beta", pvalViz = TRUE, xlab = "Vessels exp2 16S (relative abundance)")

# Look at CellScanner data for exp2
sorted.indices.cs=sortSamples(cs2,groups = cs2.metadata$Vessel, time = cs2.metadata$Time)

seqPCoA(cs2.proportions[,sorted.indices.cs],groups=cs2.metadata$Vessel[sorted.indices.cs], main="PCoA exp2 for CS", topTaxa=nrow(cs2.proportions), labels = cs2.metadata$Time[sorted.indices.cs], arrowFactor = arrowFactor, time=TRUE,  xlim=c(-0.45,xmax), ylim=c(-0.4,0.4), hideGroupLegend=TRUE, clusterQualityIndex="none")

# ordination plot without first seven time points
selected.indices.cs=selectDataSubSet(data=cs2.proportions[,sorted.indices.cs],metadata=cs2.metadata[sorted.indices.cs,],metadataItem = "Time", selectedMetadataValues = c(84,96,108,120,132,144,156,168,180,192,204,216,228,240,252))
seqPCoA(cs2.proportions[,sorted.indices.cs[selected.indices.cs]],groups=cs2.metadata$Vessel[sorted.indices.cs[selected.indices.cs]], main="PCoA exp2 for CS without first 7 time points", topTaxa=nrow(cs2.proportions), labels = cs2.metadata$Time[sorted.indices.cs[selected.indices.cs]], arrowFactor = arrowFactor, time=TRUE,  xlim=c(-0.5,xmax), ylim=c(-0.4,0.41), hideGroupLegend=TRUE, clusterQualityIndex="none")

reorder2=c(4,1,2,3) # bring CS species in the same order as in 16S exp2 data and exclude unknowns
cs2.data = plotGroupRepTS(cs2.proportions[reorder2,],groups=cs2.metadata$Vessel,time=cs2.metadata$Time, ylab=ylab, title=title, returnData = TRUE)
cv.cs2.values=cs2.data[,4]/cs2.data[,3] # standard deviation divided by mean
wilcox.test(cv.cs2.values,cv.16S2.values, paired=TRUE) # paired Wilcoxon p-value is significant
df3=cbind(cv.16S2.values,cv.cs2.values)
unit1="16Sexp2"
unit2="CSexp2"
colnames(df3)=c(unit1,unit2)
df3.melted=melt(df3)
combinations=list()
combinations[[paste(unit1,unit2,sep="")]]=c(unit1,unit2) # a single combination is tested, multiple-testing correction not required
ggplot(df3.melted, aes(Var2, value)) +geom_violin()+ ggpubr::stat_compare_means(comparisons=combinations, method="wilcox.test", p.adjust.method = "none") + xlab("Relative abundances")+ylab("Coefficient of variation across vessels per time point") # wilcox test here is not paired
# shift is also significant when including unknowns



