rm(list=ls())

########################## set working directory ####################################################################

setwd('/home/gang/Documents/R data analysis Treg/bead.solution/SnAct')

########### functions ##############################################################

# function to split majority protein IDs, store in a list
mprotein.ID.list=function(mprotein.ID)
{
  ans=list()
  for (i in 1:length(mprotein.ID))
  {ans[[i]]=unlist(strsplit(mprotein.ID[i],';'))}
  return(ans)
}

# function to get the alignment.index for majority ID alignment
alignment.index=function(mproteinID.1, mproteinID.2)
{
  temp.1=mprotein.ID.list(mproteinID.1)
  temp.2=mprotein.ID.list(mproteinID.2)
  ans=matrix(nc=2,nr=length(mproteinID.1))
  colnames(ans)=c('index.1','index.2')
  for (i in 1:length(mproteinID.1))
  {
    for (j in 1:length(mproteinID.2))
    {
      if (length(intersect(temp.1[[i]],temp.2[[j]]))>0)
      {
        ans[i,]=c(i,j)
      }
    }
  }
  return(ans)
}

# find N-glycosylation annotation in uniprot, returns a logical value
uniprot.Ngly=function(uniprot)
{
  ans=c()
  for (i in 1:nrow(uniprot))
  {
    temp=unlist(strsplit(as.vector(uniprot$Glycosylation[i]),' '))
    if ('N-linked' %in% temp)
    {
      ans=c(ans,TRUE)
    }
    else {ans=c(ans,FALSE)}
  }
  return(ans)
}


imputation=function(x=vector.for.input, down.shift, width) # funciton for imputation of a vector
{
  a=x[!(is.na(x))] #obtain the valid data
  b=rnorm(n=length(x[is.na(x)]),mean=mean(a)-down.shift*sd(a),sd=width*sd(a)) # produce the random number
  x[is.na(x)]=b # replace NA values with the random numbers
  return(x)
}


# function to do paired t test in the counter.receptor.LFQ data.frame, return p value
ttest.p=function(counter.receptor.LFQ)
{
  ans=c()
  for (i in 1:nrow(counter.receptor.LFQ))
  {
    p=t.test(as.numeric(counter.receptor.LFQ[i,c('bead.R97A.1','bead.R97A.3',
                                      'solution.R97A.1','solution.R97A.2','solution.R97A.3')]),
             as.numeric(counter.receptor.LFQ[i,c('bead.Sn.1','bead.Sn.3',
                                      'solution.Sn.1','solution.Sn.2','solution.Sn.3')]),
             paired=TRUE)
    
    ans=c(ans,p[[3]])
  }
  return(ans)
}

# function to do paired t test in the counter.receptor.LFQ data.frame, return the difference of the means
ttest.df=function(counter.receptor.LFQ)
{
  ans=c()
  for (i in 1:nrow(counter.receptor.LFQ))
  {
    p=t.test(as.numeric(counter.receptor.LFQ[i,c('bead.R97A.1','bead.R97A.3',
                                                 'solution.R97A.1','solution.R97A.2','solution.R97A.3')]),
             as.numeric(counter.receptor.LFQ[i,c('bead.Sn.1','bead.Sn.3',
                                                 'solution.Sn.1','solution.Sn.2','solution.Sn.3')]),
             paired=TRUE)
    
    ans=c(ans,-p[[5]])
  }
  return(ans)
}

########### log2 LFQ intensities for bead proximity labeling #########################################################  

# load bead proximity labeling raw data
raw.data=read.delim('/home/gang/Documents/R data analysis Treg/Raw data/Bead proteinGroups.txt')
raw.data=subset(raw.data,raw.data$Only.identified.by.site!='+')
raw.data=subset(raw.data,raw.data$Reverse!='+')
raw.data=subset(raw.data,raw.data$Potential.contaminant!='+')

if (length(levels(factor(raw.data$Majority.protein.IDs)))!=nrow(raw.data))
{stop('abnormal mproteinID found')}
if (NA %in% unlist(strsplit(as.vector(raw.data$Majority.protein.IDs),';')))
{stop('NA found in mproteinID')}

# obtain SnAct LFQ rows
LFQ=subset.data.frame(raw.data, select = c(LFQ.intensity.E1,LFQ.intensity.E2,LFQ.intensity.E3,
                                           LFQ.intensity.F1,LFQ.intensity.F2,LFQ.intensity.F3))
row.names(LFQ)=as.vector(raw.data$Majority.protein.IDs)

LFQ=log2(LFQ)

# remove -Inf (log2(0)) values of log2(LFQ) from positive experiments
LFQ=subset.data.frame(LFQ,rowSums(LFQ[4:6])!=-Inf)

# imputation separately for each column
LFQ[LFQ==-Inf]=NA # replace -Inf with NA

# imputation
LFQ[1]=imputation(LFQ[1],down.shift=1.8, width=0.3) # imputation for neg.1
LFQ[2]=imputation(LFQ[2],down.shift=1.8, width=0.3) # imputation for neg.2
LFQ[3]=imputation(LFQ[3],down.shift=1.8, width=0.3) # imputation for neg.3

bead.LFQ=LFQ

########### log2 LFQ intensities for solution proximity labeling #########################################################  

# load solution proximity labeling raw data
raw.data=read.delim('/home/gang/Documents/R data analysis Treg/Raw data/Solution proteinGroups.txt')
raw.data=subset(raw.data,raw.data$Only.identified.by.site!='+')
raw.data=subset(raw.data,raw.data$Reverse!='+')
raw.data=subset(raw.data,raw.data$Potential.contaminant!='+')

if (length(levels(factor(raw.data$Majority.protein.IDs)))!=nrow(raw.data))
{stop('abnormal mproteinID found')}
if (NA %in% unlist(strsplit(as.vector(raw.data$Majority.protein.IDs),';')))
{stop('NA found in mproteinID')}

# obtain SnAct LFQ rows
LFQ=subset.data.frame(raw.data, select = c(LFQ.intensity.E1,LFQ.intensity.E2,LFQ.intensity.E3,
                                           LFQ.intensity.F1,LFQ.intensity.F2,LFQ.intensity.F3))
row.names(LFQ)=as.vector(raw.data$Majority.protein.IDs)

LFQ=log2(LFQ)

# remove -Inf (log2(0)) values of log2(LFQ) from positive experiments
LFQ=subset.data.frame(LFQ,rowSums(LFQ[4:6])!=-Inf)

# imputation separately for each column
LFQ[LFQ==-Inf]=NA # replace -Inf with NA

# imputation
LFQ[1]=imputation(LFQ[1],down.shift=1.8, width=0.3) # imputation for neg.1
LFQ[2]=imputation(LFQ[2],down.shift=1.8, width=0.3) # imputation for neg.2
LFQ[3]=imputation(LFQ[3],down.shift=1.8, width=0.3) # imputation for neg.3

solution.LFQ=LFQ

########### combine the bead and solution proxmity labeling data #############################

# match the row names of bead.LFQ and solution.LFQ
alignment = alignment.index(row.names(bead.LFQ),row.names(solution.LFQ))
alignment=alignment[!is.na(rowSums(alignment)),]

bead.LFQ=bead.LFQ[alignment[,1],]
solution.LFQ=solution.LFQ[alignment[,2],]

# create a new data.frame, use row names of solution.LFQ as row names
bead.solution.LFQ=data.frame(row.names = row.names(solution.LFQ))

# Add gene names
row.names(raw.data)=as.vector(raw.data$Majority.protein.IDs)
bead.solution.LFQ['Gene.names']=raw.data[row.names(bead.solution.LFQ),'Gene.names']

bead.solution.LFQ[c('bead.R97A.1','bead.R97A.2','bead.R97A.3',
                    'bead.Sn.1','bead.Sn.2','bead.Sn.3')]=bead.LFQ

bead.solution.LFQ[c('solution.R97A.1','solution.R97A.2','solution.R97A.3',
                    'solution.Sn.1','solution.Sn.2','solution.Sn.3')]=solution.LFQ

########### paired ttest ############################################

LFQ.diff =ttest.df(bead.solution.LFQ)
LFQ.diff=as.vector(LFQ.diff)
bead.solution.LFQ['LFQ.diff']=LFQ.diff

LFQ.p =ttest.p(bead.solution.LFQ)
bead.solution.LFQ['LFQ.p']=LFQ.p

########### annotate counter receptor rows ######################################

# add a column with all rows annotated as non counter receptor
bead.solution.LFQ['counter.receptor.identification']='non counter recetor'

# p<=0.05, diff>1
temp = subset.data.frame(bead.solution.LFQ, LFQ.diff>1 & LFQ.p<=0.05)
bead.solution.LFQ[row.names(temp),'counter.receptor.identification']= 'non glycosylated counter receptor'

########### annotate glycosylated counter receptor rows ##################################

# load Uniprot lib.txt
uniprot=read.delim('/home/gang/Documents/R data analysis Treg/Raw data/Uniprot lib.txt')

# obtain subset of uniprot with N-glycosylation annotation
uniprot=subset.data.frame(uniprot,Glycosylation!='')
uniprot=subset.data.frame(uniprot,uniprot.Ngly(uniprot))

# annotate glycosylated counter receptor
temp=subset.data.frame(bead.solution.LFQ, counter.receptor.identification=='non glycosylated counter receptor')
for (i in 1:nrow(temp))
{
  if (length(intersect(unlist(strsplit(row.names(temp)[i],';')),as.vector(uniprot$Entry))) >0)
  {bead.solution.LFQ[row.names(temp)[i],"counter.receptor.identification"]='glycosylated counter receptor'}
}

########## use the boders set on proximity labeling tail for further filtration ################

# add a column with all row annoated with FALSE
bead.solution.LFQ['within.proximity.labeling.tail']=FALSE

# glycosylated counter receptor subset
temp=subset.data.frame(bead.solution.LFQ, counter.receptor.identification=='glycosylated counter receptor')

# border for bead pair 1, set at 0.5
temp=subset.data.frame(temp, bead.Sn.1-bead.R97A.1>0.5)

# border for bead pair 3, set at 0.75
temp=subset.data.frame(temp,bead.Sn.3-bead.R97A.3>0.75)

# border for solution pair 1, set at 1
temp=subset.data.frame(temp,solution.Sn.1-solution.R97A.1>1)

# border for solution pair 2, set at 1
temp=subset.data.frame(temp,solution.Sn.2-solution.R97A.2>1)

# border for solution pair 3, set at 1
temp=subset.data.frame(temp,solution.Sn.3-solution.R97A.3>1)

# give true value for the glycosylated counter receptors within the proximity labeling tail
bead.solution.LFQ[row.names(temp),"within.proximity.labeling.tail"]=TRUE
