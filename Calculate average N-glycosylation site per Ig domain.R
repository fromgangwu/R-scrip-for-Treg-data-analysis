rm(list=ls())

###################### set working directory ##################################################################

setwd('/home/gang/Documents/R data analysis Treg/bead.solution/SnAct')

###################### functions ##############################################################################

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

# find Ig domain annotation in uniprot, returns a logical value
uniprot.Ig=function(uniprot)
{
  ans=c()
  for (i in 1:nrow(uniprot))
  {
    temp=unlist(strsplit(as.vector(uniprot$Domain..FT.[i]),' '))
    if ('Ig-like' %in% temp | 'Ig-like.' %in% temp)
    {
      ans=c(ans,TRUE)
    }
    else {ans=c(ans,FALSE)}
  }
  return(ans)
}

# function to find number of Ig domains of a uniprot Entry
Ig.domain.count=function(Entry)
{
  temp=as.vector(uniprot[uniprot$Entry==Entry,'Domain..FT.'])
  temp=unlist(strsplit(temp,' '))
  if ('Ig-like' %in% temp | 'Ig-like.' %in% temp)
  {
    count=0
    for (i in 1:length(temp))
    {
      if ('Ig-like' == temp[i] | 'Ig-like.' == temp[i])
      {count=count+1}
    }
    ans=count
  }
  else {ans=0}
  return(as.numeric(ans))
}

# function to find Ig domain site of a uniprot Entry
Ig.domain.site=function(Entry)
{
  ans=c()
  temp=as.vector(uniprot[uniprot$Entry==Entry,'Domain..FT.'])
  temp=unlist(strsplit(temp,' '))
  if ('Ig-like' %in% temp | 'Ig-like.' %in% temp)
  {
    for (i in 1:length(temp))
    {
      if ('Ig-like' == temp[i] | 'Ig-like.' == temp[i])
      {ans=c(ans,temp[i-2]:temp[i-1])}
    }
  }
  else{ans=0}
  return(as.numeric(ans))
}

# function to find N-glycosylation site of a uniprot Entry
Ngly.site=function(Entry)
{
  ans=c()
  temp=as.vector(uniprot[uniprot$Entry==Entry,'Glycosylation'])
  temp=unlist(strsplit(temp,' '))
  if ('N-linked' %in% temp)
  {
    for (i in 1:length(temp))
    {
      if ('N-linked' == temp[i])
      {ans=c(ans,temp[i-1])}
    }
  }
  else{ans=0}
  return(as.numeric(ans))
}

# function to calculate N-glycans per Ig domain
Ngly.per.Ig=function(Entry)
{
  Nglysites=Ngly.site(Entry)
  Igsites=Ig.domain.site(Entry)
  Ngly.on.Ig=intersect(Nglysites,Igsites)
  Igdomains=Ig.domain.count(Entry)
  ans=length(Ngly.on.Ig)/Igdomains
  return(ans)
}


###################### obtain glycosylated Ig domain data from uniprot lib  ##################################################################

# load uniprot lib
uniprot=read.delim('/home/gang/Documents/R data analysis Treg/Raw data/Uniprot lib.txt')

# subset reviewed entry
uniprot=subset.data.frame(uniprot, Status=='reviewed')

# subset glycosylated entry
uniprot=subset.data.frame(uniprot, Glycosylation!='')
temp=uniprot.Ngly(uniprot)
uniprot=subset.data.frame(uniprot,temp)

# subset glycosylated Ig domain entry 
uniprot=subset.data.frame(uniprot,uniprot$Domain..FT.!='')
temp=uniprot.Ig(uniprot)
uniprot=subset.data.frame(uniprot, temp)


###################### intersect with proteomics data from Treg act ###########################################

# load Total glycoproteins Treg act
glycoproteins=read.csv('/home/gang/Documents/R data analysis Treg/bead.solution/Total glycoproteins Treg act',
                       row.names = 1)

# obtain majority protein IDs of glycoproteins
glycoproteins=unlist(strsplit(as.vector(glycoproteins[,1]),';'))

# intersect with uniprot entry
glycoproteins=intersect(glycoproteins, as.vector(uniprot$Entry))

# subset uniprot data based according to the intersect results
row.names(uniprot)=as.vector(uniprot$Entry)
uniprot=uniprot[glycoproteins,]

################## annotate N-glycan per Ig domain ##########################################################
temp=vector()
for (i in 1:nrow(uniprot))
{
  temp[i]=Ngly.per.Ig(as.vector(uniprot$Entry[i]))
}

uniprot['Nglycan.per.Ig']=temp

################## counter receptor annotation #############################################################

# load counter receptor data 
counter.receptor = read.csv('Counter receptors',row.names = 1)

# subset glycosylated counter receptors within proximity labeling tail
counter.receptor = subset.data.frame(counter.receptor, within.proximity.labeling.tail==TRUE)

# obtain counter receptor uniprot entry names
counter.receptor = row.names(counter.receptor)
counter.receptor = unlist(strsplit(counter.receptor,';'))

# subset counter receptor entry names which are in uniprot N-glycan.per.Ig data
counter.receptor = intersect(counter.receptor,as.vector(uniprot$Entry))

# counter receptor annotation
uniprot['counter.receptor']=FALSE
uniprot[counter.receptor,'counter.receptor']=TRUE

################# remove IGHG1, Sn, SigE and CD22 ##########################################################

# CD22
uniprot=subset.data.frame(uniprot,Entry!='P35329')

# IGHG1
uniprot=subset.data.frame(uniprot,Entry!='P01857')

# SigE
uniprot=subset.data.frame(uniprot,Entry!='Q91Y57')

# Sn
uniprot=subset.data.frame(uniprot, Entry!='Q62230')

########## annotate copy number per cell  ###########################################################

# load data, choose Treg act column log2
raw.data=read.csv('/home/gang/Documents/R data analysis Treg/Raw data/Copy number per cell.csv')

copy.per.cell=raw.data[5:7]
row.names(copy.per.cell)=as.vector(raw.data$Majority.protein.IDs)

copy.per.cell=log2(copy.per.cell)

# keep rows without 0 values
temp=copy.per.cell!=-Inf
temp=rowSums(temp)
temp=temp>=3

copy.per.cell=subset.data.frame(copy.per.cell,temp)

# add row means of copy number per cell
temp=rowMeans(copy.per.cell)
copy.per.cell['mean']=temp

# add row sd of copy number per cell
temp=vector()
for (i in 1:nrow(copy.per.cell))
{
  temp[i]=sd(as.vector(copy.per.cell[i,1:3]))
}
copy.per.cell['standard.deviation']=temp

# annotate uniprot data with mean and sd of copy number per cell
uniprot['copy.number.per.cell.mean']=NA
uniprot['copy.number.per.cell.sd']=NA

temp=unlist(strsplit(row.names(copy.per.cell),';'))
for (i in 1:nrow(uniprot))
{
  if (row.names(uniprot)[i] %in% temp)
  {
    for (j in 1:nrow(copy.per.cell))
      if (row.names(uniprot)[i] %in% unlist(strsplit(row.names(copy.per.cell)[j],';')))
      {
        uniprot[i,'copy.number.per.cell.mean']=as.numeric(copy.per.cell$mean[j])
        uniprot[i,'copy.number.per.cell.sd']=as.numeric(copy.per.cell$standard.deviation[j])
      }
  }
}

######### draw graphs ########################################################################

library(ggplot2)
# histogram counter receptors
original=ggplot(data=subset.data.frame(uniprot, counter.receptor), aes(Nglycan.per.Ig))

original + geom_histogram(fill='grey',color='black', breaks=seq(0,5,by=1)) + 
           scale_x_continuous(breaks = seq(-6,10,by=1)) +
           theme(text = element_text(size=10), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5, size=10))

ggsave('Ngly.per.Ig.counter.receptor', device = 'tiff', dpi='retina', width = 5, height = 3, units = 'in')

# histogram non counter receptors
original=ggplot(data=subset.data.frame(uniprot, !counter.receptor), aes(Nglycan.per.Ig))

original + geom_histogram(fill='grey',color='black', breaks=seq(0,5,by=1)) + 
  scale_x_continuous(breaks = seq(-6,10,by=1)) +
  theme(text = element_text(size=10), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5, size=10))

ggsave('Ngly.per.Ig.non.counter.receptor', device = 'tiff', dpi='retina', width = 5, height = 3, units = 'in')

# boxplot 
temp=uniprot
temp$counter.receptor[temp$counter.receptor==TRUE]='counter receptor'
temp$counter.receptor[temp$counter.receptor==FALSE]='non counter receptor'

original=ggplot(temp, aes(x=counter.receptor, y=Nglycan.per.Ig))

original + geom_boxplot() + 
  theme(text = element_text(size=10), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5, size=10)) + 
  labs(x='',y='N-glycan per Ig')+annotate('text',x=1.5, y=4.75, label='p=0.01')

ggsave('Ngly.per.Ig.box.plot', device = 'tiff', dpi='retina', width = 4, height = 4, units = 'in')


