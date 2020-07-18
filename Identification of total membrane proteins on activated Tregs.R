rm(list=ls())

###################### functions ##############################################################################

# find plasma membrane annotation in uniprot, returns a logical value
uniprot.pm=function(uniprot)
{
  ans=vector()
  for (i in 1:nrow(uniprot))
  {
    ans[i]=FALSE
    temp=unlist(strsplit(as.vector(uniprot$Gene.ontology..cellular.component.[i]),' '))
    for (j in 1:length(temp))
    {
      if (temp[j]=='plasma' & temp[j+1]=='membrane')
      {
        ans[i]=TRUE
        break
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


###################### obtain membrane protein data from uniprot lib  ##################################################################

# load uniprot lib
uniprot=read.delim('C:/Users/Gang Wu/Dropbox/R data analysis of Treg/R data analysis Treg/Raw data/Uniprot lib.txt')

# subset reviewed entry
uniprot=subset.data.frame(uniprot, Status=='reviewed')

# subset menbrane protein entry
uniprot=subset.data.frame(uniprot, Gene.ontology..cellular.component.!='')
temp=uniprot.pm(uniprot)
uniprot=subset.data.frame(uniprot,temp)

################# remove HRP IGHG1, Sn, SigE and CD22 ##########################################################

# HRP
uniprot=subset.data.frame(uniprot,Entry!='P00433')

# CD22
uniprot=subset.data.frame(uniprot,Entry!='P35329')

# IGHG1
uniprot=subset.data.frame(uniprot,Entry!='P01857')

# SigE
uniprot=subset.data.frame(uniprot,Entry!='Q91Y57')

# Sn
uniprot=subset.data.frame(uniprot, Entry!='Q62230')

################# membrane proteins on bead proximity labeling ############################################

# load the data
raw.data=read.delim('C:/Users/Gang Wu/Dropbox/R data analysis of Treg/R data analysis Treg/Raw data/Bead proteinGroups.txt')

# obtain the intensity colomns for Treg Act
raw.data=raw.data[c("Majority.protein.IDs","Intensity.F1", "Intensity.F2", "Intensity.F3",
                                           "Intensity.G1", "Intensity.G2", "Intensity.G3",
                                           "Intensity.H1", "Intensity.H2", "Intensity.H3")]

# subset rows with at least 2 valid values for Intensity.F
temp=raw.data[c("Intensity.F1", "Intensity.F2", "Intensity.F3")]
temp=temp!=0
temp=rowSums(temp)
temp=temp>=2
raw.data=subset.data.frame(raw.data,temp)

# subset rows with at least 2 valid values for Intensity.G
temp=raw.data[c("Intensity.G1", "Intensity.G2", "Intensity.G3")]
temp=temp!=0
temp=rowSums(temp)
temp=temp>=2
raw.data=subset.data.frame(raw.data,temp)

# subset rows with at least 2 valid values for Intensity.H
temp=raw.data[c("Intensity.H1", "Intensity.H2", "Intensity.H3")]
temp=temp!=0
temp=rowSums(temp)
temp=temp>=2
raw.data=subset.data.frame(raw.data,temp)

# obtain the membrane majority protein IDs
mproteinID=as.vector(raw.data$Majority.protein.IDs)

temp=c()
for (i in 1:length(mproteinID))
{
  if (length(intersect(unlist(strsplit(mproteinID[i],';')),as.vector(uniprot$Entry)))>0)
  {
    temp=c(temp,mproteinID[i])
  }
}

bead=temp

################# membrane proteins in solution proximity labeling ############################################

# load the data
raw.data=read.delim('C:/Users/Gang Wu/Dropbox/R data analysis of Treg/R data analysis Treg/Raw data/Solution proteinGroups.txt')

# obtain the intensity colomns for Treg Act
raw.data=raw.data[c("Majority.protein.IDs","Intensity.F1", "Intensity.F2", "Intensity.F3",
                                           "Intensity.G1", "Intensity.G2", "Intensity.G3",
                                           "Intensity.H1", "Intensity.H2", "Intensity.H3")]

# subset rows with at least 2 valid values for Intensity.F
temp=raw.data[c("Intensity.F1", "Intensity.F2", "Intensity.F3")]
temp=temp!=0
temp=rowSums(temp)
temp=temp>=2
raw.data=subset.data.frame(raw.data,temp)

# subset rows with at least 2 valid values for Intensity.G
temp=raw.data[c("Intensity.G1", "Intensity.G2", "Intensity.G3")]
temp=temp!=0
temp=rowSums(temp)
temp=temp>=2
raw.data=subset.data.frame(raw.data,temp)

# subset rows with at least 2 valid values for Intensity.H
temp=raw.data[c("Intensity.H1", "Intensity.H2", "Intensity.H3")]
temp=temp!=0
temp=rowSums(temp)
temp=temp>=2
raw.data=subset.data.frame(raw.data,temp)

# obtain the membrane majority protein IDs
mproteinID=as.vector(raw.data$Majority.protein.IDs)

temp=c()
for (i in 1:length(mproteinID))
{
  if (length(intersect(unlist(strsplit(mproteinID[i],';')),as.vector(uniprot$Entry)))>0)
  {
    temp=c(temp,mproteinID[i])
  }
}

solution=temp

################# membrane proteins in total membrane protein #############################################

# load the data
raw.data=read.delim('C:/Users/Gang Wu/Dropbox/R data analysis of Treg/R data analysis Treg/Raw data/Membrane proteinGroups.txt')

# obtain the intensity colomns for Treg rest
raw.data=raw.data[c("Majority.protein.IDs","Intensity.Act.1", "Intensity.Act.2", "Intensity.Act.3")]

# subset rows with at least 2 valid values for Intensity
temp=raw.data[c("Intensity.Act.1", "Intensity.Act.2", "Intensity.Act.3")]
temp=temp!=0
temp=rowSums(temp)
temp=temp>=2
raw.data=subset.data.frame(raw.data,temp)

# obtain the membrane majority protein IDs
mproteinID=as.vector(raw.data$Majority.protein.IDs)

temp=c()
for (i in 1:length(mproteinID))
{
  if (length(intersect(unlist(strsplit(mproteinID[i],';')),as.vector(uniprot$Entry)))>0)
  {
    temp=c(temp,mproteinID[i])
  }
}

membrane=temp

################# membrane proteins in copy number per cell data ##########################################

# load the data
raw.data=read.csv('C:/Users/Gang Wu/Dropbox/R data analysis of Treg/R data analysis Treg/Raw data/Copy number per cell.csv')

# obtain the intensity colomns for Treg rest
raw.data=raw.data[c("Majority.protein.IDs",'Copy.number_Act_.0','Copy.number_Act_.1','Copy.number_Act_.2')]

# subset rows with at least 2 valid values for Intensity
temp=raw.data[c('Copy.number_Act_.0','Copy.number_Act_.1','Copy.number_Act_.2')]
temp=temp!=0
temp=rowSums(temp)
temp=temp>=2
raw.data=subset.data.frame(raw.data,temp)

# obtain the membrane majority protein IDs
mproteinID=as.vector(raw.data$Majority.protein.IDs)

temp=c()
for (i in 1:length(mproteinID))
{
  if (length(intersect(unlist(strsplit(mproteinID[i],';')),as.vector(uniprot$Entry)))>0)
  {
    temp=c(temp,mproteinID[i])
  }
}

cnpc=temp

################# combine, de-redundancy ######################################################################

memproteinID=c(membrane,cnpc)
memproteinID=levels(factor(memproteinID))

temp=c()
for (i in 1:length(memproteinID)-1)
{
  for (j in i+1:length(memproteinID))
  {
    if (length(intersect(unlist(strsplit(memproteinID[i],';')),
                         unlist(strsplit(memproteinID[j],';'))))>0)
    {temp=c(temp,memproteinID[j])}
  }
}
temp=levels(factor(temp))

memproteinID=setdiff(memproteinID,temp)

################ obtain glycomemproteinID ######################################################################

# subset glycosylated rows in uniprot
temp=uniprot.Ngly((uniprot))
uniprot=subset.data.frame(uniprot,temp)

# obtain elements in memproteinID which corresponding to glycoproteins

temp=c()
for (i in 1:length(memproteinID))
{
  if (length(intersect(unlist(strsplit(memproteinID[i],';')),as.vector(uniprot$Entry)))>0)
  {
    temp=c(temp,mproteinID[i])
  }
}

glycomemproteinID=temp

