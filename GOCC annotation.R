rm(list = ls()) # remove all objects from the environment

########## set working directory ###################################################

# setwd("~/Documents/R data analysis Treg/bead.solution/SnAct")

########## functions ##########################################################################

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

# find cytoplasm annotation in uniprot, returns a logical value
uniprot.cp=function(uniprot)
{
  ans=vector()
  for (i in 1:nrow(uniprot))
  {
    ans[i]=FALSE
    temp=unlist(strsplit(as.vector(uniprot$Gene.ontology..cellular.component.[i]),' '))
    if ('cytoplasm' %in% temp)
    {
      ans[i]=TRUE
    }
  }
  return(ans)
}

# find nucleus annotation in uniprot,returns a logical value
uniprot.nuc=function(uniprot)
{
  ans=vector()
  for (i in 1:nrow(uniprot))
  {
    ans[i]=FALSE
    temp=unlist(strsplit(as.vector(uniprot$Gene.ontology..cellular.component.[i]),' '))
    if ('nucleus' %in% temp)
    {
      ans[i]=TRUE
    }
  }
  return(ans)
}


# imputation
imputation=function(x=vector.for.input, down.shift, width) # funciton for imputation of a vector
{
  a=x[!(is.na(x))] #obtain the valid data
  b=rnorm(n=length(x[is.na(x)]),mean=mean(a)-down.shift*sd(a),sd=width*sd(a)) # produce the random number
  x[is.na(x)]=b # replace NA values with the random numbers
  return(x)
}

########## load counter receptor data ##############################################

counter.receptor=read.csv('Counter receptors',row.names=1)

########## obtain subset of uniprot with GOCC annotatioin ###############

# load uniprot data
uniprot=read.delim('~/Documents/R data analysis Treg/Raw data/Uniprot lib.txt')
uniprot.temp=subset.data.frame(uniprot, Gene.ontology..cellular.component.!='')

# obtain subset of uniprot with plasma membrane annotation
plasma.membrane=subset.data.frame(uniprot.temp, uniprot.pm(uniprot.temp)==TRUE | Entry=='P00433')

# obtain subset of uniprot with cytoplasm annotation
cytoplasm=subset.data.frame(uniprot.temp, uniprot.cp(uniprot.temp)==TRUE)

# obtain subset of uniprot with neucleus annotation
nucleus=subset.data.frame(uniprot.temp, uniprot.nuc(uniprot.temp)==TRUE)

########## build logical data frame ############################################### 
plasma.membrane.logic=vector()
for (i in 1:nrow(counter.receptor))
{
  if (length(intersect(unlist(strsplit(row.names(counter.receptor)[i],';')),as.vector(plasma.membrane$Entry)))>0)
  {plasma.membrane.logic[i]=TRUE}
  else {plasma.membrane.logic[i]=FALSE}
}

cytoplasm.logic=vector()
for (i in 1:nrow(counter.receptor))
{
  if (length(intersect(unlist(strsplit(row.names(counter.receptor)[i],';')),as.vector(cytoplasm$Entry)))>0)
  {cytoplasm.logic[i]=TRUE}
  else {cytoplasm.logic[i]=FALSE}
}

nucleus.logic=vector()
for (i in 1:nrow(counter.receptor))
{
  if (length(intersect(unlist(strsplit(row.names(counter.receptor)[i],';')),as.vector(nucleus$Entry)))>0)
  {nucleus.logic[i]=TRUE}
  else {nucleus.logic[i]=FALSE}
}

# build logical data frame
GOCC=data.frame(plasma.membrane.logic, cytoplasm.logic, nucleus.logic)
other.logic=rowSums(GOCC)==0
GOCC['other.logic']=other.logic

########## GOCC annotation ##########################################################
GOCC.annotation=GOCC$other.logic
GOCC.annotation[GOCC$other.logic]='other'
GOCC.annotation[GOCC$plasma.membrane.logic]='plasma membrane'
GOCC.annotation[(GOCC$cytoplasm.logic-GOCC$plasma.membrane.logic)==1]='cytoplasm'
GOCC.annotation[(GOCC$nucleus.logic-GOCC$plasma.membrane.logic-GOCC$cytoplasm.logic)==1]='nucleus'

counter.receptor['GOCC annotation']=GOCC.annotation
