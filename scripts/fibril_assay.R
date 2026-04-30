# Function to calculate the yield of different fibril topologies. This particular
# fibril fraction function only counts nodes that are part of the interior of a
# fibrillar section of a graph. Although it leads to undercounting on the ends,
# its advantage is that it biases the genetic algorithm to produce longer stretches of fibril.
#
# For documentation on using orca package, see paper: https://www.jstatsoft.org/article/view/v071i10/0 
# Also see ergm graphlets package:
# Omer N. Yaveroglu, Sean M. Fitzhugh, Maciej Kurant, Athina Markopoulou, Carter T. Butts,
# Natasa Przulj. 2013 ergm.graphlets: A Package for ERG Modeling Based on Graphlet Statistics
# http://CRAN.R-project.org/package=ergm.graphlets.
#
# Local system requirements at the time of development required use of orca, but 
# future updates will most likely utilize graphlets. 

library(orca)

fibril_assay<-function(net){
  edgelist<-as.edgelist(net)
  #Get the orbit degrees that we'll need for classification
  od<-as.data.frame(count5(edgelist))
  #  test for 1-ribbon:
  isNode.in.1ribbon<-function(odRow){
    rowTestResult <-((odRow["O15"]==2) & (odRow["O16"]==2) & (odRow["O8"]==0) & (odRow["O7"]==0))
    rowTestResult
  }
  #  test for 2-ribbon:
  isNode.in.2ribbon<-function(odRow){
    rowTestResult <-((odRow["O15"]==12) & (odRow["O16"]==12) & (odRow["O37"]==4))
    rowTestResult
  }
  #  test for 1,2 2-ribbon:
  isNode.in.12.2ribbon<-function(odRow){
    role1 <-(odRow["O59"]==1 & odRow["O60"]==1 & odRow["O13"]==1) 
    role3 <-(odRow["O45"]==2 & odRow["O46"]==2 & odRow["O47"]==2 & odRow["O48"]==4 & odRow["O61"]==1 & odRow["O60"]==2 & odRow["O59"]==2) 
    rowTestResult <- (role1 | role3)
    rowTestResult
  }
  #  tests for 3-prism
  isNode.in.3prism<-function(odRow){
    rowTestResult <-((odRow["O53"]==4) & (odRow["O52"]==2) & (odRow["O51"]==4) & (odRow["O38"]==4) & (odRow["O37"]==8) & (odRow["O36"]==4) &
                       (odRow["O35"]==4) & (odRow["O33"]==1))
    rowTestResult
  }
  #  test for double 1,2 2-ribbon:
  isNode.in.double.12.2ribbon<-function(odRow){
    rowTestResult<-(((odRow["O67"]==4)&(odRow["O66"]==4)&(odRow["O65"]==2))|((odRow["O67"]==0)&(odRow["O66"]==2)&(odRow["O65"]==1)))
  }
  topology.types <- c("1-ribbon", "2-ribbon", "1,2 2-ribbon", "double 1,2 2-ribbon", "3-prism")
  membership <-as.data.frame(matrix(0, nrow=net$gal$n,ncol=length(topology.types)))
  colnames(membership) <- topology.types
  for (i in 1:nrow(od)) {
    if(isNode.in.1ribbon(od[i,])) membership[i, "1-ribbon"] = 1
    if(isNode.in.2ribbon(od[i,])) membership[i, "2-ribbon"] = 1
    if(isNode.in.12.2ribbon(od[i,])) membership[i, "1,2 2-ribbon"] = 1
    if(isNode.in.double.12.2ribbon(od[i,])) membership[i, "double 1,2 2-ribbon"] = 1
    if(isNode.in.3prism(od[i,])) membership[i, "3-prism"] = 1
  }
  return(membership)
}

