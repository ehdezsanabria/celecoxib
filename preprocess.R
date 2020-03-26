library(stringr)

# Provide your data path
data <- read.csv2("four.csv")

# If you want Genus instead of family,
#change accordingly in data$ and ask Ruben for further assistance 
data.sorted <- data[with(data, order(str_trim(data$Genus))),]

# Provide number of samples
number_samples = 56

#### Provide starting column of samples ####
start_col = 2

#### row number of last sample !!!!
end_col = number_samples+start_col-1

#### str_trim function removes all spaces to avoid classifying same OTU 
###as diff Genera ####
Genus <-unique(str_trim(data.sorted$Genus))
positions <- matrix(nrow=length(Genus), ncol=3)

#### Find the start and end position of all Genera in the dataframe ####
for(i in 1:length(Genus)){
  positions[i,1] = Genus[i]
  temp = which(sapply(data.sorted$Genus, 
                      function(x) any(str_trim(x) == str_trim(paste(Genus[i])))))
  positions[i,2] = temp[1]
  positions[i,3] = temp[(length(temp))] 
}
Samples.0 = colnames(data[,start_col:(end_col)])


#### Add all the rows (from start to end) of the same Genus together ####
#### end_col + as.integer() function to force integer values of positions-> this I added!!!!
data.pool<-data.frame()
for(i in 1: length(Genus)){
  temp2 = colSums(data.sorted[as.integer(positions[i,2]):as.integer(positions[i,3]),
                              (start_col:(end_col))])
  data.pool = rbind(data.pool,temp2)
}
data.pool = cbind(data.pool, Genus)
Samples = c(Samples.0, "Genus")
colnames(data.pool) = Samples

write.csv2(file="pooled4.csv",data.pool)

####BIODIVERSITY INDICES####

library(vegan)
library(xtable)

#load data for all samples
mydata <- read.csv2('indices.csv',
                    stringsAsFactors = FALSE)

#set up data as data matrix
counts <- as.matrix(mydata[-1])

#alpha diversity indices
shannon <- diversity(counts,index="shannon", MARGIN=1, base=exp(1))
simpson <- diversity(counts,index="simpson", MARGIN=1, base=exp(1))
inversesimp <- diversity(counts,index="inv", MARGIN=1, base=exp(1))

#species richness and eveness indices
totalspecies<-specnumber(counts)
Pielou <- shannon/log(totalspecies)
alpha <- fisher.alpha(counts)

#Export indices to Excel

library(xlsx)
write.xlsx(Pielou, file="indicesCX.xlsx",
           sheetName="Pielou", append=FALSE)
write.xlsx(shannon, file="indicesCX.xlsx", sheetName="Shannon", 
           append=TRUE)
write.xlsx(simpson, file="indicesCX.xlsx", sheetName="Simpson", 
           append=TRUE)
write.xlsx(alpha, file="indicesCX.xlsx", sheetName="alpha", 
           append=TRUE)
write.xlsx(inversesimp, file="indicesCX.xlsx", sheetName="InverseSimpson", 
           append=TRUE)
write.xlsx(totalspecies, file="indicesCX.xlsx", sheetName="TotalSpecies", 
           append=TRUE)

