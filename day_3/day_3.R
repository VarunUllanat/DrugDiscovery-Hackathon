BiocManager::install("STROMA4")
BiocManager::install("breastCancerMAINZ")
library(breastCancerMAINZ)
library("STROMA4")
data(mainz, package='breastCancerMAINZ')
head(fData(mainz)[, "Gene.symbol", drop=FALSE])
just.stromal.properties <- assign.properties(ESet=mainz, geneID.column="Gene.symbol",genelists="Stroma4", n=10, mc.cores=1)


property.cols <- paste0(c('T', 'B', 'D', 'E'), '.stroma.property')
patient.subtypes <- pData(just.stromal.properties)[, property.cols]

for(i in c('T', 'B', 'D', 'E')){
  patient.subtypes[, paste0(i, '.stroma.property')] <-paste0(i, '-', patient.subtypes[, paste0(i, '.stroma.property')])
}

check_score = function(cond, w=1){
  if (cond == "low"){
    cond_score = w*0
  }
  if (cond == "high"){
    cond_score = w*1
  }
  if (cond == "intermediate"){
    cond_score = w*0.5
  }
  return(cond_score)
}

get_score = function(index, patient_df = patient.subtypes){
  T_property = check_score(strsplit(patient_df$T.stroma.property[index],"-")[[1]][2])
  B_property = check_score(strsplit(patient_df$B.stroma.property[index],"-")[[1]][2])
  D_property = check_score(strsplit(patient_df$D.stroma.property[index],"-")[[1]][2])
  E_property = check_score(strsplit(patient_df$E.stroma.property[index],"-")[[1]][2])
  total_score = T_property + B_property + D_property + E_property
  return(total_score)
}
scores = sapply(1:200, get_score)

hist(scores, xlab = "Scores", main = "Histogram of test scores for Mainz dataset (n=200)")