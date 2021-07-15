library(rcdk)
library(rJava)
library(chemometrics)
library(fingerprint)
library(cluster)
targets_1 = read.table("C:/drug discovery hackathon/day_4/membrane_receptor-substances.txt", sep= "\t", as.is = c(1,2))
targets_1_list = targets_1[0:1000,"V2"]

target.mol = parse.smiles(targets_1_list)[[1]]

query.smiles = c('NC(=N)NCc1cccc([131I])c1')
query.mol =  parse.smiles(query.smiles)[[1]]

query.fp <- get.fingerprint(query.mol, type='circular')
target.fps <- lapply(target.mols, get.fingerprint, type='circular')
sims <- data.frame(sim=do.call(rbind, lapply(target.fps,
                                             fingerprint::distance,
                                             fp2=query.fp, method='tanimoto')))
subset(sims, sim >= 0.3)



fps = lapply(c(query.mol, target.mols), get.fingerprint, type = 'extended')
fp.sim <- fingerprint::fp.sim.matrix(fps, method='tanimoto')
fp.dist <- 1 - fp.sim

