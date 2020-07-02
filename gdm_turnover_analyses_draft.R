#calculate turnover and nestedness
library(betapart)
library(phyloseq)
library(ape)
library(dplyr)

tree<-read.tree("~/Desktop/CB_time_series/psabs_al_nonoise.nwk")
temp<-read.table("~/Desktop/CB_time_series/CB_time_series_temp_dist.txt")
temp1<-temp[2:211,]
row.names(temp1)<-temp1$V1
temp1$V1<-NULL
colnames(temp1)<-row.names(temp1)

ps_arc<-subset_taxa(ps_nonoise,Kingdom=="k__Archaea")#starting with ps_nonoise: relative abundances, singletons and doubletons removed
out<-taxa_names(ps_arc)
rooted<-root(tree,outgroup=out,resolve.root=TRUE)#root tree with archea as outgroup

tab_nonoise<-as.matrix(as.data.frame((otu_table(ps_nonoise))))#rel. abundances
core<-betapart.core.abund(tab_nonoise)
pair<-beta.pair.abund(core)
multi<-beta.multi.abund(core)

tab_nonoise[tab_nonoise>0]<-1#presence/absence
corePhylo<-phylo.betapart.core(tab_nonoise,rooted)
pairPhylo<-phylo.beta.pair(corePhylo)
multiPhylo<-phylo.beta.multi(corePhylo)

tab_rare<-as.matrix(as.data.frame(otu_table(p_rare)))#rarified absolute counts, singletons and doubletons removed
coreRare<-betapart.core.abund(tab_rare)
pairRare<-beta.pair.abund(coreRare)
multiRare<-beta.multi.abund(coreRare)

temp2<-temp1[row.names(temp1) %in% row.names(tab_nonoise),]
temp3<-temp2[,colnames(temp2) %in% row.names(tab_nonoise)]
temp_dis<-as.dist(temp3)

modSor<-decay.model(pairPhylo$phylo.beta.sor,temp_dis,model.type="exponential",y.type="dissimilarities",perm=100)
modSIM<-decay.model(pairPhylo$phylo.beta.sim,temp_dis,model.type="exponential",y.type="dissimilarities",perm=100)
modSNE<-decay.model(pairPhylo$phylo.beta.sne,temp_dis,model.type="exponential",y.type="dissimilarities",perm=100)

plot.decay(modSor,col="black",remove.dots=TRUE)
plot.decay(modSIM,col="red",remove.dots=TRUE,add=TRUE)
plot.decay(modSNE,col="blue",remove.dots=TRUE,add=TRUE)

#spatial distance matrix
library(geosphere)
p1<-sample_data(ps_nonoise)$Longitude
p2<-sample_data(ps_nonoise)$Latitude
pp<-cbind(p1,p2)
spat_dis<-distm(pp,fun=distVincentyEllipsoid)
row.names(spat_dis)<-sample_names(ps_nonoise)
colnames(spat_dis)<-sample_names(ps_nonoise)
spat_dis<-as.dist(spat_dis)

modSorspace<-decay.model(pairPhylo$phylo.beta.sor,spat_dis,model.type="exponential",y.type="dissimilarities",perm=100)
modSIMspace<-decay.model(pairPhylo$phylo.beta.sim,spat_dis,model.type="exponential",y.type="dissimilarities",perm=100)
modSNEspace<-decay.model(pairPhylo$phylo.beta.sne,spat_dis,model.type="exponential",y.type="dissimilarities",perm=100)

plot.decay(modSorspace,col="black",remove.dots=TRUE)
plot.decay(modSIMspace,col="red",remove.dots=TRUE,add=TRUE)
plot.decay(modSNEspace,col="blue",remove.dots=TRUE,add=TRUE)

#GDM modelling
library(gdm)
#use ps_prev fot this analysis
tab<-as.data.frame(otu_table(ps_prev))#only ASVs present in at least 10% of samples
tab$site<-row.names(tab)
meta_m<-sample_data(ps_prev)
meta_m$Year<-as.numeric(meta_m$Year)
meta_m$Depth<-as.numeric(meta_m$Depth)
meta_m$Latitude<-as.numeric(meta_m$Latitude)
meta_m$Year<-as.numeric(meta_m$Year)
meta_m$Bacteria_mL<-as.numeric(meta_m$Bacteria_mL)
meta_m$Phyto_mL<-as.numeric(meta_m$Phyto_mL)
meta_m$Picophyto_mL<-as.numeric(meta_m$Picophyto_mL)
meta_m$Nanophyto_mL<-as.numeric(meta_m$Nanophyto_mL)
meta_m$Temperature<-as.numeric(meta_m$Temperature)
meta_m$Salinity<-as.numeric(meta_m$Salinity)
meta_m$Nitrate<-as.numeric(meta_m$Nitrate)
meta_m$Silicate<-as.numeric(meta_m$Silicate)
meta_m$Phosphate<-as.numeric(meta_m$Phosphate)
meta_m$Chlorophyll<-as.numeric(meta_m$Chlorophyll)
meta_m<-as.data.frame(as.matrix(meta_m))
env_var<-meta_m[c(1,3,4,6,7,8,12,13,14,15,16)]

env_vars<-na.omit(env_var)
env_vars$Year<-as.numeric(env_vars$Year)
env_vars$Depth<-as.numeric(env_vars$Depth)
env_vars$Bacteria_mL<-as.numeric(env_vars$Bacteria_mL)
env_vars$Temperature<-as.numeric(env_vars$Temperature)
env_vars$Salinity<-as.numeric(env_vars$Salinity)
env_vars$Nitrate<-as.numeric(env_vars$Nitrate)
env_vars$Silicate<-as.numeric(env_vars$Silicate)
env_vars$Phosphate<-as.numeric(env_vars$Phosphate)
envz<-env_vars %>% mutate_at(c(2,3,6,7,8,9,10,11),funs(scale(.)))#z scale the data
names(envz)<-c("site","Year","Depth","Long","Lat","Bacteria_ml","Temperature","Salinity","Nitrate","Silicate","Phosphate")
tab2<-tab[row.names(tab) %in% envz$site,]

SP<-formatsitepair(tab2,1,siteColumn="site",predData=envz,abundance=TRUE,XColumn="Long",YColumn="Lat")
gdm_mod<-gdm(SP,geo=TRUE)
summary.gdm(gdm_mod)
plot(gdm_mod,plot.layout=c(4,4))