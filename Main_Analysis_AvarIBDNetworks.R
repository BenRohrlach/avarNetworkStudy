# Network analysis of Avar Cemetries 
# Created by: AB Rohrlach
# Created on: 11-05-2023

#### packages and custom functions ----
if(!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               ggplot2,
               spaa,
               network,
               ggnetwork,
               ergm,
               GERGM,
               philentropy,
               stringi)

source(PATH_TO_FUNCTIONS_SCRIPT)

#### Individual-based regression ----
#### Read in and clean data ----

avr <- read_delim(PATH_TO_METADATA,delim=';') %>%
  janitor::clean_names() %>%
  dplyr::mutate(yhap_short=ifelse(sex=='M',yhap_short,'-')) %>%
  dplyr::filter(!is.na(belt)) %>%
  dplyr::filter(grepl('Avar',period)) %>%
  dplyr::mutate(belt_or_coatclasp=paste0(ifelse((belt=='YES')|(coat_clasp=='YES'),'Wearing Belt Set\nor Coast Clasp','No Belt Set\nor Coast Clasp')))


avr.ibd <- read_delim(PATH_TO_IBD_DATA,delim='\t') %>%
  janitor::clean_names() %>%
  dplyr::filter((iid1%in%avr$iid)&(iid2%in%avr$iid))  %>%
  dplyr::mutate(rel=ifelse((n_ibd_12>=2)&(n_ibd_16>=1),1,0))

varOI <- 'n_ibd_16'
varCO <- 0
Dij <- df2dist(avr.ibd$iid1,avr.ibd$iid2,dplyr::pull(avr.ibd,varOI))
Bij <- df2dist(avr.ibd$iid1,avr.ibd$iid2,dplyr::pull(avr.ibd,rel))
network.metadata <- tibble(iid=row.names(Dij)) %>%
  dplyr::left_join(avr,by=c('iid'='iid'))

Nij <- (Bij>0&!is.infinite(Bij)) %>%
  network(directed=F) %>%
  network::set.vertex.attribute('sex',network.metadata$sex) %>%
  network::set.vertex.attribute('site',network.metadata$site) %>%
  network::set.vertex.attribute('ornament',network.metadata$belt_or_coatclasp) %>%
  network::set.edge.value("IBD",distRules(Dij,(Dij>=varCO&!is.infinite(Dij))))

ergm.obj <- makeModelList('Nij',
                          varStr=c('sex','ornament'),
                          cnst="nodematch('site')")
ergm.list <- ergm.obj$models

breaks <- 1:2
n_breaks <- length(breaks)
labels <- c('','','No Belt Set\nor Coast Clasp','Wearing Belt Set\nor Coast Clasp')
shapes <- c(rep(21, n_breaks), rep(22, n_breaks))
breaks2 <- rep(breaks, 2)

set.seed(123)
network.gg <- ggplot(data = Nij, aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_blank()+
  geom_edges(aes(col=IBD)) +
  geom_nodes(aes(size=as.factor(ornament)%>%as.numeric(),fill=site,shape=sex)) +
  scale_shape_manual(values=rep(21:25,100),
                     name='Genetic\nSex')+
  guides(fill = guide_legend(override.aes = list(shape=21,size=3)),
         shape = guide_legend(override.aes = list(size=3)))+
  scale_colour_gradient(low = "grey",high = "black")+
  scale_size_continuous(breaks = breaks2, labels = labels,name='Ornament',range=c(2,4),
                        guide = guide_legend(order = 2, nrow = 2, byrow = FALSE,
                                             override.aes = list(shape = shapes),
                                             label.vjust = -.5)) +
  scale_fill_manual(values=c('#E95032','#05582A','#941443'),
                    name='Site')+
  theme(legend.position='left')

network.aic <- ERGMS.aicPlot(ergm.list,
                             form_v=ergm.obj$formulae %>% 
                               gsub('Nij ~ edges\\+|Nij ~ edges \\+ ','',.) %>%
                               gsub('node','',.) %>%
                               gsub('factor','',.) %>%
                               gsub("'",'',.),
                             ang=15,tsize=8,printRes=T)

best.model <- ergm.list[[bestModelErgm(ergm.list)]]
summary(best.model)
gof.fit <- gof(best.model,unconditional=F)
plot(gof.fit)

#### Site-based analysis ----

site <- read_delim(PATH_TO_SITE_DATA,delim=',') %>%
  janitor::clean_names() 

PATH_TO_SITE_META_DATA <- '~/Dropbox/Project_Avar_Networks/Data/Forupload/SITE_META_DATA.csv'
site.avr <- read_delim(PATH_TO_SITE_META_DATA,delim=';') %>%
  janitor::clean_names() %>%
  dplyr::filter(grepl('Avar',period))  %>%
  dplyr::select(iid,site)

site.ibd.before <- read_delim(PATH_TO_IBD_DATA,delim='\t') %>%
  janitor::clean_names()  %>%
  dplyr::left_join(site.avr,by=c('iid1'='iid')) %>%
  dplyr::rename(site1=site) %>%
  dplyr::left_join(site.avr,by=c('iid2'='iid')) %>%
  dplyr::rename(site2=site) %>%
  dplyr::filter(site1%in%site$site,site2%in%site$site) %>%
  dplyr::mutate(valid=(n_ibd_12>=2)&(n_ibd_16>=1)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(sites=paste0(sort(c(site1,site2)),collapse='/')) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(sites) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(site1!=site2,n>=0)

site.ibd <- site.ibd.before %>%
  dplyr::group_by(sites) %>%
  dplyr::summarise(n=n(),across(max_ibd:n_ibd_20, mean)) %>%
  dplyr::mutate(site1=gsub('/.*','',sites),
                site2=gsub('.*/','',sites))

gij <- site %>%
  dplyr::rowwise() %>%
  dplyr::filter(any(grepl(site,site.ibd))) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(site) %>%
  as.data.frame() %>%
  dplyr::select(-ancestry) %>%
  column_to_rownames(var="site") %>%
  as.matrix() %>%
  distance(method="euclidean",
           use.row.names=TRUE)

responses <- c("max_ibd","sum_ibd_8","n_ibd_8","sum_ibd_12","n_ibd_12","sum_ibd_16","n_ibd_16","sum_ibd_20","n_ibd_20")
dij_8 <- df2dist(site.ibd$site1,site.ibd$site2,dplyr::pull(site.ibd,n_ibd_8),default=0) # Change IBD measure here ----
dij_12 <- df2dist(site.ibd$site1,site.ibd$site2,dplyr::pull(site.ibd,n_ibd_12),default=0) # Change IBD measure here ----
dij_16 <- df2dist(site.ibd$site1,site.ibd$site2,dplyr::pull(site.ibd,n_ibd_16),default=0) # Change IBD measure here ----
dij_20 <- df2dist(site.ibd$site1,site.ibd$site2,dplyr::pull(site.ibd,n_ibd_20),default=0) # Change IBD measure here ----

site.md <- tibble(site=row.names(dij_8)) %>%
  dplyr::left_join(site,by=c('site'='site'))

#### Set nij ----
nij <- (dij_8>0&!is.infinite(dij_8)) %>%
  network(directed=F) %>%
  network::set.vertex.attribute('site',site.md$site) %>%
  network::set.vertex.attribute('ancestry',site.md$ancestry) %>%
  network::set.vertex.attribute('lat',site.md$lat) %>%
  network::set.vertex.attribute('lon',site.md$lon) %>%
  network::set.edge.value("IBD",dij_8) %>%
  network::set.edge.value("Geo",gij)
ggc <- ggplot(data = nij, aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_blank()+
  geom_edges(aes(col=IBD))+
  geom_nodes(aes(fill=site,shape=site),size=7.5)
site.matrix <- site%>%bind_cols(ggplot_build(ggc)[[1]][[2]])%>%column_to_rownames(var="site")%>%dplyr::select(x,y)%>%as.matrix()

nij_8 <- (dij_8>0&!is.infinite(dij_8)) %>%
  network(directed=F) %>%
  network::set.vertex.attribute('site',site.md$site) %>%
  network::set.vertex.attribute('ancestry',site.md$ancestry) %>%
  network::set.vertex.attribute('lat',site.md$lat) %>%
  network::set.vertex.attribute('lon',site.md$lon) %>%
  set.edge.value("IBD",dij_8) %>%
  set.edge.value("Geo",gij)
tij_8 <- as.matrix(nij_8)
nij_12 <- (dij_12>0&!is.infinite(dij_12)) %>%
  network(directed=F) %>%
  network::set.vertex.attribute('site',site.md$site) %>%
  network::set.vertex.attribute('ancestry',site.md$ancestry) %>%
  network::set.vertex.attribute('lat',site.md$lat) %>%
  network::set.vertex.attribute('lon',site.md$lon) %>%
  set.edge.value("IBD",dij_12) %>%
  set.edge.value("Geo",gij)
tij_12 <- as.matrix(nij_12)
nij_16 <- (dij_16>0&!is.infinite(dij_16)) %>%
  network(directed=F) %>%
  network::set.vertex.attribute('site',site.md$site) %>%
  network::set.vertex.attribute('ancestry',site.md$ancestry) %>%
  network::set.vertex.attribute('lat',site.md$lat) %>%
  network::set.vertex.attribute('lon',site.md$lon) %>%
  set.edge.value("IBD",dij_16) %>%
  set.edge.value("Geo",gij)
tij_16 <- as.matrix(nij_16)
nij_20 <- (dij_20>0&!is.infinite(dij_20)) %>%
  network(directed=F) %>%
  network::set.vertex.attribute('site',site.md$site) %>%
  network::set.vertex.attribute('ancestry',site.md$ancestry) %>%
  network::set.vertex.attribute('lat',site.md$lat) %>%
  network::set.vertex.attribute('lon',site.md$lon) %>%
  set.edge.value("IBD",dij_20) %>%
  set.edge.value("Geo",gij)
tij_20 <- as.matrix(nij_20)
nij_max <- (dij_max>0&!is.infinite(dij_max)) %>%
  network(directed=F) %>%
  network::set.vertex.attribute('site',site.md$site) %>%
  network::set.vertex.attribute('ancestry',site.md$ancestry) %>%
  network::set.vertex.attribute('lat',site.md$lat) %>%
  network::set.vertex.attribute('lon',site.md$lon) %>%
  set.edge.value("IBD",dij_max) %>%
  set.edge.value("Geo",gij)
tij_max <- as.matrix(nij_max)

ggdefault <- ggplot(data = nij_8, aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_blank()+
  geom_edges(aes(col=IBD))+
  geom_nodes(aes(fill=site,shape=site),size=7.5)+
  scale_shape_manual(values=rep(21:25,100))+
  scale_size_manual(values=c(2,4))+
  scale_colour_gradient(low = "grey",high = "black",guide='none')+
  theme(legend.position='left')


gg8 <- ggplot(data = ggnetwork(nij_8,layout=site.matrix), aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_blank()+
  geom_edges(aes(col=IBD))+
  geom_nodes(aes(fill=site,shape=site),size=7.5)+
  scale_shape_manual(values=rep(21:25,100),
                     guide='none')+
  scale_size_manual(values=c(2,4))+
  scale_colour_gradient(low = "grey",high = "black")+
  scale_fill_discrete(guide='none')+
  theme(legend.position='left')+
  ggtitle("n_IBD_8")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
gg12 <- ggplot(data = ggnetwork(nij_12,layout=site.matrix), aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_blank()+
  geom_edges(aes(col=IBD))+
  geom_nodes(aes(fill=site,shape=site),size=7.5)+
  scale_shape_manual(values=rep(21:25,100),
                     guide='none')+
  scale_size_manual(values=c(2,4))+
  scale_colour_gradient(low = "grey",high = "black")+
  scale_fill_discrete(guide='none')+
  theme(legend.position='left')+
  ggtitle("n_IBD_12")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
gg16 <- ggplot(data = ggnetwork(nij_16,layout=site.matrix), aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_blank()+
  geom_edges(aes(col=IBD))+
  geom_nodes(aes(fill=site,shape=site),size=7.5)+
  scale_shape_manual(values=rep(21:25,100),
                     guide='none')+
  scale_size_manual(values=c(2,4))+
  scale_colour_gradient(low = "grey",high = "black")+
  scale_fill_discrete(guide='none')+
  theme(legend.position='left')+
  ggtitle("n_IBD_16")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
gg20 <- ggplot(data = ggnetwork(nij_20,layout=site.matrix), aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_blank()+
  geom_edges(aes(col=IBD))+
  geom_nodes(aes(fill=site,shape=site),size=7.5)+
  scale_shape_manual(values=rep(21:25,100),
                     guide='none')+
  scale_size_manual(values=c(2,4))+
  scale_colour_gradient(low = "grey",high = "black")+
  scale_fill_discrete(guide='none')+
  theme(legend.position='left')+
  ggtitle("n_IBD_20")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

ggmax <- ggplot(data = ggnetwork(nij_max,layout=site.matrix), aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_blank()+
  geom_edges(aes(col=IBD))+
  geom_nodes(aes(fill=site,shape=site),size=7.5)+
  scale_shape_manual(values=rep(21:25,100),
                     guide='none')+
  scale_size_manual(values=c(2,4))+
  scale_colour_gradient(low = "grey",high = "black")+
  scale_fill_discrete(guide='none')+
  theme(legend.position='left')+
  ggtitle("max_IBD")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

set.seed(12345)
input_str <- 'nij_max'
modList1 <- makeModelList(input_str,varStr=c('ancestry'),cnst="edgecov(gij)")
modList2 <- makeModelList(input_str,varStr=c('ancestry'))
modList <- list(models=c(modList1$models,modList2$models),
            formulae=c(modList1$formulae,modList2$formulae))
ERGMS.aicPlot(modList$models,form_v=modList$formulae,ang=15,tsize=6,printRes=T)
modList$models[[bestModelErgm(modList$models)]] %>% summary()
site.erg_8 <- ergm(nij_8 ~ edges + nodefactor('ancestry'))
site.erg_12 <- ergm(nij_12 ~ edges + nodefactor('ancestry'))
site.erg_16 <- ergm(nij_16 ~ edges + nodefactor('ancestry'))
site.erg_20 <- ergm(nij_20 ~ edges + nodefactor('ancestry'))
site.erg_max <- ergm(nij_max ~ edges + nodefactor('ancestry'))

mij <- site.md %>%
  as.data.frame() %>%
  column_to_rownames(var="site") 
formula_8 <- dij_8 ~ edges + nodemix("ancestry",base="european")
formula_12 <- dij_12 ~ edges + nodemix("ancestry",base="european")
formula_16 <- dij_16 ~ edges + nodemix("ancestry",base="european")
formula_20 <- dij_20 ~ edges + nodemix("ancestry",base="european")

test_8 <- gergm(formula_8,
                network_is_directed=FALSE,
                covariate_data = mij,
                number_of_networks_to_simulate = 1000000,
                thin = 1/1000,
                proposal_variance = 0.005,
                MCMC_burnin = 10000,
                seed = 456,
                convergence_tolerance = 0.05,
                generate_plots=F)
test_12 <- gergm(formula_12,
                 network_is_directed=FALSE,
                 covariate_data = mij,
                 number_of_networks_to_simulate = 1000000,
                 thin = 1/1000,
                 proposal_variance = 0.005,
                 MCMC_burnin = 10000,
                 seed = 456,
                 convergence_tolerance = 0.05,
                 generate_plots=F)
test_16 <- gergm(formula_16,
                 network_is_directed=FALSE,
                 covariate_data = mij,
                 number_of_networks_to_simulate = 1000000,
                 thin = 1/1000,
                 proposal_variance = 0.005,
                 MCMC_burnin = 10000,
                 seed = 456,
                 convergence_tolerance = 0.05,
                 generate_plots=F)
test_20 <- gergm(formula_20,
                 network_is_directed=FALSE,
                 covariate_data = mij,
                 number_of_networks_to_simulate = 1000000,
                 thin = 1/1000,
                 proposal_variance = 0.005,
                 MCMC_burnin = 10000,
                 seed = 456,
                 convergence_tolerance = 0.05,
                 generate_plots=F)


# Generate Estimate Plot
coef.gergm(test_8)
coef.gergm(test_12)
coef.gergm(test_16)
coef.gergm(test_20)

test2_8 <- conditional_edge_prediction(
  GERGM_Object = test_8,
  number_of_networks_to_simulate = 100,
  thin = 1,
  proposal_variance = 0.05,
  MCMC_burnin = 1000,
  seed = 123)
test2_12 <- conditional_edge_prediction(
  GERGM_Object = test_12,
  number_of_networks_to_simulate = 100,
  thin = 1,
  proposal_variance = 0.05,
  MCMC_burnin = 1000,
  seed = 123)
test2_16 <- conditional_edge_prediction(
  GERGM_Object = test_16,
  number_of_networks_to_simulate = 100,
  thin = 1,
  proposal_variance = 0.05,
  MCMC_burnin = 1000,
  seed = 123)
test2_20 <- conditional_edge_prediction(
  GERGM_Object = test_20,
  number_of_networks_to_simulate = 100,
  thin = 1,
  proposal_variance = 0.05,
  MCMC_burnin = 1000,
  seed = 123)

MSE_results_8 <- conditional_edge_prediction_MSE(test2_8)
MSE_results_12 <- conditional_edge_prediction_MSE(test2_12)
MSE_results_16 <- conditional_edge_prediction_MSE(test2_16)
MSE_results_20 <- conditional_edge_prediction_MSE(test2_20)


