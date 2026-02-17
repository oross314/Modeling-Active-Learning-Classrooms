library(here) 
library(igraph) 
library(ggraph) 
library(tidygraph) 
library(netrankr)
library(pals)
library(ergm)
library(tergm)
library(statnet)
library(statnet.common)
library(network)
library(sna)
library(ggplot2)
library(igraph) 
library(network) 
library(tidyselect)
library(tidyverse) 
library(dplyr)
library(stringr)
library(networkDynamic)
library(ggpubr)
library(readr)
library(broom)

### Function to plot pre and post network for each course
get.network<- function(method, course){
  nodes<-read.csv(paste0("Roster_",method,"_Course",course,".csv"))
  pre.edges<-read.csv(paste0("PRE_Edgelist_DeIdentified_",method,"_Course",course,".csv"))
  pre.net <- graph_from_data_frame(d=pre.edges, vertices=nodes, directed=F)
  pre.net <- simplify(pre.net) 
  
  post.edges<-read.csv(paste0("POST_Edgelist_DeIdentified_",method,"_Course",course,".csv"))
  post.net <- graph_from_data_frame(d=post.edges, vertices=nodes, directed=F)
  post.net <- simplify(post.net)
  
  par(mfrow=c(1,2))
  par(mar = c(4, 4, 1, 1)) 
  
  l <- layout_with_fr(pre.net)
  plot(pre.net, edge.arrow.size=0,vertex.label=NA,vertex.color="red", edge.color="black",vertex.size = 4, layout=l,edge.curved=F)
  
  l <- layout_with_fr(post.net)
  plot(post.net, edge.arrow.size=0,vertex.label=NA,vertex.color="red", edge.color="black",vertex.size = 4, layout=l,edge.curved=F)
}

get.network("ISLE","1")
get.network("ISLE","2")
get.network("ISLE","3")

get.network("PeerInstruction","1")
get.network("PeerInstruction","2")
get.network("PeerInstruction","3")
get.network("PeerInstruction","4")
get.network("PeerInstruction","5")
get.network("PeerInstruction","6")

get.network("Tutorials","1")
get.network("Tutorials","2")
get.network("Tutorials","3")
get.network("Tutorials","4")
get.network("Tutorials","5")

get.network("SCALEUP","1")
get.network("SCALEUP","2")
get.network("SCALEUP","3")
get.network("SCALEUP","4")
get.network("SCALEUP","5")

dev.off()

### Function to run TERGM for each course, including goodness-of-fit plots
run_tergm<-function(method,course){
  nodes<-read.csv(paste0("Roster_",method,"_Course",course,".csv"))
  n<-length(nodes$ID)
  pre.edges<-read.csv(paste0("PRE_Edgelist_DeIdentified_",method,"_Course",course,".csv"))
  pre.edges<-pre.edges[pre.edges$from%in%nodes$ID,]
  pre.edges<-pre.edges[pre.edges$to%in%nodes$ID,]
  post.edges<-read.csv(paste0("POST_Edgelist_DeIdentified_",method,"_Course",course,".csv"))
  post.edges<-post.edges[post.edges$from%in%nodes$ID,]
  post.edges<-post.edges[post.edges$to%in%nodes$ID,]
  
  pre.net <- graph_from_data_frame(d=pre.edges, vertices=nodes, directed=F)
  pre.net <- simplify(pre.net) 
  
  post.net <- graph_from_data_frame(d=post.edges, vertices=nodes, directed=F)
  post.net <- simplify(post.net) 
  
  pre.net <- intergraph::asNetwork(pre.net)
  post.net <- intergraph::asNetwork(post.net)
  
  networks <- networkDynamic(network.list = list(pre.net,post.net))
  
  set.seed(9)
  stergm <- tergm(networks ~ Form(~edges+
                                    gwdsp(decay=0.5,fixed=TRUE)+
                                    gwesp(decay=0.5,fixed=TRUE)+
                                    gwdegree(decay=0.5,fixed=TRUE)),estimate = "CMLE", times = 0:1)
  
  model.df <- tidy(stergm)
  
  coef<-model.df %>% 
    mutate(or = exp(estimate),  # Odds ratio
           var.diag = diag(vcov(stergm)),  # Variance of each coefficient
           or.se = sqrt(or^2 * var.diag))
  
  par(mfrow=c(1,2))
  par(mar = c(4, 4, 1, 1)) 
  plot(gof(stergm,GOF =~ degree + espartners - model),main=NULL,cex.axis=1)
  
  coef <- coef %>% mutate(Course=course,Method=method)
  coef<-cbind(coef,c("Edges","Chains","Triangles","DegreeDist"))
  
  colnames(coef)<-c("RTerm","Estimate","Std. Error","MCMC%","Z Value","P Value","OddsRatio","Var","OddsRatioSE","Course","Method","Term")
  
  return(coef)
}

### Run and store TERGM results

mat = matrix(ncol = 12, nrow = 0)
tergm_results = data.frame(mat)

tergm_results<-rbind(tergm_results,run_tergm("ISLE","1"))
tergm_results<-rbind(tergm_results,run_tergm("ISLE","2"))
tergm_results<-rbind(tergm_results,run_tergm("ISLE","3"))

tergm_results<-rbind(tergm_results,run_tergm("PeerInstruction","1"))
tergm_results<-rbind(tergm_results,run_tergm("PeerInstruction","2")) 
tergm_results<-rbind(tergm_results,run_tergm("PeerInstruction","3"))
tergm_results<-rbind(tergm_results,run_tergm("PeerInstruction","4"))
tergm_results<-rbind(tergm_results,run_tergm("PeerInstruction","5"))
tergm_results<-rbind(tergm_results,run_tergm("PeerInstruction","6"))

tergm_results<-rbind(tergm_results,run_tergm("Tutorials","1")) 
tergm_results<-rbind(tergm_results,run_tergm("Tutorials","2")) 
tergm_results<-rbind(tergm_results,run_tergm("Tutorials","3"))
tergm_results<-rbind(tergm_results,run_tergm("Tutorials","4"))
tergm_results<-rbind(tergm_results,run_tergm("Tutorials","5"))

tergm_results<-rbind(tergm_results,run_tergm("SCALEUP","1")) 
tergm_results<-rbind(tergm_results,run_tergm("SCALEUP","2")) 
tergm_results<-rbind(tergm_results,run_tergm("SCALEUP","3"))
tergm_results<-rbind(tergm_results,run_tergm("SCALEUP","4"))
tergm_results<-rbind(tergm_results,run_tergm("SCALEUP","5"))

tergm_results <- tergm_results[,c(6,7,9:12)]
colnames(tergm_results) <- c("PVal","OddsRatio","OddsRatio_Error","Course","Method","Term")

dev.off()

### Plot TERGM results

tergm_results <- tergm_results %>% mutate(Sig=ifelse(tergm_results$PVal<0.05,"Significant (p<0.05)","Not Significant"))

tergm_results[tergm_results=="DegreeDist"]<-"Degree Distribution"
tergm_results[tergm_results=="PeerInstruction"]<-"Peer Instruction"
tergm_results[tergm_results=="SCALEUP"]<-"SCALE-UP"
tergm_results$Term<-factor(tergm_results$Term,levels=c("Edges","Chains","Triangles","Degree Distribution"))
tergm_results$Method<- factor(tergm_results$Method,levels=c("ISLE","Peer Instruction","Tutorials","SCALE-UP"))

set.seed(32)
p<-ggplot() +
  geom_hline(yintercept=1, color = "black")+
  geom_boxplot(data=tergm_results[tergm_results$Term=="Chains"|tergm_results$Term=="Triangles"|tergm_results$Term=="Degree Distribution",],aes(x=Term, y=OddsRatio,fill=factor(Method)),size=0.2,alpha=1,outlier.shape=NA,colour="gray40",position=position_dodge(width=1))+
  geom_point(data=tergm_results[tergm_results$Term=="Chains"|tergm_results$Term=="Triangles"|tergm_results$Term=="Degree Distribution",],aes(x=Term, y=OddsRatio,colour=factor(Method),shape=Sig,group=factor(Method)),position=position_jitterdodge(dodge.width=1,jitter.width=0.1),size=2,alpha=1)+
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.y = element_line(color = "black", size = 0.5),axis.ticks.x=element_line(color = "black", size = 0.5),axis.title = element_text(size=9,color='black'),axis.text = element_text(size=8,color='black'),axis.text.x = element_text(size=8,angle=0),panel.spacing = unit(.15, "lines")) + theme(strip.background = element_rect(colour='black', size=0.5),strip.text=element_text(size=8,color='black'))+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +theme(legend.position = c(0.21,0.78),legend.title=element_blank(),legend.text = element_text(size=7,color="black"),legend.box.background = element_rect(colour = "black",fill="white"))+ylab("Odds Ratio")+scale_colour_manual(values=c("#56b3e9","#cc79a7","#e69f00","#009e74"))+scale_shape_manual(values=c(16,17))+xlab("Network Variable")+scale_x_discrete(labels = function(x) str_wrap(x, width = 9))+scale_y_continuous(breaks=c(-2,1,0,1,2,3,4,5,6,7,8))+ theme(axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),legend.key.size = unit(-1, 'lines'),legend.spacing.y = unit(-0.2, "cm")) +scale_fill_manual(values=c("white","white","white","white"))

p

ggarrange(p) %>%
  ggexport(filename = "TERGMs.pdf",width=3.5,height=2.5,units="in")
