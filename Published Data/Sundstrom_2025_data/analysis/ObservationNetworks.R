library(igraph) 
library(dplyr) 
library(ggplot2)
library(ggnewscale)
library(gridExtra)
library(grid)
library(tidyselect)
library(tidyverse) 
library(effsize)
library(metafor)
library(ggpubr)
library(ggsignif)
library(reshape2)
library(scales)


# Functions ---------------------------------------------------------------

create_weighted_transition_network <- function(df) {
  # Initialize list to store edges
  edge_list <- list()

  for (i in 1:(nrow(df) - 1)) {
    active_cols_row1 <- names(df)[which(df[i, ] == 1)]
    active_cols_row2 <- names(df)[which(df[i + 1, ] == 1)]
    
    # Exclude columns from row2 that were already active in row1
    new_active_cols_row2 <- setdiff(active_cols_row2, active_cols_row1)
    
    if (length(active_cols_row1) > 0 && length(new_active_cols_row2) > 0) {
      combinations <- expand.grid(active_cols_row1, new_active_cols_row2)
      edge_list[[i]] <- as.matrix(combinations)
    }
  }
  
  # Combine all combinations into one matrix
  all_edges <- do.call(rbind, edge_list)
  edge_df <- as.data.frame(all_edges, stringsAsFactors = FALSE)
  names(edge_df) <- c("from", "to")
  
 # edge_df<-edge_df[edge_df$from!=edge_df$to,]
  
  # Count edge frequencies and normalize by number of rows
  edge_weights <- edge_df %>%
    group_by(from, to) %>%
    summarise(weight = n() / (nrow(df)), .groups = "drop")
  
  # Find code frequencies to determine node sizes
  proportions <- colSums(df == 1, na.rm = TRUE) / nrow(df)
  
  # Create the result data frame
  result <- data.frame(
    Column = names(proportions),
    ProportionOfOnes = proportions,
    row.names = NULL
  )
  
  # Create igraph object
  g <- graph_from_data_frame(edge_weights, directed = TRUE, vertices = colnames(df))
  
  # Assign weights
  E(g)$weight <- edge_weights$weight
  V(g)$freq <- result$ProportionOfOnes
  
  
  V(g)$label.family <- "sans"
  
  return(g)
}

aggregate_transition_network <- function(Instructor) {
  
  dat<-read.csv(paste0(Instructor,"_Class1.csv"))
  
  dat<-dat[9:length(dat$X),c(15:21,23:24)]
  colnames(dat) <- c("Lec","RtW","FUp","PQ","CQ","AnQ","MG","D/V","Adm")
  dat <- dat[dat$L!="L",]
  dat <- dat[dat$L!="1. Students doing",]
  dat<-dat[!apply(dat == "", 1, all),]
  
  dat2<-read.csv(paste0(Instructor,"_Class2.csv"))
  
  dat2<-dat2[9:length(dat2$X),c(15:21,23:24)]
  colnames(dat2) <- c("Lec","RtW","FUp","PQ","CQ","AnQ","MG","D/V","Adm")
  dat2 <- dat2[dat2$L!="L",]
  dat2 <- dat2[dat2$L!="1. Students doing",]
  dat2<-dat2[!apply(dat2 == "", 1, all),]
  
  dat3<-read.csv(paste0(Instructor,"_Class3.csv"))
  
  dat3<-dat3[9:length(dat3$X),c(15:21,23:24)]
  colnames(dat3) <- c("Lec","RtW","FUp","PQ","CQ","AnQ","MG","D/V","Adm")
  dat3 <- dat3[dat3$L!="L",]
  dat3 <- dat3[dat3$L!="1. Students doing",]
  dat3<-dat3[!apply(dat3 == "", 1, all),]
  
  dat_all<-rbind(dat, dat2, dat3)
  
  g<-create_weighted_transition_network(dat_all)
  coords_custom <- layout_in_circle(g, order = c("Lec","RtW","CQ","PQ","MG","FUp","AnQ","D/V","Adm"))
  l <- layout_with_fr(g)
  
  color_palette <- colorRampPalette(c("white", "#cbc9e2"))
  node_colors <- color_palette(100)[as.numeric(cut(V(g)$freq, breaks=100))]
  
  plot(create_weighted_transition_network(dat_all),vertex.color=node_colors,edge.color="gray20", layout=coords_custom,edge.curved=0.1,vertex.size=45,vertex.size2=20,edge.width=(E(create_weighted_transition_network(dat_all))$weight*30),vertex.label.cex=1,vertex.label.color="black",vertex.shape="rectangle",edge.arrow.size=0.5)
  
  return(g)
}

cosine_similarity <- function(g1, g2) {
  # Convert graphs to adjacency matrices (unweighted or weighted)
  adj1 <- get.adjacency(g1, attr = "weight", sparse = FALSE)
  adj2 <- get.adjacency(g2, attr = "weight", sparse = FALSE)
  
  # Flatten the adjacency matrices into vectors
  adj1_vector <- as.vector(adj1)
  adj2_vector <- as.vector(adj2)
  
  # Compute cosine similarity
  dot_product <- sum(adj1_vector * adj2_vector)
  magnitude1 <- sqrt(sum(adj1_vector^2))
  magnitude2 <- sqrt(sum(adj2_vector^2))
  
  # Return the cosine similarity
  return(dot_product / (magnitude1 * magnitude2))
}

weighted_jaccard_directed <- function(g1, g2) {
  
  # Helper function to extract edge weights using (source, target) format
  get_edge_weights <- function(g) {
    edges <- ends(g, E(g), names = TRUE)  # Get source and target node names
    weights <- E(g)$weight
    
    # Create a unique key for each edge as "source_target"
    edge_keys <- apply(edges, 1, paste, collapse = "_")
    names(weights) <- edge_keys
    
    return(weights)
  }
  
  # Extract edge weight mappings
  w1 <- get_edge_weights(g1)
  w2 <- get_edge_weights(g2)
  
  # Union of all edge keys
  all_keys <- union(names(w1), names(w2))
  
  # Align weights across both graphs, fill missing with 0
  w1_full <- w1[all_keys]; w1_full[is.na(w1_full)] <- 0
  w2_full <- w2[all_keys]; w2_full[is.na(w2_full)] <- 0
  
  # Calculate weighted Jaccard components
  intersection <- sum(pmin(w1_full, w2_full))
  union <- sum(pmax(w1_full, w2_full))
  
  if (union == 0) return(0)  # Avoid division by zero
  return(intersection / union)
}

pairwise_function <- function(vec, func, fill_label = "Value") {
  stopifnot(is.vector(vec), is.function(func))
  
  n <- length(vec)
  result_matrix <- matrix(NA, nrow = n, ncol = n)
  
  # Apply the function to each pair
  for (i in 1:n) {
    for (j in 1:n) {
      result_matrix[i, j] <- func(aggregate_transition_network(vec[i]), aggregate_transition_network(vec[j]))
    }
  }
  
  return(result_matrix)
}

combine_upper_lower <- function(mat1, mat2) {
  n <- nrow(mat1)
  result <- matrix(NA, n, n)
  
  # Fill lower triangle (including diagonal) from mat1
  result[lower.tri(result, diag = TRUE)] <- mat1[lower.tri(mat1, diag = TRUE)]
  
  # Fill upper triangle (excluding diagonal) from mat2
  result[upper.tri(result, diag = FALSE)] <- mat2[upper.tri(mat2, diag = FALSE)]
  
  return(result)
}

backboneNetwork<-function(g,alpha,evalFunc){
  #Returns a backbone network based on LANS. 
  A<-get.adjacency(g,attr="weight")
  A<-as.matrix(A)

  p<-A/rowSums(A)
  F_hat<-function(Q){
    x<-vector()
    for(j in 1:length(Q)){
      x[j]<-length(which(Q!=0 & Q<=Q[j]))/length(which(Q>0))
    }
    return(x)
  }
  sigMatrix<-matrix(nrow = length(V(g)), ncol=length(V(g)))
  for(i in 1:length(V(g))){
    sigMatrix[i,]<-F_hat(p[i,])
  }
  sigMatrix2<-sigMatrix >= 1 - alpha
  
  mode(sigMatrix2)<-"numeric"
  sigMatrix2[is.na(sigMatrix2)] <- 0
  B<-sigMatrix2*A
  
  if(evalFunc==1){
    #hard
    rownames(B) <- colnames(B) <- V(g)$name
    h <- graph.adjacency(B, mode = "min", weighted = TRUE)
  }
  else{
    
    #soft
    rownames(B) <- colnames(B) <- V(g)$name
    h <- graph.adjacency(B, mode = "max", weighted = TRUE)
  }
  
  #h<-as.undirected(h, mode = c("collapse"),edge.attr.comb = "min")		
  
  return(h)
}

segregation<-function(memb,x){
  #calculate the weighted cross entropy, the segregation,
  #as described in Bruun and Bearden (2014), DOI: 10.1371/journal.pone.0112775
  df<-data.frame(memb,as.numeric(as.factor(x)))
  par<-as.numeric(table(x))
  tot<-sum(par)
  q<-par/tot
  N<-length(unique(memb)) 
  n<-length(memb) 
  S<-vector()
  for(i in 1:N){
    m<-sum(table(df[which(df[,1]==i),2]))
    p<-table(df[which(df[,1]==i),2])/sum(table(df[which(df[,1]==i),2]))
    pp<-as.numeric(p)
    pq<-q[as.numeric(names(p))]
    weight<-m/n
    S[i]<-weight*sum(pp*log2(pp/pq),na.rm=T)
  }
  S<-S
  return(S)
  
}

entropy<-function(x){
  par<-as.numeric(table(x))
  tot<-sum(par)
  q<-par/tot
  H<--sum(q*log2(q))
  return(H)
}

segMax<-function(memb,x, mode,iter){
  #memb is a membership vector and x is a corresponding attribute vector
  m<-length(unique(memb))
  Seg<-entropy(x)
  df<-data.frame(memb,x)
  
  Frame<-matrix(0,nrow=m,ncol=iter)
  #Frame2<-matrix(0,nrow=m,ncol=iter)
  for (i in 1:iter){
    y<-sample(x,length(memb),replace=F)
    Frame[,i]<-segregation(memb,y)
  } 
  
  if (mode==1){
    tFrame<-as.data.frame(t(Frame))
    mFrame<-sapply(tFrame,mean)
    mFrame<-as.vector(mFrame)
    sdFrame<-sapply(tFrame,sd)
    sdFrame<-as.vector(sdFrame)
    Z<-(Seg-mFrame)/(sdFrame)
    
  }
  else {
    ranMix<-colSums(Frame)
    ranMix<-ranMix/length(x)
    gMix<-sum(Seg)/length(x)
    Z<-(gMix-mean(ranMix))/(sd(ranMix))
  }
  return(Z)
}

resampleX<-function(memb,x,mode,iter){
  #memb is a membership vector and x is a corresponding attribute vector
  m<-length(unique(memb))
  
  Seg<-segregation(memb,x)
  df<-data.frame(memb,x)
  
  Frame<-matrix(0,nrow=m,ncol=iter)
  #Frame2<-matrix(0,nrow=m,ncol=iter)
  for (i in 1:iter){
    y<-sample(x,length(memb),replace=F)
    Frame[,i]<-segregation(memb,y)
  } 
  
  if (mode==1){
    #groupwise segregation Z-scores
    tFrame<-as.data.frame(t(Frame))
    mFrame<-sapply(tFrame,mean)
    mFrame<-as.vector(mFrame)
    sdFrame<-sapply(tFrame,sd)
    sdFrame<-as.vector(sdFrame)
    Z<-(Seg-mFrame)/(sdFrame)
    
  }
  else {
    #total segregation as used in paper
    ranMix<-colSums(Frame)
    ranMix<-ranMix/length(x)
    gMix<-sum(Seg)/length(x)
    Z<-(gMix-mean(ranMix))/(sd(ranMix))
  }
  return(Z)
}

analyze.concept.inventory<- function(rosterfile,method,instructor){
  data<-read.csv(rosterfile)
  n.enrolled<-length(data$ID)
  data<-na.omit(data)
  n.matched.respondents<-length(data$ID)
  percent.matched.respondents<-length(data$ID)/n.enrolled
  df<-(2*n.matched.respondents)-2
  d<-effsize::cohen.d(data$PostScore,data$PreScore,paired=TRUE,hedges.correction = TRUE)$estimate
  d.lower<-effsize::cohen.d(data$PostScore,data$PreScore,paired=TRUE,hedges.correction = TRUE)$conf.int[1]
  d.upper<-effsize::cohen.d(data$PostScore,data$PreScore,paired=TRUE,hedges.correction = TRUE)$conf.int[2]
  firstterm<-(n.matched.respondents*2/(n.matched.respondents^2))+(d^2/(2*df))
  secondterm<-n.matched.respondents*2/df
  var<-firstterm*secondterm
  
  return(c(method,instructor,n.matched.respondents,percent.matched.respondents,d,d.lower,d.upper,var))
}


# Plot all observation networks and import attributes ----------------------------------------------------


# Plot aggregate networks for SCALE-UP instructors
par(mfrow=c(2,4),mar=c(1,1,1,1))
for (i in c(1:7)) {
  aggregate_transition_network(paste0("SCALEUP_Course",i))
}
dev.off()

# Plot aggregate networks for Tutorials instructors
par(mfrow=c(2,4),mar=c(1,1,1,1))
for (i in c(2:9)) {
  aggregate_transition_network(paste0("Tutorials_Course",i))
}
dev.off()

# Plot aggregate networks for Peer Instruction instructors
par(mfrow=c(2,5),mar=c(1,1,1,1))
for (i in c(1:9)) {
  aggregate_transition_network(paste0("PeerInstruction_Course",i))
}
dev.off()

# Plot aggregate networks for ISLE instructors
par(mfrow=c(2,3),mar=c(1,1,1,1))
for (i in c(1:6)) {
  aggregate_transition_network(paste0("ISLE_Course",i))
}
dev.off()

attributes<-read.csv("Course_Attributes.csv")


# Calculate similarities -------------------------------------------


all_classes<-attributes$Instructor


Cosine<-as.data.frame(pairwise_function(all_classes, cosine_similarity, fill_label = "Cosine Similarity"))

diag(Cosine) <- 0

Jaccard<-as.data.frame(pairwise_function(all_classes, weighted_jaccard_directed, fill_label = "Jaccard Index"))

diag(Jaccard) <- 0


# Compare Jaccard and Cosine ----------------------------------------------


Combined<-combine_upper_lower(as.matrix(Jaccard),as.matrix(Cosine))

df <- melt(Combined)
colnames(df) <- c("Row", "Col", "Value")

df$triangle <- with(df, ifelse(Row > Col, "lower",
                               ifelse(Row < Col, "upper", "diag")))

compare<-ggplot() +
  geom_tile(data = subset(df, triangle == "lower"),
            aes(x = Col, y = Row, fill = Value))  +
  scale_fill_gradient2(low = "white", high = "#018571",limits = c(0,1),name="Jaccard \nindex", breaks=c(0,0.5,1)) +
  ggnewscale::new_scale_fill() +geom_tile(data = subset(df, triangle == "upper"),
                                          aes(x = Col, y = Row, fill = Value)) +
  scale_fill_gradient(low = "white", high = "#a6611a",limits = c(0,1),name="Cosine \nsimilarity", breaks=c(0,0.5,1)) +
  scale_y_reverse() +  
  labs(x = "Instructor", y = "Instructor") +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks=element_blank(),axis.title = element_text(size=10,color='black'),axis.text = element_blank(),axis.title.x = element_text(margin = margin(t = -10)),axis.title.y = element_text(margin = margin(r = -10)),legend.text = element_text(size = 7),legend.title=element_text(size = 7),legend.margin = margin(0, 2, 0, -10))+theme(legend.key.height = unit(0.4, "cm"),legend.key.width = unit(0.4, "cm"))

compare

ggarrange(compare) %>%
  ggexport(filename = "CosineVsJaccard.pdf",width=3.5,height=3,units="in")



# Create similarity network ------------------------------------------------

g<-graph_from_adjacency_matrix(as.matrix(Cosine),mode = c("min"),weighted=TRUE)


# Sparsify similarity network ------------------------------------------------


# Choosing alpha level

plot_alpha<-data.frame(alpha = numeric(), edges = numeric(),modularity = numeric(),n_modules = numeric(), giantcomponent = numeric(), duotrio = numeric(),stringsAsFactors = FALSE)

for (a in seq(0,0.2,0.0001)){
  g_sparsified<-backboneNetwork(g,alpha=a,evalFunc=0)
  edges<-length(E(g_sparsified))
  mod<-cluster_infomap(g_sparsified)$mod
  nmod<-length(communities(cluster_infomap(g_sparsified)))
  giant<-max(table(membership(cluster_infomap(g_sparsified))))/ vcount(g_sparsified)
  duotrio_prop<-length(which(membership(cluster_infomap(g_sparsified)) %in% which(sizes(cluster_infomap(g_sparsified)) %in% c(2, 3)))) / vcount(g_sparsified)
  vec<-c(a, edges, mod,nmod,giant,duotrio_prop)
  plot_alpha[nrow(plot_alpha)+1, ] <- vec
}

plot_alpha_long <- pivot_longer(plot_alpha, cols = c("edges", "modularity","n_modules","giantcomponent","duotrio"),
                        names_to = "Measure", values_to = "Value")

mod<-ggplot(plot_alpha, aes(x=alpha, y=modularity)) +
  geom_point(size=2,aes(alpha=0.4)) +theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.y = element_line(color = "black", size = 0.5),axis.ticks.x=element_line(color = "black", size = 0.5),axis.title.y = element_text(size=10,color='black'),axis.title.x=element_blank(),axis.text = element_text(size=9,color='black'),axis.text.x = element_blank()) + theme(strip.background = element_rect(colour='black', size=1),strip.text=element_text(size=9,color='black'))+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + theme(legend.position = "none",legend.title=element_blank(),legend.text=element_text(size=9,color='black'),axis.title.y = element_text(face = "italic"))+ylab("\nQ")+scale_y_continuous(labels = label_number(accuracy = 0.1))

n_mod<-ggplot(plot_alpha, aes(x=alpha, y=n_modules)) +
  geom_point(size=2,aes(alpha=0.4)) +theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.y = element_line(color = "black", size = 0.5),axis.ticks.x=element_line(color = "black", size = 0.5),axis.title.y = element_text(size=10,color='black'),axis.title.x=element_blank(),axis.text = element_text(size=9,color='black'),axis.text.x = element_blank()) + theme(strip.background = element_rect(colour='black', size=1),strip.text=element_text(size=9,color='black'))+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + theme(legend.position = "none",legend.title=element_blank(),legend.text=element_text(size=9,color='black'))+ylab("\nNumber of clusters")+scale_y_continuous(labels = label_number(accuracy = 0.1))

f_dt<-ggplot(plot_alpha, aes(x=alpha, y=duotrio)) +
  geom_point(size=2,aes(alpha=0.4)) +theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.y = element_line(color = "black", size = 0.5),axis.ticks.x=element_line(color = "black", size = 0.5),axis.title.y = element_text(size=10,color='black'),axis.title.x=element_blank(),axis.text = element_text(size=9,color='black'),axis.text.x = element_blank()) + theme(strip.background = element_rect(colour='black', size=1),strip.text=element_text(size=9,color='black'))+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + theme(legend.position = "none",legend.title=element_blank(),legend.text=element_text(size=9,color='black'))+ylab("Fraction of nodes\n in duos or trios")+
  scale_y_continuous(breaks = c(0, 0.1, 0.2),labels = c(0.0, 0.1, 0.2))


giant<-ggplot(plot_alpha, aes(x=alpha, y=giantcomponent)) +
  geom_point(size=2,aes(alpha=0.4)) +theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.y = element_line(color = "black", size = 0.5),axis.ticks.x=element_line(color = "black", size = 0.5),axis.title = element_text(size=10,color='black'),axis.text = element_text(size=9,color='black'),axis.text.x = element_text(size=9,color='black',angle=0,hjust=0.5,vjust=0.15)) + theme(strip.background = element_rect(colour='black', size=1),strip.text=element_text(size=9,color='black'))+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + theme(legend.position = "none",legend.title=element_blank(),legend.text=element_text(size=9,color='black'))+ylab("Fraction of nodes\nin largest cluster")+xlab(expression(alpha))+scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4),labels = c(0.1, 0.2, 0.3, 0.4))

ggarrange(grid.arrange(mod, n_mod, f_dt, giant, ncol = 1))

ggarrange(grid.arrange(mod, n_mod, f_dt, giant, ncol = 1)) %>%
  ggexport(filename = "Alpha.pdf",width=3.5,height=6,units="in")


g_sparsified<-backboneNetwork(g,alpha=0.05,evalFunc=0)



# Infomap clustering ------------------------------------------------------


cl <- cluster_infomap(g_sparsified)

V(g_sparsified)$cluster <- membership(cl)

cluster_colors<-c("#eaaf80","#9ad1f2","#56b3e9","#458fba","#dd7e33")


V(g_sparsified)$method <- attributes$Implementation

set.seed(72)

pdf("InfoMap.pdf", width=6, height=10)  
plot(g_sparsified,
     vertex.color = cluster_colors[V(g_sparsified)$cluster],
     vertex.label = V(g_sparsified)$method,
     vertex.label.color = "black",
     vertex.size = 14,
     edge.width = E(g_sparsified)$weight*2,
     vertex.shape="circle",
     edge.color="black",
     vertex.label.family = "sans",
     vertex.label.cex=0.48,
     layout = layout_with_fr)
dev.off()

# Assign instructors to clusters and plot cluster characteristics ------------------------------------------

# Plot node centralities by cluster

assignment<-data.frame(Instructor=all_classes,Cluster=cluster_infomap(g_sparsified)$membership)

network_characteristics<-data.frame()

for (i in all_classes) {
  net<- aggregate_transition_network(i)
 # instrengths<-as.data.frame(strength(net, mode="in"))
 #  outstrengths<-as.data.frame(strength(net, mode="out"))
  indegrees<-as.data.frame(degree(net, mode="in"))
  outdegrees<-as.data.frame(degree(net, mode="out"))
  betweenness<-as.data.frame(igraph::betweenness(net, directed=TRUE))
  combine<-cbind(indegrees,outdegrees,betweenness)
  combine$Code<-rownames(combine)
  rownames(combine) <- NULL
  colnames(combine) <- c("Indegree","Outdegree","Betweenness","Code")
  combine$Instructor<-i
  network_characteristics<-rbind(network_characteristics, combine)
}

network_characteristics <- network_characteristics %>%
  left_join(assignment, by = "Instructor")

network_characteristics_means<-network_characteristics %>% group_by(Cluster,Code) %>% summarise(across(c("Indegree","Outdegree","Betweenness"), mean, na.rm = TRUE))

network_characteristics_means<-pivot_longer(network_characteristics_means,cols=c("Indegree","Outdegree","Betweenness"),names_to="Centrality",values_to="Value")

network_characteristics_means$Centrality<-factor(network_characteristics_means$Centrality,levels=c("Indegree","Outdegree","Betweenness"))

network_characteristics_means$Cluster <- factor(network_characteristics_means$Cluster, levels = c(5,1,2,3,4))

network_characteristics_means$Code <- factor(network_characteristics_means$Code, levels = c("Lec","RtW","CQ","PQ","MG","FUp","AnQ","D/V","Adm"))

centralities<-ggplot(data=network_characteristics_means,aes(x=Code,y=Value,color=as.factor(Cluster)))+geom_point(size=2.5)+geom_line(aes(group=Cluster))+facet_wrap(~Centrality,ncol=1,scales="free_y",strip.position = "right")+theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.y = element_line(color = "black", size = 0.5),axis.ticks.x=element_line(color = "black", size = 0.5),axis.title = element_text(size=10,color='black'),axis.text = element_text(size=9,color='black'),axis.text.x = element_text(size=9,color='black',angle=0)) + theme(strip.background = element_rect(colour='black', size=0.5),strip.text=element_text(size=9,color='black'))+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + theme(legend.position = "bottom",legend.title=element_blank(),legend.text=element_text(size=9,color='black'),legend.margin = margin(t = -10))+scale_color_manual(values=cluster_colors<-c("#dd7e33","#eaaf80","#9ad1f2","#56b3e9","#458fba"), labels=c("Clicker lecture","Dialogic clicker lecture","Dialogic lecture with short groupwork activities","Short groupwork activities","Long groupwork activities"))+ylab("Average value")+xlab("COPUS code")+guides(color = guide_legend(nrow = 5))+ theme(legend.key.spacing.y = unit(-7, "pt"))

centralities

ggsave(centralities,filename = "ClusterCentralities.pdf",width=5,height=4.5,units="in",device = cairo_pdf)

# Network level features

network_features<-data.frame()

for (i in all_classes) {
  net<- aggregate_transition_network(i)
  edges<-ecount(net)
  transitivity<-transitivity(net)
  mu<-length(which(is.mutual(net)))/edges
  diam<-diameter(net)
  total_weight <- sum(E(net)$weight)
  combine<-c(edges, transitivity, mu, diam, total_weight, i)
  network_features<-rbind(network_features, combine)
}

colnames(network_features) <- c("Edges","Transitivity","Mutuality","Diameter","TotalStrength","Instructor")


network_features <- network_features %>%
  left_join(assignment, by = "Instructor")

network_features$Edges<-as.numeric(network_features$Edges)
network_features$Transitivity<-as.numeric(network_features$Transitivity)
network_features$Mutuality<-as.numeric(network_features$Mutuality)
network_features$Diameter<-as.numeric(network_features$Diameter)
network_features$TotalStrength<-as.numeric(network_features$TotalStrength)

network_features_means<-network_features %>% group_by(Cluster) %>% summarise(across(c("Edges","Transitivity","Mutuality","Diameter","TotalStrength"), mean, na.rm = TRUE))

network_features_means


# Plot observation networks by cluster ------------------------------------

par(mfrow=c(2,4),mar=c(1,1,1,1))
for (i in assignment$Instructor[assignment$Cluster==1]) {
  aggregate_transition_network(paste0(i))
}
dev.off()

par(mfrow=c(2,4),mar=c(1,1,1,1))
for (i in assignment$Instructor[assignment$Cluster==2]) {
  aggregate_transition_network(paste0(i))
}
dev.off()

par(mfrow=c(2,4),mar=c(1,1,1,1))
for (i in assignment$Instructor[assignment$Cluster==3]) {
  aggregate_transition_network(paste0(i))
}
dev.off()

par(mfrow=c(2,4),mar=c(1,1,1,1))
for (i in assignment$Instructor[assignment$Cluster==4]) {
  aggregate_transition_network(paste0(i))
}
dev.off()

par(mfrow=c(2,4),mar=c(1,1,1,1))
for (i in assignment$Instructor[assignment$Cluster==5]) {
  aggregate_transition_network(paste0(i))
}
dev.off()


# Choose one per cluster for paper
pdf("Cluster1_Example.pdf", width=6, height=6)  # Set size as needed
aggregate_transition_network("PeerInstruction_Course3")
dev.off()

pdf("Cluster2_Example.pdf", width=6, height=6)  # Set size as needed
aggregate_transition_network("ISLE_Course3")
dev.off()

pdf("Cluster3_Example.pdf", width=6, height=6)  # Set size as needed
aggregate_transition_network("Tutorials_Course2")
dev.off()

pdf("Cluster4_Example.pdf", width=6, height=6)  # Set size as needed
aggregate_transition_network("SCALEUP_Course3")
dev.off()

pdf("Cluster5_Example.pdf", width=6, height=6)  # Set size as needed
aggregate_transition_network("Tutorials_Course4")
dev.off()

# Segregation measures ----------------------------------------------------
attributes<-attributes %>% mutate(Cluster=cluster_infomap(g_sparsified)$membership)

median(attributes$ClassSize)

class_hist<-ggplot(attributes, aes(x=ClassSize)) + 
  geom_histogram(color="black", fill="grey", bins=18) +
  geom_vline(xintercept = 47, color = "black",size=1.5) +
  theme_minimal() +
  labs(x = "Number of enrolled students", y = "Number of courses")+theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.y = element_line(color = "black", size = 0.5),axis.ticks.x=element_line(color = "black", size = 0.5),axis.title = element_text(size=10,color='black'),axis.text = element_text(size=9,color='black'),axis.text.x = element_text(size=9,color='black')) + theme(strip.background = element_rect(colour='black', size=1),strip.text=element_text(size=9,color='black'))+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + theme(legend.position = "bottom",legend.title=element_blank(),legend.text=element_text(size=9,color='black'))+scale_x_continuous(breaks = c(0, 100, 200, 300, 400, 500, 600))+scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14))

class_hist

ggsave(class_hist,filename = "ClassSizeHistogram.pdf",width=3.5,height=2,units="in",device = cairo_pdf)

attributes <- attributes %>% mutate(ClassSize2=ifelse(attributes$ClassSize<45, "Small","Large"))

attributes <- attributes %>% mutate(Research2=ifelse(attributes$Research=="PUI"|attributes$Research=="RCU", "Other","R1/R2"))

Z<-vector()
Z[1]<-resampleX(attributes$Cluster,attributes$LPA_Profile,2,10000)
Z[2]<-resampleX(attributes$Cluster,attributes$ClassSize2,2,10000)
Z[3]<-resampleX(attributes$Cluster,attributes$ALMethod,2,10000)
Z[4]<-resampleX(attributes$Cluster,attributes$Discipline,2,10000)
Z[5]<-resampleX(attributes$Cluster,attributes$PublicPrivate,2,10000)
Z[6]<-resampleX(attributes$Cluster,attributes$Research2,2,10000)
Z[7]<-resampleX(attributes$Cluster,attributes$PhDGranting,2,10000)

x<-as.table(Z)
names(x)<-c("LPA profile","Class size\n(<45 or >45)","Active learning\nmethod","Discipline (physics\nor astronomy)","Public or\nprivate","Research designation\n(R1/R2 or Other)","PhD-granting")
df <- data.frame(
  Variable = factor(names(x), levels = names(x)),  # maintain order
  Zscore = as.numeric(x)
)

# Plot
seg<-ggplot(df, aes(x = Variable, y = Zscore)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.96, ymax = 1.96,
           alpha = 0.5, fill = "grey80") +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "black") +
  theme_minimal() +
  labs(x = "Course or institution attribute", y = "Segregation Z-score") +
  ylim(-2.5, 7)+theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.y = element_line(color = "black", size = 0.5),axis.ticks.x=element_line(color = "black", size = 0.5),axis.title = element_text(size=10,color='black'),axis.text = element_text(size=9,color='black'),axis.text.x = element_text(size=9,color='black',angle=90,hjust=1,vjust=0.5)) + theme(strip.background = element_rect(colour='black', size=1),strip.text=element_text(size=9,color='black'))+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + theme(legend.position = "bottom",legend.title=element_blank(),legend.text=element_text(size=9,color='black'))

seg

ggsave(seg,filename = "Segregation.pdf",width=3.7,height=3.7,units="in",device = cairo_pdf)



tbl <- table(attributes$Cluster, attributes$ClassSize2)
df <- as.data.frame(tbl)
colnames(df) <- c("Cluster", "ClassSize2", "Count")
df_fraction <- df %>%
  group_by(Cluster) %>%
  mutate(Fraction = Count / sum(Count))


df_fraction$Cluster<-factor(df_fraction$Cluster,labels=c("Dialogic clicker\nlecture","Dialogic lecture with\nshort groupwork\nactivities","Short groupwork\nactivities","Long groupwork\nactivities","Clicker lecture"))
df_fraction$Cluster<-factor(df_fraction$Cluster,levels=c("Clicker lecture","Dialogic clicker\nlecture","Dialogic lecture with\nshort groupwork\nactivities","Short groupwork\nactivities","Long groupwork\nactivities"))
df_fraction$ClassSize2<-factor(df_fraction$ClassSize2,labels=c("Large (>45)","Small (<45)"))
df_fraction$ClassSize2<-factor(df_fraction$ClassSize2,levels=c("Small (<45)","Large (>45)"))


seg_class<-ggplot(df_fraction, aes(x = Cluster, y = Fraction, fill = ClassSize2)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Proportion")+xlab("Instruction type") +theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.y = element_line(color = "black", size = 0.5),axis.ticks.x=element_line(color = "black", size = 0.5),axis.title = element_text(size=10,color='black'),axis.text = element_text(size=9,color='black'),axis.text.x = element_text(size=9,color='black',angle=90,hjust=1,vjust=0.5)) + theme(strip.background = element_rect(colour='black', size=1),strip.text=element_text(size=9,color='black'))+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + theme(legend.position = "right",legend.title=element_text(size=8,color='black'),legend.text=element_text(size=7,color='black'),legend.margin = margin(0, 0, 0, 0))+scale_fill_manual(values=c("#80cfba","#009e74"),name="Class size")+guides(fill = guide_legend(keywidth = 0.8, keyheight = 0.8))

seg_class

ggsave(seg_class,filename = "Segregation_ClassSize.pdf",width=2.5,height=3.5,units="in",device = cairo_pdf)
  
  
tbl2 <- table(attributes$Cluster, attributes$LPA_Profile)
df <- as.data.frame(tbl2)
colnames(df) <- c("Cluster", "LPA_Profile", "Count")
df_fraction <- df %>%
  group_by(Cluster) %>%
  mutate(Fraction = Count / sum(Count))


df_fraction$Cluster<-factor(df_fraction$Cluster,labels=c("Dialogic clicker\nlecture","Dialogic lecture with\nshort groupwork\nactivities","Short groupwork\nactivities","Long groupwork\nactivities","Clicker lecture"))
df_fraction$Cluster<-factor(df_fraction$Cluster,levels=c("Clicker lecture","Dialogic clicker\nlecture","Dialogic lecture with\nshort groupwork\nactivities","Short groupwork\nactivities","Long groupwork\nactivities"))

df_fraction$LPA_Profile<-factor(df_fraction$LPA_Profile,labels=c("Lecture","Clickers","Worksheets","Other\ngroupwork"))
df_fraction$LPA_Profile<-factor(df_fraction$LPA_Profile,levels=c("Lecture","Clickers","Worksheets","Other\ngroupwork"))

seg_lpa<-ggplot(df_fraction, aes(x = Cluster, y = Fraction, fill = LPA_Profile)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Proportion")+xlab("Instruction type") +theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.y = element_line(color = "black", size = 0.5),axis.ticks.x=element_line(color = "black", size = 0.5),axis.title = element_text(size=10,color='black'),axis.text = element_text(size=9,color='black'),axis.text.x = element_text(size=9,color='black',angle=90,hjust=1,vjust=0.5)) + theme(strip.background = element_rect(colour='black', size=1),strip.text=element_text(size=9,color='black'))+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + theme(legend.position = "right",legend.text=element_text(size=7,color='black'),legend.title=element_text(size=8,color='black'),axis.text.x = element_blank(),axis.title.x = element_blank(),legend.margin = margin(0, 0, 0, 0))+
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))+scale_fill_manual(values=c("#0072b2","#d55e00","#cac14a","#007757"),name = "LPA profile")+guides(fill = guide_legend(keywidth = 0.8, keyheight = 0.8))

seg_lpa

ggsave(seg_lpa,filename = "Segregation_LPA.pdf",width=3.5,height=4,units="in",device = cairo_pdf)

ggarrange(grid.arrange(seg_lpa, seg_class, ncol = 1,heights=c(1,1.9))) %>%
  ggexport(filename = "Segregation_Distributions.pdf",width=4.5,height=3.7,units="in")


# Student learning --------------------------------------------------------

mat = matrix(ncol = 8, nrow = 0)
hedgesg = data.frame(mat)

hedgesg<-rbind(hedgesg,c("ISLE","1",NA,NA,NA,NA,NA,NA))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_ISLE_Course2.csv","ISLE","2"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_ISLE_Course3.csv","ISLE","3"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_ISLE_Course4.csv","ISLE","4"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_ISLE_Course5.csv","ISLE","5"))
hedgesg<-rbind(hedgesg,c("ISLE","6",NA,NA,NA,NA,NA,NA))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course1.csv","PeerInstruction","1"))
hedgesg<-rbind(hedgesg,c("PeerInstruction","2",NA,NA,NA,NA,NA,NA))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course3.csv","PeerInstruction","3"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course4.csv","PeerInstruction","4"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course5.csv","PeerInstruction","5"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course6.csv","PeerInstruction","6"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course7.csv","PeerInstruction","7"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course8.csv","PeerInstruction","8"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course9.csv","PeerInstruction","9"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_Tutorials_Course2.csv","Tutorials","2"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_Tutorials_Course3.csv","Tutorials","3"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_Tutorials_Course4.csv","Tutorials","4"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_Tutorials_Course5.csv","Tutorials","5"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_Tutorials_Course6.csv","Tutorials","6"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_Tutorials_Course7.csv","Tutorials","7"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_Tutorials_Course8.csv","Tutorials","8"))
hedgesg<-rbind(hedgesg,c("Tutorials","9",NA,NA,NA,NA,NA,NA))
hedgesg<-rbind(hedgesg,c("SCALEUP","1",NA,NA,NA,NA,NA,NA))
hedgesg<-rbind(hedgesg,c("SCALEUP","2",NA,NA,NA,NA,NA,NA))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_SCALEUP_Course3.csv","SCALEUP","3"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_SCALEUP_Course4.csv","SCALEUP","4"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_SCALEUP_Course5.csv","SCALEUP","5"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_SCALEUP_Course6.csv","SCALEUP","6"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_SCALEUP_Course7.csv","SCALEUP","7"))


colnames(hedgesg)<-c("Pedagogy","Instructor","NMatched","Percent.Matched","HedgesG","LowerCI","UpperCI","Variance")

membership_learning<-cbind(hedgesg,attributes)

membership_learning$Cluster<-as.factor(membership_learning$Cluster)


rma(yi = as.numeric(HedgesG), vi = as.numeric(Variance), mods = ~ Cluster, data = membership_learning, method="REML")

pooled_eff <- rma(yi = as.numeric(HedgesG), vi = as.numeric(Variance), mods = ~ Cluster-1, data = membership_learning, method="REML")

pooled_plot<-as.data.frame(pooled_eff$beta,row.names = FALSE)
pooled_plot<-cbind(pooled_plot,pooled_eff$ci.lb,pooled_eff$ci.ub,pooled_eff$zval,pooled_eff$pval)
pooled_plot<-cbind(pooled_plot,c("Dialogic clicker\nlecture","Dialogic lecture with\nshort groupwork\nactivities","Short groupwork\nactivities","Long groupwork\nactivities","Clicker lecture"))
colnames(pooled_plot)<-c("g","g_lower","g_upper","zvalue","pvalue","method")
pooled_plot$method<-factor(pooled_plot$method,levels=c("Clicker lecture","Dialogic clicker\nlecture","Dialogic lecture with\nshort groupwork\nactivities","Short groupwork\nactivities","Long groupwork\nactivities"))



learn<-ggplot() +
  geom_hline(yintercept=0,color = "black")+
  geom_pointrange(data=pooled_plot,aes(x=method, y=g,ymin = g_lower, ymax = g_upper,colour=method),size=0.65)+
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.y = element_line(color = "black", linewidth = 0.5),axis.ticks.x=element_line(color = "black", linewidth = 0.5),axis.title = element_text(size=10,color='black'),axis.text = element_text(size=9,color='black')) + theme(strip.background = element_rect(colour='gray20', size=1),strip.text=element_text(size=9,color='black'))+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + theme(legend.position = "none",legend.title=element_blank(),legend.text=element_text(size=8,color='black'),axis.text.x = element_text(size=9,color='black',angle=90,hjust=1,vjust=0.5))+ylab("Effect size (Hedges' g)")+scale_colour_manual(values=c("#dd7e33","#eaaf80","#9ad1f2","#56b3e9","#458fba"))+xlab("Instruction type")


learn

ggsave(learn,filename = "Cluster_EffectSizes.pdf",width=4,height=3.25,units="in",device = cairo_pdf)

