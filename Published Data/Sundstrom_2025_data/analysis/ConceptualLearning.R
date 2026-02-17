library(ggplot2)
library(tidyselect)
library(tidyverse) 
library(dplyr)
library(effsize)
library(metafor)
library(ggpubr)
library(ggsignif)

### Function to calculate Hedge's g for one course
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


# Run function on all 25 courses with sufficient data
mat = matrix(ncol = 8, nrow = 0)
hedgesg = data.frame(mat)

hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_ISLE_Course2.csv","ISLE","2"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_ISLE_Course3.csv","ISLE","3"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_ISLE_Course4.csv","ISLE","4"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_ISLE_Course5.csv","ISLE","5"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course1.csv","PeerInstruction","1"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course3.csv","PeerInstruction","3"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course4.csv","PeerInstruction","4"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course5.csv","PeerInstruction","5"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course6.csv","PeerInstruction","6"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course7.csv","PeerInstruction","7"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course8.csv","PeerInstruction","8"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_PeerInstruction_Course9.csv","PeerInstruction","9"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_SCALEUP_Course3.csv","SCALEUP","3"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_SCALEUP_Course4.csv","SCALEUP","4"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_SCALEUP_Course5.csv","SCALEUP","5"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_SCALEUP_Course6.csv","SCALEUP","6"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_SCALEUP_Course7.csv","SCALEUP","7"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_Tutorials_Course2.csv","Tutorials","2"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_Tutorials_Course3.csv","Tutorials","3"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_Tutorials_Course4.csv","Tutorials","4"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_Tutorials_Course5.csv","Tutorials","5"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_Tutorials_Course6.csv","Tutorials","6"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_Tutorials_Course7.csv","Tutorials","7"))
hedgesg<-rbind(hedgesg,analyze.concept.inventory("Roster_Tutorials_Course8.csv","Tutorials","8"))

colnames(hedgesg)<-c("Method","Course","NMatched","Percent.Matched","HedgesG","LowerCI","UpperCI","Variance")

###  Conduct heterogeneity test

hedgesg$Method<-factor(hedgesg$Method,levels=c("SCALEUP","ISLE","PeerInstruction","Tutorials"))
rma(yi = as.numeric(HedgesG), vi = as.numeric(Variance), mods = ~ Method, data = hedgesg, method="REML")

### Calculate and plot pooled effect sizes for each method

pooled_eff <- rma(yi = as.numeric(HedgesG), vi = as.numeric(Variance), mods = ~ Method-1, data = hedgesg, method="REML")

pooled_plot<-as.data.frame(pooled_eff$beta,row.names = FALSE)
pooled_plot<-cbind(pooled_plot,pooled_eff$ci.lb,pooled_eff$ci.ub,pooled_eff$zval,pooled_eff$pval)
pooled_plot<-cbind(pooled_plot,c("SCALEUP","ISLE","PeerInstruction","Tutorials"))
colnames(pooled_plot)<-c("g","g_lower","g_upper","zvalue","pvalue","method")

pooled_plot$method[pooled_plot$method=="SCALEUP"]<-"SCALE-UP"
pooled_plot$method[pooled_plot$method=="PeerInstruction"]<-"Peer Instruction"

pooled_plot$method<-factor(pooled_plot$method,levels=c("ISLE","Peer Instruction","Tutorials","SCALE-UP"))

SignificanceDF <- data.frame(start = c("ISLE","Peer Instruction"), end = c("SCALE-UP","SCALE-UP"),
                             y = c(2.1,1.8),
                             label = c("2.25\U003C3","2.54\U003C3")) #values are from "zval" column of heterogeneity test above

ci<-ggplot() +
  geom_hline(yintercept=0,color = "black")+
  geom_pointrange(data=pooled_plot,aes(x=method, y=g,ymin = g_lower, ymax = g_upper,colour=method),size=0.65)+
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.y = element_line(color = "black", linewidth = 0.5),axis.ticks.x=element_line(color = "black", linewidth = 0.5),axis.title = element_text(size=9,color='black'),axis.text = element_text(size=8,color='black')) + theme(strip.background = element_rect(colour='gray20', linewidth=1),strip.text=element_text(size=8,color='black'))+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + theme(legend.position = "none",legend.title=element_blank(),legend.text=element_text(size=8,color='black'))+ylab("Effect Size (Hedge's g)")+scale_colour_manual(values=c("#56b3e9","#cc79a7","#e69f00","#009e74"))+xlab("Active Learning Method")+scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+ylim(-0.1,2.2)+geom_signif(
    data = SignificanceDF,
    aes(xmin = start, xmax = end, annotations = label, y_position = y),
    textsize = 2.5, vjust = -0,size=0.2,
    manual = TRUE,orientation="x")

ggsave(ci,filename = "CI.pdf",width=3.5,height=2.5,units="in",device = cairo_pdf)
