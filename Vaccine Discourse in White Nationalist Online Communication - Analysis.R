#####################################################################
#  Author: Dror Walter                                              #
#                                                                   #
#  Supplementary code for                                           #
#  Vaccine discourse in white nationalist online communication:     #
#    A mixed-methods computational approach                         #
#  Social Science and Medicine 2022                                 #
#  (https://doi.org/10.1016/j.socscimed.2022.114859)                #
#####################################################################

#### imports ####
options(stringsAsFactors = F)
options(scipen = 999)

library(lmtest)
library(tidytext)
library(tidyverse)
library(cld2)
library(stringi)
library(stringr)
library(wordcloud)
library(ggplot2)
library(lubridate)
library(quanteda)
library(ldatuning)
library(topicmodels)
library(xlsx)
library(textcat)
library(parallel)
library(doParallel)
library(tidyr)
library(tidytext)
library(dplyr)
library(igraph)

#### Data Import ####

searchwords<-c("vaccin")
vacccorpus<-SFdata[grepl(searchwords, tolower(data$text)),]

saveRDS(vacccorpus,"vaccCORP.rds")
vaccorp<-readRDS("vaccCORP.rds")

data<-vaccorp

#### Data Prep #####
# Indexing
data$index<-seq(1,nrow(data))

# Removing short and duplicate texts (to return later)
removed_short<-subset(data,nchar(as.character(data$text))<6) 
data2<-subset(data,!nchar(as.character(data$text))<6)
removed_df<-data2[duplicated(data2$text),]
data3 <- data2[!duplicated(data2$text),]

# Corpus object
mycorpus <- corpus(data3)

# cleaning and removing common/sparse tokens
stopwords_and_single<-c(stopwords("en"),"Originally","originally","posted","Posted",'by',"amp","pic.twitter.com","http", "Quote", "quote","wn","re", LETTERS, letters)
dfm_counts <- dfm(mycorpus,tolower = TRUE, remove_punct = TRUE,remove_numbers=TRUE, 
                  remove = stopwords_and_single,stem = FALSE,
                  remove_separators=TRUE) 
docnames(dfm_counts)<-dfm_counts@docvars$index

dfm_counts2<-dfm_trim(dfm_counts, max_docfreq = 0.90, min_docfreq=0.001,docfreq_type="prop")

#converting to LDA object
dtm_lda <- convert(dfm_counts2, to = "topicmodels",docvars = dfm_counts2@docvars)

full_data<-dtm_lda
rm(dtm_lda)
n <- nrow(full_data)

# quick toy model to check everything is fine before investing in searchk
TOY.lda.100 <- LDA(full_data, k = 100, method = "Gibbs",
                   control = list(alpha=0.2,seed=3522)) 


#### Modelling ####
# 5-fold Cross-validation:

print(Sys.time())
MainresultDF<-data.frame(k=c(1),perplexity=c(1),myalpha=c("x"))
MainresultDF<-MainresultDF[-1,]
candidate_alpha<- c(0.01, 0.05, 0.1, 0.2, 0.5) # we choose variaty of alphas
candidate_k <- c(2,seq(5,100,by=5))
  
for (eachalpha in candidate_alpha) { 
  print ("now running ALPHA:")
  print (eachalpha)
  print(Sys.time())
  #----------------5-fold cross-validation, different numbers of topics----------------
  cluster <- makeCluster(detectCores(logical = TRUE) - 1) # leave one CPU spare...
  registerDoParallel(cluster)
  
  clusterEvalQ(cluster, {
    library(topicmodels)
  })
  
  folds <- 5
  splitfolds <- sample(1:folds, n, replace = TRUE)

  clusterExport(cluster, c("full_data", "splitfolds", "folds", "candidate_k"))
  
  # we parallelize by the different number of topics.  A processor is allocated a value
  # of k, and does the cross-validation serially.  This is because it is assumed there
  # are more candidate values of k than there are cross-validation folds, hence it
  # will be more efficient to parallelise
  system.time({
    results <- foreach(j = 1:length(candidate_k), .combine = rbind) %dopar%{
      k <- candidate_k[j]
      print(k)
      results_1k <- matrix(0, nrow = folds, ncol = 2)
      colnames(results_1k) <- c("k", "perplexity")
      for(i in 1:folds){
        train_set <- full_data[splitfolds != i , ]
        valid_set <- full_data[splitfolds == i, ]
        
        fitted <- LDA(train_set, k = k, method = "Gibbs",
                      #control = list(alpha=eachalpha/k,burnin = burnin, iter = iter, keep = keep) )
                      control = list(alpha=eachalpha) )
        
        results_1k[i,] <- c(k, perplexity(fitted, newdata = valid_set))
      }
      return(results_1k)
    }
  })
  stopCluster(cluster)
  
  results_df <- as.data.frame(results)
  results_df$myalpha<-as.character(eachalpha)
  MainresultDF<-rbind(MainresultDF,results_df)
}

save.image("VACC search k - fullalpha 2-60.RData")

print ("DONE!!!")
print(Sys.time())

# Examining Results
MainresultDF$kalpha=paste0(as.character(MainresultDF$k),MainresultDF$myalpha) 

p<-ggplot(MainresultDF, aes(x = k, y = perplexity))
Alpha<-MainresultDF$myalpha
p+geom_point(aes(color=Alpha),size=0.1)+geom_smooth(se = FALSE, aes(color=Alpha))+
  ggtitle("5-fold cross-validation of topic modelling (5% of data)",
          "(ie five different models fit for each candidate number of topics)") +
  labs(x = "Candidate number of topics", y = "Perplexity when fitting the trained model to the hold-out set")


# Using the elbow method to identify most efficient K value
MainresultDF_MYALPHA<-MainresultDF[MainresultDF$myalpha==0.2,]
cars.spl <- with(MainresultDF_MYALPHA, smooth.spline(k, perplexity, df = 3))
with(cars, predict(cars.spl, x = MainresultDF_MYALPHA$k, deriv = 2))

plot(with(cars, predict(cars.spl, x = MainresultDF_MYALPHA$k, deriv = 2)), type = "l")

# Modeling based on optimal alpha and K
lda.35.02 <- LDA(full_data, k = 35, method = "Gibbs",
                 control = list(alpha=0.2,seed=8375)) 

save.image("final vax with model pilot.RData")


#### ANALYZING the model ####
# extracting excel matrices for topic interpretation
LDAlist<-list(lda.35.02)
datacolnum=1 #text column

for (eachLDA in LDAlist)  {
  LDAfit<-eachLDA
  mybeta<-data.frame(LDAfit@beta)
  colnames(mybeta)<-LDAfit@terms
  mybeta<-t(mybeta)
  colnames(mybeta)<-seq(1:ncol(mybeta))
  mybeta=exp(mybeta)
  
  # Cycle and print top words for each topic
  nwords=50
  
  topwords <- mybeta[1:nwords,]
  for (i in 1:LDAfit@k) {
    tempframe <- mybeta[order(-mybeta[,i]),]
    tempframe <- tempframe[1:nwords,]
    tempvec<-as.vector(rownames(tempframe))
    topwords[,i]<-tempvec
  }
  
  rownames(topwords)<-c(1:nwords)
  
  kalpha<-paste0(as.character(LDAfit@k),"_",gsub("\\.","",as.character(LDAfit@alpha)))
  write.xlsx(topwords, paste0(kalpha,"_ALTVAX_Topwords.xlsx"))
  
  # Cycle and print top FREX (Unique) words for each topic
  # get the beta
  mybeta<-data.frame(LDAfit@beta)
  colnames(mybeta)<-LDAfit@terms
  mybeta<-t(mybeta)
  colnames(mybeta)<-seq(1:ncol(mybeta))
  mybeta=exp(mybeta)
  
  # apply formula below
  # 1/(w/(bword/sumbrow)+(1-w)/(bword)) for each cell
  myw=0.3
  word_beta_sums<-rowSums(mybeta)
  my_beta_for_frex<-mybeta
  for (m in 1:ncol(my_beta_for_frex)) {
    for (n in 1:nrow(my_beta_for_frex)) {
      my_beta_for_frex[n,m]<-1/(myw/(my_beta_for_frex[n,m]/word_beta_sums[n])+((1-myw)/my_beta_for_frex[n,m]))
    }
    print (m)
  }
  # print top 50 frex:
  nwords=50
  
  topwords <- my_beta_for_frex[1:nwords,]
  for (i in 1:LDAfit@k) {
    tempframe <- my_beta_for_frex[order(-my_beta_for_frex[,i]),]
    tempframe <- tempframe[1:nwords,]
    tempvec<-as.vector(rownames(tempframe))
    topwords[,i]<-tempvec
  }
  
  rownames(topwords)<-c(1:nwords)
  
  kalpha<-paste0(as.character(LDAfit@k),"_",gsub("\\.","",as.character(LDAfit@alpha)))
  write.xlsx(topwords,paste0(kalpha,"_ALTVAX.xlsx"))
  
  #Print top texts --->
  data33<-data3
  deleted_lda_texts<-(setdiff(as.character(LDAfit@documents),as.character(data3$index)))
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  data33<-data33[data33$index %!in% deleted_lda_texts,]
  metadf<-data33
  
  print(nrow(metadf))
  print(nrow(LDAfit@gamma))
  
  meta_theta_df<-cbind(metadf[datacolnum],LDAfit@gamma)
  
  ntext=50
  
  toptexts <- mybeta[1:ntext,]
  for (i in 1:LDAfit@k) {
    print(i)
    tempframe <- meta_theta_df[order(-meta_theta_df[,i+1]),]
    tempframe <- tempframe[1:ntext,]
    tempvec<-as.vector(tempframe[,1])
    toptexts[,i]<-tempvec
  }
  
  rownames(toptexts)<-c(1:ntext)
  
  kalpha<-paste0(as.character(LDAfit@k),"_",gsub("\\.","",as.character(LDAfit@alpha)))
  write.xlsx(toptexts, paste0(kalpha,"_ALTVAX.xlsx"))
  
}

#### ANTMN ####

# Assign Topic Names

topic_names<-c("1 Stray Pets", "2 Jews and Media", "3 ZOG and 2nd American Revolution", "4 DELETE: MIXED", 
               "5 Conspiracies - Economy and Technology", "6 Holocaust Denial", "7 Circumcision and Homosexuality", 
               "8 Childrens", "9 Demographic Changes", "10 Genetics and Mutations", "11 Jews and Diseases", "12 Healthcare System", 
               "13 Outbreaks", "14 DELETE: Links of websites","15 Research","16 Geopolitics and Rights", 
               "17 Autism", "18 White Pride and Ancestry", "19 Biological Warfare", "20 Africa", 
               "21 Whites vs Blacks", "22 Toxicity and Side-Effects", "23 Conspiracies Elites", "24 Side Effects Narratives", 
               "25 Black Invention Myths", "26 DELETE: SF Operations", "27 Reaction to Apology to the Black Race", 
               "28 Conspiracies - Pedophilia","29 HPV","30 Flu", "31 White Inventions", "32 WN Constitution", "33 White Inventions", 
               "34 DELETE: MIXED", "35 Food and Health")
# Calculate topic size (for topics-as-nodes)
topicsize<-colMeans(meta_theta_df[,6:(5+LDAfit@k)])

# ANTMN function
network_from_LDA<-function(LDAobject,deleted_topics=c(),topic_names=c(),save_filename="",topic_size=c(),bbone=FALSE) {
  # Importing needed packages
  require(lsa) # for cosine similarity calculation
  require(dplyr) # general utility
  require(igraph) # for graph/network managment and output
  require(corpustools)
  
  print("Importing model")
  
  # first extract the theta matrix form the topicmodel object
  theta<-LDAobject@gamma
  # adding names for culumns based on k
  colnames(theta)<-c(1:LDAobject@k)
  
  # claculate the adjacency matrix using cosine similarity on the theta matrix
  mycosine<-cosine(as.matrix(theta))
  colnames(mycosine)<-colnames(theta)
  rownames(mycosine)<-colnames(theta)
  
  # Convert to network - undirected, weighted, no diagonal
  
  print("Creating graph")
  
  topmodnet<-graph.adjacency(mycosine,mode="undirected",weighted=T,diag=F,add.colnames="label") # Assign colnames
  # add topicnames as name attribute of node - importend from prepare meta data in previous lines
  if (length(topic_names)>0) {
    print("Topic names added")
    V(topmodnet)$name<-topic_names
  } 
  # add sizes if passed to funciton
  if (length(topic_size)>0) {
    print("Topic sizes added")
    V(topmodnet)$topic_size<-topic_size
  }
  newg<-topmodnet
  
  # delete 'garbage' topics
  if (length(deleted_topics)>0) {
    print("Deleting requested topics")
    
    newg<-delete_vertices(topmodnet, deleted_topics)
  }
  
  # Backbone
  if (bbone==TRUE) {
    print("Backboning")
    
    nnodesBASE<-length(V(newg))
    for (bbonelvl in rev(seq(0,1,by=0.05))) {
      #print (bbonelvl)
      nnodes<-length(V(backbone_filter(newg,alpha=bbonelvl)))
      if(nnodes>=nnodesBASE) {
        bbonelvl=bbonelvl
        #  print ("great")
      }
      else{break}
      oldbbone<-bbonelvl
    }
    
    newg<-backbone_filter(newg,alpha=oldbbone)
    
  }
  
  # run community detection and attach as node attribute
  print("Calculating communities")
  
  mylouvain<-(cluster_louvain(newg)) 
  mywalktrap<-(cluster_walktrap(newg)) 
  myspinglass<-(cluster_spinglass(newg)) 
  myfastgreed<-(cluster_fast_greedy(newg)) 
  myeigen<-(cluster_leading_eigen(newg)) 
  
  V(newg)$louvain<-mylouvain$membership 
  V(newg)$walktrap<-mywalktrap$membership 
  V(newg)$spinglass<-myspinglass$membership 
  V(newg)$fastgreed<-myfastgreed$membership 
  V(newg)$eigen<-myeigen$membership 
  
  # if filename is passsed - saving object to graphml object. Can be opened with Gephi.
  if (nchar(save_filename)>0) {
    print("Writing graph")
    write.graph(newg,paste0(save_filename,".graphml"),format="graphml")
  }
  
  # graph is returned as object
  return(newg)
}

# Apply ANTMN function

mynewnet<-network_from_LDA(LDAobject=LDAfit,
                           topic_names=topic_names,
                           topic_size=topicsize,
                           save_filename="altvax_pilotnet_35",
                           bbone=TRUE)

mynewnet<-network_from_LDA(LDAobject=LDAfit,
                           topic_names=topic_names,
                           topic_size=topicsize,
                           deleted_topics=(grep("DELETE:",topic_names)),
                           save_filename="altvax_pilotnet_35_nojunk",
                           bbone=TRUE)


# CREATING META THETA DF for further analysis


LDAfit<-LDAfit

data33<-data3
data33$index<-as.character(data33$index)
deleted_lda_texts<-(setdiff(data33$index, LDAfit@documents))
'%!in%' <- function(x,y)!('%in%'(x,y))
data33<-data33[data33$index %!in% deleted_lda_texts,]
metadf<-data33
meta_theta_df<-cbind(metadf,LDAfit@gamma)
removed_df2<-inner_join(removed_df,meta_theta_df,by="text")
removed_df2<-removed_df2[,-c(5:7)]
colnames(removed_df2)<-gsub("\\.x","",colnames(removed_df2))
removed_df2$index<-as.character(removed_df2$index)
meta_theta_df2<-bind_rows(meta_theta_df,removed_df2)
meta_theta_df<-meta_theta_df2
rm(meta_theta_df2)

colnames(meta_theta_df)[6:(LDAfit@k+5)]<-paste0("X",colnames(meta_theta_df)[6:(LDAfit@k+5)])

meta_theta_df_comm<-meta_theta_df

# Clacuting commulative Theme loadings from Walktrap
meta_theta_df_comm$Conspiracies<-rowSums(meta_theta_df_comm[,(as.numeric(V(mynewnet)$label[which(V(mynewnet)$walktrap == 1)]))+5]) # Conspiracies salmon 19.35
meta_theta_df_comm$Race<-rowSums(meta_theta_df_comm[,(as.numeric(V(mynewnet)$label[which(V(mynewnet)$walktrap == 2)]))+5]) # Race Blue 29.03
meta_theta_df_comm$Science<-rowSums(meta_theta_df_comm[,(as.numeric(V(mynewnet)$label[which(V(mynewnet)$walktrap == 3)]))+5]) # Science gray 41.94
meta_theta_df_comm$White_innovation<-rowSums(meta_theta_df_comm[,(as.numeric(V(mynewnet)$label[which(V(mynewnet)$walktrap == 4)]))+5]) # white inov freen 9.68

# Aggregate by Month for efficient representation
meta_theta_df_comm$dateMONTH <- format(as.Date(meta_theta_df_comm$date), "%Y-%m")



#### Hand Coding ####
# Example Samping:
 set.seed(489374)
 data3_150docs<-sample((1:nrow(data3)),150,replace=FALSE)
 data3_150docs<-data3[data3_150docs,]
 
 xlsx::write.xlsx(data3_150docs,"data3_150docs_PILOT.xlsx")
 
 set.seed(9023)
 data3_100docsFINALREL<-sample((1:nrow(data3)),100,replace=FALSE)
 data3_100docsFINALREL<-data3[data3_100docsFINALREL,]
 
 xlsx::write.xlsx(data3_100docsFINALREL,"data3_100docs_FINALREL.xlsx")
 

# Add hand coded data  
ManualsToImport<-paste0("ManualCoding/",list.files("ManualCoding/"))
rawcodings<-list()
for (file in 1:length(ManualsToImport)) {
  rawcodings[[file]]<-xlsx::read.xlsx(ManualsToImport[file],1)
}

rawcodings<-do.call(rbind,rawcodings)

meta_theta_df_comm_with_hand<-meta_theta_df_comm
meta_theta_df_comm_with_hand$dateYEAR<-str_split(meta_theta_df_comm_with_hand$dateMONTH,"-")%>%  map_chr(c(1))
meta_theta_df_comm_with_hand$indexSF<-rownames(meta_theta_df_comm)

meta_theta_df_comm_with_hand<-dplyr::right_join(meta_theta_df_comm_with_hand,rawcodings,by=c("indexSF"))

meta_theta_df_comm_with_hand$dateYEAR<-str_split(meta_theta_df_comm_with_hand$dateMONTH,"-")%>%  map_chr(c(1))

meta_theta_df_comm_with_hand$Anti_Vaccine_12<-abs(meta_theta_df_comm_with_hand$Anti_Vaccine_12-2)

save.image("altvax after pilot sampling and after manual.RData")


#### Figure 4 pro vs antivac ####
meta_theta_df_comm_with_hand$Anti_Vaccine_12FACT<-factor(meta_theta_df_comm_with_hand$Anti_Vaccine_12,levels=c(0,1),labels=c("Pro/Neutral","Anti"))

myplot <- ggplot(meta_theta_df_comm_with_hand, aes(Anti_Vaccine_12FACT)) + 
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("% of Posts")+
  xlab("Stance Towards Vaccines")+
  theme_bw()+
  ggtitle("Percent of Posts by Vaccination Sentiment")

myplot

#### figure 5 how many coded per year ####
toplot<-meta_theta_df_comm_with_hand %>% dplyr::select("Anti_Vaccine_12",date)

toplot$dateYEAR <- format(as.Date(toplot$date), "%Y")

toplot<-toplot %>% dplyr::group_by(dateYEAR) %>% dplyr::summarise(ANTI=mean(Anti_Vaccine_12),freq=n())
toplot$dateYEAR<-as.Date(toplot$dateYEAR,"%Y")

Panti_time_year<-ggplot(toplot,aes(x=dateYEAR,y=ANTI))+
  geom_col()+
  geom_smooth(method="lm",se=FALSE)+
  geom_text(data=toplot,aes(x=dateYEAR,y=ANTI,label=paste0("n=",freq)),vjust=-0.55)+
  theme(axis.text.x = element_text(angle=90, vjust=0.1, hjust=0,size=10))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")+
  ggtitle("Share of Anti-Vaccination Posts and Number of Coded Posts Per Year")+
  labs(x = "Year", 
       y = "Share of Anti-Vaccination Posts")+
  theme_bw()

#### figure 6 - share of manual categories in content ####

handcat_perc_rank<-data.frame(category=colnames(meta_theta_df_comm_with_hand[,52:63]),perc=colSums(meta_theta_df_comm_with_hand[,52:63])/nrow(meta_theta_df_comm_with_hand))

category_labels<-c("Individual Freedom (1)", "Elite Conspiracies (2)", "Non-White Conspiracies (3)", "NW Spread Diseases (4)", "White Narcissicm (5)", "Efficacy (6)", "Alternative Medicine (7)", 
                   "Vaccine Harms (8)", "Debunking Anti-Vax (9)","Anti-Vax Delegitimize WN (10)","Personal Narratives (11)","Anti-Vaccination (12)")

category_labels<-c("Individual Freedom", "Elite Conspiracies", "Non-White Conspiracies", "NW Spread Diseases", "White Narcissicm", "Efficacy", "Alternative Medicine", 
                   "Vaccine Harms", "Debunking Anti-Vax","Anti-Vax Delegitimize WN","Personal Narratives","Anti-Vaccination")

ggplot(handcat_perc_rank, aes(reorder(category, -perc), perc))+
  geom_col()+ theme_bw()+ 
  ggtitle("Share of Content Categories in Manually Coded Content")+
  labs(x = "Content Categories", 
       y = "Share of Posts")+
  theme(axis.text.x = element_text(angle=90, vjust=0.1, hjust=0,size=10))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")+
  scale_x_discrete(labels = category_labels[order(-handcat_perc_rank$perc)])


#### Figure 7 - spread of antipro for handcoded categories ####
graphlist<-list()
for (i in (1:11)) {
  #print(i)
  mycol<-colnames(meta_theta_df_comm_with_hand)[51+i]
  #print(mycol)
  toplot<-meta_theta_df_comm_with_hand %>% dplyr::select("Anti_Vaccine_12FACT",mycol)
  toplot<-toplot[toplot[,2]==1,]
  toplot<-toplot %>% dplyr::group_by(Anti_Vaccine_12FACT) %>% dplyr::summarise(freq=n())
  toplot$share<-toplot$freq/(sum(toplot$freq))
  if (nrow(toplot)==1) {
    if (toplot$Anti_Vaccine_12FACT[1]=="Anti") {
    toplot[2,]<-toplot[1,]
    toplot[1,1]<-"Pro/Neutral"
    toplot[1,2]<-0
    toplot[1,3]<-0
  } else {
    toplot[2,]<-toplot[1,]
    toplot[2,1]<-"Anti"
    toplot[2,2]<-0
    toplot[2,3]<-0
  }
  } 
  graphlist[[i]]<-ggplot(toplot,aes(x=Anti_Vaccine_12FACT,y=share))+
    geom_col()+ theme_bw()+
    ggtitle(category_labels[i])+
    ylab("# Posts")+
    xlab("Anti Vaccination Stance")
}

# Print ALL
ggpubr::ggarrange(plotlist=graphlist,nrow=2,ncol=6)

# Print only Staple vaccine argimentation (KATA,2010)
ggpubr::ggarrange(graphlist[[1]],graphlist[[2]],graphlist[[6]],graphlist[[7]],graphlist[[8]],graphlist[[9]],graphlist[[11]],nrow=2,ncol=4)

# Print only WN argumentation
ggpubr::ggarrange(graphlist[[3]],graphlist[[4]],graphlist[[5]],graphlist[[10]],nrow=2,ncol=2)

save.image("altvax after pilot sampling and after manual.RData")

#### ANALYZING THEMES AND TOPCICS OVER TIME #### 
post_freq<- meta_theta_df_comm %>% group_by(date) %>% 
  summarise(postfreq=n(),Conspiracies=mean(Conspiracies),Race=mean(Race),Science=mean(Science),White_innovation=mean(White_innovation))

post_freqMONTH<- meta_theta_df_comm %>% group_by(dateMONTH) %>% 
  summarise(postfreq=n(),Conspiracies=mean(Conspiracies),Race=mean(Race),Science=mean(Science),White_innovation=mean(White_innovation))

mybreaks<-c("","","","","2002",
            "","","","","","","","","","","","2003",
            "","","","","","","","","","","","2004",
            "","","","","","","","","","","","2005",
            "","","","","","","","","","","","2006",
            "","","","","","","","","","","","2007",
            "","","","","","","","","","","","2008",
            "","","","","","","","","","","","2009",
            "","","","","","","","","","","","2010",
            "","","","","","","","","","","","2011",
            "","","","","","","","","","","","2012",
            "","","","","","","","","","","","2013",
            "","","","","","","","","","","","2014",
            "","","","","","","","","","","","2015",
            "","","","","","","","","","","","2016",
            "","","","","","","","","","","","2017","","")

ggplot(post_freqMONTH, aes(x=dateMONTH))+geom_col(aes(y=postfreq),width=0.8)+
  theme_bw()+
  ggtitle("Frequency of Vaccine-Related Posts")+
  ylab("Post Count")+
  xlab("Date")+
  geom_smooth(aes(y=postfreq))+
  #scale_x_discrete(labels = NULL)+
  scale_x_discrete(labels = mybreaks)+
  #scale_x_date(date_breaks="Year",expand = c(0, 0), labels = date_format("%Y-%m-%d"))+
  # ylim(0,0.5)+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=0,size=10))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

data_long <- gather(post_freqMONTH, Theme, Loading, Conspiracies:White_innovation, factor_key=FALSE)

theme_month<-ggplot(data_long,aes(x=dateMONTH))+
  #geom_smooth(aes(y=Loading,color=Theme),method="loess",span=0.01,se=FALSE)+
  geom_col(aes(y=Loading,fill=Theme),position = "fill",width = 1)+
  #scale_x_date(date_breaks="day",expand = c(0, 0))+
  xlab("Date")+
  ylab("Theme Loading")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=0,size=10))+
  theme(legend.position="none")+
  ggtitle("Themes in Stormfront Vaccine-Related Discussions (Per Month)")+
  theme(plot.title = element_text(hjust = 0.5))+
  #labs(fill = "THEMES: ")+
  #scale_x_continuous(breaks = NULL)+
  #scale_fill_discrete(name="THEMES: ",
  #                    labels=c("Leisure", "Pakistan", "International (Conflicts)","U.S. Politics","News"))
  scale_fill_manual(values=c("salmon", "skyblue4","gray", "aquamarine4"),
                    name="Themes: ",
                    labels=c("Conspiracies", "Race","Science","White_Innovation"))+
  #theme(legend.position="bottom")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(labels = mybreaks)

# Plot same over years
TEMP<-data_long

TEMP$dateYEAR <- str_split(data_long$dateMONTH,"-")%>%  map_chr(c(1))

theme_year<-ggplot(TEMP,aes(x=dateYEAR))+
  #geom_smooth(aes(y=Loading,color=Theme),method="loess",span=0.01,se=FALSE)+
  geom_col(aes(y=Loading,fill=Theme),position = "fill",width = 1)+
  #scale_x_date(date_breaks="day",expand = c(0, 0))+
  xlab("Date")+
  ylab("Theme Loading")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=0,size=10))+
  theme(legend.position="none")+
  ggtitle("Themes in Stormfront Vaccine-Related Discussions (Per Year)")+
  theme(plot.title = element_text(hjust = 0.5))+
  #labs(fill = "THEMES: ")+
  #scale_x_continuous(breaks = NULL)+
  #scale_fill_discrete(name="THEMES: ",
  #                    labels=c("Leisure", "Pakistan", "International (Conflicts)","U.S. Politics","News"))
  scale_fill_manual(values=c("salmon", "skyblue4","gray", "aquamarine4"),
                    name="Themes: ",
                    labels=c("Conspiracies", "Race","Science","White_Innovation"))+
  theme(legend.position="bottom")+
  scale_y_continuous(expand = c(0,0))


ggpubr::ggarrange(theme_month,theme_year,ncol=1)
  
spansize=0.01

ggplot(post_freq, aes(x=date))+
  geom_smooth(aes(y=Conspiracies), method="loess", color="salmon",se = FALSE, size = 1.5,span=spansize)+      
  geom_smooth(aes(y=Race), method="loess", color="skyblue4",se = FALSE,size = 1.5,span=spansize)+        
  geom_smooth(aes(y=Science), method="loess", color="gray",se = FALSE, size = 1.5,span=spansize)+      
  geom_smooth(aes(y=White_innovation), method="loess", color="aquamarine3",se = FALSE, size = 1.5,span=spansize)+ 
  theme_bw()+ 
  ggtitle("Smoothed Curves for Themes Over Time")+
  labs(x = "\nDate", 
       y = "Theme Prevalence (0-1)\n")+
  scale_x_date(date_breaks="month",expand = c(0, 0))+
  # ylim(0,0.5)+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=0,size=10))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")
