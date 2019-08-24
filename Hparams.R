#############
## IMPORTS ##
#############

library('NMF')
library('GenomicRanges')
library('GenomeInfoDb')
library('MutationalPatterns')
library('plotrix')
library("ggplot2")
library("BSgenome.Hsapiens.UCSC.hg19")
library("reshape2")
library('knitr')
library("plotly")

############
## UTILS  ##
############

"

The util functions are organized by categories:
  1) Metrics
    - cosine.similarity
    - L2
    - RSS

  2) General plot functions
    - myImagePlot
    - plot_in_3D
    - plot_a_pie

  3) Signature plot functions 
    - plot_profile_1 (intermediate function useful for plot_profile)
    - plot_profile

  4) Comparison with COSMIC
    - assign_to_cosmic
    - plot_compare_to_cosmic

  5) Analysis functions
    - matrix_reconstruction
    - plot_similarity_matrix
    - exposure_to_clusters

  6) Matrix normalization
    - column_normalization
    - row_normalization

  7) Data loading functions
    - vcf_to_mutmat
    - vcfs_to_mutmat
    - get_sigpro_signatures
    - get_sigpro_activities
    - get_sigpro_stats

  8) Bootstrapping functions
    - profile_bootstrap
    - fit_boot
    - plot_boot
    - add_poisson_noise

  9) Phylogeny function
    - 

  10) Fit to signatures
    - 

  11) Penta-nucleotide functions

"


###################################### 1) Metrics ######################################

cosine.similarity<-function(A,B){
  'Compute the cosine similarity between the vectors A and B'
  return(sum(A*B)/sqrt(sum(A^2)*sum(B^2)))
}

L2<-function(A,B){
  'Compute the L2-norm between the vectors A and B'
  return (sqrt(sum((A-B)^2)))
}

RSS=function(distrib1,distrib2){
  'Compute the RSS (residual sum of squares) between the NORMALIZED distributions distrib1 and distrib2'
  rel1=distrib1/sum(distrib1)
  rel2=distrib2/sum(distrib2)
  diff=rel1-rel2
  RSS=sum(diff^2)
  return(RSS)
}


###################################### 2) General plot functions ######################################

myImagePlot <- function(x, ...){
  
  "This function takes a matrix as input and color it.
  You can change ColorRamp and ColorLevels if you want to change the color scale."
  
  min=min(x)
  max=max(x)
  
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0.,0.9,length=256),  # Red
                    seq(0.,0.4,length=256),  # Green
                    seq(0.,0.3,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
}

plot_in_3D=function(data,x,y,z,color,colors){
  
  "
  plot_in_3D is a function to visualize data in 3D.

  Args :
  - data: dataframe
  - x,y,z : data$x -> axis you want to plot
  - color,colors: vectors of colors associated with each point
  
  Returns:
  - plot
  "
  
  return(plot_ly(data=data, x = x, y = y, z = z,color=color,colors = colors))
}

plot_a_pie=function(count,categories){
  
  "
  plot_a_pie creates a pie plot based on the counts of each category
  
  Args :
  - vector of counts (ex: c(2,54,6,78))
  - vector of categories associated with the counts
  
  Returns:
  - plot
  "
  count <- data.frame(group = categories,value = count/sum(count))
  
  bp<- ggplot(count, aes(x="", y=value, fill=group,ylab="Number of mutation / category"))+geom_bar(width = 1, stat = "identity")
  pie <- bp + coord_polar("y", start=0)
  return(pie)
  
}


###################################### 3) Signature plot functions ######################################


plot_profile_1 = function(profile,
                        profile_name = c("profile"),
                        profile_ymax = 0.2,
                        diff_ylim = c(-0.02, 0.02),
                        colors){

  "
  plot_profile plots a signature profile based its vector of 96 channels. 
  /!\ it is able only to plot ONE signature 

  Args :
  - profile: vector of SBS (96 values)
  - profile_name: name of the sample (not necessary)
  - profile_ymax: ymax of the scale
  - colors: if you want to put your own colors. By default the colors used are those of MutationalPatterns
  
  Returns:
  - plot

  Example:
  plot_profile(signature[,1])
  "
  
  # if colors parameter not provided, set to default colors
  if(missing(colors)){colors = c("#2EBAED", "#000000", "#DE1C14","#D4D2D2", "#ADCC54", "#F0D0CE")}
  s1_relative = profile / sum(profile)
  
  x = cbind(s1_relative)
  colnames(x) = c(profile_name)
  
  substitutions = c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')
  index = c(rep(1,1,16), rep(2,1,16), rep(3,1,16),
            rep(4,1,16), rep(5,1,16), rep(6,1,16))
  
  # Context
  C_TRIPLETS=c("ACA", "ACC", "ACG", "ACT","CCA", "CCC", "CCG", "CCT","GCA", "GCC", "GCG", "GCT","TCA", "TCC", "TCG", "TCT")
  T_TRIPLETS=c("ATA", "ATC", "ATG", "ATT","CTA", "CTC", "CTG", "CTT","GTA", "GTC", "GTG", "GTT","TTA", "TTC", "TTG", "TTT")
  context = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))
  
  # Replace mutated base with dot
  substring(context,2,2) = "."
  
  # Construct dataframe for plotting
  df = data.frame(substitution = substitutions[index], context = context)
  rownames(x) = NULL
  df2 = cbind(df, as.data.frame(x))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  
  # These variables will be available at run-time, but not at compile-time.
  # To avoid compiling trouble, we initialize them to NULL.
  value = NULL
  substitution = NULL
  Sample = NULL
  Contribution = NULL
  Signature = NULL
  
  # Add dummy non_visible data points to force y axis limits per facet
  df4 = data.frame(substitution = rep("C>A", 4),
                   context = rep("A.A",4),
                   variable = c(profile_name),
                   value = c(profile_ymax,profile_ymax))
  
  plot = ggplot(data=df3, aes(x=context,
                              y=value,
                              fill=substitution,
                              width=1)) +
    geom_bar(stat="identity",
             position = "identity",
             colour="black", size=.2) +
    geom_point(data = df4, aes(x = context,
                               y = value), alpha = 0) +
    scale_fill_manual(values=colors) +
    facet_grid(variable ~ substitution, scales = "free_y") +
    ylab("Relative contribution") +
    # ylim(-yrange, yrange) +
    # no legend
    guides(fill=FALSE) +
    # white background
    theme_bw() +
    # format text
    theme(axis.title.y=element_text(size=12,vjust=1),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=12),
          axis.text.x=element_text(size=5,angle=90,vjust=0.4),
          strip.text.x=element_text(size=14),
          strip.text.y=element_text(size=14),
          panel.grid.major.x = element_blank(),
          panel.spacing.x = unit(0, "lines"))
  
  return(plot)
}


plot_profile=function(profiles,profile_name="none"){

  "
  Takes a set of profiles and plot the profile for which the name is provided.  
  If the set of profiles contains only one signature, the function will plot this signature.
  
  Args :
  - profiles: set of profiles [96,n_signatures]
  - profile_name: name of the signature to plot /!\ it has to be the name of a column
  
  Returns:
  - the plot of the signature
  "
  
  # Makes sure that the input profiles matrix has the right profile
  
  if(length(profiles)==96){
    return(plot_profile_1(profiles))
  }
  
  if (nrow(profiles)!=96){
    stop("The dimension has to be [nrow=96,ncol=number_of_signatures]")
  }
  
  if(profile_name=="none"){
    if(length(profiles)!=96){
      stop("You forgot to enter profile_name")
    }
  }
  
  if(sum(1*(colnames(profiles)==profile_name))==0){
    stop("There is no column named profile_name")
  }
  
  else{
    ind=which(colnames(profiles)==profile_name)
    return(plot_profile_1(profiles[,profile_name],profile_name = profile_name))
    
  }

  
  
}



###################################### 4) Comparison with COSMIC ######################################



# Assign a set of signatures to to the most similar signature in COSMIC

assign_to_cosmic=function(signatures,cosmic_signatures="none"){
  
  "
  assign_to_cosmic compares a set of signatures to the COSMIC signatures according to the cosine similarity.
  It prints the most similar COSMIC signatures and return the similarity matrix.

  Args :
  - signatures: set of signatures to compare to COSMIC signatures [nrow=96,ncol=number_of_signatures]
  - cosmic_signatures: matrix containing the COSMIC signatures [nrow=96,ncol=number_of_COSMIC_signatures] (if not provided, take the August 2019 set of COSMIC sigs)
  
  Returns:
  - prints the most similar COSMIC signature for every input signature
  - returns the similarity matrix
  
  Example:
  assign_to_cosmic(set_of_extracted_signatures)

  "
  
  if(cosmic_signatures=="none"){
    cosmic_signatures<-read.csv(file="cosmic_signatures.csv",header=TRUE)
    cosmic_signatures = as.matrix(cosmic_signatures[,3:ncol(cosmic_signatures)])
  }
  
  if(length(signatures)==96){
    A=rep(0,ncol(cosmic_signatures))
    for(i in seq(1,ncol(cosmic_signatures))){
      A[i]=cosine.similarity(signatures,cosmic_signatures[,i])
    }

    print(paste("The closest COSMIC signature is ",colnames(cosmic_signatures)[which.max(A)]," with a CS= ",max(A)))
    
    return(A)
  }
  
  A<-matrix(0L,nrow=ncol(signatures),ncol=ncol(cosmic_signatures))
  
  for(i in 1:(ncol(signatures))){
    for(j in (1:ncol(cosmic_signatures))){
      A[i,j]=cosine.similarity(signatures[,i],cosmic_signatures[,j])
    }
  }
  
  rownames(A)=colnames(signatures)
  colnames(A)=colnames(cosmic_signatures)
  
  for (i in 1:(ncol(signatures))){
    print(paste("The closest COSMIC signature for ",colnames(signatures)[i]," is ",colnames(cosmic_signatures)[which.max(A[i,])]," with a CS= ",max(A[i,])))
  }
  
  return(A)
}

plot_compare_to_cosmic=function(signatures,cosmic_signatures="none"){
  
  "
  plot_compare_to_cosmic is equivalent to assign_to_cosmic but the similarity matrix is colored
  "
  
  if(cosmic_signatures=="none"){
    cosmic_signatures<-read.csv(file="cosmic_signatures.csv",header=TRUE)
    cosmic_signatures = as.matrix(cosmic_signatures[,3:ncol(cosmic_signatures)])
  }
  
  A<-matrix(0L,nrow=ncol(signatures),ncol=ncol(cosmic_signatures))
  
  rownames(A)=colnames(signatures)
  colnames(A)=colnames(cosmic_signatures)
  
  for(i in 1:(ncol(signatures))){
    for(j in (1:ncol(cosmic_signatures))){
      A[i,j]=cosine.similarity(signatures[,i],cosmic_signatures[,j])
    }
  }
  
  return(myImagePlot(A**4))
  
}


###################################### 5) Analysis functions ######################################

matrix_reconstruction=function(exposure,signatures){
 
  "
  Compute the matrix product exposure*signatures = reconstructed matrix
  
  Args:
  - exposure: [nrow=n_signatures,ncol=n_samples]
  - signatures: [nrow=96,ncol=nsignatures]
  
  Returns:
  - reconstructed matrix
  
  "
  
  if(nrow(exposure)!=ncol(signatures)){
    stop(print("The size of the matrix are not compatible"))
  }
  
  n_signatures=nrow(exposure)
  n_samples=ncol(exposure)
  n_channels=nrow(signatures)
  
  reconstructed_matrix=matrix(0,nrow=n_channels,ncol=n_samples)
  
  for(i in seq(1,n_samples)){
    reconstructed_sample=rep(0,96)
    for(j in seq(1,n_signatures)){
      reconstructed_sample=reconstructed_sample+exposure[j,i]*signatures[,j]
    }
    reconstructed_matrix[,i]=reconstructed_sample
  }
  
  return(reconstructed_matrix)
}

plot_similarity_matrix=function(signatures_1,signatures_2){
  
  "
  plot_similarity_matrix plots the similarity matrix between signatures_1 and signatures_2 according to the cosine similarity
  
  Args:
  - signatures_1: first set of signatures [nrow=96,ncol=number_of_signatures]
  - signatures_2: second set of signatures [nrow=96,ncol=number_of_signatures]

  Returns:
  - Plot of the colored similarity matrix

  "
  
  cos_sim<-matrix(0L,nrow=ncol(signatures_1),ncol=ncol(signatures_2))
  
  for(i in 1:ncol(signatures_1)){
    for(j in 1:ncol(signatures_2)){
      cos_sim[i,j]=(cosine.similarity(signatures_1[,i],signatures_2[,j]))**4
    }
  }
  
  myImagePlot(cos_sim)
  
}

exposure_to_clusters=function(exposure){
  
  "
  exposure_to_clusters takes the matrix of exposure, performs a HC and plot the obtained clusters given by the HC
  
  Args:
  - exposure: matrix of signature activities by samples  [nrow=number_of_signatures,ncol=number_of_samples]
  
  Returns:
  - Plot the dendrogram, the given clusters (cosine sim) and the contribution of the signatures
  - the permutation to apply in order to gather samples by cluster
  
  "
  
  # Normalization
  exposure=row_normalization(exposure)
  
  # Similarity matrix (according to the cosine similarity)
  cos_sim_mat=matrix(0L,ncol=nrow(exposure),nrow=nrow(exposure))
  
  for(i in seq(1,nrow(exposure))){
    for(j in seq(1,nrow(exposure))){
      cos_sim_mat[i,j]=cosine.similarity(exposure[i,],exposure[j,])
    }
  }
  
  # Hierarchical clustering
  
  dd <- dist(scale(exposure), method = "euclidean")
  hc <- hclust(dd, method = "ward.D2")
  dend <- as.dendrogram(hc)
  dend_data <- dendro_data(dend, type = "rectangle")
  
  # Permutation of the samples in the order given by the HC
  reo=c()
  for(name in dend_data$labels$label){
    reo=c(reo,which(rownames(exposure)==name))
  }
  
  #Plots
  plot(as.phylo(hc), cex = 0.6, label.offset = 0.5)
  myImagePlot(cos_sim_mat[reo,reo]**4,min=0.,max=1)
  myImagePlot(exposure[reo,],min=0,max=1)
  
  return(reo)
  
}


###################################### 6) Matrix normalization ######################################

column_normalization=function(matrix){
  
  "
  column_normalization normalizes the columns of a matrix

  Args:
  - matrix: a matrix
  
  Returns:
  - the input matrix normalized by column
  
  "

  for (i in seq(1,ncol(matrix))){
    matrix[,i]=matrix[,i]/sum(matrix[,i])
  }
  
  return(matrix)
}

row_normalization=function(matrix){
  
  "
  row_normalization normalizes the rows of a matrix
  
  Args:
  - matrix: a matrix
  
  Returns:
  - the input matrix normalized by column
  
  "
  
  for (i in seq(1,nrow(matrix))){
    matrix[i,]=matrix[i,]/sum(matrix[i,])
  }
  
  return(matrix)
}



###################################### 7) Data loading functions ######################################

vcf_to_mutmat=function(PATH,name){
  
  "
  Takes the path of a vcf file and convert it to a count matrix
  
  Args:
  - PATH: string of the path of the vcf file
  - name: string of the name of the vcf file
  
  Returns:
  - count matrix dim=[96,1]
  
  "

  ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
  vcf_file=PATH
  sample=name
  vcf=read_vcfs_as_granges(vcf_file,sample,ref_genome)
  mut_mat <- mut_matrix(vcf_list = vcf, ref_genome = ref_genome)
  return(as.numeric(t(mut_mat)))
}

vcfs_to_mutmat=function(FOLDER_PATH,names){
  
  "
  Takes the path of a vcf file and convert it to a count matrix
  
  Args:
  - FOLDER_PATH: path of the folder in which all the vcfs are gathered
  - names: vector of the names of the vcf files
  
  Returns:
  - count matrix dim=[96,n_vcf_files]
  
  "
  
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
  vcf_files = list.files(path=FOLDER_PATH, full.names = T)
  sample=names
  vcf=read_vcfs_as_granges(vcf_files,sample,ref_genome)
  mut_mat <- mut_matrix(vcf_list = vcf, ref_genome = ref_genome)
  return(mut_mat)
}

get_sigpro_signatures=function(FOLDER_PATH,K,type="96"){
  
  "
  Takes the path of a vcf file and convert it to a count matrix
  
  Args:
  - FOLDER_PATH: folder of the SigProfiler extraction
  - K: the value of K for which we want to get the results
  
  Returns:
  - the signatures from the extraction [96,n_signatures]
  
  "
  
  SName=paste(FOLDER_PATH,"SBS",type,"/All_solutions/SBS",type,"_",K,"_Signatures/SBS",type,"_S",K,"_Signatures.txt",sep="")
  
  signatures=read.delim(SName,header=T)
  rownames(signatures)=signatures[,1]
  signatures=signatures[,2:(1+K)]
  
  return(signatures)
  

}

get_sigpro_activities=function(FOLDER_PATH,K,type="96"){
  
  "
  Takes the path of a vcf file and convert it to a count matrix
  
  Args:
  - FOLDER_PATH: folder of the SigProfiler extraction
  - K: the value of K for which we want to get the results
  
  Returns:
  - the activities of each signatures / sample [n_samples,n_signatures+1] (a column for the sample names)
  
  "
  
  AName=paste(FOLDER_PATH,"SBS",type,"/All_solutions/SBS",type,"_",K,"_Signatures/SBS",type,"_S",K,"_Activities.txt",sep="")
  
  activities=read.delim(AName,header=T)
  rownames(activities)=activities[,1]
  activities=activities[,2:(1+K)]
  
  return(activities)
  
  
}

get_sigpro_stats=function(FOLDER_PATH,K,type="96"){
  
  "
  Takes the path of a vcf file and convert it to a count matrix
  
  Args:
  - FOLDER_PATH: folder of the SigProfiler extraction
  - K: the value of K for which we want to get the results
  
  Returns:
  - the reconstruction stats for each sample [n_samples,...]
  
  "
  
  StName=paste(FOLDER_PATH,"SBS",type,"/All_solutions/SBS",type,"_",K,"_Signatures/SBS",type,"_S",K,"_Samples_stats.txt",sep="")
  stats=read.delim(StName,header=T)
  rownames(stats)=stats[,1]
  stats=stats[,2:ncol(stats)]
  
  return(stats)
  
}



###################################### 8) Bootstrapping functions ######################################

profile_bootstrap=function(profile){
  
  "
  Takes a profile and resample it -> returns the bootstrapped profile
  
  Args:
  - profile: 96-channels profile [96,1]
  
  Returns:
  - the bootstrapped profile dim=[96,1]
  
  " 
  
  profile_norm=profile/sum(profile)
  
  boot=sample(96, sum(profile), prob = profile_norm, replace = T)
  
  new_profile=rep(0,96)
  
  for(i in boot){
    new_profile[i]=new_profile[i]+1
  }
  
  return(new_profile) 
}


fit_boot=function(sample,signatures,sample_name="Signature distribution",nboot=1000){
  
  "
  Takes a sample, resample it nboot times and fit each boostrapped sample with the set of signatures
  
  Args:
  - sample: 96-channels profile
  - signatures: signatures that will be used for the fitting [96,n_sigs]
  - nboot: number of bootstrap resampling
  
  Returns:
  - a list (contribution,reconstructed) of the contribution for each bootstrapped sample
  
  " 
  
  samples=matrix(0L,ncol=nboot,nrow=96)
  
  samples[,1]=sample
  for(i in seq(2,nboot)){
    samples[,i]=profile_bootstrap(sample)
  } 
  
  boot_1=fit_to_signatures(samples,as.matrix(signatures))
  
  d=as.data.frame(t(boot_1$contribution))
  
  return(boot_1)
  
}

plot_boot=function(sample,signatures,sample_name="Signature distribution",nboot=1000){
  
  "
  Takes a sample, resample it nboot times and fit each boostrapped sample with the set of signatures and plot the distribution of activities
  
  Args:
  - sample: 96-channels profile
  - signatures: signatures that will be used for the fitting [96,n_sigs]
  - nboot: number of bootstrap resampling
  
  Returns:
  - the plot of the distribution of the activities for each signature
  
  " 
  
  boot_1=fit_boot(sample=sample,signatures=signatures,sample_name=sample_name,nboot=nboot)
  
  d=as.data.frame(t(boot_1$contribution))
  
  return(ggplot(data=melt(d))+geom_boxplot(aes(x=variable,y=value))+ggtitle(sample_name)+theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  
}

sieve_to_bootstrap=function(sieve,threshold=0.85,samples,signatures,sample_name="Signature distribution",nboot=1000){
  ind_sigs=which(sieve[,colnames(sieve)==sample_name]>threshold)
  sample=samples[,which(colnames(samples)==sample_name)]
  return(plot_boot(sample=sample,signatures[,ind_sigs],sample_name=sample_name,nboot=1000))
}

add_poisson_noise=function(profile){
  for(i in seq(1,length(profile))){
    profile[i]=rpois(1, lambda = profile[i])
  }
  
  return(profile)
}

###################################### 9) Phylogeny ######################################


#PARTIE CALCUL DE COORDONNEES

calculate_coordinates=function(K,signatures,cosine_cutoff=0.5){
  
  # INPUT
  #    K: final number of extraction
  #    signatures: /!\ has to be the list of the signatures for the different values of K : list(K1signatures,K2signatures...)
  #    cosine_cutoff: cutoff value from which we consider that 2 signatures are the same
  #    /!\ The plotting function is very sensitive to this value -> you have to try values between 0.95 and 0.99
  #    Rk: the form of Ksignatures has to be [96,number_of_signatures]
  
  #Table of coordinates
  X=matrix(0L, nrow=K, ncol = K)
  Y=matrix(0L, nrow=K, ncol = K)
  
  #Pas
  #pas=16*c(128,64,32,16,8,4,2)
  
  pas=c()
  for(i in seq(1,K)){
    pas=c(2**i,pas)
  }
  pas=16*pas
  
  
  # Keep the info to create the futur branches
  split_inds=c(1)
  target_inds=matrix(0L, nrow=2, ncol = K)
  assignation_inds=matrix(0L, nrow=K, ncol = K)
  #Init
  target_inds[1,1]=1
  target_inds[2,1]=2
  for(i in seq(1,K)){
    for(j in seq(1,i)){
      X[j,i]=5*(i-1)
    }
  }
  
  Y[1:2,2]=c(pas[2],-pas[2])
  
  
  
  
  for(Ki in seq(2,K-1)){
    
    assignation=rep(0,Ki)
    target=rep(0,Ki+1)
    sigs1=signatures[[Ki]]
    sigs2=signatures[[Ki+1]]
    
    #New line
    maxs=c()
    #
    
    for(i in seq(1,Ki)){
      
      cos=c()
      for(j in seq(1,Ki+1)){
        cos=c(cos,cosine.similarity(sigs1[,i],sigs2[,j]))
      }
      
      if(max(cos)>cosine_cutoff){
        assignation[i]=which.max(cos)
        target[which.max(cos)]=1
      }
      
      # New
      maxs=c(maxs,max(cos))
      #
    }
    
    # OLD VERSION
    #ind_splitting=which(assignation==0)
    #inds_target=which(target==0)
    #
    
    # NEW VERSION: avoid the bug when all the max(cos)>cutoff
    
    print(sum(1*(maxs>cosine_cutoff))==length(maxs))
    
    if(sum(1*(maxs>cosine_cutoff))==length(maxs)){
      assignation[which.min(maxs)]=0
      cos=c()
      for(j in seq(1,Ki+1)){
        cos=c(cos,cosine.similarity(sigs1[,which.min(maxs)],sigs2[,j]))
      }
      
      target[which.max(cos)]=0
    }

    ind_splitting=which(assignation==0)
    inds_target=which(target==0)
    
    #
    
    #Checking
    if(is.na(inds_target[1])){
      inds_target[1]=0
    }
    
    if(is.na(inds_target[2])){
      inds_target[2]=0
    }
    
    print(ind_splitting)
    if(is.na(ind_splitting)){
      ind_splitting=0
    }
    
    
    #For next
    split_inds=c(split_inds,ind_splitting)
    target_inds[1,Ki]=inds_target[1]
    target_inds[2,Ki]=inds_target[2]
    assignation_inds[1:length(assignation),Ki]=assignation
    
    
    # Parcours de assignation pour checker les signatures statiques
    for(ind in seq(1,length(assignation))){
      if(assignation[ind]>0){
        Y[assignation[ind],i+1]=Y[ind,Ki]
      }
      else{
        Y[inds_target[1],Ki+1]=Y[ind,Ki]-pas[Ki+1]
        Y[inds_target[2],Ki+1]=Y[ind,Ki]+pas[Ki+1]
      }
      
    }
  }
  
  return(list(X,Y,split_inds,target_inds,assignation_inds))
}

calculate_paths=function(K,X,Y,split_inds,target_inds,assignation_inds){
  
  paths=matrix(0L, nrow=K, ncol = K)
  
  paths[1,1]=1
  first_new=2
  
  for(i in seq(1,K-1)){
    split_ind=split_inds[i]
    target_ind=target_inds[,i]
    for(j in seq(1,i)){
      if(paths[j,i]==split_ind){
        paths[first_new,]=paths[j,]
        paths[j,i+1]=target_ind[1]
        paths[first_new,i+1]=target_ind[2]
        first_new=first_new+1
      }
      else{
        paths[j,i+1]=assignation_inds[paths[j,i],i]
      }
    }
  }
  
  return(paths)
}

calculates_paths_coordinates=function(K,paths,X,Y){
  X_plot=matrix(0L, nrow=K, ncol = K)
  Y_plot=matrix(0L, nrow=K, ncol = K)
  
  for(i in seq(1,K)){
    for(j in seq(1,K)){
      if(paths[i,j]!=0){
        X_plot[i,j]=X[paths[i,j],j]
        Y_plot[i,j]=Y[paths[i,j],j]
      }
    }
  }
  
  return(list(X_plot,Y_plot))
}

assign_one_to_cosmic=function(profile,cosmic){
  cos=rep(0,ncol(cosmic))
  for(i in seq(1,ncol(cosmic))){
    cos[i]=cosine.similarity(profile,cosmic[,i])
  }
  return(list(colnames(cosmic)[which.max(cos)],max(cos)))
}

convert_paths_to_letters=function(paths){
  paths_sigs=matrix(0L, nrow=nrow(paths), ncol = ncol(paths))
  for(i in seq(1,nrow(paths))){
    for(j in seq(1,ncol(paths))){
      paths_sigs[i,j]=paste("Sig",LETTERS[paths[i,j]],sep=" ")
    }
  }
  return(paths_sigs)
}

convert_paths_to_cosmic=function(paths,signatures,cosmic_signatures){
  paths_cosmic=matrix(0L, nrow=nrow(paths), ncol = ncol(paths))
  for(i in seq(1,ncol(paths))){
    for(j in seq(1,nrow(paths))){
      if(!(paths_cosmic[j,i]!=0)){
        paths_cosmic[j,i]=assign_one_to_cosmic(signatures[[i]][,paths[j,i]],cosmic_signatures)[[1]]
      }
    }
  }
  return(paths_cosmic)
}


plot_phylogeny=function(signatures,
                        K,
                        cosine_cutoff=0.5,
                        cosmic_assignation=F,
                        cosmic_signatures){
  
  ## CALCUL PART
  
  # Calculate the coordinate
  coo=calculate_coordinates(K=K,signatures = signatures,cosine_cutoff=cosine_cutoff)
  
  X=coo[[1]]
  Y=coo[[2]]
  split_inds=coo[[3]]
  target_inds=coo[[4]]
  assignation_inds=coo[[5]]
  
  # Create the matrix of paths: initial inds -> target inds
  paths=calculate_paths(K,X,Y,split_inds,target_inds,assignation_inds)
  
  # Calculate the coordinates that will be plot
  ret=calculates_paths_coordinates(K,paths,X,Y)
  
  X_plot=ret[[1]]
  Y_plot=ret[[2]]
  
  ## PLOT PART
  
  #Annotation
  paths_sigs=convert_paths_to_letters(paths)
  paths_cosmic=convert_paths_to_cosmic(paths,signatures,cosmic_signatures=cosmic_signatures)
  
  x.melted <- melt(as.data.frame(t(X_plot)))
  y.melted <- melt(as.data.frame(t(Y_plot)))
  melt=cbind(x.melted,y.melted[,2],as.vector(t(paths_sigs)),as.vector(t(paths_cosmic)))
  
  # Create the legend
  
  legend_x=c()
  legend_y_min=c()
  legend_y_max=c()
  legend_text=c()

  
  for(i in seq(1,K)){
    legend_x=c(legend_x,5*(i-1))
    legend_y_min=c(legend_y_min,1.2*min(Y_plot))
    legend_y_max=c(legend_y_max,max(Y_plot))
    legend_text=c(legend_text,paste("K=",toString(i)))
  }
  
  # Plot the phylogeny
  
  if(cosmic_assignation){
    pl=ggplot(data=melt)+geom_line(aes(x = value, y = y.melted[, 2], group = variable))+geom_point(aes(x = value, y = y.melted[, 2], group =variable))+geom_label(x=melt[,2],y=melt[,3]+0,label=melt[,5],colour="#661400",fontface = "bold",fill="#ffd6cc")
  }
  
  else{
    pl=ggplot(data=melt)+geom_line(aes(x = value, y = y.melted[, 2], group = variable))+geom_point(aes(x = value, y = y.melted[, 2], group =variable))+geom_label(x=melt[,2],y=melt[,3]+0,label=melt[,4],colour="#661400",fontface = "bold",fill="#ffd6cc")
  }
  
  pl=pl+annotate("text",x=legend_x,y=legend_y_min,label=legend_text,color="#cc0000")+annotate("segment", x = legend_x+2.5, xend = legend_x+2.5, y = legend_y_min, yend = legend_y_max,
                                                                                          colour = "#3399ff",alpha=0.4)
  
  pl=pl+ylim(1.2*c(min(Y_plot)-10,max(Y_plot)+5))+theme(legend.title = element_blank(),panel.background=element_rect(fill="white"),axis.title = element_blank(),axis.text = element_blank())
  
  pl

}


###################################### 10) Perso fit to signatures ######################################

remove_signatures=function(samples,signatures,set_of_signature,cutoff=0.01){
  
  #BASELINE SCORE
  boot_1=fit_to_signatures(samples,as.matrix(signatures[,set_of_signature]))
  baseline=RSS(rowMeans(samples),rowMeans(boot_1$reconstructed))
  
  new_best_score=0
  while(new_best_score<cutoff){
    #This condition is to prevent from a bug -> fitting with no signature (2 can be replaced by 1)
    if(length(set_of_signature)>2){
      #REMOVE STEPS
      RSSs=c()
      for(i in seq(1,length(set_of_signature))){
        boot_1=fit_to_signatures(samples,as.matrix(signatures[,set_of_signature[-i]]))
        RSSs=c(RSSs,RSS(rowMeans(samples),rowMeans(boot_1$reconstructed)))
      }
      
      new_best_score=min((RSSs-baseline)/baseline)
      
      if(new_best_score<cutoff){
        set_of_signature=set_of_signature[-which.min(RSSs)]
      }
    }
    else{
      new_best_score=1
    }
  }
  
  return(set_of_signature)
}

remove_specific_signatures=function(samples,signatures,set_of_signature,cutoff=0.01){
  
  #BASELINE SCORE
  boot_1=fit_to_signatures(samples,as.matrix(signatures[,set_of_signature]))
  baseline=RSS(rowMeans(samples),rowMeans(boot_1$reconstructed))
  
  new_best_score=0
  
  continue=T
  while(continue){
    #This condition is to prevent from a bug -> fitting with no signature (2 can be replaced by 1)
    continue=F
    if(length(set_of_signature)>2){
      #REMOVE STEPS
      RSSs=c()
      stds=c()
      for(i in seq(1,length(set_of_signature))){
        boot_1=fit_to_signatures(samples,as.matrix(signatures[,set_of_signature[-i]]))
        RSSs=c(RSSs,RSS(rowMeans(samples),rowMeans(boot_1$reconstructed)))
        stds=c(stds,std.error(signatures[,set_of_signature[i]]))
      }
      
      ind=which.min(RSSs*stds)
      
      new_best_score=(RSSs[ind]-baseline)/baseline
      
      
      if(new_best_score<cutoff){
        set_of_signature=set_of_signature[-ind]
        continue=T
      }
    }
  }
  
  return(set_of_signature)
}

add_signatures=function(samples,signatures,set_of_signature,cutoff_add=0.2,cutoff_cs=0.01){
  
  #BASELINE SCORE
  boot_1=fit_to_signatures(samples,as.matrix(signatures[,set_of_signature]))
  baseline=RSS(rowMeans(samples),rowMeans(boot_1$reconstructed))
  baseline_cs=cosine.similarity(rowMeans(samples),rowMeans(boot_1$reconstructed))
  
  #initialization of lists
  common_sigs=intersect(set_of_signature,seq(1,ncol(signatures)))
  sigs_to_add=seq(1,ncol(signatures))[-common_sigs]
  
  # Initialization of scores to enter the while loop
  new_best_score=cutoff_add+0.01
  new_best_cs=cutoff_cs+0.001
  
  if(length(sigs_to_add)<1){
    new_best_score=0
  }
  
  while(new_best_score>cutoff_add){
    while(new_best_cs>cutoff_cs){
      #ADDING STEPS
      RSSs=c()
      CSs=c()
      if(length(sigs_to_add)>0){
        for(i in seq(1,length(sigs_to_add))){
          boot_1=fit_to_signatures(samples,as.matrix(signatures[,c(set_of_signature,sigs_to_add[i])]))
          RSSs=c(RSSs,RSS(rowMeans(samples),rowMeans(boot_1$reconstructed)))
          CSs=c(CSs,cosine.similarity(rowMeans(samples),rowMeans(boot_1$reconstructed)))
        }
        new_best_score=abs(min(RSSs)-baseline)/baseline
        new_best_cs=abs(max(CSs)-baseline_cs)/baseline_cs
      }
      
      if(new_best_score>cutoff_add){
        if(new_best_cs>cutoff_cs){
          #print(colnames(signatures)[sigs_to_add[which.min(RSSs)]])
          set_of_signature=c(set_of_signature,sigs_to_add[which.min(RSSs)])
          sigs_to_add=sigs_to_add[-which.min(RSSs)]
          baseline=min(RSSs)
          baseline_cs=min(CSs)
        }
        else{
          new_best_score=0
        }
      }
      else{
        new_best_cs=0
      }
      
      if(length(sigs_to_add)<1){
        new_best_score=0
      }
    }
    
  }
  
  return(set_of_signature)
  
}

remove_small_exposure=function(samples,
                               signatures,
                               set_of_signature,
                               cutoff_exposure=0.05){
  
  
  #Total number of mutations in the sample
  total_num_mut=mean(colSums(samples))
  
  # loop_bool is a boolean that checks if we have to continue removing signatures
  loop_bool=T
  
  while(loop_bool){
    loop_bool=F
    
    if(length(set_of_signature)>1){
      
      # We fit or samples
      fit=fit_to_signatures(samples,signatures[,set_of_signature])
      
      #Compute the mean exposure by signature
      meanExpo=rowMeans(fit$contribution)
      
      # We remove the small activities one by one starting with the signature with the smallest activity
      min=min(meanExpo)
      
      # if the smallest activity accounts for less than cutoff_exposure*the total number of mutations in the sample, we remove it
      
      if(min<total_num_mut*cutoff_exposure){
        set_of_signature=set_of_signature[-which.min(meanExpo)]
        loop_bool=T
      } 
      
    }
    
  }
  
  return(set_of_signature)
  
}

sparse_filter=function(mut_mat,
                       signatures,
                       ind_prior_knowledge,
                       nboot=2,
                       initial_cutoff=0.2,
                       cutoff=0.02,
                       cutoff_add=0.1,
                       cutoff_cs=0.01){
  
  crible=matrix(0L,ncol=ncol(mut_mat),nrow=ncol(signatures))
  
  for(i in seq(1,ncol(mut_mat))){
    #print(paste("Sample",toString(i)))
    
    #Bootstrapping
    "
    sample=mut_mat[,i]
    samples=matrix(0L,ncol=nboot,nrow=96)
    samples[,1]=sample
    for(j in seq(2,nboot)){
    samples[,j]=profile_bootstrap(sample)
    } 
    "
    
    # With Poisson noise
    sample=mut_mat[,i]
    samples=matrix(0L,ncol=nboot,nrow=96)
    for(j in seq(1,nboot)){
      samples[,j]=add_poisson_noise(sample)
    }
    
    ##### FIRST LIST -> WITH FLAT PROFILES
    
    #print("prior knowledge")
    #print(ind_prior_knowledge)
    
    # Adding the most active signatures when we fit with the entire set of signatures
    fit_add=fit_to_signatures(samples,as.matrix(signatures))
    set_of_sig_to_add=as.numeric(which(rowMeans(fit_add$contribution)>initial_cutoff*sum(sample)))
    
    set_of_signature=union(ind_prior_knowledge,set_of_sig_to_add)
    
    #print("Set with obvious signature")
    if(i==164){
      print(set_of_signature)
    }
    
    
    #Removing as many signatures as we can
    #print("Start removing signatures")
    final_set_of_signature=remove_specific_signatures(samples,signatures,set_of_signature,cutoff=200*cutoff)
    #print(final_set_of_signature)
    
    #Adding as many signatures as we can
    #print("Start adding signatures")
    #cutoff=0.2*cutoff
    final_set_of_signature=add_signatures(samples,signatures,final_set_of_signature,cutoff_add=0.5*cutoff_add,cutoff_cs=0.2*cutoff_cs)
    
    if(i==164){
      print(final_set_of_signature)
    }
    
    # ITE 2
    #Removing as many signatures as we can
    #print("Start removing signatures")
    final_set_of_signature=remove_signatures(samples,signatures,final_set_of_signature,cutoff=cutoff)
    
    #Adding as many signatures as we can
    #print("Start adding signatures")
    final_set_of_signature=add_signatures(samples,signatures,final_set_of_signature,cutoff_add=cutoff_add*0.2)
    
    if(i==164){
      print(final_set_of_signature)
    }
    
    #We remove the signatures with small activities
    #print("Start removing signatures with small activities")
    #cutoff_expo=0.05
    final_set_of_signature=remove_small_exposure(samples,signatures,final_set_of_signature,cutoff_exposure=0.05)
    
    crible[final_set_of_signature,i]=1
    
  }
  
  return(crible)
}

perso_fit_to_signatures=function(mut_mat,
                                 signatures,
                                 ind_prior_knowledge,
                                 n_ite=10,
                                 prob_cutoff=0.3,
                                 nboot_train=2,
                                 nboot_fit=100,
                                 initial_cutoff=0.2,
                                 cutoff=0.01,
                                 cutoff_add=0.1){
  
  #Select the set of signatures by sample using several layers 
  crible=sparse_filter(mut_mat,signatures,ind_prior_knowledge,nboot=nboot_train,initial_cutoff,cutoff,cutoff_add)
  
  for(i in seq(1,n_ite-1)){
    print(print(paste("##### Vote",toString(i))))
    crible=crible+sparse_filter(mut_mat,signatures,ind_prior_knowledge,nboot=nboot_train,initial_cutoff,cutoff,cutoff_add)
  }
  
  return(crible/n_ite)
  
  '  
  crible=1*((crible/n_ite)>prob_cutoff)
  
  
  final_expo=0*crible
  deviation_expo=0*crible
  
  
  
  
  # Bootstrap version
  for(i in seq(1,ncol(mut_mat))){
  ind_sg=which(crible[,i]==1)
  fboot=fit_boot(mut_mat[,i],signatures[,ind_sg],sample_name=colnames(mut_mat[i]),nboot=nboot_fit)
  
  sample_expo=rowMeans(fboot$contribution)
  
  sample_std=rep(0,nrow(fboot$contribution))
  for(k in seq(1,nrow(fboot$contribution))){
  sample_std[k]=std.error(fboot$contribution[k,])
  }
  
  for(j in seq(1,length(sample_expo))){
  final_expo[ind_sg[j],i]=sample_expo[j]
  deviation_expo[ind_sg[j],i]=sample_std[j]
  }
  }
  
  if(length(colnames(mut_mat))==ncol(final_expo)){
  colnames(final_expo)=colnames(mut_mat)
  }
  
  if(length(colnames(signatures))==nrow(final_expo)){
  rownames(final_expo)=colnames(signatures)
  }
  
  
  return(list(final_expo,deviation_expo))
  '
}

sieve_to_exposure=function(samples,signatures,sieve,threshold=0.85){
  sieve=(1*(sieve>threshold))
  
  final_expo=0*sieve
  
  # Bootstrap version
  for(i in seq(1,ncol(samples))){
    
    ind_sg=which(sieve[,i]==1)
    
    if(length(ind_sg)==1){
      expo=fit_to_signatures(cbind(samples[,i],samples[,i]),cbind(signatures[,ind_sg],rep(0,96)))
    }
    
    else{
      expo=fit_to_signatures(cbind(samples[,i],samples[,i]),signatures[,ind_sg])
    }
    
    for(j in seq(1,length(ind_sg))){
      final_expo[ind_sg[j],i]=expo$contribution[j,1]
    }
  }
  
  if(length(colnames(samples))==ncol(final_expo)){
    colnames(final_expo)=colnames(samples)
  }
  
  if(length(colnames(signatures))==nrow(final_expo)){
    rownames(final_expo)=colnames(signatures)
  }
  
  return(final_expo)
  
}



###################################### 11) Penta-nucleotide functions ######################################

penta_to_96=function(penta_profiles){

  mut_types=c("A\\[C>A\\]A","A\\[C>A\\]C","A\\[C>A\\]G","A\\[C>A\\]T","C\\[C>A\\]A","C\\[C>A\\]C","C\\[C>A\\]G","C\\[C>A\\]T",
              "G\\[C>A\\]A","G\\[C>A\\]C","G\\[C>A\\]G","G\\[C>A\\]T","T\\[C>A\\]A","T\\[C>A\\]C","T\\[C>A\\]G","T\\[C>A\\]T",
              "A\\[C>G\\]A","A\\[C>G\\]C","A\\[C>G\\]G","A\\[C>G\\]T","C\\[C>G\\]A","C\\[C>G\\]C","C\\[C>G\\]G","C\\[C>G\\]T",
              "G\\[C>G\\]A","G\\[C>G\\]C","G\\[C>G\\]G","G\\[C>G\\]T","T\\[C>G\\]A","T\\[C>G\\]C","T\\[C>G\\]G","T\\[C>G\\]T",
              "A\\[C>T\\]A","A\\[C>T\\]C","A\\[C>T\\]G","A\\[C>T\\]T","C\\[C>T\\]A","C\\[C>T\\]C","C\\[C>T\\]G","C\\[C>T\\]T",
              "G\\[C>T\\]A","G\\[C>T\\]C","G\\[C>T\\]G","G\\[C>T\\]T","T\\[C>T\\]A","T\\[C>T\\]C","T\\[C>T\\]G","T\\[C>T\\]T",
              "A\\[T>A\\]A","A\\[T>A\\]C","A\\[T>A\\]G","A\\[T>A\\]T","C\\[T>A\\]A","C\\[T>A\\]C","C\\[T>A\\]G","C\\[T>A\\]T",
              "G\\[T>A\\]A","G\\[T>A\\]C","G\\[T>A\\]G","G\\[T>A\\]T","T\\[T>A\\]A","T\\[T>A\\]C","T\\[T>A\\]G","T\\[T>A\\]T",
              "A\\[T>C\\]A","A\\[T>C\\]C","A\\[T>C\\]G","A\\[T>C\\]T","C\\[T>C\\]A","C\\[T>C\\]C","C\\[T>C\\]G","C\\[T>C\\]T",
              "G\\[T>C\\]A","G\\[T>C\\]C","G\\[T>C\\]G","G\\[T>C\\]T","T\\[T>C\\]A","T\\[T>C\\]C","T\\[T>C\\]G","T\\[T>C\\]T",
              "A\\[T>G\\]A","A\\[T>G\\]C","A\\[T>G\\]G","A\\[T>G\\]T","C\\[T>G\\]A","C\\[T>G\\]C","C\\[T>G\\]G","C\\[T>G\\]T",
              "G\\[T>G\\]A","G\\[T>G\\]C","G\\[T>G\\]G","G\\[T>G\\]T","T\\[T>G\\]A","T\\[T>G\\]C","T\\[T>G\\]G","T\\[T>G\\]T")
  
  
  penta_to_96_sigs=matrix(0,ncol=ncol(penta_profiles),nrow=96)
  colnames(penta_to_96_sigs)=colnames(penta_profiles)
  rownames(penta_to_96_sigs)=c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T",
                                 "A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T",
                                 "A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T",
                                 "A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T",
                                 "A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T",
                                 "A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
    
    
  for(i in seq(1,length(mut_types))){
    penta_to_96_sigs[i,]=colSums(penta_profiles[grep(mut_types[i],rownames(penta_profiles)),])
  }
  
  return(penta_to_96_sigs)
}
  
