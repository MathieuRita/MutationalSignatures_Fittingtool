add_signatures=function(samples,
                        signatures,
                        set_of_signature,
                        cutoff_add=0.2,
                        cutoff_cs=0.01)
  
  {
  
  "
  Takes the current set of signatures and tests how increased the score is by adding a signature to the current of signatures.
  If adding a signature increases the score more than (1+cutoff_add)*current_score, the signature is added to the current set of signatures
  
  Args:
  - samples: matrix of samples you want to fit [nrow=96,ncol=num_of_samples]
  - signatures: matrix of all the signatures [nrow=96,ncol=num_of_signatures]
  - set_of_signature: vector with the indices of the current signatures
  - cutoff_add: cutoff to choose whether you add a signature to the set of signatures (based on the RSS)
  - cutoff_cs: cutoff to choose whether you add a signature to the set of signatures (based on the cosine similarity)
  
  Returns:
  - set_of_signature: the updated set of signatures
  "
  
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