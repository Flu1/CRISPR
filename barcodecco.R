barcodecco=function(sam,barcodelist)
{
  #####This program will analyse barcodes for Ecco.  It requires a sam file of sequences and the list of barcodes that you are interested in.
#####You will need to switch th matchseq and the refseq to look at ANP32A alleles as opposed to ANP32B alleles
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("ShortRead", version = "3.8")
  BiocManager::install("msa", version = "3.8")
  BiocManager::install("seqinr", version = "3.8")
  library(ShortRead)
  library(msa)
  library(seqinr)
  
  ####Change these bits for different barcodes/sequencing
  
  ###Pick the right seq (A or B) here:
  
  #matchseq="AGTGACGGAGTGACTGA"   #This is the A seq just before the barcode
  matchseq="TGTCTTGGACAATTGCA"   #This is the B seq just before the barcode
 
  #refseq="AGTGACGGAGTGACTGACTGCTTCCTGCGTCTGTCTCCCTGTTCTCTGAGGAGTCCTCTAACACCCTCCCCCTGCTTCAGCCCCATTCCTCAGCCCCGCTCACACGGTCGGGCCTTCTTTCAGGTGAAAGAACTTGTCCTGGACAACAGTCGGTCGAATGAAGGCAAACTCGAAGGCCTCACAGATGAATTTGAAGAACTGGAATTCTTAAGTACAATCAACGTAGGCCTCACCTC" #NNNNNNNNNNGCATNNNNNNTGTCTTGGACAATTGCAAATCAAATGATGGAAAAATTGAGGGCTTAACAGCTGAATTTGTGAACTTAGAGTTCCTCAGTTTAATAAATGTAGGCTTGATCTCAGTTTCAAATCTCCCCAAGCTGCCTAAATTGAAAAAGGTAAGTGCTTTTTCTTTAACAGTAAAAGAGAACGATCCTGGGAAGGGAAAATGTATGATTTTACCTGTAAGGAAGCACTTAGTGTAGCAGAAAGCACATGG"
  refseq="TGTCTTGGACAATTGCAAATCAAATGATGGAAAAATTGAGGGCTTAACAGCTGAATTTGTGAACTTAGAGTTCCTCAGTTTAATAAATGTAGGCTTGATCTCAGTTTCAAATCTCCCCAAGCTGCCTAAATTGAAAAAGGTAAGTGCTTTTTCTTTAACAGTAAAAGAGAACGATCCTGGGAAGGGAAAATGTATGATTTTACCTGTAAGGAAGCACTTAGTGTAGCAGAAAGCACATGG"
  cutoff=0.025  #The cutoff for determining a real error
  
  ##########
  sam1=readLines(sam)
  seqs=array(dim=c(length(sam1)))
  print("This pairs sequences.  Does not sort A and B yet.  It can split barcodes, correct errors and identify alleles.")
  for (a in 4:length(sam1))
  {
    seqs[a]=strsplit(sam1[a],"[\t]")[[1]][10]  #actual read
  }
  seqs=seqs[-c(1:3)]
  
  #Remove bad barcodes if you want here
  seqswithbs=agrep(matchseq,seqs,0.2)     #Which seqs have barcodes?
  barcodes=array(dim=length(seqswithbs))
  seqsclean=barcodes
  #Extract barcodes
  for (a in 1:length(seqswithbs))
  {
    b=aregexec(matchseq,seqs[seqswithbs[a]],0.2)   #Where is the bit just before the barcode
    barcodes[a]=substr(seqs[seqswithbs[a]],b[[1]]-6,b[[1]]-1) #extracts barcode 
    seqsclean[a]=substr(seqs[seqswithbs[a]],b[[1]],nchar(seqs[seqswithbs[a]]))   #extracts sequence without barcode
  }
  print(paste("Number of barcodes is ",length(barcodes)))
    
  
  #How many unique barcodes
  ubarcodes=data.frame(table(barcodes)) 
  print(ubarcodes)
  print(paste("Number of unique barcodes is ",dim(ubarcodes)[1]))
  print(ubarcodes[which(ubarcodes$barcodes%in%barcodelist),])

  
  
  #Which sequence is associated with which barcode
  seqs3=list()
  length(seqs3)=length(barcodelist)
  for(a in 1:length(barcodelist))
  {
    seqs3[[a]]=which(barcodes==barcodelist[a])
  }
  
  for(aa in 1:length(barcodelist))
  {
    seqcutoff=cutoff*length(seqs3[[aa]])    
    seqsdf=data.frame(table(seqsclean[seqs3[[aa]]]))
    names(seqsdf)[1]=c("seqs")
    seqsperb=as.character(seqsdf$seqs)
    seqsperbn=seqsdf$Freq
    print(paste("Finding errors now for barcode",barcodelist[aa],"barcode number",aa,"- There were",length(seqsperb),"unique sequences and",length(seqs3[[aa]]),"total sequences."))
    
    #Find all the errors
    seqef=array(0,dim=length(seqsperb))
    for(a in 1:length(seqsperb))
    {
      if(grepl(seqsperb[a],refseq))  #identify if sequence is error free
      {
        seqef[a]=1
      }
      
    }
    print(paste(sum(seqef),"out of",length(seqef),"sequences were error free",sep=" "))
    
    
    if(sum(seqef)==length(seqef))
    {
      print("No errors.  Program is about to crash.  Sorry.")
    }
    seqef2=which(seqef==0)
    seq3=DNAStringSet(c(refseq,seqsperb[which(seqef==0)]))  #Change strings to a DNA String set for msa (consider adding names?)
    
    errnums=0
    listerrs=array(dim=c(300*length(seqef2),6)) #Adjust this number if it's too big or you run out of space for errors!
    
    for(qq in 1:length(seqef2))
    {
      zz <- file("all.Rout", open = "wt")  #This bit just suppresses the output I hope
      sink(zz)
      sink(zz, type = "message") 
      
      
      seq4=msa(seq3[c(1,qq+1)],order="input")  #Key step.  This does the alignment.
      sink(type = "message")  #Closing all the sinks and connections.  It works occasionally and now I am fretting
      sink()
      closeAllConnections()
      
      seq4=as.matrix(seq4)  #Change the alignment to a matrix
      
      ###Find all the errors!!!!!###
      
      temperrs=which(seq4[2,]!=seq4[1,])
      
      for(b in 1:length(temperrs))
      {
        listerrs[errnums+b,1]=seq4[1,temperrs[b]] #Ref base
        listerrs[errnums+b,2]=seq4[2,temperrs[b]] #Mut base
        listerrs[errnums+b,3]=temperrs[b]         #Base location in alignment NB! Insertions mean can't yet compare between alignments
        listerrs[errnums+b,4]=qq #Which sequence it came from in the alignment
        listerrs[errnums+b,5]= seqsperbn[seqef2[qq]]
        listerrs[errnums+b,6]=paste0(c(listerrs[errnums+b,1:3]),collapse="")
      }
      errnums=errnums+length(temperrs)
      
    }
    listerrs=listerrs[1:errnums,]
    plot(listerrs[,4],listerrs[,3],ylim=c(1, nchar(refseq)),col='red')
    print(paste("There were",dim(listerrs)[1],"errors",sep=" "))
    
    #Identify which errors occur above a threshold
    errsUni=data.frame(table(listerrs[,6]))
    names(errsUni)[1]=c("errs")
    for(a in 1:dim(errsUni)[1])
    {
      errsUni$Freq[a]=sum(as.integer(listerrs[which(listerrs[,6]==errsUni$errs[a]),5])) #calculate the number of times each error occurs
    }
    notrealerrs=which(errsUni$Freq<seqcutoff)  #Errors that occur below the cutoff
    preseqn=length(union(listerrs[,4],listerrs[,4]))  #Finds out how many haplotypes before removing errors
    oldlisterrs=listerrs   #Keeps the pre-error removal list
    if(dim(errsUni)[1]==length(notrealerrs))
    {
      print("All errors below threshhold. Assume all Wildtype")
      next()
    }
    if(length(notrealerrs)>0.5)
    {
      listerrs=listerrs[-which(listerrs[,6]%in%errsUni$errs[notrealerrs]),]
    }
    
    nseqs=data.frame(table(listerrs[,4]))      #How many haplotypes
    names(nseqs)[1]=c("seqn")
    haplotypes=array(dim=c(dim(nseqs)[1],2))
    for(a in 1:dim(nseqs)[1])
    {
      haplotypes[a,1]=paste0(c(listerrs[which(listerrs[,4]==nseqs$seqn[a]),6]),collapse=" ")   #makes the haplotypes
      haplotypes[a,2]=listerrs[which(listerrs[,4]==nseqs$seqn[a])[1],5]
    }
    uhaplotypes=data.frame(table(haplotypes[,1]))   #finds the unique haplotypes
    names(uhaplotypes)[1]=c("haps")
    nhaps=dim(uhaplotypes)[1]
    for(a in 1:nhaps)
    {
      uhaplotypes$Freq[a]=sum(as.integer(haplotypes[which(haplotypes[,1]==uhaplotypes$haps[a]),2]))  #Adjusts the frequency
    }
    points(listerrs[,4],listerrs[,3])
    print(paste("There were",nhaps,"alleles in this sample."))
    
    #This bit just recreates the sequence
    alleles=array(refseq,dim=nhaps)
    allelesn=array(0,dim=nhaps)
    allelesn=as.integer(uhaplotypes$Freq)
    for(a in 1:nhaps)
    {
      temphap=strsplit(as.character(uhaplotypes$haps[a])," ")[[1]]
      for(a2 in 1:length(temphap))
      {
        if(substr(temphap[a2],1,1)!="-")
        {
          substr(alleles[a],substr(temphap[a2],3,nchar(temphap[a2])),substr(temphap[a2],3,nchar(temphap[a2])))=substr(temphap[a2],2,2) #replace things in the reference sequence
        }
        if(substr(temphap[a2],1,1)=="-") #Is it an insertion?
        {
          alleles[a]=paste(c(substr(alleles[a], 1, as.integer(substr(temphap[a2],3,nchar(temphap[a2])))-1), substr(alleles[a], substr(temphap[a2],3,nchar(temphap[a2])),nchar(alleles[a]))), collapse=substr(temphap[a2],2,2)) #This deals with insertions.  May break if there is an insertion at the beginning of the sequence?
        }
      }
    }
  
    allelesabovecutoff=which(allelesn>seqcutoff)
    print(paste("but there were",length(allelesabovecutoff),"alleles above the cutoff of",seqcutoff,"sequences.  Cutoff=",cutoff,"Sequences=",length(seqs3[[aa]])))
    print(alleles[allelesabovecutoff])
    print(allelesn[allelesabovecutoff])
    wildtypeseqs=0
    morewt=setdiff(oldlisterrs[,4],listerrs[,4])
    if (sum(seqef)>0)
    {
      wildtypeseqs=sum(seqsperbn[which(seqef==1)])
      if(length(morewt)==0)
      {
        print(paste("Wildtype alleles were detected -",wildtypeseqs,"total reads"))
        if(wildtypeseqs<seqcutoff)
        {
          print("The number of wildtype reads was below the cutoff and are not likely to be significant.")
        }
      }
    }
    
    wtn=0
    if(length(morewt)>0)
    {
      for(a in 1:length(morewt))
      {
        wtn=wtn+as.integer(oldlisterrs[which(oldlisterrs[,4]==morewt[a])[1],5])
      }
      print(paste("Wildtype alleles detected -",wildtypeseqs+wtn,"total reads"))
      if((wildtypeseqs+wtn)<seqcutoff)
      {
        print("The number of wildtype reads was below the cutoff and are not likely to be significant.")
      }
    }
    dat4filename=paste(strsplit(sam,split="[.]")[[1]][1],barcodelist[aa],"alleles above cutoff.fasta",sep="-")
    
    for (a in 1:length(allelesabovecutoff))
    {
      if(a==1)
      {
        write.fasta(alleles[allelesabovecutoff[a]],paste(strsplit(sam,split="[.]")[[1]][1],barcodelist[aa],a,allelesn[allelesabovecutoff[a]],sep="-"),dat4filename,open = "w", nbchar = 60, as.string = FALSE)
      }
      if(a!=1)
      {
        write.fasta(alleles[allelesabovecutoff[a]],paste(strsplit(sam,split="[.]")[[1]][1],barcodelist[aa],a,allelesn[allelesabovecutoff[a]],sep="-"),dat4filename,open = "a", nbchar = 60, as.string = FALSE)
      }
    }
    
  }
}

