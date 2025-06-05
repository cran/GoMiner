# section A contains core GoMiner functions
#'
#' runGoMinerExamples
#' 
#' @import minimalistGODB
#' @import randomGODB
#' @import HGNChelper
#' @import stats
#' @importFrom gplots heatmap.2
#' @import grDevices
#' @import utils
#' @import vprint
#' 
#' @description driver to run GoMiner under several randomization procedures
#' 
#' @param title character string descriptive title
#' @param dir character string full pathname to the directory acting result repository
#' @param sampleList character list of gene names
#' @param GOGOA3 return value of subsetGOGOA()
#' @param ontology character string c("molecular_function", "cellular_component", "biological_process")
#' @param enrichThresh  numerical acceptance threshold for enrichment
#' @param countThresh numerical acceptance threshold for gene count
#' @param pvalThresh numerical acceptance threshold for pval
#' @param fdrThresh numerical acceptance threshold for fdr
#' @param nrand numeric number of randomizations to compute FDR
#' @param mn integer param passed to trimGOGOA3, min size threshold for a category
#' @param mx integer param passed to trimGOGOA3, max size threshold for a category
#' @param verbose integer vector representing classes 
#' 
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
#' ontology<-"biological_process"
#' t<-sort(table(GOGOA3$ontologies[[ontology]][,"HGNC"]),decreasing=TRUE)
#' dir<-tempdir()
#' 
#' sampleList<-names(t)[1:50]
#' title<-"hi_hitters"
#' hh<-runGoMinerExamples(title,dir,sampleList,GOGOA3,ontology,nrand=5)
#' 
#' sampleList<-names(t)[1001:1050]
#' title<-"hi_hitters5"
#' hh<-runGoMinerExamples(title,dir,sampleList,GOGOA3,ontology,nrand=5)
#' 
#' sampleList<-cluster52
#' title<-"cluster52"
#' hh<-runGoMinerExamples(title,dir,sampleList,GOGOA3,ontology,nrand=5)
#' }
#' 
#' @return returns a list containing the return value of GoMiner()
#' 
#' @export
runGoMinerExamples<-
  function(title=NULL,dir,sampleList,GOGOA3,ontology,enrichThresh=2,
           countThresh=5,pvalThresh=0.10,fdrThresh=0.10,nrand=2,mn=2,mx=200,verbose=1) {
    
    stamp<-gsub(":","_",format(Sys.time(), "%a_%b_%d_%Y_%X"))
    if(is.null(title))
      title<-"runGoMinerExamplesResults"
    subd<-sprintf("%s/%s_%s",dir,title,stamp)
    dir.create(subd)
    
    l<-list()
    
    l1<-GoMiner(title=sprintf("%s_1",title),subd,sampleList,GOGOA3,ontology,enrichThresh,
                countThresh,pvalThresh,fdrThresh,nrand,mn,mx,opt=0,verbose)
    
    l2<-GoMiner(title=sprintf("%s_2",title),subd,sampleList,GOGOA3,ontology,enrichThresh,
                countThresh,pvalThresh,fdrThresh,nrand,mn,mx,opt=1,verbose)
    
    l3<-GoMiner(title=sprintf("%s_3",title),subd,sampleList,randomGODB(GOGOA3,FALSE),ontology,enrichThresh,
                countThresh,pvalThresh,fdrThresh,nrand,mn,mx,opt=0,verbose)
    
    l$l1<-l1
    l$l2<-l2
    l$l3<-l3
    
    return(l)
  }

#' GoMiner
#' 
#' @description driver to generate heatmap
#'
#' @param title character string descriptive title
#' @param dir character string full pathname to the directory acting result repository
#' @param sampleList character list of gene names
#' @param GOGOA3 return value of subsetGOGOA()
#' @param ontology character string c("molecular_function", "cellular_component", "biological_process")
#' @param enrichThresh  numerical acceptance threshold for enrichment
#' @param countThresh numerical acceptance threshold for gene count
#' @param pvalThresh numerical acceptance threshold for pval
#' @param fdrThresh numerical acceptance threshold for fdr
#' @param nrand numeric number of randomizations to compute FDR
#' @param mn integer param passed to trimGOGOA3, min size threshold for a category
#' @param mx integer param passed to trimGOGOA3, max size threshold for a category
#' @param opt integer 0:1 parameter used to select randomization method
#' @param verbose integer vector representing classes 
#'
#' @details
#' modes of FDR estimation:
#' opt=0 use original database with randomized geneLists
#' opt=1 use original geneList with internally scrambled genes databases
#'  (uses randomGODB())
#'  
#'  databases that can be used with the real geneList:
#'  these are explicitly passed as parameter to GoMiner()
#'  (1) original GOGOA3
#'  (2) randomized version of GOGOSA GOGOA3R<-randomGODB(GOGOA3)
#'  (3) database containing a subset of the big hitters genes (randomGODB2driver())
#'    attempts to compensate for the over-annotation of some genes,
#'    that might lead to false positive
#'    if gene G has a lot of mappings to categories, randomly sample G/category
#'    pairs to retain a reasonable number of them.
#'    e.g., reduce G from 100 category mappings to 7 category mappings, by omitting
#'    93 of the mappings G/category mappings
#'
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
#' l<-GoMiner("Cluster52",tempdir(),cluster52,
#'  GOGOA3=GOGOA3,ontology="biological_process",enrichThresh=2,
#'  countThresh=5,pvalThresh=0.10,fdrThresh=0.10,nrand=2,mn=2,mx=200,opt=0,verbose=1)
#'  
#'  # try out yeast database!
#'  load("/Users/barryzeeberg/personal/GODB_RDATA/sgd/GOGOA3_sgd.RData")
#'  # make sure this is in fact the database for the desired species
#'  GOGOA3$species
#'  # use database to find genes mapping to an interesting category
#'  cat<-"GO_0042149__cellular_response_to_glucose_starvation"
#'  w<-which(GOGOA3$ontologies[["biological_process"]][,"GO_NAME"]==cat)
#'  geneList<-GOGOA3$ontologies[["biological_process"]][w,"HGNC"]
#'  l<-GoMiner("YEAST",tempdir(),geneList,
#'   GOGOA3,ontology="biological_process",enrichThresh=2,
#'   countThresh=3,pvalThresh=0.10,fdrThresh=0.10,nrand=2,mn=2,mx=200,opt=0)
#' }
#'
#' @return returns a matrix suitable to generate a heatmap
#' 
#' @export
GoMiner<-
  function(title=NULL,dir,sampleList,GOGOA3,ontology,enrichThresh=2,
           countThresh=5,pvalThresh=0.10,fdrThresh=0.10,nrand=100,mn=2,mx=200,opt,verbose=1) {
    
    if(opt<0 | opt>1)
      stop(print(sprintf("Incorrect opt %d, must be 0 or 1",opt)))
    
    stamp<-gsub(":","_",format(Sys.time(), "%a_%b_%d_%Y_%X"))
    if(is.null(title))
      title<-"GMresults"
    subd<-sprintf("%s/%s_%s",dir,title,stamp)
    dir.create(subd)
    
    args<-sprintf("%s/args.txt",subd)
    write(sprintf("Output Results Directory: %s",dir),file=args,append=FALSE)
    write(sprintf("Ontology: %s",ontology),file=args,append=TRUE)
    write(sprintf("enrichThresh: %f",enrichThresh),file=args,append=TRUE)
    write(sprintf("countThresh: %f",countThresh),file=args,append=TRUE)
    write(sprintf("pvalThresh: %f",pvalThresh),file=args,append=TRUE)
    write(sprintf("fdrThresh: %f",fdrThresh),file=args,append=TRUE)
    write(sprintf("nrand: %d",nrand),file=args,append=TRUE)
    write(sprintf("mn: %d",mn),file=args,append=TRUE)
    write(sprintf("mx: %d",mx),file=args,append=TRUE)
    write(sprintf("opt: %d",opt),file=args,append=TRUE)
    
    write(sprintf("DATABASE:"),file=args,append=TRUE)
    write(sprintf("%s",GOGOA3$ontologies[[ontology]][1:3,]),file=args,append=TRUE)
    
    write(sprintf("sampleList:"),file=args,append=TRUE)
    write(sprintf("%s",sampleList),file=args,append=TRUE)
    
    metadata<-sprintf("%s/metadata.txt",subd)
    
    
    write(sprintf("Descriptive Title is %s",title),file=metadata,append=FALSE)
    message(sprintf("Gominer Results Directory is %s",subd))
    write(sprintf("Gominer Results Directory is %s",subd),file=metadata,append=TRUE)
    
    
    vprint(-1,verbose,"geneListDistHitters before preprocessDB")
    geneListDistHitters(sampleList,GOGOA3,ontology,FALSE)
    pp<-preprocessDB(sampleList,GOGOA3,ontology,mn,mx,thresh=.5,verbose)
    vprint(-1,verbose,c("names pp",names(pp)))
    sampleList<-pp$sampleList
    GOGOA3<-pp$GOGOA3
    vprint(-1,verbose,"geneListDistHitters after preprocessDB")
    geneListDistHitters(sampleList,GOGOA3,ontology,FALSE)
    
    DB<-GOGOA3$ontologies[[ontology]]
    
    tableSample3<-GOtable3(sampleList,DB)
    popList<-GOGOA3$genes[[ontology]]
    tablePop3<-GOtable3(popList,DB)
    
    #x_tableSample3<-tableSample3
    #save(x_tableSample3,file="data/x_tableSample3.RData")
    #x_tablePop3<-tablePop3
    #save(x_tablePop3,file="data/x_tablePop3.RData")
    
    m<-GOenrich3(tableSample3,tablePop3)
    hyper<-GOhypergeometric3(tableSample3,tablePop3)
    sink(sprintf("%s/%s",dir,"sink.txt"),append=TRUE)
    print(c(title,"GOhypergeometric3 real geneList:"),quote=FALSE)
    print(hyper[1:10],quote=FALSE)
    print(c(log10(hyper[1]),log10(hyper[10])),quote=FALSE)
    sink()
    
    plot.ecdf(hyper)
    
    #x_sampleList1<-sampleList
    #save(x_sampleList1,file="data/x_sampleList1.RData")
    
    #x_tablePop31<-tablePop3
    #save(x_tablePop31,file="data/x_tablePop31.RData")
    
    #x_hyper1<-hyper
    #save(x_hyper1,file="data/x_hyper1.RData")
    
    fdr<-FDR(sampleList,tablePop3,hyper,GOGOA3,nrand,ontology,dir,opt)
    save(fdr,file=sprintf("%s/%s",subd,"fdr.RData"))
    
    #x_m<-m
    #save(x_m,file="data/x_m.RData")
    #x_fdr<-fdr
    #save(x_fdr,file="data/x_fdr.RData")
    thresh<-GOthresh(m,fdr$sampleFDR,
                     enrichThresh=enrichThresh,countThresh=countThresh,pvalThresh=pvalThresh,fdrThresh=fdrThresh)
    
    #x_thresh<-thresh
    #save(x_thresh,file="data/x_thresh.RData")
    
    gce<-merge(tableSample3$gce,thresh,by.x="GO_NAME",by.y="Row.names")
    gce<-gce[order(gce[,"FDR"]),]
    
    heatmap<-GOheatmap(sampleList,DB,thresh,fdrThresh,verbose)
    
    if(is.null(heatmap))
      return(NULL)
    sink(sprintf("%s/%s",dir,"sink.txt"),append=TRUE)
    print(c("NUMBER OF SIGNIFICANT ONTOLOGY CATEGORIES: ",nrow(heatmap)),quote=FALSE)
    print(c("NUMBER OF SIGNIFICANT ONTOLOGY GENES: ",ncol(heatmap)),quote=FALSE)
    sink()
    write(sprintf("NUMBER OF SIGNIFICANT ONTOLOGY CATEGORIES: %d",nrow(heatmap)),file=metadata,append=TRUE)
    write(sprintf("NUMBER OF SIGNIFICANT ONTOLOGY GENES: %d",ncol(heatmap)),file=metadata,append=TRUE)
    if(nrow(heatmap)<2)
      return(NULL)
    
    svgWidth<-2*(4.75 + ncol(heatmap)*.0792)*30/20.75
    vprint(-1,verbose,c("SVG WIDTH: ",svgWidth))
    write(sprintf("SVG WIDTH: %f",svgWidth),file=metadata,append=TRUE)
    svgHeight<-2*(1.25 + nrow(heatmap)*.0403)*20/12.00
    vprint(-1,verbose,c("SVG HEIGHT: ",svgHeight))
    write(sprintf("SVG HEIGHT: %f",svgHeight),file=metadata,append=TRUE)
    
    file<-sprintf("%s/GoMiner_%d_%d.svg",subd,nrow(heatmap),ncol(heatmap))
    message(sprintf("SVG Heatmap is stored as file://%s",file))
    write(sprintf("SVG Heatmap is stored as %s",file),file=metadata,append=TRUE)
    
    svg(filename=file,width=svgWidth,height=svgHeight)
    hm<-heatmap.2(heatmap,col = heat.colors(n=100,rev=FALSE),trace="none",lhei=c(1,15),lwid=c(1,15),key=FALSE,margins = c(5, 50))
    #hm<-heatmap.2(heatmap,col = heat.colors(n=100,rev=FALSE),trace="none",
    #      lmat = rbind(c(0,3),c(2,1),c(0,4)),lhei = c(1.5,4,1),lwid = c(1.5,4),
    #     key=TRUE,keysize=1)
    dev.off()
    write.table(gce,file=sprintf("%s/gce_%d_%d.txt",subd,nrow(heatmap),ncol(heatmap)),quote=FALSE,sep="\t",col.names=NA)
    write.table(tableSample3$t,file=sprintf("%s/thresh_%d_%d.txt",subd,nrow(heatmap),ncol(heatmap)),quote=FALSE,sep="\t",col.names=NA)
    write.table(heatmap[rev(hm$rowInd),hm$colInd],file=sprintf("%s/heatmap_%d_%d.txt",subd,nrow(heatmap),ncol(heatmap)),quote=FALSE,sep="\t",col.names=NA)
    
    l<-list()
    l$thresh<-thresh
    l$gce<-gce
    
    return(l)
  }

#' GOtable3
#' 
#' @description tabulate number of geneList mappings to GO categories
#' 
#' @param hgncList character list of gene names
#' @param DB selected ontology branch of return value of subsetGOGOA
#' 
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
#' x<-GOtable3(cluster52,GOGOA3$ontologies[["biological_process"]])
#' }
#'
#' @return returns a list whose components are c("DB","table","ngenes")
#'  where 'DB' is the GO DB subsetted to the desired ONTOLOGY,
#'  and 'table' is tabulation of number of occurrences of each GO
#'  category name within the desired ONTOLOGY,
#'  and ngenes is the total number of hgncList genes mapping to GOGOA 
#' 
#' @export
GOtable3<-
  function(hgncList,DB) {
    l<-list()
    w<-which(DB[,"HGNC"] %in% hgncList)
    t<-table(DB[w,"GO_NAME"])
    l$DB<-DB[w,]
    mappingGenes<-unique(DB[w,"HGNC"])
    l$table<-t
    l$ngenes<-length(mappingGenes)
    l$gce<-DB[w,c("HGNC","GO_NAME")]
    
    return(l)
  }

#' GOenrich3
#' 
#' @description compute the gene enrichment in a GO category
#' 
#' @param tableSample3 sample return value of GOtable3()
#' @param tablePop3 population return value of GOtable3()
#' 
#' @examples
#' m<-GOenrich3(x_tableSample3,x_tablePop3)
#'
#' @return returns a matrix with columns c("SAMPLE","POP","ENRICHMENT")
#'
#' @export
GOenrich3<-
  function(tableSample3,tablePop3) {
    m<-matrix(ncol=3,nrow=length(tableSample3$table))
    colnames(m)<-c("SAMPLE","POP","ENRICHMENT")
    rownames(m)<-names(tableSample3$table)
    
    for(cat in names(tableSample3$table)) {
      m[cat,"SAMPLE"]<-tableSample3$table[cat]
      m[cat,"POP"]<-tablePop3$table[cat]
      m[cat,"ENRICHMENT"]<-(tableSample3$table[cat]/tableSample3$ngenes)/(tablePop3$table[cat]/tablePop3$ngenes)
    }
    return(m)
  }

#' FDR
#' 
#' @description compute the false discovery rate (FDR) of the hypergeometric
#' p values of genes mapping to gene ontology (GO) categories
#'
#' @param sampleList character vector of user-supplied genes of interest
#' @param tablePop3 return value of GOtable3()
#' @param hyper return value of GOhypergeometric3()
#' @param GOGOA3 return value of subsetGOGOA()
#' @param nrand integer number of randomizations
#' @param ontology c("molecular_function","cellular_component","biological_process")
#' @param subd character string pathname for directory containing sink.txt
#' @param opt integer 0:1 parameter used to determine randomization method
#' 
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
#' fdr<-FDR(x_sampleList1,x_tablePop31,x_hyper1,GOGOA3,3,"biological_process",tempdir(),0)
#' }
#' 
#' @return returns a list with FDR information
#' 
#' @export
FDR<-
  function(sampleList,tablePop3,hyper,GOGOA3,nrand,ontology,subd,opt=0) {
    
    
    if(opt<0 | opt >1)
      stop(sprintf("FDR opt must either 0 or 1, not %d",opt))
    
    l<-list()
    ngenes<-length(sampleList)
    
    # rcpd() is a function that interpolates observed p values
    rcpd<-RCPD(GOGOA3,tablePop3,sampleList,nrand,ontology,hyper,subd,opt)
    
    sampleFDR<-matrix(nrow=length(names(hyper)),ncol=2)
    rownames(sampleFDR)<-names(hyper)
    colnames(sampleFDR)<-c("OBSERVED_LOG10(P)","FDR")
    for(cat in names(hyper)) {
      x<-log10(hyper[cat])
      sampleFDR[cat,"OBSERVED_LOG10(P)"]<-x
      if(is.null(rcpd)) # if no randoms met criteria, then just set FDR to 0
        sampleFDR[cat,"FDR"]<-0
      else
        sampleFDR[cat,"FDR"]<-rcpd(x)
    }
    
    l$rcpd<-rcpd
    l$sampleFDR<-sampleFDR
    
    return(l)
  }

#' GOhypergeometric
#' 
#' @description compute the hypergeometric p value for gene enrichment in a GO category
#' 
#' @param tableSample3 sample return value of GOtable3()
#' @param tablePop3 population return value of GOtable3()
#' 
#' @examples
#' hyper<-GOhypergeometric3(x_tableSample3,x_tablePop3)
#'
#' @return returns a matrix with columns c("x","m","n","k","p")
#' 
#' @export
GOhypergeometric3<-
  function(tableSample3,tablePop3) {
    # https://www.geeksforgeeks.org/hypergeometric-distribution-in-r-programming/
    # x: number of genes in the cat in sample
    # m: number of genes in the cat in population
    # n: number of genes not in the category in population
    # k: total number of genes in the sample
    
    m<-matrix(nrow=length(tableSample3$table),ncol=5)
    rownames(m)<-names(tableSample3$table)
    colnames(m)<-c("x","m","n","k","p")
    for(cat in names(tableSample3$table)) {
      m[cat,"x"]<-tableSample3$table[cat]
      m[cat,"m"]<-tablePop3$table[cat]
      m[cat,"n"]<-tablePop3$ngenes-tablePop3$table[cat]
      m[cat,"k"]<-tableSample3$ngenes
      m[cat,"p"]<-dhyper(m[cat,"x"],m[cat,"m"],m[cat,"n"],m[cat,"k"])
    }
    mm<-m[order(m[,"p"]),"p"]
    return(m[order(m[,"p"]),"p"])
  }

#' GOthresh
#' 
#' @description retrieve lines of m that meet both enrichThresh and countThresh
#' 
#' @param m return value of GOenrich3()
#' @param sampleFDR component of return value of RCPD()
#' @param enrichThresh  numerical acceptance threshold for enrichment
#' @param countThresh numerical acceptance threshold for gene count
#' @param pvalThresh numerical acceptance threshold for pval
#' @param fdrThresh numerical acceptance threshold for fdr
#' 
#' @examples
#' thresh<-GOthresh(x_m,x_fdr$sampleFDR,enrichThresh=2,countThresh=2,pvalThresh=0.1,fdrThresh=0.100)
#'
#' @return returns a subset of matrix (m joined with fdr$sampleFDR) with entries meeting all thresholds
#' 
#' @export
GOthresh<-
  function(m,sampleFDR,enrichThresh,countThresh,pvalThresh,fdrThresh) {
    mm<-merge(m,sampleFDR,by=0)
    
    w1<-which(mm[,"ENRICHMENT"]>=enrichThresh)
    w2<-which(mm[,"SAMPLE"]>=countThresh)
    w3<-which(mm[,"OBSERVED_LOG10(P)"]<=log10(pvalThresh))
    w4<-which(mm[,"FDR"]<=fdrThresh)
    w12<-intersect(w1,w2)
    w123<-intersect(w12,w3)
    w1234<-intersect(w123,w4)
    
    return(mm[w1234,])
  }

#' GOheatmap
#' 
#' @description generate a matrix to be used as input to a heat map
#'
#' @param sampleList character list of gene names
#' @param x DB component of return value of GOtable3()
#' @param thresh output of GOthresh()
#' @param fdrThresh numeric value of FDR acceptance threshold
#' @param verbose integer vector representing classes
#'
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
#' heatmap<-GOheatmap(cluster52,GOGOA3$ontologies[["biological_process"]],x_thresh,verbose=1)
#' }
#'
#' @return returns a matrix to be used as input to a heat map
#' 
#' @export
GOheatmap<-
  function(sampleList,x,thresh,fdrThresh=.105, verbose) {
    genes<-vector("character")
    
    for(cat in thresh[,"Row.names"]) {
      w<-which(x[,"GO_NAME"]==cat)
      genes<-unique(c(genes,x[w,"HGNC"])) # genes that map to a 'good' category
    }
    
    vprint(-1,verbose,c("TOTAL NUMBER OF ONTOLOGY GENES MAPPING TO A SIGNIFICANT ONTOLOGY CATEGORY: ",length(genes)))
    vprint(-1,verbose,c("NUMBER OF GENES IN INPUT SAMPLE GENELIST: ",length(sampleList)))
    
    genes<-intersect(genes,sampleList)
    
    vprint(-1,verbose,c("NUMBER OF GENES IN INPUT SAMPLE GENELIST MAPPING TO A SIGNIFICANT ONTOLOGY CATEGORY: ",length(genes)))
    
    if(length(genes)==0)
      return(NULL)
    
    # columns of m are genes
    # rows of m are categories
    m<-matrix(fdrThresh,nrow=nrow(thresh),ncol=length(genes))
    rownames(m)<-thresh[,"Row.names"]
    colnames(m)<-genes
    # value m[cat,gene] is 0 if gene does not map to cat
    # otherwise m[cat,gene] is FDR[cat]
    for(r in 1:nrow(thresh)) {
      cat<-thresh[r,"Row.names"]
      w<-which(x[,"GO_NAME"]==cat)
      genes<-unique(x[w,"HGNC"])
      genes<-intersect(genes,sampleList)
      
      m[cat,genes]<-thresh[r,"FDR"]
    }
    
    return(m)
  }

# section B contains randomization functions

#' RCPD
#' 
#' @description prepare a cpd of p values from randomized gene sets
#' 
#' @param GOGOA3 return value of subsetGOGOA()
#' @param tablePop return value of GOtable3()
#' @param geneList character vector lisgt of genes to randomize
#' @param nrand integer number of randomizations
#' @param ontology c("molecular_function","cellular_component","biological_process")
#' @param hyper return value of GOhypergeometric3() from real (nonrandom) data
#' @param subd character string pathname for directory containing sink.txt
#' @param opt integer 0:1 parameter used to select randomization method
#' 
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
#' rcpd<-RCPD(GOGOA3,x_tablePop31,10,3,"biological_process",x_hyper1,tempdir(),0)
#' }
#'
#' @details the cpd of the randomizations is to be used for estimating
#' the false discovery rate (FDR) of the real sampled genes
#'
#' @return returns a histogram of log10(p)
#' 
#' @export
RCPD<-
  function(GOGOA3,tablePop,geneList,nrand,ontology,hyper,subd,opt) {
    
    if(opt!=0 & opt !=1)
      stop(sprintf("RCPD opt must be 0 or 1, not %d",opt))
    
    p<-vector("numeric",0)
    lh<-length(hyper)  
    for(i in 1:nrand) {
      if(opt==0) { # use original database with randomized geneLists
        sampleList<-randSubsetGeneList(GOGOA3$genes[[ontology]],length(geneList))
        tableSample<-GOtable3(sampleList,GOGOA3$ontologies[[ontology]])
      }
      if(opt==1) { # use original geneList with internally scrambled genes database
        GOGOA3R<-randomGODB(GOGOA3,FALSE)
        tableSample<-GOtable3(sampleList,GOGOA3R$ontologies[[ontology]])
      }
      
      p0<-GOhypergeometric3(tableSample,tablePop)
      sink(sprintf("%s/%s",subd,"sink.txt"),append=TRUE)
      print(c("RCPD",i),quote=FALSE)
      print(p0[1:10],quote=FALSE)
      print(c(log10(p0[1]),log10(p0[10])),quote=FALSE)
      sink()
      
      l<-min(lh,length(p0))
      if(l>2)
        plot(hyper[1:l],p0[1:l])
      p<-c(p,p0)
      
    }
    
    #plot.ecdf(ecdf(log10(p)))
    if(length(p)>2)
      return(ecdf(log10(p)))
    return(NULL)
  }

#' randSubsetGeneList
#' 
#' @description retrieve n unique random genes
#' 
#' @param geneList character vector geneList
#' @param ngenes integer desired number of random genes
#' 
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
#' genes<-randSubsetGeneList(GOGOA3$genes[["biological_process"]],20)
#' }
#'  
#' @return returns a character vector of genes
#' 
#' @export
randSubsetGeneList<-
  function(geneList,ngenes) {
    return(geneList[sample(length(geneList),ngenes)])
  }

# section C contains preprocessing functions

#' preprocessDB
#'
#' @description driver to perform several preprocessing steps:
#'  quick peek
#'  trim small and large categories
#'  is the database for human species
#'  validate validated HGNC symbols in sampleList
#'  determine up to date (ie, contains GOGOA3$species) or legacy version of human database
#'  
#' @param sampleList character list of gene names
#' @param GOGOA3 return value of subsetGOGOA()
#' @param ontology character string c("molecular_function", "cellular_component", "biological_process")
#' @param mn integer param passed to trimGOGOA3, min size threshold for a category
#' @param mx integer param passed to trimGOGOA3, max size threshold for a category
#' @param thresh numerical paramter passed to checkGeneListVsDB()
#' @param verbose integer vector representing classes
#'
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
#' pp<-preprocessDB(cluster52,GOGOA3,"biological_process",20,200,0.5,3)
#' }
#' 
#' @return returns a list whose components are a trimmed version of GOGOA3 and (for human)
#' a sampleList with validated HGNC symbols
#' 
#' @export
preprocessDB<-
  function(sampleList,GOGOA3,ontology,mn,mx,thresh,verbose) {
    
    l<-list()
    
    vprint(-1,verbose,"GOMINER DATABASE QUICK PEEK:")
    vprint(-1,verbose,GOGOA3$ontologies[[ontology]][1:5,])
    
    if(!is.null(mn) & !is.null(mx))
      GOGOA3<-hitterBeforeAfterDriver(GOGOA3,mn,mx,verbose)
    
    #####sampleList<-intersect(unique(sampleList),GOGOA3$genes[[ontology]])
    
    if(human(GOGOA3))
      sampleList<-validHGNCSymbols(sampleList)$geneList
    
    sampleList<-intersect(unique(sampleList),GOGOA3$genes[[ontology]])
    
    checkGeneListVsDB(sampleList,ontology,GOGOA3,thresh=0.5,verbose)
    
    l$GOGOA3<-GOGOA3
    l$sampleList<-sampleList
    
    return(l)
  }

#' hitterBeforeAfterDriver
#' 
#' @description driver to invoke hitters2() and trimGOGOA3()
#'  
#' @param GOGOA3 return value of minimalistGODB::buildGODatabase()
#' @param mn integer minimum category size
#' @param mx integer maximum category size
#' @param verbose integer vector representing classes
#' 
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # This example is given in full detail in the package vignette.
#' # You can generate GOGOA3.RData using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO
#' dir<-"/Users/barryzeeberg/personal/GODB_RDATA/goa_human/"
#' load(sprintf("%s/%s",dir,"GOGOA3_goa_human.RData"))
#' geneList<-GOGOA3$ontologies[["biological_process"]][1:10,"HGNC"]
#' GOGOA3tr<-hitterBeforeAfterDriver(GOGOA3,mn=20,mx=200,1)
#' }
#' 
#' @return returns the return value of trimGOGOA3()
#' 
#' @export
hitterBeforeAfterDriver<-
  function(GOGOA3,mn=20,mx=200,verbose) {
    vprint(1,verbose,"before trimming")
    hitters2(GOGOA3,verbose)
    GOGOA3tr<-trimGOGOA3(GOGOA3,mn,mx,verbose)
    vprint(1,verbose,"after trimming")
    hitters2(GOGOA3tr,verbose)
    
    return(GOGOA3tr)
  }

#' hitters2
#' 
#' @description determine the number of mappings for the top several genes
#'  
#' @param GOGOA3 return value of minimalistGODB::buildGODatabase()
#' @param verbose integer vector representing classes
#' 
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # This example is given in full detail in the package vignette.
#' # You can generate GOGOA3.RData using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO
#' dir<-"/Users/barryzeeberg/personal/GODB_RDATA/goa_human/"
#' load(sprintf("%s/%s",dir,"GOGOA3_goa_human.RData"))
#' geneList<-GOGOA3$ontologies[["biological_process"]][1:10,"HGNC"]
#' hitters2(GOGOA3,1)
#' }
#' 
#' @return returns no value, but has side effect of printing information
#' 
#' @export
hitters2<-
  function(GOGOA3,verbose=1) {
    ontologies<-names(GOGOA3$ontologies)
    #for(ontology in ontologies) {
    for(ontology in "biological_process") {
      t<-sort(table(GOGOA3$ontologies[[ontology]][,"HGNC"]),decreasing=TRUE)
      vprint(1,verbose,ontology)
      vprint(1,verbose,dim(GOGOA3$ontologies[[ontology]]))
      vprint(1,verbose,t[1:3])
      
      t<-sort(table(GOGOA3$ontologies[[ontology]][,"GO"]),decreasing=FALSE)
      vprint(1,verbose,t[1:3])
      t<-sort(table(GOGOA3$ontologies[[ontology]][,"GO"]),decreasing=TRUE)
      vprint(1,verbose,t[1:3])
    }
  }

#' trimGOGOA3
#'
#' @description remove categories from GOGOA3 that are too small or too large
#' 
#' @param GOGOA3 return value of subsetGOGOA()
#' @param mn integer min size threshold for a category
#' @param mx integer max size threshold for a category
#' @param verbose integer vector representing classes
#' 
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # This example is given in full detail in the package vignette.
#' # You can generate GOGOA3.RData using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' 
#' GOGO3tr<-trimGOGOA3(GOGOA3,mn=2,mx=200,1)
#' }
#' 
#' @details
#' If a category is too small, it is unreliable for statistical evaluation
#' Also, in the extreme case of size = 1, then that category is essentially
#' equivalent to a gene rather than a category. Same is partially true for size = 2.
#' If a category is too large, it is too generic to be useful for categorization.
#' Finally, by trimming the database, analyses will run faster.
#' 
#' @return returns trimmed version of GOGOA3
#' 
#' @export
trimGOGOA3<-
  function(GOGOA3,mn,mx,verbose) {
    if(!("tcats" %in% names(GOGOA3$stats)))
      stop("trimGOGOA3(): GOGOA3$stats$tcats not present. please use updated version of GOGOA3")
    onts<-names(GOGOA3$stats$tcats)
    for(ont in onts) {
      w<-which((GOGOA3$stats$tcats[[ont]][,1] >= mn) & (GOGOA3$stats$tcats[[ont]][,1] <= mx))
      good<-rownames(GOGOA3$stats$tcats[[ont]])[w]
      vprint(-1,verbose,c("GOOD",good[1:3]))
      vprint(-1,verbose,c(ont,nrow(GOGOA3$ontologies[[ont]]),length(w)))
      w1<-which(GOGOA3$ontologies[[ont]][,"GO_NAME"] %in% good)
      vprint(-1,verbose,length(w1))
      GOGOA3$ontologies[[ont]]<-GOGOA3$ontologies[[ont]][w1,]
      vprint(-1,verbose,c(ont,nrow(GOGOA3$ontologies[[ont]]),length(w)))
    }
    return(GOGOA3)
  }

#' human
#' 
#' @description determine if database represents human species
#' 
#' @param GOGOA3 return value of subsetGOGOA()
#' @param verbose integer vector representing classes
#' 
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
#' hum<-human(GOGOA3)
#' 
#' load("/Users/barryzeeberg/personal/GODB_RDATA/sgd/GOGOA3_sgd.RData")
#' hum<-human(XENOPUS,1)
#' }
#' 
#' @return returns Boolean TRUE if species is human
#'  
#' @export
human<-
  function(GOGOA3,verbose=TRUE) {
    w<-which(names(GOGOA3)=="species")
    if(length(w)==0) { # original version of GOGOA3 has only human genes
      message("You are using the older version of the human database")
      message("Next time you might want to replace it with the updated version")
      message("Good news, you can now use alternative databases like 'XENOPUS' to study different species!")
      message("You can generate it using the updated version of the package 'minimalistGODB'")
      message(" or retrieve it from https://github.com/barryzee/GO/tree/main/databases")
      
      return(TRUE)
    }
    
    if(substr(GOGOA3$species,1,9)=="goa_human") {
      if(verbose)
        message("Congratulations you are using the most up to date version of the human database")
      return(TRUE)
    }
    
    if(verbose)
      message(sprintf("Congratulations you are using the most up to date version of the %s database",GOGOA3$species))
    
    return(FALSE)
  }

#' validHGNCSymbols
#' 
#' @description convert outdated HGNC symbols to current HGNC symbols
#' 
#' @param geneList character vector of HGNC symbols
#' 
#' @details removes NA and /// from output of checkGeneSymbols()
#' 
#' @examples
#' geneList<-c("FN1", "tp53", "UNKNOWNGENE","7-Sep",
#'  "9/7", "1-Mar", "Oct4", "4-Oct","OCT4-PG4", "C19ORF71",
#'   "C19orf71")
#' l<-validHGNCSymbols(geneList)
#'
#' @return returns list of mapping table and vector of current HGNC symbols
#' 
#' @export
validHGNCSymbols<-
  function(geneList) {
    l<-list()
    
    x<-suppressMessages(suppressWarnings(checkGeneSymbols(geneList)))
    l$map<-x
    l$geneList<-setdiff(unique(x[,"Suggested.Symbol"]),NA) # exclude NA
    l$geneList<-unlist(strsplit(l$geneList," /// ")) # split 'gene1 /// gene2'
    
    return(l)
  }

#' checkGeneListVsDB
#' 
#' @description determine if gene list and database contain compatible identifiers 
#' 
#' @param geneList character list of gene names
#' @param GOGOA3 return value of subsetGOGOA()
#' @param ontology character string c("molecular_function", "cellular_component", "biological_process")
#' @param thresh numeric acceptance threshold for fraction of gene list matching database identifiers
#' @param verbose integer vector representing classes
#' 
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
#' checkGeneListVsDB(geneList=cluster52,ontology="biological_process",
#'  GOGOA3,thresh=0.5,verbose=TRUE)
#' 
#' # supposed to generate error message
#' load("/Users/barryzeeberg/personal/GODB_RDATA/sgd/GOGOA3_sgd.RData")
#' checkGeneListVsDB(geneList=xenopusGenes,ontology="biological_process",
#'  GOGOA3,thresh=0.5,verbose=TRUE)
#' }
#' 
#' @return returns no value, but may have side effect of aborting the computation
#' 
#' @export
checkGeneListVsDB<-
  function(geneList,ontology,GOGOA3,thresh=0.5,verbose=FALSE) {
    u<-unique(geneList)
    lu<-length(u)
    x<-intersect(GOGOA3$genes[[ontology]],u)
    lx<-length(x)
    diff1<-setdiff(GOGOA3$genes[[ontology]],u)
    diff2<-setdiff(u,GOGOA3$genes[[ontology]])
    
    match<-length(x)/length(u)
    
    if(match<thresh)
      warning("SIGNIFICANT MISMATCH BETWEEN IDENTIFIERS IN GENELIST AND DATABASE")
    if((match<thresh)) {
      print(sprintf("MATCH IS %f",match))
      print(sprintf("GENE LIST = %d INTERSECT = %d DIFF1 = %d DIFF2 = %d",lu,lx,length(diff1),length(diff2)))
      print("EXAMPLE IDENTIFIERS IN DATABASE:")
      print(diff1[1:5])
      print("EXAMPLE IDENTIFIERS IN GENELIST:")
      print(diff2[1:5])
    }
    if(match<thresh) 
      stop("ABORTING BECAUSE OF IDENTIFIER MISMATCH!!")
  }
