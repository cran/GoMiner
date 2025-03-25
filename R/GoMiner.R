#' validHGNCSymbols
#'
#' @import minimalistGODB
#' @import HGNChelper
#' @import stats
#' @importFrom gplots heatmap.2
#' @import grDevices
#' @import utils
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
    x<-suppressMessages(suppressWarnings(checkGeneSymbols(geneList,expand.ambiguous = FALSE)))
    l$map<-x
    l$geneList<-setdiff(unique(x[,"Suggested.Symbol"]),NA) # exclude NA
    l$geneList<-unlist(strsplit(l$geneList," /// ")) # split 'gene1 /// gene2'
    
    return(l)
  }

#' randSubsetGeneList
#' 
#' @description retrieve n unique random genes
#' 
#' @param geneList character vector geneList
#' @param ngenes integer desired number of random genes
#' 
#' @examples
#' #load("data/GOGOA3small.RData")
#' genes<-randSubsetGeneList(GOGOA3small$genes[["biological_process"]],20)
#'  
#' @return returns a character vector of genes
#' 
#' @export
randSubsetGeneList<-
  function(geneList,ngenes) {
    return(geneList[sample(length(geneList),ngenes)])
  }

#' GOtable3
#' 
#' @description tabulate number of geneList mappings to GO categories
#' 
#' @param hgncList character list of gene names
#' @param DB selected ontology branch of return value of subsetGOGOA
#' 
#' @examples
#' #load("data/GOGOA3small.RData")
#' DB<-GOGOA3small$ontologies[["biological_process"]]
#' 
#' # housekeeping genes downloaded from https://housekeeping.unicamp.br/?download
#' #load("data/Housekeeping_Genes.RData")
#' hgncList<-Housekeeping_Genes[,"Gene.name"]
#' x<-GOtable3(hgncList,DB)
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
    hgncList<-validHGNCSymbols(hgncList)$geneList
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
#' #load("data/x_tableSample3.RData")
#' #load("data/x_tablePop3.RData")
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

#' GOhypergeometric
#' 
#' @description compute the hypergeometric p value for gene enrichment in a GO category
#' 
#' @param tableSample3 sample return value of GOtable3()
#' @param tablePop3 population return value of GOtable3()
#' 
#' @examples
#' #load("data/x_tableSample3.RData")
#' #load("data/x_tablePop3.RData")
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
    return(m[order(m[,"p"]),"p"])
  }

#' FDR
#' 
#' @description compute the false discovery rate (FDR) of the hypergeometric
#' p values of genes mapping to gene ontology (GO) categories
#'
#' @param sampleList character vector of user-supplied genes of interest
#' @param GOGOA3 return value of subsetGOGOA()
#' @param nrand integer number of randomizations
#' @param ONT c("molecular_function","cellular_component","biological_process") 
#' 
#' @examples
#' #load("data/GOGOA3small.RData")
#' sampleList<-randSubsetGeneList(GOGOA3small$genes[["biological_process"]],10)
#' fdr<-FDR(sampleList,GOGOA3small,nrand=100,"biological_process")
#' 
#' @return returns a list with FDR information
#' 
#' @export
FDR<-
  function(sampleList,GOGOA3,nrand,ONT) {
    l<-list()
    DB<-GOGOA3$ontologies[[ONT]]
    sampleList<-validHGNCSymbols(sampleList)$geneList
    sampleList<-intersect(unique(sampleList),GOGOA3$genes[[ONT]])
    ngenes<-length(sampleList)
    tableSample3<-GOtable3(sampleList,DB)
    
    popList<-GOGOA3$genes[[ONT]]
    tablePop3<-GOtable3(popList,DB)
    
    hyper<-GOhypergeometric3(tableSample3,tablePop3)
    
    # rcpd() is a function that interpolates observed p values
    rcpd<-RCPD(GOGOA3,ngenes,nrand,ONT) 
    
    sampleFDR<-matrix(nrow=length(names(hyper)),ncol=2)
    rownames(sampleFDR)<-names(hyper)
    colnames(sampleFDR)<-c("OBSERVED_LOG10(P)","FDR")
    for(cat in names(hyper)) {
      x<-log10(hyper[cat])
      sampleFDR[cat,"OBSERVED_LOG10(P)"]<-x
      sampleFDR[cat,"FDR"]<-rcpd(x)
    }

    l$rcpd<-rcpd
    l$sampleFDR<-sampleFDR
    
    return(l)
  }

#' GOthresh
#' 
#' @description retrieve lines of m that meet both enrichThresh and countThresh
#' 
#' @param m return value of GOenrich3()
#' @param sampleFDR component of return value of RCPD()
#' @param enrichThresh  numerical acceptance threshold for enrichment
#' @param countThresh numerical acceptance threshold for gene count
#' @param fdrThresh numerical acceptance threshold for fdr
#' @examples
#' #load("data/x_m.RData")
#' #load("data/x_fdr.RData")
#' thresh<-GOthresh(x_m,x_fdr$sampleFDR,enrichThresh=2,countThresh=2,fdrThresh=0.100)
#'
#' @return returns a subset of matrix (m joined with fdr$sampleFDR) with entries meeting all thresholds
#' 
#' @export
GOthresh<-
  function(m,sampleFDR,enrichThresh,countThresh,fdrThresh) {
    mm<-merge(m,sampleFDR,by=0)
    w1<-which(mm[,"ENRICHMENT"]>=enrichThresh)
    w2<-which(mm[,"SAMPLE"]>=countThresh)
    w3<-which(mm[,"FDR"]<=fdrThresh)
    w12<-intersect(w1,w2)
    w123<-intersect(w12,w3)
    
    return(mm[w123,])
  }

#' RCPD
#' 
#' @description prepare a cpd of p values from randomized gene sets
#' 
#' @param GOGOA3 return value of subsetGOGOA()
#' @param ngenes integer number of genes to randomize
#' @param nrand integer number of randomizations
#' @param ONT c("molecular_function","cellular_component","biological_process")
#' 
#' @examples
#' #load("data/GOGOA3small.RData")
#' rcpd<-RCPD(GOGOA3small,ngenes=100,nrand=10,ONT="biological_process")
#'
#' @details the cpd of the randomizations is to be used for estimating
#' the false discovery rate (FDR) of the real sampled genes
#'
#' @return returns a histogram of log10(p)
#' 
#' @export
RCPD<-
  function(GOGOA3,ngenes,nrand,ONT) {
    p<-vector("numeric",0)
   
    DB<-GOGOA3$ontologies[[ONT]]
    
    popList<-GOGOA3$genes[[ONT]]
    tablePop<-GOtable3(popList,DB)
    
    for(i in 1:nrand) {
      sampleList<-randSubsetGeneList(GOGOA3$genes[[ONT]],ngenes)
      tableSample<-GOtable3(sampleList,DB)
      p<-c(p,GOhypergeometric3(tableSample,tablePop))
    }
    
    return(ecdf(log10(p)))
  }

#' GOheatmap
#' 
#' @description generate a matrix to be used as input to a heat map
#'
#' @param sampleList character list of gene names
#' @param x DB component of return value of GOtable3()
#' @param thresh output of GOthresh()
#' @param fdrThresh numeric value of FDR acceptance threshold
#'
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO
#' #load("~/GODB_RDATA/GOGOA3.RData")
#' ONT<-"biological_process"
#' DB<-GOGOA3$ontologies[[ONT]]
#'
#' #load("data/cluster52.RData")
#' sampleList<-cluster52
#'
#' #load("data/x_thresh.RData")
#' heatmap<-GOheatmap(sampleList,DB,x_thresh)
#' }
#'
#' @return returns a matrix to be used as input to a heat map
#' 
#' @export
GOheatmap<-
  function(sampleList,x,thresh,fdrThresh=.105) {
    genes<-vector("character")
    for(cat in thresh[,"Row.names"]) {
      w<-which(x[,"GO_NAME"]==cat)
      genes<-unique(c(genes,x[w,"HGNC"])) # genes that map to a 'good' category
    }
    
    message(c("TOTAL NUMBER OF ONTOLOGY GENES MAPPING TO A SIGNIFICANT ONTOLOGY CATEGORY: ",length(genes)))
    
    message(c("NUMBER OF GENES IN INPUT SAMPLE GENELIST: ",length(sampleList)))
    
    genes<-intersect(genes,sampleList)
    
    message(c("NUMBER OF GENES IN INPUT SAMPLE GENELIST MAPPING TO A SIGNIFICANT ONTOLOGY CATEGORY: ",length(genes)))
    
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

#' GoMiner
#' 
#' @description driver to generate heatmap
#'
#' @param title character string descriptive title
#' @param dir character string full pathname to the directory acting result repository
#' @param sampleList character list of gene names
#' @param GOGOA3 return value of subsetGOGOA()
#' @param ONT character string c("molecular_function", "cellular_component", "biological_process")
#' @param enrichThresh  numerical acceptance threshold for enrichment
#' @param countThresh numerical acceptance threshold for gene count
#' @param fdrThresh numerical acceptance threshold for fdr
#' @param nrand numeric number of randomizations to compute FDR
#'
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO
#' load("~/GODB_RDATA/GOGOA3.RData")
#' load("data/cluster52.RData")
#' l<-GoMiner("Cluster52",tempdir(),cluster52,
#'  GOGOA3,ONT="biological_process",enrichThresh=2,
#'  countThresh=5,fdrThresh=0.10,nrand=10)
#' }
#'
#' @return returns a matrix suitable to generate a heatmap
#' 
#' @export
GoMiner<-
  function(title=NULL,dir,sampleList,GOGOA3,ONT,enrichThresh=2,
            countThresh=5,fdrThresh=0.10,nrand=100) {
    l<-list()
    
    DB<-GOGOA3$ontologies[[ONT]]
  
    stamp<-gsub(":","_",format(Sys.time(), "%a_%b_%d_%Y_%X"))
    if(is.null(title))
      title<-"GMresults"
    subd<-sprintf("%s/%s_%s",dir,title,stamp)
    dir.create(subd)
    
    args<-sprintf("%s/args.txt",subd)
    write(sprintf("Output Results Directory: %s",dir),file=args,append=FALSE)
    write(sprintf("Ontology: %s",ONT),file=args,append=TRUE)
    write(sprintf("enrichThresh: %f",enrichThresh),file=args,append=TRUE)
    write(sprintf("countThresh: %f",countThresh),file=args,append=TRUE)
    write(sprintf("fdrThresh: %f",fdrThresh),file=args,append=TRUE)
    write(sprintf("nrand: %d",nrand),file=args,append=TRUE)
    
    write(sprintf("sampleList:"),file=args,append=TRUE)
    write(sprintf("%s",sampleList),file=args,append=TRUE)
    
    metadata<-sprintf("%s/metadata.txt",subd)
    
    
    write(sprintf("Descriptive Title is %s",title),file=metadata,append=FALSE)
    message(sprintf("Gominer Results Directory is %s",subd))
    write(sprintf("Gominer Results Directory is %s",subd),file=metadata,append=TRUE)
    
    tableSample3<-GOtable3(sampleList,DB)
    popList<-GOGOA3$genes[[ONT]]
    tablePop3<-GOtable3(popList,DB)
 
    #x_tableSample3<-tableSample3
    #save(x_tableSample3,file="data/x_tableSample3.RData")
    #x_tablePop3<-tablePop3
    #save(x_tablePop3,file="data/x_tablePop3.RData")
    m<-GOenrich3(tableSample3,tablePop3)
     
    fdr<-FDR(sampleList,GOGOA3,nrand=100,ONT)
 
    #x_m<-m
    #save(x_m,file="data/x_m.RData")
    #x_fdr<-fdr
    #save(x_fdr,file="data/x_fdr.RData")
    thresh<-GOthresh(m,fdr$sampleFDR,
      enrichThresh=enrichThresh,countThresh=countThresh,fdrThresh=fdrThresh)
    
    #x_thresh<-thresh
    #save(x_thresh,file="data/x_thresh.RData")
    
    gce<-merge(tableSample3$gce,thresh,by.x="GO_NAME",by.y="Row.names")
    gce<-gce[order(gce[,"FDR"]),]
    
    heatmap<-GOheatmap(sampleList,DB,thresh,fdrThresh)

    message(c("NUMBER OF SIGNIFICANT ONTOLOGY CATEGORIES: ",nrow(heatmap)))
    write(sprintf("NUMBER OF SIGNIFICANT ONTOLOGY CATEGORIES: %d",nrow(heatmap)),file=metadata,append=TRUE)
    
    svgWidth<-2*(4.75 + ncol(heatmap)*.0792)*30/20.75
    message(c("SVG WIDTH: ",svgWidth))
    write(sprintf("SVG WIDTH: %f",svgWidth),file=metadata,append=TRUE)
    svgHeight<-2*(1.25 + nrow(heatmap)*.0403)*20/12.00
    message(c("SVG HEIGHT: ",svgHeight))
    write(sprintf("SVG HEIGHT: %f",svgHeight),file=metadata,append=TRUE)
    
    file<-sprintf("%s/GoMiner_%d_%d.svg",subd,nrow(heatmap),ncol(heatmap))
    message(sprintf("SVG Heatmap is stored as %s",file))
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
    
    l$thresh<-thresh
    l$gce<-gce
    
    return(l)
  }