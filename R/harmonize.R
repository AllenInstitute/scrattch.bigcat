library(Matrix)
#library(ggplot2)
library(matrixStats)

jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
blue.red <-colorRampPalette(c("blue", "white", "red"))


#' Sample sets lists
#'
#' @param cells.list 
#' @param cl.list 
#' @param cl.sample.size 
#' @param sample.size 
#'
#' @return
#' @export
#'
#' @examples
sample_sets_list <- function(cells.list, cl.list, cl.sample.size=100, sample.size=5000)
  {
    for(x in names(cells.list)){
      if(length(cells.list[[x]]) > sample.size){
        if(is.null(cl.list[[x]])){
          cells.list[[x]] = sample(cells.list[[x]], sample.size)
        }
        else{
          tmp.cl = cl.list[[x]][cells.list[[x]]]
          if(is.factor(tmp.cl)){
            tmp.cl = droplevels(tmp.cl)
          }
          cl.size = table(tmp.cl)
          cells.list[[x]] = sample_cells(tmp.cl, max(cl.sample.size,round(sample.size/sum(cl.size >= 4))))
        }
      }
    }
    return(cells.list)
  }




sample_cl_dat <- function(comb.dat, sets, cl, cl.sample.size=200)
  {
    dat.list = sapply(sets, function(set){
      print(set)
      select.cells = intersect(row.names(comb.dat$meta.df)[comb.dat$meta.df$platform==set], names(cl))
      if(length(select.cells)==0){
        return(NULL)
      }
      tmp.cl = cl[select.cells]
      if(is.factor(tmp.cl)){
        tmp.cl = droplevels(tmp.cl)
      }
      select.cells = sample_cells(tmp.cl,cl.sample.size)
      print(length(select.cells))
      get_logNormal(comb.dat$dat.list[[set]], select.cells)
    },simplify=F)
    dat.list = dat.list[!sapply(dat.list,is.null)]
    return(dat.list)
  }

knn_combine <- function(result.1, result.2)
{
  knn.index = rbind(result.1[[1]], result.2[[1]])
  knn.distance = rbind(result.1[[2]], result.2[[2]])
  return(list(knn.index, knn.distance))
}

#' get knn batch
#'
#' @param dat 
#' @param ref.dat 
#' @param k 
#' @param method 
#' @param dim 
#' @param batch.size 
#' @param mc.cores 
#'
#' @return
#' @export
#'
#' @examples
get_knn_batch <- function(dat, ref.dat, k, method="cor", dim=NULL, batch.size, mc.cores=1,return.distance=FALSE,transposed=TRUE, index=NULL,clear.index=FALSE, ntrees=50) 
  {
    library(BiocNeighbors)
    if(return.distance){
      fun = "knn_combine"
    }
    else{
      fun = "rbind"
    }
    if(is.null(index) & method %in% c("Annoy.Euclidean", "Annoy.Cosine", "cor")) {
      if(transposed){
        map.ref.dat = Matrix::t(ref.dat)
      }
      else{
        map.ref.dat = ref.dat
      }
      if (method == "cor") {
        map.ref.dat = map.ref.dat - rowMeans(map.ref.dat)
        map.ref.dat = l2norm(map.ref.dat, by = "row")
      }
      if (method == "Annoy.Cosine") {
        map.ref.dat = l2norm(map.ref.dat, by = "row")
      }
      index = buildAnnoy(map.ref.dat, ntrees = ntrees)
      rm(map.ref.dat)
      #gc()
    }    
    if(transposed){      
      results <- batch_process(x=1:ncol(dat), batch.size=batch.size, mc.cores=mc.cores, .combine=fun, FUN=function(bin){
        get_knn(dat=dat[row.names(ref.dat),bin,drop=F], ref.dat=ref.dat, k=k, method=method, dim=dim,return.distance=return.distance, transposed=transposed,index=index,ntrees=ntrees)
      })
    }
    else{
      results <- batch_process(x=1:nrow(dat), batch.size=batch.size, mc.cores=mc.cores, .combine=fun, FUN=function(bin){
        get_knn(dat=dat[bin,colnames(ref.dat),drop=F], ref.dat=ref.dat, k=k, method=method, dim=dim,return.distance=return.distance, transposed=transposed,index,ntrees=ntrees)
      })
    }
    if(clear.index){
      cleanAnnoyIndex(index)
    }
    else{
      if(!return.distance){
        results = list(knn.index=results)
      }
      results$index=index
    }
    return(results)
  }

knn_cor <- function(ref.dat, query.dat, k = 15)
  {
                                        
    sim = cor(as.matrix(query.dat), as.matrix(ref.dat))
    sim[is.na(sim)] = 0
    return(sim_knn(sim, k=k))
  }


#' Get KNN
#'
#' @param dat 
#' @param ref.dat 
#' @param k 
#' @param method 
#' @param dim 
#'
#' @return
#' @export
#'
                                        #' @examples

get_knn <- function(dat, ref.dat, k, method ="cor", dim=NULL,index=NULL, build.index=FALSE, transposed=TRUE, return.distance=FALSE, ntrees=100)
  {
    if(transposed){
      cell.id = colnames(dat)
    }
    else{
      cell.id= row.names(dat)
    }
    
    if(transposed){
      if(is.null(index)){
        ref.dat = Matrix::t(ref.dat)
      }
      dat = Matrix::t(dat)
    }
    if(method %in% c("Euclidean","Cosine")){
      build.index=FALSE
      index=NULL
    }
    if(method=="RANN"){
      knn.result = RANN::nn2(ref.dat, dat, k=k)
    }
    else if(method %in% c("Annoy.Euclidean", "Annoy.Cosine","cor","Cosine")){
      library(BiocNeighbors)
      if(is.null(index)){
        if(method=="cor"){
          ref.dat = ref.dat - Matrix::rowMeans(ref.dat)
          ref.dat = l2norm(ref.dat,by = "row")
        }
        if (method %in% c("Annoy.Cosine","Cosine")){
          ref.dat = l2norm(ref.dat,by = "row")
        }
        if(build.index){
          index= buildAnnoy(ref.dat, ntrees=ntrees)
        }
      }
      if (method %in% c("Annoy.Cosine","Cosine")){
        dat = l2norm(dat,by="row")
      } 
      if (method %in% c("Annoy.cor","cor")){
        dat = dat - Matrix::rowMeans(dat)
        dat = l2norm(dat,by = "row")
      }
      if (method %in% c("Annoy.Cosine","Annoy.Euclidean","Annoy.cor")){
        knn.result = queryAnnoy(X= ref.dat, query=dat, k=k, precomputed = index)
      }
      else if(method %in% c("Cosine", "Euclidean","cor")){
        knn.result = queryKNN(X= ref.dat, query=dat, k=k)
      }
    }
    else if(method == "CCA"){
      mat3 = crossprod(ref.dat, dat)
      cca.svd <- irlba(mat3, dim=dim)
      ref.dat = cca.svd$u
      dat = cca.svd$v      
      knn.result =get_knn(dat, ref.dat, method="cor")
    }
    else{
      stop(paste(method, "method unknown"))
    }    
    knn.index= knn.result[[1]]
    knn.distance = knn.result[[2]]
    row.names(knn.index) = row.names(knn.distance)=cell.id
    rm(dat)
    #gc()
    if(!return.distance){
      return(knn.index)
    }
    else{
      list(knn.index=knn.index, knn.distance=knn.distance)
    }
  }




cleanAnnoyIndex <- function(index)
  {
    unlink(index@path)
    rm(index)
  }

knn_combine <- function(result.1, result.2)
{
  knn.index = rbind(result.1[[1]], result.2[[1]])
  knn.distance = rbind(result.1[[2]], result.2[[2]])
  return(list(knn.index=knn.index, knn.distance=knn.distance))
}



###comb.dat include the following elements
###dat.list a list of data matrix
###ref.de.param.list the DE gene criteria for each reference dataset (optional)
###meta.df merged meta data for all datasets. 
###cl.list clusters for each dataset (optional)
###cl.df.list cluster annotations for each dataset (optional) 


prepare_harmonize_big <- function(dat.list, meta.df=NULL, cl.list=NULL, cl.df.list = NULL, de.param.list=NULL, de.genes.list=NULL, rename=TRUE)
  {
    common.genes = dat.list[[1]]$row_id
    if(length(dat.list)>1){
      for(x in 2:length(dat.list)){
        common.genes= intersect(common.genes, dat.list[[x]]$row_id)
      }
    }
    if(rename){
      for(x in names(dat.list)){
        dat.list[[x]]$col_id = paste(x, dat.list[[x]]$col_id, sep=".")
      }
      if(!is.null(cl.list)){
        for(x in names(cl.list)){
          names(cl.list[[x]]) = paste(x, names(cl.list[[x]]), sep=".")
        }
      }
    }
    
    platform = do.call("c",lapply(names(dat.list), function(p){
      dat = dat.list[[p]]
      setNames(rep(p, length(dat$col_id)), dat$col_id)
    }))
    #gene.counts <- do.call("c",lapply(names(dat.list), function(p){
    #  dat = dat.list[[p]]
    #  setNames(bg_colSums(dat > 0), colnames(dat))
    #}))
    df = data.frame(platform)
    if(!is.null(meta.df)){
      common.cells = intersect(row.names(meta.df), row.names(df))
      meta.df = cbind(meta.df[common.cells,,drop=F], df[common.cells,,drop=F])
    }
    else{
      meta.df = df
    }
    meta.df$platform = factor(meta.df$platform)
    all.cells = unlist(lapply(dat.list, function(x)x$col_id))
    comb.dat = list(dat.list=dat.list, meta.df = meta.df, cl.list=cl.list, cl.df.list = cl.df.list, de.genes.list = de.genes.list, de.param.list= de.param.list, common.genes=common.genes, all.cells= all.cells, type="big")
  }




select_joint_genes <-  function(comb.dat, ref.dat.list, select.cells = comb.dat$all.cells, maxGenes=2000, vg.padj.th=0.5, max.dim=20, top.n=100,rm.eigen=NULL, rm.th=rep(0.7,ncol(rm.eigen)),use.common=TRUE)
  {
    require(matrixStats)
    require(data.table)
    require(dplyr)
    select.genes = lapply(names(ref.dat.list), function(ref.set){
      ref.dat = ref.dat.list[[ref.set]]
      ref.cells=colnames(ref.dat)
      cat(ref.set, length(ref.cells),"\n")
      tmp.dat = ref.dat
      tmp.dat@x = 2^tmp.dat@x - 1
      vg = find_vg(tmp.dat)      
      rm(tmp.dat)
      gc()
      vg = vg %>% arrange(loess.padj,-abs(loess.z))
      if(nrow(vg) > maxGenes){
        select.genes = with(vg, gene[which(loess.padj < vg.padj.th | dispersion >3)])
        if(length(select.genes) < 5){
          return(NULL)
        }
        select.genes = head(select.genes, maxGenes+1000)
        rd = rd_PCA(norm.dat=ref.dat,select.genes, ref.cells, max.pca = max.dim,method="elbow")
        if(is.null(rd)){
          return(NULL)
        }
        if(!is.null(rm.eigen)){
          rd.dat = filter_RD(rd$rd.dat, rm.eigen, rm.th)
        }
        else{
          rd.dat = rd$rd.dat
        }
        if(is.null(rd.dat)){
          return(NULL)
        }
        rot = t(rd$pca$rotation[,colnames(rd.dat)])
        if(is.null(rot)){
          return(NULL)
        }
        rot.scaled = (rot  - rowMeans(rot))/rowSds(rot)
        gene.rank = t(apply(-abs(rot.scaled), 1, rank))
        select = gene.rank <= top.n & abs(rot.scaled ) > 2
        select.genes = colnames(select)[colSums(select)>0]
        if(length(select.genes) < 3){
          return(NULL)
        }
        vg=vg %>% filter(gene %in% select.genes)
      }
      vg$rank = 1:nrow(vg)
      vg
    })
    library(data.table)
    gene.score = rbindlist(select.genes)
    if(nrow(gene.score)==0){
      return(NULL)
    }
    gene.score = gene.score %>% group_by(gene) %>% summarize(gene.score = sum(maxGenes +1000 - rank)) %>% arrange(-gene.score)    
    select.genes= head(gene.score$gene, maxGenes)
    #gg.cons = gene_gene_cor_conservation(ref.dat.list, select.genes, select.cells)
    #select.genes = row.names(gg.cons)[gg.cons > conservation.th]
    if(use.common){
      common.genes=Reduce("intersect", lapply(ref.dat.list, row.names))
      select.genes= intersect(select.genes,common.genes)
    }
    return(select.genes)
  }

build_annoy_index <- function(ref.dat, knn.method, transposed=TRUE,ntrees=50)
  {
    library(BiocNeighbors)
    if(transposed){
      knn.ref.dat = t(ref.dat)
    }
    if(knn.method=="cor"){
      knn.ref.dat = knn.ref.dat - rowMeans(knn.ref.dat)
      knn.ref.dat = l2norm(knn.ref.dat, by = "row")
    }
    else if(knn.method=="Annoy.Cosine"){
      knn.ref.dat = l2norm(knn.ref.dat, by = "row")
    }          
    index = buildAnnoy(knn.ref.dat,ntrees=ntrees)
    rm(knn.ref.dat)
    gc()
    return(index)
  }

compute_knn <- function(comb.dat, select.genes, ref.list, ref.dat.list=NULL, select.sets=names(comb.dat$dat.list), select.cells=comb.dat$all.cells, k=15, cross.knn.method=c("Annoy.Cosine","cor"), self.knn.method=c("Annoy.Cosine"), block.size=10000, mc.cores=1)
  {    
    cat("Number of select genes", length(select.genes), "\n")
    cat("Get knn\n")
    dat.list = comb.dat$dat.list
    if(is.null(ref.dat.list)){
      for(ref.set in names(ref.list)){
        if(length(ref.list[[ref.set]]) <= k) {
          ##Not enough reference points to compute k
          next
        }
        dat = dat.list[[ref.set]]
        ref.cells = ref.list[[ref.set]]
        if(comb.dat$type =="mem"){
          ref.dat.list[[ref.set]] = dat[,ref.cells,drop=FALSE]
        }
        else{
          ref.dat.list[[ref.set]] = get_logNormal(dat, ref.cells, mc.cores=mc.cores)
        }
      }
    }
    knn.list = list()
    for(ref.set in names(ref.list)){
      cat("Ref set ", ref.set, "\n")
      ###index is the index of knn from all the cells              
      k.tmp = k
      if(length(ref.list[[ref.set]]) <= k*2) {
        k.tmp = round(k/2)
      }
      ref.cells = ref.list[[ref.set]]
      index=NULL
      ref.dat = ref.dat.list[[ref.set]]
      ref.dat = ref.dat[select.genes[select.genes %in% row.names(ref.dat)],ref.cells,,drop=FALSE]
      idx = match(ref.cells, comb.dat$all.cells)      
      tmp.knn.list =list()
      for(set in c(ref.set, setdiff(select.sets,ref.set))){
        dat = dat.list[[set]]
        cat("Set ", set, "\n")
        if(comb.dat$type=="mem"){
          map.cells=  intersect(colnames(dat), select.cells)          
          map.select.genes = select.genes[select.genes %in% intersect(row.names(ref.dat),row.names(dat))]
        }
        else{
          map.cells=  dat$col_id[dat$col_id %in% select.cells]
          map.select.genes = select.genes[select.genes %in% intersect(row.names(ref.dat),dat$row_id)]
        }        
        if(length(map.cells)==0){
          next
        }
        cat("map cells",length(map.cells),"\n")
        map.index=NULL
        if(!is.null(index)){
          if(length(map.select.genes)==ncol(index)){
            map.index=index
          }          
        }
        if(set==ref.set){
          method = self.knn.method
        }
        else{
          method = cross.knn.method
        }
        tmp.ref.dat = ref.dat
        if(length(map.select.genes)!=nrow(ref.dat)){
          tmp.ref.dat = tmp.ref.dat[map.select.genes,,drop=FALSE]
        }
        if(all(map.cells %in% colnames(ref.dat.list[[set]]))){
          dat = ref.dat.list[[set]][map.select.genes,map.cells,drop=FALSE]          
          knn.result=get_knn_batch(dat=dat, ref.dat = tmp.ref.dat, k=k.tmp, method = method, batch.size = block.size, mc.cores=mc.cores, index=map.index, transposed=TRUE)
        }
        else{
          if(comb.dat$type=="mem"){
            dat = comb.dat$dat.list[[set]][map.select.genes,map.cells,drop=FALSE]          
            knn.result=get_knn_batch(dat=dat, ref.dat = tmp.ref.dat, k=k.tmp, method = method, batch.size = block.size, mc.cores=mc.cores, index=map.index, transposed=TRUE)
          }
          else{
            knn.result=get_knn_batch_big(big.dat=dat.list[[set]], ref.dat = tmp.ref.dat,select.cells=map.cells, k=k.tmp, method = cross.knn.method, block.size = block.size, mc.cores=mc.cores, index=map.index)
          }
        }
        if(is.null(index)){
          index = knn.result$ref.index
        }
        else{
          cleanAnnoyIndex(knn.result$ref.index)
        }
        #if(!is.null(comb.dat$cl.list)){  
        #  test.knn = test_knn(knn, comb.dat$cl.list[[set]], colnames(ref.dat), comb.dat$cl.list[[ref.set]])          
        #  if(!is.null(test.knn)){
        #    cat("Knn", set, ref.set, method, "cl.score", test.knn$cl.score, "cell.score", test.knn$cell.score,"\n")
        #  }
        #}
        knn = knn.result[[1]]
        knn = matrix(idx[knn], nrow=nrow(knn), dimnames=list(row.names(knn), NULL))
        tmp.knn.list[[set]]=knn
      }
      knn = do.call("rbind",tmp.knn.list)
      if(!is.null(index)){
        cleanAnnoyIndex(index)
      }
      knn.list[[ref.set]] = knn[select.cells,]      
    }
    knn.comb = do.call("cbind",knn.list)
    return(knn.comb)
  }

  

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param comb.dat 
##' @param ref.sets 
##' @param select.sets 
##' @param merge.sets 
##' @param select.cells 
##' @param select.genes 
##' @param cross.knn.method 
##' @param self.knn.method 
##' @param k 
##' @param sample.size 
##' @param cl.sample.size 
##' @param block.size 
##' @param verbose 
##' @param mc.cores 
##' @param rm.eigen 
##' @param rm.th 
##' @param ... 
##' @return 
##' @author Zizhen Yao
knn_joint <- function(comb.dat, joint.rd.dat=NULL, ref.sets=names(comb.dat$dat.list), select.sets= names(comb.dat$dat.list), merge.sets=ref.sets, select.cells=comb.dat$all.cells, select.genes=NULL,cross.knn.method="Annoy.Cosine", self.knn.method = "Annoy.Euclidean", method="leiden", k=15,  sample.size = 50000, jaccard.sampleSize = 200000, cl.sample.size = 100, block.size = 10000, verbose=TRUE,mc.cores=1,rm.eigen=NULL, rm.th=0.7,max.dim=20,rm.genes=NULL,...)
{
  select.cells = comb.dat$all.cells[comb.dat$all.cells %in% select.cells]
  cat("Number of select cells", length(select.cells), "\n")
  cells.list = split(select.cells, comb.dat$meta.df[select.cells, "platform"])[select.sets]
  ref.list = cells.list[ref.sets]
  ref.sets = ref.sets[sapply(ref.list,length) >= sapply(comb.dat$de.param.list[ref.sets], function(x)x$min.cells)]
  if(length(ref.sets)==0){
    return(NULL)
  }
  ref.list = ref.list[ref.sets]
  if(!is.null(joint.rd.dat)){
      ref.cells=unlist(ref.list)
      knn.comb=get_knn(dat=joint.rd.dat[select.cells,], ref.dat = joint.rd.dat[ref.cells,], k=k, method = self.knn.method, transposed=FALSE)
      result = knn_jaccard_clust(knn.comb,prune=1/(ncol(knn.comb)-1), method=method)      
    }
  else{
    if(length(select.cells) < block.size){
      mc.cores=1
    }
    if(length(select.cells) >= sample.size){
      ref.list =  sample_sets_list(cells.list[ref.sets], comb.dat$cl.list[ref.sets], sample.size=sample.size, cl.sample.size = cl.sample.size)
      ref.list = sapply(ref.list, function(x){
        ref.cells=sample(x, min(sample.size, length(x)))
        ref.cells = comb.dat$all.cells[comb.dat$all.cells %in% ref.cells]
      },simplify=F)
    }
    ###Select genes for joint analysis
    cat("Get ref.dat.list\n")
    if(comb.dat$type=="mem"){
      ref.dat.list = sapply(ref.sets, function(ref.set){
        ref.dat= comb.dat$dat.list[[ref.set]]
        ref.dat= ref.dat[!row.names(ref.dat) %in% rm.genes, ref.list[[ref.set]],drop=FALSE]
      },simplify=F)
    }
    else{
      ref.dat.list = sapply(ref.sets, function(ref.set){
        ref.dat=get_logNormal(comb.dat$dat.list[[ref.set]], ref.list[[ref.set]])
        ref.dat[!row.names(ref.dat) %in% rm.genes,,drop=FALSE]
      },simplify=F)
    }
    if(is.null(select.genes)){
      select.genes = select_joint_genes(comb.dat, ref.dat.list = ref.dat.list,select.cells=select.cells, max.dim= max.dim,rm.eigen=rm.eigen, rm.th=rm.th)
    }
    if(length(select.genes) < 20){
      return(NULL)
    }
    ref.dat.list = sapply(ref.dat.list, function(ref.dat){
      tmp.genes = intersect(select.genes, row.names(ref.dat))
      tmp = Matrix::colSums(ref.dat[tmp.genes,]) == 0
      if(sum(tmp)>0){
        ref.dat[,!tmp,drop=F]
      }
      else{
        ref.dat
      }
    },simplify=F)
    ref.list = sapply(ref.dat.list, colnames, simplify=FALSE)
    cat("Get knn\n")
    knn.comb= compute_knn(comb.dat, select.genes=select.genes, ref.list=ref.list, ref.dat.list= ref.dat.list, select.sets=select.sets, select.cells=select.cells, k=k, cross.knn.method=cross.knn.method, self.knn.method=self.knn.method, block.size=block.size, mc.cores=mc.cores)
    if(is.null(knn.comb)){
      return(NULL)
    }
    #########
    sampled.cells = unlist(cells.list)
    if(length(sampled.cells)> sample.size){
      tmp.list =  sample_sets_list(cells.list, comb.dat$cl.list[select.sets], sample.size=sample.size, cl.sample.size = cl.sample.size)
      sampled.cells = unlist(tmp.list)
      if(length(sampled.cells)>jaccard.sampleSize){
        sampled.cells = sample(sampled.cells, jaccard.sampleSize)
      }
      sampled.cells=union(sampled.cells, unlist(ref.list))
    }     
    result = knn_jaccard_clust(knn.comb[sampled.cells,],prune=1/(ncol(knn.comb)-1), method=method)
  }
  result$knn = knn.comb
  cl = result$cl
  result$ref.list = ref.list
  if(length(cl) < nrow(result$knn)){
    diff.cells = setdiff(row.names(result$knn), names(cl))
    pred.df = predict_knn(result$knn[diff.cells,,drop=F], comb.dat$all.cells, cl, mc.cores=mc.cores)$pred.df
    pred.cl= setNames(pred.df$pred.cl, row.names(pred.df))
    cl = c(setNames(as.character(cl), names(cl)), setNames(as.character(pred.cl), names(pred.cl)))
  }
  cl.platform.counts = table(comb.dat$meta.df[names(cl), "platform"],cl)
  #print(cl.platform.counts)
  ###If a cluster is not present in reference sets, split the cells based on imputed cluster based on cells in reference set.
  ref.de.param.list = comb.dat$de.param.list[ref.sets]
  cl.min.cells = sapply(ref.de.param.list, function(x)x$min.cells)
  cl.big= cl.platform.counts[ref.sets,,drop=F] >= cl.min.cells
  bad.cl = colnames(cl.big)[colSums(cl.big) ==0]
  cl.big = setdiff(colnames(cl.big), bad.cl)
  if(length(cl.big)<=1){
    return(NULL)
  }
  if(length(bad.cl) > 0){
    print("Bad.cl")
    print(bad.cl)
    tmp.cells = names(cl)[cl %in% bad.cl]
    pred.prob = predict_knn(knn.comb[tmp.cells,,drop=F], comb.dat$all.cells, cl)$pred.prob
    pred.prob$freq=NULL
    pred.prob = pred.prob %>% group_by(query) %>% filter(!nn.cl %in% bad.cl) %>% mutate(freq = n/sum(n))
    pred.df = pred.prob %>% group_by(query) %>% summarize(pred.cl = nn.cl[which.max(freq)], pred.score = max(freq))
    cl[pred.df$query]= pred.df$pred.cl
  }
  cl= merge_cl_multiple(comb.dat=comb.dat,merge.sets=merge.sets, cl=cl, anchor.genes=select.genes,mc.cores=mc.cores)  
  if(length(unique(cl))<=1){
    return(NULL)
  }
  print(table(cl))
  result$cl = cl
  result$select.genes= select.genes
  result$ref.de.param.list = ref.de.param.list
  return(result)
}
 
harmonize <- function(comb.dat, prefix, overwrite=TRUE, dir="./",...)
  {
    fn = file.path(dir, paste0(prefix, ".rda"))
    print(fn)
    if(!overwrite){
      if(file.exists(fn)){
        load(fn)
        return(result)
      }
    }
    result = knn_joint(comb.dat, ...)
    save(result, file=fn)
    if(is.null(result)){
      return(NULL)
    }
    print("Cluster size")
    print(table(result$cl))
    #g = plot_cl_meta_barplot(result$cl, meta.df[names(result$cl), "platform"])
    #g = g + theme(axis.text.x = element_text(angle=45,hjust=1, vjust=1))
    #ggsave(paste0(prefix, ".platform.barplot.pdf"),g,height=5, width=12)
    #plot_confusion(result$cl, prefix,comb.dat)
    return(result)
  }
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param comb.dat 
##' @param select.cells 
##' @param prefix 
##' @param result 
##' @param ... 
##' @return 
##' @author Zizhen Yao
i_harmonize<- function(comb.dat, select.cells=comb.dat$all.cells, ref.sets=names(comb.dat$dat.list), select.sets=names(comb.dat$dat.list), prefix="", joint.rd.dat=NULL, result=NULL, overwrite=TRUE, sample.size = 50000, dir="./",rm.genes=NULL,split.size = 100, ...)
  {
    #attach(comb.dat)
    if(is.null(result)){
      result = harmonize(comb.dat=comb.dat, select.cells=select.cells, ref.sets=ref.sets, joint.rd.dat=joint.rd.dat, select.sets=select.sets,prefix=prefix, overwrite=overwrite,sample.size=sample.size, dir=dir,rm.genes=rm.genes,...)
    }
    if(is.null(result)){
      return(NULL)
    }
    all.results= list(result)
    names(all.results) = prefix
    cl = result$cl
    if(length(cl) < split.size){
      return(all.results)
    }
    for(i in as.character(sort(unique(result$cl)))){
      tmp.prefix=paste(prefix, i,sep=".")
      print(tmp.prefix)
      select.cells= names(cl)[cl == i]
      platform.size = table(comb.dat$meta.df[select.cells, "platform"])      
      print(platform.size)
      pass.th = sapply(select.sets, function(set)platform.size[[set]] >= comb.dat$de.param.list[[set]]$min.cells)
      pass.th2 = sapply(ref.sets, function(set)platform.size[[set]] >= comb.dat$de.param.list[[set]]$min.cells*2)
      if(sum(pass.th) > 1 & sum(pass.th[ref.sets]) == length(ref.sets) & sum(pass.th2) >= 1){
        computed=FALSE
        fn = file.path(dir, paste0(tmp.prefix, ".rda"))
        print(fn)
        result=NULL
        if(!overwrite){
          if(file.exists(fn)){
            load(fn)
            computed=TRUE
          }
        }
        if(comb.dat$type!="big"){
          print("mem")
          tmp.result = i_harmonize(comb.dat, select.cells=select.cells, ref.sets=ref.sets, joint.rd.dat=joint.rd.dat, select.sets=select.sets,prefix=tmp.prefix, overwrite=overwrite, sample.size=sample.size,dir=dir,result=result,...)
        }
        else{
          if(!computed & length(select.cells) < sample.size){
            new.comb.dat = comb.dat
            dat.list = sample_cl_dat(comb.dat, sets=select.sets, cl=cl[select.cells],sample.size)
            if(!is.null(rm.genes)){
              dat.list = rm_genes(dat.list, rm.genes)
            }
            new.comb.dat$dat.list = dat.list
            new.comb.dat$type="mem"
            print("mem")
            tmp.result = i_harmonize(comb.dat=new.comb.dat, select.cells=select.cells, ref.sets=ref.sets, prefix=tmp.prefix, joint.rd.dat=joint.rd.dat, overwrite=overwrite, sample.size=sample.size,dir=dir,result=result,...)
            rm(new.comb.dat)
            #gc()
          }
          else{
            print("big")
            tmp.result = i_harmonize(comb.dat, select.cells=select.cells, ref.sets=ref.sets, prefix=tmp.prefix, joint.rd.dat = joint.rd.dat, overwrite=overwrite, sample.size=sample.size, dir=dir,result=result,rm.genes=rm.genes,...)
          }
        }
        if(!is.null(tmp.result)){
          all.results[names(tmp.result)] = tmp.result
        }
      }
    }
    return(all.results)
  }

rm_genes <- function(dat.list, rm.genes)
  {
    for(x in names(dat.list)){
      dat = dat.list[[x]]
      rm = row.names(dat) %in% rm.genes
      if(any(rm)){
        dat = dat[!rm,,drop=FALSE]
        dat.list[[x]] = dat
      }
    }
    return(dat.list)
  }


get_de_result_big <- function(comb.dat, cl, cl.bin, sets=names(comb.dat$dat.list),cl.stats.list=NULL, de.d = "de_parquet",de.sum.d = "de_summary",pairs.fn="pairs.parquet",block.size=10000,mc.cores=25,...)
  {
    de.result=list()
    if(is.null(cl.stats.list)){
      cl.stats.list = get_cl_stats_list(comb.dat, merge.sets=sets, cl=cl,mc.cores=mc.cores)
    }
    cl.means.list = cl.stats.list$cl.means.list
    cl.present.list = cl.stats.list$cl.present.list
    cl.sqr.means.list = cl.stats.list$cl.sqr.means.list
    
    if(file.exists(pairs.fn)){
      pairs.df=read_parquet(pairs.fn)
    }
    else{
      pairs.df = as.data.frame(create_pairs(unique(cl)))
      pairs.df$pair_id = 1:nrow(pairs.df)
      pairs.df$pair = row.names(pairs.df)
      pairs.df$pair_bin = ceiling(pairs.df$pair_id/block.size)
      write_parquet(pairs.df, pairs.fn)
    }
    dir.create(de.d)
    dir.create(de.sum.d)
    
    for(x in sets){
      print(x)
      tmp.cl = cl[names(cl) %in% (comb.dat$dat.list[[x]]$col_id) & cl %in% colnames(cl.means.list[[x]])]
      if(is.factor(tmp.cl)){
        tmp.cl= as.factor(tmp.cl)
      }
      tmp.de.d = file.path(de.d, paste0("set=",x))
      tmp.de.sum.d = file.path(de.sum.d, paste0("set=",x))
      de.result[[x]] = de_selected_pairs(norm.dat = NULL, cl = tmp.cl, cl.bin=cl.bin, pairs = pairs.df, de.param = comb.dat$de.param.list[[x]], cl.means= cl.means.list[[x]], cl.present=cl.present.list[[x]], cl.sqr.means=cl.sqr.means.list[[x]], mc.cores=mc.cores, out.dir=tmp.de.d, summary.dir = tmp.de.sum.d, return.summary=TRUE,...)
    }
    return(de.result)
  }


comb_de_result_big <- function(ds, cl.means.list, cl.bin, sets=names(cl.means.list), max.num=1000, lfc.conservation.th=0.7, mc.cores=20,out.dir="comb_de_parquet",overwrite=FALSE)
{
  library(data.table)  
  library(dplyr)
  require(doMC)
  require(foreach)
  registerDoMC(cores=mc.cores)
  if(!dir.exists(out.dir)){
    dir.create(out.dir)
  }
  cl.means.list = sapply(cl.means.list, as.matrix, simplify=FALSE)
  all.cl = unique(unlist(lapply(cl.means.list,colnames)))
  all.gene = unique(unlist(lapply(cl.means.list,row.names)))
  all.bins = unique(cl.bin$bin)
  tmp=foreach(bin1 = all.bins,.combine="c")%:%
    foreach(bin2 = all.bins,.combine="c")%dopar% {
      if(!overwrite & dir.exists(file.path(out.dir, paste0("bin.x=",bin1), paste0("bin.y=",bin2)))){
        return(NULL)
      }        
      de.df = droplevels(ds %>% filter(bin.x==bin1 & bin.y==bin2) %>% collect())     
      de.genes = de.df %>% group_by(gene,P1,P2) %>% summarize(num=n(),logPval=mean(logPval))
      pairs = de.df %>% select(pair, P1, P2, bin.x, bin.y) %>% distinct()
      lfc = list()
      for(set in sets){
        cl.means= cl.means.list[[set]]
        missing.cl = setdiff(all.cl,colnames(cl.means))
        missing.gene = setdiff(all.gene,row.names(cl.means))
        select = with(de.genes,!(P1 %in% missing.cl | P2 %in% missing.cl| gene %in% missing.gene))
        exp1=get_pair_matrix(cl.means, de.genes$gene[select], de.genes$P1[select])
        exp2=get_pair_matrix(cl.means, de.genes$gene[select], de.genes$P2[select])
        lfc[[set]] = rep(NA, nrow(de.genes))
        lfc[[set]][select] = exp1 - exp2
      }
      lfc = do.call("cbind",lfc)
      denom = rowSums(!is.na(lfc))
      up.ratio = rowSums(lfc > 1, na.rm=TRUE)/denom      
      lfc.mean = rowMeans(lfc, na.rm=TRUE)    
      de.genes = cbind(de.genes, data.frame(up.ratio, lfc=lfc.mean))
###dplyr filter is slow for some reason
      de.genes = with(de.genes,de.genes[up.ratio > lfc.conservation.th ,,drop=FALSE])
      de.genes = de.genes %>% arrange(P1,P2, -num, -abs(lfc)) %>% group_by(P1,P2) %>% mutate(rank=1:n())
      de.genes = de.genes %>% left_join(pairs, by=c("P1","P2"))
      write_dataset(de.genes,  out.dir, partition=c("bin.x","bin.y"))
    }
}

 
get_de_result_recursive <- function(comb.dat, all.results, sets=names(comb.dat$dat.list),ref.dat.list, max.cl.size = 300, ...)
  {
    #impute.dat.list <<- list()
    for(x in names(all.results)){
      print(x)
      result = all.results[[x]]
      cl = result$cl
      cl = cl[names(cl) %in% colnames(ref.dat.list)]
      cl = cl[sample_cells(cl, 300)]      
      de.result = get_de_result(ref.dat.list, comb.dat$de.param.list, cl = cl)
      cl.means.list = get_cl_means_list(ref.dat.list, comb.dat$de.param.list, cl=cl, sets=sets)
      de.result$comb.de.genes = comb_de_result(de.result$de.genes.list, cl.means.list = cl.means.list, common.genes=comb.dat$common.genes, ...)
      all.results[[x]]$de.result = de.result
    }
    return(all.results)
  }
        


impute_knn_cross <- function(comb.dat, split.results, impute.dat.list, ref.sets, select.cells, init=TRUE)
  {
###cross-modality Imputation based on nearest neighbors in each iteraction of clustering using anchoring genes or genes shown to be differentiall expressed.
    select.genes = row.names(impute.dat.list[[1]])
    for(x in names(split.results)){
      print(x)
      result = split.results[[x]]
      cl = result$cl
      knn = result$knn    
      for(ref.set in ref.sets){
        if(ref.set %in% names(result$ref.list)){
          tmp.cells = row.names(result$knn)
          query.cells = intersect(tmp.cells[comb.dat$meta.df[tmp.cells,"platform"] != ref.set], select.cells)
          if(init & any(!query.cells %in% row.names(impute.dat.list[[ref.set]]))){
            impute.genes = select.genes
          }
          else{
            impute.genes=intersect(select.genes,c(result$select.markers, result$select.genes))
          }
          select.cols = comb.dat$meta.df[comb.dat$all.cells[result$knn[1,]],"platform"] == ref.set
          if(sum(select.cols)==0){
            next
          }
          if(length(query.cells)==0){
            next
          }          
          select.knn = result$knn[query.cells,select.cols,drop=F]
          dat = impute.dat.list[[ref.set]]
          gene.id = match(impute.genes, row.names(dat))
          cell.id = match(query.cells, colnames(dat))
          reference.id = match(comb.dat$all.cells, colnames(dat))
          ImputeKnn(select.knn, reference.id, cell.id, gene_idx=gene.id, 
                      dat=dat,impute_dat=dat, w_mat_= NULL,
                      transpose_input=FALSE, transpose_output=FALSE)                  
        }
      }
    }
    return(impute.dat.list=impute.dat.list)
  }




#### assume within data modality have been performed
####
impute_knn_global_big<- function(comb.dat, split.results, select.genes, select.cells, ref.dat.list, ref.sets=names(ref.dat.list), sets=comb.dat$sets, rm.eigen=NULL, rm.th=0.7, verbose=FALSE,mc.cores=5, org.rd.dat.list=NULL, blocksize=100000)
  {

    knn.list <- list()    
    impute.dat.list <- list()
    
    ###Impute the reference dataset in the original space globally
    if(is.null(org.rd.dat.list)){
      org.rd.dat.list <- list()
      for(x in ref.sets){
          print(x)
          ref.dat = ref.dat.list[[x]][select.genes, ]
          tmp.cells = intersect(select.cells, comb.dat$dat.list[[x]]$col_id)
          rd.result <- rd_PCA_big(comb.dat$dat.list[[x]], ref.dat, select.cells=tmp.cells, max.dim=100, th=0.5, mc.cores=mc.cores,method="elbow",verbose=verbose)
          rd.dat = rd.result$rd.dat
          if(!is.null(rm.eigen)){
            rd.dat = filter_RD(rd.dat, rm.eigen, rm.th, verbose=verbose)
          }
          ref.cells=colnames(ref.dat)
          knn = get_knn_batch(rd.dat, rd.dat[ref.cells,], method="Annoy.Euclidean", mc.cores=mc.cores, batch.size=50000,k=k, transposed=FALSE)                  
          reference.id = 1:length(ref.cells)
          gene.id = 1:length(select.genes)
          impute.big.dat = create_big.dat(select.genes, row.names(rd.dat))
          bin = ceiling(seq_len(nrow(rd.dat))/blocksize)
          for(i in 1:max(bin)){
            select= bin==x
            tmp.cells = row.names(rd.dat)[select]
            cell.id = 1:length(tmp.cells)
            impute.dat = matrix(0, nrow=length(select.genes), ncol=length(cell.id))
            dimnames(impute.dat) = list(select.genes, tmp.cells)
            dat = as.matrix(ref.dat)
            ImputeKnn(knn[select,], reference.id, cell.id, gene.id, dat=dat, impute.dat,
                      w_mat_ = NULL,transpose_input=FALSE, transpose_output=FALSE)
            impute.big.dat$FBM[,select] = impute.dat            
          }          
          impute.dat.list[[x]] = impute.big.dat
        }
    }
    impute.dat.list = impute_knn_cross(comb.dat, split.results, impute.dat.list, ref.sets, init=TRUE)
  }



gene_gene_cor_conservation <- function(dat.list, select.genes, select.cells,pairs=NULL)
  {
    sets = names(dat.list)
    gene.cor.list = sapply(sets, function(set){
      print(set)
      dat = dat.list[[set]]
      gene.cor = cor(t(as.matrix(dat[select.genes,intersect(colnames(dat),select.cells)])))
      gene.cor[is.na(gene.cor)] = 0
      gene.cor
    },simplify=F)
    if(is.null(pairs)){
      n.sets = length(sets)	
      pairs = cbind(rep(sets, rep(n.sets,n.sets)), rep(sets, n.sets))
      pairs = pairs[pairs[,1]<pairs[,2],,drop=F]
    }
    gene.cor.mat= sapply(1:nrow(pairs), function(i){
      p = pairs[i,]
      print(p)
      pair_cor(gene.cor.list[[p[1]]], gene.cor.list[[p[2]]])
    })
    colnames(gene.cor.mat) = paste0(pairs[,1],":",pairs[,2])
    return(gene.cor.mat)
  }

#' Sim knn
#'
#' @param sim 
#' @param k 
#'
#' @return
#' @export
#'
#' @examples
sim_knn <- function(sim, k=15)
{
  require(matrixStats)
  th =  rowOrderStats(as.matrix(sim), which=ncol(sim)-k+1)
  select = sim >= th
  knn.index = t(apply(select, 1, function(x)head(which(x),k)))
  if(k==1){
    knn.index= matrix(knn.index, ncol=1)
  }
  knn.distance = do.call("rbind",lapply(1:nrow(sim), function(i) (1- sim[i,,drop=F])[knn.index[i,,drop=F]]))
  return(list(knn.index, knn.distance))
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param knn.index 
##' @param bin_size 
##' @param prune 
##' @return 
##' @author Zizhen Yao
jaccard_big <-  function(knn.index, bin_size =50000, prune=0)
{
  library(arrow)
  library(dplyr)
  library(Matrix)
  k = ncol(knn.index)
  total = nrow(knn.index)
  id=row.names(knn.index)
  max.index=max(knn.index)
  ## common values:
  bins = ceiling((1:nrow(knn.index))/bin_size)
  nbin = max(bins)
  fn.list = c()
  for(b1 in 1:nbin){
    for(b2 in b1:nbin){
      i1 = which(bins==b1)
      j1 = as.vector(knn.index[i1,])
      ii1 = rep(i1, k)
      knn.mat1 = sparseMatrix(i = j1, j=ii1, x=1, dims=c(max.index, total))
      i2 = which(bins==b2)
      j2 = as.vector(knn.index[i2,])
      ii2 = rep(i2, k)
      knn.mat2 = sparseMatrix(i = j2, j=ii2, x=1, dims=c(max.index, total))
      A <-  Matrix::crossprod(knn.mat1, knn.mat2)
      A@x = A@x / (2*k - A@x)
      #A should be a sym matrix
      A <- as(A, "TsparseMatrix")
      if(nbin>1){
        df = data.frame(from=A@i+1, to=A@j+1, weight=A@x)
        df = df %>% filter(weight > prune & from <= to)
        fn = tempfile()
        write_parquet(df, sink=fn)
        fn.list = c(fn.list, fn)
      }
      #rm(A)
      #rm(df)
      #gc()
    }
  }
  if(nbin>1){
    df = open_dataset(fn.list) %>% collect()
    unlink(fn.list)
    A = sparseMatrix(i=df[[1]],j=df[[2]],x=df[[3]], repr="T",symmetric=TRUE)    
  }
  return(A)
}

#' KNN Jaccard Louvain
#'
#' @param knn.index 
#'
#' @return
#' @export
#'
#' @examples
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param knn.index 
##' @param method 
##' @param prune 
##' @return 
##' @author Zizhen Yao
knn_jaccard_clust <- function(knn.index, method=c("leiden","louvain"),prune=0.05,return.graph=FALSE,...)
  {
    require(igraph)
    cat("Get jaccard\n")
    #sim=knn_jaccard(knn.index,...)
    #sim = ComputeSNN(knn.index,prune=prune)
    sim=jaccard_big(knn.index, prune=prune)
    rownames(sim) = colnames(sim) = row.names(knn.index)

    gr <- igraph::graph.adjacency(sim, mode = "upper", 
                                  weighted = TRUE)
    if(method[1]=="louvain"){
      cat("Louvain clustering\n")
      result <- igraph::cluster_louvain(gr,...)
    }
    else{
      cat("Leiden clustering\n")
      library(leidenAlg)
      result <- leiden.community(gr,...)      
    }
    result$cl=membership(result)
    rm(sim)
    if(return.graph){
      result$gr = gr
    }
    else{
      rm(gr)
    }
    #gc()
    return(result)
  }


#' combine cl
#'
#' @param all.results 
#'
#' @return
#' @export
#'
#' @examples
combine_cl <- function(all.results)
  {
    cl = all.results[[1]]$cl
    cl = setNames(as.integer(cl),names(cl))
    markers=all.results[[1]]$markers
    n.cl = max(cl)
    for(i in 2:length(all.results)){
      #if(is.null(all.results[[i]]$cl) | length(unique(all.results[[i]]$cl)) < 2) next
      if(is.null(all.results[[i]]$cl)) next
      new.cl = all.results[[i]]$cl
      new.cl = setNames(as.integer(new.cl)+ n.cl,names(new.cl))
      cl[names(new.cl)] = new.cl
      n.cl = max(cl)
      cat(names(all.results)[i], n.cl, "\n")
    }
    return(cl)
  }



#' Batch process
#'
#' @param x 
#' @param batch.size 
#' @param FUN 
#' @param mc.cores 
#' @param .combine 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
batch_process <- function(x, batch.size, FUN, mc.cores=1, .combine="c",bins=NULL,...)
  {
    require(foreach)
    require(doMC)
    if(is.null(bins)){
      bins = split(x, floor((1:length(x))/batch.size))
    }
    mc.cores=min(mc.cores, length(bins))
    registerDoMC(cores=mc.cores)
    results= foreach(i=1:length(bins), .combine=.combine) %dopar% {
      FUN(bins[[i]],...)
    }
    return(results)
  }




#' Title
#'
#' @param cl.means.list 
#'
#' @return
#' @export
#'
#' @examples
get_gene_cl_correlation <- function(cl.means.list)
  {
    sets=names(cl.means.list)
    gene.cl.cor = list()
    for(i in 1:(length(cl.means.list)-1)){
      for(j in (i+1):length(cl.means.list)){
        pair= paste(sets[i], sets[j], sep=":")
        common.cl = intersect(colnames(cl.means.list[[i]]), colnames(cl.means.list[[j]]))
        common.genes = intersect(row.names(cl.means.list[[i]]), row.names(cl.means.list[[j]]))
        gene.cor =  pair_cor(cl.means.list[[i]][common.genes,common.cl],cl.means.list[[j]][common.genes,common.cl])
        gene.cl.cor[[pair]] = gene.cor
      }
    }
    return(gene.cl.cor)
  }


#' Title
#'
#' @param consensus.cl 
#' @param prefix 
#' @param comb.dat 
#' @param consensus.cl.df 
#' @param do.droplevels 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_confusion <- function(consensus.cl, prefix, comb.dat,consensus.cl.df = NULL, sets=names(comb.dat$cl.list), do.droplevels = FALSE,cex.x=6, cex.y=6,...)
{
  g.list=list()
  for(x in sets){
    if(sum(names(comb.dat$cl.list[[x]]) %in% names(consensus.cl)) > 0){
      if(is.null(consensus.cl.df)){
        g = compare_annotate(consensus.cl, comb.dat$cl.list[[x]], comb.dat$cl.df.list[[x]], rename = FALSE, do.droplevels=do.droplevels,cex.x=cex.x, cex.y=cex.y)$g
        g = g + xlab("consensus cluster") + ylab(x)
      }
      else{
        g = compare_annotate(comb.dat$cl.list[[x]], consensus.cl, consensus.cl.df, rename = FALSE, do.droplevels=do.droplevels)$g
        g = g + ylab("consensus cluster") + xlab(x)
      }
      g.list[[x]] <- g
      ggsave(paste(prefix, x, "pdf", sep="."), g,...)
    }
  }
  return(g.list)
}



