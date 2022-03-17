rd_PCA_big <- function(big.dat, dat, select.cells=big.dat$col_id, max.dim=10, th=2, verbose=TRUE, mc.cores=1,method="zscore")
{
  system.time({tmp = get_PCA(dat, max.pca=max.dim, verbose=verbose,method=method,th=th)})
  if(is.null(tmp)){
    return(NULL)
  }
  rot = tmp$rot
  pca = tmp$pca
  rd.dat = tmp$rd.dat
  if(ncol(dat)< length(select.cells)){
    if(verbose){
      print("project")
    }
    system.time({rd.dat = big_project(big.dat, select.cells, rot, mc.cores=mc.cores)})
  }
  rm(rot)
  rm(dat)
  return(list(rd.dat=rd.dat, pca=pca))
}


big_project <- function(big.dat, select.cells, rot, mc.cores=1,...)
  {
    require(Matrix)
    my_project <- function(big.dat, cols){
      fn = tempfile()
      dat = get_logNormal(big.dat, cols)[row.names(rot),]
      dat = dat[row.names(rot),]
      dat = Matrix::crossprod(dat, rot)
      df= as.data.frame(as.matrix(dat))
      df$gene = row.names(df)      
      write_parquet(df, sink=fn)
      rm(df)
      rm(dat)
      gc()
      fn
    }
    rd.dat.fn = big_dat_apply(big.dat, select.cells, FUN=my_project, .combine="c",  mc.cores=mc.cores,...)
    library(arrow)
    df = open_dataset(rd.dat.fn) %>% collect()
    rd.dat = as.matrix(df[,1:(ncol(df)-1)])
    row.names(rd.dat)=df$gene
    unlink(rd.dat.fn)
    return(rd.dat)
  }



#' One round of clustering in the iteractive clustering pipeline 
#'
#' @param norm.dat normalized expression data matrix in log transform, using genes as rows, and cells and columns. Users can use log2(FPKM+1) or log2(CPM+1).
#' @param select.cells The cells to be clustered. Default: columns of norm.dat
#' @param counts Raw gene counts. Default NULL, inferred from norm.dat.
#' @param method Clustering method. It can be "louvain", "hclust" and "kmeans". Default "louvain"
#' @param vg.padj.th High variance gene adjusted pvalue cut off. Default 0.5.
#' @param dim.method Dimension reduction techniques. Current options include "pca" and "WGCNA". Default "pca"
#' @param max.dim The number of top dimensions retained. Default 20. Since clustering is performed iteratively, not all relevant dimensions need to be captured in one iterations. 
#' @param rm.eigen The reduced dimensions that need to be masked and removed. Default NULL.  
#' @param rm.th The cutoff for correlation between reduced dimensions and rm.eigen. Reduced dimensions with correlatin with any rm.eigen vectors are not used for clustering. Default 0.7
#' @param de.param The differential gene expression threshold. See de_param() function for details. 
#' @param min.genes The minimal number of high variance and differentially expressed genes genes. Default 5. 
#' @param type Can either be "undirectional" or "directional". If "undirectional", the differential gene threshold de.param is applied to combined up-regulated and down-regulated genes, if "directional", then the differential gene threshold is applied to both up-regulated and down-regulated genes. 
#' @param maxGenes Only used when dim.method=="WGCNA". The maximum number of genes to calculate gene modules. 
#' @param sampleSize The number of sampled cells to compute reduced dimensions.
#' @param max.cl.size Sampled cluster size. This is to speed up limma DE gene calculation. Instead of using all cells, we randomly sampled max.cl.size number of cells for testing DE genes.    
#' @param prefix Used to keep track of intermediate results in "verbose" mode. Default NULL.
#' @param verbose Default FALSE
#'
#' @return Clustering result is returned as a list with two elements: 
#'         cl: cluster membership for each cell
#'         markers: top markers that seperate clusters     
#'         
onestep_clust_big<- function(big.dat, 
                             select.cells= big.dat$col_id,
                             genes.allowed = big.dat$row_id,
                             counts = NULL,
                             method = c("louvain","leiden","ward.D", "kmeans"),
                             vg.padj.th = 0.5,                             
                             dim.method = c("pca","WGCNA"), 
                             max.dim = 20,
                             mc.cores=20,
                             rm.eigen = NULL, 
                             rm.th = 0.7, 
                             de.param = de_param(),
                             merge.type = c("undirectional", "directional"), 
                             maxGenes = 3000,
                             sampleSize = 50000,
                             jaccard.sampleSize = 300000,
                             max.cl.size = 300,
                             k.nn=15,
                             prefix = NULL, 
                             verbose = FALSE)
                            
  {
    library(matrixStats)    
    method=method[1]
    merge.type=merge.type[1]

    
    sampled = sample(select.cells, min(sampleSize, length(select.cells)))
    if(length(sampled) > length(select.cells)/2){
      sampled = select.cells
    }
    if(is.null(counts) & sum(!sampled %in% colnames(counts))>0){
      counts = get_counts(big.dat, sampled,sparse=TRUE,mc.cores=mc.cores)
    }
    else{
      counts = counts[,sampled]
    }
    if(verbose){
      print("Find high variance genes")
    }    
    system.time({vg = findVG(counts)})
    select.genes = with(vg, row.names(vg)[which(loess.padj < vg.padj.th | dispersion >2.5)])
    select.genes = intersect(select.genes, genes.allowed)
    if(length(select.genes) < de.param$min.genes){
      return(NULL)
    }
    select.genes = head(select.genes[order(vg[select.genes, "loess.padj"],-vg[select.genes, "z"])],maxGenes)
    if(verbose){
      print("Reduce dimensions")
    }
    counts = logCPM(counts)[select.genes,]
    rd.result = rd_PCA_big(big.dat, dat=counts, select.cells, max.dim=max.dim, verbose=verbose,mc.cores=mc.cores,method="elbow")
    if(is.null(rd.result)){
      return(NULL)
    }
    rd.dat = rd.result$rd.dat
    rm(counts)
    if(!is.null(rm.eigen)){
      
      rd.dat <- filter_RD(rd.dat, rm.eigen, rm.th, verbose=verbose)
    }
    if(is.null(rd.dat)||ncol(rd.dat)==0){
      return(NULL)
    }
    if(verbose){
      print("Clustering")
    }
    max.cl = ncol(rd.dat) *2 + 1
    if(method %in% c("louvain","leiden")){
      k = pmin(k.nn, round(nrow(rd.dat)/2))
      if(verbose){
        print("Compute KNN")
      }
      jaccard.sampled=row.names(rd.dat)
      if(length(jaccard.sampled)> jaccard.sampleSize){        
        jaccard.sampled = sample(jaccard.sampled, jaccard.sampleSize)
      }
      ref.rd.dat = rd.dat[jaccard.sampled,]
      library(BiocNeighbors)
      index = buildAnnoy(ref.rd.dat, distance ="Euclidean", transposed = FALSE, ntrees=100)
      
      knn.result=get_knn_batch(dat=rd.dat, ref.dat = ref.rd.dat, k=k, method = "Annoy.Euclidean", mc.cores=mc.cores, batch.size = 10000,transposed=FALSE,index=index, return.distance=TRUE)
      knn.index=knn.result[[1]]
      knn.dist=knn.result[[2]]
      if(verbose){
        print("Jaccard Cluster")
      }
      tmp=knn_jaccard_clust(knn.index[jaccard.sampled,], method=method,prune=0.05)
      #tmp = jaccard_louvain(rd.dat, k)
      if(is.null(tmp)){
        return(NULL)
      }
      cl = tmp$cl
      if(length(cl) < nrow(knn.index)){
        diff.cells = setdiff(row.names(knn.index), names(cl))
        pred.df = predict_knn(knn.index[diff.cells,,drop=F], jaccard.sampled, cl, mc.cores=mc.cores)$pred.df
        pred.cl= setNames(pred.df$pred.cl, row.names(pred.df))
        cl = c(setNames(as.character(cl), names(cl)), setNames(as.character(pred.cl), names(pred.cl)))
      }
      if(length(unique(cl))>max.cl){
        tmp.means = get_cl_means(rd.dat, cl)
        tmp.hc = hclust(dist(t(tmp.means)), method="average")
        tmp.cl= cutree(tmp.hc, pmin(max.cl, length(unique(cl))))
        cl = setNames(tmp.cl[as.character(cl)], names(cl))
      }
    }
    else if(method=="ward.D"){
      hc = hclust(dist(rd.dat),method="ward.D")
      #print("Cluster cells")
      cl = cutree(hc, max.cl)
    }
    else if(method=="kmeans"){
      cl = kmeans(rd.dat, max.cl)$cluster
    }
    else{
      stop(paste("Unknown clustering method", method))
    }
    #sampled.cells = sample_cells(cl, max.cl.size)
    #norm.dat = get_logNormal(big.dat, sampled.cells)
    #merge.result=merge_cl(norm.dat, cl=cl, rd.dat=rd.dat, merge.type=merge.type, de.param=de.param, max.cl.size=max.cl.size, verbose=verbose)
    if(verbose){
      print("merge cluster")
    }
    merge.result=merge_cl_big(big.dat, cl=cl, rd.dat=rd.dat, merge.type=merge.type, de.param=de.param, max.cl.size=max.cl.size, verbose=verbose, return.markers=TRUE,mc.cores=mc.cores)
    #rm(norm.dat)
    gc()
    if(is.null(merge.result)){
      return(NULL)
    }
    cl = merge.result$cl
    de.genes = merge.result$de.genes
    markers= merge.result$markers
    result=list(cl=cl, markers=markers)
    rm(rd.dat.t)
    rm(merge.result)
    rm(counts)
    gc()
    return(result)
  }



#' Iterative clustering algorithm for single cell RNAseq dataset
#'
#' @param norm.dat normalized expression data matrix in log transform, using genes as rows, and cells and columns. Users can use log2(FPKM+1) or log2(CPM+1)
#' @param select.cells The cells to be clustered
#' @param prefix The character string to indicate current iteration.
#' @param split.size The minimal cluster size for further splitting
#' @param result The current clustering result as basis for further splitting.
#' @param method Clustering method. It can be "auto", "louvain", "hclust"
#' @param ... Other parameters passed to method `onestep_clust()`
#'
#' @return Clustering result is returned as a list with two elements: 
#'         cl: cluster membership for each cell
#'         markers: top markers that seperate clusters     
#'         
#' @examples clust.result <- iter_clust(norm.dat)
#'           clust.result <- iter_clust(norm.dat, de.param = de_param(q1.th = 0.5, de.score.th = 100))
iter_clust_big<- function(big.dat=NULL,
                          select.cells = big.dat$col_id,
                          prefix = NULL, 
                          split.size = 10, 
                          result = NULL,
                          method = "auto",
                          counts = NULL,
                          sampleSize = 10000,
                          mc.cores=10,
                          overwrite=TRUE,
                          verbose=FALSE,
                          ...)
  {
    if(!is.null(prefix)) {
      cat(prefix, length(select.cells),"\n")
    }
    
    if(method == "auto"){
      if(length(select.cells) > 2000){
        select.method="louvain"
      }
      else{
        select.method="ward.D"
      }
    }
    else{
      select.method=method
    }

    if(is.null(result)){
      outfile=paste0(prefix, ".rda")
      if(file.exists(outfile) & !overwrite){
        load(outfile)       
      }
      else{
        result=onestep_clust_big(big.dat=big.dat, select.cells=select.cells, prefix=prefix,method=select.method, counts=counts, sampleSize= sampleSize,mc.cores=mc.cores,verbose=verbose,...)
        if(verbose){
          save(result, file=outfile)
        }
        gc()
      }
      if(is.null(result)){
        return(NULL)
      }
    }
    cl = result$cl[select.cells]
    gene.mod = result$gene.mod
    markers=result$markers
    cl = setNames(as.integer(cl),names(cl))
    new.cl =cl
    cl.size = table(cl)
    to.split = names(cl.size)[cl.size >=split.size]
    if(length(to.split)>0){
      n.cl = 1
      for(x in sort(unique(cl))){
        tmp.cells = names(cl)[cl==x]
        if(!x %in% to.split){
          new.cl[tmp.cells]=n.cl
        }
        else{
          tmp.prefix = paste(prefix, x, sep=".")
          if(length(tmp.cells) < 50000){
            norm.dat = get_logNormal(big.dat, tmp.cells)
            tmp.result=iter_clust(norm.dat=norm.dat, select.cells=tmp.cells, prefix=tmp.prefix,split.size=split.size,method= method, counts=counts, sampleSize=sampleSize, overwrite=overwrite,verbose=verbose,...)
            rm(norm.dat)
          }
          else{
            tmp.result=iter_clust_big(big.dat=big.dat, select.cells=tmp.cells, prefix=tmp.prefix,split.size=split.size,method= method, counts=counts, sampleSize=sampleSize, overwrite=overwrite,mc.cores=mc.cores,verbose=verbose,...)
          }
          gc()
          if(is.null(tmp.result)){
            new.cl[tmp.cells]=n.cl
          }
          else{
            tmp.cl = tmp.result$cl
            if(length(unique(tmp.cl)>1)){
              new.cl[names(tmp.cl)] = n.cl + as.integer(tmp.cl)
              markers=union(markers, tmp.result$markers)
            }
          }
        }
        n.cl = max(new.cl)+1
      }
      cl = new.cl
    }
    result=list(cl=cl, markers=markers)
    return(result)
  }



iter_clust_merge_big <- function(big.dat, select.cells=big.dat$col_id, merge.type="undirectional", de.param = de_param(), max.cl.size = 300,...)
{
  result <- iter_clust_big(big.dat=big.dat, select.cells=select.cells, de.param = de.param, merge.type=merge.type, ...)
  cl = result$cl
  markers = result$markers
  tmp.cells = sample_cells(cl, max.cl.size)
  norm.dat = get_logNormal(big.dat, tmp.cells)
  merge.result= merge_cl(norm.dat=norm.dat, cl=cl, rd.dat.t=norm.dat[markers,], de.param = de.param, merge.type=merge.type, return.markers=FALSE)
  return(merge.result)
}



get_cols_delayedArray <- function(big.dat_delayedArray, cols, keep.col=TRUE, sparse=TRUE)
{
  library(Matrix)
  if(is.character(cols)){
    id = match(cols, colnames(big.dat_delayedArray))
  }
  else{
    id = cols
  }
  ord = order(id)
  
  mat = big.dat_delayedArray[,id[ord],drop=F]
  if(keep.col){
    org.order = (1:length(id))[order(ord)]    
    mat = mat[,org.order,drop=F]
    colnames(mat) = colnames(big.dat_delayedArray)[id]      
  }
  else{
    colnames(mat) = colnames(big.dat_delayedArray)[id[ord]]      
  }
  if(sparse){
    mat = Matrix(mat,sparse=TRUE)
  }
  rownames(mat) = rownames(big.dat_delayedArray)
  return(mat)
}



get_cl_stats_fbm <- function(big.dat, cl, max.cells=100000, stats=c("means"))
  {
    if(!is.factor(cl)){
      cl = as.factor(cl)
    }
    cl.size = table(cl)
    cl.bins=round(cumsum(cl.size) / max.cells)
    cl.bins=split(names(cl.size), cl.bins)
    tmp.results=sapply(cl.bins, function(select.cl){
      tmp.cl= droplevels(cl[cl %in% select.cl])      
      dat = get_logNormal(big.dat, names(tmp.cl))
      result = sapply(stats, function(x){
          get_cl_stats(dat, cl=tmp.cl,stats=x)
      }, simplify=F)
    },simplify=F)
    cl.results = sapply(stats, function(x){
      do.call("cbind",sapply(tmp.results, function(result) result[[x]], simplify=F))
    },simplify=F)
    return(cl.results)    
  }

merge_cl_big <- function(big.dat,
                    cl, 
                    rd.dat=NULL,
                    rd.dat.t = NULL,
                    de.param = de_param(), 
                    merge.type = c("undirectional","directional"), 
                    max.cl.size = 300,
                    de.genes = NULL, 
                    return.markers = FALSE,
                    verbose = 0,
                    k=4,
                    mc.cores=1)
  {
    de.method = "fast_limma"
    if(!is.integer(cl)){
      cl = setNames(as.integer(as.character(cl)), names(cl))
    }
    merge.type=merge.type[1]
    de.df=list()
    pairs=NULL
    if(!is.null(de.genes)){
      pairs=do.call("rbind",strsplit(names(de.genes), "_"))
      row.names(pairs)=names(de.genes)
    }
     ###Merge small clusters with the closest neighbors first.
    if(!is.null(rd.dat)){
      cl.rd = as.data.frame(get_cl_means(rd.dat,cl[names(cl) %in% row.names(rd.dat)]))
    }
    else{
      cl.rd = as.data.frame(get_cl_means(rd.dat.t,cl[names(cl) %in% colnames(rd.dat.t)]))
    }
    cl.size = table(cl)
    while(TRUE){
      if(length(cl.size)==1){
        break
      }
      cl.small =  names(cl.size)[cl.size < de.param$min.cells]
      ###if all clusters are small, not need for further split. 
      if(length(cl.small)==length(cl.size)){
        return(NULL)
      }
      if(length(cl.small)==0){
        break
      }      
      merge.pairs = get_knn_pairs(cl.rd[,!colnames(cl.rd) %in% cl.small, drop=F], cl.rd[,cl.small,drop=F], k=1)
      x = merge.pairs[1,1]
      y=  merge.pairs[1,2]
      if(verbose > 0){
        cat("Merge: ", x,y, "sim:", merge.pairs[1,"sim"],"\n")
      }
      cl[cl==y]= x
      p = as.character(c(x,y))
      cl.rd[[p[1]]] = get_weighted_means(as.matrix(cl.rd[,p]), as.vector(cl.size[p]))
      cl.rd[[p[2]]] = NULL
      cl.size[p[1]] = sum(cl.size[p])
      cl.size = cl.size[names(cl.size)!=p[2]]
    }
    tmp.cl = cl[names(cl) %in% big.dat$col_id]
    tmp = get_cl_stats_big(big.dat, cl, max.cl.size=max.cl.size, stats=c("means","present","sqr_means"),mc.cores=mc.cores)
    cl.means = as.data.frame(tmp$means)
    cl.present = as.data.frame(tmp$present)    
    cl.sqr.means = as.data.frame(tmp$sqr_means)
    
    while(length(unique(cl)) > 1){
      merge.pairs = get_knn_pairs(cl.rd, cl.rd, k=k)
      ###Determine the de score for these pairs
      if(nrow(merge.pairs)==0){
        break
      }
      
      #####get DE genes for new pairs
      new.pairs = setdiff(row.names(merge.pairs),names(de.genes))
      if(verbose > 0){
        cat("Compute DE genes\n")
      }
      tmp.pairs = merge.pairs[new.pairs,,drop=FALSE] 
      de.result = de_selected_pairs(norm.dat=NULL, cl=cl, pairs=tmp.pairs, de.param= de.param, method=de.method, cl.means=cl.means, cl.present=cl.present, cl.sqr.means=cl.sqr.means,mc.cores=mc.cores, blocksize=1000)
      tmp.de.genes = de.result$de.genes
      de.genes[names(tmp.de.genes)] = tmp.de.genes
      pairs = get_pairs(names(de.genes))
      
      tmp.pairs= intersect(names(de.genes), row.names(merge.pairs))
      sc = sapply(de.genes[tmp.pairs], function(x){
        if(length(x)>0){x$score}
        else{0}
      })
      sc = sort(sc)
                                        #print(head(sc,10))      
      to.merge = sapply(names(sc), function(p){
        to.merge = test_merge(de.genes[[p]], de.param, merge.type=merge.type)
      })
      if(sum(to.merge)==0){
        break
      }
      sc = sc[to.merge]
      to.merge= merge.pairs[names(sc),,drop=FALSE]
      to.merge$sc = sc          
      
      merged =c()
      ###The first pair in to.merge always merge. For the remaining pairs, if both clusters have already enough cells,
      ###or independent of previus merging, then they can be directly merged as well, without re-assessing DE genes. 
      for(i in 1:nrow(to.merge)){
        p = c(to.merge[i,1], to.merge[i,2])
        if(i == 1 | sc[i] < de.param$de.score.th /2  & length(intersect(p, merged))==0){
          cl[cl==p[2]] = p[1]
          
          p = as.character(p)
          if(verbose > 0){
            cat("Merge ",p[1], p[2], to.merge[i,"sc"], to.merge[i, "sim"], cl.size[p[1]],"cells", cl.size[p[2]],"cells", "\n")
          }
          
          cl.rd[[p[1]]] = get_weighted_means(as.matrix(cl.rd[,p]), cl.size[p])
          cl.rd[[p[2]]] = NULL          
          cl.means[[p[1]]] = get_weighted_means(as.matrix(cl.means[, p]), cl.size[p])
          cl.means[[p[2]]] = NULL
          cl.present[[p[1]]] = get_weighted_means(as.matrix(cl.present[, p]), cl.size[p])
          cl.present[[p[2]]] = NULL
          cl.sqr.means[[p[1]]] = get_weighted_means(as.matrix(cl.sqr.means[, p]), cl.size[p])
          cl.sqr.means[[p[2]]] = NULL
          cl.size[p[1]] = sum(cl.size[p])
          cl.size = cl.size[names(cl.size)!=p[2]]
          rm.pairs = row.names(pairs)[pairs[,1]%in% p | pairs[,2]%in% p]
          de.genes = de.genes[setdiff(names(de.genes),rm.pairs)]
          merged = c(merged,p)
        }
      }
    }
    if(length(unique(cl))<2){
      return(NULL)
    }
    if(verbose > 0){
      print(table(cl))
    }
    markers = NULL
    if(return.markers){
      if(!is.null(max.cl.size)){
        sampled.cells = sample_cells(cl[names(cl) %in% big.dat$col_id],  max.cl.size)
        tmp.cl= cl[sampled.cells]
      }
      else{
        tmp.cl= cl
      }
      de.genes = de_all_pairs(norm.dat=NULL, cl=tmp.cl, de.param=de.param, cl.means=cl.means, cl.present=cl.present, cl.sqr.means=cl.sqr.means, mc.cores=mc.cores)    
      markers = select_markers(norm.dat=NULL, cl, de.genes=de.genes, n.markers=50, mc.cores=mc.cores)$markers
    }
    sc = sapply(de.genes, function(x){
      if(length(x)>0){x$score}
      else{0}
    })
    return(list(cl=cl, de.genes=de.genes,sc=sc, markers=markers))
  }




  
