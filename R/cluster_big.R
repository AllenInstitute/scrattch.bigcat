library(bigstatsr)
convert_big.dat <- function(mat, logNormal=TRUE, backingfile=file.path(getwd(), "fbm"),...)
  {
    library(bigstatsr)
    m = FBM(nrow=nrow(mat),ncol=ncol(mat),backingfile=backingfile, ...)
    ind_nozero <- which(mat != 0, arr.ind = TRUE)
    m[ind_nozero] <- mat[ind_nozero]
    big.dat = list(fbm=m, row_id = row.names(mat), col_id = colnames(mat))
    big.dat$logNormal = logNormal
    return(big.dat)
  }

append_big.dat <- function(big.dat, mat)
  {
    offset = big.dat$fbm$ncol
    ind_nozero <- which(mat != 0, arr.ind = TRUE)
    tmp = ind_nozero
    tmp[,2] = tmp[,2] + offset
    big.dat$fbm$add_columns(ncol(mat))
    big.dat$fbm[tmp] <- mat[ind_nozero]
    big.dat$col_id = c(big.dat$col_id, colnames(mat))
    return(big.dat)
  }


convert_mat_list_big.dat <- function(mat.list, backingfile=NULL, ...)
  {
    if(is.null(backingfile)){
      backingfile=file.path(getwd(), paste0(Sys.Date(),"_fbm"))
    }
    big.dat =  convert_big.dat(mat.list[[1]],backingfile=backingfile, ...)
    for(i in 2:length(mat.list)){
      big.dat = append_big.dat(big.dat, mat.list[[i]])
    }
    return(big.dat)
  }

reorder_big.dat <- function(big.dat, new.cols,backingfile=NULL)
  {
    if(is.character(new.cols)){
      cols = match(new.cols, big.dat$col_id)
    }
    else{
      cols = new.cols
    }
    if(is.null(backingfile)){
      backingfile=file.path(getwd(), paste0(Sys.Date(),"_fbm"))
    }
    big.dat$fbm = big_copy(big.dat$fbm, ind.col=cols, backingfile=backingfile)
    big.dat$col_id = big.dat$col_id[cols]
    return(big.dat)
  }

convert_mat_list_big.dat_old<- function(mat.list, ...)
  {
    cn = unlist(lapply(mat.list, colnames))
    ncols = length(cn)
    nrows = nrow(mat.list[[1]])
    m = FBM(nrow=nrows,ncol=ncols, ...)
    big.dat = list(fbm=m, row_id = row.names(mat.list[[1]]) , col_id = cn)
    add.start = 1
    for(i in 1:length(mat.list)){
      mat = mat.list[[i]]
      add.end = add.start + ncol(mat) - 1
      big.dat$fbm[,add.start:add.end] = mat
      add.start = add.end + 1
    }
    return(big.dat)
  }

big_dat_apply <- function(big.dat, cols, p.FUN, .combine="c",  mcores=1, block.size = 10000,...)
  {
    require(foreach)
    require(doParallel)
    if (mcores == 1) {
      registerDoSEQ()
    }
    else {
      cl <- makeForkCluster(mcores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl), add = TRUE)
    }
    if(is.character(cols)){
      cols = match(cols, big.dat$col_id)
    }
    cols = sort(cols)    
    bins = split(cols, ceiling(1:length(cols)/block.size))
    res = foreach(bin = bins, .combine=.combine) %dopar% p.FUN(big.dat, bin, ...)
  }


get_cols <- function(big.dat, cols, keep.col=TRUE, sparse=TRUE)
  {
    library(Matrix)
    if(is.character(cols)){
      id = match(cols, big.dat$col_id)
    }
    else{
      id = cols
    }
    ord = order(id)

    mat = big.dat$fbm[,id[ord],drop=F]
    if(keep.col){
      org.order = (1:length(id))[order(ord)]    
      mat = mat[,org.order,drop=F]
      colnames(mat) = big.dat$col_id[id]      
    }
    else{
      colnames(mat) = big.dat$col_id[id[ord]]      
    }
    if(sparse){
      mat = Matrix(mat,sparse=sparse)
    }
    rownames(mat) = big.dat$row_id
    return(mat)
  }

get_logNormal <- function(big.dat, cols, select.genes=NULL, ...)
  {
    mat = get_cols(big.dat, cols,  ...)
    if(big.dat$logNormal){
      norm.dat = mat
    } 
    else{
      norm.dat=logCPM(mat)
    }
    if(!is.null(select.genes)){
      mat = mat[select.genes,,drop=F]
    }
    return(mat)
  }


get_counts <- function(big.dat, cols, ...)
  {
    mat = get_cols(big.dat, cols,...)
    if(big.dat$logNormal){
      mat@x = 2^mat@x -1
    }
    mat
  }



rd_PCA_big <- function(big.dat, dat, select.cells, max.dim=10, th=2, verbose=TRUE, mcores=1,method="zscore")
{
  system.time({tmp = get_PCA(dat, max.pca=max.dim, verbose=verbose,method=method,th=th)})
  if(is.null(tmp)){
    return(NULL)
  }
  rot = tmp$rot
  rd.dat = tmp$rd.dat
  pca = tmp$pca
  if(ncol(dat)< length(select.cells)){
    if(verbose){
      print("project")
    }
    system.time({rd.dat = big_project(big.dat, select.cells, rot, mcores=mcores)})
  }
  rm(rot)
  rm(dat)
  return(list(rd.dat=rd.dat, pca=pca))
}

big_project <- function(big.dat, select.cells, rot, mcores=1,...)
  {
    require(Matrix)
    my_project <- function(big.dat, cols, rot){
      dat = get_logNormal(big.dat, cols,...)[row.names(rot),]
      dat = Matrix::crossprod(dat, rot)      
      row.names(dat) = row.names(dat)
      gc()
      dat
    }
    rd.dat= big_dat_apply(big.dat, select.cells, p.FUN=my_project, .combine="rbind",  mcores=mcores, rot=rot)
    rd.dat = as.matrix(rd.dat)
    return(rd.dat)
  }


big_project2 <- function(big.dat, select.cells, rot, mcores=1)
  {
    cols=match(select.cells, big.dat$col_id)
    cols = sort(cols)
    rows = match(row.names(rot), big.dat$row_id)
    rd.dat=big_cprodMat(big.dat$fbm, rot, ind.row=rows, ind.col=cols, ncores=mcores)
    gc()
    rd.dat = as.matrix(rd.dat)
    colnames(rd.dat) = colnames(rot)
    row.names(rd.dat) = big.dat$col_id[cols]
    return(rd.dat[select.cells,])
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
                             counts = NULL,
                             method = c("louvain","ward.D", "kmeans"), 
                             vg.padj.th = 0.5, 
                             dim.method = c("pca","WGCNA"), 
                             max.dim = 20, 
                             rm.eigen = NULL, 
                             rm.th = 0.7, 
                             de.param = de_param(),
                             merge.type = c("undirectional", "directional"), 
                             maxGenes = 3000,
                             sampleSize = 20000,
                             max.cl.size = 300,
                             k.nn=15,
                             prefix = NULL, 
                             verbose = FALSE, overwrite=FALSE)
                            
  {
    library(matrixStats)
    outfile=paste0(prefix, ".rda")
    if(file.exists(outfile) & !overwrite){
      load(outfile)
      return(result)
    }
    
    method=method[1]
    merge.type=merge.type[1]

    
    sampled = sample(select.cells, min(sampleSize, length(select.cells)))
    if(length(sampled) > length(select.cells)/2){
      sampled = select.cells
    }
    if(is.null(counts) & sum(!sampled %in% colnames(counts))>0){
      system.time({counts = get_counts(big.dat, sampled,sparse=TRUE)})
    }
    else{
      counts = counts[,sampled]
    }
    if(verbose){
      print("Find high variance genes")
    }    
    system.time({vg = findVG(counts)})
    select.genes = with(vg, row.names(vg)[which(loess.padj < vg.padj.th | dispersion >3)])
    if(length(select.genes) < de.param$min.genes){
      return(NULL)
    }
    select.genes = head(select.genes[order(vg[select.genes, "loess.padj"],-vg[select.genes, "z"])],maxGenes)
    if(verbose){
      print("Reduce dimensions")
    }
    counts = logCPM(counts)[select.genes,]
    rd.result = rd_PCA_big(big.dat, dat=counts, select.cells, max.dim=max.dim, verbose=verbose)
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
    if(method=="louvain"){
      k = pmin(k.nn, round(nrow(rd.dat)/2))
      tmp = jaccard_louvain(rd.dat, k)
      if(is.null(tmp)){
        return(NULL)
      }
      cl = tmp$cl
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
    sampled.cells = sample_cells(cl, max.cl.size)
    norm.dat = get_logNormal(big.dat, sampled.cells)
    merge.result=merge_cl(norm.dat, cl=cl, rd.dat=rd.dat, merge.type=merge.type, de.param=de.param, max.cl.size=max.cl.size, verbose=verbose)
    rm(norm.dat)
    gc()
    if(is.null(merge.result)){
      return(NULL)
    }
    cl = merge.result$cl
    de.genes = merge.result$de.genes
    markers= merge.result$markers
    result=list(cl=cl, markers=markers)

    if(verbose){
      save(result, file=outfile)
    }
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
    counts = NULL
    if(is.null(result)){
      #if(is.null(counts)){
      #  sampled = sample(select.cells, min(sampleSize, length(select.cells)))
      #  if(length(sampled) > length(select.cells)/2){
      #    sampled = select.cells
      #  }
      #  cat("Allocate count matrix",length(sampled),"\n")
      #  counts = get_counts(big.dat, sampled,sparse=TRUE)
      #}          
      result=onestep_clust_big(big.dat=big.dat, select.cells=select.cells, prefix=prefix,method=select.method, counts=counts, sampleSize= sampleSize,...)      
      #if(sum(!select.cells %in% colnames(counts))>0){
      #  print("free count matrix")
      #  rm(counts)
      #  gc()
      #  counts=NULL
      #}
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
          tmp.result=iter_clust_big(big.dat=big.dat, select.cells=tmp.cells, prefix=tmp.prefix,split.size=split.size,method= method, counts=counts, sampleSize=sampleSize, ...)
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

