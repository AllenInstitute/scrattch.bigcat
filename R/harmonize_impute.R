
#' Title
#'
#' @param knn.dist 
#' @param scale 
#' @param exclude.th 
#'
#' @return
#' @export
#'
#' @examples
get_knn_weight <- function(knn.dist, scale=0.2, exclude.th = 0.0001)
  {
    w = exp(-knn.dist*scale)
    if(exclude.th >= 0){
      w[knn.dist < exclude.th] = 0
    }
    return(w)
  }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param knn.idx 
##' @param reference 
##' @param cl 
##' @return 
##' @author Zizhen Yao
predict_knn <- function(knn.idx, reference, cl, mc.cores=1)
  {
    library(matrixStats)
    library(data.table)
    library(dplyr)
    library(parallel)
    defaultW <- getOption("warn")
    options(warn = -1)
    query = row.names(knn.idx)
    if(nrow(knn.idx) < 10000){
      mc.cores=1
    }
    results = parallel::pvec(query,function(x){
      df = data.table(nn=as.vector(knn.idx[x,,drop=F]), query=rep(x, ncol(knn.idx)))
      df = df %>% filter(!is.na(nn))
      df$nn.cl = cl[reference[df$nn]]
      tb = df %>% group_by(query, nn.cl) %>% summarise(n=n())%>% mutate(freq = n/sum(n))
      pred.df = tb %>% group_by(query) %>% summarize(pred.cl = nn.cl[which.max(freq)], pred.score = max(freq))
      list(list(pred.df = pred.df, tb=tb))
    },mc.cores=mc.cores)
    pred.df = do.call("rbind",lapply(results, function(x)x[[1]]))
    pred.prob =  do.call("rbind",lapply(results, function(x)x[[2]]))    
    pred.df = as.data.frame(pred.df)
    row.names(pred.df) = pred.df$query
    pred.df$query=NULL
    options(warn = defaultW)
    return(list(pred.df=pred.df, pred.prob = pred.prob))
  }



#' Title
#'
#' @param knn.idx 
#' @param reference 
#' @param dat 
#' @param knn.dist 
#' @param w 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
impute_knn <- function(knn.idx, reference, dat, knn.dist=NULL, w=NULL, transpose_input=FALSE, transpose_output= FALSE)
  {
    
    if(transpose_input){
      reference.id = match(reference, row.names(dat))
      genes = colnames(dat)
    }
    else{
      reference.id = match(reference, colnames(dat))
      genes=row.names(dat)
    }
    
    cell.id = 1:nrow(knn.idx)
    gene.id = 1:length(genes)
    if(transpose_output){        
      impute.dat = matrix(0,length(cell.id), length(genes))
      row.names(impute.dat) = row.names(knn.idx)
      colnames(impute.dat) = genes
    }
    else{
      impute.dat = matrix(0, length(genes),length(cell.id))
      row.names(impute.dat) = genes
      colnames(impute.dat) =  row.names(knn.idx)
    }
    if(!is.null(w) & !is.numeric(w) & !is.null(knn.dist)){
      w= knn.dist
      w = w - rowMaxs(w) - 0.001
      w = -w/rowMaxs(abs(w))
      w = w/rowSums(w)      
    }
    ImputeKnn(knn.idx, reference.id, cell.id, gene.id,
              dat=dat,impute.dat, w_mat_=w,
              transpose_input=transpose_input, transpose_output=transpose_output)    
    impute.dat
  }




impute_knn_global <- function(comb.dat, split.results, select.genes, select.cells, ref.list, sets=names(comb.dat$dat.list), max.dim=100, k=15, th=0.5, rm.eigen=NULL,rm.th=0.65,method="zscore",mc.cores=1,verbose=FALSE)
  {
    library(matrixStats)
    org.rd.dat.list <- list()
    knn.list <- list()
    impute.dat.list <- list()
    ###Impute the reference dataset in the original space globally
    for(x in names(ref.list))
      {
        print(x)
        tmp.cells= select.cells[comb.dat$meta.df[select.cells,"platform"]==x]
        ref.cells = ref.list[[x]]
        
        rd.result <- rd_PCA(comb.dat$dat.list[[x]], select.genes, select.cells=tmp.cells, sampled.cells = ref.cells, max.pca =max.dim, th=th, method=method,mc.cores=mc.cores,verbose=verbose)
        org.rd.dat.list[[x]] = rd.result
        rd.result = org.rd.dat.list[[x]]
        if(!is.null(rm.eigen)){
          rd.dat  = filter_RD(rd.result$rd.dat, rm.eigen, rm.th,verbose=verbose)
        }
        #print(ncol(rd.dat))        
        knn = get_knn_batch(rd.dat, rd.dat[ref.cells,], method="Annoy.Euclidean", mc.cores=mc.cores, batch.size=50000,k=k,transposed=FALSE)
        dat = as.matrix(comb.dat$dat.list[[x]][select.genes,ref.cells])        
        reference.id = 1:length(ref.cells)
        cell.id = match(row.names(rd.dat), select.cells)                
        gene.id = 1:length(select.genes)
        impute.dat = matrix(0, nrow=length(select.genes), ncol=length(select.cells))
        dimnames(impute.dat) = list(select.genes, select.cells)
        ImputeKnn(knn, reference.id, cell.id, gene.id, dat=dat, impute.dat, w_mat_ = NULL,
                  transpose_input=FALSE, transpose_output=FALSE);
        impute.dat.list[[x]] = impute.dat            
      }
    
    ###cross-modality Imputation based on nearest neighbors in each iteraction of clustering using anchoring genes or genes shown to be differentiall expressed. 
    for(x in names(split.results)){
      print(x)
      result = split.results[[x]]
      if(x == names(split.results)[1]){
        impute.genes = select.genes
      }
      else{
        impute.genes=intersect(select.genes,c(result$markers, result$select.genes))
      }
      cl = result$cl
      knn = result$knn
      for(ref.set in intersect(names(result$ref.list),names(ref.list))){
        tmp.cells = row.names(knn)
        query.cells = intersect(tmp.cells[comb.dat$meta.df[tmp.cells,"platform"] != ref.set], select.cells)
        select.cols = comb.dat$meta.df[comb.dat$all.cells[knn[1,]],"platform"] == ref.set
        if(sum(select.cols)==0){
          next
        }
        if(length(query.cells)==0){
          next
        }

        select.knn = knn[query.cells,select.cols,drop=F]
        dat = impute.dat.list[[ref.set]]
        gene.id = match(impute.genes, row.names(dat))
        cell.id = match(query.cells, colnames(dat))
        reference.id = match(comb.dat$all.cells, colnames(dat))
        ImputeKnn(select.knn, reference.id, cell.id, gene_idx=gene.id, 
                  dat=dat,impute_dat=dat, w_mat_= NULL,
                  transpose_input=FALSE, transpose_output=FALSE)        
      }
    }
    return(list(knn.list =knn.list, org.rd.dat.list = org.rd.dat.list,impute.dat.list=impute.dat.list, ref.list=ref.list))
  }


prepare_impute_genes <- function(ds, cl.bin, all.group, top.n=10)
  {
    all.g = sort(unique(cl.group))
    for(g in all.g){
      tmp.cl=names(cl.group)[cl.group==g]
      tmp.bin = cl.bin %>% filter(cl %in% tmp.cl) %>% pull(bin) %>% unique
      select.genes = ds %>% filter(bin.x %in% tmp.bin & bin.y %in% tmp.bin & P1 %in% tmp.cl & P2 %in% tmp.cl & rank <= top.n) %>% select(gene) %>% distinct() %>% collect() %>% pull(gene)      
      markers = ds %>% filter(bin.x %in% tmp.bin & bin.y %in% tmp.bin & P1 %in% tmp.cl & P2 %in% tmp.cl & rank <= 20) %>% select(gene) %>% distinct() %>% collect() %>% pull(gene)
    }    
    markers = intersect(markers, select.markers)
    split.results[[as.character(g)]] = list(cl=droplevels(cl[cl %in% tmp.cl]), select.genes=select.genes, markers=markers)
  }

get_impute_knn <- function(comb.dat, cl, cl.group, ref.set)
  {
    for(g in names(cl.group)){
      tmp.cl=cl[cl %nii% cl.group[[g]]]
      select.genes = split.results[[g]]$select.genes
      select.cells = names(tmp.cl)
      cells.list = split(select.cells, comb.dat$meta.df[select.cells, "platform"])
      ref.cells = sample_cells(droplevels(tmp.cl[cells.list[[ref.set]]]), 100)
      idx = match(ref.cells, comb.dat$all.cells)
      map.cells= unlist(cells.list[!names(cells.list)==ref.set])
      if(length(map.cells)==0){
        split.results[[g]]=NULL
        next
      }
  ref.dat = get_logNormal(comb.dat$dat.list[[ref.set]], ref.cells, select.genes)
  knn.list = list()
  for(x in map.set){
    tmp.cells= intersect(map.cells, comb.dat$dat.list[[x]]$col_id)
    knn=get_knn_batch_big(comb.dat$dat.list[[x]], ref.dat = ref.dat, select.cells=tmp.cells, k=15, method = "Annoy.Cosine", mc.cores=10, clear.index=TRUE)
    knn.list[[x]] = matrix(idx[knn], nrow=nrow(knn), dimnames=list(row.names(knn), NULL))
  }
  knn = do.call("rbind",knn.list)
  split.results[[g]]$knn = knn
  ref.list = list(ref.cells)
  names(ref.list)=ref.set
  split.results[[g]]$ref.list=ref.list

    }
  }


impute_knn_global<- function(comb.dat,joint_cl, select.genes, select.cells, ref.list, sets=names(comb.dat$dat.list), max.dim=100, k=15, th=0.5, rm.eigen=NULL,rm.th=0.65,method="zscore",mc.cores=1,verbose=FALSE)
  {
    library(matrixStats)
    org.rd.dat.list <- list()
    knn.list <- list()
    impute.dat.list <- list()
    ###Impute the reference dataset in the original space globally
    for(x in names(ref.list))
      {
        print(x)
        tmp.cells= select.cells[comb.dat$meta.df[select.cells,"platform"]==x]
        ref.cells = ref.list[[x]]
        
        rd.result <- rd_PCA(comb.dat$dat.list[[x]], select.genes, select.cells=tmp.cells, sampled.cells = ref.cells, max.pca =max.dim, th=th, method=method,mc.cores=mc.cores,verbose=verbose)
        org.rd.dat.list[[x]] = rd.result
        rd.result = org.rd.dat.list[[x]]
        if(!is.null(rm.eigen)){
          rd.dat  = filter_RD(rd.result$rd.dat, rm.eigen, rm.th,verbose=verbose)
        }
        #print(ncol(rd.dat))        
        knn = get_knn_batch(rd.dat, rd.dat[ref.cells,], method="Annoy.Euclidean", mc.cores=mc.cores, batch.size=50000,k=k,transposed=FALSE)
        dat = as.matrix(comb.dat$dat.list[[x]][select.genes,ref.cells])        
        reference.id = 1:length(ref.cells)
        cell.id = match(row.names(rd.dat), select.cells)                
        gene.id = 1:length(select.genes)
        impute.dat = matrix(0, nrow=length(select.genes), ncol=length(select.cells))
        dimnames(impute.dat) = list(select.genes, select.cells)
        ImputeKnn(knn, reference.id, cell.id, gene.id, dat=dat, impute.dat, w_mat_ = NULL,
                  transpose_input=FALSE, transpose_output=FALSE);
        impute.dat.list[[x]] = impute.dat            
      }
    
    ###cross-modality Imputation based on nearest neighbors in each iteraction of clustering using anchoring genes or genes shown to be differentiall expressed. 
    for(x in names(split.results)){
      print(x)
      result = split.results[[x]]
      if(x == names(split.results)[1]){
        impute.genes = select.genes
      }
      else{
        impute.genes=intersect(select.genes,c(result$markers, result$select.genes))
      }
      cl = result$cl
      knn = result$knn
      for(ref.set in intersect(names(result$ref.list),names(ref.list))){
        tmp.cells = row.names(knn)
        query.cells = intersect(tmp.cells[comb.dat$meta.df[tmp.cells,"platform"] != ref.set], select.cells)
        select.cols = comb.dat$meta.df[comb.dat$all.cells[knn[1,]],"platform"] == ref.set
        if(sum(select.cols)==0){
          next
        }
        if(length(query.cells)==0){
          next
        }

        select.knn = knn[query.cells,select.cols,drop=F]
        dat = impute.dat.list[[ref.set]]
        gene.id = match(impute.genes, row.names(dat))
        cell.id = match(query.cells, colnames(dat))
        reference.id = match(comb.dat$all.cells, colnames(dat))
        ImputeKnn(select.knn, reference.id, cell.id, gene_idx=gene.id, 
                  dat=dat,impute_dat=dat, w_mat_= NULL,
                  transpose_input=FALSE, transpose_output=FALSE)        
      }
    }
    return(list(knn.list =knn.list, org.rd.dat.list = org.rd.dat.list,impute.dat.list=impute.dat.list, ref.list=ref.list))
  }

