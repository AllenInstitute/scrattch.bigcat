impute_dat_big <- function(impute.dat.big, ref.dat, knn, ref.cells, select.genes, block.size=100000, g.block.size=1000)
  {
    used.ref.cells = ref.cells[unique(as.vector(knn))]
    if(!is(ref.dat,"Matrix") & !is.matrix(ref.dat)){
      ref.dat = get_cols(ref.dat, cols=used.ref.cells, rows=select.genes,sparse=FALSE)
      if(is.numeric(ref.dat)){
        ref.dat = matrix(ref.dat, nrow=length(select.genes),dimnames=list(select.genes, used.ref.cells))
      }
    }
    ##Sort select.genes and knn based on their order in impute.dat.big
    select.genes = impute.dat.big$row_id[impute.dat.big$row_id %in% select.genes]
    query.cells = row.names(knn)
    query.cells = impute.dat.big$col_id[impute.dat.big$col_id %in% query.cells]
    knn= knn[query.cells,]
    bins = ceiling(seq_len(nrow(knn))/block.size)    
    g.bins = split(select.genes, ceiling(1:length(select.genes)/g.block.size))
    for(j in 1:length(g.bins)){
      tmp.genes=g.bins[[j]]
      select.ref.dat = as.matrix(ref.dat[tmp.genes,,drop=FALSE])
      for(i in 1:max(bins)){
        cat(i,j,"\n")
        tmp.cells = row.names(knn)[bins==i]
        impute.dat = impute_knn(knn[tmp.cells, ,drop=FALSE], ref.cells, dat = select.ref.dat, transpose_input = FALSE)
        set_cols_fbm(impute.dat.big,impute.dat)
      }
    }    
  }

impute_knn_global_big <- function(comb.dat, split.results, select.genes, select.cells, ref.list, sets=names(comb.dat$dat.list), org.rd.dat.list=NULL, max.dim=100, k=15, th=0.5, rm.eigen=NULL,rm.th=0.65,method="zscore",mc.cores=1,verbose=FALSE,impute.dat.list=NULL)
  {
    library(matrixStats)
    knn.list <- list()
    if(is.null(impute.dat.list)){
      impute.dat.list <- list()
    }
    if(is.null(org.rd.dat.list)){
      org.rd.dat.list <- list()
###Impute the reference dataset in the original space globally
      for(x in names(ref.list))
        {
          print(x)
          tmp.cells= select.cells[comb.dat$meta.df[select.cells,"platform"]==x]
          ref.cells = ref.list[[x]]
          if(comb.dat$type=="mem"){
            rd.dat <- rd_PCA(comb.dat$dat.list[[x]], select.genes, select.cells=tmp.cells, sampled.cells = ref.cells, max.pca =max.dim, th=th, method=method,mc.cores=mc.cores,verbose=verbose)$rd.dat
          }
          else{
            rd.dat <- rd_PCA_big(comb.dat$dat.list[[x]], select.genes, select.cells=tmp.cells, sampled.cells = ref.cells, max.pca =max.dim, th=th, method=method,mc.cores=mc.cores,verbose=verbose)$rd.dat
          }          
          if(!is.null(rm.eigen)){
            rd.dat  = filter_RD(rd.dat, rm.eigen, rm.th,verbose=verbose)
          }
          org.rd.dat.list[[x]] = rd.result
        }
    }
    for(ref.set in names(ref.list))
      {
        rd.dat = org.rd.dat.list[[ref.set]]
        rd.dat = rd.dat[row.names(rd.dat) %in% select.cells,]
        ref.cells = ref.list[[ref.set]]
        if(is.null(impute.dat.list[[ref.set]])){          
          impute.dat.list[[ref.set]] = create_big.dat_fbm(col.id=select.cells, row.id=select.genes,backingfile=paste0("impute_data_",ref.set))
        }
        knn = get_knn_batch(rd.dat, rd.dat[ref.cells,], method="Annoy.Euclidean", mc.cores=mc.cores, batch.size=50000,k=k,transposed=FALSE,clear.index=TRUE)        
        impute.dat.big = impute.dat.list[[ref.set]]
        impute_dat_big(impute.dat.big, comb.dat$dat.list[[x]], knn, ref.cells, select.genes)
      }
    ###cross-modality Imputation based on nearest neighbors in each iteraction of clustering using anchoring genes or genes shown to be differentiall expressed.
    for(x in names(split.results)){      
      result = split.results[[x]]
      impute.genes = intersect(c(result$impute.genes,result$knn.genes), select.genes)
      cat("split group",x,length(impute.genes),"\n")
      cl = result$cl
      for(ref.set in names(ref.list)){
        impute.dat.big=impute.dat.big.list[[ref.set]]
        ref.cells=result$ref.list[[ref.set]]                      
        tmp.cells = row.names(knn)
        query.cells = intersect(tmp.cells[comb.dat$meta.df[tmp.cells,"platform"] != ref.set], select.cells)
        select.cols = comb.dat$meta.df[comb.dat$all.cells[knn[1,]],"platform"] == ref.set
        if(sum(select.cols)==0){
          next
        }
        if(length(query.cells)==0){
          next
        }
        knn = result$knn        
        select.knn = knn[query.cells,select.cols,drop=F]
        impute_dat_big(impute.dat.big, big.dat=impute.dat.big, knn=knn, ref.cells=ref.cells,select.genes=impute.genes)                           
      }
    }       
    return(impute.dat.list)
  }

impute_cross_knn_big <- function(split.results, ref.dat, query.dat, query.cells, impute.genes = split.results[[1]]$impute.genes, prefix=format(Sys.time(), '%Y_%m_%d.%H.%M'),k=15,method = "Annoy.Cosine", mc.cores=10, clear.index=TRUE, impute.dat.big=NULL)
  {
    if(is.null(impute.dat.big)){
      impute.dat.big = create_big.dat_fbm(col.id=query.cells, row.id=impute.genes,backingfile=paste0("impute_data_",prefix))
    }       
    for(g in names(split.results)){
      result = split.results[[g]]
      tmp.cl=result$cl
      select.impute.genes = intersect(result$impute.genes,impute.genes)
      knn.genes = result$knn.genes
      if(length(knn.genes)<5){
        next
      }
      cat("split group",g,length(select.impute.genes),"\n")
      select.query.cells= intersect(names(tmp.cl), query.cells)
      if(length(select.query.cells)==0){
        next
      }
      if(is(ref.dat,"Matrix")|is.matrix(ref.dat)){
        select.ref.cells = intersect(names(tmp.cl), colnames(ref.dat))
        select.ref.cells = sample_cells(tmp.cl[select.ref.cells], 100)
        select.ref.dat=ref.dat[knn.genes,select.ref.cells,drop=FALSE]
      }
      else{
        select.ref.cells = intersect(names(tmp.cl), ref.dat$col_id)
        select.ref.cells = sample_cells(tmp.cl[select.ref.cells], 100)        
        select.ref.dat = get_logNormal(ref.dat, select.ref.cells, knn.genes,sparse=FALSE)
      }
      if(is(query.dat,"Matrix")){
        knn = get_knn_batch(query.dat[,select.query.cells], ref.dat = select.ref.dat,  k=k, method = method, mc.cores=mc.cores, clear.index=clear.index)
      }
      else{
        knn=get_knn_batch_big(query.dat, ref.dat = select.ref.dat, select.cells=select.query.cells, k=k, method = method, mc.cores=mc.cores, clear.index=clear.index)
      }      
      split.results[[g]]$knn = knn
      print("impute")
      impute_dat_big(impute.dat.big, ref.dat=ref.dat, knn=knn, ref.cells=select.ref.cells,select.genes=select.impute.genes)
    }
    return(list(impute.dat.big=impute.dat.big,split.results=split.results))
  }

