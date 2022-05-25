recompute_knn <- function(comb.dat, cl, ref.sets, select.sets, cl.group.list, cl.group.markers,sample.size=80000, knn.method="Annoy.Cosine")
  {
    knn.result.list = list()    
    for(l in names(cl.group.list)){
      cl.group = cl.group.list[[l]]
      if(length(cl.group)<=1){
        next
      }
      markers = cl.group.markers[[l]]      
      select.cl = cl[cl %in% cl.group]      
      group.cl = setNames(cl.group[as.character(select.cl)], names(select.cl))
      tmp.ref.list = sapply(ref.sets, function(x){
        tmp.cl= droplevels(select.cl[names(select.cl) %in% comb.dat$dat.list$col_id])
        sample_cells(select.cl, 100)
      })
      tmp.ref.list = sapply(tmp.ref.list, function(x){
        sample(x, min(length(x),sample.size))
      })
      select.cells=names(select.cl)
      knn.comb = compute_knn(comb.dat, select.genes=markers, ref.list=tmp.ref.list, select.sets=select.sets, select.cells=select.cells,cross.knn.method=cross.knn.method)
      knn.result.list[[l]] = list(select.genes=markers, knn=knn.comb)      
    }
    return(knn.result.list)
  }

impute_knn_big <- function(comb.dat, knn.result.list, select.genes, select.cells, ref.list, sets=names(comb.dat$dat.list), org.rd.dat.list=NULL, max.dim=100, k=15, th=0.5, rm.eigen=NULL,rm.th=0.65,method="zscore",mc.cores=1,verbose=FALSE)
  {
    library(matrixStats)
    knn.list <- list()
    impute.dat.list <- list()    
    if(is.null(org.rd.dat.list)){
      org.rd.dat.list <- list()
###Impute the reference dataset in the original space globally
      for(x in names(ref.list))
        {
          print(x)
          tmp.cells= select.cells[comb.dat$meta.df[select.cells,"platform"]==x]
          ref.cells = ref.list[[x]]
          if(comb.dat$type=="mem"){
            rd.result <- rd_PCA(comb.dat$dat.list[[x]], select.genes, select.cells=tmp.cells, sampled.cells = ref.cells, max.pca =max.dim, th=th, method=method,mc.cores=mc.cores,verbose=verbose)
          }
          else{
            rd.result <- rd_PCA_big(comb.dat$dat.list[[x]], select.genes, select.cells=tmp.cells, sampled.cells = ref.cells, max.pca =max.dim, th=th, method=method,mc.cores=mc.cores,verbose=verbose)
          }
          
          rd.result = org.rd.dat.list[[x]]
          if(!is.null(rm.eigen)){
            rd.result4rd.dat  = filter_RD(rd.result$rd.dat, rm.eigen, rm.th,verbose=verbose)
          }
          org.rd.dat.list[[x]] = rd.result
        }
    }
    for(x in names(ref.list))
      {
        ref.cells = ref.list[[x]]
        impute.dat.list[[x]] = init_mat_parqeut_dense(paste0("impute_dat_",x), col.id=select.cells, row.id=select.genes, val=0)        
        #print(ncol(rd.dat))        
        knn = get_knn_batch(rd.dat, rd.dat[ref.cells,], method="Annoy.Euclidean", mc.cores=mc.cores, batch.size=50000,k=k,transposed=FALSE)
        
        dat = as.matrix(get_cols(comb.dat$dat.list[[x]], select.cells=ref.cells, rows=select.gene))
        reference.id = 1:length(ref.cells)       
        gene.id = 1:length(select.genes)
        col.df = impute.dat.list[[x]]$col.df 
        col.df = col.df %>% filter(col_names %in% select.cells)
        foreach(c.id %in% unique(col.df$col_bin)) %dopar% {
          tmp.col.df = col.df %>% filter(col_bin==c.id)
          tmp.cells = tmp.col.df$col_names
          impute.dat = impute_knn(knn[tmp.cells,], reference.id, dat)
          update_mat_parquet_dense(impute.dat.list[[x]], impute.dat,mc.cores=1)
        }        
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
        big.dat = impute.dat.list[[ref.set]]
        dat = as.matrix(get_cols(big.dat, rows=impute.genes, cols=result$ref.list[[ref.set]]))
        gene.id = match(impute.genes, row.names(dat))
        cell.id = match(query.cells, colnames(dat))
        reference.id = match(comb.dat$all.cells, colnames(dat))        
        impute.dat = impute_knn(knn[tmp.cells,], reference.id, dat)
        update_mat_parquet_dense(impute.dat.list[[x]], impute.dat,mc.cores=1)
      }
    }
    return(list(knn.list =knn.list, org.rd.dat.list = org.rd.dat.list,impute.dat.list=impute.dat.list, ref.list=ref.list))
  }

