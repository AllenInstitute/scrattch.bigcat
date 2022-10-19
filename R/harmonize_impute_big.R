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

impute_dat_big <- function(impute.dat.big, big.dat, knn, ref.cells, select.genes, block.size=100000, g.block.size=1000)
  {
    used.ref.cells = ref.cells[unique(as.vector(knn))]
    dat = get_cols(big.dat, cols=used.ref.cells, rows=select.genes)
    bins = ceiling(seq_len(nrow(knn))/block.size)
    g.bins = split(select.genes, ceiling(1:length(select.genes)/g.block.size))
    for(j in 1:length(g.bins)){
      tmp.genes=g.bins[[j]]
      ref.dat = as.matrix(dat[tmp.genes,])
      for(i in 1:max(bins)){
        cat(i,j,"\n")
        tmp.cells = row.names(knn)[bins==i]
        impute.dat = impute_knn(knn[tmp.cells, ,drop=FALSE], ref.cells, dat = ref.dat, transpose_input = FALSE)
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
    for(x in names(ref.list))
      {
        rd.dat = org.rd.dat.list[[x]]
        rd.dat = rd.dat[row.names(rd.dat) %in% select.cells,]
        ref.cells = ref.list[[x]]
        if(is.null(impute.dat.list[[x]])){          
          impute.dat.list[[x]] = create_big.dat_fbm(col.id=select.cells, row.id=select.genes,backingfile=paste0("impute_data_",x))
        }
        knn = get_knn_batch(rd.dat, rd.dat[ref.cells,], method="Annoy.Euclidean", mc.cores=mc.cores, batch.size=50000,k=k,transposed=FALSE,clear.index=TRUE)        
        impute.dat.big = impute.dat.list[[x]]
        impute_dat_big(impute.dat.big, comb.dat$dat.list[[x]], knn, ref.cells, select.genes)
###cross-modality Imputation based on nearest neighbors in each iteraction of clustering using anchoring genes or genes shown to be differentiall expressed. 
        for(x in names(split.results)){      
          result = split.results[[x]]
          impute.genes = intersect(c(result$markers,result$select.genes), select.genes)
          cat("split group",x,length(impute.genes),"\n")
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
            impute_dat_big(impute.dat.big, big.dat=impute.dat.big, knn=knn, ref.cells=comb.dat$all.cells, select.genes=impute.genes)        
          }
        }
      }
    
    return(impute.dat.list)
  }

