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

impute_knn_big <- function(comb.dat, knn.result.list, select.genes, select.cells, ref.sets=names(comb.dat$dat.list), select.sets=names(comb.dat$dat.list), max.dim=100, k=15, th=0.5, rm.eigen=NULL,rm.th=0.65,method="zscore",mc.cores=1,verbose=FALSE)
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

