library(Matrix)
#library(ggplot2)
library(matrixStats)

jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
blue.red <-colorRampPalette(c("blue", "white", "red"))



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
      get_logNormal(comb.dat$dat.list[[set]], select.cells, select.genes=comb.dat$common.genes)
    },simplify=F)
    dat.list = dat.list[!sapply(dat.list,is.null)]
    return(dat.list)
  }


get_cells_logNormal <- function(comb.dat, cells)
  {
    dat = matrix(0, nrow=length(comb.dat$common.genes),ncol=length(cells), dimnames=list(comb.dat$common.genes, cells))
    for(set in names(comb.dat$dat.list)){
      select.cells = intersect(cells, comb.dat$dat.list[[set]]$col_id)
      if(length(select.cells)>0){
         tmp.dat = get_logNormal(comb.dat$dat.list[[set]], select.cells, select.genes=comb.dat$common.genes)
         dat[,colnames(tmp.dat)] = as.matrix(tmp.dat)
       }
    }
    return(dat)
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
    for(x in 2:length(dat.list)){
      common.genes= intersect(common.genes, dat.list[[x]]$row_id)
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
    comb.dat = list(dat.list=dat.list, meta.df = meta.df, cl.list=cl.list, cl.df.list = cl.df.list, de.genes.list = de.genes.list, de.param.list= de.param.list, common.genes=common.genes, all.cells= all.cells)
  }




select_joint_genes_big <-  function(comb.dat, ref.dat.list, select.cells = comb.dat$all.cells, maxGenes=2000, vg.padj.th=0.5, max.dim=20,use.markers=TRUE, top.n=100,rm.eigen=NULL, conservation.th = 0.5,rm.th=rep(0.7,ncol(rm.eigen)))
  {
    require(matrixStats)
    select.genes = lapply(names(ref.dat.list), function(ref.set){
      ref.dat = ref.dat.list[[ref.set]]
      ref.cells=colnames(ref.dat)
      cat(ref.set, length(ref.cells),"\n")
      tmp.cells=  intersect(select.cells, ref.cells)
###if cluster membership is available, use cluster DE genes
      if(use.markers & !is.null(comb.dat$de.genes.list[[ref.set]])){
        cl = droplevels(comb.dat$cl.list[[ref.set]][tmp.cells])
        cl.size = table(cl)
        cl = droplevels(cl[cl %in% names(cl.size)[cl.size > comb.dat$de.param.list[[ref.set]]$min.cells]])
        if(length(levels(cl)) <= 1){
          return(NULL)
        }
        de.genes = comb.dat$de.genes.list[[ref.set]]
        print(length(de.genes))     
        select.genes = display_cl(cl, norm.dat=ref.dat, max.cl.size = 200, n.markers=20, de.genes= de.genes)$markers
        select.genes = intersect(select.genes, comb.dat$common.genes)
      }
####if cluster membership is not available, use high variance genes and genes with top PCA loading
      else{
        tmp.dat = ref.dat
        tmp.dat@x = 2^tmp.dat@x - 1
        vg = find_vg(tmp.dat)
        rm(tmp.dat)
        gc()
        select.genes = row.names(vg)[which(vg$loess.padj < vg.padj.th | vg$dispersion >3)]
        if(length(select.genes) < 5){
          return(NULL)
        }
        select.genes = head(select.genes[order(vg[select.genes, "padj"],-vg[select.genes, "z"])],maxGenes)
        rd = rd_PCA(norm.dat=ref.dat,select.genes, ref.cells, max.pca = max.dim)
        if(is.null(rd)){
          return(NULL)
        }
        rd.dat = rd$rd.dat
        rot = t(rd$pca$rotation[,1:ncol(rd$rd.dat)])
        if(!is.null(rm.eigen)){
          rm.cor=cor(rd.dat, rm.eigen[row.names(rd.dat),])
          rm.cor[is.na(rm.cor)]=0
          select = colSums(t(abs(rm.cor)) >= rm.th) ==0
          print("Select PCA")
          print(table(select))
          if(sum(select)==0){
            return(NULL)
          }
          rot = rot[select,,drop=FALSE]
        }
        if(is.null(rot)){
          return(NULL)
        }
        rot.scaled = (rot  - rowMeans(rot))/rowSds(rot)
        gene.rank = t(apply(-abs(rot.scaled), 1, rank))
        select = gene.rank <= top.n & abs(rot.scaled ) > 2
        select.genes = colnames(select)[colSums(select)>0]
      }
    })
    gene.score = table(unlist(select.genes))
    if(length(gene.score)==0){
      return(NULL)
    }
    select.genes= names(head(sort(gene.score, decreasing=T), maxGenes))
    #gg.cons = gene_gene_cor_conservation(ref.dat.list, select.genes, select.cells)
    #select.genes = row.names(gg.cons)[gg.cons > conservation.th]
    return(select.genes)
  }


get_knn_batch_big <- function(big.dat, ref.dat, select.cells, k, method="cor", dim=NULL, batch.size, mc.cores=1,...)
  {
    results <- batch_process(x=select.cells, batch.size=batch.size, mc.cores=mc.cores, .combine="rbind", FUN=function(bin){
      dat = get_logNormal(big.dat, bin, keep.col=FALSE, sparse=FALSE)[row.names(ref.dat),,drop=F]
      knn=get_knn(dat=dat, ref.dat=ref.dat, k=k, method=method, dim=dim,...)
      rm(dat)
      gc()
      knn
    })
    return(results)
  }



compute_knn_big <- function(comb.dat, select.genes, ref.list, ref.dat.list=NULL, select.sets=names(comb.dat$dat.list), select.cells=comb.dat$all.cells, k=15, cross.knn.method=c("Annoy.Cosine","cor"), self.knn.method=c("Annoy.Euclidean","RANN"), batch.size=10000, mc.cores=1, rm.eigen=NULL, rm.th=0.7, max.dim=20)
  {    
    cat("Number of select genes", length(select.genes), "\n")
    cat("Get knn\n")
    dat.list = comb.dat$dat.list
###index is the index of knn from all the cells
    knn.list = list()
    for(ref.set in names(ref.list)){
      cat("Ref ", ref.set, "\n")
      if(length(ref.list[[ref.set]]) <= k) {
        ##Not enough reference points to compute k
        next
      }
      k.tmp = k
      if(length(ref.list[[ref.set]]) <= k*2) {
        k.tmp = round(k/2)
      }
      big.dat = dat.list[[ref.set]]
      map.cells=  big.dat$col_id[big.dat$col_id %in% select.cells]
      if(length(map.cells)==0){
        next
      }
      ref.cells = ref.list[[ref.set]]
      if(is.null(ref.dat.list)){
        ref.dat = get_logNormal(big.dat, ref.cells, sparse=FALSE, keep.col=FALSE)[select.genes, ,drop=F]}
      else{
        ref.dat = ref.dat.list[[ref.set]][select.genes,]
      }      
      tmp.cores = mc.cores
      if(length(map.cells)< batch.size){
        tmp.cores = 1
      }
      if(length(map.cells) == length(ref.cells)){
        rd.dat = rd_PCA(ref.dat, select.cells=map.cells, th=1, max.pca=max.dim)$rd.dat
      }
      else{
        rd.dat = rd_PCA_big(big.dat=dat.list[[ref.set]],dat = ref.dat, select.cells=map.cells, max.dim = max.dim, th=1, mc.cores=tmp.cores)$rd.dat
      }
      if(!is.null(rm.eigen)){
        rd.dat = filter_RD(rd.dat, rm.eigen, rm.th=rm.th)
      }
      ref.rd.dat = rd.dat[ref.cells,,drop=F]
      if(is.null(ref.rd.dat)){
        next
      }
      idx = match(ref.cells, comb.dat$all.cells)
      index = NULL
      if(length(select.cells) >50000 & self.knn.method %in% c("Annoy.Euclidean")){
        require(BiocNeighbors)
        index = buildAnnoy(ref.rd.dat, distance ="Euclidean", transposed = FALSE)
      }
      knn=get_knn_batch(dat=rd.dat, ref.dat = ref.rd.dat, k=k.tmp, method = self.knn.method, batch.size = batch.size, mc.cores=tmp.cores, index=index, transposed=FALSE)
      if(!is.null(index)){
        cleanAnnoyIndex(index)
      }      
      knn = matrix(idx[knn], nrow=nrow(knn), dimnames=list(row.names(knn), NULL))
      self.knn = knn
      index = NULL
      if(cross.knn.method  %in% c("Annoy.Euclidean", "Annoy.Cosine")){
        if(cross.knn.method=="Annoy.Cosine"){
          distance = "Cosine"          
        }
        else{
          distance = "Euclidean"
        }
        if(length(select.cells)>50000){
          index = buildAnnoy(ref.dat, distance =distance, transposed = TRUE)
        }
      }
      knn =do.call("rbind", lapply(setdiff(select.sets,ref.set), function(set){
        cat("Set ", set, "\n")
        map.cells=  dat.list[[set]]$col_id[dat.list[[set]]$col_id %in% select.cells]
        if(length(map.cells)==0){
          return(NULL)
        }
        tmp.cores = mc.cores
        if(length(map.cells)< batch.size){
          tmp.cores = 1
        }
        if(length(map.cells)==length(ref.list[[set]])){
          dat = ref.dat.list[[set]][select.genes,map.cells]          
          knn=get_knn_batch(dat=dat, ref.dat = ref.dat, k=k.tmp, method = cross.knn.method, batch.size = batch.size, mc.cores=tmp.cores, index=index, transposed=TRUE)
        }
        else{
          knn=get_knn_batch_big(big.dat=dat.list[[set]], ref.dat = ref.dat,select.cells=map.cells, k=k.tmp, method = cross.knn.method, batch.size = batch.size, mc.cores=tmp.cores, index=index, transposed=TRUE)
        }
        #if(!is.null(comb.dat$cl.list)){  
        #  test.knn = test_knn(knn, comb.dat$cl.list[[set]], colnames(ref.dat), comb.dat$cl.list[[ref.set]])          
        #  if(!is.null(test.knn)){
        #    cat("Knn", set, ref.set, method, "cl.score", test.knn$cl.score, "cell.score", test.knn$cell.score,"\n")
        #  }
        #}
        knn = matrix(idx[knn], nrow=nrow(knn), dimnames=list(row.names(knn), NULL))        
      }))
      if(!is.null(index)){
        cleanAnnoyIndex(index)
      }
      knn = rbind(self.knn, knn)
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
##' @param batch.size 
##' @param verbose 
##' @param mc.cores 
##' @param rm.eigen 
##' @param rm.th 
##' @param ... 
##' @return 
##' @author Zizhen Yao
knn_joint_big <- function(comb.dat, ref.sets=names(comb.dat$dat.list), select.sets= names(comb.dat$dat.list), merge.sets=ref.sets, select.cells=comb.dat$all.cells, select.genes=NULL,cross.knn.method="Annoy.Cosine", self.knn.method = "Annoy.Euclidean", method="leiden", k=15,  sample.size = 50000, cl.sample.size = 100, batch.size = 10000, verbose=TRUE,mc.cores=1,rm.eigen=NULL, rm.th=0.7,max.dim=20,...)
{
  if(length(select.cells) < batch.size){
    mc.cores=1
  }
  select.cells = comb.dat$all.cells[comb.dat$all.cells %in% select.cells]
  cat("Number of select cells", length(select.cells), "\n")
  cells.list = split(select.cells, comb.dat$meta.df[select.cells, "platform"])[select.sets]
  ref.list =  sample_sets_list(cells.list[ref.sets], comb.dat$cl.list[ref.sets], sample.size=sample.size, cl.sample.size = cl.sample.size)
  ref.list = sapply(ref.list, function(x){
    ref.cells=sample(x, min(sample.size, length(x)))
    ref.cells = comb.dat$all.cells[comb.dat$all.cells %in% ref.cells]
  },simplify=F)
  lapply(ref.list, function(x)print(length(x)))
  ref.sets = ref.sets[sapply(ref.list,length) >= sapply(comb.dat$de.param.list[ref.sets], function(x)x$min.cells)]
  if(length(ref.sets)==0){
    return(NULL)
  }
  ref.list = ref.list[ref.sets]
  ###Select genes for joint analysis
  cat("Get ref.dat.list\n")
  ref.dat.list = sapply(ref.sets, function(ref.set){
    get_logNormal(comb.dat$dat.list[[ref.set]], ref.list[[ref.set]], select.genes=comb.dat$common.genes)
  },simplify=F)
  if(is.null(select.genes)){
    select.genes = select_joint_genes_big(comb.dat, ref.dat.list = ref.dat.list,select.cells=select.cells, max.dim= max.dim,...)
  }
  if(length(select.genes) < 20){
    return(NULL)
  }
  ref.dat.list = sapply(ref.dat.list, function(ref.dat){
    tmp = Matrix::colSums(ref.dat[select.genes,]) == 0
    if(sum(tmp)>0){
      ref.dat[,!tmp,drop=F]
    }
    else{
      ref.dat
    }
  },simplify=F)
  ref.list = sapply(ref.dat.list, colnames, simplify=FALSE)
  cat("Get knn\n")
  knn.comb= compute_knn_big(comb.dat, select.genes=select.genes, ref.list=ref.list, ref.dat.list= ref.dat.list, select.sets=select.sets, select.cells=select.cells, k=k, cross.knn.method=cross.knn.method, self.knn.method=self.knn.method, batch.size=batch.size, mc.cores=mc.cores, rm.eigen=rm.eigen, rm.th=rm.th, max.dim=max.dim)
  if(is.null(knn.comb)){
    return(NULL)
  }
  #########
  #save(knn.comb, file="knn.comb.rda")
  sampled.cells = unlist(cells.list)
  if(length(sampled.cells)> sample.size){
    tmp.list =  sample_sets_list(cells.list, comb.dat$cl.list[select.sets], sample.size=sample.size, cl.sample.size = cl.sample.size)
    sampled.cells = unlist(tmp.list)
    if(length(sampled.cells)>100000){
      sampled.cells = sample(sampled.cells, 100000)
    }
    sampled.cells=union(sampled.cells, unlist(ref.list))
  }     
  result = knn_jaccard_clust(knn.comb[sampled.cells,],prune=1/ncol(knn.comb), method=method)
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
  print(cl.platform.counts)
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
  merge.dat.list = sapply(merge.sets, function(x){
    tmp.cells = with(comb.dat,intersect(row.names(meta.df)[meta.df$platform==x], names(cl)))
    if(length(tmp.cells)==0){
      return(NULL)
    }
    if(length(tmp.cells) == length(ref.list[[x]])){
      return(ref.dat.list[[x]])
    }
    sampled.cells = sample_cells(cl[tmp.cells],200)
    get_logNormal(comb.dat$dat.list[[x]], sampled.cells, select.genes=comb.dat$common.genes, keep.col=FALSE, sparse=TRUE)
  },simplify=F)
  cl= merge_cl_multiple(comb.dat=comb.dat, merge.dat.list=merge.dat.list, cl=cl, anchor.genes=select.genes)
  if(length(unique(cl))<=1){
    return(NULL)
  }
  print(table(cl))
  result$cl = cl
  result$select.genes= select.genes
  result$ref.de.param.list = ref.de.param.list
  return(result)
}



 

harmonize_big <- function(comb.dat, prefix, overwrite=TRUE, dir="./",...)
  {
    fn = file.path(dir, paste0(prefix, ".rda"))
    print(fn)
    if(!overwrite){
      if(file.exists(fn)){
        load(fn)
        return(result)
      }

    }
    result = knn_joint_big(comb.dat, ...)
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
i_harmonize_big<- function(comb.dat, select.cells=comb.dat$all.cells, ref.sets=names(comb.dat$dat.list), prefix="", result=NULL, overwrite=TRUE, sample.size = 50000, dir="./",...)
  {
    #attach(comb.dat)
    if(is.null(result)){
      if(comb.dat$type=="big"){
        print("big")
        result = harmonize_big(comb.dat=comb.dat, select.cells=select.cells, ref.sets=ref.sets, prefix=prefix, overwrite=overwrite,sample.size=sample.size, dir=dir,...)
      }
      else{
        print("mem")
         result = harmonize(comb.dat=comb.dat, select.cells=select.cells, ref.sets=ref.sets, prefix=prefix, overwrite=overwrite, sample.size=sample.size,dir=dir,...)
      }
    }
    if(is.null(result)){
      return(NULL)
    }
    all.results= list(result)
    names(all.results) = prefix
    cl = result$cl
    for(i in as.character(sort(unique(result$cl)))){
      tmp.prefix=paste(prefix, i,sep=".")
      print(tmp.prefix)
      select.cells= names(cl)[cl == i]
      platform.size = table(comb.dat$meta.df[select.cells, "platform"])      
      print(platform.size)
      pass.th = sapply(sets, function(set)platform.size[[set]] >= comb.dat$de.param.list[[set]]$min.cells)
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
          tmp.result = i_harmonize(comb.dat, select.cells=select.cells, ref.sets=ref.sets, prefix=tmp.prefix, overwrite=overwrite, sample.size=sample.size,dir=dir,result=result,...)
        }
        else{
          if(!computed & length(select.cells) < sample.size){
            new.comb.dat = comb.dat
            dat.list = sample_cl_dat(comb.dat, sets, cl=cl[select.cells],sample.size)
            new.comb.dat$dat.list = dat.list
            new.comb.dat$type="mem"
            print("mem")
            tmp.result = i_harmonize(comb.dat=new.comb.dat, select.cells=select.cells, ref.sets=ref.sets, prefix=tmp.prefix, overwrite=overwrite, sample.size=sample.size,dir=dir,result=result,...)
            rm(new.comb.dat)
            gc()
          }
          else{
            print("big")
            tmp.result = i_harmonize_big(comb.dat, select.cells=select.cells, ref.sets=ref.sets, prefix=tmp.prefix, overwrite=overwrite, sample.size=sample.size, dir=dir,result=result,...)
          }
        }
        if(!is.null(tmp.result)){
          all.results[names(tmp.result)] = tmp.result
        }
      }
    }
    return(all.results)
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
        




#### assume within data modality have been performed
####
impute_knn_global_big<- function(comb.dat, split.results, select.genes, select.cells, ref.dat.list, ref.sets=names(ref.dat.list), sets=comb.dat$sets, rm.eigen=NULL, rm.th=0.7, verbose=FALSE,mc.cores=5, org.rd.dat.list=NULL)
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
          org.rd.dat.list[[x]] = rd.result
        }
    }
    
    for(x in ref.sets){
        rd.dat  = org.rd.dat.list[[x]]$rd.dat
        ref.cells=colnames(ref.dat.list[[x]])
        if(!is.null(rm.eigen)){
          rd.dat = filter_RD(rd.dat, rm.eigen, rm.th, verbose=verbose)
        }
        if(verbose){
          print(ncol(rd.dat))
        }
        knn.result <- RANN::nn2(data=rd.dat[ref.cells,], query=rd.dat, k=15)
        knn <- knn.result[[1]]
        row.names(knn) = row.names(rd.dat)    
        knn.list[[x]]=knn
        impute.dat.list[[x]] <- impute_knn(knn, ref.cells, as.matrix(t(ref.dat.list[[x]][select.genes,ref.cells])))
      }
    ###cross-modality Imputation based on nearest neighbors in each iteraction of clustering using anchoring genes or genes shown to be differentiall expressed.

    for(x in names(split.results)){
      print(x)
      result = split.results[[x]]
      cl = result$cl
      knn = result$knn
      
      for(ref.set in ref.sets){
        if(ref.set %in% names(result$ref.list)){
          tmp.cells = row.names(result$knn)
          add.cells=FALSE
          query.cells = intersect(tmp.cells[comb.dat$meta.df[tmp.cells,"platform"] != ref.set], select.cells)
          if(any(!query.cells %in% row.names(impute.dat.list[[ref.set]]))){
            add.cells=TRUE
            impute.genes = select.genes
          }
          else{
            impute.genes=intersect(select.genes,c(result$select.markers, result$select.genes))
          }
          select.cols = comb.dat$meta.df[comb.dat$all.cells[result$knn[1,]],"platform"] == ref.set
          if(sum(select.cols)==0){
            next
          }
          else{
            ref.cells = intersect(comb.dat$all.cells[unique(as.vector(knn[, select.cols]))],select.cells)            
            select.knn = result$knn[query.cells,select.cols]
            impute.dat = impute_knn(select.knn, comb.dat$all.cells, impute.dat.list[[ref.set]][ref.cells,impute.genes])
          }
          if(!add.cells){
            impute.dat.list[[ref.set]][query.cells, impute.genes] <- impute.dat
          }
          else{
            impute.dat.list[[ref.set]] <- rbind(impute.dat.list[[ref.set]],impute.dat)
          }
          rm(impute.dat)
          gc()
        }
      }
    }
    return(list(knn.list =knn.list, org.rd.dat.list = org.rd.dat.list,impute.dat.list=impute.dat.list))
  }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param comb.dat 
##' @param all.results 
##' @param select.genes 
##' @param select.cells 
##' @param ref.sets 
##' @return 
##' @author Zizhen Yao
impute_knn_cross <- function(comb.dat, all.results, select.genes, select.cells, ref.sets=ref.sets)
  {
    #impute.dat.list <<- list()
    return(impute.dat.list)
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

