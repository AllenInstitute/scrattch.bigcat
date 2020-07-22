iter_consensus_clust_big <- function(cl.list, co.ratio=NULL,  cl.mat=NULL, big.dat=NULL,  select.cells=names(cl.list[[1]]), diff.th=0.25, prefix=NULL, method=c("auto", "louvain","ward.D"), verbose=FALSE, de.param = de.param, max.cl.size = 200, result=NULL, split.size = de.param$min.cells*2,merge.type=c("undirectional", "directional"))
{
  method=method[1]
  require(igraph)
  if(verbose){
    cat(prefix, length(select.cells),"\n")
  }
  if(!is.null(result)){
    markers=result$markers
    cl = setNames(as.integer(as.character(result$cl)),names(result$cl))
    cell.cl.co.ratio= get_cell.cl.co.ratio(cl, co.ratio= co.ratio, cl.mat=cl.mat[,names(cl)])
  }
  else{
    markers=NULL
    if(length(select.cells)  < split.size){
      return(NULL)
    }
    co.ratio.sampled = FALSE
    if(is.null(co.ratio)){
      cl.size = table(cl.list[[1]][select.cells])
      graph.size= sum(cl.size^2)
      if(graph.size > 10^8){
        co.ratio.sampled=TRUE
        tmp.cl.list = lapply(cl.list, function(cl)cl[select.cells])
        sampled.cells = sample_cl_list(tmp.cl.list, max.cl.size=max.cl.size)
        
        cl.size = table(cl.list[[1]][sampled.cells])
        graph.size= sum(cl.size^2)
        if(graph.size > 10^8){
          sampled.cells = sample_cells(cl.list[[1]][sampled.cells], max.cl.size)
        }
      }
      else{
        sampled.cells=select.cells
      }
      co.ratio = Matrix::crossprod(cl.mat[,sampled.cells])
      co.ratio@x = co.ratio@x/length(cl.list)
    }
    if(method=="auto"){
      if (length(select.cells)> 3000){
        select.method = "louvain"
      }
      else{
        select.method="ward.D"
      }
    }
    else{
      select.method = method
    }
    if(select.method=="ward.D"){
      if(!is.matrix(co.ratio)){
        tmp.co.ratio = as.matrix(co.ratio[select.cells, select.cells])
      }
      else{
        tmp.co.ratio = co.ratio
      }
      tmp.cl = init_cut(tmp.co.ratio, select.cells, cl.list, min.cells= de.param$min.cells, th = diff.th,method=select.method)
      rm(tmp.co.ratio)
      if(is.null(tmp.cl)){
        return(NULL)
      }  
    }
    else{###louvain
      if(co.ratio.sampled | ncol(co.ratio)==length(select.cells)){
        adj.mat = co.ratio
      }
      else{
        adj.mat = co.ratio[select.cells, select.cells]
      }
      gr = graph.adjacency(adj.mat, mode="undirected",weighted=TRUE)
      comm= cluster_louvain(gr)
      rm(gr)

      if(pass_louvain(modularity(comm), adj.mat)){
        tmp.cl = setNames(comm$membership,colnames(adj.mat))
        rm(comm)
        gc()
        if(length(unique(tmp.cl))==1){
          return(NULL)
        }
       
      }
      else{
        return(NULL)
      }
      rm(adj.mat)
      gc()
    }
    if(verbose){
      print(table(tmp.cl))
    }
    tmp.cl=merge_cl_by_co(tmp.cl, co.ratio=co.ratio, cl.mat=cl.mat[,names(tmp.cl)],diff.th)
    cell.cl.co.ratio= get_cell.cl.co.ratio(tmp.cl, co.ratio= co.ratio, cl.mat=cl.mat[,names(tmp.cl)])
    tmp.cells = sample_cells(tmp.cl, max.cl.size)
    norm.dat = get_logNormal(big.dat, tmp.cells)
    tmp= merge_cl(norm.dat=norm.dat, cl=tmp.cl, rd.dat.t=t(cell.cl.co.ratio), verbose=verbose,  de.param = de.param, max.cl.size= max.cl.size)
    rm(norm.dat)
    rm(tmp.cells)
    gc()
    
    if(is.null(tmp) | !is.list(tmp)) return(NULL)
    if (length(unique(tmp$cl))==1) return(NULL)
    markers=tmp$markers
    tmp.cl= tmp$cl
    tmp.cl = setNames(as.integer(tmp.cl),names(tmp.cl))
    rm(tmp)
    gc()
    if(length(unique(tmp.cl))==1) {
      return(NULL)
    }
    if(co.ratio.sampled){
      co.ratio = NULL
    }
    
    cell.cl.co.ratio= get_cell.cl.co.ratio(tmp.cl, co.ratio= co.ratio, cl.mat=cl.mat[,select.cells])
    if(length(tmp.cl) < length(select.cells)){
      cl = setNames(as.integer(colnames(cell.cl.co.ratio)[apply(cell.cl.co.ratio, 1, which.max)]), row.names(cell.cl.co.ratio))
    }
    else{
      cl = tmp.cl
    }
    rm(tmp.cl)
    gc()
    if(verbose){
      cat("Total:", length(cl), "\n")
      print(table(cl))
    }
  }
  cl.size = sort(table(cl))
  uncertain.cells=sapply(names(cl.size), function(i){
    tmp.cells=names(cl)[cl==i]
    sum(cell.cl.co.ratio[tmp.cells, as.character(i)] < 1 - diff.th)
  })
  rm(cell.cl.co.ratio)
  gc()
  if(verbose){
    print("Uncertain cells")
    print(uncertain.cells)
  }
  select.cl = names(uncertain.cells)[uncertain.cells > de.param$min.cells]
  if(is.null(select.cl)){
    return(list(cl=cl, markers=markers))
  }
  n.cl=max(cl)
  new.cl=cl
  for(i in select.cl){
    tmp.prefix= paste0(prefix, ".", i)
    tmp.cells=names(cl)[cl==i]
    result= iter_consensus_clust_big(cl.list=cl.list, co.ratio=co.ratio, cl.mat = cl.mat, big.dat = big.dat, select.cells=tmp.cells, prefix=tmp.prefix,  diff.th =diff.th, method=method, de.param = de.param, verbose=verbose, max.cl.size=max.cl.size, merge.type=merge.type)
    if(is.null(result)){
      next
    }
    tmp.cl= result$cl
    new.cl[names(tmp.cl)] = tmp.cl + n.cl
    n.cl = max(new.cl)
    markers=union(markers, result$markers)
    rm(result)
    gc()
  }
  cl=new.cl
  cl = setNames(as.integer(as.factor(cl)), names(cl))
  return(list(cl=cl, markers=markers))
}

map_result <- function(big.dat, cl.dat, test.cells, ncores=1, block.size=20000)
  {
    map <- function(big.dat, cols, cl.dat)
      {
        source("~/zizhen/My_R/scrattch.bigcat/R/cluster_big.R")
        test.dat = get_logNormal(big.dat, cols,sparse=FALSE)[row.names(cl.dat),,drop=F]
        test.cl.cor <- cor(test.dat, cl.dat)
        test.cl.cor[is.na(test.cl.cor)] <- 0
        max.cl.cor <- apply(test.cl.cor, 1, which.max)
        test.cl <- setNames(colnames(test.cl.cor)[max.cl.cor], colnames(test.dat))
        return(test.cl)
      }
    test.cl = big_dat_apply(big.dat, cols=test.cells, map, .combine="c",  ncores=ncores, block.size = block.size, cl.dat = cl.dat)
  }

run_consensus_clust_big <- function(big.dat=NULL, select.cells=big.dat$col_id, niter=100, iter=1:niter, sample.frac=0.8,  mc.cores=1, co.result = NULL,  output_dir="subsample_result",de.param=de_param(),merge.type=c("undirectional","directional"), override=FALSE, init.result=NULL, method="auto",max.cl.size = 200,...)
{
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  all.cells=select.cells
  if(!is.null(init.result)){
    all.cells= intersect(all.cells, names(init.result$cl))
  }
  if(is.null(co.result)){
    for(i in iter){
      outfile= file.path(output_dir, paste0("result.",i,".rda"))
      if(file.exists(outfile)& !override){
        next
      }
      if(!is.null(sample.frac)){
        select.cells=sample(all.cells, round(length(all.cells)*sample.frac))
      }
      else{
        select.cells = all.cells
      }
      save(select.cells, file=file.path(output_dir, paste0("cells.",i,".rda")))
    }
    
    run <- function(i,all.cells, ...){
      library(scrattch.hicat)
      source("~/zizhen/My_R/scrattch.bigcat/R/cluster_big.R")
      source("~/zizhen/My_R/scrattch.bigcat/R/consensusCluster_big.R")      
      prefix = paste("iter",i,sep=".")
      print(prefix)
      outfile= file.path(output_dir, paste0("result.",i,".rda"))
      if(file.exists(outfile)& !override){
        return(NULL)
      }
      load(file.path(output_dir, paste0("cells.",i,".rda")))
      result <- iter_clust_big(big.dat = big.dat,  select.cells=select.cells, prefix=prefix, de.param = de.param, merge.type=merge.type, result=init.result, ...)
      #save(result, file=outfile)
      #load(outfile)
      ###Test on remaining cells
      test.cells = setdiff(all.cells, names(result$cl))
      if(length(test.cells)>0){
        cl= result$cl              
        train.cells = sample_cells(cl, max.cl.size)
        train.norm.dat =  get_logNormal(big.dat, train.cells)[result$markers,,drop=F]
        cl.dat <- get_cl_means(train.norm.dat, cl[train.cells])
        test.cl = map_result(big.dat, cl.dat, test.cells)
        result$test.cl=test.cl
      }
      save(result, file=outfile)
    }
    if (mc.cores==1){
      sapply(iter, function(i){run(i,all.cells,...)})
    }
    else{
      require(foreach)
      require(doParallel)
      cl <- makeCluster(mc.cores)
      registerDoParallel(cl)
      foreach(i=iter, .combine='c') %dopar% run(i,all.cells,...)
      stopCluster(cl)
    }  
    result.files= file.path(output_dir, dir(output_dir, pattern="result.*.rda"))
    co.result <- collect_subsample_cl_matrix_big(big.dat,result.files,all.cells)
  }
  load(result.files[[1]])
  cl.size = table(result$cl)

  consensus.result = iter_consensus_clust_big(cl.list=co.result$cl.list, cl.mat = co.result$cl.mat,  big.dat=big.dat, select.cells=all.cells, de.param = de.param, merge.type=merge.type, method=method, result=init.result)
  refine.result = refine_cl(consensus.result$cl, cl.mat = co.result$cl.mat, tol.th=0.01, confusion.th=0.6, min.cells= de.param$min.cells)
  markers = consensus.result$markers
  
  #result <- iter_clust_big(big.dat = big.dat, select.cells=all.cells, de.param = de.param, merge.type=merge.type, result=init.result, binSize=binSize,...)
  #cl=merge_cl_by_co(result$cl, co.ratio=NULL, cl.mat=co.result$cl.mat, diff.th=0.25)
  #refine.result = refine_cl(cl, cl.mat = co.result$cl.mat, tol.th=0.01, confusion.th=0.6, min.cells=de.param$min.cells)
  #markers=result$markers      

  cl = refine.result$cl
  tmp.cells = sample_cells(cl, max.cl.size)
  norm.dat = get_logNormal(big.dat, tmp.cells)
  merge.result= merge_cl(norm.dat=norm.dat, cl=cl, rd.dat.t=norm.dat[markers,], de.param = de.param, merge.type=merge.type, return.markers=FALSE)
  return(list(co.result=co.result, cl.result=merge.result))
}


collect_subsample_cl_matrix_big <- function(big.dat,result.files,all.cells, max.cl.size=NULL,mc.cores=1)
{
  select.cells=c()
  run <- function(f){
    library(scrattch.hicat)
    #source("~/zizhen/My_R/scrattch.hicat_big/R/cluster.R")
    print(f)
    tmp=load(f)
    all.cl = with(result, setNames(as.character(cl),names(cl)))
    if(is.null(result$test.cl)){
      test.cells = setdiff(all.cells, names(result$cl))
      if(length(test.cells)>0){
        train.cells = sample_cells(result$cl, 100)
        train.norm.dat =  get_logNormal(big.dat, train.cells)[result$markers,]
        cl.dat <- get_cl_means(train.norm.dat, result$cl[train.cells])   
        result$test.cl = map_result(big.dat, cl.dat, test.cells)
        save(result, file=f)
        all.cl = with(result, c(all.cl, setNames(as.character(test.cl), names(test.cl))))
      }     
    }
    return(all.cl)
  }
  if (mc.cores==1){
    cl.list=sapply(result.files, function(f){run(f)},simplify=F)
  }
  else{
    require(foreach)
    require(doParallel)
    cl <- makeCluster(mc.cores)
    registerDoParallel(cl)
    cl.list= foreach(f=result.files, .combine='list', .multicombine=TRUE) %dopar% run(f)
    names(cl.list) = result.files
    stopCluster(cl)
  }
  save(cl.list, file="cl.list.rda")
  if(!is.null(max.cl.size)){
    select.cells= sample_cl_list(cl.list, max.cl.size=max.cl.size)
  }
  else{
    select.cells= all.cells
  }
  cl.mat = compile_cl_mat(cl.list, select.cells)
  return(list(cl.list=cl.list, cl.mat = cl.mat))
}

