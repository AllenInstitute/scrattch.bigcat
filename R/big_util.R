create_big.dat_fbm <- function(row.id, col.id, logNormal=TRUE, backingfile=file.path(getwd(), "fbm"),...)
  {
    library(bigstatsr)
    m = FBM(nrow=length(row.id),ncol=length(col.id),backingfile=backingfile, type="float", ...)
    big.dat = list(fbm=m, row_id = row.id, col_id = col.id)
    big.dat$logNormal = logNormal
    return(big.dat)
  }

convert_big.dat_fbm <- function(mat, logNormal=TRUE, backingfile=file.path(getwd(), "fbm"),...)
  {
    library(bigstatsr)
    m = FBM(nrow=nrow(mat),ncol=ncol(mat),backingfile=backingfile, type="float",...)
    ind_nozero <- Matrix::which(mat != 0L, arr.ind = TRUE)
    m[ind_nozero] <- mat[ind_nozero]
    big.dat = list(fbm=m, row_id = row.names(mat), col_id = colnames(mat))
    big.dat$logNormal = logNormal
    return(big.dat)
  }

###Only works for fbm
append_big.dat_fbm <- function(big.dat, mat)
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

###Only works for fbm
convert_mat_list_big.dat_fbm <- function(mat.fn, backingfile=NULL, ...)
  {
    if(is.null(backingfile)){
      backingfile=file.path(getwd(), paste0(Sys.Date(),"_fbm"))
    }
    big.dat=NULL
    for(i in 1:length(mat.fn)){
      load(mat.fn[i])      
      if(i==1){
        big.dat =  convert_big.dat(mat.list[[1]],backingfile=backingfile, ...)
      }
      else{
        big.dat = append_big.dat(big.dat, mat)
      }
    }
    return(big.dat)
  }



###Only works for fbm
reorder_big.dat_fbm <- function(big.dat, new.cols,backingfile=NULL)
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

get_cols <-  function(big.dat, cols, rows=big.dat$row_id, ...)
  {
    if(!is.list(big.dat)){
      mat = big.dat[rows, cols]
    }
    else if(big.dat$type=="fbm"){
      mat=get_cols_fbm(big.dat, cols, rows=rows, ...)
    }
    else if(big.dat$type=="parquet"){
      mat = get_cols_parquet(big.dat, cols, rows=rows,...)
    }
    else if(big.dat$type=="parquet_dens"){
      mat = get_cols_parquet_dense(big.dat, cols, rows=rows,...)
    }
    return(mat)
  }

get_cols_fbm<- function(big.dat, cols, rows=big.dat$row_id, keep.col=TRUE, sparse=TRUE)
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
    is.ordered=all(ord == 1:length(ord))    
    if(keep.col & !is.ordered){
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
    if(!identical(rows, big.dat$row_id)){
      mat = mat[rows,]
    }
    return(mat)
  }


get_logNormal <- function(big.dat, cols, rows=NULL, ...)
  {
    if(big.dat$logNormal){
      mat = get_cols(big.dat, cols, rows, ...)    
      norm.dat = mat
    } 
    else{
      mat = get_cols(big.dat, cols, ...)    
      norm.dat=logCPM(mat)
      if(!is.null(rows)){
        norm.dat = norm.dat[rows,,drop=F]
      }
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


init_mat_parqeut_dense <- function(dir, col.id, row.id, val=0, col.bin.size=50000, row.bin.size= 500, mc.cores=10)
  {
    library(arrow)
    if(!dir.create(dir)){
      stop("Create a new directory")
    }    
    row.df = data.frame(row_id=1:length(row.id), row_name=row.id)
    row.df$row_bin = ceiling((1:nrow(row.df))/row.bin.size)
    row.df = row.df %>% group_by(row_bin) %>% mutate(row_bin_id=row_id - min(row_id)+1) %>% ungroup
    
    col.df = data.frame(col_id=1:length(col.id), col_name=col.id)
    col.df$col_bin = ceiling((1:nrow(col.df))/col.bin.size)
    col.df = col.df %>% group_by(col_bin) %>% mutate(col_bin_id=col_id - min(col_id)+1) %>% ungroup

    write_parquet(row.df, sink=file.path(dir,"row.parquet"))
    write_parquet(col.df, sink=file.path(dir,"col.parquet"))
    
    mat = matrix(val, nrow=row.bin.size, ncol=col.bin.size)
    mat.df = as.data.frame(mat)        
    for(a=unique(col.df$col_bin)){
      d = file.path(dir,a)
      dir.create(d)
      for(b=unique(row.df$row_bin)){
        dir.create(file.path(dir,a,b))
        write_parquet(mat.df, file.path(d, "data.parquet"))        
      }
    }
    big.dat = list(type="parquet.dense", parquet.dir=dir,row_id = row.id, col_id = col.id,row.df=row.df, col.df = col.df)
    return(big.dat)
  }        

update_mat_parqeut_dense <- function(big.dat, mat, mc.cores=10)
  {
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    library(arrow)
    library(dplyr)
    row.df = read_parquet(file.path(big.dat$parquet.dir, "row.parquet")) %>% filter(row_names %in% row.names(mat))
    col.df = read_parquet(file.path(big.dat$parquet.dir, "col.parquet")) %>% filter(col_names %in% col.names(mat))
    foreach(a=unique(col.df$col_bin)) %:%
      foreach(b=unique(row.df$row_bin)) %dopar% {
        fn = file.path(big.dat$parquet.dir, a,b, "data.parquet")
        org.mat =  read_parquet(fn)
        tmp.col.df = col.df %>% filter(col_bin==a)
        tmp.row.df = row.df %>% filter(row_bin==b)        
        tmp.mat = mat[tmp.row.df$row_names, tmp.col.df$col_names]              
        org.mat[tmp.row.df$row_id, tmp.col.df$col_id] = tmp.mat
        write_parquet(org.mat, fn)
      }    
  }


get_cols_parquet_dense <- function(big.dat, cols, rows, mc.cores=5)
  {
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    library(arrow)
    library(dplyr)

    if(is.null(big.dat$row.df)){
      row.df = read_parquet(file.path(big.dat$parquet.dir, "row.parquet"))
    }
    row.df = row.df %>% filter(row_names %in% rows)
    if(is.null(big.dat$col.idf)){
      col.df = read_parquet(file.path(big.dat$parquet.dir, "col.parquet"))
    }
    col.df = col.df %>% filter(col_names %in% cols)    
    mat=foreach(a=unique(col.df$col_bin), .combine="cbindlist") {
      tmp.col.df = col.df %>% filter(col_bin==a)      
      mat = foreach(b=unique(row.df$row_bin),.combine="rbindlist") %dopar% {
        tmp.row.df = row.df %>% filter(row_bin==b)
        fn = file.path(big.dat$parquet.dir, a,b, "data.parquet")
        mat =  read_parquet(fn)
        mat=mat[tmp.row.df$row_id, tmp.col.df$col_id]
        row.names(mat) = tmp.row.df$row_names
        colnames(mat) = tmp.coldf$col_names
        list(mat)
      }
      list(mat)
    }
    mat = as.matrix(mat)
    return(mat)    
  }



convert_big.dat_parqeut <- function(mat, dir="./",parquet.dir=file.path(dir,"norm.dat_parquet"), col.fn=file.path(dir,"samples.parquet"), row.fn=file.path(dir,"gene.parquet"),col.bin.size=50000, row.bin.size=500,logNormal=TRUE,mc.cores=10)
  {
    if(!dir.exists(dir)){
      dir.create(dir)
    }
    if(!dir.exists(parquet.dir)){
      dir.create(parquet.dir)
    }
    col.df = data.table(col_name=colnames(mat))
    col.df$col_id = 1:nrow(col.df)
    col.df$col_bin = ceiling((1:nrow(col.df))/col.bin.size)
    write_parquet(col.df, col.fn)

    
    row.id = row.names(mat.list[[1]])
    row.df = data.frame(row_id=1:length(row.id), row_name=row.id)
    row.df$row_bin = ceiling((1:nrow(row.df))/row.bin.size)
    write_parquet(row.df, row.fn)

    
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    tmp <- foreach(bin = 1:max(col.df$col_bin), .combine="c") %dopar% {            
      col.id = col.df %>% filter(col_bin==bin) %>% pull(col_id)
      tmp.mat = mat[,col.id]
      i = tmp.mat@i
      x = tmp.mat@x
      p = tmp.mat@p      
      j = p[-1] - p[-length(p)]
      j = rep(col.id-1, j)
      df = data.frame(i=i, j=j, x=x)
      df$row_bin = ceiling((i+1)/row.bin.size)
      path = file.path(parquet.dir, paste0("col_bin=",bin))
      dir.create(path)
      write_dataset(df, path, partition=c("row_bin"))      
    }
    big.dat = list(type = "parquet", row_id = row.df$row_name, col_id = col.df$col_name, logNormal=logNormal,col.fn = col.fn, row.fn = row.fn, parquet.dir = parquet.dir)    
    return(big.dat)    
  }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param mat 
##' @param big.dat 
##' @param col.bin.size 
##' @param row.bin.size 
##' @return 
##' @author Zizhen Yao
append_big.dat_parquet <- function(mat, big.dat, col.bin.size=50000, row.bin.size=500, mc.cores=10)
  {
    parquet.dir = big.dat$parquet.dir
    col.df = read_parquet(big.dat$col.fn)
    max.bin = max(col.df$col_bin)
    max.id = max(col.df$col_id)
    
    new.col.df = data.table(col_name=colnames(mat),col_id = 1:ncol(mat))
    new.col.df$col_bin = ceiling((1:nrow(new.col.df))/col.bin.size) + max.bin    
    row.df = read_parquet(big.dat$row.fn)
    row.df$i = row.df$row_id - 1
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    tmp <- foreach(bin = unique(new.col.df$col_bin), .combine="c") %dopar% {            
      col.id = col.df %>% filter(col_bin==bin) %>% pull(col_id)
      tmp.mat = mat[,col.id]
      i = tmp.mat@i
      x = tmp.mat@x
      p = tmp.mat@p      
      j = p[-1] - p[-length(p)]
      j = rep(col.id-1, j) + max.id
      df = data.frame(i=i, j=j, x=x)
      df = df %>% left_join(row.df[,c("i","row_bin")],by="i")
      path = file.path(parquet.dir, paste0("col_bin=",bin))
      dir.create(path)
      print(path)
      write_dataset(df, path, partition=c("row_bin"))      
    }
    new.col.df$col_id = new.col.df$col_id + max.id
    col.df = rbindlist(list(col.df, new.col.df))    
    write_parquet(col.df, sink=big.dat$col.fn)
    big.dat$col_id = c(big.dat$col_id, new.col.df$col_name)    
    return(big.dat)    
  }

reorder_big.dat_parquet <- function(big.dat, cols, parquet.dir, col.parquet.fn, col.bin.size=50000,mc.cores=10)
  {
    big.dat = init_big.dat_parquet(big.dat)
    row.df = big.dat$row.df
    row.df$i = row.df$row_id - 1
    new.col.df = big.dat$col.df
    new.col.df$new_col_id = match(new.col.df$col_name, cols)
    new.col.df = new.col.df %>% filter(!is.na(new_col_id)) %>% arrange(new_col_id)
    new.col.df$new_col_bin = ceiling(new.col.df$new_col_id/col.bin.size)
    dir.create(parquet.dir)
    nbin = max(new.col.df$new_col_bin)
    tmp <- foreach(bin = 1:nbin, .combine="c") %dopar% {
      tmp.col.df = new.col.df %>% filter(new_col_bin==bin) %>% select(col_name,new_col_id)      
      tmp.mat = get_cols(big.dat, tmp.col.df$col_name)
      i = tmp.mat@i
      x = tmp.mat@x
      p = tmp.mat@p
      j = p[-1] - p[-length(p)]
      j = rep(tmp.col.df$new_col_id-1, j)
      df = data.frame(i=i, j=j, x=x)
      df = df %>% left_join(row.df[,c("i","row_bin")],by="i")
      path = file.path(parquet.dir, paste0("col_bin=",bin))
      dir.create(path)
      print(path)
      write_dataset(df, path, partition=c("row_bin"))
    }
    new.col.df = new.col.df[,c("new_col_id","col_name","new_col_bin")]    
    colnames(new.col.df)=c("col_id","col_name","col_bin")    
    write_parquet(new.col.df, sink=col.parquet.fn)
    big.dat$col.df = new.col.df
    big.dat$col.fn = col.parquet.fn
    big.dat$col_id = cols
    big.dat$parquet.dir = parquet.dir
    big.dat$ds = open_dataset(big.dat$parquet.dir)
    return(big.dat)    
  }

fbm_to_parquet <- function(big.dat, cols=big.dat$col_id, rows=big.dat$row_id, dir="./",parquet.dir=file.path(dir,"norm.dat_parquet"), col.fn=file.path(dir,"samples.parquet"), row.fn=file.path(dir,"gene.parquet"), row.bin.size = 500,    col.bin.size = 50000)
  {
    if(!dir.exists(dir)){
      dir.create(dir)
    }
    if(!dir.exists(parquet.dir)){
      dir.create(parquet.dir)
    }
    if(is.character(cols)){
      col.id = sort(match(cols, big.dat$col_id))
    }
    else{
      col.id = sort(cols)
    }
    if(is.character(rows)){
      row.id = sort(match(rows, big.dat$row_id))
    }
    else{
      row.id = sort(rows)
    }
    col.bin.size = max(length(col.id)/1000, col.bin.size)
    col.df = data.frame(col_id=col.id, col_name=big.dat$col_id[col.id])
    col.df$col_bin = ceiling((1:nrow(col.df))/col.bin.size)
    write_parquet(col.df, col.fn)

    row.df = data.frame(row_id=row.id, row_name=big.dat$row_id[row.id])
    row.df$row_bin = ceiling((1:nrow(row.df))/row.bin.size)
    write_parquet(row.df, row.fn)

    for(bin in 1:max(col.df$col_bin)){
      col.id = col.df %>% filter(col_bin==bin) %>% pull(col_id)
      cat("bin", bin, length(col.id),"\n")
      mat = get_cols(big.dat, col.id,sparse=TRUE)
      i = mat@i
      x = mat@x
      p = mat@p      
      j = p[-1] - p[-length(p)]
      j = rep(col.id-1, j)
      df = data.frame(i=i, j=j, x=x)
      df$row_bin = ceiling((i+1)/row.bin.size)
      path = file.path(parquet.dir, paste0("col_bin=",bin))
      dir.create(path)
      write_dataset(df, path, partition=c("row_bin"))      
    }
    new.big.dat = big.dat
    new.big.dat$type = "parquet"
    new.big.dat$col.fn = col.fn
    new.big.dat$row.fn = row.fn
    new.big.dat$parquet.dir = parquet.dir
    new.big.dat$fbm = NULL
    return(new.big.dat)    
  }


filter_cols <- function(col.df, cols)
  {
    if(is.character(cols)){      
      col.df = col.df %>% filter(col_name %in% cols)
      select.j = sort(col.df$col_id - 1)
      col.bin = col.df %>% pull(col_bin) %>% unique 
    }
    else{
      col.df = col.df %>% filter(col_id %in% cols)
      select.j = sort(col.df$col_id - 1)
      col.bin = col.df %>% pull(col_bin) %>% unique
      cols = data.table(col_id=cols) %>% left_join(col.df[,c("col_id","col_name")],by="col_id") %>% pull(col_name) 
    }
    return(list(cols=cols, col.bin=col.bin, select.j=select.j))
  }


init_big.dat_parquet <- function(big.dat.parquet)
  {
    library(arrow)
    library(data.table)
    dir = big.dat.parquet$parquet.dir
    col.fn = big.dat.parquet$col.fn
    row.fn = big.dat.parquet$row.fn
    
    ds = open_dataset(dir)
    col.df = read_parquet(col.fn)
    row.df = read_parquet(row.fn)
    big.dat.parquet$ds = ds
    big.dat.parquet$col.df = col.df
    big.dat.parquet$row.df = row.df    
    return(big.dat.parquet)
  }

get_cols_parquet <- function(big.dat.parquet, cols, rows=NULL,keep.col=FALSE, sparse=TRUE, mc.cores=5)
  {
    library(data.table)
    library(arrow)
    library(dplyr)
    if(is.null(big.dat.parquet$ds)){
      big.dat.parquet = init_big.dat_parquet(big.dat.parquet)
    }
    ds = big.dat.parquet$ds
    col.df = big.dat.parquet$col.df
    row.df = big.dat.parquet$row.df
    
    if(is.character(cols)){      
      col.df = col.df %>% filter(col_name %in% cols)       
    }
    else{
      col.df = col.df %>% filter(col_id %in% cols)
      cols = data.table(col_id=cols) %>% left_join(col.df[,c("col_id","col_name")],by="col_id") %>% pull(col_name)    
    }
    select.j = sort(col.df$col_id - 1)
    col.bin = col.df %>% pull(col_bin) %>% unique 

    #col.df$new_j = 1:nrow(col.df)
    col.df$new_j = match(col.df$col_name, cols)
    col.id = rep(0, max(col.df$col_id))
    col.id[col.df$col_id] = col.df$new_j

    if(!is.null(rows)){      
      if(is.character(rows)){      
        row.df = row.df %>% filter(row_name %in% rows)
      }
      else{
        row.df = row.df %>% filter(row_id %in% rows)
        rows = data.table(row_id=rows) %>% left_join(row.df[,c("row_id","row_name"),by="row_id"]) %>% pull(row_name)
      }
    }
    else{
      rows=big.dat.parquet$row_id
    }
    select.i = sort(row.df$row_id) - 1
    row.bin = row.df %>% pull(row_bin) %>% unique
    
    #row.df$new_i = 1:nrow(row.df)
    row.df$new_i= match(row.df$row_name, rows)
    row.id = rep(0, max(row.df$row_id))
    row.id[row.df$row_id] = row.df$new_i

    if(as.numeric(length(select.j))*length(select.i) > 2^32){
      stop("Matrix too big")
    }

    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    tmp <- foreach(c.id = col.bin, .combine="c") %:%      
      foreach(r.id =row.bin, .combine="c") %dopar% {
        mat.df = ds %>% filter(r.id==row_bin & c.id==col_bin  & j %in% select.j & i %in% select.i) %>% select(i,j,x) %>% collect()
        list(mat.df)
      }
    
    mat.df = rbindlist(tmp)
    #mat.df = ds %>% filter(row_bin %in% row.bin & col_bin %in% col.bin & j %in% select.j & i %in% select.i) %>% select(i,j,x) %>% collect()
    
    library(Matrix)
    mat = sparseMatrix(i = row.id[mat.df$i+1], j=col.id[mat.df$j+1], x=mat.df$x, dims=c(nrow(row.df),nrow(col.df)))
    #colnames(mat) = col.df$col_name
    #row.names(mat) = row.df$row_name
    colnames(mat) = cols
    row.names(mat) = rows
    if(keep.col){
      mat = mat[,cols]
    }
    if(!sparse){
      mat = as.matrix(mat)
    }
    return(mat)
  }




get_cl_stats_parquet <- function(big.dat.parquet, cl, mc.cores=20,stats=c("means","present","sqr_means"), return.matrix=TRUE)
  {
    library(data.table)
    library(arrow)
    library(dplyr)
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)

    if(is.null(big.dat.parquet$ds)){
      mat.dir = big.dat.parquet$parquet.dir
      ds = open_dataset(mat.dir)
      col.fn = big.dat.parquet$col.fn
      row.fn = big.dat.parquet$row.fn
      col.df = read_parquet(col.fn)
      row.df = read_parquet(row.fn)
    }
    else{
      ds = big.dat.parquet$ds
      row.df = big.dat.parquet$row.df
      col.df = big.dat.parquet$col.df
    }
    cl.anno = data.table(cl=as.factor(cl), col_name=names(cl)) %>% left_join(col.df, by="col_name")
    cl.anno$j = cl.anno$col_id - 1
    cl.size = cl.anno %>% group_by(cl) %>% summarize(size=n(),.groups = 'drop')
    cl.l = levels(cl.anno$cl)
    ncl = length(cl.l)
    col.bin = cl.anno %>% select(col_bin) %>% distinct() %>% pull(col_bin)
    row.bin = row.df %>% pull(row_bin) %>% unique
    tmp.dir = tempdir()
    if(dir.exists(tmp.dir)){
      unlink(tmp.dir,recursive=TRUE, force=TRUE)
    }
    dir.create(tmp.dir)
    
    ###for each cluster, find markers that discriminate it from other types
    tmp <-  foreach(r.id = row.bin, .combine="c") %:%      
      foreach(c.id =col.bin, .combine="c") %dopar% {        
        tmp=ds %>% filter(col_bin== c.id & row_bin == r.id) %>% collect() %>%
          right_join(cl.anno %>% filter(col_bin==c.id) %>% select(j, cl),by="j") %>%
            group_by(cl,i) %>% filter(!is.na(i))
        cl.stats=NULL
        if("means" %in% stats){
          cl.stats <- tmp %>% summarize(cl.sum=sum(x),.groups = 'drop')
        }
        if("present" %in% stats){
          tmp1=tmp %>% summarize(cl.present.sum=sum(x>0),.groups = 'drop')
          if(is.null(cl.stats)){
            cl.stats <- tmp1
          }
          else{
            cl.stats <- cl.stats %>% left_join(tmp1, by=c("cl","i"))
          }
        }
        if("sqr_means" %in% stats){          
          tmp1 = tmp %>% summarize(cl.sqr.sum=sum(x^2),.groups = 'drop')
          if(is.null(cl.stats)){
            cl.stats <- tmp1
          }
          else{
            cl.stats <- cl.stats %>% left_join(tmp1, by=c("cl","i"))
          }
        }
        #cl.sum=sum(x),cl.present.sum = sum(x > 0), cl.sqr.sum=sum(x^2))
        idx=paste0(r.id, "_", c.id)
        write_parquet(as.data.frame(cl.stats), sink=file.path(tmp.dir, paste0(idx,".parquet")))
      }
    #print(tmp.dir)
    fn = dir(tmp.dir)
    rm(cl.anno)
    gc()
    result.df <-  foreach(r.id=row.bin, .combine="c") %dopar% {
      select.fn  = file.path(tmp.dir,grep(paste0("^",r.id,"_"), fn,value=T))
      df = open_dataset(select.fn) %>% collect() %>% group_by(cl, i)
      cl.stats = NULL
      if("means" %in% stats){
        cl.stats <- df %>% summarize(sum= sum(cl.sum),.groups = 'drop')        
      }
      if("present" %in% stats){
        tmp <- df %>% summarize(present.sum= sum(cl.present.sum),.groups = 'drop')
        if(is.null(cl.stats)){
          cl.stats <- tmp          
        }
        else{
          cl.stats = cl.stats %>% left_join(tmp,by=c("cl","i"))
        }
      }
      if("sqr_means" %in% stats){
        tmp <- df %>% summarize(sqr_sum= sum(cl.sqr.sum),.groups = 'drop')
        if(is.null(cl.stats)){
          cl.stats <- tmp          
        }
        else{
          cl.stats = cl.stats %>% left_join(tmp, by=c("cl","i"))
        }
      }      
      list(cl.stats)
    }
    
    cl.stats = rbindlist(result.df)
    cl.stats = cl.stats %>% left_join(cl.size, by="cl")    
    if("means" %in% stats){
      cl.stats <- cl.stats %>% mutate(means= sum/size)
    }
    if("present" %in% stats){
      cl.stats <- cl.stats %>% mutate(present = present.sum/size)
    }
    if("sqr_means" %in% stats){
      cl.stats <- cl.stats %>% mutate(sqr_means = sqr_sum/size)
    }
    
    if(return.matrix){
      ngenes = nrow(row.df)
      coor =  (as.integer(cl.stats$cl) - 1)* ngenes  + cl.stats$i + 1
      mat.list= list()
      for(x in stats){
        mat = matrix(0, nrow=ngenes,ncol=ncl)
        row.names(mat) = big.dat.parquet$row_id
        colnames(mat) = cl.l
        mat[coor] = cl.stats[[x]]
        mat = as.matrix(mat, ncol=ncl)
        mat.list[[x]] = mat
      }
      return(mat.list)
    }
    return(cl.stats)    
  }

get_cl_stats_big <- function(big.dat, cl, max.cl.size=200,stats=c("means"),...)
  {
    sampled.cells = sample_cells(cl, max.cl.size)
    cl = cl[sampled.cells]
    if(big.dat$type=="fbm"){
      get_cl_stats_fbm(big.dat, cl, stats=stats,...)
    }
    else{
      get_cl_stats_parquet(big.dat, cl, stats=stats,...)
    }
  }


get_big_bins <- function(big.dat, cols, block.size=10000)
  {
    if(is.character(cols)){
      cols = match(cols, big.dat$col_id)
    }
    cols = sort(cols)    
    bins=NULL
    if(big.dat$type=="parquet"){
      library(arrow)
      library(dplyr)
      if(is.null(big.dat$col.df)){
        col.df = read_parquet(big.dat$col.fn)
      }
      else{
        col.df = big.dat$col.df
      }
      col.df = col.df %>% filter(col_id %in% cols)
      col.df$new_bin = ceiling(1:nrow(col.df)/block.size)
      match.bin=col.df %>% group_by(col_bin) %>% summarize(match_bin=min(new_bin))
      col.df = col.df %>% left_join(match.bin, by="col_bin")
      bins = with(col.df, split(col_id, match_bin))
    }
    else{
      bins = split(cols, ceiling(1:length(cols)/block.size))
    }
  }


big_dat_apply <- function(big.dat, cols, FUN, .combine="c",  mc.cores=1, block.size=10000,...)
  {    
    require(foreach)    
    library(parallel)    
    require(doMC)
    bins = get_big_bins(big.dat, cols, block.size=block.size)
    #print(length(bins))
    mc.cores = min(mc.cores, length(bins))
    registerDoMC(cores=mc.cores)
    res = foreach(bin = bins, .combine=.combine) %dopar% FUN(big.dat, bin, ...)
  }


build_train_index <- function(cl.dat, method= c("Annoy.Cosine","cor","Annoy.Euclidean"),fn=tempfile(fileext=".idx"))
  {
    library(BiocNeighbors)
    method = method[1]
    ref.dat = Matrix::t(cl.dat)
    if(method=="cor"){
      ref.dat = ref.dat - rowMeans(ref.dat)
      ref.dat = l2norm(ref.dat,by = "row")
    }
    if (method=="Annoy.Cosine"){
      ref.dat = l2norm(ref.dat,by = "row")
    }
    index= buildAnnoy(ref.dat, fname=fn)
    return(index)    
  }

get_knn_batch_big <- function(big.dat, ref.dat, select.cells,block.size=10000, mc.cores,k, method="cor", dim=NULL, return.distance=FALSE, index=NULL, clear.index=FALSE, ntrees=50,transposed=TRUE)     
  {
    
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
      gc()
    }    
    result = big_dat_apply(big.dat, cols=select.cells, .combine=fun, mc.cores=mc.cores, block.size=block.size, FUN = function(big.dat, bin,...){
      {
        dat = get_logNormal(big.dat, bin, rows=row.names(ref.dat))
        knn=get_knn(dat=dat, ref.dat=ref.dat, k=k, method=method, dim=dim,return.distance=return.distance,index=index,...)
        rm(dat)
        gc()
        knn
      }
    })    
    if(clear.index){
      cleanAnnoyIndex(index)
    }
    else{
      if(!is.list(result)){
        result = list(result)
      }
      result$ref.index=index
    }
    return(result)
  }

map_cells_knn_big <- function(big.dat, cl.dat, select.cells, train.index=NULL, method = c("Annoy.Cosine","cor"), block.size=10000, mc.cores=10)
  {
    library(bigstatsr)
    cl.knn =  get_knn_batch_big(big.dat, cl.dat, select.cells=select.cells, k=1, index=train.index, method=method, transposed=TRUE, block.size=block.size, mc.cores=mc.cores,return.distance=TRUE)
    knn.index = cl.knn[[1]]
    knn.dist = cl.knn[[2]]
    map.df = data.frame(sample_id=row.names(knn.index), cl = colnames(cl.dat)[knn.index], dist = knn.dist)
    return(map.df)
  }

    
