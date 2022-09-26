create_big.dat_fbm <- function(row.id, col.id, logNormal=TRUE, backingfile=file.path(getwd(), "fbm"),...)
  {
    library(bigstatsr)
    m = FBM(nrow=length(row.id),ncol=length(col.id),backingfile=backingfile, type="float", ...)
    big.dat = list(fbm=m, row_id = row.id, col_id = col.id)
    big.dat$logNormal = logNormal
    big.dat$type="fbm"
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
    else if(big.dat$type=="parquet.dense"){
      mat = get_cols_parquet_dense(big.dat, cols, rows=rows,...)
    }
    return(mat)
  }

get_cols_fbm<- function(big.dat, cols, rows=big.dat$row_id, keep.col=TRUE, sparse=TRUE, mc.cores=0)
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

set_cols_fbm<- function(big.dat, mat)
  {
    library(Matrix)    
    cols = match(colnames(mat), big.dat$col_id)
    rows = match(row.names(mat), big.dat$row_id)
    big.dat$fbm[rows,cols]=mat
    return(NULL)
  }


get_logNormal <- function(big.dat, cols, rows=big.dat$row_id, ...)
  {
    if(big.dat$logNormal){
      mat = get_cols(big.dat, cols, rows, ...)    
      norm.dat = mat
    } 
    else{
      mat = get_cols(big.dat, cols, ...)    
      norm.dat=logCPM(mat)
    }
    if(!is.null(rows)){
      norm.dat = norm.dat[rows,,drop=F]
    }
    return(norm.dat)
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
    x = rep(val, row.bin.size*col.bin.size)
    fn = file.path(dir, "x.parquet")
    write_parquet(as.data.table(x), sink=fn)
    registerDoMC(cores=mc.cores)
    foreach(a = unique(col.df$col_bin)) %dopar%{
      d= file.path(dir,a)
      dir.create(d)
      for(b in unique(row.df$row_bin)){
        d2= file.path(d,b)
        print(d2)
        dir.create(d2)
        file.copy(fn, file.path(d2, "x.parquet"))
      }
    }
    file.remove(fn)
    big.dat = list(type="parquet.dense", parquet.dir=dir,row_id = row.id, col_id = col.id,row.df=row.df, col.df = col.df, row.bin.size=row.bin.size, col.bin.size = col.bin.size)
    return(big.dat)
  }        

update_mat_parquet_dense <- function(big.dat, mat, mc.cores=10)
  {
    library(parallel)    
    require(doMC)
    require(foreach)
    library(arrow)
    library(dplyr)
    row.df = read_parquet(file.path(big.dat$parquet.dir, "row.parquet"))
    row.df = row.df %>% filter(row_name %in% row.names(mat))
    col.df = read_parquet(file.path(big.dat$parquet.dir, "col.parquet"))
    col.df = col.df %>% filter(col_name %in% colnames(mat))
    mc.cores = min(mc.cores, length(unique(col.df$col_bin))*length(unique(row.df$row_bin)))
    registerDoMC(cores=mc.cores)
    foreach(a=unique(col.df$col_bin)) %:%
      foreach(b=unique(row.df$row_bin)) %dopar% {
        fn = file.path(big.dat$parquet.dir, a,b, "x.parquet")
        print(fn)
        org.mat =  read_parquet(fn)[["x"]]
        org.mat = matrix(org.mat, big.dat$row.bin.size, big.dat$col.bin.size)        
        tmp.col.df = col.df %>% filter(col_bin==a)
        tmp.row.df = row.df %>% filter(row_bin==b)        
        tmp.mat = mat[tmp.row.df$row_name, tmp.col.df$col_name]              
        org.mat[tmp.row.df$row_bin_id, tmp.col.df$col_bin_id] = tmp.mat
        df = data.table(x=as.vector(org.mat))
        write_parquet(df, fn)
        return(NULL)
      }
    return(NULL)
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
    else{
      row.df=big.dat$row.df
    }
    row.df = row.df %>% filter(row_name %in% rows)
    if(is.null(big.dat$col.df)){
      col.df = read_parquet(file.path(big.dat$parquet.dir, "col.parquet"))
    }
    else{
      col.df = big.dat$col.df
    }
    col.df = col.df %>% filter(col_name %in% cols)    
    mat=foreach(a=unique(col.df$col_bin), .combine="cbind") %dopar% {
      tmp.col.df = col.df %>% filter(col_bin==a)
      result = matrix(0, nrow=nrow(row.df), ncol=nrow(tmp.col.df),dimnames=list(row.df$row_name, tmp.col.df$col_name))
      for(b in unique(row.df$row_bin)){
        tmp.row.df = row.df %>% filter(row_bin==b)
        fn = file.path(big.dat$parquet.dir, a,b, "x.parquet")
        print(fn)
        mat = read_parquet(fn)$x
        mat = matrix(mat, nrow=big.dat$row.bin.size, ncol=big.dat$col.bin.size)
        result[tmp.row.df$row_name, tmp.col.df$col_name]=mat[tmp.row.df$row_bin_id, tmp.col.df$col_bin_id]
      }
      result
    }
    mat = mat[rows, cols]
    return(mat)    
  }



convert_big.dat_parquet <- function(mat, dir=getwd(),parquet.dir=file.path(dir,"norm.dat_parquet"), col.fn=file.path(dir,"samples.parquet"), row.fn=file.path(dir,"gene.parquet"),col.bin.size=50000, row.bin.size=500,logNormal=TRUE,mc.cores=10)
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

    
    #row.id = row.names(mat.list[[1]])
    row.id = rownames(mat)
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
    big.dat = list(type = "parquet", row_id = row.df$row_name, col_id = col.df$col_name, logNormal=logNormal,col.fn = col.fn, row.fn = row.fn, parquet.dir = parquet.dir,row.df=row.df, col.df = col.df)    
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
    if(is.null(big.dat$col.df)){
      col.df = read_parquet(big.dat$col.fn)
    }
    else{
      col.df = big.dat$col.df
    }
    max.bin = max(col.df$col_bin)
    max.id = max(col.df$col_id)    
    new.col.df = data.table(col_name=colnames(mat),col_id = 1:ncol(mat))
    new.col.df$col_bin = ceiling((1:nrow(new.col.df))/col.bin.size) + max.bin
    if(is.null(big.dat$row.df)){
      row.df = read_parquet(big.dat$row.fn)
    }
    else{
      row.df = big.dat$row.df
    }
    row.df$i = row.df$row_id - 1
    library(parallel)    
    require(doMC)
    require(foreach)
    bins = unique(new.col.df$col_bin)
    registerDoMC(cores= min(mc.cores, length(bins)))
    tmp <- foreach(bin = bins, .combine="c") %dopar% {            
      col.id = new.col.df %>% filter(col_bin==bin) %>% pull(col_id)
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
    big.dat$col.df = col.df
    return(big.dat)    
  }

convert_mat_list_big.dat_parquet <- function(mat.fn.list, dir=getwd(),parquet.dir=file.path(dir,"norm.dat_parquet"), col.fn=file.path(dir,"samples.parquet"), row.fn=file.path(dir,"gene.parquet"),col.bin.size=50000, row.bin.size=500,logNormal=TRUE,mc.cores=10)
  {
    load(mat.fn.list[[1]])
    cat(mat.fn.list[[1]],ncol(mat),"cells\n")
    big.dat = convert_big.dat_parquet(mat, dir=dir, parquet.dir = parquet.dir, col.fn=col.fn, row.fn=row.fn, col.bin.size=col.bin.size, row.bin.size=row.bin.size, logNormal=logNormal, mc.cores=mc.cores)
    if(length(mat.fn.list)>1){
      for(i in 2:length(mat.fn.list)){        
        load(mat.fn.list[[i]])
        cat(mat.fn.list[[i]],ncol(mat),"cells\n")
        big.dat=append_big.dat_parquet(mat, big.dat, col.bin.size=col.bin.size, row.bin.size=row.bin.size, mc.cores=mc.cores)
      }
    }
    return(big.dat)        
  }


reorder_big.dat_parquet <- function(big.dat, cols, dir=getwd(),parquet.dir=file.path(dir,"norm.dat_parquet"), col.fn=file.path(dir,"samples.parquet"), row.fn=file.path(dir,"gene.parquet"), col.bin.size=50000,mc.cores=10)
  {
    big.dat = init_big.dat_parquet(big.dat)
    if(!dir.exists(dir)){
      dir.create(dir)
    }
    if(!dir.exists(parquet.dir)){
      dir.create(parquet.dir)
    }
    file.copy(big.dat$row.fn, row.fn)
    row.df = big.dat$row.df
    row.df$i = row.df$row_id - 1
    new.col.df = big.dat$col.df
    new.col.df$new_col_id = match(new.col.df$col_name, cols)
    new.col.df = new.col.df %>% filter(!is.na(new_col_id)) %>% arrange(new_col_id)
    new.col.df$new_col_bin = ceiling(new.col.df$new_col_id/col.bin.size)
    
    dir.create(parquet.dir)
    nbin = max(new.col.df$new_col_bin)
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores= min(mc.cores, nbin))    
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
    write_parquet(new.col.df, sink=col.fn)
    big.dat$col.df = new.col.df    
    big.dat$col.fn = col.fn
    big.dat$row.fn = row.fn
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
    if(!dir.exists(dir)){
      stop(paste("Data directory ", dir, "does not exists\n"))
    }
    ds = open_dataset(dir)
    if(!file.exists(col.fn)){
      stop(paste("col.fn", col.fn, "does not exists\n"))
    }
    col.df = read_parquet(col.fn)
    if(!file.exists(row.fn)){
      stop(paste("row.fn", row.fn, "does not exists\n"))
    }
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
        rows = data.table(row_id=rows) %>% left_join(row.df[,c("row_id","row_name")],by="row_id")%>% pull(row_name)
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
    registerDoMC(cores=min(mc.cores, length(col.bin)))
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


big.dat_logNormal_parquet <- function(big.dat, parquet.dir, denom=10^6,mc.cores=10)
  {
    if(big.dat$logNormal){
      return(big.dat)
    }
    dir.create(parquet.dir)
    big.dat.new = big.dat
    big.dat = init_big.dat_parquet(big.dat)    
    ds = big.dat$ds
    library(parallel)    
    require(doMC)
    require(foreach)
    bins = big.dat$col.df %>% pull(col_bin) %>% unique
    registerDoMC(cores=min(mc.cores, length(bins)))
    tmp <- foreach(bin= bins, .combine="c") %dopar% {      
      map.df = ds %>% filter(col_bin ==bin) %>% collect()
      sf=map.df %>% group_by(j) %>% summarize(sf=sum(x)/denom)      
      map.df = map.df %>% left_join(sf) %>% mutate(x.norm = log2(x/sf+1))
      tmp.map.df = map.df %>% mutate(x=x.norm) %>% select(i,j,x,"row_bin")
      path = file.path(parquet.dir, paste0("col_bin=", bin))
      dir.create(path)
      print(path)
      write_dataset(tmp.map.df, path, partition = c("row_bin"))
      NULL
    }
    big.dat.new$parquet.dir = parquet.dir
    big.dat.new$ds = NULL
    big.dat.new$logNormal=TRUE
    return(big.dat.new)
  }


get_cl_stats_fbm <- function(big.dat, cl, max.cells=100000, stats=c("means"),genes.allowed=NULL,mc.cores=1)
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
      if(!is.null(genes.allowed)){
        select = row.names(dat) %in% genes.allowed
        dat = dat[select,,drop=FALSE]
      }
      result = sapply(stats, function(x){
          get_cl_stats(dat, cl=tmp.cl,stats=x)
      }, simplify=F)
    },simplify=F)
    cl.results = sapply(stats, function(x){
      do.call("cbind",sapply(tmp.results, function(result) result[[x]], simplify=F))
    },simplify=F)
    return(cl.results)    
  }

get_cl_stats_parquet <- function(big.dat.parquet, cl, mc.cores=20,stats=c("means","present","sqr_means"), genes.allowed=NULL, return.matrix=TRUE)
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
      if(!is.null(genes.allowed)){
        select = big.dat.parquet$row_id %in% genes.allowed
      }

      for(x in stats){
        mat = matrix(0, nrow=ngenes,ncol=ncl)
        row.names(mat) = big.dat.parquet$row_id
        colnames(mat) = cl.l
        mat[coor] = cl.stats[[x]]
        mat = as.matrix(mat, ncol=ncl)
        if(!is.null(genes.allowed)){
          mat = mat[select,,drop=FALSE]
        }
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
    cat("Total bins",length(bins),"\n")
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

get_knn_batch_big <- function(big.dat, ref.dat, select.cells,block.size=10000, mc.cores,k, method="cor", dim=NULL, return.distance=FALSE, index=NULL, clear.index=FALSE, ntrees=50)
  {
    library(BiocNeighbors)    
    if(return.distance){
      fun = "knn_combine"
    }
    else{
      fun = "rbind"
    }
    if(is.null(index) & method %in% c("Annoy.Euclidean", "Annoy.Cosine", "cor")) {     
      map.ref.dat = Matrix::t(ref.dat)
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
        knn=get_knn(dat=dat, ref.dat=ref.dat, k=k, method=method, dim=dim,return.distance=return.distance,index=index,transposed=TRUE,...)
        rm(dat)
        gc()
        knn
      }
    })    
    if(clear.index){
      cleanAnnoyIndex(index)
    }
    else{
      if (!return.distance) {
        result = list(knn.index = result)
      }
      result$index = index
    }
    return(result)
  }

map_cells_knn_big <- function(big.dat, cl.dat, select.cells, train.index=NULL, method = c("Annoy.Cosine","cor"), block.size=10000, mc.cores=10,clear.index=is.null(train.index))
  {
    cl.knn =  get_knn_batch_big(big.dat, cl.dat, select.cells=select.cells, k=1, index=train.index, method=method, block.size=block.size, mc.cores=mc.cores,return.distance=TRUE)
    knn.index = cl.knn[[1]]
    knn.dist = cl.knn[[2]]
    map.df = data.frame(sample_id=row.names(knn.index), cl = colnames(cl.dat)[knn.index], dist = knn.dist)
    return(map.df)
  }


big.dat.parquet.to.h5 <- function(big.dat, h5.file, type="float",mc.cores=10)
  {
    library(rhdf5)
    library(dplyr)
    library(arrow)
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    
    h5createFile(h5.file)
    big.dat = init_big.dat_parquet(big.dat)
    col.df= big.dat$col.df
    row.df = big.dat$row.df
    c.bin.size = col.df %>% pull(col_bin) %>% table %>% max
    r.bin.size = row.df %>% pull(row_bin) %>% table %>% max
    group="X"
    h5createDataset(h5.file,"X",c(nrow(row.df), nrow(col.df)),storage.mode=type,chunk=c(r.bin.size, c.bin.size))
    
    for(c.id in unique(big.dat$col.df$col_bin)) {
      cols = col.df %>% filter(col_bin==c.id) %>% pull(col_id)
      foreach(r.id = unique(big.dat$row.df$row_bin)) %dopar% {        
        rows = row.df %>% filter(row_bin==r.id) %>% pull(row_id)
        mat = as.matrix(get_cols(big.dat, cols=cols, rows=rows))
        h5write(mat, file=h5.file, name="X",index=list(rows, cols))
        cat("Finish", c.id, r.id, "\n")
      }
    }    
  }

display_cl_big <- function(big.dat, cl, ds, cl.bin, max.cl.size=200,markers=NULL,...)
  {
    tmp.cells= sample_cells(cl, max.cl.size)
    norm.dat = get_logNormal(big.dat, tmp.cells)
    if(is.null(markers)){
      markers = select_markers_ds(ds, cl.bin=cl.bin, select.cl = unique(cl))
    }
    display_cl(cl=cl[tmp.cells], norm.dat, markers=markers,...)
  }


cell_cl_cor <- function(big.dat, cl, cl.dat, cl.bin.size=50)
  {
    library(parallel)    
    require(doMC)
    require(foreach)
    foreach(a = unique(col.df$col_bin)) %dopar%{
      d= file.path(dir,a)
    }
  }
