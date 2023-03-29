get_gene_score_ds <- function(ds, to.add, genes, cl.bin, de=NULL, max.num=1000,mc.cores=20)
  {
    require(doMC)
    registerDoMC(cores=mc.cores)
    if(!is.null(de)){
      tmp.de = de %>% right_join(to.add[,c("P1","P2")])
      gene.score = tmp.de %>% group_by(gene) %>% summarize(score = sum(as.numeric(max.num- rank))) %>% filter(gene %in% genes) %>% arrange(-score)   
    }
    else{
      to.add = to.add %>% left_join(cl.bin,by=c("P1"="cl")) %>% left_join(cl.bin,by=c("P2"="cl")) 
      cl.bin.x = to.add %>% pull(bin.x) %>% unique
      cl.bin.y = to.add %>% pull(bin.y) %>% unique  
      tmp=foreach::foreach(bin1=cl.bin.x,.combine="c")%:%
        foreach::foreach(bin2=cl.bin.y,.combine="c")%dopar% {
          tmp.pairs = to.add %>%  filter(bin.x == bin1 & bin.y ==bin.y)
          tmp.de = ds %>% filter(bin.x == bin1 & bin.y ==bin.y & gene %in% genes & rank < max.num & P1 %in% tmp.pairs$P1 & P2 %in% tmp.pairs$P2) %>% collect
          tmp.de = tmp.de %>% right_join(tmp.pairs,by=c("P1","P2"))
          gene.score = tmp.de %>% group_by(gene) %>% summarize(rank.sum = sum(max.num- rank))
          list(gene.score)
        }
      gene.score=rbindlist(tmp)
      gene.score = gene.score %>% group_by(gene) %>% summarize(score=sum(as.numeric(rank.sum))) %>% arrange(-score)
    }
    return(gene.score)
  }


update_gene_score_ds <- function(gene.score, ds, to.remove, cl.bin, de=NULL, max.num=1000,mc.cores=20)
  {
    require(doMC)
    registerDoMC(cores=mc.cores)
    to.remove = to.remove %>% left_join(cl.bin,by=c("P1"="cl")) %>% left_join(cl.bin,by=c("P2"="cl"))
    cl.x = to.remove %>% pull(P1) %>%unique
    cl.y = to.remove %>% pull(P2) %>%unique
    
    cl.bin.x = to.remove %>% pull(bin.x) %>% unique
    cl.bin.y = to.remove %>% pull(bin.y) %>% unique
    if(!is.null(de)){
      tmp.de = de %>% right_join(to.remove)
      rm.gene.score = tmp.de %>% group_by(gene) %>% summarize(rm.score = sum(max.num- rank))
    }
    else{
      if(length(cl.bin.x)*length(cl.bin.y) * nrow(gene.score) < 10^6){
        de=  ds %>% filter(bin.x %in% cl.bin.x & bin.y %in% cl.bin.y & gene %in% gene.score$gene & rank < max.num & P1 %in% cl.x & P2 %in% cl.y) %>% collect
        tmp.de = de %>% right_join(to.remove)
        rm.gene.score = tmp.de %>% group_by(gene) %>% summarize(rm.score = sum(max.num- rank))        
      }
      else{
        tmp=foreach::foreach(bin1=cl.bin.x,.combine="c")%:%
          foreach::foreach(bin2=cl.bin.y,.combine="c")%dopar% {
            tmp.pairs = to.remove %>%  filter(bin.x == bin1 & bin.y ==bin2)
            tmp.de = ds %>% filter(bin.x == bin1 & bin.y ==bin2 & gene %in% gene.score$gene & rank < max.num & P1 %in% tmp.pairs$P1 & P2 %in% tmp.pairs$P2) %>% collect
            tmp.de = tmp.de %>% right_join(tmp.pairs)
            gene.score = tmp.de %>% group_by(gene) %>% summarize(rank.sum = sum(max.num- rank))
            list(gene.score)
          }
        rm.gene.score=rbindlist(tmp)      
        if(is.null(rm.gene.score)|nrow(rm.gene.score)==0){
          return(gene.score)
        }
        rm.gene.score = rm.gene.score %>% group_by(gene) %>% summarize(rm.score=sum(as.numeric(rank.sum)))
      }
    }
    tmp = gene.score %>% left_join(rm.gene.score) %>% mutate(rm.score = ifelse(is.na(rm.score), 0, rm.score)) %>% mutate(new.score = score - rm.score) 
    tmp = tmp%>% select(gene, new.score) %>% rename(score=new.score) %>% filter(score > 0) %>% arrange(-score)
    return(tmp)
  }


check_pairs_lfc <- function(to.add, genes, cl.means, lfc.th=2, mc.cores=1)
  {
    require(doMC)
    registerDoMC(cores=min(mc.cores,length(genes)))
    to.add$P1 = as.character(to.add$P1)
    to.add$P2 = as.character(to.add$P2)
    cl.means=as.matrix(cl.means)
    de.checked=foreach::foreach(g=genes, .combine="cbind")%dopar% {
      d1=cl.means[g,to.add$P1]
      d2 = cl.means[g, to.add$P2]
      lfc = d1 -d2
      lfc > lfc.th
    }
    de.checked = as.matrix(de.checked, nrow=nrow(to.add))
    to.add$checked = rowSums(de.checked)
    de.checked = to.add %>% filter(checked>0)
    return(de.checked)
  }

check_pairs_ds <- function(de.dir, to.add, genes,cl.bin, de=NULL, mc.cores=10,max.num=1000, ds=NULL)
  {
    require(doMC)
    registerDoMC(cores=mc.cores)
    to.add = to.add[,c("P1","P2")]
    to.add = to.add %>% left_join(cl.bin,by=c("P1"="cl")) %>% left_join(cl.bin,by=c("P2"="cl"))
    select.P1 = to.add %>% pull(P1) %>% unique
    select.P2 = to.add %>% pull(P2) %>% unique    
    if(!is.null(de)){
      de.checked = to.add %>% left_join(de %>% filter(rank < max.num & gene %in% genes)) %>% group_by(P1,P2) %>% summarize(checked=sum(!is.na(gene)))
    }
    else{
      cl.bin.x = to.add %>% pull(bin.x) %>% unique
      cl.bin.y = to.add %>% pull(bin.y) %>% unique
      if(length(cl.bin.x)*length(cl.bin.y)*length(genes) < 10^5){
        if(is.null(ds)){
          ds = open_dataset(de.dir)
        }
        de = ds %>% filter(bin.x  %in% cl.bin.x & bin.y %in% cl.bin.y & gene %in% genes & P1 %in% select.P1 & P2 %in% select.P2 & rank < max.num ) %>% collect
        de.checked = to.add %>% left_join(de %>% filter(rank < max.num & gene %in% genes)) %>% group_by(P1,P2) %>% summarize(checked=sum(!is.na(gene)))
      }
      else{
        de.checked=foreach::foreach(bin1=cl.bin.x,.combine="c")%:%
          foreach::foreach(bin2=cl.bin.y,.combine="c")%dopar% {
            cat("bin", bin1, bin2, "\n")
            tmp.pairs = to.add %>%  filter(bin.x == bin1 & bin.y ==bin2)
            ###
            #tmp.de =ds %>% filter(bin.x == bin1 & bin.y ==bin.y) %>% collect %>% filter(gene %in% genes & rank < max.num)
            #Directly read from parquet file is faster
            d = file.path(de.dir, paste0("bin.x=",bin1), paste0("bin.y=",bin2))
            fn = dir(d, pattern="parquet")
            tmp.de = read_parquet(file.path(d, fn))
            tmp.de = tmp.de %>% filter(gene %in% genes & rank < max.num)
            tmp.de = tmp.de %>% right_join(tmp.pairs)
            de.checked = tmp.de %>% group_by(P1,P2) %>% summarize(checked=n())
            list(de.checked)
          }      
        de.checked = rbindlist(de.checked)
      }
    }
    return(de.checked)
  }




select_markers_ds <- function(ds, cl.bin, select.cl=NULL, top.n=20,mc.cores=10)
  {
    if(!is.null(select.cl)){
      cl.bin = cl.bin %>% filter(cl %in% select.cl)
    }
    select.bin = cl.bin %>% pull(bin) %>% unique
    mc.cores=min(mc.cores, length(select.bin))
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)    
    tmp=foreach::foreach(bin1=select.bin,.combine="c")%:%
      foreach::foreach(bin2=select.bin,.combine="c")%dopar% {
        de = ds %>% filter(bin.x %in% bin1 & bin.y %in% bin2)
        if(!is.null(select.cl)){
          de = de %>% filter(P1 %in% select.cl & P2 %in% select.cl)            
        }
        de %>% filter(rank <= top.n) %>% pull(gene) %>% unique
      }
    select.markers=unique(tmp)
    return(select.markers)
  }


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param de.dir 
##' @param add.num 
##' @param genes 
##' @param cl.bin 
##' @param de 
##' @param mc.cores 
##' @param max.genes 
##' @param cl.means 
##' @param lfc.th 
##' @return 
##' @author Zizhen Yao
select_markers_pair_direction_ds <- function(de.dir, add.num, genes, cl.bin, de=NULL, mc.cores=mc.cores,max.genes=1000,cl.means=NULL,lfc.th=2, ds=NULL)
  {
    if(is.null(ds)){
      ds = open_dataset(de.dir)
    }
    select.genes=c()
    if(!is.null(de)){
      de = de %>% filter(gene %in% genes)
    }
    add.num = add.num %>% left_join(cl.bin,by=c("P1"="cl")) %>% left_join(cl.bin,by=c("P2"="cl"))
    gene.score= get_gene_score_ds(ds, to.add=add.num, genes=genes,cl.bin=cl.bin, de=de, mc.cores=mc.cores)
    while(nrow(add.num)>0 & length(genes)>0 & length(select.genes)< max.genes){      
      if(is.null(gene.score) | nrow(gene.score)==0){
        break
      }
      g = as.character(gene.score$gene[1])
      print(g)      
      if(is.null(cl.means)){
        new.checked = check_pairs_ds(de.dir, to.add=add.num %>% select(P1,P2), genes=g,cl.bin=cl.bin,de=de, mc.cores=mc.cores)
      }
      else{
        new.checked = check_pairs_lfc(to.add=add.num %>% select(P1,P2), genes=g,cl.means=cl.means, lfc.th=lfc.th)
      }
      if(is.null(new.checked)){
        break
      }
      if(nrow(new.checked)==0){
        break
      }
      add.num$checked=NULL
      add.num = add.num %>% left_join(new.checked, by=c("P1","P2"))      
      add.num = add.num %>% mutate(checked=ifelse(is.na(checked),0,checked)) %>% mutate(num = num-checked)
      if(nrow(add.num)==0){
        break
      }
      to.remove = add.num %>% filter(num <=0) %>% select(P1,P2)
      add.num = add.num %>% filter(num > 0)
      genes = setdiff(genes,g)
      select.genes=c(select.genes,g)
      gene.score = update_gene_score_ds(gene.score, ds=ds, to.remove=to.remove, cl.bin=cl.bin, de=de,mc.cores=mc.cores)
      gene.score = gene.score %>% filter(gene %in% genes)
      if(!is.null(de)){
        de = de %>% filter(gene!=g)
        if(is.null(de)| nrow(de)==0){
          break
        }
      }
    }
    return(list(select.genes=select.genes, de=de))
  }

#' Title
#'
#' @param cl 
#' @param g1 
#' @param g2 
#' @param ds
#' @param top.n 
#' @param max.num 
#' @param n.markers 
#' @param up.gene.score 
#' @param down.gene.score 
#'
#' @return
#' @export
#'
#' @examples
# Always directory
select_markers_pair_group_top_ds <- function(g1,g2,ds, genes, cl.bin, select.sign="up", n.markers=20,mc.cores=1, ...)
{
  require(matrixStats)
  require(data.table)
  require(arrow)
  require(dplyr)

  require(doMC)
  registerDoMC(cores=mc.cores)
  up.to.add = down.to.add=NULL
  up.genes=down.genes=NULL
  if("up" %in% select.sign){
    up.to.add = as.data.frame(create_pairs(g1, g2, direction="unidirectional"))
    up.gene.score  = get_gene_score_ds(ds, to.add=up.to.add, genes=genes, cl.bin=cl.bin,...)
    up.genes = head(up.gene.score$gene, n.markers)
  }
  if("down" %in% select.sign){
    down.to.add = as.data.frame(create_pairs(g2, g1,direction="unidirectional"))
    down.gene.score  = get_gene_score_ds(ds, to.add=down.to.add, genes=genes, cl.bin=cl.bin, ...)
    down.genes = head(down.gene.score$gene, n.markers)
  }
  to.add=rbind(up.to.add, down.to.add)
  return(list(up.genes=up.genes, down.genes=down.genes,to.add=to.add))
}

#############
select_markers_pair_group_ds <- function(g1,g2,de.dir, genes, cl.bin, n.markers=20,select.sign=c("up","down"),max.genes=50,default.markers=NULL, cl.means=NULL,lfc.th=2, mc.cores=1, ds=NULL)
  {
    if(is.null(ds)){
      ds = open_dataset(de.dir)
    }
    if(is.null(default.markers)){
      result = select_markers_pair_group_top_ds(g1,g2,ds=ds, genes=genes, cl.bin=cl.bin, n.markers=n.markers,select.sign=select.sign,mc.cores=mc.cores)
      markers=c(result$up.genes, result$down.genes)
      add.num = result$to.add
    }
    else{
      add.num=NULL
      if("up" %in% select.sign){
        add.num = as.data.frame(create_pairs(g1, g2, direction="unidirectional"))
      }
      if("down" %in% select.sign){
        add.num = rbind(add.num,as.data.frame(create_pairs(g2, g1, direction="unidirectional")))
      }
      markers = default.markers
    }
    add.num$num = n.markers
    if(is.null(cl.means)){
      new.checked = check_pairs_ds(de.dir, to.add=add.num %>% select(P1,P2),genes=markers,cl.bin=cl.bin,mc.cores=mc.cores)
    }
    else{
      new.checked = check_pairs_lfc(to.add=add.num %>% select(P1,P2),genes=markers,cl.means=cl.means, lfc.th=lfc.th,mc.cores=mc.cores)
    }
    add.num = add.num %>% left_join(new.checked, by=c("P1","P2"))      
    add.num = add.num %>% mutate(checked=ifelse(is.na(checked),0,checked)) %>% mutate(num = num-checked)            
    add.num = add.num %>% filter(num > 0)
    max.genes= max.genes - length(markers)
    if(nrow(add.num)>0 & max.genes > 0){
      add.num$checked=NULL
      genes = setdiff(genes,markers)
      more.markers <- select_markers_pair_direction_ds(de.dir, add.num=add.num, genes=genes,cl.bin=cl.bin,max.genes=max.genes,cl.means=cl.means, lfc.th=lfc.th, mc.cores=mc.cores, ds=ds)
      more.markers <- more.markers$select.genes
      markers= c(markers, more.markers)
    }
    return(markers)
  }

select_N_markers_ds<- function(de.dir, select.cl=NULL,pair.num=1, add.num=NULL, genes, cl.bin, default.markers=NULL,...)
  {
    if(is.null(add.num)){
      add.num = as.data.frame(create_pairs(select.cl, direction="directional"))
      add.num$num = pair.num
    }
    if(!is.null(default.markers)){      
      de.checked.num = check_pairs_ds(de.dir, add.num, genes=default.markers, cl.bin=cl.bin,...)
      add.num = add.num %>% left_join(de.checked.num,by=c("P1","P2"))
      add.num = add.num %>% mutate(checked=ifelse(is.na(checked),0,checked)) %>% mutate(num = num-checked)
      add.add = add.num %>% filter(num > 0)
      genes = setdiff(genes, default.markers)
      add.num$checked=NULL
    }
    markers <- select_markers_pair_direction_ds(de.dir, add.num=add.num, genes=genes,cl.bin=cl.bin,...)$select.genes    
  }



#' Title
#'
#' @param de.genes 
#' @param cl 
#' @param n.markers 
#' @param default.markers 
#' @param rm.genes 
#' @param up.gene.score 
#' @param down.gene.score 
#'
#' @return
#' @export
#'
#' @examples
select_pos_markers_ds<- function(de.dir, cl, select.cl, genes, cl.bin, n.markers=1,  mc.cores=1,out.dir="cl.markers",...)
  {
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    ds = open_dataset(de.dir)
    if (!dir.exists(out.dir)) {
        dir.create(out.dir)
    }
    
    ###for each cluster, find markers that discriminate it from other types
    cl.markers <- foreach(x=select.cl, .combine="c") %dopar% {
      print(x)
      g1=x
      g2 = setdiff(cl, x)
      #select.de = ds %>% filter(bin.x==cl.bin.x & P1==g1 & P2 %in% g2) %>% collect()
      markers <- select_markers_pair_group_ds(g1,g2, de.dir=de.dir,  genes=genes, cl.bin=cl.bin, n.markers=n.markers,select.sign="up",...)
      save(markers, file=file.path(out.dir, paste0(x, ".markers.rda")))
      tmp=list(markers)
      names(tmp)=x
      tmp
    }
    return(cl.markers)
  }

#marker.result <- select_markers_pair_group_ds(g1,g2, ds=ds, de=select.de, genes=genes, cl.bin=cl.bin, n.markers=n.markers,select.sign="up")

select_top_pos_markers_ds<- function(ds, cl, select.cl, genes, cl.bin, n.markers=3, mc.cores=10,...)
  {
    library(parallel)
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    
    ###for each cluster, find markers that discriminate it from other types

    cl.markers <- foreach(x=select.cl, .combine="c") %dopar% {
      print(x)
      g1=x
      g2 = setdiff(cl, x)
      cl.bin.x = cl.bin %>% filter(cl==g1) %>% pull(bin)      
      select.de = ds %>% filter(bin.x==cl.bin.x & P1==g1 & P2 %in% g2) %>% collect()
      markers= select_markers_pair_group_top_ds(g1,g2,ds, genes=genes,cl.bin=cl.bin, de=select.de,n.markers=n.markers,select.sign="up",...)$up.genes
      tmp=list(markers)
      names(tmp)=x
      tmp
    }
    return(cl.markers)
  }



###select_markers_groups 
#' Title
#'
#' @param de.genes 
#' @param cl.group Assignment of clusters to groups cluster as names, and group id as values.
#'
#' @return
#' @export
#'
#' @examples


select_markers_groups_top_ds <- function(ds, cl.group, select.groups=names(cl.group), n.markers=3,mc.cores=1,...)
  {

    library(parallel)
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    
    all.cl = unlist(cl.group)
    group.markers <- foreach(x=select.groups, .combine="c") %dopar% {
      print(x)
      g1 = cl.group[[x]]
      g2 = setdiff(all.cl, g1)              
      markers=select_markers_pair_group_top_ds(g1,g2,ds=ds, select.sign="up",n.markers=n.markers, ...)$up.genes
      list(markers)
    }
    names(group.markers) = select.groups
    return(group.markers)
  }

#cl.group is a data.frame with column "cl", and "group"
select_markers_groups <- function(de.dir, cl.group, genes, cl.bin, select.groups= unique(cl.group$group), n.markers=20,mc.cores=1,...)
  {
    ds = open_dataset(de.dir)
    cl.group$cl = as.character(cl.group$cl)
    group_pair=create_pairs(unique(cl.group$group))
    library(parallel)
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)

    group.markers <- foreach(i=1:nrow(group_pair), .combine="c") %dopar% {
      x= group_pair[i,1]
      y= group_pair[i,2]
      cat(x,y,"\n")
      g1 = cl.group %>% filter(group==x) %>% pull(cl)
      g2 = cl.group %>% filter(group==y) %>% pull(cl)
      result=select_markers_pair_group_top_ds(g1,g2,ds=ds, n.markers=n.markers, genes=genes,cl.bin=cl.bin,select.sign=c("up","down"),...)
      list(c(result$up.genes, result$down.genes))
    }
    group.markers=unique(unlist(group.markers))   
    pairs = as.data.frame(create_pairs(cl.group$cl,direction="directional"))
    pairs$pair = row.names(pairs)
    pairs = pairs %>% left_join(cl.group, by=c("P1"="cl"))
    pairs = pairs %>% left_join(cl.group, by=c("P2"="cl"))
    pairs = pairs %>% filter(group.x!=group.y) 
   
    de.checked.num = check_pairs_ds(de.dir, pairs[,c("P1","P2")], genes=group.markers,cl.bin=cl.bin, mc.cores=mc.cores)    
    add.numm = pairs %>% left_join(de.checked.num) %>% mutate(num=n.markers - checked)
    
    select.markers=group.markers
    genes = setdiff(genes, select.markers)
    more.markers <- select_markers_pair_direction_ds(de.dir, add.num=to.add,  genes=genes, cl.bin=cl.bin, mc.cores=mc.cores,...)    
    select.markers = c(select.markers, more.markers$markers)
  }

  
get_de_genes <- function(ds, cl.bin, cl1, cl2)
  {
    bin1 = cl.bin %>% filter(cl==cl1) %>% pull(bin)
    bin2 = cl.bin %>% filter(cl==cl2) %>% pull(bin)
    select.genes = ds %>% filter(bin.x==bin1 & P1==cl1 & bin.y==bin2 & P2==cl2 | bin.x==bin2 & P2==cl1 & bin.y==bin1 & P1==cl2) %>% collect()
    return(select.genes)
  }

score_group_markers_fast <- function(cl.dat, g1, g2, diff.th=2)
  {
    g1 = intersect(g1, colnames(cl.dat))
    g2 = intersect(g2, colnames(cl.dat))
    gene.m1=rowMeans(cl.dat[,g1,drop=FALSE])
    gene.m2=rowMeans(cl.dat[,g2,drop=FALSE])
    th  = gene.m2 + diff.th
    fg = rowMeans(cl.dat[,g1,drop=FALSE] > th)
    bg = rowMeans(cl.dat[,g2,drop=FALSE] > diff.th)
    r = sort(fg/(bg+0.0005), decreasing=TRUE)
    return(r)
  }
