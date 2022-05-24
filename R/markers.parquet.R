get_gene_score <- function(ds, to.add, genes, de=NULL, max.num=1000,all.pairs=NULL,mc.cores=20)
  {
    require(doMC)
    registerDoMC(cores=mc.cores)
    up.pair = to.add %>% filter(sign=="up") %>% pull(pair)
    down.pair = to.add %>% filter(sign=="down") %>% pull(pair)
    if(is.null(all.pairs)){
      all.bins = ds %>% pull(pair_bin) %>% unique()
      tmp <- foreach(x=all.bins, .combine="c") %dopar% {
        pair = ds %>% filter(pair_bin ==x) %>% select(pair_bin, pair) %>% collect() %>% distinct()
        list(pair)
      }
      all.pairs=rbindlist(tmp)
    }
    select.pair.bin = all.pairs %>% filter(pair %in% c(up.pair, down.pair)) %>% pull(pair_bin) %>% unique
    bin = max(ceiling(length(select.pair.bin)/(2*mc.cores)),2)
    bins = split(select.pair.bin, ceiling((1:length(select.pair.bin))/bin))
    gene.score <- foreach(i=1:length(bins), .combine="c") %dopar% {
      x=bins[[i]]
      if(is.null(de)){
        tmp.de = ds %>% filter(pair_bin %in% x) %>%
          filter(gene %in% genes & rank < max.num & ((sign=="up" & pair %in% up.pair) | (sign=="down" & pair %in% down.pair))) %>% collect()
      }
      else{
        tmp.de = de %>% filter(pair_bin %in% x)%>%
          filter(gene %in% genes & rank < max.num) %>% 
          filter(((sign=="up" & pair %in% up.pair) | (sign=="down" & pair %in% down.pair)))
      }
      
      gene.score = tmp.de %>% group_by(gene) %>% summarize(rank.sum = sum(max.num- rank))
      list(gene.score)
    }
    gene.score=rbindlist(gene.score)
    gene.score = gene.score %>% group_by(gene) %>% summarize(score=sum(as.numeric(rank.sum))) %>% arrange(-score)
    return(gene.score)
  }

check_pairs <- function(ds, to.add, genes,de=NULL, all.pairs=NULL, mc.cores=10)
  {
    require(doMC)
    registerDoMC(cores=mc.cores)
    up.pair = to.add %>% filter(sign=="up") %>% pull(pair)
    down.pair = to.add %>% filter(sign=="down") %>% pull(pair)
    if(is.null(all.pairs)){
      all.bins = ds %>% pull(pair_bin) %>% unique
      tmp <- foreach(x=all.bins, .combine="c") %dopar% {        
        pair = ds %>% filter(pair_bin ==x) %>% select(pair_bin, pair) %>% collect() %>% distinct()
        list(pair)
      }
      all.pairs=rbindlist(tmp)
    }
    select.pair.bin = all.pairs %>% filter(pair %in% c(up.pair, down.pair)) %>% pull(pair_bin) %>% unique
    bin = max(ceiling(length(select.pair.bin)/(2*mc.cores)),2)
    bins = split(select.pair.bin, ceiling((1:length(select.pair.bin))/bin))
    de.checked <- foreach(i=1:length(bins), .combine="c") %dopar% {
      x = bins[[i]]
      if(is.null(de)){
        de.checked = ds %>% filter(pair_bin %in% x & gene %in% genes) %>% 
          filter(((sign=="up" & pair %in% up.pair) | (sign=="down" & pair %in% down.pair))) %>%
            group_by(pair,sign) %>% collect()
      }
      else{
        de.checked = de %>% filter(pair_bin %in% x & gene %in% genes) %>% 
          filter(((sign=="up" & pair %in% up.pair) | (sign=="down" & pair %in% down.pair))) %>%
            group_by(pair,sign)
      }
      de.checked = de.checked %>%summarize(checked = n())      
      list(de.checked)
    }
    de.checked = rbindlist(de.checked)        
  }

select_markers <- function(ds, pairs = NULL, top.n=20)
  {
    if(is.null(pairs)){
      select.markers = ds %>% filter(rank < top.n) %>% pull(gene)  %>% unique
    }
    else{
      select.markers = ds %>% filter(pair %in% pairs & rank < top.n) %>% pull(gene) %>% unique
    }
  }


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param ds 
##' @param add.num 
##' @param genes.allowed 
##' @param de 
##' @param max.num 
##' @param cl.means 
##' @param all.pairs 
##' @param mc.cores 
##' @return 
##' @author Zizhen Yao
select_markers_pair_direction <- function(ds, add.num, genes.allowed, de=NULL, max.num=1000, cl.means=NULL, all.pairs=NULL,mc.cores=10)
  {
    select.genes=c()    
    while(nrow(add.num)>0 & length(genes.allowed)>0){
      gene.score= get_gene_score(ds, to.add=add.num, genes=genes.allowed,de=de, max.num=max.num, all.pairs=all.pairs,mc.cores=mc.cores)
      if(is.null(gene.score) | nrow(gene.score)==0){
        break
      }
      print(head(gene.score))
      g = as.character(gene.score$gene[1])
      print(g)
      new.checked = check_pairs(ds, add.num,de=de, genes=g,all.pairs=all.pairs,mc.cores=mc.cores)
      add.num$checked=NULL
      add.num = add.num %>% left_join(new.checked, by=c("sign","pair"))
      add.num = add.num %>% mutate(checked=ifelse(is.na(checked),0,checked)) %>% mutate(num = num-checked)      
      to.remove.up = add.num %>% filter(num <= 0 & sign=="up") %>% pull(pair)
      to.remove.down = add.num %>% filter(num <= 0 & sign=="down") %>% pull(pair)
      add.num = add.num %>% filter(num > 0)    
      genes.allowed = setdiff(genes.allowed,g)
      select.genes=c(select.genes,g)
      if(!is.null(de)){
        de = de %>% filter(gene!=g & ((sign=="up" & !pair %in% to.remove.up)| (sign=="down"& !pair %in% to.remove.down)))
      }
    }
    return(list(markers=select.genes, de=de))
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
select_markers_pair_group_top<- function(cl, g1,g2,ds, genes, de= NULL, select.sign=c("up","down"),max.num=1000,all.pairs=NULL, n.markers=20)
{
  require(matrixStats)
  require(data.table)
  require(arrow)
  require(dplyr)
  pairs.df = as.data.frame(create_pairs(g1,g2))
  pairs.df$pair =row.names(pairs.df)
  to.add = pairs.df
  to.add$sign = "up"
  to.add = rbindlist(list(to.add, to.add %>% mutate(sign="down")))
  flip.sign = c("up"="down","down"="up")
  to.add = to.add %>% mutate(group.sign = ifelse(P1 %in% g1,sign, flip.sign[sign]))
  to.add = to.add %>% filter(group.sign %in% select.sign)
  if(is.null(de)){
    if(!is.null(all.pairs)){
      select.pair.bin = all.pairs %>% filter(pair %in% to.add$pair) %>% pull(pair_bin) %>% unique
      ds = ds %>% filter(pair_bin %in% select.pair.bin)
    }
    de = ds %>% filter(pair %in% to.add$pair & gene %in% genes & rank < max.num) %>% collect() %>% left_join(to.add, by=c("pair","sign")) %>% filter(!is.na(group.sign))
  } 
  else{
    de = de %>% filter(pair %in% to.add$pair & gene %in% genes & rank < max.num) %>% collect() %>% left_join(to.add, by=c("pair","sign")) %>% filter(!is.na(group.sign))
  }
  gene.score = de %>% group_by(group.sign, gene) %>% summarize(score = sum(max.num - rank))  %>% arrange(-score)
  up.genes = head(gene.score %>% filter(group.sign=="up") %>% pull(gene), n.markers)
  down.genes = head(gene.score %>% filter(group.sign=="down") %>% pull(gene), n.markers)
  return(list(up.genes=up.genes, down.genes=down.genes,de=de,to.add=to.add))
}

select_markers_pair_group<- function(cl, g1,g2,ds, de=NULL, genes, select.sign=c("up", "down"), max.num=1000,n.markers=20,all.pairs=NULL)
  {
    tmp = select_markers_pair_group_top(cl, g1,g2,ds=ds,de=de, genes=genes, select.sign=select.sign, max.num=max.num, n.markers=n.markers,all.pairs=all.pairs)
    de=tmp$de
    default.markers=c(tmp$up.genes,tmp$down.genes)
    add.num = tmp$to.add
    add.num$num = n.markers
    result = select_N_markers(ds, pairs=NULL, genes=genes, add.num=add.num, de= de, default.markers=default.markers,all.pairs=all.pairs)
    result$markers= c(default.markers, result$markers)
    return(result)
  }

select_N_markers <- function(ds, pairs, genes, pair.num=1, add.num=NULL, de=NULL, default.markers=NULL,all.pairs=NULL,mc.cores=10)
  {
    if(is.null(add.num)){
      add.up=data.frame(pair=pairs, sign="up", num=pair.num)
      add.down= data.frame(pair=pairs, sign="down", num=pair.num)
      add.num = rbind(add.up, add.down)
    }
    to.add = add.num
    if(!is.null(all.pairs)){
      select.pair.bin = all.pairs %>% filter(pair %in% to.add$pair) %>% pull(pair_bin) %>% unique
      ds = ds %>% filter(pair_bin %in% select.pair.bin)
    }    
    if(!is.null(default.markers)){      
      de.checked.num = check_pairs(ds, add.num, genes=default.markers, all.pairs=all.pairs, mc.cores=mc.cores,de=de)
      add.num = add.num %>% left_join(de.checked.num,by=c("pair","sign"))
      add.num = add.num %>% mutate(checked=ifelse(is.na(checked),0,checked)) %>% mutate(num = num-checked)
      to.add = add.num %>% filter(num > 0)
      genes = setdiff(genes, default.markers)      
    }
    to.add = to.add %>% select(pair, sign, num)
    if(!is.null(de)){
      de = de %>% filter(pair %in% to.add$pair & gene %in% genes)
    }
    marker.select.result <- select_markers_pair_direction(ds, to.add, de=de, genes.allowed=genes,all.pairs=all.pairs,mc.cores=mc.cores)
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
select_pos_markers <- function(ds,  cl, select.cl, genes, de=NULL, n.markers=1,  mc.cores=1, max.num=1000, all.pairs=NULL)
  {
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)

    ###for each cluster, find markers that discriminate it from other types

    cl.markers <- foreach(x=select.cl, .combine="c") %dopar% {
    #cl.markers <- list()
    #for(x in select.cl){
      print(x)
      g1=x
      g2 = setdiff(cl, x)
      select.pairs = as.data.frame(create_pairs(g1,g2))
      select.pairs$pair = row.names(select.pairs)
      if(is.null(all.pairs)){
        select.pairs = select.pairs %>% left_join(all.pairs)
      }
      if(is.null(de)){
        if(is.null(all.pairs)){
          ds = ds %>% filter(pair %in% select.pairs$pair_bin)
        }
        select.de = ds %>% filter(pair %in% select.pairs$pair) %>% collect()
      }
      else{
        select.de = de %>% filter(pair %in% select.pairs$pair) %>% collect()
      }
      marker.result <- select_markers_pair_group(cl, g1,g2, ds=ds, de=select.de, genes=genes, select.sign="up", max.num=max.num, all.pairs=all.pairs, n.marker=n.markers)
      list(marker.result$markers)
    }
    return(cl.markers)
  }


select_top_pos_markers <- function(ds, cl, select.cl, genes, n.markers=3, mc.cores=10, max.num=1000,all.pairs=NULL)
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
      markers= select_markers_pair_group_top(cl, g1,g2,ds, genes=genes, select.sign="up",max.num=max.num,n.markers=n.markers,all.pairs=NULL)$up.genes
      list(markers)
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


select_markers_groups_top <- function(ds, cl.group, select.groups=names(cl.group), n.markers=3,mc.cores=1,...)
  {

    library(parallel)
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    
    ###for each cluster, find markers that discriminate it from other types
    all.cl = unlist(cl.group)
    group.markers <- foreach(x=select.groups, .combine="c") %dopar% {
    #group.markers=list()
    #for(x in select.groups){
      print(x)
      g1 = cl.group[[x]]
      g2 = setdiff(all.cl, g1)              
      markers=select_markers_pair_group_top(cl=all.cl, g1,g2,ds=ds, select.sign="up",n.markers=n.markers, ...)$up.genes
      list(markers)
    }
    names(group.markers) = select.groups
    return(group.markers)
  }

#cl.group is a data.frame with column "cl", and "group"
select_markers_groups <- function(ds, cl.group, genes, select.groups= unique(cl.group$group), n.markers=20,mc.cores=1,all.pairs=NULL)
  {

    library(parallel)
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    cl.group$cl = as.character(cl.group$cl)
    group_pair=create_pairs(unique(cl.group$group))    
    
    group.markers <- foreach(i=1:nrow(group_pair), .combine="c") %dopar% {
      x= group_pair[i,1]
      y= group_pair[i,2]
      cat(x,y,"\n")
      g1 = cl.group %>% filter(group==x) %>% pull(cl)
      g2 = cl.group %>% filter(group==y) %>% pull(cl)
      result=select_markers_pair_group_top(cl=all.cl, g1,g2,ds=ds, n.markers=n.markers, genes=genes,all.pairs=all.pairs)
      list(c(result$up.genes, result$down.genes))
    }
    group.markers=unique(unlist(group.markers))
   
    pairs = as.data.frame(create_pairs(cl.group$cl))
    pairs$pair = row.names(pairs)
    pairs = pairs %>% left_join(cl.group, by=c("P1"="cl"))
    pairs = pairs %>% left_join(cl.group, by=c("P2"="cl"))
    pairs = pairs %>% filter(group.x!=group.y)
    
    tmp.pairs = pairs[,1:3]
    to.add = rbindlist(list(tmp.pairs %>% mutate(sign="down"),tmp.pairs %>% mutate(sign="up"))) %>% mutate(num=n.markers)
    de.checked.num = check_pairs(ds, to.add, genes=group.markers, all.pairs=all.pairs, mc.cores=20)    
    to.add = to.add %>% left_join(de.checked.num,by=c("pair","sign"))
    to.add = to.add %>% mutate(checked=ifelse(is.na(checked),0,checked)) %>% mutate(num = num-checked) %>% filter(num > 0)
    
    tmp.add = to.add %>% filter(num >= n.markers/2)
    more.markers = ds %>% filter(pair %in% tmp.add$pair & rank < n.markers) %>% collect() %>% right_join(tmp.add) %>% pull(gene) %>% unique    
    more.markers= setdiff(more.markers, group.markers)
    de.checked.num = check_pairs(ds, to.add, genes=more.markers, all.pairs=all.pairs, mc.cores=20)    
    to.add = to.add[,1:5] %>% left_join(de.checked.num,by=c("pair","sign"))
    to.add = to.add %>% mutate(checked=ifelse(is.na(checked),0,checked)) %>% mutate(num = num-checked) %>% filter(num > 0)
    select.markers=union(more.markers, group.markers)
    de = ds %>% filter(pair %in% to.add$pair & !gene %in% select.markers) %>% collect()    
    more.markers2 <- select_markers_pair_direction(ds, to.add, de=de, genes.allowed=unique(de$gene) ,all.pairs=all.pairs)    
    markers = head(c(select.markers, more.markers2$markers), 3000)    
  }

  
