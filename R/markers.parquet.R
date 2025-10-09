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
          tmp.pairs = to.add %>%  filter(bin.x == bin1 & bin.y ==bin2)
          if(nrow(tmp.pairs)==0){
            return(NULL)
          }
          tmp.de = ds %>% filter(bin.x == bin1 & bin.y ==bin2 & gene %in% genes & rank < max.num & P1 %in% tmp.pairs$P1 & P2 %in% tmp.pairs$P2) %>% collect
          tmp.de = tmp.de %>% right_join(tmp.pairs,by=c("P1","P2"))
          if(nrow(tmp.de)==0){
            return(NULL)
          }
          tmp.de = tmp.de %>% filter(!is.na(gene))
          if(nrow(tmp.de)==0){
            return(NULL)
          }
          cat(bin1, bin2, "\n")
          gene.score = tmp.de %>% group_by(gene) %>% summarize(rank.sum = sum(max.num- rank))
          list(gene.score)
        }
      gene.score=rbindlist(tmp)
      if(is.null(gene.score)){
        return(NULL)
      }
      if(nrow(gene.score)==0){
        return(NULL)
      }
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
    tmp = tmp%>% select(gene, new.score) %>% rename("new.score"="score") %>% filter(score > 0) %>% arrange(-score)
    return(tmp)
  }


check_pairs_lfc <- function(to.add, genes, cl.means, lfc.th=2, mc.cores=1)
  {

    require(doMC)
    registerDoMC(cores=min(mc.cores,length(genes)))
    to.add$P1 = as.character(to.add$P1)
    to.add$P2 = as.character(to.add$P2)
    cl.means=as.matrix(cl.means)

    g.bin.size = ceiling(length(genes)/mc.cores)
    gene.bin = data.frame(gene=genes, bin = ceiling((1:length(genes))/g.bin.size))
    de.checked=foreach::foreach(b=unique(gene.bin$bin), .combine="c")%dopar% {
      checked.sum=0
      v = gene.bin %>% filter(bin==b) %>% pull(gene)
      for(g in v){
        d1 = cl.means[g,to.add$P1]
        d2 = cl.means[g,to.add$P2]
        lfc = d1 -d2
        checked.sum=  checked.sum + as.integer(as.vector(lfc > lfc.th))
      }
      list(checked.sum)
    }
    checked.sum = Reduce("+",de.checked)    
    to.add$checked = checked.sum    
    return(to.add)
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
##' @title select_markers_pair_direction_ds
##' @param de.dir 
##' @param add.num 
##' @param genes 
##' @param cl.bin 
##' @param de 
##' @param mc.cores 
##' @param max.genes 
##' @param cl.means 
##' @param lfc.th 
##' @return value.  
##' @author Zizhen Yao
select_markers_pair_direction_ds <- function(de.dir, add.num, genes, cl.bin, de=NULL, mc.cores=20,max.genes=1000,cl.means=NULL,lfc.th=2, ds=NULL)
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
      if(is.null(gene.score)){
        break
      }
      if (nrow(gene.score)==0){
        break
      }
      g = as.character(gene.score$gene[1])
      if(is.na(g)){
        break
      }
      print(g)      
      if(is.null(cl.means)){
        new.checked = check_pairs_ds(de.dir, to.add=add.num %>% select(P1,P2), genes=g,cl.bin=cl.bin,de=de, ds=ds,mc.cores=mc.cores)
      }
      else{
        new.checked = check_pairs_lfc(to.add=add.num %>% select(P1,P2), genes=g,cl.means=cl.means, lfc.th=lfc.th,mc.cores=mc.cores)
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
      to.remove = add.num %>% filter(num <=0) %>% select(P1,P2)
      add.num = add.num %>% filter(num > 0)
      genes = setdiff(genes,g)
      select.genes=c(select.genes,g)
      if(nrow(add.num)==0){
        break
      }
      #gene.score = update_gene_score_ds(gene.score, ds=ds, to.remove=to.remove, cl.bin=cl.bin, de=de,mc.cores=mc.cores)
      #gene.score = gene.score %>% filter(gene %in% genes)
      if(nrow(to.remove)> nrow(add.num) |nrow(add.num)<1000){
        gene.score=get_gene_score_ds(ds, to.add=add.num, genes=genes,cl.bin=cl.bin, de=de, mc.cores=mc.cores)
      }
      else{
        rm.gene.score = get_gene_score_ds(ds, to.add=to.remove, genes=genes,cl.bin=cl.bin, de=de, mc.cores=mc.cores)
        if(is.null(rm.gene.score)){
          gene.score = gene.score %>% filter(gene !=g)
          next
        }
        rm.gene.score=rm.gene.score %>% rename(rm.score=score)
        tmp = gene.score %>% filter(!gene ==g) %>% left_join(rm.gene.score) %>% mutate(rm.score = ifelse(is.na(rm.score), 0, rm.score)) %>% mutate(new.score = score - rm.score) 
        tmp = tmp%>% select(gene, new.score) %>% rename(score=new.score) %>% filter(score > 0) %>% arrange(-score)
        gene.score=tmp
      }
      if(!is.null(de)){
        de = de %>% filter(gene!=g)
        if(is.null(de)| nrow(de)==0){
          break
        }
      }
    }
    return(list(select.genes=select.genes, de=de))
  }

#' Title: select_markers_pair_group_top_ds
#'
#' @param cl cl
#' @param g1 g1
#' @param g2 g2
#' @param ds ds
#' @param top.n top.n
#' @param max.num max.num
#' @param n.markers n.markers
#' @param up.gene.score up.gene.score
#' @param down.gene.score down.gene.score 
#'
#' @return value. Note # Always directory
#' @export

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
    cat("mc.cores=", mc.cores, "\n")
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
    if(is.null(markers)){
      return(markers)
    }
    add.num$num = n.markers
    if(is.null(cl.means)){
      new.checked = check_pairs_ds(de.dir, to.add=add.num %>% select(P1,P2),genes=markers,cl.bin=cl.bin,ds=ds,mc.cores=mc.cores)
    }
    else{
      new.checked = check_pairs_lfc(to.add=add.num %>% select(P1,P2),genes=markers,cl.means=cl.means, lfc.th=lfc.th,mc.cores=mc.cores)
    }
    add.num = add.num %>% left_join(new.checked, by=c("P1","P2"))      
    add.num = add.num %>% mutate(checked=ifelse(is.na(checked),0,checked)) %>% mutate(num = num-checked)            
    add.num = add.num %>% filter(num > 0)
    max.genes= max.genes - length(markers)
    if(nrow(add.num)>0 & max.genes > 0 & length(genes)>1){
      add.num$checked=NULL
      genes = setdiff(genes,markers)
      more.markers <- select_markers_pair_direction_ds(de.dir, add.num=add.num, genes=genes,cl.bin=cl.bin,max.genes=max.genes,cl.means=cl.means, lfc.th=lfc.th, mc.cores=mc.cores, ds=ds)
      more.markers <- more.markers$select.genes
      markers= c(markers, more.markers)
    }
    return(markers)
  }

select_N_markers_ds<- function(de.dir, select.cl=NULL,pair.num=1, add.num=NULL, genes, cl.bin, default.markers=NULL,cl.means=NULL,...)
  {
    if(is.null(add.num)){
      add.num = as.data.frame(create_pairs(select.cl, direction="directional"))
      add.num$num = pair.num
    }
    if(!is.null(default.markers)){
      if(is.null(cl.means)){
        de.checked.num = check_pairs_ds(de.dir, add.num, genes=default.markers, cl.bin=cl.bin,...)
      }
      else{
        de.checked.num = check_pairs_lfc(to.add=add.num %>% select(P1,P2), genes=default.markers,cl.means=cl.means,...)
      }
      add.num = add.num %>% left_join(de.checked.num,by=c("P1","P2"))
      add.num = add.num %>% mutate(checked=ifelse(is.na(checked),0,checked)) %>% mutate(num = num-checked)
      add.num = add.num %>% filter(num > 0)
      genes = setdiff(genes, default.markers)
      add.num$checked=NULL
    }
    print(dim(add.num))
    if(nrow(add.num)==0){
      return(NULL)
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
#' @return value. 
#' @export

select_pos_markers_ds<- function(de.dir, cl, select.cl, genes, cl.bin, ds=NULL,n.markers=1, out.dir="cl.markers",  mc.cores=1, overwrite=TRUE, ...)
  {
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    if(is.null(ds)){
      ds = open_dataset(de.dir)
    }
    if (!dir.exists(out.dir)) {
      dir.create(out.dir)
    }
    ###for each cluster, find markers that discriminate it from other types
    cl.markers <- foreach(x=select.cl, .combine="c") %dopar% {
      fn = file.path(out.dir,gsub("/", "", paste0(x, ".markers.rda")))
      if(file.exists(fn) & !overwrite){
        load(fn)
        tmp=list(markers)
        names(tmp)=x
        tmp
        return(tmp)
      }
      print(x)
      g1=x
      g2 = setdiff(cl, x)
      #select.de = ds %>% filter(bin.x==cl.bin.x & P1==g1 & P2 %in% g2) %>% collect()
      markers <- select_markers_pair_group_ds(g1,g2, de.dir=de.dir,  ds=ds, genes=genes, cl.bin=cl.bin, n.markers=n.markers,select.sign="up",...)
      save(markers, file=fn)
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
#' @return value. 
#' @export



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
select_markers_groups <- function(de.dir, cl.group, genes, cl.bin, select.groups= unique(cl.group$group), lfc.th=2,cl.means=NULL, n.markers=20,ds=NULL, mc.cores=1,...)
  {
    if(is.null(ds)){
      ds = open_dataset(de.dir)
    }
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
    
    if(is.null(cl.means)){
      de.checked.num = check_pairs_ds(de.dir, pairs[,c("P1","P2")], genes=group.markers,cl.bin=cl.bin, ds=ds, mc.cores=mc.cores)          
    }
    else{
      de.checked.num = check_pairs_lfc(to.add=pairs[,c("P1","P2")], genes=group.markers,cl.means=cl.means, lfc.th=lfc.th,mc.cores=mc.cores)
      }


    add.num = pairs %>% left_join(de.checked.num) %>% mutate(num=n.markers - checked)
    
    select.markers=group.markers
    genes = setdiff(genes, select.markers)
    more.markers <- select_markers_pair_direction_ds(de.dir, add.num=add.num,  genes=genes, cl.bin=cl.bin, ds=ds, mc.cores=mc.cores,...)    
    select.markers = c(select.markers, more.markers$markers)
  }

  
get_de_genes <- function(ds, cl.bin, cl1, cl2)
  {
    bin1 = cl.bin %>% filter(cl==cl1) %>% pull(bin)
    bin2 = cl.bin %>% filter(cl==cl2) %>% pull(bin)
    select.genes = ds %>% filter(bin.x==bin1 & P1==cl1 & bin.y==bin2 & P2==cl2 | bin.x==bin2 & P2==cl1 & bin.y==bin1 & P1==cl2) %>% collect()
    return(select.genes)
  }


select_markers_pair_group_lc <- function(g1,g2,cl.dat, n.markers=20, genes=genes,select.sign=c("up","down"), base=2)
  {
    g1 = intersect(g1, colnames(cl.dat))
    g2 = intersect(g2, colnames(cl.dat))
    gene.m1=rowMeans(cl.dat[,g1,drop=FALSE])
    gene.m2=rowMeans(cl.dat[,g2,drop=FALSE])
    select.markers=c()
    if("up" %in% select.sign){
      r = sort(gene.m1/(gene.m1 + gene.m2 + base),decreasing=TRUE)
      r = r[r > 0.6]
      select.markers=c(select.markers,head(names(r), n.markers))
    }
    if("down" %in% select.sign){
      r = sort(gene.m2/(gene.m1 + gene.m2 + base),decreasing=TRUE)
      r = r[r > 0.6]
      select.markers=c(select.markers,head(names(r), n.markers))
    }
    return(select.markers)
  }
