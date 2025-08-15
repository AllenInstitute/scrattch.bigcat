#' Title
#'
#' @param cl.rd 
#'
#' @return
#' @export
#'
#' @examples
get_cl_sim <- function(cl.rd)
{
  if(ncol(cl.rd)>2 & nrow(cl.rd) > 2){
    sim=cor(cl.rd)
  }
  else{
    cl.diff=as.matrix(dist(t(cl.rd)))
    sim = 1 - cl.diff/max(cl.diff)
  }
  return(sim)
}


get_cl_sim_multiple <- function(cl.rd.list, FUN =pmax)
  {
    all.cl = unique(unlist(lapply(cl.rd.list, colnames)))
    cl.sim = matrix(-1, nrow=length(all.cl),ncol=length(all.cl))
    colnames(cl.sim) = row.names(cl.sim) = all.cl
    cl.rd.list = cl.rd.list[!sapply(cl.rd.list, is.null)]
    for(cl.rd in cl.rd.list){
      if(nrow(cl.rd) <= 1){
        next
      }
      tmp.sim = get_cl_sim(cl.rd)
      cl.sim[row.names(tmp.sim),colnames(tmp.sim)] = FUN(tmp.sim, cl.sim[row.names(tmp.sim),colnames(tmp.sim)])
    }
    diag(cl.sim)=1
    return(cl.sim)
  }


##' ##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param cl.rd.list 
##' @param cl 
##' @param comb.dat 
##' @return 
##' @author Zizhen Yao
combine_cl_sim <- function(cl.rd.list, cl, comb.dat)
  {
    select.cl.rd.list = sapply(names(cl.rd.list),function(set){
      cl.rd=cl.rd.list[[set]]
      if(comb.dat$type=="mem"){
        cl.size = table(cl[names(cl) %in% colnames(comb.dat$dat.list[[set]])])
      }
      else{
        cl.size = table(cl[names(cl) %in% comb.dat$dat.list[[set]]$col_id])
      }
      select.cl = names(cl.size)[cl.size >= comb.dat$de.param.list[[set]]$min.cells]
      if(length(select.cl)==0){
        return(NULL)
      }      
      cl.rd[,select.cl,drop=F]    
    },simplify=FALSE)
    cl.sim = get_cl_sim_multiple(select.cl.rd.list)
  }

                                        #' Title
#'
#' @param dat.list 
#' @param de.param.list 
#' @param common.genes
#' @param cl 
#' @param pairs 
#' @param cl.means.list 
#' @param cl.present.list 
#' @param lfc.conservation.th 
#' @param de.genes.list 
#' @param max.cl.size 
#'
#' @return
#' @export
#'
#' @examples
de_genes_pairs_multiple <- function(dat.list, de.param.list,cl, pairs, cl.means.list, cl.present.list, cl.sqr.means.list, lfc.conservation.th=0.6, de.genes.list=NULL, max.cl.size=200, mc.cores=1, method="fast_limma")
  {

    if(is.null(de.genes.list)){
      de.genes.list = sapply(names(de.param.list), function(x)list())
    }
    
    cl.size.platform=list()
    new.de.genes.list=list()
    for(x in names(de.param.list)){
      dat = dat.list[[x]]
      if(is.list(dat)){
        tmp.cl = cl[names(cl) %in% dat$col_id]
      }
      else{
        tmp.cl = cl[names(cl) %in% colnames(dat)]
      }
      cl.size = table(tmp.cl)
      cl.size.platform[[x]]= cl.size
      select.cl = names(cl.size)[cl.size >= de.param.list[[x]]$min.cells]
      if(length(select.cl) < 2){
        next
      }
      tmp.cl = tmp.cl[tmp.cl %in% select.cl]
      if(is.factor(tmp.cl)){
        tmp.cl = droplevels(tmp.cl)
      }
      tmp.pairs= pairs[pairs[,1] %in% select.cl & pairs[,2] %in% select.cl,,drop=FALSE]
      if(nrow(tmp.pairs)==0){
        next
      }
      new.de.genes.list[[x]] = de_selected_pairs(norm.dat=NULL, cl= tmp.cl, pairs = tmp.pairs, de.param=de.param.list[[x]], method=method,cl.means=cl.means.list[[x]], cl.present = cl.present.list[[x]], cl.sqr.means=cl.sqr.means.list[[x]],mc.cores=mc.cores)$de.genes
    }
    require(foreach)
    ###find the conserved DE genes
    pairs = as.data.frame(pairs)
    pairs$pair = row.names(pairs)
    pairs$pair_bin = ceiling((1:nrow(pairs))/5000)

    mc.cores = min(mc.cores, max(pairs$pair_bin))
    registerDoMC(cores=mc.cores)
    top.n=1000
    cl.means.list = sapply(cl.means.list, as.matrix, simplify=FALSE)
    tmp = foreach::foreach(pid = unique(pairs$pair_bin),.combine="c") %dopar% {
      tmp=lapply(names(new.de.genes.list),function(set){
        de.genes = new.de.genes.list[[set]]
        select.pairs = pairs %>% filter(pair_bin==pid & pair %in% names(de.genes)) %>% pull(pair)
        de.df = lapply(select.pairs, function(p){
          if(is.null(de.genes[[p]])|de.genes[[p]]$num==0){
            return(NULL)
          }
          up = de.genes[[p]]$up.genes
          down = de.genes[[p]]$down.genes
          up = head(up, top.n)
          down = head(down, top.n)
          pair = strsplit(p, "_")[[1]]
          p1 = pair[1]
          p2 = pair[2]
          gene = c(names(up),names(down))
          logPval = c(up, down)
          sign = rep(c("up","down"),c(length(up),length(down)))
          df = data.frame(gene=gene, logPval=logPval,sign=factor(sign,c("up","down")))
          df$pair = p
          df$set = set
          df       
        })
        de.df = rbindlist(de.df)
      })
      de.df = rbindlist(tmp)
      if(nrow(de.df)==0){
        return(de.genes.list)
      }
      de.gene.df = de.df %>% select(pair, gene, sign) %>% distinct()
      de.gene.df = de.gene.df %>% left_join(pairs[,c("pair","P1","P2")],by="pair")
      de.gene.set.df = lapply(names(de.genes.list), function(set){
        tmp.genes = row.names(cl.means.list[[set]])
        tmp=de.gene.df %>% mutate(set=set)%>% filter(gene %in% tmp.genes)
        exp1 = get_pair_matrix(cl.means.list[[set]], tmp$gene, tmp$P1)
        exp2 = get_pair_matrix(cl.means.list[[set]], tmp$gene, tmp$P2)
        tmp = tmp %>% mutate(lfc = exp1-exp2) %>% filter(!is.na(lfc))
      })
      de.gene.set.df = rbindlist(de.gene.set.df)
      de.gene.df = de.gene.set.df %>% group_by(pair,gene,sign) %>% summarize(set.num=n(), lfc.num = sum(sign=="up"& lfc > 1 | sign=="down" & lfc < -1),.groups="drop")      
      de.gene.df = de.gene.df %>% filter(lfc.num >= lfc.conservation.th * set.num)
      de.df = de.df %>% right_join(de.gene.df,by=c("gene","sign","pair"))
      if(nrow(de.df)==0){
        return(de.genes.list)
      }      
      tmp.de.genes.list = with(de.df, tapply(1:nrow(de.df), list(set, pair), function(x){
        tmp.df = de.df[x,]
        up.genes = with(tmp.df %>% filter(sign=="up"), setNames(logPval, gene))
        down.genes = with(tmp.df %>% filter(sign=="down"), setNames(logPval, gene))        
        tmp = up.genes
        tmp[tmp > 20] = 20
        up.score <- sum(tmp)
        tmp = down.genes
        tmp[tmp > 20] = 20
        down.score <- sum(tmp)    
        list(
             up.genes=up.genes,
             down.genes=down.genes,
             up.score = up.score,
             down.score = down.score,
             score = up.score + down.score,
             up.num = length(up.genes),
             down.num = length(down.genes),
             num = length(up.genes) + length(down.genes)
             )        
      },simplify=FALSE))
      list(tmp.de.genes.list)
    }    
    for(set in names(de.genes.list)){      
      tmp.list = do.call("c",sapply(tmp, function(x){
        if(set %in% row.names(x)){
          y=x[set,]
          names(y)=colnames(x)
          y
        }
        else{
          NULL
        }
      },simplify=FALSE))
      for(p in setdiff(pairs$pair, c(names(de.genes.list[[set]]), names(tmp.list)))){
        tmp.list[[p]] = null_de()
      }
      de.genes.list[[set]] = c(de.genes.list[[set]], tmp.list)
    }
    return(de.genes.list)
  }    

de_genes_multiple <- function(comb.dat, cl, merge.sets=names(comb.dat$dat.list),mc.cores=10,...)
  {
    cl.stats.list = get_cl_stats_list(comb.dat, merge.sets, cl, max.cl.size=300,mc.cores=mc.cores)
    cl.means.list = cl.stats.list$cl.means.list
    cl.present.list = cl.stats.list$cl.present.list
    cl.sqr.means.list = cl.stats.list$cl.sqr.means.list
    de_genes_pairs_multiple(comb.dat$dat.list, comb.dat$de.param.list,cl=cl, pairs=create_pairs(unique(cl)), cl.means.list, cl.present.list, cl.sqr.means.list, de.genes.list=NULL,mc.cores=mc.cores,...)    
  }

  
get_cl_stats_list <- function(comb.dat, merge.sets, cl, max.cl.size=300,mc.cores=10, use.min.cells=TRUE)
  {
    cl.stats.list = list()
    meta.df=comb.dat$meta.df
    for(set in merge.sets){
      dat = comb.dat$dat.list[[set]]
      de.param = comb.dat$de.param.list[[set]]
      tmp.cl = cl[names(cl) %in% row.names(meta.df)[meta.df$platform==set]]
      tmp.size = table(tmp.cl)
      if(use.min.cells){
        tmp.cl= tmp.cl[tmp.cl %in% names(tmp.size)[tmp.size >= de.param$min.cells]]
      }
      if(is.factor(tmp.cl)){
        tmp.cl=droplevels(tmp.cl)
      }
      tmp.cells = sample_cells(tmp.cl, max.cl.size)
      if(length(tmp.cells)==0){
        cl.stats.list[[set]]=NULL
        next
      }
      tmp.cl = tmp.cl[tmp.cells]
      
      if(comb.dat$type=="mem"){
        tmp=sapply(c("means","present","sqr_means"), function(x){
          get_cl_stats(dat, cl=tmp.cl, max.cl.size = max.cl.size, stats=x)
        },simplify=F)
      }
      else{      
        if(length(tmp.cl)< 50000){
          dat = get_logNormal(dat, names(tmp.cl))
          tmp=sapply(c("means","present","sqr_means"), function(x){
            get_cl_stats(dat, cl=tmp.cl, max.cl.size = max.cl.size, stats=x)
          },simplify=F)
        }
        else{
          tmp=get_cl_stats_big(dat, tmp.cl, max.cl.size = max.cl.size, stats=c("means","present","sqr_means"),mc.cores=mc.cores)
        }
      }
      cl.stats.list[[set]]=tmp
    }
    cl.means.list  = sapply(cl.stats.list, function(x) as.data.frame(x$means),simplify=F)
    cl.present.list  = sapply(cl.stats.list, function(x) as.data.frame(x$present),simplify=F)
    cl.sqr.means.list  = sapply(cl.stats.list, function(x) as.data.frame(x$sqr_means),simplify=F)
    return(list(cl.means.list=cl.means.list, cl.present.list=cl.present.list, cl.sqr.means.list=cl.sqr.means.list))
  }

####Change criteria. If one of the platform shows significant DE genes, and the other platform show consistent fold change, keep the clusters seperate. 

#' Title
#'
#' @param comb.dat 
#' @param merge.dat.list 
#' @param cl 
#' @param anchor.genes 
#' @param verbose 
#' @param pairBatch 
#' @param de.genes.list 
#' @param lfc.conservation.th 
#' @param merge.type 
#'
#' @return
#' @export
#'
#' @examples
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param comb.dat 
##' @param merge.sets 
##' @param cl 
##' @param anchor.genes 
##' @param verbose 
##' @param pairBatch 
##' @param lfc.conservation.th 
##' @param merge.type 
##' @param de.method 
##' @param max.cl.size 
##' @param compare.k 
##' @return 
##' @author Zizhen Yao
merge_cl_multiple <- function(comb.dat, merge.sets, cl, anchor.genes=NULL, genes.allowed=NULL, joint.rd.dat=NULL, verbose=TRUE, lfc.conservation.th=0.7, merge.type="undirectional", de.method="fast_limma",max.cl.size=300, compare.k=4, mc.cores=5, pairBatch=100, cl.stats.list=NULL)
{
  #print("merge_cl_multiple")
  cl = setNames(as.character(cl),names(cl))

  merge_x_y <- function(x, y)
  {
    cat("merge", x, y, "\n")
    if(length(unique(cl))==1){
      return(NULL)
    }
    cl.size = sort(table(cl[cl %in% c(x,y)]),decreasing=T)
    p = names(cl.size)
    x = p[1]
    y = p[2]
    for(set in merge.sets){
      if(comb.dat$type=="mem"){
        tmp.cells = intersect(names(cl)[cl %in% c(x,y)], colnames(comb.dat$dat.list[[set]]))
      }
      else{
        tmp.cells = intersect(names(cl)[cl %in% c(x,y)], comb.dat$dat.list[[set]]$col_id)
      }
      tmp.cl.size = table(cl[tmp.cells])
      if(length(tmp.cl.size)>1){
        cl.means.list[[set]][[x]] = get_weighted_means(as.matrix(cl.means.list[[set]][, p]), tmp.cl.size[p])
        cl.sqr.means.list[[set]][[x]] = get_weighted_means(as.matrix(cl.sqr.means.list[[set]][, p]), tmp.cl.size[p])
        cl.present.list[[set]][[x]] = get_weighted_means(as.matrix(cl.present.list[[set]][, p]), tmp.cl.size[p])
        if(is.null(joint.rd.dat)){
          cl.rd.list[[set]][[x]] = get_weighted_means(as.matrix(cl.rd.list[[set]][,p]), cl.size[p])
        }
      }
      else if (y %in% names(tmp.cl.size)){
        cl.means.list[[set]][[x]] = cl.means.list[[set]][[y]] 
        cl.sqr.means.list[[set]][[x]] = cl.sqr.means.list[[set]][[y]]
        cl.present.list[[set]][[x]] = cl.present.list[[set]][[y]]
        if(is.null(joint.rd.dat)){
          cl.rd.list[[set]][[x]] = cl.rd.list[[set]][y]
        }
      }
      cl.means.list[[set]][[y]] = NULL
      cl.sqr.means.list[[set]][[y]] = NULL
      cl.present.list[[set]][[y]] = NULL
      if(is.null(joint.rd.dat)){
        cl.rd.list[[set]][[y]] = NULL
      }
    }
    cl[cl == y] = x
    if(!is.null(joint.rd.dat)){
      cl.rd[[x]] = get_weighted_means(as.matrix(cl.rd[,p]), cl.size[p])
      cl.rd[[y]] =  NULL
      cl.rd.list=NULL
    }
    return(list(cl=cl, cl.rd.list=cl.rd.list, cl.rd=cl.rd, cl.means.list=cl.means.list, cl.present.list = cl.present.list, cl.sqr.means.list = cl.sqr.means.list))
  }
  
  add_pairs_de_genes <- function(de.genes.list, new.pairs)
    {
      de.genes.list <- de_genes_pairs_multiple(comb.dat$dat.list[merge.sets], de.param.list= comb.dat$de.param.list[merge.sets], cl=cl, de.genes.list=de.genes.list, pairs=new.pairs, cl.means.list = cl.means.list, cl.present.list = cl.present.list, cl.sqr.means.list=cl.sqr.means.list, lfc.conservation.th=lfc.conservation.th, max.cl.size=max.cl.size,mc.cores=mc.cores)
      return(de.genes.list)
    }
  
  rm_pairs_de_genes <- function(de.genes.list, rm.pairs)
    {
      for(x in names(de.genes.list)){
        de.genes.list[[x]] = de.genes.list[[x]][setdiff(names(de.genes.list[[x]]), rm.pairs)]
      }
      return(de.genes.list)
    }
  
  test_merge_multiple <-  function(de.genes.list, merge.type="undirectional")
    {
      merge.pairs <- unique(unlist(lapply(de.genes.list, names)))      
      to.merge.df <- do.call("rbind", lapply(names(de.genes.list), function(x){
        to.merge <- sapply(merge.pairs, function(p){
          test_merge(de.genes.list[[x]][[p]], merge.de.param.list[[x]], merge.type=merge.type)
        })
        sc <- sapply(merge.pairs, function(p){de.genes.list[[x]][[p]]$score})
        sc[sapply(sc, is.null)] = 0
        sc = unlist(sc)
        df=data.frame(to.merge, sc, pair = merge.pairs,stringsAsFactors=FALSE)                     
      }))
      not.merged = with(to.merge.df, tapply(!to.merge, pair, sum))
      not.merged.pair = names(not.merged[not.merged > 0])
      to.merge.df = to.merge.df[!to.merge.df$pair %in% not.merged.pair,,drop=F]
      to.merge.sc =with(to.merge.df, tapply(sc, pair, max))
      return(sort(to.merge.sc))
    }
  
  merge.de.param.list = comb.dat$de.param.list[merge.sets]
  de.score.th = mean(sapply(merge.de.param.list, function(x)x$de.score.th))
  cl.platform.counts = table(comb.dat$meta.df[names(cl), "platform"],cl)[merge.sets,,drop=F]
  tmp = table(comb.dat$meta.df$platform)[merge.sets]
  cl.platform = cl.platform.counts / as.vector(tmp) 
  cl.platform = t(t(cl.platform) / colSums(cl.platform))  
  cl.min.cells = sapply(merge.de.param.list, function(x)x$min.cells)
  cl.big= cl.platform.counts >= cl.min.cells[rownames(cl.platform.counts)]
  cl.small = colnames(cl.big)[colSums(cl.big) == 0]
  cl.big =  colnames(cl.big)[colSums(cl.big) > 0]

  if(length(cl.big)==0){
    return(NULL)
  }

  ###Merge small clusters based on KNN prediction
  cl.small.cells= names(cl)[cl %in% cl.small]
  cl.big.cells= names(cl)[cl %in% cl.big]
  if(is.null(cl.stats.list)){
    cl.stats.list = get_cl_stats_list(comb.dat, merge.sets, cl,mc.cores=mc.cores, use.min.cells=FALSE)
  }
  cl.means.list      = cl.stats.list$cl.means.list
  cl.present.list    = cl.stats.list$cl.present.list
  cl.sqr.means.list  = cl.stats.list$cl.sqr.means.list
  if(!is.null(genes.allowed)){
    for(set in merge.sets){
      select.genes = intersect(row.names(cl.means.list[[set]]), genes.allowed)
      cl.means.list[[set]] = cl.means.list[[set]][select.genes, ]
      cl.present.list[[set]] = cl.present.list[[set]][select.genes, ]
      cl.sqr.means.list[[set]] = cl.present.list[[set]][select.genes, ]
    }
  }

  if(length(cl.small)>0){
    ######if there is a joint space, use it to compute cluster centroid for mapping
    ######else use marker gene expressiona as cluster centroid for mapping within each modality. 
    if(!is.null(joint.rd.dat)){
      cl.dat = get_cl_means(joint.rd.dat, cl)
      query.dat = t(joint.rd.dat[cl.small.cells,,drop=FALSE])
      map.df = map_cells_knn(query.dat, cl.dat[,cl.big,drop=FALSE])
      new.cl = setNames(map.df$cl, map.df$sample_id)      
    }
    else{       
      cl.small.cells.byplatform = split(cl.small.cells, as.character(comb.dat$meta.df[cl.small.cells, "platform"]))
      cl.big.cells.byplatform = split(cl.big.cells, comb.dat$meta.df[cl.big.cells, "platform"])
      unresolved.cells=c()
      new.cl = c()
      for(set in names(cl.small.cells.byplatform)){
        query.cells =cl.small.cells.byplatform[[set]]
        if(length(query.cells)==0){
          next
        }
        else{
          cl.dat = cl.means.list[[set]]
        }
        cl.dat = cl.dat[,colnames(cl.dat)%in% cl.big]
        if(is.null(cl.dat)|ncol(cl.dat)==0){
          unresolved.cells = c(unresolved.cells, query.cells)
          next
        }
        cl.dat = cl.dat[intersect(anchor.genes,row.names(cl.dat)),]
   OA     if(comb.dat$dat.list[[set]]$type=="mem"){
          dat = comb.dat$dat.list[[set]][row.names(cl.dat),query.cells]
          map.df = map_cells_knn(dat, cl.dat, mc.cores=mc.cores)
        }
        else{
          if(length(query.cells)<10000){
            dat = get_logNormal(comb.dat$dat.list[[set]], cols=query.cells, rows=row.names(cl.dat))
            map.df = map_cells_knn(dat, cl.dat, mc.cores=mc.cores)
          }
          else{
            map.df = map_cells_knn_big(comb.dat$dat.list[[set]], cl.dat, select.cells = query.cells, mc.cores=mc.cores)
          }
        }
        map.cl = setNames(map.df$cl, map.df$sample_id)      
        new.cl = c(new.cl, map.cl)
      }
      ####Used to assign all cells in a small clusters to the best matching big cluster. Now allow assign individual cells in a small cluster independently. 
      #if(length(unresolved.cells)>0){
      #  tb = table(cl[names(new.cl)], new.cl)
      #  map.cl = setNames(colnames(tb)[apply(tb, 1, which.max)], row.names(tb))
      #  new.cl[unresolved.cells] = map.cl[as.character(cl[unresolved.cells])]
      #}
    }
    ### Update cluster statistics based on updated cluster membership
    new.cl.stats.list   = get_cl_stats_list(comb.dat, merge.sets, new.cl, use.min.cells=FALSE)
    new.cl.means.list   = new.cl.stats.list$cl.means.list
    new.cl.present.list = new.cl.stats.list$cl.present.list
    new.cl.sqr.means.list = new.cl.stats.list$cl.sqr.means.list
    new.cl.platform.counts = table(comb.dat$meta.df[names(new.cl), "platform"],new.cl)[merge.sets,,drop=F]    
    ###update cl.stats
    for(set in names(new.cl.means.list)){
      tmp.cl.size = cbind(new.cl.platform.counts[set,], cl.platform.counts[set,colnames(new.cl.platform.counts)])
      for(p in row.names(tmp.cl.size)){
        if(new.cl.platform.counts[set,p] > 0){
          tmp.dat= cbind(new.cl.means.list[[set]][,p],cl.means.list[[set]][,p])
          cl.means.list[[set]][[p]] = get_weighted_means(tmp.dat, tmp.cl.size[p,])
          tmp.dat= cbind(new.cl.sqr.means.list[[set]][,p], cl.sqr.means.list[[set]][,p])
          cl.sqr.means.list[[set]][[p]] = get_weighted_means(tmp.dat, tmp.cl.size[p,])
          tmp.dat= cbind(new.cl.present.list[[set]][,p],cl.present.list[[set]][,p])
          cl.present.list[[set]][[p]] = get_weighted_means(tmp.dat, tmp.cl.size[p,])
        }
      }
    }
    cl[names(new.cl)]=new.cl    
  }
  ####Deprecated code for computating cl statistics from scratch. 
  #cl.stats.list = get_cl_stats_list(comb.dat, merge.sets, cl)
  #cl.means.list      = cl.stats.list$cl.means.list
  #cl.present.list    = cl.stats.list$cl.present.list
  #cl.sqr.means.list  = cl.stats.list$cl.sqr.means.list
  
  if(is.null(joint.rd.dat)){
    cl.rd.list = sapply(cl.means.list, function(x){
      x[row.names(x) %in% anchor.genes,,drop=F]
    },simplify=FALSE)
  }
  else{
    cl.rd = as.data.frame(get_cl_means(joint.rd.dat, cl))
  }
  
  de.pairs = NULL
  de.genes.list = sapply(merge.sets, function(x)list(),simplify=F)
  while(length(unique(cl)) > 1) {
    if(!is.null(joint.rd.dat)){
      cl.sim=get_cl_sim(cl.rd)
    }
    else{
      cl.sim=combine_cl_sim(cl.rd.list, cl, comb.dat)
    }
    if (length(cl.sim)==0) return(NULL)
    ###Find pairs of nearest neighbrs as candidates for merging.
    k.tmp = pmin(compare.k, ncol(cl.sim))
    nn=colnames(cl.sim)[sim_knn(cl.sim, k= k.tmp)[[1]]]
    
    merge.pairs = data.frame(P1=rep(row.names(cl.sim), length(k.tmp)), P2=nn,stringsAsFactors=FALSE)
    merge.pairs = merge.pairs[merge.pairs[,1]!=merge.pairs[,2],]
    merge.pairs$sim = get_pair_matrix(cl.sim, merge.pairs$P1, merge.pairs$P2)
    ##filter pairs that do not have overlapping modality.
    cl.size.platform = table(cl, as.character(comb.dat$meta.df[names(cl), "platform"]))
    merge.pairs.select = sapply(colnames(cl.size.platform), function(set){
      pmin(cl.size.platform[as.character(merge.pairs$P1), set], cl.size.platform[as.character(merge.pairs$P2), set]) >= comb.dat$de.param.list[[set]]$min.cells
    })
    merge.pairs.select = apply(merge.pairs.select, 1, any)
    merge.pairs = merge.pairs[merge.pairs.select,,drop=FALSE]
    tmp1 = pmin(merge.pairs[, 1], merge.pairs[, 2])
    tmp2 = pmax(merge.pairs[, 1], merge.pairs[, 2])
    merge.pairs$P1 = tmp1
    merge.pairs$P2 = tmp2
    
    p = paste(merge.pairs[, 1], merge.pairs[, 2], sep = "_")
    merge.pairs = merge.pairs[!duplicated(p), , drop = F]
    row.names(merge.pairs) = p[!duplicated(p)]
    merge.pairs = merge.pairs[order(merge.pairs$sim, decreasing = T), ,drop=F]
    merge.pairs = merge.pairs[!row.names(merge.pairs) %in% row.names(de.pairs),,drop=F]
    if(nrow(merge.pairs)==0){
      break
    }
    
    while(nrow(merge.pairs) > 0){
      cat("Merge.pairs", nrow(merge.pairs), "\n")
      new.pairs = head(row.names(merge.pairs),pairBatch)
      select.merge.pairs=merge.pairs[new.pairs, ,drop=F]
      if(is.null(de.pairs)){
        de.pairs = select.merge.pairs
      }else{
        de.pairs = rbind(select.merge.pairs,de.pairs)
      }
      de.genes.list = add_pairs_de_genes(de.genes.list, select.merge.pairs)
      if(is.null(de.genes.list) || sum(sapply(de.genes.list, length))==0){
        return(NULL)
      }
      merge.sc = test_merge_multiple(de.genes.list, merge.type=merge.type)
      if(length(merge.sc)>0){
        break
      }
      merge.pairs = merge.pairs[!row.names(merge.pairs)%in% new.pairs,]
    }
    if(length(merge.sc)==0){
      break
    }
    merged = c()   
    for(i in 1:length(merge.sc)){
      tmp = names(merge.sc[i])
      if(!tmp %in% row.names(de.pairs)){
        next
      }       
      p = de.pairs[tmp,]
      x = p[1,1]
      y = p[1,2]
      sim = p[,3]
      if (i == 1 | merge.sc[i] < de.score.th/2 & sum(c(x,y) %in% merged) == 0){              
        if (verbose > 0) {
          cat("Merge ", x, y, merge.sc[i], sim, length(which(cl == x)),  "cells", length(which(cl == y)), "cells", "\n")
        }         
        update.result=merge_x_y(x=x, y=y)
        if(is.null(update.result)){
          return(NULL)
        }
        cl = update.result$cl
        cl.rd.list = update.result$cl.rd.list
        cl.rd = update.result$cl.rd
        cl.means.list = update.result$cl.means.list
        cl.present.list = update.result$cl.present.list
        cl.sqr.means.list = update.result$cl.sqr.means.list
        rm.pairs = de.pairs[, 1] %in% c(x,y) | de.pairs[, 2] %in% c(x,y)
        de.genes.list = rm_pairs_de_genes(de.genes.list, row.names(de.pairs)[rm.pairs])
        de.pairs = de.pairs[!rm.pairs,,drop=F]        
        merged = c(merged, c(x,y))
      }
    }
  }
  if (length(unique(cl)) < 2) {
    return(NULL)
  }
  return(cl)
}

