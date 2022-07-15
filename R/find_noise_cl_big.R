library(arrow)

find_triplets_big <- function(de.summary="de_summary", all.pairs, select.cl = NULL, min.up.num=30, max.down.num=10, min.de.num=50)
  {
    ds= open_dataset(de.summary, partition ="pair_bin")
    select.pair = ds %>% filter(up.num > min.up.num & down.num < max.down.num & up.num  - down.num > min.up.num | down.num > min.up.num & up.num < max.down.num & down.num - up.num > min.up.num) %>% collect()
    select.pair = select.pair %>% filter(pair %in% all.pairs$pair)
    pairs = get_pairs(select.pair$pair)    
    pairs$pair = row.names(pairs)
    
    select.pair = select.pair %>% left_join(pairs)
    if(!is.null(select.cl)){
      select.pair = select.pair %>% filter(P1 %in% select.cl & P2 %in% select.cl)
    }
    select.pair = select.pair %>% mutate(cl.up= ifelse(up.num > down.num, P1, P2), cl.down =ifelse(up.num > down.num, P2, P1), up.num.new = ifelse(up.num > down.num, up.num, down.num), down.num.new = ifelse(up.num > down.num, down.num, up.num))
    
    tmp=table(select.pair$cl.up)
    select.cl = names(tmp)[tmp > 1]
    select.pair = select.pair %>% filter(cl.up %in% select.cl)
    
    tmp = select.pair[,c("cl.up","cl.down","up.num.new","down.num.new")]
    triplets = tmp %>% left_join(tmp, by="cl.up")
    triplets = triplets %>% filter(cl.down.x!=cl.down.y)
    triplets = triplets %>% mutate(P1 = ifelse(cl.down.x < cl.down.y, cl.down.x, cl.down.y),
      P2 = ifelse(cl.down.x < cl.down.y, cl.down.y, cl.down.x))
    triplets = triplets %>% left_join(all.pairs)
    diff.pair = ds %>% filter(up.num > min.de.num & down.num > min.de.num) %>% pull(pair)
    triplets = triplets %>% filter(pair %in% diff.pair) 
    triplets = triplets %>% mutate(
      pair1 = ifelse(cl.up < cl.down.x, paste0(cl.up,"_", cl.down.x), paste0(cl.down.x,"_", cl.up)),
      pair2 = ifelse(cl.up < cl.down.y, paste0(cl.up,"_", cl.down.y), paste0(cl.down.y,"_", cl.up)))
    tmp.cols=c("up.num.new.x","up.num.new.x","down.num.new.x","down.num.new.y")
    colnames(triplets)[match(tmp.cols,colnames(triplets))] = gsub("\\.new", "", tmp.cols)
    triplets = triplets %>% arrange(cl.up, down.num.x + down.num.y)
    return(triplets)
  }



get_de_pair_big<- function(de.df, cl1, cl2)
  {
    tmp.pair = create_pairs(c(cl1, cl2))
    de = de.df %>% filter(pair %in% row.names(tmp.pair))
    #flip the sign if cl1!=P1
    if(cl1!= tmp.pair[,1]){
      de = de %>% mutate(sign = ifelse(sign=="up","down","up"))
    }
    de$P1 = cl1
    de$P2 = cl2
    de
  }


check_triplet_big<- function(de.df, triplet,top.n=50)
  {
    cl = tirplet$cl.up
    cl1 = triplet$cl.down.x
    cl2 = triplet$cl.down.y        
    up.genes = with(de.df %>% filter(P1==cl1 & P2 ==cl2 & rank <= top.n), setNames(logPval, gene))
    down.genes = with(de.df %>% filter(P1==cl2 & P2 == cl1 & rank <=top.n), setNames(logPval, gene))
    up.genes.score= get_de_truncate_score_sum(up.genes)
    down.genes.score = get_de_truncate_score_sum(down.genes)
    
    tmp1.de = get_de_pair_big(de.df,triplet$cl.up,triplet$P1)
    tmp2.de = get_de_pair_big(de.df,triplet$cl.up,triplet$P2)


    tmp.genes=de.df %>% filter(P1==cl & P2==cl2) %>% pull(gene)
    olap.up.genes1 = intersect(tmp.genes, names(up.genes))
    olap.up.num1 = length(olap.up.genes1)
    olap.up.score1 = get_de_truncate_score_sum(up.genes[olap.up.genes1])   
    olap.up.ratio1 = olap.up.score1 / up.genes.score
    
    tmp.genes=de.df %>% filter(P1==cl & P2==cl1) %>% pull(gene)
    olap.down.genes1 = intersect(tmp.genes, names(down.genes))    
    olap.down.num1 = length(olap.down.genes1)
    olap.down.score1 = get_de_truncate_score_sum(down.genes[olap.down.genes1])
    olap.down.ratio1 = olap.down.score1 / down.genes.score
    
    up.genes2 =  with(de.df %>% filter(P1=cl1 & P2=cl& rank <= top.n), setNames(logPval, gene))
    up.genes.score2 = get_de_truncate_score_sum(up.genes2)
    olap.up.genes2 = intersect(names(up.genes2),names(up.genes))
    olap.up.num2 = length(olap.up.genes2)
    olap.up.score2 = get_de_truncate_score_sum(up.genes2[olap.up.genes2])
    olap.up.ratio2 = olap.up.score2 /up.genes.score2
    
    up.genes2 =  with(de.df %>% filter(P1=c2 & P2=cl& rank <= top.n), setNames(logPval, gene))
    down.genes.score2 = get_de_truncate_score_sum(down.genes2)
    olap.down.genes2 = intersect(names(down.genes2),names(down.genes))
    olap.down.num2 = length(olap.down.genes2)
    olap.down.score2 = get_de_truncate_score_sum(down.genes2[olap.down.genes2])
    olap.down.ratio2 = olap.down.score2 /down.genes.score2

    olap.num=c(olap.up.num1, olap.down.num1, olap.up.num2, olap.down.num2)
    olap.ratio = c(olap.up.ratio1, olap.down.ratio1, olap.up.ratio2, olap.down.ratio2)
    olap.score = c(olap.up.score1, olap.down.score1, olap.up.score2, olap.down.score2)
    names(olap.num) = paste0("olap.num.",c("up.1","down.1","up.2","down.2"))
    names(olap.ratio) = paste0("olap.ratio.",c("up.1","down.1","up.2","down.2"))
    score = sum(olap.score) / sum(c(up.genes.score, down.genes.score, up.genes.score2, down.genes.score2))
    
    result = list(
      cl = cl,
      cl1= cl1,
      cl2= cl2,
      up.num = length(up.genes),
      down.num = length(down.genes),
      score = score      
      )
    result = c(result, olap.ratio, olap.num)
    return(result)
  }


find_doublets_all_big <- function(de.dir, summary.dir = NULL, triplets=NULL, cl.bin, mc.cores=30,score.th=0.8, olap.th=1.6,out.dir="doublets_result",overwrite=TRUE,...)
  {
    require(parallel)
    require(doMC)
    require(foreach)
    library(data.table)
    library(dplyr)
    if(!dir.exists(out.dir)){
      dir.create(out.dir)
    }
    ds = open_dataset(de.dir, partition="pair_bin")
    if(is.null(triplets)){
      triplets=find_triplets_big(de.summary=summary.dir, cl.bin=cl.bin,...)
    }
    tmp = triplets %>% select("pair", "cl.up") %>% group_by(cl.up) %>% collect() %>% summarize(size=n())
    candidates = tmp %>% arrange(-size) %>% pull(cl.up)
    registerDoMC(cores=min(mc.cores,length(candidates)))
    mcoptions <- list(preschedule = FALSE)
    result.df=foreach::foreach(x=candidates,.options.multicore = mcoptions)%dopar% {
      fn = file.path(out.dir, paste0(x, ".data.parquet"))
      if(!overwrite & file.exists(fn)){
        result.df = read_parquet(fn)
        return(result.df)
      }
      tmp = triplets %>% filter(cl.up==x) %>% collect()
      cat(x, "triplets:", nrow(tmp),"\n")
      result.list= list()
      for(i in 1:nrow(tmp)){
        triplet = unlist(tmp[i,c("cl.up","cl.down.x","cl.down.y")])
        triplet.bin = cl.bin %>% filter(cl %in% triplet) %>% pull(bin)%in% unique
        de.df = ds %>% filter(bin.x %in% triplet.bin & bin.y %in% triplet.bin & P1 %in% triplet & P2 %in% triplet) %>% collect()
        result = check_triplet_big(de.df, triplet)
        result.list = c(result.list, list(result))
        if(result$score > score.th & result$olap.ratio.up.1 + result$olap.ratio.down.1 > olap.th){
          print(result)
          break
        }
      }
      fields = names(result.list[[1]])      
      result.df = do.call("data.frame",sapply(fields, function(f){
        sapply(result.list, function(x)x[[f]])
      },simplify=FALSE))
      write_parquet(result.df, sink=fn)
      return(NULL)
    }
    ds = open_dataset(out.dir)
    return(ds)
  }


find_low_quality_big <- function(ds, low.th=2,pairs)
  {
    library(arrow)
    library(dplyr)
    df = ds %>% filter(up.num < low.th | down.num < low.th) %>% collect
    df = df %>% left_join(pairs)
    df = df %>% mutate(cl=ifelse(up.num < low.th,P2, P1))
    df = df %>% mutate(cl.low=ifelse(up.num < low.th,P1, P2))
    return(df)
  }


plot_doublet_big <- function(big.dat, ds, cl, doublet.df, ...)
  {
    for(i in 1:nrow(doublet.df)){                                  
      x = as.character(doublet.df[i, "cl"])
      y = as.character(doublet.df[i, "cl1"])
      z = as.character(doublet.df[i, "cl2"])
      tmp.cl = cl[cl %in% c(x, y, z)]
      tmp.cl = setNames(factor(as.character(tmp.cl), c(x,y,z)), names(tmp.cl))
      tmp.cells= sample_cells(tmp.cl, 300)
      norm.dat = get_logNormal(big.dat, tmp.cells)
      tmp.pairs = create_pairs(unique(tmp.cl))
      markers = select_markers_ds(ds, pairs= row.names(tmp.pairs))
      tmp=display_cl(cl=tmp.cl[tmp.cells], norm.dat, prefix=paste0("doublet.",paste(levels(tmp.cl), collapse="_")), max.cl.size=100, markers=markers,...)
    }
  }

