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
    de = de.df %>% filter(pair == triplet$pair)    
    up.genes = with(de %>% filter(sign=="up" & rank <= top.n), setNames(logPval, gene))
    down.genes = with(de %>% filter(sign=="down" & rank <=top.n), setNames(logPval, gene))
    up.genes.score= get_de_truncate_score_sum(up.genes)
    down.genes.score = get_de_truncate_score_sum(down.genes)
    
    tmp1.de = get_de_pair_big(de.df,triplet$cl.up,triplet$P1)
    tmp2.de = get_de_pair_big(de.df,triplet$cl.up,triplet$P2)


    tmp.genes=tmp2.de %>% filter(sign=="up") %>% pull(gene)
    olap.up.genes1 = intersect(tmp.genes, names(up.genes))
    olap.up.num1 = length(olap.up.genes1)
    olap.up.score1 = get_de_truncate_score_sum(up.genes[olap.up.genes1])   
    olap.up.ratio1 = olap.up.score1 / up.genes.score

    tmp.genes=tmp1.de %>% filter(sign=="up") %>% pull(gene)
    olap.down.genes1 = intersect(tmp.genes, names(down.genes))    
    olap.down.num1 = length(olap.down.genes1)
    olap.down.score1 = get_de_truncate_score_sum(down.genes[olap.down.genes1])
    olap.down.ratio1 = olap.down.score1 / down.genes.score
    
    up.genes2 =  with(tmp1.de %>% filter(sign=="down" & rank <= top.n), setNames(logPval, gene))
    up.genes.score2 = get_de_truncate_score_sum(up.genes2)
    olap.up.genes2 = intersect(names(up.genes2),names(up.genes))
    olap.up.num2 = length(olap.up.genes2)
    olap.up.score2 = get_de_truncate_score_sum(up.genes2[olap.up.genes2])
    olap.up.ratio2 = olap.up.score2 /up.genes.score2
    
    
    down.genes2 = with(tmp2.de %>% filter(sign=="down" & rank <= top.n), setNames(logPval, gene))
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
      cl = triplet[,"cl.up"],
      cl1= triplet[,"P1"],
      cl2= triplet[,"P2"],
      up.num = length(up.genes),
      down.num = length(down.genes),
      score = score      
    )
    result = c(result, olap.ratio, olap.num)
    return(result)
  }


find_doublets_all_big <- function(de.dir, summary.dir = NULL, triplets.fn=NULL, all.pairs, mc.cores=40,score.th=0.8, olap.th=1.6,out.dir="doublets_result",overwrite=TRUE,...)
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
    if(!is.null(triplets.fn)){
      triplets = open_dataset(triplets.fn)
    }
    else{
      triplets=find_triplets_big(de.summary=summary.dir, all.pairs, ...)
    }
    tmp = triplets %>% select("pair", "cl.up") %>% group_by(cl.up) %>% collect() %>% summarize(size=n())
    candidates = tmp %>% arrange(-size) %>% pull(cl.up)
    registerDoMC(cores=min(mc.cores,length(candidates)))
    mcoptions <- list(preschedule = FALSE)
    result.df=foreach::foreach(x=candidates,.combine="rbindlist",.options.multicore = mcoptions)%dopar% {
      fn = file.path(out.dir, paste0(x, ".data.parquet"))
      if(!overwrite & file.exists(fn)){
        result.df = read_parquet(fn)
        return(result.df)
      }
      tmp = triplets %>% filter(cl.up==x) %>% collect()
      cat(x, "triplets:", nrow(tmp),"\n")
      result.list= list()
      for(x in 1:nrow(tmp)){
        triplet = tmp[x,]
        pairs =  c(triplet$pair, triplet$pair1, triplet$pair2)
        tmp.pair.bin=all.pairs %>% filter(pair %in% pairs) %>% pull(pair_bin) %>% unique
        de.df = ds %>% filter(pair_bin %in% tmp.pair.bin & pair %in% pairs) %>% collect()
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
      return(result.df)
    }
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
