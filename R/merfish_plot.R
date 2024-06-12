
plot_sections <- function(merfish.anno, 
                          meta.f, 
                          meta.col, 
                          select.cells, 
                          coord.to.use = "corrected",
                          negative.y = T,
                          frac.overlap = 0.05,
                          sections, 
                          ncol=length(sections), 
                          layout=NULL, 
                          plot.option=c("all","left","right","folded","crop"), 
                          label.top=4, 
                          fg.alpha=1, 
                          fg.cex=1,
                          combine.plots=TRUE,
                          frac.crop = c(0,0,0,0), #crop x(l), x(r), y(t), y(b)
                          label.cex = 1,
                          plot.ccf = NULL,
                          ccf.column = "CCF_acronym",
                          ...)
{
  x_coord <- paste0(coord.to.use,"_x")
  y_coord <- paste0(coord.to.use,"_y")
  
  rd.dat = as.data.frame(merfish.anno[,c(x_coord,y_coord)])  
  if(isTRUE(negative.y)){
    rd.dat[,2] = - rd.dat[,2]
  }
  
  #xlims <- c(min(rd.dat[,1]), max(rd.dat[,1]))
  #ylims <- c(min(rd.dat[,2]), max(rd.dat[,2]))
  
  
  row.names(rd.dat) = merfish.anno$cell_id
  meta = setNames(merfish.anno[[meta.f]], merfish.anno$cell_id)
  meta.section = merfish.anno %>% filter(cell_id %in% select.cells & section %in% sections) %>% select(one_of(c(meta.f, "section"))) %>% group_by_at(vars(one_of(meta.f)))  %>% summarise(select.section = names(which.max(table(section))))
  select.meta.section = setNames(meta.section[["select.section"]],meta.section[[meta.f]])
  
  g.list=list()
  for(sec in sections){
    # print(sec)
    all.cells = intersect(merfish.anno %>% filter(section==sec) %>% pull(cell_id),row.names(rd.dat))      
    select.rd.dat = rd.dat[all.cells,]      
    if(plot.option=="folded"){
      select.rd.dat[,1] =abs(select.rd.dat[,1])
    }      else if(plot.option=="left"){
      #cutoff = rearrr::centroid(select.rd.dat[,x_coord]) + (rearrr::centroid(select.rd.dat[,x_coord])*frac.overlap)
      cutoff = mean(select.rd.dat[,x_coord]) + mean(select.rd.dat[,x_coord])*frac.overlap
      select.rd.dat = select.rd.dat[select.rd.dat[,1]<=cutoff,,drop=F]
    } else if(plot.option=="right"){
      cutoff = mean(select.rd.dat[,x_coord]) - mean(select.rd.dat[,x_coord])*frac.overlap
      select.rd.dat = select.rd.dat[select.rd.dat[,1]>= cutoff,,drop=F]
    } else if(plot.option == "crop"){
      
      min.x = min(select.rd.dat[,x_coord])
      max.x = max(select.rd.dat[,x_coord])
      cutoff.xl = min.x + ((max.x-min.x)*frac.crop[1])
      cutoff.xr = min.x + ((max.x-min.x)*(1-frac.crop[2]) )
      
      select.rd.dat = select.rd.dat[select.rd.dat[,x_coord]>=cutoff.xl,,drop=F]
      select.rd.dat = select.rd.dat[select.rd.dat[,x_coord]<=cutoff.xr,,drop=F]
      
      min.y = min(select.rd.dat[,y_coord])
      max.y = max(select.rd.dat[,y_coord])       
      cutoff.yb = min.y + ((max.y-min.y)*frac.crop[4])
      cutoff.yt = min.y + ((max.y-min.y)*(1-frac.crop[3]) )
      
      select.rd.dat = select.rd.dat[select.rd.dat[,y_coord]>=cutoff.yb,,drop=F]
      select.rd.dat = select.rd.dat[select.rd.dat[,y_coord]<=cutoff.yt,,drop=F]
      
    }
    fg.cells = intersect(select.cells, row.names(select.rd.dat))
    cl.size = table(meta[fg.cells])      
    if(length(fg.alpha)>1){
      select.fg.alpha = fg.alpha[fg.cells]        
    }      else{
      select.fg.alpha= fg.alpha
    }
    if(length(fg.cex)>1){
      select.fg.cex = fg.cex[fg.cells]        
    }      else{
      select.fg.cex = fg.cex
    }
    label.meta = names(select.meta.section)[select.meta.section==sec]      
    label.meta = union(label.meta, names(tail(sort(cl.size[cl.size > 0]), label.top)))
    g.list[[sec]]=plot_RD_highlight(select.rd.dat,
                                    meta=meta,
                                    meta.col=meta.col,
                                    label.meta = label.meta,
                                    fg.cells = fg.cells,
                                    fg.alpha=select.fg.alpha, 
                                    fg.cex=select.fg.cex, 
                                    label.cex=label.cex,
                                    #xlims = xlims,
                                    #ylims = ylims,
                                    ...
    )      
  }
  
  if(!is.null(plot.ccf)){
    require(ggforce)
    for(sect in names(g.list)){
      print(sect)
      
      sel.sect <- sect
      hull.df <- merfish.anno %>% 
        filter(section %in% sect)
      
      if(isTRUE(negative.y)){
        hull.df[,y_coord] = - hull.df[,y_coord]
      }
      
      plot.meta <- g.list[[sect]][["data"]]
      
      plot.meta <- plot.meta %>% 
        left_join(hull.df[,c(x_coord,y_coord,"CCF_bin_x", "CCF_bin_y",ccf.column)], by=c("Dim1"=x_coord, "Dim2"=y_coord) ) 
      
      g.list[[sect]] <- g.list[[sect]] + #theme_bw() + 
        ggnewscale::new_scale_color() +
        geom_mark_hull(data = plot.meta %>% ungroup(), inherit.aes = F,
                       aes( x=Dim1, y=Dim2, filter = get(ccf.column) %in% plot.ccf),
                       lwd = 0.1,
                       concavity = 5,
                       expand=0,
                       radius=0
        )
    }
  }
  
  if(combine.plots){
    multiplot(plotlist= g.list, cols = ncol, layout=layout,byrow=TRUE)
    #gridExtra::grid.arrange(grobs=g.list,ncol=ncol,top=textGrob(plotTitle))
  }
  else{
    return(g.list)
  }
}


plot_gene_sections <- function(merfish.anno, big.dat=NULL, gene.exp=NULL, gene,  sections, ncol=length(sections), colorset= c("#ffffd9","#c7e9b4","#41b6c4","#1d91c0","#225ea8","#081d58"),plot.option=c("left","right","both"),...)
{
  rd.dat = merfish.anno[,c("global_x","global_y")]
  row.names(rd.dat) = merfish.anno$cell_id
  if(is.null(gene.exp)){
    gene.exp = get_cols(big.dat, rows=gene, cols = big.dat$col_id,sparse=FALSE)
  }
  meta.bin  = setNames(cut(gene.exp[gene,], 100), colnames(gene.exp))
  meta.col = setNames(colorRampPalette(colorset)(100), levels(meta.bin))
  g.list=list()
  for(s in sections){      
    all.cells = merfish.anno %>% filter(section==s) %>% pull(cell_id)
    select.rd.dat = rd.dat[all.cells,]
    if(plot.option=="left"){
      select.rd.dat = select.rd.dat[select.rd.dat[,1]<=0.5,,drop=F]
    }else if(plot.option=="right"){
      select.rd.dat = select.rd.dat[select.rd.dat[,1]>=-0.5,,drop=F]
    }
    all.cells = row.names(select.rd.dat)
    fg.cells = all.cells[gene.exp[gene, all.cells] > 0]
    g.list[[paste0(s,gene)]] =plot_RD_highlight(select.rd.dat, meta= meta.bin[all.cells], meta.col=meta.col, show.legend="none",raster=TRUE, fg.cells=fg.cells,...)
  }
  return(g.list)
}




