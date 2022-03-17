find_low_quality_ds <- function(ds, low.th=2)
  {
    library(arrow)
    library(dplyr)
    df = ds %>% filter(up.num < low.th | down.num < low.th) %>% collect
    pairs = get_pairs(df$pair)
    pairs$pair = row.names(pairs)
    df = df %>% left_join(pairs)
    df = df %>% mutate(cl=ifelse(up.num < low.th,P2, P1))
    df = df %>% mutate(cl.low=ifelse(up.num < low.th,P1, P2))
    return(df)
  }
