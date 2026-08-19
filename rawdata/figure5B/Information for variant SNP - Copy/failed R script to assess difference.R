library(tidyverse)
library(GenomicRanges)
library(magrittr)
library(knitr)

#https://jmonlong.github.io/Hippocamplus/2017/09/19/mummerplots-with-ggplot2/
#Tutorial for how I did this! 

setwd("C:/Users/cf466/OneDrive - University of Exeter/Documents/Documents - CF laptop/PhD/Phage Cocktails/Section 2 - Optimisation/Figure 1 - Clarification/genomic differences between CPL272 and var272")

readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

mumgp = readDelta("CPL_vs_var_differences.delta")

mumgp %>% head %>% kable

filterMum <- function(df, minl=1000, flanks=1e4){
  coord = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) %>%
    summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
    ungroup %>% arrange(desc(rs)) %>%
    mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
  merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
    mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}

mumgp.filt = filterMum(mumgp, minl=1e4)
#mumgp.filt %>% head %>% kable

ggplot(mumgp.filt, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(qid~., scales='free', space='free', switch='y') +
  theme_bw() + theme(strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.99,.01), legend.justification=c(1,0),
                     strip.background=element_blank(),
                     axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab('PipetteDweller Sequence (65.8 kb)') + ylab('PipetteDweller Variant Sequence') + scale_colour_brewer(palette='Set1')

mumgp %<>% mutate(similarity=1-error/abs(qe-qs))
mumgp.filt %<>% mutate(similarity=1-error/abs(qe-qs))

ggplot(mumgp, aes(x=rs, xend=re, y=similarity, yend=similarity)) + geom_segment() +
  theme_bw() + xlab('reference sequence') + ylab('similarity') + ggtitle('All contigs') +
  ylim(0,1)

