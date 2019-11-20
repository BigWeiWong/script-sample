## ---------------------------
## Script name: sample.R
## Purpose of script:
## Author: Qiongzi Qiu
## Date Created: 2019-01-01
## Copyright (c) Qiongzi Qiu, 2019
## Email: 3110102477@zju.edu.cn
## ---------------------------
## Notes: blahblah
## ---------------------------

setwd("work dictionary")

library(ggplot2)
library(reshape2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

source("some package")

options(stringsAsFactors = FALSE)

##### boxplot sample #####
infile_list = c('replace',
              'replace',
              'replace')
for(infile in infile_list){
  dat <- read.table(infile, header=T, sep='\t')
  mutation_type = strsplit(infile, '\\.')[[1]][2]
  dat['replace'] <- ifelse(dat$KO_non_pval_adj<0.05 & dat$KO_non_t_statistic>0 & !is.na(dat$KO_non_t_statistic), TRUE, FALSE)
  dat_melt = melt(dat[,c('replace', 'replace', 'replace', 'replace')], id.vars='replace')
  dat_melt[dat_melt=='replace'] <- gsub('_', ' ', mutation_type)
  dat_melt[dat_melt=='replace'] <- 'replace'
  my_comparisons = list(c('replace', 'replace'), c('replace', 'replace'), c('replace', 'replace'))
  
  boxplot_tmp = ggboxplot(dat_melt, 
                          x='variable', y='value',
                          xlab=NULL, ylab='replace') +
    stat_compare_means(comparisons = my_comparisons,
                       method.args=list(alternative='greater')) +
    facet_wrap( ~ replace, scales="fixed")
  outfile = paste0(mutation_type, '_replace.pdf')
  ggsave(outfile, boxplot_tmp)
}

##### heatmap sample #####
infile_list = c('replace',
                'replace',
                'replace',
                'replace')
for(infile in infile_list){
  dat = read.table(infile, header=T, sep='\t')
  dat[,'replace'] = dat[,'replace']/(dat[,'replace'] + dat[,'replace'])
  dat[is.na(dat[,'replace']),'replace'] = 'replace'
  for(type_tmp in c(TRUE, FALSE)){
    if(type_tmp==TRUE){
      dat_tmp = dat[dat$KO_non_up>0, ]
      outfile = gsub('replace', 'replace', infile)
    }else{
      dat_tmp = dat[dat$KO_non_prop==1, ]
      outfile = gsub('replace', 'replace', infile)
    }
    
    dat_reshape = reshape(dat_tmp[,c('replace', 'replace', 'replace')], idvar='replace', timevar='replace', direction='wide')
    rownames(dat_reshape) = dat_reshape[,1]
    dat_reshape = t(dat_reshape[,-1])
    count_row = as.numeric(apply(is.na(dat_reshape)==F, 1, sum))
    count_col = as.numeric(apply(is.na(dat_reshape)==F, 2, sum))
    dat_reshape = dat_reshape[count_row>0, count_col>0]
    dat_reshape = t(apply(dat_reshape, 1, as.numeric))

    dat_reshape[dat_reshape==1]=2
    dat_reshape[dat_reshape<1]=1
    dat_reshape[is.na(dat_reshape)]=0
    count_row = as.matrix(cbind(as.numeric(apply(dat_reshape==2, 1, sum)),
                      as.numeric(apply(dat_reshape==1, 1, sum))))
    count_col = as.matrix(cbind(as.numeric(apply(dat_reshape==2, 2, sum)),
                      as.numeric(apply(dat_reshape==1, 2, sum))))
    
    annotation_row = rowAnnotation(Count=row_anno_barplot(count_row,
                                   gp = gpar(fill = c("black", "grey"))))
    annotation_col = HeatmapAnnotation(Count = anno_barplot(count_col,
                                       gp = gpar(fill = c("black", "grey"),
                                                 col = c(NA, NA))))
    
    pdf(outfile, width=6, height=6)
    p = Heatmap(as.matrix(dat_reshape), 
                col = c('white', 'grey', 'black'),
            show_column_names = F, na_col='white',
            row_order = order(rowSums(count_row), decreasing = T),
            column_order = order(rowSums(count_col), decreasing = T),
            top_annotation = annotation_col,
            right_annotation = annotation_row,
            heatmap_legend_param = list(
              title = "", at = c(2, 1, 0), border = "black",  ncol = 3,
              labels = c("replace", "replace", "replace")
            ))
    draw(p, heatmap_legend_side = "bottom")
    dev.off()
  }
}

