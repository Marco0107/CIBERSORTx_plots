## Functions to plot CIBERSORTx results, still in development
## By Marco Pretti on February 28, 2021
## Last update on Sep 18, 2021
## Contact: marco-mp@hotmail.com or marco.pretti@inca.gov.br

plot.ciber.heat <- function(ciber.obj, ann_info, needHelp=F, is.Absolute=T, sample.column=1){
  
  if(needHelp==T){
    print("ciber.obj (required) is a df or tibble with cell types and sampleID in columns\n")
    print("ann_info is a data.frame with samples as rownames and infos to annotateas columns")
    print("is.Absolute assumes you ran CIBERSORTx in Absolute mode and normalize sample wise")
    print("sample.column define in which column is your sample ID")
  }
  require(data.table)
  require(RColorBrewer)
  require(dplyr)
  require(ggdendro)
  require(gridExtra)
  require(grid)
  require(ggplot2)
  require(cowplot)
  set.seed(123)

  # default values
  ciber.obj = as.data.frame(ciber.obj)
  sample.column=colnames(ciber.obj)[sample.column]
  value.columns <- colnames(ciber.obj)[!colnames(ciber.obj) %in% sample.column]
  
  ann.label1 <- colnames(ann_info)[1]
  if(ncol(ann_info)>1){ann.label2 <- colnames(ann_info)[2]} # to be tested
  ann_info[[sample.column]] <- rownames(ann_info)

  size.axis.Y = 12 # size do axis.text.y, mudar aqui para alterar no dend e barplot
  dend.Top = 50 # valor inversamente proporcional a altura do dendo, limite superior
  dend.Bot = 0 # valor do limite inferior da altura do dendo, alterar se for incluir multiplas labels
  height = F # T, F, plota valores de height pro dendo

  tema=list(theme_bw(),theme(axis.title=element_blank(), 
                             axis.text.y = element_text(size = size.axis.Y, vjust = 1), ### nao alterar
                             axis.text.x = element_blank(), #axis.text.y = element_blank(),
                             legend.position = "none",
                             plot.margin = unit(c(0,0,0,0),"mm")),
            scale_x_discrete(expand = c(0,0)))
  
  # normalize data
  if(is.Absolute==T){
    sum.ciber <- rowSums(ciber.obj %>% dplyr::select(-sample.column))
    ciber.obj[,value.columns] = ciber.obj[,value.columns]/sum.ciber*100
    scaleFUN <- function(x) sprintf("%.f", x*100) # , labels = scaleFUN
    scaleFUN.dend <- function(x) substr(round(x*10000, digits = 3),1,4) # , labels = scaleFUN
  }else{
    scaleFUN <- function(x) sprintf("%.f", x) # , labels = scaleFUN
    scaleFUN.dend <- function(x) round(x, digits = 1) # , labels = scaleFUN
  }
  
  # dendogram
  hc_complete <- hclust(as.dist(1-cor(t(ciber.obj[,value.columns]))), method="ward.D")
  ord <- hc_complete$order ## to change clustering order
  ciber.obj[[sample.column]] <- factor(ciber.obj[[sample.column]], levels = ciber.obj[[sample.column]][ord]) ## to change clustering order
  dendogram <- as.dendrogram(hc_complete)
  ddata <- ggdendro::dendro_data(dendogram, type = "rectangle")
  
  ## IMPROVE HERE - dend alignement left and right
  dend = ggplot(ggdendro::segment(ddata))+ 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
    theme(axis.title = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(size = size.axis.Y, vjust = -0.2), axis.ticks.x = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(dend.Top,0.5,dend.Bot,0.5),"mm"))+
    scale_y_continuous(labels = scaleFUN.dend, expand = c(0,0), position = "left")+
    scale_x_continuous(expand = c(0,0))
  if(height == F){dend = dend+theme(axis.text.y = element_text(color = "white"), axis.ticks.y = element_blank())}
  
  ## barplot
  colors <- c(brewer.pal(12,"Paired"), brewer.pal(8,"Dark2")) # suited for 20 cell types max
  p1=suppressWarnings(ciber.obj %>% melt() %>% 
                        ggplot(aes_string(x = sample.column, y = 'value', fill = 'variable'))+
                        geom_bar(width = 1, stat = "identity")+
                        labs(fill = "Subpopulation", y = "Relative percent")+
                        ylab("teste")+tema+
                        scale_y_continuous(expand = c(0,0), position = "left")+ ### nao alterar
                        guides(fill = guide_legend(ncol = 2))+
                        scale_fill_manual(values = colors))
  
  # barplot legend
  legend <- cowplot::get_legend(p1+theme(legend.text=element_text(size=12), legend.position="right", legend.key.size = unit(0.8,"line")))
  
  
  ## ann_label
  plot.ann_label <- function(data, ann_label, palette){
    # check size of palette and group
    if(brewer.pal.info[palette,][1] < length(unique(data[[ann_label]]))){
      print("Provided more labels than colors in palette, Ramping palette instead")
      palette = colorRampPalette(brewer.pal(name=palette, n = 8))(length(unique(data[[ann_label]]))) # 8 cobre todas
      color=scale_fill_manual(values = palette)
    }else{color <- scale_fill_brewer(palette = palette)}
    
    ## to match clustering order
    data[[sample.column]] <- factor(data[[sample.column]], levels = ciber.obj[[sample.column]][ord])
    
    ggplot(data = data, aes_string(x=sample.column, y='value', fill=ann_label))+
      geom_col(width = 1,position = "fill")+tema+
      theme(axis.text.x = element_blank(), axis.ticks = element_blank(),axis.text.y = element_text(colour = "white"))+
      scale_y_continuous(expand=c(0,0), position = "left", labels = scaleFUN)+color
    }
  
  if(length(ann_info) >2){
      p.ann1=plot.ann_label(data = suppressWarnings(ciber.obj %>% melt()) %>% inner_join(ann_info), ann_label = ann.label1, palette = "Set1")
      p.ann2=plot.ann_label(data = suppressWarnings(ciber.obj %>% melt()) %>% inner_join(ann_info), ann_label = ann.label2, palette = "Dark2")
      p.ann <- cowplot::plot_grid(p.ann2, p.ann1, ncol = 1)
      ann_legend1 = cowplot::get_legend(p.ann1+theme(legend.text=element_text(size=12), legend.position="right", legend.key.size = unit(0.8,"line"))+guides(fill=guide_legend(ncol=2)))
      ann_legend2 = cowplot::get_legend(p.ann2+theme(legend.text=element_text(size=12), legend.position="right", legend.key.size = unit(0.8,"line"))+guides(fill=guide_legend(ncol=2)))
      legends <- cowplot::plot_grid(ann_legend1, ann_legend2, legend, ncol=1, align = "v", rel_heights = c(2,2,2))
      
  }else{
    p.ann=plot.ann_label(data = suppressWarnings(ciber.obj %>% melt()) %>% inner_join(ann_info), ann_label = ann.label1, palette = "Set3")
    ann_legend <- cowplot::get_legend(p.ann+theme(legend.text=element_text(size=12), legend.position="right", legend.key.size = unit(0.8,"line"))+guides(fill=guide_legend(ncol=2)))
    legends <- cowplot::plot_grid(ann_legend, NULL, legend, ncol=1, align = "v", rel_heights = c(1,0,1))
  }

  ## Plot
  plots <- cowplot::plot_grid(dend, p.ann, p1, ncol=1, rel_heights = c(4,0.3,6))
  cowplot::plot_grid(plots, legends, rel_heights = c(2,1), ncol = 2, align = "hv")
}


plot.ciber.box <- function(ciber.obj, ann_info, needHelp=F, is.Absolute=T, sample.column=1,method="wilcox",
                           p.adj.paired.method="BH", is.paired=F, comparisons=NULL, exclude.zeros=TRUE,
                           exclude.manual=NULL){
  if(needHelp==T){
    print("ciber.obj (required) is a df or tibble with cell types and sampleID in columns\n")
    print("ann_info is a data.frame with samples as rownames and infos to annotateas columns")
    print("is.Absolute assumes you ran CIBERSORTx in Absolute mode and normalize sample wise")
    print("sample.column define in which column is your sample ID")
    print("p.adjust.method: see methods available at t_test and wilcox_test")
  }
  require(RColorBrewer)
  require(dplyr)
  require(gridExtra)
  require(grid)
  require(ggplot2)
  require(rstatix)
  require(data.table)
  
  # default values
  ciber.obj = as.data.frame(ciber.obj)
  sample.column=colnames(ciber.obj)[sample.column]
  value.columns <- colnames(ciber.obj)[!colnames(ciber.obj) %in% sample.column]
  
  ann.label1 <- colnames(ann_info)[1]
  ann_info[[sample.column]] <- rownames(ann_info)
  
  # comparisons
  if(is.null(comparisons) & is.paired==T){
    if(length(unique(ann_info[[ann.label1]])) > 2)
    print("You must provide the comparison list for a paired test with more than 2 levels !")
    break
  }
  
  if(is.null(comparisons) & is.paired==F){
    comparisons = unique(ann_info[[ann.label1]])
    split.comp <- function(x){ # combine comparisons two by two
      x = combn(x, 2)
      split(x, rep(1:ncol(x), each = nrow(x)))
    }
    comparisons = split.comp(x = comparisons)
    }

  # Filter data excluding zero values
  if(exclude.zeros==TRUE){
    ind = sapply(ciber.obj[value.columns], function(x) sum(x>0)) > nrow(ciber.obj)/2
    exclude = ind[!ind]
    print(paste("Removing cell types with > 50% zero observations:", names(exclude)))
    ciber.obj = ciber.obj %>% dplyr::select(-names(exclude))
  }
  
  # Manually exclude from boxplot
  if(!is.null(exclude.manual)){
    keep = colnames(ciber.obj)[!colnames(ciber.obj) %in% exclude.manual]
    ciber.obj = ciber.obj %>% dplyr::select(keep)
  }
  
  # calculate p-value
  print(paste("Using",method,"method to calculate difference between groups"))
  stat.test = suppressWarnings(ciber.obj %>% melt()) %>% inner_join(ann_info) %>%
    dplyr::rename(var='variable') %>%
    group_by(var)

  
  # Rename group during mean test
  colnames(stat.test)[ncol(stat.test)] <- "Group"
  
  if(method=="wilcox" | method=="wilcoxon"){
    stat.test = wilcox_test(stat.test, value ~ Group, paired = is.paired, comparisons = comparisons,
                            p.adjust.method = p.adj.paired.method) %>% 
      add_y_position() %>% mutate(variable=var, Group=group1)

    print(comparisons)
    }else if(method=="t.test" | method=="student"| method=="t"){
      stat.test=t_test(data = stat.test, value ~ Group, paired = is.paired, comparisons = comparisons,
                       p.adjust.method = p.adj.paired.method) %>% 
        add_y_position() %>% mutate(variable=var, Group=group1)
    }
    colnames(stat.test)[ncol(stat.test)] <- ann.label1

#  if(is.paired==T){
    stat.test=stat.test %>%
      dplyr::filter(p.adj.signif != "ns") 
#  }else{
  #   stat.test=stat.test %>%
  #     dplyr::filter(p != "ns") %>% dplyr::rename(p.corrected = p)
  # }
  
  # boxplot - ## ADICIONAR ncol pro facet AQUI ########
  boxplot <- suppressWarnings(ciber.obj %>% melt()) %>% inner_join(ann_info) %>%
    ggplot(aes_string(x=ann.label1, y='value', fill=ann.label1))+geom_boxplot()+
    facet_wrap(.~variable, scales="free_y", ncol=3, shrink = T)+ylab("Amount of cell type")+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    scale_fill_brewer(palette = "Set3")
  
  if(nrow(stat.test) > 0){
    boxplot=boxplot+
      geom_text(data = stat.test, aes(y=`y.position`*0.55, label=paste("p=",`p.adj`), x=group1), hjust=0)+ ## corrigir o p
      geom_segment(data = stat.test, aes(x = group1, xend=group2, y = `y.position`*0.5, yend=`y.position`*0.5))
  }
  return(boxplot)
}


## Correlation matrix

plot.ciber.corr <- function(ciber.obj, ann_info=NA, needHelp=F, is.Absolute=T, sample.column=1, cor.type="spearman", exclude_zeros=TRUE){
  if(needHelp==T){
    print("ciber.obj (required) is a df or tibble with cell types and sampleID in columns\n")
    print("ann_info is a data.frame with samples as rownames and infos to annotateas columns")
    print("is.Absolute assumes you ran CIBERSORTx in Absolute mode and normalize sample wise")
    print("sample.column define in which column is your sample ID")
    print("cor.type: see rcorr for 'types' of correlation supported")
  }
  require(RColorBrewer)
  require(dplyr)
  require(cowplot)
  require(pheatmap)
  require(corrplot)
  require(Hmisc)
  require(ggplot2)
  
  # default values
  set.seed(123)
  ciber.obj = as.data.frame(ciber.obj)
  sample.column=colnames(ciber.obj)[sample.column]
  value.columns <- colnames(ciber.obj)[!colnames(ciber.obj) %in% sample.column]
  
  corr.matrix <- function(x){
    
    if(exclude_zeros==TRUE){
      # Filter data excluding zero values
      ind = sapply(x[value.columns], function(x) sum(x>0)) > nrow(x)/2
      exclude = ind[!ind]
      if (length(exclude >0)) {
        print(paste("Removing cell types with > 50% zero observations:", names(exclude)))
        x = x %>% dplyr::select(-names(exclude))
      }
    }else{
      if(cor.type!='spearman'){
        print("Attention!!\nDon't include zeros when performing Pearson correlation")
      }
    }
    
    # Correlation matrix
    ciber.cor <- rcorr(as.matrix(x %>% dplyr::select(-all_of(sample.column))), type = cor.type)
    ciber.cor$P[ciber.cor$P > 0.05] <- 0
    ciber.cor$r[ciber.cor$P > 0.05] <- 0
    ciber.cor$r[is.na(ciber.cor$r)] <- 0
    return(ciber.cor)
  }
  
  # Mode: comparing groups
  if(is.data.frame(ann_info)){
    print("Mode: Comparing groups")
    ann.label1 <- colnames(ann_info)[1]
    ann_info[[sample.column]] <- rownames(ann_info)
    
    # Create and Assign different matrices for each group
    lista <- lapply(unique(ann_info[[1]]), function(y){
      ann_info.sub = ann_info[ann_info[,1] %in% y,]
      ciber.obj.sub <- inner_join(ciber.obj, ann_info.sub)
      return(corr.matrix(x = ciber.obj.sub %>% dplyr::select(-ann.label1)))
    }); names(lista) <- unique(ann_info[[1]])
    
    # Plot
    lista = suppressWarnings(melt(lapply(lista, '[', 1)) %>%
     ggplot(aes(y=Var1, x=value, fill=L1))+geom_point(shape=23, color="black")+facet_wrap(.~Var2)+
      geom_vline(xintercept = 0, linetype = "dashed")+
      ylab("Cell types")+xlab(paste(cor.type, "correlation coeff"))+
      scale_fill_brewer(palette = "Set3")+
      theme(axis.title = element_text(size=12), legend.title = element_blank(),
            strip.text = element_text(size=10)))
    
    return(cowplot::plot_grid(lista))
  }

  # Mode: Correlation matrix
  if (!is.data.frame(ann_info)) {
    print("Mode: Correlation matrix")
    ciber.cor <- corr.matrix(ciber.obj)
    pheatmap(ciber.cor$r, scale = "none", angle_col = 45, breaks = seq(-1,1,by=0.02),
                       color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100))    
  }
  
}
