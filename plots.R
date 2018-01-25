library(viridis)
library(tidyverse)
library(ggbeeswarm)
library(reshape2)
library(ggdendro)

magma <- c('#424949', '#F9A825', '#CD3F0A', '#8A023C')
redtogreen <- c('#550527', '#A10702', '#F44708', '#F9A113', '#599124')
greentored <- c('#599124' , '#F9A113', '#F44708', '#A10702', '#550527')
material <- c('#15889C', '#ED496F', '#8E1382', '#C62828', '#FFB300', '#43A047', '#FF6F00')
vibrant <- c('#95190C', '#FFB300', '#1F7530', '#086788', '#6300B5', '#EF8737', '#550527')
blass <- c('#e67e22', '#c0392b', '#6c3483', '#2471a3', '#229954', '#d4ac0d')
magenta <- c('#D96DED', '#6300B5', '#8E1382', '#ED496F', '#FF9D7C')
bluelime <- c('#012824', '#15889C', '#7BCFC0', '#88E200', '#028F03')
bricksky <- c('#0B3954', '#28AAE1', '#FDC120', '#D94551', '#55545C')
colorful <- c('#9E0031', '#FFBB00', '#198749', '#2D62A3', '#8E1382', '#5A0002')
flat <- c('#471b9a' , '#ff6f00', '#1565c0', '#43a047', '#ffb300', '#c62828')

three <- c('#16A085', '#2980B9', '#8E44AD')

many <- c('#Ae0031', '#ffbb00', '#198749', '#2d62a3', '#573794','#Ff7708', '#15889C', '#ED49B1', '#8E13A2', '#ead9d5', '#666686', '#CD6155', '#F7DC6F', '#7DCEA0', '#85C1E9', '#EB984E', '#03cd4A', '#CDDC39', '#e0594b', '#c76cde', '#24B177', '#8D6E63', '#486ff7', '#6300b5', '#88e200', '#012824', '#0d3290', '#ead9d5', '#a347fb', '#54fc7a', '#eb1388', '#b0978d', '#fe52cf', '#83f1f6', '#f1f847', '#2b1dfc', '#6c6f15', '#6ca05c', '#7788cd', '#f502f3', '#0dc290', '#fa0e03', '#3caa0a', '#befc8d', '#08f8eb', '#b1cd3f', '#d6a5fa', '#ce606c', '#ab1eba', '#6ecc9f', '#054ddc', '#486ff7', '#854f49', '#f22B21', '#3a0e43', '#225805', '#37d160', '#e4b974', '#a8bade', '#47edd1', '#f47a92', '#c76cde', '#9106eb', '#81aa20', '#d7fdfd', '#5deb2e', '#f82745', '#6435e0', '#027ffe', '#8e3101', '#16f648', '#1c15bc', '#8be46e', '#8d6fa0', '#e68fc6', '#058ca9', '#9e018a', '#bdfd0b', '#b22760', '#2bf49f', '#cb9348', '#9d8303', '#c251a1', '#46adaf', '#a3e3af', '#22bb34', '#6ea3fa', '#260374', '#1c3854', '#405d37', '#c21df3', '#fcea92', '#537f88', '#fd4c18', '#f2d71e', '#fd4c7a')

insert_minor <- function(major_labs, n_minor) {
  labs <- c( sapply( major_labs, function(x) c(x, rep("", n_minor))))
  labs[1:(length(labs)-n_minor)]
}

color_gradient <- function(start='white', end='#550527', n=100){
  return(colorRampPalette(c(start, end))(n))
}

color_gradient_scrambled <- function(start='white', end='#550527', n=100){
  return(sample(colorRampPalette(c(start, end))(n)))
}

bar_plot <- function(data, x, y, fill, xlab, ylab, flab, title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, fill=fill)) + 
    ggtitle(title) + 
    geom_bar(stat = 'summary', fun.y = "mean") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    labs(x=xlab, y=ylab, fill=flab)
    if (is.null(colorscheme)){
      p <- p + scale_fill_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_fill_manual(values=colorscheme)
    }
  return(p)
}

beeswarm_plot <- function(data, x, y, color, alpha=1, xlab, ylab, clab, title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, color=color, alpha=alpha)) + 
    ggtitle(title) + 
    geom_beeswarm() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    labs(x=xlab, y=ylab, color=clab)
    if (is.null(colorscheme)){
      p <- p + scale_fill_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_fill_manual(values=colorscheme)
    }
  return(p)
}

randombee_plot <- function(data, x, y, color, alpha=1, box=FALSE, median=FALSE, xlab, ylab, clab, title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, color=color)) + 
    ggtitle(title) + 
    geom_quasirandom(dodge.width=0.8, alpha=alpha) + 
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    theme_bw() + 
    labs(x=xlab, y=ylab, color=clab)
    if (is.null(colorscheme)){
      p <- p + scale_color_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_color_manual(values=colorscheme)
    }
    if (box){
      p <- p + stat_boxplot(geom='boxplot', alpha=0.5, outlier.shape=NA, position='dodge', width=0.8)
    }
    if (median){
      p <- p + stat_summary(fun.y='median', fun.ymin='median', fun.ymax='median', geom='crossbar', position='dodge', width=0.8)
    }
  return(p)
}

bar_plot_dodge <- function(data, x, y, fill=NA, xlab='x', ylab='y', flab='fill', title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, fill=fill)) + 
    geom_bar(stat = 'summary', fun.y = "mean", position = 'dodge') + 
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    labs(x=xlab, y=ylab, fill=flab)
    if (is.null(colorscheme)){
      p <- p + scale_fill_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_fill_manual(values=colorscheme)
    }
  return(p)
}

bar_plot_ident <- function(data, x, y, fill, xlab='x', ylab='y', title='', colorscheme=NULL, flegend=TRUE, alegend=FALSE, alpha=NA){
  p <- ggplot(data=data, aes(x=x, y=y, fill=fill, alpha=alpha)) + 
    geom_bar(stat = 'identity') + 
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    scale_alpha_discrete(range = c(1, 0.3)) + 
    labs(x=xlab, y=ylab) + 
    theme(legend.position="none")
    if (is.null(colorscheme)){
      p <- p + scale_fill_viridis(discrete=T, option='viridis', na.value='grey')
    } else {
      p <- p + scale_fill_manual(values=colorscheme, na.value='grey')
    }
    if (!(flegend)){
      p <- p + guides(fill=FALSE)
    }
    if (!(alegend)){
      p <- p + guides(alpha=FALSE)
    }
  return(p)
}

bar_plot_count <- function(data, x, fill=NA, xlab=NA, ylab='count', title='', colorscheme=NULL, flegend=TRUE, alegend=FALSE, alpha=NA){
  p <- ggplot(data=data, aes(x=x, fill=fill)) + 
    geom_bar() + 
    ggtitle(title) +
    theme_bw() + 
    scale_alpha_discrete(range = c(1, 0.3)) + 
    labs(x=xlab, y=ylab) + 
    theme(
      axis.text.x=element_text(angle = 90, hjust = 1)
      )
    if (is.null(colorscheme)){
      p <- p + scale_fill_viridis(discrete=T, option='viridis', na.value='grey')
    } else {
      p <- p + scale_fill_manual(values=colorscheme, na.value='grey')
    }
    if (!(flegend)){
      p <- p + guides(fill=FALSE)
    }
    if (!(alegend)){
      p <- p + guides(alpha=FALSE)
    }
  return(p)
}

violin_plot <- function(data, x, y, fill, alpha=1, xlab='x', ylab='y', flab='fill', legend=TRUE, ylim=c(0, NA), title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, fill=fill), alpha=alpha) + 
    geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) + 
    # stat_summary(fun.y="mean", geom="crossbar", position='dodge', size=1) +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() +
    scale_alpha_discrete(range = c(1, 0.6)) +  
    scale_y_continuous(limits = ylim) +  
    labs(x=xlab, y=ylab, fill=flab)
    if (is.null(colorscheme)){
      p <- p + scale_fill_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_fill_manual(values=colorscheme)
    }
    if (!(legend)){
      p <- p + theme(legend.position="none")
    }
  return(p)
}

box_plot <- function(data, x, y, fill, alpha=NA, xlab='x', ylab='y', flab='fill', flegend=TRUE, alegend=FALSE, out_size=1.5, ylim=c(0, NA), title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, fill=fill, alpha=alpha)) + 
    geom_boxplot(outlier.size=out_size) + 
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() +
    scale_alpha_discrete(range = c(1, 0.6)) +  
    scale_y_continuous(limits = ylim) +  
    labs(x=xlab, y=ylab, fill=flab)
    if (is.null(colorscheme)){
      p <- p + scale_fill_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_fill_manual(values=colorscheme)
    }
    if (!(flegend)){
      p <- p + guides(fill=FALSE)
    }
    if (!(alegend)){
      p <- p + guides(alpha=FALSE)
    }
  return(p)
}

logbox_plot <- function(data, x, y, fill, alpha=1, xlab='x', ylab='y', flab='fill', flegend=TRUE, alegend=FALSE, out_size=1.5, ylim=c(NA, 1), blank_x=FALSE, title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, fill=fill, alpha=alpha)) + 
    geom_boxplot(outlier.size=out_size) + 
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() +
    scale_alpha_discrete(range = c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3)) +  
    scale_y_log10(limits = ylim) +
    labs(x=xlab, y=ylab, fill=flab)
    if (is.null(colorscheme)){
      p <- p + scale_fill_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_fill_manual(values=colorscheme)
    }
    if (!(flegend)){
      p <- p + guides(fill=FALSE)
    }
    if (!(alegend)){
      p <- p + guides(alpha=FALSE)
    }
    if (blank_x){
      p <- p + theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
    }
  return(p)
}

loghex_plot <- function(data, x, y, color='#424949', xlab='x', ylab='y', clab='color', slab='shape', title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y), shape=shape) + 
    geom_hex() + 
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_x_log10() +
    theme_bw() +  
    labs(x=xlab, y=ylab, color=clab)
    if (is.null(colorscheme)){
      p <- p + scale_fill_viridis(option='viridis')
    } else {
      p <- p + scale_fill_manual(values=colorscheme)
    }
  return(p)  
}

qq_plot <- function(data, sample, shape='a', col=sample, xlab='x', ylab='y',slab='shape', clab='color', title='', colorscheme=NULL){
  p <-ggplot(data=data, aes(col=col, shape=shape)) + 
    ggtitle(title) +
    geom_qq(aes(sample = genome_fraction), distribution = stats::qnorm) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    labs(x=xlab, y=ylab, col=clab)
    if (is.null(colorscheme)){
      p <- p + scale_color_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_color_manual(values=colorscheme)
    }
  return(p)
}

qq_plot_ss <- function(x=c(), y=c(), size=1, variable, split, ylab='y', xlab='x', clab='c', legend_title='Legend', legend_labels=c(), dist, title='', colorscheme=NULL){
  x_sep <- x %>% select_(variable, split) %>% group_by_(split) %>% do(data_frame(quantile(.[[variable]], dist))) 
  y_sep <- y %>% select_(variable, split) %>% group_by_(split) %>% do(data_frame(quantile(.[[variable]], dist)))
  xy_tbl <- bind_cols(x_sep, y_sep)[,c(1, 2, 4)] %>% print()
  names(xy_tbl) <- c('c', 'x', 'y')
 
  p <- ggplot(xy_tbl) +
    geom_point(aes(x = x, y = y, color = c), size=size) + 
    ggtitle(title) +
    theme_bw() +
    geom_abline(slope=1, linetype=2) +
    labs(x = xlab, y = ylab, color = clab)
    if (is.null(colorscheme)){
      p <- p + scale_color_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_color_manual(values=colorscheme)
    }
  return(p)
}

ecdf_plot <- function(data, x, line='1', col='#424949', xlab='x', ylab='y',llab='line', clab='color', title='', colorscheme=NULL){
  p <-ggplot(data=data, aes(x=x, col=col), linetype=line) + 
    ggtitle(title) +
    stat_ecdf(size=1) + 
    theme_bw() + 
    labs(x=ylab, y=xlab, col=clab, line=llab) + 
    coord_flip()
    if (is.null(colorscheme)){
      p <- p + scale_color_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_color_manual(values=colorscheme)
    }
  return(p)
}

dot_plot <- function(data, x, y, color, xlab='x', ylab='y', clab='color', alpha=1, size=2, title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, color=color)) + 
    geom_point(size=size, alpha=alpha) + 
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
    theme_bw() +  
    labs(x=xlab, y=ylab, color=clab)
    if (is.null(colorscheme)){
      p <- p + scale_color_viridis(discrete=T, option='viridis', na.value='grey')
    } else {
      p <- p + scale_color_manual(values=colorscheme, na.value='grey')
    }
  return(p)  
}

dot_plot_smooth <- function(data, x, y, color=NULL, xlab='x', ylab='y', clab='color', alpha=1, se=FALSE, span=0.8, size=2, title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, color=color)) + 
    geom_point(size=size, alpha=alpha) + 
    geom_smooth(span=span, se=se) +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
    theme_bw() +  
    labs(x=xlab, y=ylab, color=clab)
    if (is.null(colorscheme)){
      p <- p + scale_color_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_color_manual(values=colorscheme)
    }
  return(p)  
}

log_plot <- function(data, x, y, color='#424949', shape='x', xlab='x', ylab='y', clab='color', slab='shape', alpha=0.8, size=2, title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, color=color), shape=shape) + 
    geom_point(size=size, alpha=alpha) + 
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_x_log10() +
    theme_bw() +  
    labs(x=xlab, y=ylab, color=clab)
    if (is.null(colorscheme)){
      p <- p + scale_color_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_color_manual(values=colorscheme)
    }
  return(p)  
}


loglog_plot <- function(data, x, y, color, alpha=1, xlab, ylab, clab, colorscheme=NULL) {
  p <- ggplot(data=data, aes(x=x, y=y, color=color)) + 
    geom_point(size=1, alpha=alpha) + 
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    scale_y_log10(breaks=c(1, 0.1, 0.01, 0.001, 0.0001, 1e-04, 1e-05), limits=c(1e-05, 1)) +
    scale_x_log10(breaks=c(1, 0.1, 0.01, 0.001, 0.0001, 1e-04, 1e-05), limits=c(1e-05, 1)) +
    labs(x=xlab, y=ylab, color=clab)
    if (is.null(colorscheme)){
      p <- p + scale_color_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_color_manual(values=colorscheme)
    }
  return(p)  
}

gather_matrix <- function(data, key, value){
  reshape2:::melt.matrix(data, varnames=key, value.name=value)
}

heatmap_plot <- function(data, x, y, fill, xlab='x', ylab='y', flab='fill', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, fill=fill)) +
    geom_tile(color='white') + 
    theme_bw()
    if (is.null(colorscheme)){
      p <- p + scale_fill_viridis(option='magma')
    } else {
      p <- p + scale_fill_gradiant(low=colorscheme[1], high=colorscheme[2])
    }
  return(p)
}