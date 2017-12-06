library(viridis)
library(tidyverse)
library(ggbeeswarm)

magma <- c('#424949', '#F9A825', '#CD3F0A', '#8A023C')
redtogreen <- c('#550527', '#A10702', '#F44708', '#F9A113', '#599124')
greentored <- c('#599124' , '#F9A113', '#F44708', '#A10702', '#550527')
material <- c('#15889C', '#ED496F', '#8E1382', '#C62828', '#FFB300', '#FF6F00', '#43A047')
vibrant <- c('#95190C', '#FFB300', '#1F7530', '#086788', '#6300B5', '#EF8737', '#550527')
blass <- c('#e67e22', '#c0392b', '#6c3483', '#2471a3', '#229954', '#d4ac0d')
magenta <- c('#D96DED', '#6300B5', '#8E1382', '#ED496F', '#FF9D7C')
bluelime <- c('#012824', '#15889C', '#7BCFC0', '#88E200', '#028F03')
bricksky <- c('#0B3954', '#28AAE1', '#FDC120', '#D94551', '#55545C')
magenta <- c('#D96DED', '#6300B5', '#8E1382', '#ED496F', '#FF9D7C')
colorful <- c('#9E0031', '#FFBB00', '#198749', '#2D62A3', '#8E1382', '#5A0002')
flat <- c('#471b9a' , '#ff6f00', '#1565c0', '#43a047', '#ffb300', '#c62828')

many <- c('#9e0031', '#ffbb00', '#198749', '#2d62a3', '#5a0002', '#e67e22', '#c0392b', '#6c3483', '#247177', '#d4ac0d', '#28aae1', '#6300b5', '#88e200', '#012824', '#0d3290', '#ead9d5', '#a347fb', '#54fc7a', '#eb1388', '#03dd14', '#b0978d', '#fe52cf', '#83f1f6', '#f1f847', '#2b1dfc', '#6c6f15', '#6ca05c', '#7788cd', '#f502f3', '#0dc290', '#fa0e03', '#3caa0a', '#befc8d', '#08f8eb', '#b1cd3f', '#d6a5fa', '#ce606c', '#ab1eba', '#6ecc9f', '#054ddc', '#486ff7', '#854f49', '#3a0e43', '#225805', '#37d160', '#e4b974', '#a8bade', '#47edd1', '#f47a92', '#c76cde', '#9106eb', '#81aa20', '#d7fdfd', '#5deb2e', '#f82745', '#6435e0', '#027ffe', '#8e3101', '#16f648', '#1c15bc', '#8be46e', '#8d6fa0', '#e68fc6', '#058ca9', '#9e018a', '#bdfd0b', '#b22760', '#2bf49f', '#cb9348', '#9d8303', '#c251a1', '#46adaf', '#770764', '#a3e3af', '#22bb34', '#6ea3fa', '#260374', '#1c3854', '#405d37', '#c21df3', '#fcea92', '#537f88', '#fd4c18', '#f2d71e', '#fd4c7a')

barPlot <- function(data, x, y, fill, xlab, ylab, flab, title='', colorscheme=NULL){
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

beeswarmPlot <- function(data, x, y, color, alpha=1, xlab, ylab, clab, title='', colorscheme=NULL){
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

randombeePlot <- function(data, x, y, color, alpha=1, box=FALSE, xlab, ylab, clab, title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, color=color)) + 
    ggtitle(title) + 
    geom_quasirandom(dodge.width=0.8, alpha=alpha) + 
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    theme_bw() + 
    labs(x=xlab, y=ylab, color=clab)
    if (is.null(colorscheme)){
      p <- p + scale_fill_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_fill_manual(values=colorscheme)
    }
    if (box){
      p <- p + stat_boxplot(geom='boxplot', alpha=0.5, outlier.shape=NA, position='dodge', width=0.8)
    }
  return(p)
}

barPlotDodge <- function(data, x, y, fill=NA, xlab, ylab, flab=NA, title='', colorscheme=NULL){
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

barPlotCount <- function(data, x, y, fill, xlab, ylab, title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, fill=fill)) + 
    geom_bar(stat = 'identity') + 
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    labs(x=xlab, y=ylab) + 
    theme(legend.position="none")
    if (is.null(colorscheme)){
      p <- p + scale_fill_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_fill_manual(values=colorscheme)
    }
  return(p)
}

boxPlot <- function(data, x, y, fill, xlab='x', ylab='y', flab='fill', title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, fill=fill)) + 
    geom_boxplot() + 
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

logboxPlot <- function(data, x, y, fill, xlab='x', ylab='y', flab='fill', title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, fill=fill)) + 
    geom_boxplot() + 
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    scale_y_log10() +
    labs(x=xlab, y=ylab, fill=flab)
    if (is.null(colorscheme)){
      p <- p + scale_fill_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_fill_manual(values=colorscheme)
    }
  return(p)
}

qqPlot <- function(data, sample, shape='a', col=sample, xlab='x', ylab='y',slab='shape', clab='color', title='', colorscheme=NULL){
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

qqPlotSS <- function(x=c(), y=c(), colors=c(), variable, split, ylab='y', xlab='x', clab='c', legend_title='Legend', legend_labels=c(), dist, title='', colorscheme=NULL){
  x_sep <- x %>% select_(variable, split) %>% group_by_(split) %>% do(data_frame(quantile(.[[variable]], dist))) 
  y_sep <- y %>% select_(variable, split) %>% group_by_(split) %>% do(data_frame(quantile(.[[variable]], dist)))
  xy_tbl <- bind_cols(x_sep, y_sep)[,c(1, 2, 4)] %>% print()
  names(xy_tbl) <- c('c', 'x', 'y')
 
  p <- ggplot(xy_tbl) +
    geom_point(aes(x = x, y = y, color = c)) + 
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

ecdfPlot <- function(data, x, line='1', col='#424949', xlab='x', ylab='y',llab='line', clab='color', title='', colorscheme=NULL){
  p <-ggplot(data=data, aes(x=x, col=col, linetype=line)) + 
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

dotPlot <- function(data, x, y, color, xlab='x', ylab='y', clab='color', alpha=1, size=2, title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, color=color)) + 
    geom_point(size=size, alpha=alpha) + 
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

dotPlotSmooth <- function(data, x, y, color=NULL, xlab='x', ylab='y', clab='color', alpha=1, se=FALSE, span=0.8, size=2, title='', colorscheme=NULL){
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

logPlot <- function(data, x, y, color='#424949', shape='x', xlab='x', ylab='y', clab='color', slab='shape', alpha=0.8, size=2, title='', colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, color=color, shape=shape)) + 
    geom_point(size=size, alpha=alpha) + 
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_x_log10() +
    theme_bw() +  
    labs(x=xlab, y=ylab, color=clab, shape=slab)
    if (is.null(colorscheme)){
      p <- p + scale_color_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_color_manual(values=colorscheme)
    }
  return(p)  
}

loglogPlot <- function(data, x, y, color, xlab, ylab, clab, title='', xlims=c(-5, 0), ylims=c(-5, 0), colorscheme=NULL){
  p <- ggplot(data=data, aes(x=x, y=y, color=color)) + 
    geom_point(size=2, alpha=1) + 
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    scale_y_log10() +
    scale_x_log10() +
    scale_x_continuous(limits = xlims) + 
    scale_y_continuous(limits = ylims) +  
    labs(x=xlab, y=ylab, color=clab)
    if (is.null(colorscheme)){
      p <- p + scale_color_viridis(discrete=T, option='viridis')
    } else {
      p <- p + scale_color_manual(values=colorscheme)
    }
  return(p)  
}