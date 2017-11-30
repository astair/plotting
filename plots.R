library(viridis)

magma <- c('#424949', '#F9A825', '#CD3F0A', '#8A023C')

barPlot <- function(data, x, y, fill, xlab, ylab){
  p <- ggplot(data=data, aes(x=x, y=y, fill=fill)) + 
    geom_bar(stat = 'summary', fun.y = "mean") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    scale_fill_viridis(discrete=T, option='viridis') +
    labs(x=xlab, y=ylab) + 
    theme(legend.position="none")
  return(p)
}

barPlotDodge <- function(data, x, y, fill=NA, xlab, ylab, flab=NA){
  p <- ggplot(data=data, aes(x=x, y=y, fill=fill)) + 
    geom_bar(stat = 'summary', fun.y = "mean", position = 'dodge') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    scale_fill_viridis(discrete=T, option='viridis') +
    labs(x=xlab, y=ylab, fill=flab)
  return(p)
}

barPlotCount <- function(data, x, y, fill, xlab, ylab){
  p <- ggplot(data=data, aes(x=x, y=y, fill=fill)) + 
    geom_bar(stat = 'identity') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    scale_fill_viridis(discrete=T, option='viridis') +
    labs(x=xlab, y=ylab) + 
    theme(legend.position="none")
  return(p)
}

boxPlot <- function(data, x, y, fill, xlab='x', ylab='y', flab='fill'){
  p <- ggplot(data=data, aes(x=x, y=log10(y), fill=fill)) + 
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    scale_fill_viridis(discrete=T, option='viridis') +
    labs(x=xlab, y=ylab, fill=flab)
  return(p)
}

boxPlot2 <- function(data, x, y, fill, xlab='x', ylab='y', flab='fill'){
  p <- ggplot(data=data, aes(x=x, y=log10(y), fill=fill)) + 
    geom_boxplot(notch=F) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    # scale_fill_manual(values=magma[2:4])
    scale_fill_viridis(discrete=T, option='viridis') +
    labs(x=xlab, y=ylab, fill=flab)
  return(p)
}

qqPlot <- function(data, sample, shape='a', col=sample, xlab='x', ylab='y',slab='shape', clab='color'){
  p <-ggplot(data=data, aes(col=col, shape=shape)) + 
    geom_qq(aes(sample = genome_fraction), distribution = stats::qnorm) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    scale_color_viridis(discrete=T, option='viridis') +
    labs(x=xlab, y=ylab, col=clab)
  return(p)
}

qqPlotSS <- function(x=c(), y=c(), colors=c(), variable, split, ylab='y', xlab='x', legend_title='Legend', legend_labels=c(), dist){
  x_sep <- x %>% select_(variable, split) %>% group_by_(split) %>% do(data_frame(quantile(.[[variable]], dist))) 
  y_sep <- y %>% select_(variable, split) %>% group_by_(split) %>% do(data_frame(quantile(.[[variable]], dist)))
  xy_tbl <- bind_cols(x_sep, y_sep)[,c(1, 2, 4)] %>% print()
  names(xy_tbl) <- c('c', 'x', 'y')
 
  p <- ggplot(xy_tbl) +
    geom_point(aes(x = x, y = y, color = c)) + 
    theme_bw() +
    geom_abline(slope=1, linetype=2) +
    # scale_color_manual(values=colors) + 
    scale_color_viridis(discrete=T, option='viridis') +
    labs(x = xlab, y = ylab)
  return(p)
}

ecdfPlot <- function(data, x, line='1', col='#424949', xlab='x', ylab='y',llab='line', clab='color'){
  p <-ggplot(data=data, aes(x=x, col=col, linetype=line)) + 
    stat_ecdf(size=1) + 
    theme_bw() + 
    # scale_color_manual(values=magma[2:4]) +
    scale_color_viridis(discrete=T, option='viridis') +
    labs(x=xlab, y=ylab, col=clab, line=llab) + 
    coord_flip()
  return(p)
}

dotPlot <- function(data, x, y, color, xlab, ylab, clab) {
  p <- ggplot(data=data, aes(x=x, y=y, color=color)) + 
    geom_point(size=2, alpha=0.8) + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
    theme_bw() + 
    scale_color_viridis(discrete=T, option='viridis') + 
    labs(x=xlab, y=ylab, color=clab)
  return(p)  
}

logPlot <- function(data, x, y, color='#424949', shape='x', xlab='x', ylab='y', clab='color', slab='shape', alpha=0.8, size=2) {
  p <- ggplot(data=data, aes(x=log10(x), y=y, color=color, shape=shape)) + 
    geom_point(size=size, alpha=alpha) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    scale_color_viridis(discrete=T, option='viridis') + 
    labs(x=xlab, y=ylab, color=clab, shape=slab)
  return(p)  
}

loglogPlot <- function(data, x, y, color, xlab, ylab, clab) {
  p <- ggplot(data=data, aes(x=log10(x), y=log10(y), color=color)) + 
    geom_point(size=2, alpha=1) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    scale_color_viridis(discrete=T, option='viridis') + 
    labs(x=xlab, y=ylab, color=clab)
  return(p)  
}