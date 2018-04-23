library(viridis)
library(tidyverse)
library(ggbeeswarm)
library(reshape2)
library(ggdendro)
library(ggridges)
library(ggfortify)
library(gridExtra)
library(grid)


##################
## COLORSCHEMES ##
##################

magma <- c("#424949", "#F9A825", "#CD3F0A", "#8A023C")
redtogreen <- c("#550527", "#A10702", "#F44708", "#F9A113", "#599124")
greentored <- c("#599124" , "#F9A113", "#F44708", "#A10702", "#550527")
material <- c("#15889C", "#ED496F", "#8E1382", "#C62828", "#FFB300", "#43A047", "#FF6F00")
vibrant <- c("#95190C", "#FFB300", "#1F7530", "#086788", "#6300B5", "#EF8737", "#550527")
blass <- c("#e67e22", "#c0392b", "#6c3483", "#2471a3", "#229954", "#d4ac0d")
magenta <- c("#D96DED", "#6300B5", "#8E1382", "#ED496F", "#FF9D7C")
bluelime <- c("#012824", "#15889C", "#7BCFC0", "#88E200", "#028F03")
bricksky <- c("#0B3954", "#28AAE1", "#FDC120", "#D94551", "#55545C")
colorful <- c("#9E0031", "#FFBB00", "#198749", "#2D62A3", "#8E1382", "#5A0002")
flat <- c("#471b9a" , "#ff6f00", "#1565c0", "#43a047", "#ffb300", "#c62828")
three <- c("#16A085", "#2980B9", "#8E44AD")
many <- c("#Ae0031", "#ffbb00", "#198749", "#2d62a3", "#573794","#Ff7708", "#15889C", 
    "#ED49B1", "#8E13A2", "#ead9d5", "#666686", "#CD6155", "#F7DC6F", "#7DCEA0", 
    "#85C1E9", "#EB984E", "#03cd4A", "#CDDC39", "#e0594b", "#c76cde", "#24B177", 
    "#8D6E63", "#486ff7", "#6300b5", "#88e200", "#012824", "#0d3290", "#ead9d5", 
    "#a347fb", "#54fc7a", "#eb1388", "#b0978d", "#fe52cf", "#83f1f6", "#f1f847", 
    "#2b1dfc", "#6c6f15", "#6ca05c", "#7788cd", "#f502f3", "#0dc290", "#fa0e03", 
    "#3caa0a", "#befc8d", "#08f8eb", "#b1cd3f", "#d6a5fa", "#ce606c", "#ab1eba", 
    "#6ecc9f", "#054ddc", "#486ff7", "#854f49", "#f22B21", "#3a0e43", "#225805", 
    "#37d160", "#e4b974", "#a8bade", "#47edd1", "#f47a92", "#c76cde", "#9106eb", 
    "#81aa20", "#d7fdfd", "#5deb2e", "#f82745", "#6435e0", "#027ffe", "#8e3101", 
    "#16f648", "#1c15bc", "#8be46e", "#8d6fa0", "#e68fc6", "#058ca9", "#9e018a", 
    "#bdfd0b", "#b22760", "#2bf49f", "#cb9348", "#9d8303", "#c251a1", "#46adaf", 
    "#a3e3af", "#22bb34", "#6ea3fa", "#260374", "#1c3854", "#405d37", "#c21df3", 
    "#fcea92", "#537f88", "#fd4c18", "#f2d71e", "#fd4c7a")
cluster_cols <- c("#Ae0031", "#ffbb00", "#4CAF50", 
    "#2d62a3", "#573794",
    "#Ff7708", "#15889C", "#ED49B1", "#8E13A2", 
    "#1B5E20", 
    "#626567", 
    "#f47a92", "#F7DC6F", 
    "#A6ACAF")




#######################
## UTILITY FUNCTIONS ##
#######################

insert_minor <- function(major_labs, n_minor) {
  labs <- c(sapply( major_labs, function(x) c(x, rep("", n_minor))))
  labs[1:(length(labs)-n_minor)]
}

color_gradient <- function(start="white", end="#550527", n=100){
  return(colorRampPalette(c(start, end))(n))
}

color_gradient_scrambled <- function(start="white", end="#550527", n=100){
  return(sample(colorRampPalette(c(start, end))(n)))
}

qq_ss_plot <- function(x=c(), y=c(), size=1, variable, split, legend_title="Legend", 
    legend_labels=c(), dist, colorscheme=NULL)
{
  x_sep <- x %>% 
    select_(variable, split) %>% 
    group_by_(split) %>% 
    do(data_frame(quantile(.[[variable]], dist))) 
  y_sep <- y %>% 
    select_(variable, split) %>% 
    group_by_(split) %>% 
    do(data_frame(quantile(.[[variable]], dist)))
  xy_tbl <- bind_cols(x_sep, y_sep)[,c(1, 2, 4)] %>% print()
  names(xy_tbl) <- c("c", "x", "y")
 
  p <- ggplot(xy_tbl) +
    geom_point(aes(x = x, y = y, color = c), size=size) + 
    theme_bw() +
    geom_abline(slope=1, linetype=2)
    if (is.null(colorscheme)){
      p <- p + scale_color_viridis(discrete=T, option="viridis")
    } else {
      p <- p + scale_color_manual(values=colorscheme)
    }
  return(p)
}

gather_matrix <- function(data, key, value){
  reshape2:::melt.matrix(data, varnames=key, value.name=value)
}

heatmap_plot <- function(data, x, y, fill=NA, colorscheme=NULL)
{
  p <- ggplot(data=data, aes(x=x, y=y, fill=fill)) +
    geom_tile(color="white") + 
    theme_bw()
    if (is.null(colorscheme)){
      p <- p + scale_fill_viridis(option="magma")
    } else {
      p <- p + scale_fill_gradiant(low=colorscheme[1], high=colorscheme[2])
    }
  return(p)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot = function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots = c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout = matrix(seq(1, cols * ceiling(numPlots/cols)),
            ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots == 1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                layout.pos.col = matchidx$col))
        }
    }
}

multiplot_shared_legend <- function(..., plotlist = NULL, nrow = 1, ncol = length(c(list(...), plotlist)), position = c("bottom", "right")) {
    
    plots <- c(list(...), plotlist)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
    gl <- c(gl, nrow = nrow, ncol = ncol)
    
    combined <- switch(position,
        "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
            legend,
            ncol = 1,
            heights = unit.c(unit(1, "npc") - lheight, lheight)),
        "right" = arrangeGrob(do.call(arrangeGrob, gl),
            legend,
            ncol = 2,
            widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
    grid.newpage()
    grid.draw(combined)
    
}
