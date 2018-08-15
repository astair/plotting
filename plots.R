
##########################################################
## UTILITY FUNCTIONS FOR PLOTTING AND OTHER SHENANIGANS ##
##########################################################

require(viridis)
require(tidyverse)
require(ggbeeswarm)
require(reshape2)
require(ggdendro)
require(ggridges)
require(ggfortify)
require(gridExtra)
require(grid)

theme_set(theme_bw())

as_matrix <- function(df){
    mat <- as.matrix(df[2:ncol(df)])
    rownames(mat) <- df[[1]]
    return(mat)
}

write_matrix <- function(x, path, delim="\t", rownames="X1"){
    as_tibble(as.matrix(x), rownames=rownames) %>%
        write_delim(path=path, delim=delim)
}

read_matrix <- function(path, delim="\t"){
    tbl <- read_delim(path, delim="\t")
    mat <- as.matrix(tbl[,2:length(colnames(tbl))])
    rownames(mat) <- tbl[[1]]
    return(mat)
}

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

extract_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
