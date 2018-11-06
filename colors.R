
######################################################
##### ACCEPTABLE COLORS FOR PEOPLE WITH EYEBALLS #####
######################################################

require(viridis)
theme_set(theme_bw())

redpink <- "#c2185b"
purple <- "#7b1fa2"
orange <- "#e65100"
yellow <- "#fbc02d"
green <- "#388e3c"
red <- "#b71c1c" 

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
many <- c("#Ae0031", "#ffbb00", "#2d62a3", "#198749", "#573794","#Ff7708", "#15889C",
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
cl_types_colors <- c("#Ae0031", "#ffbb00", "#198749", "#573794","#Ff7708", "#15889C",
    "#ED49B1", "#aeb6bf", "#5d6d7e")
scale_fill_relaxing <- scale_fill_gradient(low = "#fffbd5", high = "#b20a2c")
scale_color_relaxing <- scale_color_gradient(low = "#fffbd5", high = "#b20a2c")


scale_colour_discrete <- function(...) scale_colour_manual(values=many)
scale_colour_continuous <- function(...) scale_color_distiller(palette="Blues")
scale_fill_discrete <- function(...) scale_fill_manual(values=many)
scale_fill_continuous <- function(...) scale_fill_distiller(palette="Blues")