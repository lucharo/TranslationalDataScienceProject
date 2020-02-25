############################################################################
############################################################################
###                                                                      ###
###                           GGPLOT TEMPLATES                           ###
###                                                                      ###
############################################################################
############################################################################

######################################################################
##            This script contains two functions one to             ##
##  set up the palette for papers and the other for presentations,  ##
##   these can be called at any instance when plotting is desires   ##
######################################################################

library(ggplot2)
library(GGally)
if (!require(ggthemes)) install.packages("ggthemes")
library(ggthemes)

# check ggthemes at: 
# https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/
# https://cran.r-project.org/web/packages/ggthemes/ggthemes.pdf
# https://www.datanovia.com/en/blog/ggplot-themes-gallery/
# https://rpubs.com/Koundy/71792

paperPalette = function(theme.default = theme_gdocs(base_size = 10)){
  theme.toSet = theme.default + 
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.2, "cm"),
          legend.spacing = unit(0, "cm"),
          legend.title = element_text(face="italic"),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold"))
          
          
    theme_set(theme.toSet)
}

presentationPalette = function(theme.default = theme_gdocs(base_size = 20)){
  theme.toSet = theme.default + 
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.2, "cm"),
          legend.spacing = unit(0, "cm"),
          legend.title = element_text(face="italic"),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold"))
  
  
  theme_set(theme.toSet)
  
  
}

