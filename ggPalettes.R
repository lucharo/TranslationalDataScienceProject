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

paperPalette = function(theme = theme_gdocs()){
  theme_set(theme)
  
  
  
}