################################################################################
# Introduction
#
# Determine the semi-residuals from a given model.
#
################################################################################
# License
#
# Copyright (C) 2016 Steffen Ehrmann
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful , but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
# function body
semi.residual <- function(model, predictors, f = 3, levels){

  # model       model from which residuals can be derived.
  # predictors  list of variables comprised in the data-frame based on which the model has been
  #             fitted.
  # f           number of standard deviation equivalents to define extreme values.
  # levels      levels which should be emphasized with different colours.

  df <- eval(model@call[[3]])
  n_sub <- subset(df, rownames(df)%in%names(predict(model)))
  levels_name <- levels
  levels <- eval(parse(text=paste("n_sub$", levels, sep="")))
  response <- getME(model, "y")
  response_name <- colnames(model@frame)[1]
  resid <- resid(model)

  for(i in 1:length(predictors)){

    var_ext <- subset(n_sub, xtreme(eval(parse(text=paste("n_sub$", predictors[i], sep=""))), f=f))
    resid_ext <- subset(n_sub, xtreme(resid, f=f))
    names_ext <- c(rownames(var_ext), rownames(resid_ext))
    names_ext <- names_ext[!duplicated(names_ext)]
    extremes_x <- subset(n_sub, rownames(n_sub)%in%names_ext)
    not_extremes_x <- !rownames(n_sub)%in%names_ext
    extremes_y <- subset(resid, names(resid)%in%names_ext)

    x1 <- eval(parse(text=paste("n_sub$", predictors[i], sep="")))
    x2 <- subset(x1, not_extremes_x)
    resid1 <- resid
    resid2 <- subset(resid, not_extremes_x)
    levels2 <- subset(levels, not_extremes_x)

    g1 <- ggplot(dat = n_sub, aes(x = x1, y = resid1)) +
      geom_point(size = 4.2) +
      geom_point(aes(colour = levels), size = 3.5) +
      geom_smooth(method = "lm", colour = "dodgerblue3", formula = y ~ poly(x, 2), linetype = "dashed", se = F) +
      geom_smooth(method = "lm", colour = "firebrick3", linetype = "dashed", se = F) +
      geom_point(data = subset(n_sub, not_extremes_x), aes(x = x2, y = resid2, colour = levels2), size = 3.5) +
      geom_smooth(data = subset(n_sub, not_extremes_x), aes(x = x2, y = resid2), method = "lm", colour = "dodgerblue3", formula = y ~ poly(x, 2), linetype = "solid") +
      geom_smooth(data = subset(n_sub, not_extremes_x), aes(x = x2, y = resid2), method = "lm", colour = "firebrick3", linetype = "solid") +
      geom_point(dat = extremes_x, aes(x = eval(parse(text=paste("extremes_x$", predictors[i], sep=""))), y = extremes_y)) +
      xlab(predictors[i])
    g2 <- ggplot(dat = n_sub, aes(x = eval(parse(text=paste("n_sub$", predictors[i], sep=""))), y = response)) +
      geom_point(size = 4.2) +
      geom_point(aes(colour = levels), size = 3.5) +
      geom_smooth() +
      geom_smooth(method = lm, colour = "red") +
      xlab(predictors[i]) + ylab(response_name)
    mylegend<-g_legend(g1)
    grid.arrange(arrangeGrob(g1 + theme(legend.position="none"), g2 + theme(legend.position="none"), nrow=2), mylegend, ncol=2, widths = c(5, 1))
  }
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}