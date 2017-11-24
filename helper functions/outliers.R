################################################################################
# Introduction
#
# Determine the outliers of a given model based on the number of
# standard-deviation equivalents.
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
outliers <- function(model, predictor, f = 3){

  df <- eval(model@call[[3]])
  n_sub <- subset(df, rownames(df)%in%names(predict(model)))
  resid <- resid(model)

  var_ext <- subset(n_sub, xtreme(eval(parse(text=paste("n_sub$", predictor, sep=""))), f=f))
  resid_ext <- subset(n_sub, xtreme(resid, f=f))
  names_ext <- c(rownames(var_ext), rownames(resid_ext))
  names_ext <- names_ext[!duplicated(names_ext)]
  extremes_x <- subset(n_sub, rownames(n_sub)%in%names_ext)
  extremes_y <- subset(resid, names(resid)%in%names_ext)
  extremes <- list(x = extremes_x, y = extremes_y)
  return(extremes)

}
