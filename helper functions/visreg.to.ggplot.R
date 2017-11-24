################################################################################
# Introduction
#
# Transform data from visreg::visreg() to be compatible with ggplot2.
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
visreg.to.ggplot <- function(model, type = "contrast", levels = NULL, var_excl = NULL){

  # model       model from which response profiles can be derived.
  # type        type of plot to be produced. see ?visreg
  # levels      optional; levels which should be emphasized with different colours.
  # var_excl    variables to be excluded from ploted response profiles. Usually random factors.

  require(visreg)

  vis_obj <- visreg(model, type = type, plot = F)
  gg_vis_obj3 <- NULL

  if(is.null(levels)){
    levels <- rep(1, dim(model@frame)[1])
  } else{
    levels <- model@frame[, which(colnames(model@frame)==levels)]
  }

  for(i in 1:length(vis_obj)){
    if(dim(vis_obj[[i]]$res)[1]!=dim(model@frame)[1]){
      warning("there is a discrepanicy in cases between the model and the visreg visualisation attempt")
      next
    }
    if(vis_obj[[i]]$meta$x%in%var_excl) {
      next
    } else{
      gg_vis_obj <- data.frame(variable = vis_obj[[i]]$meta$x,
                               type = "fit",
                               x = as.numeric(vis_obj[[i]]$fit[which(colnames(vis_obj[[i]]$fit)==vis_obj[[i]]$meta$x)][[1]]),
                               y = vis_obj[[i]]$fit$visregFit,
                               groups = NA,
                               Upr = vis_obj[[i]]$fit$visregLwr,
                               Lwr = vis_obj[[i]]$fit$visregUpr
      )

      gg_vis_obj2 <- data.frame(variable = vis_obj[[i]]$meta$x,
                                type = "res",
                                x = as.numeric(vis_obj[[i]]$res[which(colnames(vis_obj[[i]]$res)==vis_obj[[i]]$meta$x)][[1]]),
                                y = vis_obj[[i]]$res$visregRes,
                                groups = levels,
                                Upr = NA,
                                Lwr = NA
      )

      if(is.null(gg_vis_obj3)){
        gg_vis_obj3 <- rbind(gg_vis_obj, gg_vis_obj2)
      } else{
        gg_vis_obj3 <- rbind(gg_vis_obj3, gg_vis_obj, gg_vis_obj2)
      }
    }
  }
  return(gg_vis_obj3)
}
