################################################################################
# Introduction
#
# Carry out a variance decomposition of a given model. This enables calculation
# of a relative importance value for each predictor.
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
var.decomp <- function(model, ddf = "Satterthwaite", type = 3, estimator = "etasq",
                       semipartial = FALSE, penalised = FALSE, SS_adjust = FALSE,
                       order = NULL, verbose = NULL){

  # model           specify model-object here (currently only merMod).
  # ddf             determine degrees of freedom; see ?anova.merModLmerTest.
  # type            determine type of hypothesis to be tested; see ?anova.merModLmerTest.
  # estimator       estimator for magnitude of effect size, currently either 'etasq' or
  #                 'omegasq'; see 'https://support.sas.com/documentation/cdl/en/statug/63347/HTML/default/viewer.htm#statug_glm_sect032.htm'
  #                 for SAS-approach and 'http://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/effectSize' to get ranges for
  #                 judging these values.
  # semipartial     logical; should the magnitude of each effect size also be calculated with
  #                 respect to the total variance including random effects? default = FALSE.
  # penalised       logical; should penalised (weighted) residual sum of squares be calculate
  #                 and used for inference? default = FALSE.
  # SS_adjust       logical; when using either approximation for calculating denominator
  #                 degrees of freedom, the sum of squares of the fixed variables change, in
  #                 some cases so drastically, that their values do not add up to 100% with the
  #                 residual sum of squares anymore. When setting this argument to TRUE, the
  #                 resulting sum of the fixed variable sum of squares will be scaled to the
  #                 sum of squares remaining after removing residual (and random) sum of
  #                 squares; default = FALSE
  # order           give name of the column of final output to order by and if order should be
  #                 carried out decreasing or not. (e.g. c("Pr(>F)", FALSE)).
  # verbose         create more output by specifying '1' or '2'.

  require(lme4); require(lmerTest)
  ddf <- match.arg(ddf, c("Satterthwaite", "Kenward-Roger", "lme4"))
  estimator <- match.arg(estimator, c("etasq", "omegasq"))

  if(!any(is(model)=="merMod")){
    stop("this function runs only with merMod-objects currently")
  } else{
    estim <- summary(model)$coefficients[,1][-1]
    aov <- anova(model, ddf = ddf, type = type)
    aov$coef <- estim

    N <- nobs(model)
    N_ran <- getME(model, "l_i"); names(N_ran) <- NULL
    MS_f <- aov$'Mean Sq'
    SS_f <- aov$'Sum Sq'
    if(ddf == "lme4"){
      DF_f <- aov$Df
    } else{
      DF_f <- aov$'NumDF'
    }

    Var <- as.data.frame(VarCorr(model))
    Var_res <- Var$vcov[which(Var$grp=="Residual")]
    Var_ran <- Var$vcov[which(Var$grp!="Residual")]

    if(penalised){
      # penalised weighted residual (error + random) sum of squares (pwrss)
      SS_e <- getME(model, name="devcomp")$cmp[5]; names(SS_e) <- NULL
    } else{
      # weighted residual (error + random) sum of squares (wrss)
      SS_e <- getME(model, name="devcomp")$cmp[3]; names(SS_e) <- NULL
    }

    MS_e <- SS_e / df.residual(model)

    SS_total <- sum((model.response(model.frame(model)) - mean(model.response(model.frame(model))))**2)

    SS_reg <- SS_total - SS_e
    # correct all fixed effect sum of squares by what they can maximally take up in a
    # deterministic interpretation (based on really existing sum of squares).
    SS_f_adj <- SS_reg/sum(SS_f) * SS_f

    if(!is.null(verbose)){

      if(!SS_adjust){
        stats <- c(SS_reg, sum(SS_f), SS_e, SS_total)
        stats <- data.frame(stats, stats/SS_total*100)
        colnames(stats) <- c('Sum Sq', 'Fract'); rownames(stats) <- c("Regression", "  Fixed", "Residual", "Total")
      } else{
        stats <- c(SS_reg, sum(SS_f), sum(SS_f_adj), SS_e, SS_total)
        stats <- data.frame(stats, stats/SS_total*100)
        if(ddf == "Satterthwaite"){
          colnames(stats) <- c('Sum Sq', 'Fract'); rownames(stats) <- c("Regression", "  Fix (ST)", "  Fix (adj)", "Residual", "Total")
        } else if(ddf == "Kenward-Roger"){
          colnames(stats) <- c('Sum Sq', 'Fract'); rownames(stats) <- c("Regression", "  Fix (KR)", "  Fix (adj)", "Residual", "Total")
        } else if(ddf == "lme4"){
          colnames(stats) <- c('Sum Sq', 'Fract'); rownames(stats) <- c("Regression", "  Fix (lme4)", "  Fix (adj)", "Residual", "Total")
        }
      }
    } else{
      stats <- c(SS_reg, SS_e, SS_total)
      stats <- data.frame(stats, stats/SS_total*100)
      colnames(stats) <- c('Sum Sq', 'Fract'); rownames(stats) <- c("Regression", "Residual", "Total")
    }

    # this sets all sum of squares so that they correspond to their relative value for the
    # overall model.
    if(SS_adjust){
      SS_f <- SS_f_adj
      aov$'Sum Sq adj' <- SS_f_adj
      if(ddf == "lme4"){
        aov <- aov[c(1, 2, 6, 3:5)]
      } else {
        aov <- aov[c(1, 8, 2:7)]
      }
    }

    stats <- round(stats, 2)

    rsq_fixed <- 1-(SS_e/N)/(SS_total/N)
    rsq_fixed_adj <- 1-MS_e/(SS_total/(N-1))

    # R^2 of a linear model for the response of the mixed model on the fitted values of the
    # same model
    lmfit <-  lm(model.response(model.frame(model)) ~ fitted(model))
    rsq_total <- summary(lmfit)$r.squared

    # "comparing the residual variance of the full model against the residual variance of a
    # (fixed) intercept-only null model"
    omega0sq <- 1-var(residuals(model))/(var(model.response(model.frame(model))))

    # put stuff together
    modelfit <- data.frame(rsq_fixed, rsq_fixed_adj, rsq_total, omega0sq)
    colnames(modelfit) <- c("R^2", "R^2.adj", "R^2 R~F", "omega0^2"); rownames(modelfit) <- c('')
    modelfit <- round(modelfit, 4)

    if(estimator=="etasq"){
      etasq_p <- SS_f/(SS_f + SS_e)
      p_relImp <- etasq_p/sum(etasq_p)*rsq_fixed_adj
      aov$'eta^2_part' <- etasq_p
      if(semipartial){
        etasq <- SS_f/(SS_total)
        relImp <- etasq/sum(etasq)*rsq_fixed_adj
        aov$'eta^2' <- etasq
      }
    } else if(estimator=="omegasq"){
      omegasq_p <- (SS_f - DF_f * MS_e) / (SS_f + (N - DF_f) * MS_e)
      p_relImp <- omegasq_p/sum(omegasq_p)*rsq_fixed_adj
      aov$'omega^2_part' <- omegasq_p
      if(semipartial){
        omegasq <- (SS_f - DF_f * MS_e) / (SS_total + MS_e)
        relImp <- omegasq/sum(omegasq)*rsq_fixed_adj
        aov$'omega^2' <- omegasq
      }
    }

    aov$'relImp_part' <- p_relImp
    if(semipartial){
      aov$' relImp' <- relImp
    }

    if(!is.null(order)){
      aov <- aov[order(aov[, which(colnames(aov)==order[1])], decreasing=order[2]),]
    }

  }
  object <- list()
  object$'Fixed Part' <- aov
  object$'Var. Partitioning' <- stats
  object$'modelfit' <- modelfit

  object$warning <- NULL
  if(any(Var_ran==0)){
    object$warning <- "variance of (one) of the random effects is zero."
  }
  if(round(SS_reg, 2) < round(sum(SS_f), 2)){
    object$warning <- c(object$warning, paste("'Regression sum of squares' are smaller than the sum of 'fixed variable sum of squares' based on ", ddf, "-approximation. Consider setting 'SS_adjust = TRUE'.", sep=""))
  }

  if(!is.null(verbose)){
    if(verbose==1){
      object <- object
    }
    if(verbose==2){
      object$estim <- estim
      object$'VarCorr' <- Var
    }
  } else{
    if(ddf == "lme4"){
      if(semipartial){
        object[[1]] <- object[[1]][-c(1, 7, 9)]
      } else{
        object[[1]] <- object[[1]][-c(1, 7)]
      }
    } else{
      if(semipartial){
        object[[1]] <- object[[1]][-c(4, 5, 9, 11)]
      } else{
        object[[1]] <- object[[1]][-c(4, 5, 9)]
      }
    }
  }

  return(object)

}
