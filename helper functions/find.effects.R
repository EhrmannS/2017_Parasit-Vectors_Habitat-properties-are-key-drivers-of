################################################################################
# Introduction
#
# We want to add a whole bunch of variables to an existing model and want to
# judge which set of variables would be the best subset. Ordination is
# recommended to decrease the amount of work which has to be put into this
# workflow, but since there is a way of finding covariation is given in this
# function, it is not necessary. We hence add each potential variable to the
# existing model and check the variables signal. We are interested in judging
# the # model parsimony (AIC), the significance and the estimate of the
# potential variable and how strong the potential variables covaries with the
# already existing model.
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
find.effects <- function(model=NULL, variables=NULL, degree=1, interaction=NULL, stat="F",
                         type="III", sbst=NULL, order=NULL, verbose=FALSE, print=TRUE){

  # model       model formula of any kind that can be updated via update()
  # variables   character vector of variable names, which can also be found as columns in the
  #             'data'-argument of the model. If a variable is already included in the model,
  #             it will be excluded automatically.
  # degree      give the polynom degrees for each variable that should be tested.
  # interaction give the name of a variable with which interaction terms should be tested.
  #             Interactions can so far only be tested for first order polynomials.
  # stat        specify which test-statistic to use (as to be found in ?Anova).
  # sbst        specify if a subset for one of/a combination of 'estim', 'p_val', 'gof' or
  #             'r_sq' should be specified. (e.g. 'subset="p_val < 0.05 & r_sq < 0.5')
  # order       specify if the resulting data.frame should be sorted by the estimates. 'estim'
  #             for Estimates, 'gof' for REML/AIC, 'p_val' for p-values and r_sq for tolerance.
  # verbose     logical; specify for more or less information. default = FALSE
  # print       logical; print output or not. default = FALSE

  require(car); require(lme4); require(pascal)

  if(is.null(variables)){
    stop("please specify variables which should be tested in this routine.")
  }
  if(is.null(model)){
    stop("please specify a model which should be tested in this routine.")
  }

  g <- get.args()

  var_names <- as.character(attr(attr(model@frame, "term"), "term.labels"))[!as.character(attr(attr(model@frame, "term"), "term.labels")) %in% as.character(attr(attr(model@frame, "term"), "predvars.random"))[-c(1, 2)]]

  if(!any(interaction==var_names) & !is.null(interaction)){
    stop("term to test interactions with is not yet part of the model. Please include this before testing interactions!")
  }

  estim <- NULL; p_val <- NULL; stars <- NULL; names <- NULL; gof <- NULL; r_sq <- NULL; failed <- list()
  terms <- NULL; term_names <- NULL; term_searches <- NULL
  pb <- txtProgressBar(min = 0, max = length(variables), style = 3, char=">", width=getOption("width")-14)

  for(i in 1:length(variables)){

    switch1 <- FALSE
    iterations <- degree+length(interaction)

    for(j in 1:iterations){

      # define variables (term), their "names" (term_name) and what has to be searched for in the model created output (term_search) here
      if(j==1){

        if(any(variables[i]==var_names)){
          if(verbose){
            cat(paste("\n'", variables[i], "' already available\n", sep=""))
          }
          switch1 <- T
          next
        } else {
          if(is.factor(get(paste(model@call[[3]], sep=""))[,eval(variables[i])])){
            term <- variables[i]
            term_search1 <- variables[i]
            term_search2 <- paste(variables[i], "0", sep="")
            term_name <- term
          } else{
            term <- variables[i]
            term_search1 <- variables[i]
            term_search2 <- variables[i]
            term_name <- term
          }
        }

      } else {

        if(j==iterations & !is.null(interaction)){

          if(any(paste(variables[i], ":", interaction, sep="")==var_names)) {
            if(verbose){
              cat(paste("\n'", variables[i], ":", interaction, "' already available\n", sep=""))
            }
            next
          } else if(any(paste(interaction, ":", variables[i], sep="")==var_names)){
            if(verbose){
              cat(paste("\n'", interaction, ":", variables[i], "' already available\n", sep=""))
            }
            next
          }
          if(switch1 & any(interaction==var_names)){
            if(which(var_names == variables[i]) < which(var_names == interaction)){
              term_search1 <- paste(variables[i], ":", interaction, sep="")
              term_search2 <- paste(variables[i], ":", interaction, "0", sep="")
            } else{
              term_search1 <- paste(interaction, ":", variables[i], sep="")
              term_search2 <- paste(interaction, "0", ":", variables[i], sep="")
            }
            term <- term_search1
            term_name <- term_search1
          } else{
            term_search1 <- paste(variables[i], ":", interaction, sep="")
            term_search2 <- paste(variables[i], ":", interaction, "0", sep="")
            term <- paste(variables[i], "*", interaction, sep="")
            term_name <- term_search1
          }

        } else{

          if(any(paste("I(", variables[i], "^", j, ")", sep="")==var_names)){
            if(verbose){
              cat(paste("\n'I(", variables[i], "^", j, ")' already available\n", sep=""))
            }
            next
          } else if(is.factor(get(paste(model@call[[3]], sep=""))[,eval(variables[i])])){
            if(verbose){
              cat(paste("\n'", variables[i], "' is a factor and can't be squared\n", sep=""))
            }
            next
          } else{
            term_search1 <- paste("I(", variables[i], "^", j, ")", sep="")
            term_search2 <- paste("I(", variables[i], "^", j, ")", sep="")
            term_name <- paste(variables[i], " + I(", variables[i], "^", j, ")", sep="")
            if(switch1){
              term <- term_search1
            } else {
              term <- paste(term, " + ", term_search1, sep="")
            }
          }
        }

      }

      # build new models based on the input variables
      updated <- paste(". ~ . + ", term, sep="")
      updated <- tryCatch(suppressMessages(update(model, updated)), error = function(e) cat(paste("\n", term, " created an error:", sep=""))) # still problem with chol.default(Hessian)

      # check for flawed models and exclude them (with a warning)
      if(!is.null(attr(updated@pp$X, "msgRankdrop"))){

        cat("\n'", term, "': ", attr(updated@pp$X, "msgRankdrop"), "\n", sep="")
        p_val <- c(p_val, NA)
        estim <- c(estim, NA)
        r_sq <- c(r_sq, NA)
        gof <- c(gof, NA)
        names <- c(names, NA)
        failed$rank <- c(failed$rank, term)

      } else if(!length(updated@optinfo$conv$lme4)==0){

        cat("\n'", term, "': ", paste(updated@optinfo$conv$lme4$messages, collapse = "; "), "\n", sep="")
        p_val <- c(p_val, NA)
        estim <- c(estim, NA)
        r_sq <- c(r_sq, NA)
        gof <- c(gof, NA)
        names <- c(names, NA)
        failed$conv <- c(failed$conv, term)

      } else {

        # derive p-values
        a <- Anova(updated, type=type, test.statistic = stat)
        df.Anova <- data.frame(p_val=a[,dim(a)[2]], names=rownames(a))
        p_val <- c(p_val, df.Anova$p_val[df.Anova$names==term_search1])

        # estimate and tolerance
        if(j==iterations & !is.null(interaction)){
          estim <- c(estim, NA)
          r_sq <- c(r_sq, 0)
          names <- c(names, NA)
        } else{
          s <- summary(updated)$coef
          df.summary <- data.frame(estim=s[,1], names=rownames(s))
          estim <- c(estim, df.summary$estim[df.summary$names==term_search2])
          variable_fit <- paste(variables[i], " ~ ", paste(var_names, collapse = "+"), sep="")
          variable_fit <- eval(call('lm', formula=variable_fit, data=model@call[[3]]))
          r_sq <- c(r_sq, summary(variable_fit)$r.squared)
          names <- c(names, term_name)
        }

        # goodness of fit
        if(isREML(updated)){
          gof_temp <- getME(updated, "devcomp")$cmp[7]
        } else{
          gof_temp <- AIC(updated)
        }
        gof <- c(gof, gof_temp)

      }
    }

    setTxtProgressBar(pb, i)

  }

  close(pb)

  if(!is.null(sbst)){
    sbst <- eval(parse(text=sbst))
    estim <- subset(estim, sbst)
    p_val <- subset(p_val, sbst)
    gof <- subset(gof, sbst)
    r_sq <- subset(r_sq, sbst)
    names <- subset(names, sbst)
  }


  if(length(names)==0){
    warning("none of the selected variables satisfy the specifications of this call")
  } else {
    df <- data.frame(name = names, Estimate = format(round(estim, 3)), gof = format(round(gof, 0)), rsq = format(round(r_sq, 2)), p_val = format(round(p_val, 5)))
    stars <- symnum(p_val, corr = FALSE,
                    cutpoints = c(0, .001, .01, .05, .1, 1),
                    symbols = c("***", "**", "*", ".", " "))
    df$stars <- format(noquote(stars))

    if(isREML(updated)){
      colnames(df) <- c("name", "Estimate", "REML", "R²", paste("Pr(>", stat, ")", sep=""), " ")
    } else{
      colnames(df) <- c("name", "Estimate", "AIC", "R²", paste("Pr(>", stat, ")", sep=""), " ")
    }

    if(!is.null(order)){
      if(order=="estim"){
        df <- df[order(abs(as.numeric(df[,2])), abs(as.numeric(df[,3]))),]
      } else if(order=="gof"){
        df <- df[order(abs(as.numeric(df[,3])), abs(as.numeric(df[,3]))),]
      } else if(order=="r_sq"){
        df <- df[order(abs(as.numeric(df[,4]))),]
      } else if(order=="p_val"){
        df <- df[order(abs(as.numeric(df[,5]))),]
      }
    }

    assign(paste("df_", g[[1]], sep=""), df, envir=.GlobalEnv)

    if(print){
      if(verbose){
        assign(paste("failed_", g[[1]], sep=""), failed, envir=.GlobalEnv)
        cat("\n please see ", paste("failed_", g[[1]], sep=""), " for a list of failed variables")
      }
      cat(paste("\n", as.character(eval(as.formula(as.character(g[[1]]))))[2], as.character(eval(as.formula(as.character(g[[1]]))))[1], as.character(eval(as.formula(as.character(g[[1]]))))[3], "+ ...\n\n", sep=" "))
      return(df)
    }

  }

}
