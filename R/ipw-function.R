#' Causal inference by Inverse Probability Weighting (IPW)
#'
#' This function conduct causal inference by implementing inverse probability
#' weighting (IPW). Please read
#' (https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/) before
#' applying this test.
#'
#' @param formula An object of class "formula": a symbolic description of the
#' model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param bootstrap An integer determining how many times this procedure
#' being repeated by resampling with replacement.
#' @param state An integer to set random seed for reproducible results.
#' @param verbose A logical determining whether a progress bar is shown.
#'
#' @return output A list containing the formula, exposure of interest, marginal
#' effect, 95% confidence interval (CI), significance by p-value obtained from
#' the CI (https://doi.org/10.1136/bmj.d2304), data, bootstrapping times,
#' random seed, and index for each bootstrap set.
#'
#' @keywords inverse probability weighting, IPW
#'
#' @export
#'
#' @examples
#'
#' # Load example data for formula and data
#' input=input_example()
#' formula=input$formula
#' data=input$data
#'
#' # Conduct g-formula
#' ipw(formula,data)

ipw=function(formula,data,bootstrap=30,state=33,verbose=F){

  once_ipw=function(formula,data){

    exposure=str_split(as.character(formula)[3],' \\+ ')[[1]]
    exposure=paste(
      c(exposure[1],paste(exposure[-1],collapse=' + '))
      ,collapse=' ~ '
    )
    exposure=as.formula(exposure)
    exposure_model=glm(exposure,binomial,data)
    exposure_prob=predict(exposure_model,type='response')

    exposure=str_split(as.character(formula)[3],' \\+ ')[[1]][1]
    data$ipew=1/(
      (exposure_prob*data[[exposure]])
      +((1-exposure_prob)*(1-data[[exposure]]))
    )

    data$censored_outcome=is.na(data[[as.character(formula)[2]]])
    data$censored_outcome=as.numeric(data$censored_outcome)
    censored_outcome=as.character(formula)
    censored_outcome=paste(
      c('censored_outcome',censored_outcome[3])
      ,collapse=' ~ '
    )
    censored_outcome=as.formula(censored_outcome)

    if(sum(data$censored_outcome)>0){
      censoring_model=
        suppressWarnings(
          glm(censored_outcome,binomial(link='logit'),data)
        )
      censoring_prob=predict(censoring_model,type='response')
    }else{
      censoring_prob=setNames(
        data[[as.character(formula)[2]]]*0
        ,rownames(data)
      )
    }

    data$ipcw=(
      (1-censoring_prob)*
      as.numeric(!is.na(data[[as.character(formula)[2]]]))
    )

    outcome=
      paste(
        c(as.character(formula)[2]
          ,str_split(as.character(formula)[3],' \\+ ')[[1]][1])
        ,collapse=' ~ '
      ) %>%
      as.formula()


    outcome_model=
      geeglm(
        formula=outcome
        ,family=gaussian(link='identity')
        ,data=data %>% mutate(id=seq(nrow(.)))
        ,weights=ipew*ipcw
        ,id=id
        ,corstr='exchangeable'
      )

    m_effect=outcome_model$coefficients[2]

    m_effect
  }

  set.seed(state)
  if(verbose){
    verboseapply=pblapply
  }else{
    verboseapply=lapply
  }
  bs_m_effect=
    seq(bootstrap) %>%
    verboseapply(function(X){
      index=sample(seq(nrow(data)),nrow(data),T)
      list(
        m_effect=once_ipw(formula,data[index,])
        ,index=index
      )
    })

  index=lapply(bs_m_effect,function(x)x$index)
  bs_m_effect=sapply(bs_m_effect,function(x)x$m_effect)

  m_effect=mean(bs_m_effect)+(-1:1)*qnorm(0.975)*sd(bs_m_effect)/sqrt(bootstrap)
  m_effect=setNames(m_effect[c(2,1,3)],c('mean','LB','UB'))
  z=mean(bs_m_effect)/(sd(bs_m_effect)/sqrt(bootstrap))
  p_value=exp(-0.717*z-0.416*z^2)

  output=
    list(
      formula=formula
      ,exposure=str_split(as.character(formula)[3],' \\+ ')[[1]][1]
      ,marginal_effect=m_effect[1]
      ,CI95_interval=m_effect[2:3]
      ,p_value=p_value
      ,data=data
      ,bootstrap=bootstrap
      ,state=state
      ,index=index
    )

  class(output)='ipw'

  assign(x='print.ipw',envir=baseenv(),value=function(x,digits=NULL){

    formula=
      as.character(x$formula)[2:length(as.character(x$formula))] %>%
      paste(collapse=' = ')

    if(is.null(digits)) digits=4
    marginal_effect=round(x$marginal_effect,digits)
    CI95_interval=round(x$CI95_interval,digits)
    p_value=ifelse(x$p_value<0.0001,'<0.0001',round(x$p_value,4))

    cat(paste0(
      '
    Method: Inverse Probability Weighting (IPW)

    Model: ',formula,'

    Exposure of interest: ',x$exposure,'
      Marginal effect: ',marginal_effect,'
      95% CI: ',paste(CI95_interval,collapse=' to '),'
      p-value: ',p_value,'

    Sample size: ',nrow(x$data),'
    Boostrapping times: ',x$bootstrap,'

    Callable variables in this object:
      ',paste(names(x)[1:3],collapse=', '),',
      ',paste(names(x)[4:7],collapse=', '),'
    '
    ))

  })

  output

}
