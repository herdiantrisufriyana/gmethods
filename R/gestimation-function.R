#' Causal inference by g-estimation
#'
#' This function conduct causal inference by implementing g-estimation. Please
#' read (https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/)
#' before applying this test.
#'
#' @param formula An object of class "formula": a symbolic description of the
#' model to be fitted. The exposure of interest should be plugged-in as the
#' first covariate at the right-hand size of the formula.
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
#' @keywords g-estimation
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
#' gestimation(formula,data)

gestimation=function(formula,data,bootstrap=30,state=33,verbose=F){

  once_gestimation=function(formula,data){

    data_notNA=subset(data,!is.na(data[[as.character(formula)[2]]]))

    Hpsi=function(outcome,exposure){
      function(psi) outcome-psi*exposure
    }

    exposure_formula=str_split(as.character(formula)[3],' \\+ ')[[1]]
    exposure_formula=
      paste(
        c(exposure_formula[1]
          ,paste(exposure_formula[-1],collapse=' + '))
        ,collapse=' ~ '
      )
    exposure_formula=as.formula(exposure_formula)

    absCoefHpsi=function(exposure_formula,data){
      function(psi){
        exposure_model=
          suppressWarnings(
            glm(
              formula=exposure_formula
              ,family=binomial(link='logit')
              ,data=data
            )
          )

        data$..Hpsi..=
          Hpsi(
            outcome=data[[as.character(formula)[2]]]
            ,exposure=data[[as.character(exposure_formula)[2]]]
          )(psi) %>%
          setNames('Hpsi')

        suppressWarnings(update(exposure_model,formula=~.+..Hpsi..)) %>%
          coef() %>%
          tail(1) %>%
          abs()
      }
    }

    absCoefHpsiFun=absCoefHpsi(exposure_formula,data_notNA)

    psiArgMinHpsi=
      suppressWarnings(
        optim(
          par=0
          ,fn=absCoefHpsiFun
          ,lower=-10
          ,upper=10
        )
      )

    m_effect=psiArgMinHpsi$value

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
        m_effect=once_gestimation(formula,data[index,])
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

  class(output)='gestimation'

  assign(x='print.gestimation',envir=baseenv(),value=function(x,digits=NULL){

    formula=
      as.character(x$formula)[2:length(as.character(x$formula))] %>%
      paste(collapse=' = ')

    if(is.null(digits)) digits=16
    marginal_effect=round(x$marginal_effect,digits)
    CI95_interval=round(x$CI95_interval,digits)
    p_value=ifelse(x$p_value<0.0001,'<0.0001',round(x$p_value,4))

    cat(paste0(
      '
    Method: g-estimation

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
