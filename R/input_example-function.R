#' Make an input example for gmethods package
#'
#' This function load a causal model as formula object and a data frame of
#' nhefs table, a cleaned version of NHEFS data. In Causal Inference: What If
#' book by Hernán and Robins, four NHEFS datasets are described in this book
#' (https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/). The
#' remaining datasets from that book could be retrieve from cidata R package
#' (https://github.com/malcolmbarrett/cidata) for more information.
#'
#' @return output A list of a formula and a data frame with dimension of 1629
#' rows and 10 columns.
#'
#' @keywords example data
#'
#' @export
#'
#' @examples
#'
#' # Load example data for formula and data
#' input=input_example()
#' formula=input$formula
#' data=input$data

input_example=function(){

  output=
    list(
      formula=formula(wt82_71~qsmk+sex+poly(age,2)+race+poly(smokeyrs,2))
      ,data=
        read_csv('../data/nhefs.csv',col_types=cols()) %>%
        select_at(c('wt82_71','qsmk','sex','age','race','smokeyrs'))
    )

  output

}

