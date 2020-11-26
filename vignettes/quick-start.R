# Load this package

library(gmethods)


# Load example data for formula and data

input=input_example()
formula=input$formula
data=input$data


# Conduct g-methods

results_gformula=gformula(formula,data,verbose=T)
results_ipw=ipw(formula,data,verbose=T)
results_gestimation=gestimation(formula,data,verbose=T)


# Results

## g-formula

results_gformula


## Inverse probability weighting (IPW)

results_ipw


# g-estimation

results_gestimation
