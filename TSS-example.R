library(MASS)
library(glmnet)
library(changepoint)
library(SIS)
library(ncvreg)
library(stabs)
library(car)
library(plsRglm)
library(rJava)
library(glmulti)
library(lars)
library(dplyr)


# Define Data -------------------------------------------------------------

# read in a csv file (place in the base of the repo)
model_data <- read.csv("model_data.csv")
variables_to_remove <- c("clinical_stageCT3")
endpoint_event_col <- c("mr")

model_data = model_data |> 
  dplyr::select(-all_of(variables_to_remove)) |>
  dplyr::select(all_of(endpoint_event_col), everything()) |>
  na.omit() 

# Separate predictors and outcome
y <- model_data[[endpoint_event_col]]
x <- model_data[ , !(names(model_data) %in% endpoint_event_col)]

# Define Arguments ---------------------------------------------------------------

nrepeat = 50 #number of subsamples in Stage 1 
nrepeat2 = 50 #number of subsamples in Stage 2
percent_split = 0.5 #take subsamples of half size
nsim = 1 #number of times you run the method
family="binomial" #type of outcome

set.seed(10) 

# Run Function ------------------------------------------------------------


results <- tss(
  x = x,
  y = y,
  nrepeat = 50,
  nrepeat2 = 50,
  nsim = 1,
  percent_split = .5,
  cutoff = .5, 
  family = "binomial")


print(results)

freq_TSS=paste(results$models_TSS)
as.data.frame(table(freq_TSS))

