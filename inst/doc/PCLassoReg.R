## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("PCLassoReg")

## ---- eval=FALSE--------------------------------------------------------------
#  devtools::install_github("weiliu123/PCLassoReg")

## -----------------------------------------------------------------------------
library("PCLassoReg")

## ----load data----------------------------------------------------------------
# load data
data(survivalData)
data(PCGroups)

x <- survivalData$Exp
y <- survivalData$survData

## ----view survData------------------------------------------------------------
head(survivalData$survData)

## ----get protein complexes----------------------------------------------------
# get human protein complexes
PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
                        Type = "EntrezID")

## ----Split data set-----------------------------------------------------------
set.seed(20150122)
idx.train <- sample(nrow(x), round(nrow(x)*2/3))
x.train <- x[idx.train,]
y.train <- y[idx.train,]
x.test <- x[-idx.train,]
y.test <- y[-idx.train,]

## ----fit model----------------------------------------------------------------
# fit cv.PCLasso model
cv.fit1 <- cv.PCLasso(x = x.train, y = y.train, group = PC.Human, nfolds = 5)

## ----plot norm----------------------------------------------------------------
# plot the norm of each group
plot(cv.fit1, norm = TRUE)

## ----plot coef----------------------------------------------------------------
# plot the individual coefficients
plot(cv.fit1, norm = FALSE)

## ----plot cve-----------------------------------------------------------------
# plot the cross-validation error (deviance)
plot(cv.fit1, type = "cve")

## ----lambda.min---------------------------------------------------------------
cv.fit1$cv.fit$lambda.min

## -----------------------------------------------------------------------------
# Selected protein complexes at lambda.min
sel.groups <- predict(object = cv.fit1, type="groups",
                      lambda = cv.fit1$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
sel.groups <- predict(object = cv.fit1, type="groups",
                      lambda = c(0.1, 0.05))

## -----------------------------------------------------------------------------
# The number of risk protein complexes at lambda.min
sel.ngroups <- predict(object = cv.fit1, type="ngroups",
                       lambda = cv.fit1$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
sel.ngroups <- predict(object = cv.fit1, type="ngroups",
                      lambda = c(0.1, 0.05))

## -----------------------------------------------------------------------------
# The coefficients of protein complexes at lambda.min
groups.norm <- predict(object = cv.fit1, type="coefficients",
                       lambda = cv.fit1$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
groups.norm <- predict(object = cv.fit1, type="coefficients",
                       lambda = c(0.1, 0.05))

## -----------------------------------------------------------------------------
# Selected genes/proteins at lambda.min
sel.vars <- predict(object = cv.fit1, type="vars",
                    lambda=cv.fit1$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
sel.vars <- predict(object = cv.fit1, type="vars",
                    lambda=c(0.1, 0.05))

## -----------------------------------------------------------------------------
# The number of risk genes/proteins at lambda.min
sel.nvars <- predict(object = cv.fit1, type="nvars",
                     lambda=cv.fit1$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
sel.vars <- predict(object = cv.fit1, type="nvars",
                    lambda=c(0.1, 0.05))

## -----------------------------------------------------------------------------
# Selected genes/proteins at lambda.min
sel.vars <- predict(object = cv.fit1, type="vars.unique",
                    lambda=cv.fit1$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
sel.vars <- predict(object = cv.fit1, type="vars.unique",
                    lambda=c(0.1, 0.05))
# The number of risk genes/proteins at lambda.min
sel.nvars <- predict(object = cv.fit1, type="nvars.unique",
                     lambda=cv.fit1$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
sel.vars <- predict(object = cv.fit1, type="nvars.unique",
                    lambda=c(0.1, 0.05))

## -----------------------------------------------------------------------------
# predict risk scores of samples in x.test
s <- predict(object = cv.fit1, x = x.test, type="link",
             lambda=cv.fit1$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
s <- predict(object = cv.fit1, x = x.test, type="link",
             lambda=c(0.1, 0.05))

## ----load class data----------------------------------------------------------
# load data
data(classData)
data(PCGroups)

x <- classData$Exp
y <- classData$Label

## ----get protein complexes 2--------------------------------------------------
# get human protein complexes
PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
                        Type = "GeneSymbol")

## ----Split data set 2---------------------------------------------------------
set.seed(20150122)
idx.train <- sample(nrow(x), round(nrow(x)*2/3))
x.train <- x[idx.train,]
y.train <- y[idx.train]
x.test <- x[-idx.train,]
y.test <- y[-idx.train]

## ----fit model 2--------------------------------------------------------------
cv.fit2 <- cv.PCLasso2(x = x.train, y = y.train, group = PC.Human,
                       penalty = "grLasso", family = "binomial", nfolds = 10)

## ----plot norm 2--------------------------------------------------------------
# plot the norm of each group
plot(cv.fit2, norm = TRUE)

## ----plot coef 2--------------------------------------------------------------
# plot the individual coefficients
plot(cv.fit2, norm = FALSE)

## ----plot cve 2---------------------------------------------------------------
# plot the cross-validation error (deviance)
plot(cv.fit2, type = "cve")

## ----lambda.min 2-------------------------------------------------------------
cv.fit2$cv.fit$lambda.min

## ----predict groups 2---------------------------------------------------------
# Selected protein complexes at lambda.min
sel.groups <- predict(object = cv.fit2, type="groups",
                      lambda = cv.fit2$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
sel.groups <- predict(object = cv.fit2, type="groups",
                      lambda = c(0.1, 0.05))

## ----predict ngroups 2--------------------------------------------------------
# The number of risk protein complexes at lambda.min
sel.ngroups <- predict(object = cv.fit2, type="ngroups",
                       lambda = cv.fit2$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
sel.ngroups <- predict(object = cv.fit2, type="ngroups",
                      lambda = c(0.1, 0.05))

## -----------------------------------------------------------------------------
# The coefficients of protein complexes at lambda.min
groups.norm <- predict(object = cv.fit2, type="coefficients",
                       lambda = cv.fit2$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
groups.norm <- predict(object = cv.fit2, type="coefficients",
                       lambda = c(0.1, 0.05))

## -----------------------------------------------------------------------------
# Selected genes/proteins at lambda.min
sel.vars <- predict(object = cv.fit2, type="vars",
                    lambda=cv.fit2$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
sel.vars <- predict(object = cv.fit2, type="vars",
                    lambda=c(0.1, 0.05))

## -----------------------------------------------------------------------------
# The number of risk genes/proteins at lambda.min
sel.nvars <- predict(object = cv.fit2, type="nvars",
                     lambda=cv.fit2$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
sel.vars <- predict(object = cv.fit2, type="nvars",
                    lambda=c(0.1, 0.05))

## -----------------------------------------------------------------------------
# Selected genes/proteins at lambda.min
sel.vars <- predict(object = cv.fit2, type="vars.unique",
                    lambda=cv.fit2$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
sel.vars <- predict(object = cv.fit2, type="vars.unique",
                    lambda=c(0.1, 0.05))
# The number of risk genes/proteins at lambda.min
sel.nvars <- predict(object = cv.fit2, type="nvars.unique",
                     lambda=cv.fit2$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
sel.vars <- predict(object = cv.fit2, type="nvars.unique",
                    lambda=c(0.1, 0.05))

## -----------------------------------------------------------------------------
# predict probabilities of samples in x.test
s <- predict(object = cv.fit2, x = x.test, type="response",
             lambda=cv.fit2$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
s <- predict(object = cv.fit2, x = x.test, type="response",
             lambda=c(0.1, 0.05))

## ----predict class------------------------------------------------------------
# predict class labels of samples in x.test
s <- predict(object = cv.fit2, x = x.test, type="class",
             lambda=cv.fit2$cv.fit$lambda.min)
# For values of lambda not in the sequence of fitted models, linear
# interpolation is used.
s <- predict(object = cv.fit2, x = x.test, type="class",
             lambda=c(0.1, 0.05))

## -----------------------------------------------------------------------------
# load data
data(survivalData)
data(PCGroups)

x = survivalData$Exp
y = survivalData$survData

PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
                        Type = "EntrezID")

# fit PCSCAD model
fit.PCSCAD <- PCLasso(x, y, group = PC.Human, penalty = "grSCAD", gamma = 6)

# fit PCMCP model
fit.PCMCP <- PCLasso(x, y, group = PC.Human, penalty = "grMCP", gamma = 5)

## -----------------------------------------------------------------------------
# load data
data(classData)
data(PCGroups)

x = classData$Exp
y = classData$Label

PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
                        Type = "GeneSymbol")

# fit PCSCAD model
fit.PCSCAD2 <- PCLasso2(x, y, group = PC.Human, penalty = "grSCAD",
                       family = "binomial", gamma = 10)

# fit PCMCP model
fit.PCMCP2 <- PCLasso2(x, y, group = PC.Human, penalty = "grMCP",
                      family = "binomial", gamma = 9)

