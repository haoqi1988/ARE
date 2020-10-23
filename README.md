README
================

## GitHub Documents

This package implements the algorithm developed in Nickel, Jiang and
Tresp, 2014.

## How to install

To install from GitHub, run this code

    devtools::install_github("haoqi1988/ARE")

## To prepare inputs:

This function requires 14 inputs: 4 are dimension parameters: N, number
of entities (based on data). K, number of layers in X / number of link
categories (each layer is one kind of link) (user specified). r, number
of latent attributes for each entity (arbitrary number). p, number of
the kinds of indirect relationships between entities to consider (user
specified). for example, (just for demo purposes)

    N = 100
    K = 2
    r = 3
    p = 3

4 are the names of these dims, for example

    entity_names <- paste0('Entity',1:N)
    link_names <- paste0('layer',1:K)
    attribute_names <- paste0('latentAttr', 1:r)
    indirect_link_names <- paste0('IndiLink', 1:p)

3 are hyperparameters, these values are quite arbitrary and not very
sensitive.  
lamdaW, an arbitrary weight to adjust the weight of W, bigger value
penalizes complex W. lamdaA, an arbitrary weight to adjust the weight of
A, bigger value penalizes complex A. lamdaR, an arbitrary weight to
adjust the weight of R, bigger value penalizes complex R. so for
example,

    lamdaW <- 0.1 
    lamdaA <- 1   
    lamdaR <- 0.1 

1 parameter is to control how accurate the learning needs to be.
error\_tolerance, user specified, the smaller the more accurate estimate
and longer time to converge. For large datasets, choose a bigger number
for shorter runtime. For smaller datasets, smaller values are
acceptable. for example,

    error_tolerance <- 0.33 

After specifying the parameters, the user needs to prepare the input
data. First, build X, the multilayer network, based on your empirical
data. X is the 3-d array of direct relationships between entities to be
factorized, dim(X) is N by N by k, for example, (here i am using random
data for demo
    purposes)

    X_1 <- matrix(sample(0:1, N, replace=TRUE), nrow = N, ncol = N, dimnames = list(entity_names,entity_names))
    X_2 <- matrix(sample(0:1, N, replace=TRUE), nrow = N, ncol = N, dimnames = list(entity_names,entity_names))
    X <- abind(X_1, X_2, along = 3)

Second, Build M, from empirical data, based/dependent on X, M is the 3-d
array of indirect relationship between entities, dim(M) is N by N by p,
for example, (here i am using co-1, co-2 and 1-2 for demo purposes)

    M_1 <- t(X_1) %*% X_1
    M_2 <- t(X_2) %*% X_2
    M_3 <- X_2 %*% X_1
    M <- abind(M_1, M_2, along = 3)
    M <- abind(M, M_3, along = 3)

## How to use the package

After all the inputs are specified, you can use the ARE() function to
factorize X, for example: (make sure the X,M,A,R,W and the dimnames are
consistent with the dimensions specified)

    results <- ARE(N = 100,
                   K = 2,
                   r = 3,
                   p = 3,
                   entity_names = entity_names,
                   link_names = link_names,
                   attribute_names = attribute_names,
                   indirect_link_names = indirect_link_names,
                   lamdaA = lamdaA, # weight of A
                   lamdaR = lamdaR, # weight of R
                   lamdaW = lamdaW, # weight of W
                   X = X,      # here X is a 3-d array, does not need to be a tensor
                   M = M,
                   error_tolerance = 0.33)

This way, the results list will contain A, R and W

## Plot Precision Recall Curve

The package also comes with a function that plots the within-sample
precision recall curve. For example

    x <- PlotPRCurve(A = results$A,
                R = results$R,
                W = results$W,
                X = X,
                M = M)
    x
    plot(x)

For more details, please reference the original paper, or read this
blogpost of mine: \[link here\]
