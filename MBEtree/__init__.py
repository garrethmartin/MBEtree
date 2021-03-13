#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt


def epanechnikov_density(
    X,
    sigmaopt,
    lambdaopt,
    dt,
    verbose=False,
    BallTree=True,
    chunks=1,
    weights=None,
    ):

    rhos = []
    ndim = X.shape[1]
    N = X.shape[0]

    # Calculate volume on n-sphere

    Vd = np.pi ** (ndim / 2.) / gamma(ndim / 2. + 1.)

    # Split into chunks if necessary

    for (X_i, lambda_i) in zip(np.array_split(X, chunks),
                               np.array_split(lambdaopt, chunks)):

        # Obtain Euclidean distances to all points within the kernel bandwidth for each point

        if BallTree:
            (inds, dists) = dt.query_radius(X_i, sigmaopt * lambda_i,
                    return_distance=True)
        else:
            (dist, ind) = dt.query(X_i, N,
                                   distance_upper_bound=sigmaopt
                                   * np.max(lambdaopt))
            dists = [d[np.isfinite(d)] for d in dist]
            inds = [i[np.isfinite(d)] for (i, d) in zip(ind, dist)]
        if weights is not None:
            ws = [weights[i] for i in inds]
        else:
            ws = np.zeros_like(dists)

        # For each point calculate the squared length for each nearby point

        t_dot_t = [(d / (sigmaopt * l)) ** 2 for (d, l) in zip(dists,
                   lambda_i)]

        # Calculate the density at each point using an Epanechnikok kernel

        for (t, l, w) in zip(t_dot_t, lambda_i, ws):
            Ke = (ndim + 2.) / (2 * Vd) * (1. - t[t < 1])
            if weights is not None:
                rho = 1. / N * np.sum([x * w for (x, w) in
                        zip((sigmaopt * l) ** -ndim * Ke, w[t < 1])])
            else:
                rho = 1. / N * np.sum((sigmaopt * l) ** -ndim * Ke)
            rhos.append(rho)
    return np.asarray(rhos)


def modified_breiman_density(
    X,
    verbose=False,
    BallTree=True,
    chunks=1,
    weights=None,
    alpha=None,
    ):
    '''
    purpose:
            Computes number density at each point using the modified Breiman density estimator with variable
            Epanechnikov kernel
    inputs:
            X: N-d array of positions: N-d array
            BallTree: if true use a ball tree to estimate distance between points, otherwise use KDTree. Ball tree is
            most efficient for ndim >=3. default=True
            chunks: number of chunks to break up the calculation of densities into. This can be helpful if the number
            of elements is large, especially when using the KDTree method: int, default=1
            weights: optional weighting for each point: 1-d array
            alpha: the sensitivity parameter: if none alpha=1/d
    outputs:
            rho: density at each point: 1-d array
    '''

    X = np.asarray(X)
    assert len(X.shape) == 2 and X.shape[0] > X.shape[1]
    if weights is not None:
        weights = np.asarray(weights)
        assert len(weights) == X.shape[0]
    if BallTree:
        from sklearn.neighbors import BallTree as tree
    else:
        from scipy.spatial import cKDTree as tree

    # Get the data dimensions

    d = X.shape[1]
    N = X.shape[0]
    if alpha == None:
        alpha = 1. / d

    # (0) Construct tree

    if verbose:
        print('Constructing tree...')
    dt = tree(X)

    # (1) Obtain initial pilot bandwidth using 20th and 80th percentiles

    P = np.percentile(X, [20, 80], axis=0)
    sigma = (P[1] - P[0]) / np.log(N)

    # Take minimum value of sigma to avoid over-smoothing

    sigmaopt = np.min(sigma)
    lambdaopt = np.ones(X.shape[0])

    # (2) Calculate the initial pilot density

    if verbose:
        print('Calculating pilot densities with fixed bandwidth...')
    pilot_rho = epanechnikov_density(
        X,
        sigmaopt,
        lambdaopt,
        dt,
        BallTree=BallTree,
        chunks=chunks,
        weights=weights,
        )

    # (3) Calculate the new local bandwith for each point

    if verbose:
        print('Calculating densities using variable bandwidth from pilot densities...'
              )
    g = np.exp(np.sum(np.log(pilot_rho) / N))
    lambdaopt_new = (pilot_rho / g) ** -alpha

    # (4) Use new bandwidths to calculate density

    rho = epanechnikov_density(
        X,
        sigmaopt,
        lambdaopt_new,
        dt,
        BallTree=BallTree,
        chunks=chunks,
        weights=weights,
        )

    return rho
