# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os


# A fazer:
#boxplot
# 1) alterar translate_curve para calcular param a em função do ultimo
# dos top n (fazer o a = score do ultimo ponto).
# 2) alterar classify_points para retornar C, D, Z de forma otimizada



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


def f_powerlaw(x, a, k):
    return a*(x)**(-k)

#------------------------------------------------------------------------------

def f_exp(x, a, b):
    return a*np.exp(-b*x)

#------------------------------------------------------------------------------



#==============================================================================
# FUNCTIONS
#==============================================================================



#------------------------------------------------------------------------------

def GetPercentilePoints(X, Y, window, slide, percentile, x_ini, x_end):
    '''
    Get percentile points to be used in the fitting curve,
    i.e. for each X, take only the max(abs(Y)) values.
    Obs: only points where X is within [x_ini, x_end] are returned.

    Parameters
    ----------
    X, Y:       numpy vectors.
    percentile: float, value.
                Percentile of points to be returned.
    window:     float, value.
                size of window.
    slide:      float, value.
                size of window's slides.

    Returns
    -------
    Xabs, Yabs: numpy vectors containing selected points
                by its percentile.
    '''
    Xabs, Yabs = [], []
    #
    half_window = window / 2.0
    n_slides = (x_end - x_ini) / slide
    xspace = np.linspace(x_ini + half_window, x_end - half_window,
                         n_slides, endpoint = False)
    for xp in xspace:
        indX = np.logical_and(X >= xp-half_window, X < xp+half_window)
        # if window is not empty
        if indX.any():
            Ycurr = np.abs(Y[indX])
            yp = np.percentile(Ycurr, percentile)
            #
            Xabs += [xp]
            Yabs += [yp]
    #
    Xabs = np.asarray(Xabs)
    Yabs = np.asarray(Yabs)
    #
    return Xabs, Yabs

#------------------------------------------------------------------------------

def PercentileCurveFitting(func, X, Y, window, slide, percentile, x_ini, x_end):
    '''
    Returns the Curve Fitting parameters to a given percentile of points.
    '''
    #
    x_ini = float(min(X) if x_ini is None else x_ini)
    x_end = float(max(X) if x_end is None else x_end)
    #
    Xabs, Yabs = GetPercentilePoints(X, Y, window, slide, percentile, x_ini, x_end)
    # curve fitting
    popt, pcov = curve_fit(func, Xabs, Yabs)
    #
    return popt, pcov

#------------------------------------------------------------------------------

#def _compute_norm_delta(dfScores, func, popt):
#    X = dfScores['X']
#    Y = dfScores['Y']
#    #
#    dfScores['S'] = abs(Delta_norm)
#    # Rank by X
#    dfScores['RankX'] = dfScores['X'].rank(ascending=False)
#    # Rank by S: SCORE(x) = abs(y(x) / yf(x))
#    dfScores['RankS'] = dfScores['S'].rank(ascending=False)
 


#------------------------------------------------------------------------------
	
def _read_extracted_features(file_out, new_header):
    fields = ['G/E', 'FreqC', 'FreqD']
    dfScores = pd.read_table(file_out, skiprows=5, usecols=fields)
    dfScores.columns = new_header
    return dfScores

#------------------------------------------------------------------------------

def _insert_new_features(dfScores):
    #
    C = dfScores['C']
    D = dfScores['D']
    #
    X = (C + D)
    Y = (D - C) / (C + D)
    #
    popt = PercentileCurveFitting(func, X, Y, window, slide, percentile,
                                  x_ini, x_end)[0]
    Yf = func(X, *popt)
    #
    # S: SCORE(x) = abs(y(x) / yf(x))
    Score_Delta_norm = abs(Y / Yf)
    #
    dfScores['X']  = X
    dfScores['Y']  = Y
    dfScores['Yf'] = Yf
    dfScores['S']  = Score_Delta_norm
    #
    dfScores['RankX'] = dfScores['X'].rank(ascending=False)
    dfScores['RankS'] = dfScores['S'].rank(ascending=False)

#------------------------------------------------------------------------------

def read_scores(file_out, new_header=['GENE','C','D']):
    dfScores = _read_extracted_features(file_out, new_header)
    _insert_new_features(dfScores)
    return dfScores

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# PARAMETERS
x_ini      = 5 ########
x_end      = None
percentile = 75
window     = 2
slide      = 1
func = f_powerlaw ########


