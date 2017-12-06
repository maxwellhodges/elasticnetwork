from rpy2 import robjects
from rpy2.robjects.numpy2ri import numpy2ri
from rpy2.robjects.packages import importr
import numpy as np

def quant_reg(x, y, q):
    """ Use R library quantreg to do a linear quantile fit to datas.  
    More precisely.
    
    Parameters
    ----------
    x : list/:numpy:array
      independent variable
    y : list/:numpy:array
      dependent variable
    q : float 
      quantile to fit (in [0,1])
    
    Returns
    -------
    ??
      y-values for linear fit at input x-values

    """
    r = robjects.r
    try:
        r('library(quantreg)')
    except:
        r("install.packages('quantreg')")
        r("library(quantreg)")
    r('library(splines)')
    robjects.numpy2ri.activate()
    x = robjects.FloatVector(x)
    y = robjects.FloatVector(y)
    robjects.globalenv['x'] = x
    robjects.globalenv['y'] = y
    robjects.globalenv['q'] = q
    r('qrfit <- rq(y ~ x, tau = q)')
    yfit = r('yfit <- fitted(qrfit)')
    return yfit