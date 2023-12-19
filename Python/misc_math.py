import numpy as np
from math import sqrt
from scipy.optimize import curve_fit

def mean2(x) :
    return np.sum(x) / np.size(x)

def corr2(a,b) :
    a = a - mean2(a)
    b = b - mean2(b)
    return (a*b).sum() / sqrt((a*a).sum() * (b*b).sum())

def stdDev(x, y) : 
    mean = np.sum(x*y)/np.sum(y)
    sigma = np.sqrt(np.sum(y*(x-mean)**2)/np.sum(y))
    
    return mean, sigma

def gauss(x, a, x0, sigma) :
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def gaussFit(x,y) :
    
    mean, sigma = stdDev(x, y)
    popt, pcov = curve_fit(gauss, x, y, p0=[np.max(y), mean, sigma])
    return gauss(x, *popt)

def gaussRadius(x,y,radiusType = '1/e') :
    
    if radiusType == '1/e' :
        idx = np.argwhere(y > np.max(y)/np.exp(1))
        a = np.min(x[idx])
        b = np.max(x[idx])
        rad = (b-a)/2
    if radiusType == '1/e2' :
        idx = np.argwhere(y > np.max(y)/np.exp(2))
        a = np.min(x[idx])
        b = np.max(x[idx])
        rad = (b-a)/2        
    if radiusType == 'FWHM' :
        idx = np.argwhere(y > np.max(y)/2)
        a = np.min(x[idx])
        b = np.max(x[idx])
        rad = (b-a)/2
    return rad

def sech(x) :
    return 1/np.cosh(x)