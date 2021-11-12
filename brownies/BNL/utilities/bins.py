from numpy import logspace, flip
import math


def equal_lethargy_bins(numBins, domainMin=1e-5, domainMax=20.0e6, reverse=False):
    seq = logspace(start=math.log10(domainMin), stop=math.log10(domainMax), num=numBins)
    if reverse:
        seq = flip(seq)
    return seq


