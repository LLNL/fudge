import numpy
import math
import scipy

# ---------------------------------------------------------------------------------
#
# Utilities
#
# ---------------------------------------------------------------------------------

def bounded_nonlinear_curve_fit(func, xdata, ydata,
                                yerr=None, fitparms0=None, oparms=None, oparmscov=None,
                                bounds=(-numpy.inf, numpy.inf), size=100, method=None, verbose=False):
    """
    Assume test function called like this:

        ydata = f(xdata, *fitparms, **kwrds) + eps

    But, here's the interesting thing, the other parameters (`oparms`) can be uncertain and assumed to have a
    multivariate normal PDF with covariance `oparmscov`.  If you want to use this aspect of the fitter, make sure
    `oparms` is a keyword argument of your function `func`

    Rather than think too hard, we're going to Monte-Carlo it

    :param func : callable.
                  FIXME The model function, f(x, ...). It must take the independent variable as the first argument
                        and the parameters to fit as separate remaining arguments.

    :param xdata : An M-length sequence or an (k,M)-shaped array for functions with k predictors.  The independent
                   variable where the data is measured.

    :param ydata : M-length sequence  The dependent data - nominally f(xdata, ...)

    :param yerr : None or M-length sequence or MxM array, optional.  Determines the uncertainty in ydata. If we
                  define residuals as r = ydata - f(xdata, *popt), then the interpretation of sigma depends on its
                  number of dimensions:

                    * A 1-d sigma should contain values of standard deviations of errors in ydata. In this case,
                      the optimized function is chisq = sum((r / sigma) ** 2).

                    * A 2-d sigma should contain the covariance matrix of errors in ydata. In this case, the
                      optimized function is chisq = r.T * inv(sigma) * r.

                    * None (default) is equivalent of 1-d sigma filled with ones.

    :param fitparms0 : None, scalar, or N-length sequence, optional  Initial guess for the parameters. If None,
                       then the initial values will all be 1 (if the number of parameters for the function can be
                       determined using introspection, otherwise a ValueError is raised).
    :params oparms: FIXME

    :parms oparmscov: FIXME

    :parms size: FIXME

    :param bounds : 2-tuple of array_like, optional.  Lower and upper bounds on parameters. Defaults to no bounds.
                    Each element of the tuple must be either an array with the length equal to the number of parameters,
                    or a scalar (in which case the bound is taken to be the same for all parameters.) Use np.inf with
                    an appropriate sign to disable bounds on all or some parameters.  Passed to scipy.optimize.curve_fit

    :param method : {'lm', 'trf', 'dogbox'}, optional.  Method to use for optimization. See least_squares for more
                    details.  Default is 'lm' for unconstrained problems and 'trf' if bounds are provided. The method
                    'lm' won't work when the number of observations is less than the number of variables, use 'trf' or
                    'dogbox' in this case.  Passed to scipy.optimize.curve_fit
    """

    # If no oparm covariance, we don't need to Monte-Carlo anything so it's just a
    # plain old non-linear optimization
    if oparms is None or oparmscov is None:
        raise ValueError(
            "Either oparms or oparmscov is None, don't use this routine, use scipy.optimize.curve_fit directly")

        # Since have oparm covariance, apply Monte-Carlo forward propagation
    Rs = []
    # Rcovs=[]
    for i, oparms_realization in enumerate(numpy.random.multivariate_normal(oparms, oparmscov, size=size)):
        def newfunc(xdata, *fitparms):
            return func(xdata, *fitparms, oparms=oparms_realization)

        thisR, thisRcov = scipy.optimize.curve_fit(newfunc, ydata, xdata, p0=fitparms0, sigma=yerr, bounds=bounds,
                                                   method=method)
        Rs.append(thisR)
        #    Rcovs.append(thisRcov)
        if verbose:
            print(i, thisR, thisRcov)
    # numpy.mean(Rcovs, axis=0)  # FIXME: don't think I need to save these?
    # FIXME: get funky numbers when hit fitting boundaries and/or the
    #        various matrix inverses fail 'cuz of kernels.  Need tests?
    return numpy.mean(Rs, axis=0), numpy.cov(Rs, rowvar=False)



def linear_regression(xs, ys):
    """
    Performs linear regression fit to the specified points.

    We assume `ys[i] = a + b * xs[i] + eps[i]` where eps is stochastic noise at point `xs[i]`

    :param xs: list-like (best is numpy.array) list of x values
    :param ys: lisk-like (best is numpy.array) list of y values
    :return: (a, b, stddev(b))
    """
    xave = numpy.average(xs)
    yave = numpy.average(ys)
    covxy = numpy.cov(xs, ys)
    varx = covxy[0, 0]
    vary = covxy[1, 1]
    b = covxy[0, 1] / varx
    a = yave - b * xave
    db = math.sqrt(vary/(len(xs)-2))/math.sqrt(varx)
    return a, b, db
