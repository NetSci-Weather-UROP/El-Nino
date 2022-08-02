"""
This is a modified version of utils.py that uses CuPy (GPU-accelerated w/ CUDA)
to speed up the computation. This is partly an experiment to see how well
(and if at all) it works with generating the large dataset to use for plotting
the evolution data (recreating fig 3 in the PNAS paper). evolution.py calls
year_series 12x68 times and so far this modified script has not created any 
issues - it speeds up the calculations about thricefold.

Note: 
You will need to have the CUDA toolkit installed. CuPy and Numpy are
basically interchangeable apart from converting from one format to another
as evident in the code. Functions now take CuPy inputs for arrays and
return CuPy arrays.

Inconsistencies:
This script is only up to date with the original utils.py committed on 21/07
(bcbc30870e3538ac9c38c74307f2ba6a9e335855 by GR). It could be that some of the
calcultions involved have been changed since then, in which case they are 
not reflected in this script.
"""

from calendar import c
import netCDF4 as nc
import numpy as np
import cupy as cp

def read_air_data(year, latlon=True):
    """
    Read NCEP data for specified year.

    Parameters:
    year - integer specifying the year.

    Returns:
    T - (365, 71, 144) numpy array of temperature data every day at 12am;
    lat - (71,) numpy array of the latitude coordinates;
    lon - (144,) numpy array of the longitude coordinates.
    """

    # load dataset
    ds = nc.Dataset(f"air.sig995/air.sig995.{year}.nc")

    # get temperature data
    T = cp.array(ds['air'])
    T = T[::4,1:-1,:]

    # load latitude and longitude gridpoints
    if latlon:
        lat = cp.array(ds['lat'])[1:-1]
        lon = cp.array(ds['lon'])
        return T, lat, lon

    return T


def _avoid_seasonality(T, start_year, end_year, yearly_seasonal=False):
    """
    Remove the effect of seasonality on the measurements.
    """

    # for every day
    for day in range(365):

        # compute the daily mean over all years at all gridpoints
        dmean = cp.mean(T[day::365,:,:], axis=0)

        # compute the corresponding standard deviations at all gridpoints
        dstandev = cp.std(T[day::365,:,:], axis=0)

        # update temperature data
        T[day::365,:,:] = (T[day::365,:,:] - dmean) / dstandev

    return T


def _find_inside_indices(lat, lon, x0, x1, y0, y1):
    """
    Find indices for points inside a "rectangular" region on the globe.
    """

    # without crossing the 180° longitude
    if x0 < x1:
        lon_in = [i for i, coord in enumerate(lon) if x0 < coord < x1]
    # with crossing the 180° longitude
    else:
        lon_in = [
            i for i, coord in enumerate(lon)
            if (coord < x1 or x0 < coord)
        ]

    # y0 always less than y1
    lat_in = [j for j, coord in enumerate(lat) if y0 < coord < y1]

    return cp.asarray(lat_in), cp.asarray(lon_in)


def separate_regions(T, lat, lon, x0, x1, y0, y1):
    """
    Separate nodes based on being inside or outside a rectangular region on the
    globe.

    Parameters:
    T - temperature measurements at the nodes;
    lat - 1-D array containing the lateral coordinates;
    lon - 1-D array containing the longitudinal coordinates;
    x0 - western limit of the region;
    x1 - eastern limit of the region;
    y0 - southern limit of the region;
    y1 - northern limit of the region (has to satisfy `y1` > `y0`).

    Returns:
    T_in/T_out - 2-D numpy array of in/out-temperature measurements, as well
    as longitudinal and lateral coordinates.
    """

    # find regional indices
    lat_in, lon_in = _find_inside_indices(lat, lon, x0, x1, y0, y1)
    
    # compute total number of nodes and number of inside nodes
    nodecount = cp.size(lat) * cp.size(lon)
    in_nodecount = cp.size(lat_in) * cp.size(lon_in)

    # total number of days to be examined (usually 365)
    days = cp.shape(T)[0]

    lat, lon = cp.meshgrid(lat,lon)

    # initalise gridpoint time series sets
    T_in = cp.empty([in_nodecount, days+2])
    T_out = cp.empty([nodecount-in_nodecount, days+2])

    # add time series to correct set
    inrow, outrow = 0, 0
    for i in range(cp.shape(lat)[0]):
        for j in range(cp.shape(lat)[1]):
            if (j in lat_in and i in lon_in):
                T_in[inrow,:-2] = T[:,j,i]
                T_in[inrow,-2:] = cp.asarray([i,j])
                inrow += 1
            else:
                T_out[outrow,:-2] = T[:,j,i]
                T_out[outrow,-2:] = cp.asarray([i,j])
                outrow += 1

    return T_in, T_out


def get_data(start_year, end_year):
    """
    Obtain all normalised measurements for a given period.

    Parameters:
    start_year - first year of data;
    end_year - last year of data (included).

    Returns:
    T - 3-D numpy array containing the measurements;
    lat - 1-D numpy array of lateral coordinates;
    lon - 1-D numpy array of longitudinal coordinates.
    """

    # compute T and lat, lon first (so lat and lon don't have to be loaded
    # for every year)
    T, lat, lon = read_air_data(start_year)

    # remove possible gap day
    if not start_year % 4:
        T = T[:-1,:,:]

    # iterate through years, end_year inclusive
    for year in range(start_year+1, end_year+1):

        # append yearly data
        T = cp.append(T, read_air_data(year, latlon=False), axis=0)

        # remove possible gap days
        if not year % 4:
            T = T[:-1,:,:]

    # standardise data
    T = _avoid_seasonality(T, start_year, end_year)

    return T, lat, lon


def pearson_coeffs(x, y):
    """
    Calculates pearson coefficients between two sets  of time series:
    (mean(xy) - mean(x)mean(y)) / (std(x)std(y))
    """

    # number of rows and columns
    n, m = cp.shape(x)[0], cp.shape(y)[0]

    # number of observations
    k = cp.shape(x)[1]

    # check matching shapes
    if not k == cp.shape(y)[1]:
        raise ValueError(
            f"Both x and y must have the same number of observations, not {k} "
            f"and {cp.shape(y)[1]}."
        )

    return (
        (1/k * x @ y.transpose()
        - x.mean(axis=1).reshape([n,1]) * y.mean(axis=1))
        / (x.std(axis=1).reshape([n,1]) * y.std(axis=1))
    )

def maximod(x, axis=None):
    """
    Find the maximal modulus inside an array.
    """
    pos = cp.max(x, axis=axis)
    neg = cp.min(x, axis=axis)
    return pos * (pos > -neg) + neg * (pos < -neg)

def comp_c(T_in, T_out, tau_max=200):
    """
    Compute the maximum of the correlations, the respective offset, mean and
    standard deviation over all offsets.
    """
    # rows and columns
    n, m = cp.shape(T_in)[0], cp.shape(T_out)[0]

    days = 365

    # temporary array holding corrcoefs for all offsets
    temp = cp.empty([n,m,tau_max+1])
    t_in = T_in[:,:days]

    # iterate through all offsets
    for tau in range(tau_max+1):
        if not tau % 50:
            print(chr(9608), end="")
        # set offset window for outside nodes
        t_out = T_out[:,tau:days+tau]

        # append corrcoefs for specific offset
        temp[:,:,tau] = pearson_coeffs(t_in, t_out)
    print("")
    C = cp.empty([n,m,4])

    # theta_i,j
    C[:,:,1] = cp.argmax(cp.abs(temp), axis=2)
    # C(theta_i,j)
    C[:,:,0] = maximod(temp, axis=2) * (C[:,:,1] < 151)
    # mean over all theta
    C[:,:,2] = cp.mean(temp, axis=2)
    # standard deviation over  all theta
    C[:,:,3] = cp.std(temp, axis=2)

    return C


def year_series(T, lat, lon, year, start_year=1948):
    """
    Compute `crosscorr` for a specific year window in the time series.
    """

    # compute start date (1 July) index
    start_date = (year-start_year) * 365 + 182

    # compute end date index
    end_date = start_date + 565

    # split regions for this period
    T_in, T_out = separate_regions(
        T[start_date:end_date,:,:], lat, lon, 190, 240, -5, 5
    )

    # compute various parameters
    C = comp_c(T_in, T_out)

    return C, T_in, T_out


##################### DEFUNCT

# def season_indices(start_year, end_year, gap_years = False):
#     """
#     Give the indices corresponding to new seasons in a 4-year cycle.

#     Parameters:
#     start_year - the first year considered;
#     end_year - the last year considered.

#     Returns:
#     indices - 1-D numpy array of stard dates of new seasons since
#     01/01/`start_year`, including this, as well as 31/12/`end_year`
#     """
#     gap = 4 - start_year % 4
#     if not gap % 4:
#         gap = 0
#     indices = [59 + gap*gap_years]
#     for i in range(end_year - start_year + 1):
#         for d in [92, 92, 91, 90]:
#             indices.append(indices[-1]+d)
#     indices.pop()
#     indices = np.array(indices)
#     if gap_years:
#         for year in range(end_year - start_year + 1):
#             if not (year % 4 - gap):
#                 k = year*4
#                 if not k:
#                     indices[:] += 1
#                 else:
#                     indices[k:] += 1
#     indices = np.concatenate(([0], indices, [indices[-1]+31]))
#     return indices

#####################