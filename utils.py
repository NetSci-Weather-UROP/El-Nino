import netCDF4 as nc
import numpy as np

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
    T = np.array(ds['air'])
    T = T[::4,1:-1,:]

    # load latitude and longitude gridpoints
    if latlon:
        lat = np.array(ds['lat'])[1:-1]
        lon = np.array(ds['lon'])
        return T, lat, lon

    return T


def avoid_seasonality(T):
    """
    Remove the effect of seasonality on the measurements.
    """

    # for every day
    for day in range(365):

        # compute the daily mean over all years at all gridpoints
        dmean = np.mean(T[day::365,:,:], axis=0)

        # compute the corresponding standard deviations at all gridpoints
        dstandev = np.std(T[day::365,:,:], axis=0)

        # update temperature data
        T[day::365,:,:] = (T[day::365,:,:] - dmean) / dstandev

    return T


def find_start_date(day, month):
    """
    Return the ordinal number of a date in a year.
    """
    days = np.array([0,31,28,31,30,31,30,31,31,30,31,30])
    startdate = np.sum(days[:month]) + day - 1

    return startdate


def find_inside_indices(lat, lon, x0, x1, y0, y1):
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

    return lat_in, lon_in


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
    lat_in, lon_in = find_inside_indices(lat, lon, x0, x1, y0, y1)
    
    # compute total number of nodes and number of inside nodes
    nodecount = np.size(lat) * np.size(lon)
    in_nodecount = np.size(lat_in) * np.size(lon_in)

    # total number of days to be examined (usually 365)
    days = np.shape(T)[0]

    lat, lon = np.meshgrid(lat,lon)

    # initalise gridpoint time series sets
    T_in = np.empty([in_nodecount, days+2])
    T_out = np.empty([nodecount-in_nodecount, days+2])

    # add time series to correct set
    inrow, outrow = 0, 0
    for i in range(np.shape(lat)[0]):
        for j in range(np.shape(lat)[1]):
            if (j in lat_in and i in lon_in):
                T_in[inrow,:-2] = T[:,j,i]
                T_in[inrow,-2:] = [i,j]
                inrow += 1
            else:
                T_out[outrow,:-2] = T[:,j,i]
                T_out[outrow,-2:] = [i,j]
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
        T = np.append(T, read_air_data(year, latlon=False), axis=0)

        # remove possible gap days
        if not year % 4:
            T = T[:-1,:,:]

    # standardise data
    T = avoid_seasonality(T)

    return T, lat, lon


def pearson_coeffs(x, y):
    """
    Calculates pearson coefficients between two sets  of time series:
    (mean(xy) - mean(x)mean(y)) / (std(x)std(y))
    """

    # number of rows and columns
    n, m = np.shape(x)[0], np.shape(y)[0]

    # number of observations
    k = np.shape(x)[1]

    # check matching shapes
    if not k == np.shape(y)[1]:
        raise ValueError(
            f"Both x and y must have the same number of observations, not {k} "
            f"and {np.shape(y)[1]}."
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
    pos = np.max(x, axis=axis)
    neg = np.min(x, axis=axis)
    return pos * (pos > -neg) + neg * (pos < -neg)


def comp_c(T_in, T_out, tau_max=200):
    """
    Compute the maximum of the correlations, the respective offset, mean and
    standard deviation over all offsets.
    """
    # rows and columns
    n, m = np.shape(T_in)[0], np.shape(T_out)[0]

    days = 365

    # temporary array holding corrcoefs for all offsets
    temp = np.empty([n,m,2*tau_max+1])
    t_in = T_in[:,tau_max:tau_max+days]

    # iterate through all offsets
    for tau in range(-tau_max, tau_max+1):
        if not tau % 50:
            print("Tau:", tau)
        # set offset window for outside nodes
        t_out = T_out[:,tau_max+tau:tau_max+days+tau]

        # append corrcoefs for specific offset
        temp[:,:,tau_max+tau] = pearson_coeffs(t_in, t_out)

    C = np.empty([n,m,4])

    # theta_i,j
    C[:,:,1] = np.argmax(np.abs(temp), axis=2) - tau_max
    # C(theta_i,j)
    C[:,:,0] = maximod(temp, axis=2) * (C[:,:,1] < 151) * (0 <= C[:,:,1])
    # mean over all theta
    C[:,:,2] = np.mean(temp[:,:,tau_max:tau_max+150], axis=2)
    # standard deviation over  all theta
    C[:,:,3] = np.std(temp[:,:,tau_max:tau_max+150], axis=2)

    return C


def reformat_c(C, T_out):
    """
    return M, containing local C, N, W data (plottable with latlon=True)
    """
    M = np.empty([
        np.max(T_out[:,-1]+1).astype(int), np.max(T_out[:,-2]+1).astype(int), 3
    ])

    M_N = np.sum((0 <=C[:,:,1])*(C[:,:,1] < 151), axis=0)
    M_C = np.sum(C[:,:,0], axis=0)
    M_W = np.sum((C[:,:,0]-C[:,:,2])/C[:,:,3], axis=0)

    i = 0
    for line in T_out:
        coord1, coord2 = line[-1].astype(int), line[-2].astype(int)
        M[coord1, coord2,:] = [M_N[i], M_C[i], M_W[i]]
        i += 1

    return M


def year_series(
        T, lat, lon, year, start_year=1948, window_day=1, window_month=7
    ):
    """
    Compute `crosscorr` for a specific year window in the time series.
    """

    # compute start date (1 July by default) index
    start_date = (
        (year-start_year) * 365 + find_start_date(window_day, window_month)
    )

    # compute end date index
    end_date = start_date + 565

    # split regions for this period
    T_in, T_out = separate_regions(
        T[start_date-200:end_date,:,:], lat, lon, 190, 240, -5, 5
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