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


def _avoid_seasonality(T, start_year, end_year, yearly_seasonal=False):
    """
    Removes the effect of seasonality on the measurements.
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

    lat_in, lon_in = _find_inside_indices(lat, lon, x0, x1, y0, y1)
    days = np.shape(T)[0]
    lat, lon = np.meshgrid(lat,lon)
    T_in, T_out = np.empty([57, days+2]), np.empty([10167, days+2])

    inrow, outrow = 0, 0
    for i in range(np.shape(lat)[0]):
        for j in range(np.shape(lat)[1]):
            x, y = lon[i,j], lat[i,j]
            if (j in lat_in and i in lon_in):
                T_in[inrow,:-2] = T[:,j,i]
                T_in[inrow,-2:] = [i,j]
                inrow += 1
            else:
                T_out[outrow,:-2] = T[:,j,i]
                T_out[outrow,-2:] = [i,j]
                outrow += 1

    return T_in, T_out


def get_data(start_year, end_year, gap_years=False):
    """
    Obtain all normalised measurements for a given period.

    Parameters:
    start_year - first year of data;
    end_year - last year of data (included);
    gap_years=False - include gap days.

    Returns:
    T - 3-D numpy array containing the measurements;
    lat - 1-D numpy array of lateral coordinates;
    lon - 1-D numpy array of longitudinal coordinates.
    """

    if gap_years:
        T, lat, lon = read_air_data(start_year)
        for year in range(start_year+1, end_year+1):
            T = np.append(T, read_air_data(year, latlon=False), axis=0)
    else:
        T, lat, lon = read_air_data(start_year)
        if not start_year % 4:
            T = T[:-1,:,:]
        for year in range(start_year+1, end_year+1):
            T = np.append(T, read_air_data(year, latlon=False), axis=0)
            if not year % 4:
                T = T[:-1,:,:]


    T = _avoid_seasonality(T, start_year, end_year)

    return T, lat, lon


def pearson_coeffs(x, y):
    """
    Calculates pearson coefficients between two sets  of time series:
    (mean(xy) - mean(x)mean(y)) / (std(x)std(y))
    """
    n, m = np.shape(x)[0], np.shape(y)[0]
    k = np.shape(x)[1]
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


def comp_c(T_in, T_out, tau_max=200, gap=False):
    """
    Compute the maximum of the correlations, the respective offset, mean and
    standard deviation over all offsets.
    """
    n, m = np.shape(T_in)[0], np.shape(T_out)[0]
    days = 365 + gap
    temp = np.empty([n,m,tau_max+1])
    t_in = T_in[:,:days]
    for tau in range(tau_max+1):
        if not tau % 25:
            print("Tau:", tau)
        t_out = T_out[:,tau:days+tau]
        temp[:,:,tau] = pearson_coeffs(t_in, t_out)
    C = np.empty([n,m,4])
    tempmax = np.empty([n,m,2])
    C[:,:,1] = np.argmax(abs(temp), axis=2)
    tempmax[:,:,0] = np.max(temp, axis=2) * (C[:,:,1] < 151)
    tempmax[:,:,1] = np.max(-temp, axis=2) * (C[:,:,1] < 151)
    C[:,:,0] = (
        - tempmax[:,:,1] * (tempmax[:,:,1] > tempmax[:,:,0])
        + tempmax[:,:,0] * (tempmax[:,:,1] < tempmax[:,:,0])
    )
    C[:,:,2] = np.mean(temp, axis=2)
    C[:,:,3] = np.std(temp, axis=2)
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