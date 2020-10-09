import numpy as np
import math as m

RADIUS = 8493

def resampleData(data, dimension):
    newData = np.full((dimension,data.shape[1]),0)
    nrays, nbins = data.shape
    for i in range(data.shape[1]):
        #print("{} {}".format(np.arange(nrays).shape[0],data[:,i].shape[0]))
        newData[:,i] = np.interp(np.arange(0,nrays,nrays/dimension),np.arange(nrays),data[:,i])
    return newData

def isValidVolfile(path):
    with open(path) as volfile:
        start = volfile.readline()
        if start.startswith("<volume"):
            result = True
        else:
            result = False
    return result

def polar_to_cart(polar_data, theta_step, range_step, x, y, order=3):

    from scipy.ndimage.interpolation import map_coordinates as mp

    # "x" and "y" are numpy arrays with the desired cartesian coordinates
    # we make a meshgrid with them
    X, Y = np.meshgrid(x, y)

    # Now that we have the X and Y coordinates of each point in the output plane
    # we can calculate their corresponding theta and range
    Tc = np.degrees(np.arctan2(Y, X)).ravel()
    Rc = (np.sqrt(X**2 + Y**2)).ravel()

    # Negative angles are corrected
    Tc[Tc < 0] = 360 + Tc[Tc < 0]

    # Using the known theta and range steps, the coordinates are mapped to
    # those of the data grid
    Tc = Tc / theta_step
    Rc = Rc / range_step

    # An array of polar coordinates is created stacking the previous arrays
    coords = np.vstack((Tc, Rc))

    # To avoid holes in the 360º - 0º boundary, the last column of the data
    # copied in the begining
    polar_data = np.vstack((polar_data, polar_data[-1,:]))

    # The data is mapped to the new coordinates
    # Values outside range are substituted with nans
    cart_data = mp(polar_data, coords, order=order, mode='constant', cval=np.nan)

    # The data is reshaped and returned
    return(cart_data.reshape(len(y), len(x)).T)

def getBeamHeight(angle, stoprange, rangestep):
    x = np.arange(0 + rangestep / 2, stoprange + rangestep / 2, rangestep)
    height = x.copy()
    for i in range(len(height)):
        height[i] = x[i] * m.tan(m.radians(angle)) + (-m.sqrt(RADIUS ** 2 - x[i] ** 2) + RADIUS)
    return height

def fillRange(data, bins, fillValue = np.nan):
    filler = np.full((data.shape[0], bins - data.shape[1]), fillValue)
    newData = np.concatenate((data, filler), axis=1)
    return newData