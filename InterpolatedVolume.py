from PolarVolume import PolarVolume
import numpy as np

class InterpolatedVolume:
    def __init__(self, volume: PolarVolume, hmin: float, hmax: float, vres: float):
        self.__hmin = hmin
        self.__hmax = hmax
        self.__vres = vres
        nlevs = int((hmax - hmin)/vres)
        nrays, nbins = volume.getDimensions()
        dims = (nlevs, nrays, nbins)
        self.data = np.full(dims, 0.0)
        self.__levels = []
        for i in range(nlevs):
            h = i * vres * 1000     #converting from km to meters
            self.data[i, :, :] = volume.getCappi(h,pseudo=False)
            self.__levels.append(h)
        print(self.__levels)

    def __validateProductBounds(self, bounds):
        if len(bounds) != 2:
            return False
        if bounds[0] >= bounds[1]:
            return False
        if bounds[0] < self.__hmin or bounds[1] > self.__hmax:
            return False
        return True

    def getCmax(self, bounds = None):
        hmin = self.__hmin
        hmax = self.__hmax
        if bounds is not None:
            if not self.__validateProductBounds(bounds):
                raise ValueError("Specified boundary exceeds volume dimensions")
            else:
                hmin = bounds[0]
                hmax = bounds[1]

        startIdx = next(x for x, val in enumerate(self.__levels) if val >= hmin)
        stopIdx = len(self.__levels) - next(x for x, val in enumerate(reversed(self.__levels)) if val <= hmax)
        cmax = np.max(self.data[startIdx:stopIdx], axis=0)
        return cmax
