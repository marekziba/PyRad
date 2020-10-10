from PolarVolume import PolarVolume
import numpy as np


class InterpolatedVolume:
    def __init__(self, volume: PolarVolume, hmin: float, hmax: float, vres: float):
        self.__hmin = hmin
        self.__hmax = hmax
        self.__vres = vres
        self.__volume = volume
        nlevs = int((hmax - hmin) / vres)
        nrays, nbins = volume.getDimensions()
        dims = (nlevs, nrays, nbins)
        self.data = np.full(dims, 0.0)
        self.__levels = []
        for i in range(nlevs):
            h = i * vres * 1000  # converting from km to meters
            self.data[i, :, :] = volume.getCappi(h, pseudo=False)
            self.__levels.append(h)
        self.data[np.isnan(self.data)] = -999.0

    def __validateProductBounds(self, bounds):
        if len(bounds) != 2:
            return False
        if bounds[0] >= bounds[1]:
            return False
        if bounds[0] < self.__hmin or bounds[1] > self.__hmax:
            return False
        return True

    def __getIndexes(self, bounds):
        hmin = self.__hmin * 1000
        hmax = self.__hmax * 1000
        if bounds is not None:
            if not self.__validateProductBounds(bounds):
                raise ValueError("Specified boundary exceeds volume dimensions")
            else:
                hmin = bounds[0] * 1000
                hmax = bounds[1] * 1000

        startIdx = next(x for x, val in enumerate(self.__levels) if val >= hmin)
        stopIdx = len(self.__levels) - next(x for x, val in enumerate(reversed(self.__levels)) if val <= hmax)
        return (startIdx, stopIdx)

    def getCmax(self, bounds=None, absoluteMax=True):
        startIdx, stopIdx = self.__getIndexes(bounds)
        data = self.data.copy()
        data[data == -999.0] = 0.0
        print("CMAX ({},{})".format(startIdx, stopIdx))
        if not absoluteMax:
            cmax = np.max(data[startIdx:stopIdx], axis=0)
        else:
            maxArray = np.max(data[startIdx:stopIdx], axis=0)
            minArray = np.min(data[startIdx:stopIdx], axis=0)
            cmax = np.where(maxArray < abs(minArray), minArray, maxArray)
        cmax[cmax == 0.0] = np.nan


        return cmax

    def getVIL(self, bounds=None):
        startIdx, stopIdx = self.__getIndexes(bounds)
        data = self.data.copy()
        data = data[startIdx:stopIdx]
        vil = ((10 ** (data / 10)) / 24000) ** (1 / 1.82)
        noData = ((10 ** (-999.0 / 10)) / 24000) ** (1 / 1.82)
        vil[vil == noData] = 0.0
        vil = np.sum(vil, axis=0)
        vil[vil == 0.0] = np.nan
        return vil
