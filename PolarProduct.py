from utils import polar_to_cart
import datetime as dt
import numpy as np


class PolarProduct:
    def __init__(self, data, metadata):
        self.data = data
        self.nrays = data.shape[0]
        self.nbins = data.shape[1]

        self.dtype = metadata['dtype']
        self.anglestep = metadata['anglestep']
        self.rangestep = metadata['rangestep']
        self.radar = metadata['radname']
        self.radID = metadata['radid']
        self.__date = dt.datetime.strptime(metadata['datetime'], "%H:%M:%S %Y-%m-%d")

    def resample(self, dimension):
        newData = np.full((dimension, self.data.shape[1]), 0.0)
        for i in range(self.nbins):
            newData[:, i] = np.interp(np.arange(0, self.nrays, self.nrays / dimension), np.arange(self.nrays),
                                      self.data[:, i])
        return newData

    def getDate(self, format):
        return self.__date.strftime(format)

    def toCart(self, dim, order=0):
        # dim - dimension of output cartesian grid
        assert (dim % 2 == 0)
        assert (self.anglestep is not None)
        assert (self.rangestep is not None)

        scanRange = self.nbins * self.rangestep
        diff = (2 * scanRange) / dim
        x = y = np.arange(-scanRange, scanRange, diff)
        cart_data = polar_to_cart(self.data, self.anglestep, self.rangestep, x, y, order=order)
        return cart_data

    def nanify(self):
        self.data[self.data == self.datamin] = np.nan


class PolarPPI(PolarProduct):
    def __init__(self, data, volumeMetadata, productMetadata):
        super().__init__(data, volumeMetadata)

        self.elevation = productMetadata["elevation"]
        self.tilt = productMetadata["tilt"]


class PolarCAPPI(PolarProduct):
    def __init__(self, data, volumeMetadata, productMetadata):
        super().__init__(data, volumeMetadata)

        self.height = productMetadata["height"]
        self.rmin = productMetadata["startrange"]
        self.rmax = productMetadata["stoprange"]


class PolarCMAX(PolarProduct):
    def __init__(self, data, volumeMetadata, productMetadata):
        super().__init__(data, volumeMetadata)

        self.hmin = productMetadata["hmin"]
        self.hmax = productMetadata["hmax"]


class PolarVIL(PolarProduct):
    def __init__(self, data, volumeMetadata, productMetadata):
        super().__init__(data, volumeMetadata)

        self.hmin = productMetadata["hmin"]
        self.hmax = productMetadata["hmax"]
