from abc import ABC, abstractmethod
import numpy as np
import wradlib as wrl
from PolarProduct import *
import utils as u

class PolarVolume(ABC):
    def getCappi(self, cappi_h, pseudo=True, raw=False):
        cappi = np.full((self.nrays, self.nbins), np.nan)
        rh = np.full((self.numele, self.nbins), 0.0)

        for i in range(self.numele):
            ele = float(self.ppis[i].elevation)
            rangestep = float(self.ppis[i].rangestep)
            rh[i,:] = u.getBeamHeight(ele,self.stoprange,rangestep)

        rh = (rh * 1000) - cappi_h  # <--------------- convert beam height to meters and calculate deviation from specified CAPPI height
        rh = rh.astype(int)
        abs_rh = abs(rh)
        nbin = 0
        ele = 0
        pos = self.numele - ele - 1
        volume = self.voldata.copy()
        volume[np.isnan(volume)] = 0
        while rh[pos, nbin] < 0 and nbin < 249:
            if pseudo:
                cappi[:, nbin] = volume[pos, :, nbin]
            else:
                cappi[:, nbin] = np.nan
            nbin = nbin + 1
        startrange = nbin * self.rangestep
        for ele in range(0, self.numele):
            ele = ele + 1
            pos = self.numele - ele
            while rh[pos, nbin] < 0 and nbin < 249:
                cappi[:, nbin] = (volume[pos + 1, :, nbin] * abs_rh[pos, nbin] + volume[pos, :, nbin] *
                                      abs_rh[pos + 1, nbin]) / (abs_rh[pos, nbin] + abs_rh[pos + 1, nbin])
                if nbin < 249:
                    nbin = nbin + 1
        stoprange = nbin * self.rangestep
        if nbin < 249:
            # print("nbin<249")
            for nbin in range(nbin, 250):
                if pseudo:
                    cappi[:, nbin] = volume[0, :, nbin]
                else:
                    cappi[:, nbin] = np.nan
        cappi[cappi == 0.0] = np.nan

        if raw:
            return cappi
        else:
            volumeMetadata = self.getVolumeMetadata()

            productMetadata = dict()
            productMetadata['height'] = cappi_h
            productMetadata['startrange'] = startrange
            productMetadata['stoprange'] = stoprange

            return PolarCAPPI(cappi, volumeMetadata, productMetadata)

    def getDimensions(self):
        return (self.nrays, self.nbins)

    @abstractmethod
    def getPPI(self, n, raw=False): pass

    @abstractmethod
    def date(self): pass

    @abstractmethod
    def getVolumeMetadata(self): pass

    @abstractmethod
    def _getProjection(self): pass



class Rb5Volume(PolarVolume):
    def __extractElevationData(self, elevation, fillAzimuth = False, nan = True):
        data = self.rbdict['volume']['scan']['slice'][elevation]['slicedata']['rawdata']['data']
        datadepth = float(self.rbdict['volume']['scan']['slice'][elevation]['slicedata']['rawdata']['@depth'])
        datamin = float(self.rbdict['volume']['scan']['slice'][elevation]['slicedata']['rawdata']['@min'])
        datamax = float(self.rbdict['volume']['scan']['slice'][elevation]['slicedata']['rawdata']['@max'])

        try:
            azi = self.rbdict['volume']['scan']['slice'][elevation]['slicedata']['rayinfo']['data']
            azidepth = float(self.rbdict['volume']['scan']['slice'][elevation]['slicedata']['rayinfo']['@depth'])
            azirange = float(self.rbdict['volume']['scan']['slice'][elevation]['slicedata']['rayinfo']['@rays'])
        except(TypeError):
            azi = self.rbdict['volume']['scan']['slice'][elevation]['slicedata']['rayinfo'][0]['data']
            azidepth = float(self.rbdict['volume']['scan']['slice'][elevation]['slicedata']['rayinfo'][0]['@depth'])
            azirange = float(self.rbdict['volume']['scan']['slice'][elevation]['slicedata']['rayinfo'][0]['@rays'])
        try:
            azires = float(self.rbdict['volume']['scan']['slice'][elevation]['anglestep'])
        except(KeyError):
            try:
                azires = float(self.rbdict['volume']['scan']['pargroup']['anglestep'])
            except(KeyError):
                azires = 1

        try:
            rangestep = float(self.rbdict['volume']['scan']['slice'][elevation]['rangestep'])
        except(KeyError):
            try:
                rangestep = float(self.rbdict['volume']['scan']['pargroup']['rangestep'])
            except(KeyError):
                rangestep = 1

        azi = (azi * azirange / 2 ** azidepth) * azires
        azi = azi.round().astype(int)
        nroll = azi[0] + 1

        data = data.astype(np.double)

        if data.shape[1] < 250 and fillAzimuth:
            shape = (data.shape[0], 250 - data.shape[1])
            filler = np.full(shape, 0, dtype=np.double)
            data = np.concatenate((data, filler), axis=1)

        data = np.roll(data, nroll, 0)
        data = datamin + data * (datamax - datamin) / 2 ** datadepth
        if nan:
            data[data == datamin] = float('NaN')

        # metadata = dict()
        # metadata['anglestep'] = 360 / azirange
        # metadata['rangestep'] = rangestep
        # metadata['min'] = datamin
        # metadata['max'] = datamax
        # metadata['radid'] = self.rbdict['volume']['sensorinfo']['@id']
        # metadata['radname'] = self.rbdict['volume']['sensorinfo']['@name']
        # metadata['datetime'] = self.rbdict['volume']['scan']['@time'] + " " + self.rbdict['volume']['scan']['@date']
        # metadata['dtype'] = self.rbdict['volume']['scan']['slice'][elevation]['slicedata']['rawdata']['@type']

        volumeMetadata = self.getVolumeMetadata()

        productMetadata = dict()
        productMetadata['elevation'] = self.rbdict['volume']['scan']['slice'][elevation]['posangle']
        productMetadata['tilt'] = elevation

        ppr = PolarPPI(data, volumeMetadata, productMetadata)

        return ppr

    def getPPI(self, elevation, raw=False):
        if not raw:
            return self.ppis[elevation]
        else:
            return self.voldata[elevation]

    def date(self):
        return self.rbdict['volume']['scan']['@time'] + " " + self.rbdict['volume']['scan']['@date']

    def getVolumeMetadata(self):
        metadata = dict()
        metadata['dtype'] = self.dtype
        metadata['anglestep'] = self.anglestep
        metadata['rangestep'] = self.rangestep
        metadata['radname'] = self.rbdict['volume']['sensorinfo']['@name']
        metadata['radid'] = self.sensor
        metadata['datetime'] = self.date()
        return metadata

    def _getProjection(self):
        self.lon = float(self.rbdict['volume']['sensorinfo']['lon'])
        self.lat = float(self.rbdict['volume']['sensorinfo']['lat'])
        self.proj4_def = "+proj=aeqd +lat_0=" + str(self.lat) + " +lon_0=" + str(
            self.lon) + " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

    def __init__(self, path, aziCount = 360):
        self.rbdict = wrl.io.read_rainbow(path)
        self._getProjection()
        self.ppis = []
        self.numele = int(self.rbdict['volume']['scan']['pargroup']['numele'])
        self.stoprange = float(self.rbdict['volume']['scan']['pargroup']['stoprange'])
        self.rangestep = float(self.rbdict['volume']['scan']['pargroup']['rangestep'])
        self.anglestep = float(self.rbdict['volume']['scan']['pargroup']['anglestep'])
        try:
            self.dtype = self.rbdict['volume']['scan']['slice'][0]['slicedata']['rawdata']['@type']
        except:
            self.dtype = self.rbdict['volume']['scan']['slice']['slicedata']['rawdata']['@type']
        self.sensor = self.rbdict['volume']['sensorinfo']['@id']
        #self.voldata = np.full((self.numele, aziCount, nbins))
        bins = []
        for i in range(self.numele):
            data = self.__extractElevationData(i, nan=True)
            self.ppis.append(data)
            bins.append(data.nbins)
        self.nbins = max(bins)
        self.voldata = np.full((self.numele, aziCount, self.nbins), 0.0)
        for i in range(self.numele):
            data = self.ppis[i].resample(360)
            if data.shape[1] < self.nbins:
                data = u.fillRange(data, self.nbins)
            self.voldata[i] = data
        self.nrays = aziCount

