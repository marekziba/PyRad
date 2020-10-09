from abc import ABC, abstractmethod
import numpy as np
import wradlib as wrl
from PolarProduct import PolarProduct
import utils as u

class PolarVolume(ABC):
    def getCappi(self, cappi_h, pseudo=True):
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
        for ele in range(0, self.numele):
            ele = ele + 1
            pos = self.numele - ele
            while rh[pos, nbin] < 0 and nbin < 249:
                cappi[:, nbin] = (volume[pos + 1, :, nbin] * abs_rh[pos, nbin] + volume[pos, :, nbin] *
                                      abs_rh[pos + 1, nbin]) / (abs_rh[pos, nbin] + abs_rh[pos + 1, nbin])
                if nbin < 249:
                    nbin = nbin + 1
        if nbin < 249:
            # print("nbin<249")
            for nbin in range(nbin, 250):
                if pseudo:
                    cappi[:, nbin] = volume[0, :, nbin]
                else:
                    cappi[:, nbin] = np.nan
        cappi[cappi == 0.0] = np.nan
        return cappi

    @abstractmethod
    def getPPI(self, n, raw=False): pass

    @abstractmethod
    def date(self): pass



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

        metadata = dict()
        metadata['anglestep'] = 360 / azirange
        metadata['rangestep'] = rangestep
        metadata['min'] = datamin
        metadata['max'] = datamax
        metadata['numslice'] = elevation
        metadata['radid'] = self.rbdict['volume']['sensorinfo']['@id']
        metadata['radname'] = self.rbdict['volume']['sensorinfo']['@name']
        metadata['datetime'] = self.rbdict['volume']['scan']['@time'] + " " + self.rbdict['volume']['scan']['@date']
        metadata['elevation'] = self.rbdict['volume']['scan']['slice'][elevation]['posangle']
        metadata['dtype'] = self.rbdict['volume']['scan']['slice'][elevation]['slicedata']['rawdata']['@type']

        ppr = PolarProduct(data, metadata)
        ppr.anglestep = azires
        ppr.datamin = datamin
        ppr.datamax = datamax

        return ppr

    def getPPI(self, elevation, raw=False):
        if not raw:
            return self.ppis[elevation]
        else:
            return self.voldata[elevation]

    def date(self):
        return self.rbdict['volume']['scan']['@time'] + " " + self.rbdict['volume']['scan']['@date']

    def __init__(self, path, aziCount = 360):
        self.rbdict = wrl.io.read_rainbow(path)
        self.ppis = []
        self.numele = int(self.rbdict['volume']['scan']['pargroup']['numele'])
        self.stoprange = float(self.rbdict['volume']['scan']['pargroup']['stoprange'])
        self.rangestep = float(self.rbdict['volume']['scan']['pargroup']['rangestep'])
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

