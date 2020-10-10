import matplotlib.pyplot as pl
from utils import *
from pyproj import Proj, transform

class CartesianImage:
    def __init__(self, ppi, dim=500, order=0):
        assert (dim % 2 == 0)
        assert (ppi.anglestep is not None)
        assert (ppi.rangestep is not None)

        # print("dim = {}".format(dim))
        self.dtype = ppi.dtype
        self.dim = dim

        self.scanRange = (ppi.nbins * ppi.rangestep) * 1000
        diff = (2 * self.scanRange) / dim
        self.__x = self.__y = np.arange(-self.scanRange, self.scanRange, diff)
        self.__xx, self.__yy = np.meshgrid(self.__x, self.__y)
        self.data = polar_to_cart(ppi.data, ppi.anglestep, ppi.rangestep * 1000, self.__x, self.__y, order=order)

        if self.dtype == 'dBZ':
            self.data[self.data <= 0.0] = float('NaN')

    def reproject(self, source_crs, dest_crs):
        epsg4326 = Proj('epsg:4326')
        if isinstance(source_crs, Proj):
            pass
        elif isinstance(source_crs, str):
            source_crs = Proj(source_crs)
        else:
            raise TypeError("source_crs is instance of {}, expected str or pyproj.Proj".format(type(source_crs)))

        if isinstance(dest_crs, Proj):
            pass
        elif isinstance(dest_crs, str):
            dest_crs = Proj(dest_crs)
        else:
            raise TypeError("dest_crs is instance of {}, expected str or pyproj.Proj".format(type(dest_crs)))

        self.__xx, self.__yy = transform(source_crs, dest_crs, self.__xx, self.__yy)

    def getBoundingBox(self):
        epsg4326 = Proj('epsg:4326')
        epsg3857 = Proj('epsg:3857')
        lrx, lry = transform(epsg3857, epsg4326, np.max(self.__xx), np.min(self.__yy))
        ulx, uly = transform(epsg3857, epsg4326, np.min(self.__xx), np.max(self.__yy))
        return [[ulx, uly], [lrx, lry]]

    def plot(self, outpath, scale, bounds):
        #   set the figure aspect ratio to match aspect ratio of data in destination CRS
        xdim = abs(np.min(self.__xx) - np.max(self.__xx))
        ydim = abs(np.min(self.__yy) - np.max(self.__yy))
        aspect = ydim / xdim
        fsize = self.dim / 100
        fx, fy = fsize, fsize * aspect
        fig = pl.figure(figsize=(fx, fy))

        #   set up axes
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
        ax.set_xlim(np.min(self.__xx), np.max(self.__xx))
        ax.set_ylim(np.min(self.__yy), np.max(self.__yy))
        ax.set_aspect('equal', adjustable='box')
        ax.axis('off')

        #   plot the data
        pl.pcolormesh(self.__xx, self.__yy, self.data, figure=fig, cmap=scale, vmin=bounds[0], vmax=bounds[1])
        pl.savefig(outpath, dpi=100, transparent=True)
        pl.close('all')
        del fig
        del ax