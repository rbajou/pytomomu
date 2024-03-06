#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

class Track:

    @property
    def id(self):
        """Id evt"""
        return self._id
    
    @id.setter
    def id(self, value:int):
        self._id= value

    @property
    def a(self):
        """Coefficient directeur [ax, ay]"""
        return self._a

    @property
    def tthetax(self): 
        """"Tangente of projection angle track in (0xz)"""
        return self._tthetax

    @property
    def tthetay(self): 
        """"Tangente of projection angle track in (0yz)"""
        return self._tthetay
    
    @a.setter
    def a(self, value:np.ndarray):
        self._a = value
        self._tthetax, self._tthetay = value.T

    @property
    def b(self): 
        """Ordonnée à l'origine [bx, by]"""
        return self._b

    @b.setter
    def b(self, value:np.ndarray):
        self._b = value

    @property
    def r(self): 
        """"Résidus trace-cluster"""
        return self._r

    @r.setter
    def r(self, value:np.ndarray):
        self._r = value

    def intersection(self, z:float):
        x, y = self._a * z + self._b
        xyz = np.stack((x, y, z))
        return xyz

    def __dir__(self):
        return['a', 'b', 'tan_theta', 'tan_phi', 'res']

    
    def fit_ls(self, x:np.ndarray, y:np.ndarray):

        '''
        author : Paul Lardillier
        
        Computes a and b such as y = a*x + b, according to the least squared method
        x.shape = y.shape = (nclus,)
        '''
        assert x.shape == y.shape, 'x and y must have the same shape'
        axis = 0
        vxs = np.var(x, axis=axis)
        covxy = np.mean(x*y, axis=axis) - np.mean(x, axis=axis)*np.mean(y, axis=axis)
        a = covxy/vxs
        b = np.mean(y, axis=axis) - a * np.mean(x, axis=axis)    
        yp = a * x.T + b
        r = y - yp.T

        return [a, b, r]


    def fit_minres(self, xyz:np.ndarray):
        '''
        author : Paul Lardillier

        Returns the fit when taking into account only n-1 detectors, keeping the fit with minimum residual
        x and y are (nclus,)
        ''' 
        x, y, _ = xyz.T 
     
        assert x.shape == y.shape, 'x and y must have the same shape'

        ndet = xyz.shape[0]
        nfit = ndet
     
        arr_a, arr_b, arr_r = np.zeros((2,ndet)), np.zeros((2,ndet)), np.zeros((2, ndet, ndet-1))

        af, bf, rf = -999.*np.ones(2), -999.*np.zeros(2), -999.*np.zeros((2, ndet))

        for c in range(2):
            for i in range(nfit): 
                mask = np.full(nfit, True)
                mask[i] = False
                xyz_mask = xyz[mask]
                x, y = xyz_mask[:,0,-1], xyz_mask[:,0,c] #shape (ndet-1, 1)
                nnull = ( x != -999. )
                if len(nnull[nnull==True]) >= ndet - 1: 
                    a0, b0, r0 = self.fit_ls(x[nnull], y[nnull])
                else : 
                    a0, b0, r0 = -999., -999., -999.*np.ones(nfit-1)
                arr_a[c, i] = a0
                arr_b[c, i] = b0
                arr_r[c, i] = r0

            imin = np.argmin(np.abs(arr_r[c, :]), axis=0)[0]
            #We recompute the residuals, but with all the 4 detectors
            X, Y = xyz[:,0,-1], xyz[:,0,c]
            af[c], bf[c] = arr_a[c, imin], arr_b[c, imin]
            YP = af[c]* X.T + bf[c]
            rf[c] = Y - YP.T
                
        self.a, self._b, self._r = af, bf, rf            
        
    

if __name__ == "__main__":
    pass