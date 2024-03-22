#!/usr/bin/python3
# -*- coding: utf-8 -*-

from dataclasses import dataclass
from typing import List, Union
import numpy as np

from cluster import StripCollection, ClusterCollection
from detector import DetectorEnsemble, Bench, Telescope, VetoPlane
from tracking import Track
from simu.g4hit import HitCollection


class Event:

    @property
    def id(self):
        return self._id
    
    @id.setter
    def id(self, value:int):
        self._id = value
    
    @property 
    def time(self):
        return self._time
    
    @time.setter
    def time(self, value:float):
        self._time = value

    @property 
    def strip_collection(self):
        return self._strip_collection
    
    @strip_collection.setter
    def strip_collection(self, value:StripCollection):
        self._strip_collection = value

    @property 
    def cluster_collection(self):
        return self._cluster_collection
    
    @cluster_collection.setter
    def cluster_collection(self, value:ClusterCollection):
        self._cluster_collection = value

    @property 
    def hit_collection(self):
        return self._hit_collection
    
    @hit_collection.setter
    def hit_collection(self, value:HitCollection):
        self._hit_collection = value

    @property 
    def track(self):
        return self._track
    
    @track.setter
    def track(self, value:Track):
        self._track = value

    @property
    def xyz(self): 
        return self._xyz 
    
    @xyz.setter
    def xyz(self, value:np.ndarray):
        self._xyz = value


@dataclass
class EventCollection:

    def __init__(self) :
        self.events = {} #events
        self._tracks = list()
        # self._df = pd.DataFrame(columns=['id', 'time'])

    def __setitem__(self, event:Event):
        self.events[event.id] = event
        # _dict = {'id':event.id, 'time':event.time}
        # _dftmp = pd.DataFrame(data=[_dict])
        # self._df = pd.concat([self._df.astype(_dftmp.dtypes), _dftmp.astype(self._df.dtypes)])
        
    def __getitem__(self, id:int):
        return self.events[id]

    # def to_df(self):
    #     self._df = pd.DataFrame.from_records([self.events]).iloc[0]

    # @property
    # def df(self): 
    #     return self._df
    
    @property
    def detector_ensemble(self):
        return self._ensemble
    
    @property
    def detectors(self):
        return self._detectors

    @detector_ensemble.setter
    def detector_ensemble(self, value:Union[DetectorEnsemble, Bench]):
        self._ensemble = value
        self._detectors = self._ensemble.detectors
    
    @property
    def xyz(self):
        return self._xyz
    
    @xyz.setter
    def xyz(self, value:np.ndarray):
        self._xyz = value
    
    @property
    def tracks(self):
        return self._tracks

    @tracks.setter
    def tracks(self, value:Union[np.ndarray, List]):
        self._tracks = value

    def get_xyz(self, mask:np.ndarray=None):
        '''
        Conversion from cluster position (in strips) to detector coordinates
        '''
        det_xyz = np.array([d.position for d in self._detectors])
        det_pitch = np.array([d.strip_layers[0].pitch for d in self._detectors])
        det_side = np.array([d.side for d in self._detectors])
        xyz_off = np.zeros((len(self._detectors), 3))
        xyz_off = np.array([-det_side[i]/2 + det_xyz[i]  for i, d in enumerate(self._detectors) ]) 

        lay_index = [lay.index for det in self._detectors for lay in det.strip_layers]
        nev = len(self.events)
        id0 = list(self.events.keys())[0]
        shp = self.events[id0].cluster_collection.amplitude[lay_index,:].shape
        nlay, nclus = shp[0], shp[1]
        npan = int(nlay/2)
        xmm, ymm, zmm = np.zeros((nev, npan, nclus)), np.zeros((nev, npan, nclus)), np.zeros((nev, npan, nclus))
        self.mask_xyz = np.ones((nev, npan, nclus), dtype=bool)
        det_axis = 1
        dim_array = np.ones((1,zmm.ndim),int).ravel()
        dim_array[det_axis] = -1
    
        ix_even, ix_odd = range(0, nlay, 2), range(1, nlay, 2)
       
        # m = (mask[:, ix_even] == True) & (mask[:, ix_odd] == True)
      
        for i, (id, ev) in enumerate(self.events.items()):
            cc = ev.cluster_collection
            xc, yc = cc.position[ix_even], cc.position[ix_odd] 
            xy = np.stack((xc, yc), axis=-1)
            xmm[i] = (xc.T * det_pitch + xyz_off[:,0]).T 
            ymm[i] = (yc.T * det_pitch + xyz_off[:,1]).T
            zmm[i] = (np.ones(xc.shape).T * det_xyz[:,-1]).T
            m = np.all(  xy != 0 , axis=-1 )
            if mask is not None: 
                m = m & (mask[i, ix_even] == True) & (mask[i, ix_odd] == True) 
            xmm[i, ~m], ymm[i, ~m], zmm[i, ~m] = -999., -999., -999.
            ev.xyz = np.stack((xmm[i], ymm[i], zmm[i]),  axis=-1)

        self._xyz = np.stack((xmm, ymm, zmm), axis=-1)
   

    def get_tracks(self, XYZ:np.ndarray=None, mask:np.ndarray=None) -> None:

        if mask is None: mask = np.full((XYZ.shape[0], XYZ.shape[1]), True)

        for i, (id, ev) in enumerate(self.events.items()): 
            xyz, m = XYZ[i], mask[i]           
            ev.track = Track()
            ev.track.id = id
            ev.track.fit_minres(xyz[m]) 
            if np.all(ev.track.a != -999.)==True :
                self._tracks.append(ev.track)
    
     
        
        



if __name__ == "__main__":

    pass
    