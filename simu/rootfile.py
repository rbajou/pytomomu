#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import uproot 
from pathlib import Path
from typing import Union
 
from g4hit import HitCollection
from event import Event, EventCollection


class RootFile: 

    def __init__(self, file:Union[Path, str]) -> None:
        self.file = file
        self.content = {}

    def open(self,  treename:str="Muons", **kwargs): 
        if isinstance(self.file, Path) : self.file = str(self.file)
        self.rootfile = uproot.open(self.file + f':{treename}')
        _branch_names = list(self.rootfile.keys())
        self.content =  self.rootfile.arrays(_branch_names, library='np', **kwargs)


class SimuOutput(RootFile):

    def __init__(self, file:Union[Path, str], treename:str="Muons", **kwargs) -> None: #treename = "Muons", "MetaData"
        RootFile.__init__(self, file=file)
        self.open(treename = treename, **kwargs)
        print('rootfile branches : ', list(self.content.keys()))

    def fill_event_collection(self, ndet:int=4, nh:int=8, mrotx:np.ndarray=np.diag((1,1,1)), mrotz:np.ndarray=np.diag((1,1,1))):
        '''
        Fill an event collection object keeping for the 'hit_collection' the same 'cluster_collection' dimensions for xyz and energy deposits arrays, as in real analyse data
        i.e xyz.shape = (ndet, nh).

        ndet (int) : number of micromegas detectors
        nh (int) : number of the most energetic hits to keep per detector 
        '''
        
        self.event_collection = EventCollection()
        lk = list(self.content.keys())
        nev = self.content[lk[0]].shape[0]
        k=0
        for i in range(nev):
            hit_collection = object.__new__(HitCollection)
            hit_collection.nhit = self.content[f'HitNb'][i]
            edep = self.content[f'HitEdep'][i]
            x, y, z = self.content[f'HitPosX'][i], self.content[f'HitPosY'][i], self.content[f'HitPosZ'][i]
            xyz = np.vstack((x, y, z)).T
            # xyz = np.dot(xyz, np.dot(mrotz, mrotx))
            xyz = xyz @ mrotz @ mrotx
            mat_xyz, mat_edep = np.ones((ndet, nh, 3))*-999, np.zeros((ndet, nh))
            det_ix = self.content['HitDet'][i]
            
            if len(xyz) != 0:
                for j in range(ndet) : 
                    ix = np.argwhere(det_ix==j).flatten()
                    order = np.flip(np.argsort(edep[ix])) #hits sorted in energy deposit decreasing order 
                    if len(ix) == 0: continue
                    k = nh
                    if len(ix) < nh : k = len(ix)
                    mat_xyz[j,:k] = xyz[ix][order][:nh]
                    mat_edep[j,:k] = edep[ix][order][:nh]
                    
            hit_collection.edep = mat_edep
            hit_collection.position = mat_xyz
            event = Event()
            event.id = int(self.content['Event'][i])
            event.xyz = mat_xyz
            event.hit_collection = hit_collection
            self.event_collection.__setitem__(event)
            k+=1
        

if __name__ == "__main__":
    pass