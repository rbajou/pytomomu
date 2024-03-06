#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from typing import List, Union
import numpy as np
from pathlib import Path
import uproot
import time
import re
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd

#package module(s)
from cluster import ClusterCollection
from event import Event, EventCollection
from rootfile import RootFile


class Analyse(RootFile):

    def __init__(self, file:Union[Path, str]) -> None:
        RootFile.__init__(self, file=file)
        self.open()
        print('rootfile branches : ', self._branch_names)
        try : 
            regex = r'MGv(\d)'
            matches = np.array([ re.search(regex, bn) for bn in self._branch_names])
            self.prefix_mg = matches[matches != None][0].group(0)
            
        except : 
            regex = r'INB72'
            matches = np.array([ re.search(regex, bn) for bn in self._branch_names])
            self.prefix_mg = matches[matches != None][0].group(0)
      

        self.evvar_names = [self._branch_names[i] for i, m in enumerate(matches) if m is None ]
        self.clusvar_names = [self._branch_names[i] for i, m in enumerate(matches) if m is not None ]
        print(f'MicroMegas version : {self.prefix_mg}')
        print(f"self.clusvar_names = {self.clusvar_names}")


    def fill_event_collection(self, prefix:str=None):
        
        if len(self.content) == 0 : raise ValueError("Run 'get_content(...)' first.")

        self.dict_clus_var = {k: self.content[k] for k in self._branch_names}
        lk = list(self.dict_clus_var.keys())
    
        nev = self.dict_clus_var[lk[0]].shape[0]

        if prefix is None: prefix = self.prefix_mg

        self.event_collection = EventCollection()
       
        for i in range(nev):
            cluster_collection = object.__new__(ClusterCollection)
            cluster_collection.nclus = self.dict_clus_var[f'{prefix}_NClus'][i] if f'{prefix}_NClus' in self.dict_clus_var.keys() else 1
            cluster_collection.position = self.dict_clus_var[f'{prefix}_ClusPos'][i]
            cluster_collection.amplitude = self.dict_clus_var[f'{prefix}_ClusAmpl'][i]
            cluster_collection.max_amplitude = self.dict_clus_var[f'{prefix}_ClusMaxStripAmpl'][i]
            cluster_collection.max_amplitude_strip = self.dict_clus_var[f'{prefix}_StripMaxAmpl'][i] if f'{prefix}_StripMaxAmpl' in self.dict_clus_var.keys() else 1
            cluster_collection.size = self.dict_clus_var[f'{prefix}_ClusSize'][i]
            cluster_collection.tot = self.dict_clus_var[f'{prefix}_ClusTOT'][i]
            cluster_collection.max_sample = self.dict_clus_var[f'{prefix}_ClusMaxSample'][i]
            cluster_collection.time = self.dict_clus_var[f'{prefix}_ClusT'][i]
           
            event = Event()
            event.id = int(self.dict_clus_var['evn'][i])
            event.time = self.dict_clus_var['evttime'][i] if f'evttime' in self.dict_clus_var.keys() else -999
            event.cluster_collection = cluster_collection

            self.event_collection.__setitem__(event)


class Final(RootFile):

    def __init__(self, file:Union[Path, str]) -> None:
        RootFile.__init__(self, file)
        self.open()




if __name__ == "__main__": 
    pass