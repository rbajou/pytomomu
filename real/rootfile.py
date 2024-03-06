#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import uproot 
from pathlib import Path
from typing import Union
 


class RootFile: 

    def __init__(self, file:Union[Path, str]) -> None:
        self.file = file
        self._branch_names = []
        self.content = {}

    def open(self,  treename:str="T", ): 
        if isinstance(self.file, Path) : self.file = str(self.file)
        self.rootfile = uproot.open(self.file + f':{treename}')
        self._branch_names = list(self.rootfile.keys())

    def get_content(self, n_events:int, branch_list:list=None,start_event:int=0, 
                     step_size:int=50000): #, save_npy=False, save_path=None, save_name=None):
        '''
        Import the branches with names given in branch_list as array numpy
        '''
        
        if branch_list is None: 
            branch_list = self._branch_names.copy()

        rf = self.rootfile

        n_events_tot = int(rf.num_entries)
        n_events = min(n_events, n_events_tot)
 
        del_branch =[]
        for branch in branch_list:
            bshape = list(np.array(rf[branch].array(entry_start = 0, entry_stop = 5)).shape)
            if bshape[0] == 0 : 
                del_branch.append(branch)
                continue

            bshape[0] = n_events
            # print(f"bshape = {bshape}")
            self.content[branch] = np.zeros(shape=bshape)

        for b in del_branch : 
            branch_list.remove(b)
            self._branch_names.remove(b) 
     
        index = 0
        for batch in rf.iterate(step_size=step_size, filter_name=branch_list, entry_start=start_event, entry_stop=min(start_event+n_events, n_events_tot), library='np'):
            batch_size = batch[branch_list[0]].shape[0]
            for branch in branch_list:
                self.content[branch][index:index+batch_size] = np.array(batch[branch])
            index += batch_size
            print(f"[Data Import] {index} events imported")        



def general_import(file:Union[str, Path], branch_list:list, n_events:int, treename:str="T", start_event:int=0, step_size:int=50000, save_npy=False, save_path=None, save_name=None):
    '''
    Import the branches with names given in branch_list as array numpy
    '''
    if isinstance(file, Path) : file = str(file)
    
    rf = uproot.open(file+f':{treename}')
    n_events_tot = int(rf.num_entries)
    
    n_events = min(n_events, n_events_tot)

    dict = {}
    for branch in branch_list:
        bshape = list(np.array(rf[branch].array(entry_start = 0, entry_stop = 1)).shape)
        bshape[0] = n_events
        dict[branch] = np.zeros(shape=bshape)

    index = 0
    for batch in rf.iterate(step_size=step_size, filter_name=branch_list, entry_start=start_event, entry_stop=min(start_event+n_events, n_events_tot), library='np'):
        batch_size = batch[branch_list[0]].shape[0]
        # batch = batch.reshape(s[0]*s[-1], s[1])
        for branch in branch_list:
            dict[branch][index:index+batch_size] = np.array(batch[branch])
        index += batch_size
        print(f"[Data Import] {index} events imported")

    if save_npy==True:
        for branch in branch_list:
            np.save(save_path + save_name + '_' + branch)
    return dict

if __name__ == "__main__":
    pass