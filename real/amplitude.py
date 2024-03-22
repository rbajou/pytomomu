#!/usr/bin/python3
# -*- coding: utf-8 -*-

from typing import List, Union
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

from event import EventCollection


params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'axes.labelpad':10,
         'mathtext.fontset': 'stix',
         'font.family': 'STIXGeneral',
         }
plt.rcParams.update(params)


class Amplitude: 

    def __init__(self, event_collection:EventCollection) -> None:
        self.event_collection = event_collection

    def plot(self, evtid:list, file:Union[str, Path]=None, **kwargs): 
            
        ec = self.event_collection

        lid, lev = np.array(list(ec.events.keys())), np.array(list(ec.events.values()))
        
        mask = np.full(len(lid), False)            
        for m in np.argwhere([i==lid  for i in evtid]): mask[m[-1]] = True
        lid, lev = lid[mask], lev[mask]
        lsc = np.array([ev.strip_collection for ev in lev])
        sc0 = lsc[0]
        nev, nchip, nch, nsample = len(lid), sc0.nchip, sc0.nchannel, sc0.nsample

        key_var = 'Amplitude channels'
    
        fig = plt.figure(figsize=(12,10),)
        npan = nchip
        mosaic = np.array([ f'{key_var}_{ipan}' for ipan in range(npan)])
        mosaic = mosaic.reshape((2,npan//2))
        nrow, ncol = mosaic.shape[0], mosaic.shape[1]
        kwargs_size = dict(wspace=0.2, hspace=0.2, 
                    top=0.94-0.1/(nrow+1), bottom=0.3/(nrow+1), 
                    left=0.3/(ncol+1), right=1-0.3/(ncol+1) )
        axs = fig.subplot_mosaic(mosaic, sharex=True, sharey=True, gridspec_kw=kwargs_size)

        for ipan in range(npan):
            m = f'{key_var}_{ipan}'
            ax = axs[m]
           
            for iev in range(nev): 

                for ich in range(nch):
                    # try : 
                    amplitude = lsc[iev].amplitude
                    xall, yall =  np.arange(0, nsample), amplitude[ipan, ich, :]
                    nnull = (xall != 0 ) & (yall != 0 ) 
                    x, y =  xall[nnull], yall[nnull]
                    ax.plot(x, y, **kwargs)
                    # except: 
                    #     print(f"Issue plotting {iev}, {ipan}, {ich}")
                    
            schip = f"ASIC {ipan}"
            # anchored_text = AnchoredText(schip, loc="upper right", frameon=False, prop=dict(fontsize='large', color="red"))
            # ax.add_artist(anchored_text)
            ax.set_title(f"{schip}")
            ax.set_xlim([0,50])
            ax.set_ylim([0, 4096])
            
        plt.gcf().text(x=0.05/(ncol+1), y=0.5, s="Amplitude",  
                    fontsize='x-large', 
                    rotation='vertical',
                    ha='center')       
        plt.gcf().text(x=0.5, y=0.1/(nrow+1), s="Sample",  
                    fontsize='x-large', 
                    rotation='horizontal', 
                    ha='center')
                     

            # cbar.set_label(label='entries', fontsize='large')
            
        nev = len(lev)

        label = "(evtid "+" ".join([str(i) for i in evtid])+")"
        if kwargs:
            if "label" in kwargs.keys(): 
                label = kwargs["label"] + f" : {key_var} " + label 
       
        plt.gcf().text(x=0.1, y=0.985-0.1/(nrow+1), s= label,  
                    fontsize='xx-large', 
                    rotation='horizontal', 
                    color = 'red',
                    fontweight='bold', ha='left')  

        if file : 
            fig.savefig(file, dpi=300)
            print(f'Save figure {file}')
        
        plt.close()