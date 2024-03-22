#!/usr/bin/python3
# -*- coding: utf-8 -*-

from typing import List, Union
import numpy as np
from pathlib import Path
import uproot
import time
import pickle
import re
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.offsetbox import AnchoredText


from event.event import EventCollection

params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large',
         'axes.labelpad':10,
         'mathtext.fontset': 'stix',
         'font.family': 'STIXGeneral',
         'axes.grid': False,
         }
plt.rcParams.update(params)


class HitMap: 

    def __init__(self, event_collection:EventCollection) -> None:
        self.event_collection = event_collection
     

    def plot(self, hit_pos:np.ndarray, mask:np.ndarray, file:Union[Path,str]=None, **kwargs):

        ens = self.event_collection.detector_ensemble
        lay_pitch = [lay.pitch for det in ens.detectors for lay in det.strip_layers]
        lay_nstrips = [lay.nstrips for det in ens.detectors for lay in det.strip_layers]
       
        key_var = 'HitPos'
    
        fig = plt.figure(figsize=(10,10), )
        ndet = len(ens.detectors)
        npan = ndet
        mosaic = np.array([ f'{key_var}_xy_{ipan}' for ipan in range(npan)])
        mosaic = mosaic.reshape((2,npan//2))
        nrow, ncol = mosaic.shape[0], mosaic.shape[1]
        kwargs_size = dict(wspace=0.2, hspace=0.2, 
                    top=0.94-0.1/(nrow+1), bottom=0.3/(nrow+1), 
                    left=0.3/(ncol+1), right=1-0.3/(ncol+1) )
        axs = fig.subplot_mosaic(mosaic, sharex=True, sharey=True, gridspec_kw=kwargs_size)

        kwargs_map = {}

        for ipan in range(npan):
            
            x, y, z = hit_pos[:, ipan].T
            # m = np.full(x.shape, True)
            # if mask : m = mask[:, ipan]
            px, py = lay_pitch[2*ipan], lay_pitch[2*ipan+1]
            nx, ny = lay_nstrips[2*ipan], lay_nstrips[2*ipan+1]
            nnan = ~np.isnan(x) & ~np.isnan(y) 
            (xmin, xmax), (ymin, ymax) =  (-px*nx/2, px*nx/2), (-py*ny/2,py*ny/2)
            # print(xmin, xmax, ymin, ymax)
            bins_xy, range_xy = [150,150], [[xmin, xmax], [ymin, ymax]]
            xr, yr = np.linspace(xmin, xmax, bins_xy[0]), np.linspace(ymin, ymax, bins_xy[1])
            Xr, Yr = np.meshgrid(xr, yr)
            hxy = np.histogram2d(x[nnan], y[nnan], bins=bins_xy, range=range_xy)[0]#range=range_xy
            if np.all(hxy == 0) : continue
            # mnull = hxy==0 
            # hxy[mnull] = np.nan
            m = f'{key_var}_xy_{ipan}'
            ax = axs[m]
            vmin, vmax = 1, np.nanmax(hxy)
            # print(vmin, vmax)
            # im = ax.imshow(hxy, norm=LogNorm(vmin=vmin, vmax=vmax),**kwargs_map)
            im = ax.pcolormesh(Xr, Yr, hxy.T, norm=LogNorm(vmin=vmin, vmax=vmax),**kwargs_map)
            ax.invert_yaxis()
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = fig.colorbar(im, cax=cax, extend='max')
            
            sdet = f"X{2*ipan}/Y{2*ipan+1}"
            anchored_text = AnchoredText(f"entries: {np.nansum(hxy):.2e}", loc="upper right", frameon=True, prop=dict(fontsize="medium", color="red"))
            ax.add_artist(anchored_text)
            # ax.set_title(f"$\\mu$M {sdet}")
            ax.set_title(f"{sdet}")
            
        plt.gcf().text(x=0.05/(ncol+1), y=0.5, s="Y [mm]",  
                    fontsize='x-large', 
                    rotation='vertical',
                    ha='center')       
        plt.gcf().text(x=0.5, y=0.1/(nrow+1), s="X [mm]",  
                    fontsize='x-large', 
                    rotation='horizontal', 
                    ha='center')
                    # color = 'red',
                    # fontweight='bold')  

            # cbar.set_label(label='entries', fontsize='large')
            
        label = f"{key_var}"
        if kwargs:   
            if "label" in kwargs.keys(): label = kwargs["label"] + " : " + label 
       
        plt.gcf().text(x=0.1, y=0.985-0.1/(nrow+1), s=label,  
                    fontsize='xx-large', 
                    rotation='horizontal', 
                    color = 'red',
                    fontweight='bold')  


        fig.savefig(file, dpi=300)
        print(f'Save figure {file}')
        plt.close()
    
            

