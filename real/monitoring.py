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
import glob
from datetime import datetime
import sys
import warnings
warnings.filterwarnings("ignore")
#package module(s)
from detector import INB72
from real.rootfile import RootFile



params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large',
         'axes.labelpad':10,
         'mathtext.fontset': 'stix',
         'font.family': 'STIXGeneral',
         'axes.grid' : True 
         }
plt.rcParams.update(params)

class Monitoring(RootFile):

    def __init__(self, file:Union[Path, str]) -> None:
        RootFile.__init__(self, file=file)
        self.open(treename = "HV_Mon")
        print('rootfile branches : ', self._branch_names)

if __name__ == "__main__":

    t0 = time.time()
    main_path = Path.home()
    tel = INB72
    data_path = main_path / "Data" / "tomomu" / "inb72" 
    #INB72_monitoring_20240122_12H21_HV_mon.root
    irun = sys.argv[1]
    str_run = f"run{irun}" 

    run_path = data_path / str_run

    prefix = f"INB72_2024_{str_run}"

    datestr = "20240206_09h37"
    fmon = sys.argv[2] if len(sys.argv > 2) else glob.glob(str(main_path / "Projects" / "tomomu" / "Data" / str_run / "*monitoring*"))[0] #"/local/home/rb277364/Projects/tomomu/inb72/mnt/INB72_monitoring_20240208_12H15_HV_mon.root"
    
    mon = Monitoring(file= fmon) 
    print(f"Load file {fmon}")

    nev=1e7

    mon.get_content(n_events=nev)

    timestamp = mon.content['time']

    # ampl = mon.content['ampl']
    flow = mon.content['flow']

    dout = data_path / "monitoring" 
    
    dout.mkdir(parents = True, exist_ok = True)
    lvmon = ['flowIn', 'flowOut']
    dict_monflow = {'flowIn [L/h]':flow[:,0,0], 'flowOut [L/h]': flow[:,0,1]}
    nvar = len(dict_monflow)
    nrow, ncol = nvar, 1 
    fig = plt.figure(figsize=(16,11), )
    left, right, top, bottom = 0.2/(ncol+1), 0.95-0.1/(ncol+1), 0.94-0.1/(nrow+1), 0.4/(nrow+1)
    kwargs_size = dict(wspace=0.0, hspace=0.1, 
                top=top, bottom=bottom, 
                left=left, right=right)
    gs = fig.add_gridspec(nrow, ncol, **kwargs_size )
    kwargs = {}
    step = 1
    for ivar, (key, var) in enumerate(dict_monflow.items()):
        ax = fig.add_subplot(gs[ivar, 0])
        ax.plot(timestamp[::step], var[::step], **kwargs)
        ax.set_ylabel(f'{key}')
        if ivar != nvar-1 : ax.set_xticklabels([])
        else : 
            datetime_ticks = [datetime.fromtimestamp(int(ts)).strftime('%d/%m %H:%M') for ts in ax.get_xticks()]
            ax.set_xticklabels(datetime_ticks)
            ax.set_xlabel('time')   
            for label in ax.get_xticklabels(which='major'):
                label.set(rotation=30, horizontalalignment='right')

    plt.gcf().text(x=left, y=0.985-0.1/(nrow+1), s=f"{prefix} : Monitoring gas flow",  
                fontsize='xx-large', 
                rotation='horizontal', 
                color = 'red',
                fontweight='bold')  

    foutname = f"mosaic_monitoring_flow_{prefix}.png"
    fout = dout / foutname
    fig.savefig(fout, dpi=300)
    print(f'Save figure {fout}')
    plt.close()

    
    #TENSIONS
    HVPS_5V  = mon.content['HVPS_5V']
    HVPS_15V  = mon.content['HVPS_15V']
    HVPS_3V3  = mon.content['HVPS_3V3']
    Vmon  = mon.content['V'] #(nent, 1, 5) ?
    Iset = mon.content['I']
    feedback_type = mon.content['feedback_type']
    # print("voltage var : ", HVPS_5V.shape, HVPS_15V.shape, HVPS_3V3.shape, Vmon.shape, Iset.shape, feedback_type.shape)

    dict_montens = {'HVPS_5V [V]':HVPS_5V, 'HVPS_15V [V]': HVPS_15V, 'HVPS_3V3 [V]' : HVPS_3V3, 'feedback_type': feedback_type}
    nvar = len(dict_montens)
    nrow, ncol = nvar, 1 
    fig = plt.figure(figsize=(16,11), )
    left, right, top, bottom = 0.2/(ncol+1), 0.95-0.1/(ncol+1), 0.94-0.1/(nrow+1), 0.4/(nrow+1)
    kwargs_size = dict(wspace=0.0, hspace=0.1, 
                top=top, bottom=bottom, 
                left=left, right=right)
    gs = fig.add_gridspec(nrow, ncol, **kwargs_size )
    kwargs = {}
    for ivar, (key, var) in enumerate(dict_montens.items()):
        ax = fig.add_subplot(gs[ivar, 0])
        ax.plot(timestamp[::step], var[::step], **kwargs)
        ax.set_ylabel(f'{key}')
        if ivar != nvar-1 : ax.set_xticklabels([])
        else : 
            datetime_ticks = [datetime.fromtimestamp(int(ts)).strftime('%d/%m %H:%M') for ts in ax.get_xticks()]
            ax.set_xticklabels(datetime_ticks)
            ax.set_xlabel('time')
            for label in ax.get_xticklabels(which='major'):
                label.set(rotation=30, horizontalalignment='right')
          
    plt.gcf().text(x=left, y=0.985-0.1/(nrow+1), s=f"{prefix} : HV",  
                fontsize='xx-large', 
                rotation='horizontal', 
                color = 'red',
                fontweight='bold')  

    foutname = f"mosaic_HV_{prefix}.png"
    fout = dout / foutname
    fig.savefig(fout, dpi=300)
    print(f'Save figure {fout}')
    plt.close()

    feedback_type = mon.content['feedback_type']
    dict_montens = {'Rstrip_1 [V]':Vmon[:,:,0], 'Rstrip_2 [V]':Vmon[:,:,1], 'Rstrip_3 [V]':Vmon[:,:,2], 'Drift [V]':Vmon[:,:,3], 'Rstrip4 [V]':Vmon[:,:,4]}
    nvar = len(dict_montens)
    nrow, ncol = nvar, 1 
    fig = plt.figure(figsize=(16,11), )
    left, right, top, bottom = 0.2/(ncol+1), 0.95-0.1/(ncol+1), 0.94-0.1/(nrow+1), 0.4/(nrow+1)
    kwargs_size = dict(wspace=0.0, hspace=0.1, 
                top=top, bottom=bottom, 
                left=left, right=right)
    gs = fig.add_gridspec(nrow, ncol, **kwargs_size )
    kwargs = {}
    for ivar, (key, var) in enumerate(dict_montens.items()):
        ax = fig.add_subplot(gs[ivar, 0])
        ax.plot(timestamp[::step], var[::step], **kwargs)
        ax.set_ylabel(f'{key}')
        if ivar != nvar-1 : ax.set_xticklabels([])
        else : 
            datetime_ticks = [datetime.fromtimestamp(int(ts)).strftime('%d/%m %H:%M') for ts in ax.get_xticks()]
            ax.set_xticklabels(datetime_ticks)
            ax.set_xlabel('time')
            for label in ax.get_xticklabels(which='major'):
                label.set(rotation=30, horizontalalignment='right')
          
    plt.gcf().text(x=left, y=0.985-0.1/(nrow+1), s=f"{prefix} : Vmon",  
                fontsize='xx-large', 
                rotation='horizontal', 
                color = 'red',
                fontweight='bold')  

    foutname = f"mosaic_Vmon_{prefix}.png"
    fout = dout / foutname
    fig.savefig(fout, dpi=300)
    print(f'Save figure {fout}')
    plt.close()


    #T/P
    T_PC = mon.content['T_PC']
    T_HV = mon.content['T_HV']
    # print(T_PC.shape, T_HV.shape)

    dict_montemp = {'T_PC [°C]':T_PC, 'T_HV [°C]': T_HV}
    nvar = len(dict_montemp)
    nrow, ncol = nvar, 1 
    fig = plt.figure(figsize=(16,11), )
    left, right, top, bottom = 0.2/(ncol+1), 0.95-0.1/(ncol+1), 0.94-0.1/(nrow+1), 0.4/(nrow+1)
    kwargs_size = dict(wspace=0.0, hspace=0.1, 
                top=top, bottom=bottom, 
                left=left, right=right)
    gs = fig.add_gridspec(nrow, ncol, **kwargs_size )
    kwargs = {}
    for ivar, (key, var) in enumerate(dict_montemp.items()):
        ax = fig.add_subplot(gs[ivar, 0])
        ax.plot(timestamp[::step], var[::step], **kwargs)
        ax.set_ylabel(f'{key}')
        if ivar != nvar-1 : ax.set_xticklabels([])
        else : 
            datetime_ticks = [datetime.fromtimestamp(int(ts)).strftime('%d/%m %H:%M') for ts in ax.get_xticks()]
            ax.set_xticklabels(datetime_ticks)
            ax.set_xlabel('time')
            for label in ax.get_xticklabels(which='major'):
                label.set(rotation=30, horizontalalignment='right')
          
    plt.gcf().text(x=left, y=0.985-0.1/(nrow+1), s=f"{prefix} : Monitoring temperature",  
                fontsize='xx-large', 
                rotation='horizontal', 
                color = 'red',
                fontweight='bold')  

    foutname = f"mosaic_monitoring_temp_{prefix}.png"
    fout = dout / foutname
    fig.savefig(fout, dpi=300)
    print(f'Save figure {fout}')
    plt.close()




