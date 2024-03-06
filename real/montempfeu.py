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


if __name__=="__main__":

    prefix = "INB72_2024"
    str_run = "run13"
    file = Path.home() / f"Data/tomomu/inb72/run13/TEST_{str_run}.txt"
    dout = Path(__file__).parent / 'out' / str_run 
    dout.mkdir(parents=True, exist_ok=True)
    Max_Tmp_Int, Max_Tmp_SideX, Max_Tmp_SideA = [], [], []
    dict_temp = {'Max_Tmp_Int': [], 'Max_Tmp_SideX': [], 'Max_Tmp_SideA': []}
    nl=0
    nval = 0
    with open(str(file), 'r') as f :
        lines= f.readlines()
        for l in lines:
            nl+=1
            lval = np.asarray(l.split("  "))
            if l.startswith(' Max_Tmp_Int'): 
                arg = np.argwhere(lval=='0 Avr =')[0][0]
                val = float(lval[arg+1])*0.1
                dict_temp['Max_Tmp_Int'].append(val)
                pass
            elif l.startswith(' Max_Tmp_SideX'):
                arg = np.argwhere(lval=='0 Avr =')[0][0]
                val = float(lval[arg+1])*0.1
                dict_temp['Max_Tmp_SideX'].append(val)
                pass
            elif l.startswith(' Max_Tmp_SideA'):
                arg = np.argwhere(lval=='0 Avr =')[0][0]
                val = float(lval[arg+1])*0.1
                dict_temp['Max_Tmp_SideA'].append(val)
                nval +=1
                pass
            else:
                pass
    
    tstep = 2
    timestamp = np.arange(0, nval) *tstep
    nstep = 1
    print(len(timestamp), len(dict_temp['Max_Tmp_SideA']), len(dict_temp['Max_Tmp_SideX']))
    nvar = len(dict_temp)
    nrow, ncol = nvar, 1 
    fig = plt.figure(figsize=(16,11), )
    wspace, hspace = 0., 0.1
    left, right, top, bottom = 0.2/(ncol+1), 0.95-0.1/(ncol+1), 0.94-0.1/(nrow+1), 0.4/(nrow+1)
    kwargs_size = dict(wspace=wspace, hspace=hspace, 
                top=top, bottom=bottom, 
                left=left, right=right)
    gs = fig.add_gridspec(nrow, ncol, **kwargs_size )
    kwargs = {}
    for ivar, (key, var) in enumerate(dict_temp.items()):
        ax = fig.add_subplot(gs[ivar, 0])
        ax.plot(timestamp[::nstep], var[::nstep], **kwargs)
        ax.set_ylabel(f'{key}')
        #     datetime_ticks = [datetime.fromtimestamp(int(ts)).strftime('%d/%m %H:%M') for ts in ax.get_xticks()]
        #     ax.set_xticklabels(datetime_ticks)
        xval = np.arange(min(timestamp), max(timestamp), 3600)
        ax.set_xticks(xval)
        new_labels = np.array([v/3600 for v in xval])
        ax.set_xticklabels(new_labels)
        if ivar != nvar-1 : ax.set_xticklabels([])
        else : 
            ax.set_xlabel('time [h]')
          
    plt.gcf().text(x=left, y=0.985-0.1/(nrow+1), s=f"SlowControl FEU : Temp [°C]    ({prefix}_{str_run})",  
                fontsize='xx-large', 
                rotation='horizontal', 
                color = 'red',
                fontweight='bold')  
    foutname = f"mon_temp_FEU.png"
    fout = dout / foutname
    fig.savefig(fout, dpi=300)
    print(f'Save figure {fout}')

    plt.close()