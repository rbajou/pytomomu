#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import time

params = {'legend.fontsize': 'large',
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large',
         'axes.labelpad':10,
         'mathtext.fontset': 'stix',
         'font.family': 'STIXGeneral',
         'axes.grid' : True,
         }
plt.rcParams.update(params)

from datetime import datetime, date, timezone
from argparse import ArgumentError
from pathlib import Path
from typing import Union


class EventRate:


    def __init__(self, time:np.ndarray, t0:int=0, dt_gap:int=3600):
        self.run_duration = 0    
        self.time = np.sort(time) #s
        dtime = np.diff(self.time) 
        self.t0 = t0 #unix timestamp start 
        self.run_duration = np.sum(dtime[dtime < dt_gap])  # in second
        self.mean = 0

    def __call__(self, ax, width:float=1., label:str="",  t_off:float=0., tlim=None, **kwargs):
        # if tlim is None: tlim =  ( 0, int(datetime(2032, 4, 2, hour=16,minute=00,second=00).replace(tzinfo=timezone.utc).timestamp())   )
        # t_min, t_max = tlim
        # if t_min > t_max: raise ArgumentError("t_min > t_max")
        # mask = (t_min <= self.time) & (self.time <= t_max)
        # time = self.time[mask]
        arr_time = self.time - self.time[0]
        ntbins = int(round((np.max(arr_time)-np.min(arr_time)) / width, 1))
        entries, bins = np.histogram(arr_time,   bins=ntbins)
        center = (bins[1:]+bins[:-1])/2
        width = (bins[1:]-bins[:-1])
    
        # ntimebins = int(abs(t_end - t_start)/width) #hour
        # (self.entries, self.tbin, self.patch) = ax.hist(arr_time, bins=ntimebins, edgecolor='None', label=f"{label}\nnevts={len(arr_time):1.3e}", **kwargs)

        # if self.t0 != 0:
        center += self.t0 
        center -= self.run_duration
        t_start = int(np.min(center))
        t_end = int(np.max(center))
        date_start = datetime.fromtimestamp(t_start+t_off)
        self.start = date_start
        date_start = date(date_start.year, date_start.month, date_start.day )#
        date_end = datetime.fromtimestamp(t_end+t_off)
        self.end = date_end
        date_end = date(date_end.year, date_end.month, date_end.day )#str(datetime.fromtimestamp(data_tomo[:, 1][-1]))
        self.date_start, self.date_end=  date_start, date_end
        label +=  f"\nstart:{self.start.strftime('%y/%m/%d %H:%M')}\nend:{self.end.strftime('%y/%m/%d %H:%M')}\nduration={self.run_duration/3600:.1f}h\nnevts={len(arr_time):1.3e}"
    
        ax.bar(center, entries, width, edgecolor='None', label=label, **kwargs)
        if self.t0 != 0:
            datetime_ticks = [datetime.fromtimestamp(int(ts)).strftime('%d/%m %H:%M') for ts in ax.get_xticks()]
            ax.set_xticks(ax.get_xticks())
            ax.set_xticklabels(datetime_ticks)

        ax.set_ylabel("nevt [hour$^{-1}$]")
        ax.set_xlabel("time")
        # title =  f"Event time distribution from {str(date_start)} to {str(date_end)}"
        # ax.set_title(title)
        #plt.figtext(.5,.95, title, ha='center')
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    
    def to_csv(self, file:Union[str, Path], **kwargs):
        self.df = pd.DataFrame(data={'nevt': self.entries, 'tbinc':self.bins}, **kwargs)
        self.df.to_csv(file)


if __name__ == "__main__":

    from detector import INB72
    from analyse.analyse import Analyse
    
    t0 = time.time()
    main_path = Path.home()
    tel = INB72
    # data_path = main_path / "Data" / "tomomu" / "local_mimosa" / "DATA" / "Izen" / "COLIS_2274" / "MID" 
    data_path = main_path /"Projects" / "tomomu" / "Data" / "inb72" 
    dout = Path(__file__).parent / 'out' 
    dout.mkdir(parents = True, exist_ok = True)
    
    gen = tel.name #"MGv3"
    
    prefix = f"{tel.name}_2024"
   
    l_irun = [1, 2, 3, 4, 5, 6, 7, 12, 13, 14]#np.arange(1,15) ## #[int(arg) for arg in sys.argv[1:3]]
    iana = 1
    nev_max = int(1e7)
    fig, ax = plt.subplots(figsize=(16,9), gridspec_kw={'left':0.08, 'right':0.85, 'top':0.95})
    kwargs = {"alpha":0.3}
    for irun in l_irun:
        str_run = f"run{irun}" 
        run_path = data_path / str_run
        str_ana = f"analyse{iana}"
        filename = f"{prefix}_{str_run}_{str_ana}.root"
        basename = filename.split('.')[0]
        # fana = data_path / "1_analyse" / filename
        fana = run_path  / filename
        print(f"Load file {fana}")
        if not fana.exists(): continue
        with open(str(run_path/"timestamp0.txt"), 'r') as f : trun0 = int(f.readline())

        print(f"Import nev_max = {nev_max}")
        ana = Analyse(file = fana)
        ana.open()
        print(f"TBranches to NumPy -- {time.time() - t0:.1f} s")
        ana.get_content(n_events = nev_max) 
        evttime = ana.content["evttime"]
        label = f"Run {irun}"
        print(f"Event Rate {label}")
        er = EventRate(evttime, t0=trun0)
        print(er.time[0])
        er(ax, width=3600, label=label, **kwargs) #width = size time bin width in seconds 
 
    fout = dout / f"eventrate_runs{np.min(l_irun)}_{np.max(l_irun)}.png"
    ax.legend(loc='upper left', bbox_to_anchor=(1,1))
    plt.savefig(str(fout))
    print(f"Save figure {fout}")