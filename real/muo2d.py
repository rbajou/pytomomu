#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from typing import List, Union
import numpy as np
from pathlib import Path
import uproot
import time
import pickle
import re
import matplotlib.pyplot as plt

from analyse import Analyse, Final
from cluster import ClusterCollection
from detector import DetectorEnsemble, Bench, VB, IZEN
from event import Event, EventCollection
from hitmap import HitMap


class Muo2D:

    # def __init__(self, event_collection :EventCollection) -> None:
    #     tracks = [ev.track for ev in event_collection]
    #     ttheta, tthetax = []

    def plot(self, file:Union[Path,str]=None):


        fig, ax = plt.Subplots()

        Z = np.histogram2d(x, y )

        ax.pcolormesh(X, Y, Z)

        if file:
            fig.savefig(file, dpi=300)
            print(f'Save figure {file}')
        plt.close()
    

if __name__ == "__main__":
    pass

    t0 = time.time()

    home = Path.home()
    data_path = home / "Data" / "tomomu" / "inb72" 
    ana_path = data_path / "analyse" 
    final_path = data_path / "final"  
    final_path.mkdir(parents=True, exist_ok=True)

    n_events_import = 100000

    ####### Change file name here 
    irun, iana= 1, 1 #[int(arg) for arg in sys.argv[1:3]] 
    prefix = "INB72_2024"
    str_run = f"run{irun}" # 
    str_ana = f"analyse{iana}"
    
    # filename = f"{prefix}_{str_run}_{str_ana}.root"
    filename = f"{prefix}_{str_run}_final.root"
    basename = filename.split('.')[0]
    ffinal = final_path / filename

    final = Final(ffinal)
    final.open()
    print(final.branch_names)
    
    # x, y = final
    # fig, ax = plt.subplots()