#!/usr/bin/python3
# -*- coding: utf-8 -*-

from typing import List, Union
import numpy as np
from pathlib import Path
import time
import re
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#package module(s)
from .amplitude import Amplitude
from cluster import StripCollection
from detector import INB72
from event import Event, EventCollection
from rootfile import RootFile

import sys
# Increase the recursion limit
new_limit = 3000  # Example value
sys.setrecursionlimit(new_limit)

class Signal(RootFile):

    def __init__(self, file:Union[Path, str]) -> None:
        RootFile.__init__(self, file=file) 
        self.open()
        print('rootfile branches : ', self._branch_names)
        try : 
            regex = r'MGv(\d)'
            matches = np.array([ re.search(regex, bn) for bn in self._branch_names])
            self.basename_mg = matches[matches != None][0].group(0)
            
        except : 
            regex = r'INB72'
            matches = np.array([ re.search(regex, bn) for bn in self._branch_names])
            self.basename_mg = matches[matches != None][0].group(0)
      

        self.evvar_names = [self._branch_names[i] for i, m in enumerate(matches) if m is None ]
        self.stripbranch_names = [self._branch_names[i] for i, m in enumerate(matches) if m is not None ]
        print(f'MicroMegas version : {self.basename_mg}')
        print(f"Strip branches : {self.stripbranch_names}")

    def fill_strip_collection(self, basename:str=None):
        if len(self.content) == 0 : raise ValueError("Run 'get_content(...)' first.")
  
        self.strip_content = {k: self.content[k] for k in self._branch_names}
        lk = list(self.strip_content.keys())
    
        nev = self.strip_content[lk[0]].shape[0]

        if basename is None: basename = self.basename_mg

        self.event_collection = EventCollection()
       
        for i in range(nev):
            strip_collection = object.__new__(StripCollection)
            strip_collection.amplitude = self.strip_content[f'StripAmpl_{basename}'][i, :, :, :]
            # print(strip_collection.amplitude.shape)
            # strip_collection.amplitude = strip_collection.amplitude[range(0,2)]
            # if i > 2: 
            #    print(i, self.strip_content[f'StripAmpl_{basename}'][i].shape)
            
            event = Event()
            event.id = int(self.strip_content['Nevent'][i])
            event.time = self.strip_content['evttime'][i]
            event.strip_collection = strip_collection

            self.event_collection.__setitem__(event)
    
def update(i):
    im.set_array(image_array[i])
    return im, 

if __name__ == "__main__":
  
    t0 = time.time()
    
    main_path = Path.home()
    
    tel = INB72
    # data_path = main_path / "Data" / "tomomu" / "local_mimosa" / "DATA" / "Izen" / "COLIS_2274" / "MID" 
    mnt_path = main_path / "Projects" / "tomomu" / "inb72" / "mnt"
    
    basename = "INB72_2024"
    irun = sys.argv[1]
    # basename = "IZEN_COLIS_MID"
    str_run = f"run{irun}" 
    # str_run = f"runtest_{irun}"  
    prefix = f"{basename}_{str_run}"
    filename = f"{basename}_{str_run}_signal.root"
    basename = filename.split('.')[0]

    # fana = data_path / "1_analyse" / filename
    fsig = mnt_path  / filename
    print(f"Load file {fsig}")

    dout = Path(__file__).parent / 'out' / str_run / "signal"
    dout.mkdir(parents = True, exist_ok = True)

    sig = Signal(fsig)
    nevts = int(1e3)

    sig.get_content(n_events=nevts)
    lkey, lstrip = list(sig.content.keys()), list(sig.content.values())
   
    sig.fill_strip_collection()

    ec = sig.event_collection
    ec.detector_ensemble = tel

    levtid = np.arange(501, 511)

    lfamp = []
    kwargs = {"label": prefix, 'linewidth':0.2, 'color':'blue'}
    amplitude = Amplitude(event_collection=ec)
    for evtid in levtid:
        filename = f"amplitude_channel_{ec.detector_ensemble.name}_"+ basename + f'_evt{evtid}.png'
        famp = dout / filename
        lfamp.append(famp)
        amplitude.plot(evtid=[evtid], file=famp, **kwargs)


    image_array = []
    for my_file in lfamp:
        
        image = Image.open(str(my_file))
        image_array.append(image)


  

    # Create the figure and axes objects
    fig, ax = plt.subplots(figsize=(12,10),)

    # Set the initial image
    im = ax.imshow(image_array[0], animated=True)



    # Create the animation object
    animation_fig = animation.FuncAnimation(fig, update, frames=len(image_array), interval=1000, blit=True,repeat_delay=10,)
    plt.axis('off')
    # Show the animation
    # plt.show()

    animation_fig.save(dout/"amplitude_channels.gif", dpi=300)  
    print(f"Save gif {dout/'amplitude_channels.gif'}")