'''
Adapted from Paul's script
'''

import numpy as np
import matplotlib
import tkinter
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from pathlib import Path 
import time

#tomomuLib module(s)
from config_INB72 import def_INB72

import uproot as up
from rootfile import general_import


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
filename = f"{prefix}_{str_run}_{str_ana}.root"
basename = filename.split('.')[0]
fana = ana_path / filename


#Importation of the data

#If data real measured data
n_events = int(1e4)
print(f"Import n_events = {n_events}")
gen = "INB72"
branch_list = [f'{gen}_ClusPos', f'{gen}_ClusSize', 'evn','evttime']
dat = general_import(file=fana, branch_list = branch_list, n_events=n_events)
 
print(f"TBranches to NumPy -- {time.time() - t0:.1f} s")

arr_pos = dat[f'{gen}_ClusPos']
arr_size = dat[f'{gen}_ClusSize']
event = np.asarray(dat['evn'], dtype=np.int32)
n_events = event.shape[0]


print("data imported !")

 
tel = def_INB72(arr_pos, arr_size)#, unaligned=False)




print('Now fitting the events...')
#Definition of resx and resy from the fits
[ax, bx, rx], [ay, by, ry] = tel.fits()


ind_resx = np.argmin(np.abs(rx), axis=0)
ind_resy = np.argmin(np.abs(ry), axis=0)

resx = np.choose(ind_resx, rx)
resy = np.choose(ind_resx, ry)

ntrack = np.zeros(n_events, dtype=np.int32)
quarter = np.zeros(n_events, dtype=np.int32)
etime = np.zeros(n_events, dtype=np.float64)


X0 = bx + ax * tel.plans[0].positions_center[0,2]
Y0 = by + ay * tel.plans[0].positions_center[0,2]

ttheta = ax
tthetax = ay


fout = final_path / f"{prefix}_{str_run}_final.root"
#Datas in the right type
print(f'Write file {fout}...')
resx = resx.astype(np.float64)
resy = resy.astype(np.float64)
X0 = X0.astype(np.float64)
Y0 = Y0.astype(np.float64)
ttheta = ttheta.astype(np.float64)
tthetax = tthetax.astype(np.float64)

root_out = up.recreate(str(fout))

root_out['T'] = {'event': event,
                 'etime': etime,
                 'ntrack': ntrack,
                 'quarter': quarter,
                 'resx': resx,
                 'resy': resy,
                 'X0': X0,
                 'Y0': Y0,
                 'tthetay': ttheta,
                 'tthetax': tthetax,
                 }

root_out.close()
