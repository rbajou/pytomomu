#!/usr/bin/python3
# -*- coding: utf-8 -*-

import subprocess
import numpy as np
from pathlib import Path
import time
import matplotlib.pyplot as plt
import argparse
from typing import Union
from matplotlib.cm import ScalarMappable
import matplotlib.colors as cm
import matplotlib.lines as mlines  #use for legend settings
import palettable
import random
import seaborn as sns
#package module(s)
from analyse import Analyse, Final
from detector import Bench, Telescope
from survey import set_real_survey

params = {'legend.fontsize': 'xx-large',
          'legend.title_fontsize' : 'xx-large',
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


t0 = time.time()

parser=argparse.ArgumentParser(
description='''Main script for processing real data *analysis*.root file, reconstruct muon tracks and save output in *final*.root file.''', epilog='''All is well that ends well.''')
parser.add_argument('--survey', '-s', default='izen', help="Survey name",  type=str)
parser.add_argument('--run_num', '-r', default=4, help="Run number",  type=int)
parser.add_argument('--ntracks', '-n', default=int(1e3), help="Number of event to display", type=int)
parser.add_argument('--nevt_max', '-max', default=int(1e3), help="Number of events to load", type=int)
parser.add_argument('--ana', '-a', default=1, help="Analysis file '1' or '2'", type=Union[int, str])
args=parser.parse_args()

main_path = Path.home() / "Projects" / "tomomu"
CURRENT_SURVEY = set_real_survey(args.survey)
det = CURRENT_SURVEY.detector 

if args.ntracks > args.nevt_max : args.nevt_max = args.ntracks

if isinstance(det, Telescope): tel = det
elif isinstance(det, Bench) : tel =  det.telescope[0]
else : raise TypeError("Unknown detector type")

data_path = Path(CURRENT_SURVEY.data_path) #main_path / "Data" / "inb72" 
print(f"Data path : {data_path}")

run_prefix = CURRENT_SURVEY.run_prefix 

lrun  = [args.run_num] if args.run_num else CURRENT_SURVEY.run_num 
n_events = args.nevt_max 

clusvar_prefix = CURRENT_SURVEY.clusvar_prefix
iana=  args.ana
str_ana = f"analyse{iana}"
ana_path = data_path / "1_analyse"
final_path =  data_path / "2_final" 
mon_path = data_path / "3_monitoring"
out_path = data_path / "9_out" 
out_path.mkdir(parents=True, exist_ok=True)

irun = args.run_num


str_run = f"run{irun}" 

label = f"{run_prefix}_run{irun}"
filename = f"{run_prefix}_{str_run}_{str_ana}.root"
print(f"\nProcess {filename.split('.')[0]}\t. . .")
basename = filename.split('.')[0]

rout = out_path / str_run
rout.mkdir(parents = True, exist_ok = True)


fana = ana_path  / filename
if not fana.exists():
    try : 
        short_ana1_script = Path(__file__).parents[1] / "macro" / f"shorten_analysis1.sh"
        stdout = subprocess.run(["bash", str(short_ana1_script), det.name.lower(), run_prefix, str(irun)], check=True, stdout=subprocess.PIPE).stdout.decode('utf-8')
        while True:
            line = stdout.readline()
            if not line: break
        print(stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error : {e} \n Failed to convert '{run_prefix}_{str_run}_analyse.root' to '{run_prefix}_{str_run}_{str_ana}.root'.")

print(f"\nLoad file {fana}")

lclus_var = [f'{clusvar_prefix}_ClusPos', f'{clusvar_prefix}_ClusAmpl', f'{clusvar_prefix}_ClusMaxStripAmpl', f'{clusvar_prefix}_ClusSize', f'{clusvar_prefix}_ClusTOT', f'{clusvar_prefix}_ClusMaxSample', f'{clusvar_prefix}_ClusT']
branch_list = lclus_var

ana = Analyse(file = fana)
n_events_tot = ana.rootfile.num_entries
print(f"In file n_events_tot = {n_events_tot}" )
n_events = n_events if n_events < n_events_tot else n_events_tot
print(f"Import n_events = {n_events}")
ana.open(treename="T")
print(f"TBranches to NumPy -- {time.time() - t0:.1f} s\n")

ana.get_content(n_events = n_events) 
ana.fill_event_collection()
ec = ana.event_collection
ec.detector_ensemble = det
mask = np.full(len(ec.events), True) 
# event_no = np.array([ev.id for i, (_, ev) in enumerate(ec.events.items())])
# event_time = np.array([ev.time for i, (_, ev) in enumerate(ec.events.items())])
clus_size = np.array([ev.cluster_collection.size for i, (_, ev) in enumerate(ec.events.items())])# if mask[i] == True ])
mask_clus_size = (clus_size > 1) #np.full(clus_size.shape, True)  #
ec.get_xyz(mask=mask_clus_size)

liev = list(ec.events.keys())
riev = np.random.randint(min(liev), max(liev), args.ntracks)
str_riev = f"{'_'.join([str(i) for i in np.sort(riev)])}"
# col_riev = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
            #  for i in range(args.ntracks)]
col_riev = sns.color_palette("Spectral_r", args.ntracks).as_hex()

fig, axs = plt.subplots(figsize=(12, 10), nrows=1, ncols=2)


kwargs = dict(alpha=0.2, color='greenyellow', edgecolor='none')

##colormap & colorbar
cmap = palettable.matplotlib.Viridis_20.mpl_colormap 
# vmin, vmax, n = 0, np.max(ampl_sum[:, icl]), 100
vmin, vmax, n = 0, 4800, 100 #amplitude range in ADC
range_val = np.linspace(vmin, vmax, n)
norm = cm.Normalize(vmin=vmin, vmax=vmax)
color_scale =  cmap(norm(range_val))

fig, ax = plt.subplots(figsize=(12,8))

#plot telescope or bench
if isinstance(det, Telescope): 
    det.plotXZ(ax=ax, **kwargs) 
elif isinstance(det, Bench) : 
    kwargs = dict(alpha=1., color='greenyellow')
    for d in det.detector_ensemble: 
        print(d.name)
        d.plotXZ(ax=ax, **kwargs) 
        d.plotXZ(ax=ax, **kwargs) 
        d.plotXZ(ax=ax, **kwargs) 
else : raise TypeError("Unknown detector type")

xyz_tel = ec.xyz[:, -4:, :]
ec.get_tracks(xyz_tel) 


for i, iev in enumerate(riev) : 

    print(f"Plot evt {iev}")
    # print(f'mask_clus_size = {mask_clus_size[iev-1]}, {mask_clus_size[iev-1].shape}')

    ev = ec.events[iev]
    nc = ev.cluster_collection.nclus
    # print(nc)
    xyz = ev.xyz
    cl_ampl = ev.cluster_collection.amplitude
    cl_size = ev.cluster_collection.size #in number of strips 


    zmin, zmax= det.detectors[0].position[-1] - 100, det.detectors[-1].position[-1] + 100

    npan, nclus = len(det.detectors), ec.xyz.shape[2]
    ampl = cl_ampl.reshape(npan,2,nclus)
    ampl_sum = np.sum(ampl, axis=1)
    cl_sum = np.sum(cl_size.reshape(npan,2,nclus), axis=1)
    icl = range(0,8)   

 
    m = np.all(xyz[:, icl, :] > -999., axis=2)
    arg_col =  [np.argmin(abs(range_val-v))for v in ampl_sum[:,icl][m].flatten()]   
    color_values = color_scale[arg_col]#.reshape(m.shape[0], m.shape[1], 4)


    scatter = ax.scatter(
            xyz[:, icl, 0][m],
            xyz[:, icl, 2][m],
            s=cl_sum[:, icl][m],
            c= color_values,#'limegreen',
            marker='o',
            # edgecolor='darkgrey',
    )
    #plot track :
    # if mask_3p[iev] :
    xyz_trk = np.array([ev.track.intersection(z) for z in np.linspace(zmin, zmax, 50)])
    ax.plot(xyz_trk[:,0], xyz_trk[:,-1],
        c=col_riev[i], linewidth=0.75)
    fontdict= {"size":"xx-large", "color": col_riev[i], "weight":"bold", "ha":"center"}
    ax.text(xyz_trk[-1,0],xyz_trk[-1,-1], s=iev, fontdict=fontdict)
    ax.text(xyz_trk[0,0],xyz_trk[0,-1], s=iev, fontdict=fontdict)


##cluster size legend 
leg_hdl= []
for s in [5, 10, 20]:
    hdl = mlines.Line2D([], [], marker='o', fillstyle='none', markeredgecolor='black', linestyle='None', alpha=0.3, markersize=s, label=f"{s}")
    leg_hdl.append(hdl)
ax.legend(title='Cluster size', handles = leg_hdl, loc='upper left')#, bbox_to_anchor=(0.01,1.1))

#colorbar
cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, pad=0.05, shrink=1)
cbar.set_label('Amplitude [ADC]')

det_side = tel.detectors[0].side
ax.set_xlim(-det_side, det_side)
ax.set_ylim(zmin-100, zmax+100)
ax.set_title(f"Event(s) {', '.join([str(i) for i in np.sort(riev)])}", weight='bold', )
plt.tight_layout()

ax.set_xlim(-500, 500)
ax.set_ylim(0, 4500)
fout = f"{tel.name}_xz.png"
fig.savefig(fout)
# plt.show(block=True)
print(f'Save {fout}')