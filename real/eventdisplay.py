#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from matplotlib.cm import ScalarMappable
import matplotlib.colors as cm
import matplotlib.lines as mlines  #use for legend settings
# from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import palettable
from random import sample, choice
import time

from analyse import Analyse
from detector import INB72
from utils import tools



if __name__ == "__main__":

    t0 = time.time()
    main_path = Path(__file__).parents[2]
    tel = INB72
    print(main_path)
    # data_path = main_path / "Data" / "tomomu" / "local_mimosa" / "DATA" / "Izen" / "COLIS_2274" / "MID" 
    data_path = main_path / "Data" / "inb72" 
    irun, iana= 1, 1 #[int(arg) for arg in sys.argv[1:3]]
    prefix = "INB72_2024"
    # prefix = "IZEN_COLIS_MID"
    str_run = f"run{irun}" 
    ana_path = data_path / str_run #"analyse"
    # str_run = f"runtest_{irun}"  
    str_ana = f"analyse{iana}"
    filename = f"{prefix}_{str_run}_{str_ana}.root"
    basename = filename.split('.')[0]
    # fana = data_path / "1_analyse" / filename
    fana = ana_path  / filename
    print(f"Load file {fana}")
    dout = Path(__file__).parent / 'out' / str_run / str_ana
    dout.mkdir(parents = True, exist_ok = True)
    # gen = 'MGv3'
    gen = "INB72"
    # ltime_var = ['time', 'deltat']
    lclus_var = [f'{gen}_ClusPos', f'{gen}_ClusAmpl', f'{gen}_ClusMaxStripAmpl', f'{gen}_ClusSize', f'{gen}_ClusTOT', f'{gen}_ClusMaxSample', f'{gen}_ClusT']
    # dict_time_var = {k: np.asarray(tree[k]) for k in ltime_var}
    branch_list = lclus_var
    # n_events = int(sys.argv[3]) if (len(sys.argv) >= 4) else int(1e5)
    n_events = int(1e2)
    print(f"Import n_events = {n_events}")
    ana = Analyse(file = fana)
    ana.open()
    print(f"TBranches to NumPy -- {time.time() - t0:.1f} s")
    ana.get_content(n_events = n_events) 
    ana.fill_event_collection()
    ec = ana.event_collection
    ec.detector_ensemble = tel
    ec.get_xyz()
    ec.get_tracks(ec.xyz)
    
    m0 = ec.mask_xyz[0]
    iev = 1
    ev = ec.events[iev]
    nc = ev.cluster_collection.nclus
    xyz0 = ev.xyz
    cl_ampl = ev.cluster_collection.amplitude
    cl_size = ev.cluster_collection.size #in number of strips 
    print(np.min(cl_size), np.mean(cl_size), np.max(cl_size))
    xyz_t = np.array([ev.track.intersection(z) for z in np.linspace(0, 180, 50)])

    ampl = cl_ampl.reshape(4,2,8)
    ampl_sum = np.sum(ampl, axis=1)
    cl_sum = np.sum(cl_size.reshape(4,2,8), axis=1)
    print(ampl_sum.shape)
    print(f"cluster_size shape = {cl_size.shape}")
    print(f"cl_sum shape = {cl_sum.shape}")

    icl = range(0,8)   
    # print(palettable.scientific.sequential.print_maps())
    cmap = palettable.matplotlib.Inferno_20_r.mpl_colormap 
    vmin, vmax, n = 0, np.max(ampl_sum[:, icl]), 100
    range_val = np.linspace(vmin, vmax, n)
    norm = cm.Normalize(vmin=vmin, vmax=vmax, )
    color_scale =  cmap(norm(range_val))
    print(type(norm))
    print(type(color_scale))
    m = xyz0[:, icl, 0] > 0
    arg_col =  [np.argmin(abs(range_val-v))for v in ampl_sum[:,icl][m].flatten()]   
    print("arg_col = ", len(arg_col)) 
    color_values = color_scale[arg_col]#.reshape(m.shape[0], m.shape[1], 4)
    print("color_values.shape = ",color_values.shape) 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
#     ax = Axes3D(fig)
    ec.detector_ensemble.plot3D(fig=fig, ax=ax)
    scatter = ax.scatter(
            xyz0[:, icl, 0][m],
            xyz0[:, icl, 1][m],
            xyz0[:, icl, 2][m],
            s=cl_sum[:, icl][m],
            c= color_values,#'limegreen',
            marker='o',
            # edgecolor=color_values,
    )
    leg_hdl= []
    for s in [5, 10, 20]:
        hdl = mlines.Line2D([], [], marker='o', fillstyle='none', markeredgecolor='grey', linestyle='None', alpha=0.3, markersize=s, label=f"{s}")
        leg_hdl.append(hdl)
    ax.legend(title='Cluster size', handles = leg_hdl, loc='best',fontsize='large', bbox_to_anchor=(1,1))
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, pad=0.1)
    cbar.set_label('Amplitude')
    
#     # ax.legend(handles=handles, fontsize='large',  loc=(0.5,0.85))#'upper right')#, pad=-10)
#     ax.view_init(elev=10., azim=-60)
#     ax.plot(xyz_t[:,0],xyz_t[:,1], xyz_t[:,-1],
#             c="blue", linewidth=0.75)
    plt.show()
 



    plt.close()

    # evtID_good= sample(list(df.index), N) 
    # print(f"evtID = {evtID_good}")