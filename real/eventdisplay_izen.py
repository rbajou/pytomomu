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
from detector import IZEN
from utils import tools



if __name__ == "__main__":

    t0 = time.time()
    main_path = Path(__file__).parents[2]
    tel = IZEN
    print(main_path)
    # data_path = main_path / "Data" / "tomomu" / "local_mimosa" / "DATA" / "Izen" / "COLIS_2274" / "MID" 
    data_path = main_path / "Data" / "izen" 
    irun, iana= 12, 1 #[int(arg) for arg in sys.argv[1:3]]
    prefix = "IZEN_COLIS_MID"
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
    clusvar_prefix = 'MGv3'
    lclus_var = [f'{clusvar_prefix}_ClusPos', f'{clusvar_prefix}_ClusAmpl', f'{clusvar_prefix}_ClusMaxStripAmpl', f'{clusvar_prefix}_ClusSize', f'{clusvar_prefix}_ClusTOT', f'{clusvar_prefix}_ClusMaxSample', f'{clusvar_prefix}_ClusT']
    branch_list = lclus_var
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
    
    
    iev = 3
    for iev in range(1,2):
        ev = ec.events[iev]
        nc = ev.cluster_collection.nclus
        xyz0 = ev.xyz
        cl_ampl = ev.cluster_collection.amplitude
        cl_size = ev.cluster_collection.size #in number of strips 
        ec.get_tracks(ec.xyz[:, -4:, :])
        xyz_t = np.array([ev.track.intersection(z) for z in np.linspace(0, 4500, 50)])
        # print(f"xyz = {xyz0}")
        npan, nclus = 12, 8
        ampl = cl_ampl.reshape(npan,2,nclus)
        ampl_sum = np.sum(ampl, axis=1)
        cl_sum = np.sum(cl_size.reshape(npan,2,nclus), axis=1)
        icl = range(0,8)   
        cmap = palettable.matplotlib.Viridis_20.mpl_colormap 
        vmin, vmax, n = 0, np.max(ampl_sum[:, icl]), 100
        range_val = np.linspace(vmin, vmax, n)
        norm = cm.Normalize(vmin=vmin, vmax=vmax)
        color_scale =  cmap(norm(range_val))
        m = np.all(xyz0[:, icl, :] > -999., axis=2)
        arg_col =  [np.argmin(abs(range_val-v))for v in ampl_sum[:,icl][m].flatten()]   
        color_values = color_scale[arg_col]#.reshape(m.shape[0], m.shape[1], 4)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        kwargs = dict(alpha=0.2, color='greenyellow', edgecolor='none')
        ec.detector_ensemble.plot3D(fig=fig, ax=ax, **kwargs)
        scatter = ax.scatter(
                xyz0[:, icl, 0][m],
                xyz0[:, icl, 1][m],
                xyz0[:, icl, 2][m],
                s=cl_sum[:, icl][m],
                c= color_values,#'limegreen',
                marker='o',
                # edgecolor=color_values,
        )
        ax.plot(xyz_t[:,0],xyz_t[:,1], xyz_t[:,-1],
			c="blue", linewidth=0.75)

        cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, pad=0.1)
        cbar.set_label('Amplitude')
        
        leg_hdl= []
        for s in [5, 10, 20]:
            hdl = mlines.Line2D([], [], marker='o', fillstyle='none', markeredgecolor='grey', linestyle='None', alpha=0.3, markersize=s, label=f"{s}")
            leg_hdl.append(hdl)
        ax.legend(title='Cluster size', handles = leg_hdl, loc='best',fontsize='large', bbox_to_anchor=(1,1))
        ax.set_xlim(-1000, 500)
        ax.set_ylim(-750, 500)
        ax.set_zlim(0, 4500)
        plt.show()
        plt.close()
 