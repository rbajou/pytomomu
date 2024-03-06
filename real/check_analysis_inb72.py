#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from typing import Union
from pathlib import Path
import uproot 
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.offsetbox import AnchoredText
import time
import sys

#package module(s)
from rootfile import general_import
from detector import INB72


params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'axes.labelpad':10,
         'mathtext.fontset': 'stix',
         'font.family': 'STIXGeneral',
         }
plt.rcParams.update(params)




if __name__ == "__main__": 
   
    #avoid running the script for 

    t0 = time.time()
    
    main_path = Path(__file__).parents[2]
    
    # data_path = main_path / "Data" / "tomomu" / "local_mimosa" / "DATA" / "Izen" / "COLIS_2274" / "MID" 
    data_path = main_path / "Data" / "inb72" 
    irun, iana= [int(arg) for arg in sys.argv[1:3]]
    str_run = f"run{irun}" # 
    run_path = data_path / str_run
    
    # str_run = f"runtest_{irun}"  #
    str_ana = f"analyse{iana}"
    # filename = f"IZEN_COLIS_MID_{str_run}_{str_ana}.root"
    filename = f"INB72_2024_{str_run}_{str_ana}.root"
    #issue with Data/tomomu/local_mimosa/DATA/Izen/COLIS_2274/MID/IZEN_COLIS_MID_run0_analyse.root
    #It causes the pc to suddenly shut down while converting TBranch to NumPy !!
    # fana = data_path / "1_analyse" / filename
    fana = run_path / filename
    print(f"Load file {fana}")

    dout = Path(__file__).parents[1] / 'inb72' / 'out' / f"{str_run}" / f"{str_ana}"
    dout.mkdir(parents = True, exist_ok = True)
    
    # prefix = 'MGv3'
    prefix = "INB72"
    # ltime_var = ['time', 'deltat']
    lclus_var = [f'{prefix}_ClusPos', f'{prefix}_ClusAmpl', f'{prefix}_ClusMaxStripAmpl', f'{prefix}_StripMaxAmpl', f'{prefix}_ClusSize', f'{prefix}_ClusTOT',  f'{prefix}_ClusT', f'{prefix}_ClusMaxSample']
    # dict_time_var = {k: np.asarray(tree[k]) for k in ltime_var}
    branch_list = lclus_var
    n_events = int(sys.argv[3]) if (len(sys.argv) >= 4) else int(1e6)
    print(f"Import n_events = {n_events}")
    
    dict_clus_var = general_import(file=fana, branch_list = branch_list, n_events=n_events)
   
    # print(dict_clus_var[f'{prefix}_ClusMaxSample'][:1,], dict_clus_var[f'{prefix}_ClusAmpl'][:1,])
  
    print(f"TBranches to NumPy -- {time.time() - t0:.1f} s")
    # dict_clus_var = tree.arrays(lclus_var, library='numpy')

    bins = { v : 100 for v in lclus_var}
    bins[f'{prefix}_ClusSize'] = 50
    bins[f'{prefix}_ClusTOT'] = 15
    bins[f'{prefix}_ClusMaxSample'] = 50
    hrange = { v : None for v in lclus_var}
    hrange[f'{prefix}_ClusAmpl'] = [0, 5e3] #ADC
    kwargs_hist = {}#'edgecolor':'black', 'linewidth':0.1}
    ###PLOTS
    # fig = plt.figure(1, figsize= (12,7))
    # gs = GridSpec(1, len(dict_time_var))#, left=0.02, right=0.98, wspace=0.1, hspace=0.5)
    # for i, (k, v) in enumerate(dict_time_var.items()):
    #     ax = fig.add_subplot(gs[0,i], aspect='equal')
    #     ax.hist(v, bins=bins, **kwargs_hist )
    #     ax.set_xlabel(k)
    #     ax.set_ylabel('entries')
    # fig.tight_layout()
    # foutname = "time_var_"+ filename.split('.')[0] + '.png'
    # fout = dout / foutname
    # fig.savefig(fout)
    # print(f'Save figure {fout}')
    # # plt.show()


    var0 = dict_clus_var[lclus_var[0]]
    ndim = len(var0.shape)
    if ndim == 2: 
         for k,v in dict_clus_var.items(): v = v[:, :, np.newaxis]
    nevt_tot, ndet_tot = var0.shape[0], var0.shape[1]
    print(f"clus_var.shape = {var0.shape}")
    print(f"nevt_tot, ndet_tot = {nevt_tot}, {ndet_tot}")

    # nrow, ncol = ndet_tot, len(lclus_var)
    # fig = plt.figure(figsize=(16,11), )#, constrained_layout=True)
    # mosaic = [[ f'{v}_{idet}' for v in lclus_var ] for idet in range(ndet)]
    # # print(mosaic)
    # kwargs_size = dict(wspace=0.0, hspace=0.0, 
    #             top=0.97-0.1/(nrow+1), bottom=0.6/(nrow+1), 
    #             left=0.4/(ncol+1), right=0.95-0.1/(ncol+1) )
    # axs = fig.subplot_mosaic(mosaic, sharey=True, gridspec_kw=kwargs_size)
    # for ivar, (key, var) in enumerate(dict_clus_var.items()):
    #     for idet in range(ndet) : 
    #         v = var[:,idet]    
    #         m = f'{key}_{idet}'
    #         axs[m].hist(v, bins=bins, **kwargs_hist)
    #         if ivar == 0 : 
    #             axs[m].set_ylabel(idet)
    #         if idet == 0 : 
    #             axs[m].set_title(key)
    #         axs[m].set_yticklabels([])
    # plt.gcf().text(x=0.025/(ncol+1), y=0.5, s=f"$\\mu$M",  
    #                fontsize='x-large', 
    #                rotation='vertical')       
    # foutname = "mosaic_clus_var_alldet_"+ filename.split('.')[0] + '.png'
    # fout = dout / foutname
    # fig.savefig(fout)
    # print(f'Save figure {fout}')
    # plt.close()


    
    k = f'{prefix}_ClusPos'
    clus_pos_all_det = dict_clus_var[k]
    # print(f"np.any(~np.isfinite(clus_pos_all_det)) = {np.any(np.any(~np.isfinite(clus_pos_all_det)==True), axis=-1)}")
    # print(f"np.any(np.isnan(clus_pos_all_det)) = {np.any(np.any(np.isnan(clus_pos_all_det)==True), axis=-1)}")
    is_transmu_col = np.int32(np.any(clus_pos_all_det != 0, axis=-1))
    # print("_ClusPos = ", dict_clus_var[f'{prefix}_ClusPos'][0, :, :])
    
    is_transmu = np.count_nonzero(is_transmu_col==True, axis=ndim-2) >= ndet_tot-12 
    is_transmu = np.ones(nevt_tot, dtype=bool) 

    n_transmu = len(is_transmu[is_transmu==True])
    print(f"n_transmitted_mu / nevt_tot = {n_transmu} / {nevt_tot} = {n_transmu/nevt_tot*100:.2f}% ")


 
    tel = INB72
    lay_index = [lay.index for det in tel.detectors for lay in det.strip_layers]
    det_name, det_range = tel.name, range(min(lay_index), max(lay_index)+1)
    lay_pitch = [lay.pitch for det in tel.detectors for lay in det.strip_layers]
    lay_nstrips = [lay.nstrips for det in tel.detectors for lay in det.strip_layers]
    dict_det_clus_var = {k: v[:,det_range] for k, v in dict_clus_var.items()}
    
    v = dict_det_clus_var[lclus_var[0]] 
    ndet = v.shape[1]
    nrow, ncol = ndet, len(lclus_var)
    fig = plt.figure(figsize=(16,11), )
    mosaic = [[ f'{v}_{idet}' for v in lclus_var ] for idet in range(ndet)]
    kwargs_size = dict(wspace=0.0, hspace=0.0, 
                top=0.94-0.1/(nrow+1), bottom=0.4/(nrow+1), 
                left=0.4/(ncol+1), right=0.95-0.1/(ncol+1) )
    axs = fig.subplot_mosaic(mosaic, sharey=True, gridspec_kw=kwargs_size)

    for ivar, (key, var) in enumerate(dict_det_clus_var.items()):

        for idet in range(ndet) : 
            
            v = var[is_transmu,idet]

            mask = ~np.isnan(v) & np.isfinite(v)
            # print(f"{key}_{idet}, np.all(mask==True) = ",np.all(mask==True))
            m = f'{key}_{idet}'
            ax = axs[m]
            kwargs_hist['range'] = hrange[key]
            ax.hist(v[mask], bins=bins[key], **kwargs_hist)
            ax.set_yscale('log')
            if ivar == 0 : 
                if idet%2 == 0 : coord = 'X'
                else: coord = 'Y'
                ax.set_ylabel(f'{coord}{idet}')
            if idet == 0 : ax.set_title(key)
            if idet != ndet-1 : ax.set_xticklabels([])
            # ax.set_yticklabels([])
    
    plt.gcf().text(x=0.025/(ncol+1), y=0.5, s=f"$\\mu$M",  
                fontsize='x-large', 
                rotation='vertical')       
    plt.gcf().text(x=0.5, y=0.985-0.1/(nrow+1), s=f"{det_name}",  
                fontsize='xx-large', 
                rotation='horizontal', 
                color = 'red',
                fontweight='bold')  

    foutname = f"mosaic_clus_var_{det_name}_"+ filename.split('.')[0] + '.png'
    fout = dout / foutname
    fig.savefig(fout, dpi=300)
    print(f'Save figure {fout}')
    plt.close()

    ##Occupation maps
    key_cluspos = f'{prefix}_ClusPos'
    
    fig = plt.figure(figsize=(10,10), )
    npan = ndet//2
    mosaic = np.array([ f'{key_cluspos}_xy_{ipan}' for ipan in range(npan)])
    mosaic = mosaic.reshape((2,2))
    nrow, ncol = mosaic.shape[0], mosaic.shape[1]
    kwargs_size = dict(wspace=0.2, hspace=0.1, 
                top=0.94-0.1/(nrow+1), bottom=0.4/(nrow+1), 
                left=0.3/(ncol+1), right=1-0.3/(ncol+1) )
    axs = fig.subplot_mosaic(mosaic, sharex=True, sharey=True, gridspec_kw=kwargs_size)

    kwargs_map = {}

    for ipan in range(npan):
        v = dict_det_clus_var[key_cluspos][is_transmu, 2*ipan:2*(ipan+1)]
       
        # v = np.swapaxes(v, 1, -1)
        px, py = lay_pitch[2*ipan], lay_pitch[2*ipan+1]
        nx, ny = lay_nstrips[2*ipan], lay_nstrips[2*ipan+1]
        x, y = v[:,0].flatten()*px, v[:,1].flatten()*py #conversion to cluster pos in mm
        nnan = ~np.isnan(x) & ~np.isnan(y) 
        bins_xy, range_xy = [150,150], [[0, px*nx], [0, py*ny]]
        xr, yr = np.linspace(0, nx*px, bins_xy[0]), np.linspace(0, ny*py, bins_xy[1])
        Xr, Yr = np.meshgrid(xr, yr)
        hxy = np.histogram2d(x[nnan], y[nnan], bins=bins_xy, range=range_xy)[0]#range=range_xy
        # hxy = hxy.transpose()
        m = f'{key_cluspos}_xy_{ipan}'
        ax = axs[m]
        vmin, vmax = 1, 1e2#np.nanmax(hxy))
        # print(vmin, vmax)
        # im = ax.imshow(hxy, norm=LogNorm(vmin=vmin, vmax=vmax),**kwargs_map)
        im = ax.pcolormesh(Xr,Yr, hxy, norm=LogNorm(vmin=vmin, vmax=vmax),**kwargs_map)
        ax.invert_yaxis()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(im, cax=cax, extend='max')
        
        sdet = f"X{2*ipan}/Y{2*ipan+1}"
        # anchored_text = AnchoredText(sdet, loc="upper right", frameon=False, prop=dict(fontsize='large', color="red"))
        # ax.add_artist(anchored_text)
        ax.set_title(f"$\\mu$M {sdet}")
        
        # cbar.set_label(label='entries', fontsize='large')
        

    plt.gcf().text(x=0.1, y=0.985-0.1/(nrow+1), s=f"{det_name} : {key_cluspos}",  
                fontsize='xx-large', 
                rotation='horizontal', 
                color = 'red',
                fontweight='bold')  

    foutname = f"mosaic_clus_xy_pos_{det_name}_"+ filename.split('.')[0] + '.png'
    fout = dout / foutname
    fig.savefig(fout, dpi=300)
    print(f'Save figure {fout}')
    plt.close()

