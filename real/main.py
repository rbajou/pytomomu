#!/usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import subprocess
import numpy as np
from pathlib import Path
import time
import uproot
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
import argparse
#package module(s)
# from config import  DETECTOR, DATA_PATH, RUNS
from analyse import Analyse, Final
from config.survey import CURRENT_SURVEY
from detector import Bench, Telescope, VetoPlane
from eventrate import EventRate
from hitmap import HitMap
from utils import gaus

t0 = time.time()
main_path = Path.home() / "Projects" / "tomomu"
det = CURRENT_SURVEY.detector 
if isinstance(det, Telescope): tel = det
elif isinstance(det, Bench) : tel =  det.telescope[0]
else : raise TypeError("Unknown detector type")

data_path = Path(CURRENT_SURVEY.data_path) #main_path / "Data" / "inb72" 
print(f"Data path : {data_path}")

run_prefix = CURRENT_SURVEY.run_prefix 
# prefix = "INB72_2024"

parser=argparse.ArgumentParser(
description='''.''', epilog='''All is well that ends well.''')
parser.add_argument('--run_num', '-r', default=None, help="Run number (int)",  type=int)
parser.add_argument('--nevt_max', '-n', default=10000,  type=int)
parser.add_argument('--ana', '-a', default=1,  type=int, help="'1' or '2'")

args=parser.parse_args()

lrun  = [args.run_num] if args.run_num else CURRENT_SURVEY.run_num 
n_events = args.nevt_max if args.nevt_max else int(1e4)

clusvar_prefix = CURRENT_SURVEY.clusvar_prefix
iana=  args.ana
str_ana = f"analyse{iana}"



for irun in lrun:
    str_run = f"run{irun}" 
    run_path = data_path / str_run
    label = f"{run_prefix}_run{irun}"
    # str_run = f"runtest_{irun}"  
    filename = f"{run_prefix}_{str_run}_{str_ana}.root"
    print(f"\nProcess {filename.split('.')[0]}\t. . .")
    basename = filename.split('.')[0]
    #dout = Path(__file__).parent / 'out' / str_run / str_ana
    dout = run_path / "out"
    dout.mkdir(parents = True, exist_ok = True)
    # fana = data_path / "1_analyse" / filename
    fana = run_path  / filename
    final_path =  run_path / "final"  
    final_path.mkdir(parents=True, exist_ok=True)
    ffinal = final_path / f"{run_prefix}_{str_run}_final.root"
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
    ana.open()
    print(f"TBranches to NumPy -- {time.time() - t0:.1f} s\n")

    ana.get_content(n_events = n_events) 
    ana.fill_event_collection()
    ec = ana.event_collection
    ec.detector_ensemble = det

    '''
    print("\nEvent Rate")
    fout = dout / "eventrate.png"
    fig, ax = plt.subplots(figsize=(12,8))
    evttime = ana.dict_clus_var["evttime"]
    with open(str(run_path / "timestamp0.txt"), 'r') as f : t0_run = int(f.readline())
    er = EventRate(evttime, t0=t0_run)
    print(f"Run duration : {er.run_duration:.0f} s")
    label_er = f"all"
    er(ax, width=1, label=label_er) #width = size time bin width in hour 
    ax.legend(loc='best')
    plt.gcf().text(x=0.1, y=0.95, s=f"{label} : event rate", fontsize='xx-large', rotation='horizontal', color = 'red', fontweight='bold', ha='left')  
    plt.savefig(str(fout), dpi=300)
    print(f"Save figure {fout}\n")
    # cut = np.any(ana.dict_clus_var[f'{clusvar_prefix}_ClusSize'][:,] <= 1, axis=1)
'''
    mask = np.full(len(ec.events), True) 
    event_no = np.array([ev.id for i, (_, ev) in enumerate(ec.events.items())])
    event_time = np.array([ev.time for i, (_, ev) in enumerate(ec.events.items())])
    clus_size = np.array([ev.cluster_collection.size for i, (_, ev) in enumerate(ec.events.items())])# if mask[i] == True ])
    mask_clus_size = (clus_size > 1) #np.full(clus_size.shape, True)  #

    hitmap = HitMap(event_collection=ec)
    foutname = f"mosaic_clus_xy_pos_{ec.detector_ensemble.name}_"+ basename + '.png'
    fout = dout / foutname
    clus_pos = np.array([ev.cluster_collection.position for i, (_, ev) in enumerate(ec.events.items())])
 
    ##clus_pos shape : (nev, nlay, nclus)
    kwargs = {"label":f"{run_prefix}_run{irun}"}
    hitmap.plot(clus_pos=clus_pos, mask=mask_clus_size, file=fout, **kwargs)

    # llayer_id_tel = [lay.index for det in det.detectors for lay in det.strip_layers]
    # range_layer = range( min(llayer_id_tel), max(llayer_id_tel)+1 )
    # ldet_id_tel = [det.index for det in det.detectors]
 
    # if not ffinal.exists():
    ec.get_xyz(mask=mask_clus_size)
    # print(ec.xyz.shape, ec.xyz[10,:])
    ##col0 : evt_ix, col1: det_ix, col2: clus_ix, col3: coord_ix
    XYZ = np.expand_dims(ec.xyz[:,-4:,0,:], axis=2) #expand dim to keep four axis and keep only detector panels
    # '''
    ec.get_tracks(XYZ)

    print(f"\nGet tracks ({len(ec.tracks)}) -- {time.time()-t0:.1f} s")
    dict_final = {}
    lid, lev = list(ec.events.keys()), list(ec.events.values())#[:10]
    dict_final["evn"] = list(lid)
    dict_final["evttime"] = np.array([ev.time for ev in lev])
    zdet = np.array([d.position[-1] for d in tel.detectors])
    ix_bottom = np.argmin(zdet)
    z0 = tel.detectors[ix_bottom].position[-1]
    ltrack = [ev.track for ev in lev]
    a, b = np.array([trk.a for trk in ltrack]), np.array([trk.b for trk in ltrack])
    dict_final['tthetax'], dict_final['tthetay'],  =  np.array([trk.tthetax for trk in ltrack]), np.array([trk.tthetay for trk in ltrack])
    (ax, ay), (bx, by) = a.T, b.T
    x0, y0 = ax * z0 + bx, ay * z0 + by 
    dict_final["x0"], dict_final["y0"]  = x0, y0
    r = np.array([trk.r for trk in ltrack])
    dict_final["resx"], dict_final["resy"]  = r[:, 0], r[:,1]
    root_out = uproot.recreate(str(ffinal))
    root_out['T'] = dict_final
    root_out.close()
    print(f"Save file {ffinal}\n")

    # if ffinal.exists():
    final = Final(ffinal)
    final.get_content(n_events=n_events)
    tthetax, tthetay = final.content['tthetax'], final.content['tthetay']
    mask_track = (tthetax != -999) & (tthetay!= -999) 
    # print(f"mask_track = {np.all(mask_track == True)}, {np.all(mask_track == False)}")

    x, y = tthetax[mask_track], tthetay[mask_track]
    fig, ax = plt.subplots(figsize=(10,10))
    xmin, xmax, ymin, ymax = np.nanmin(x), np.nanmax(x), np.nanmin(y), np.nanmax(y)
    bins, range = 150, [[-1, 1], [-1, 1]]#[[xmin, xmax], [ymin, ymax]]
    h2d, binsx, binsy = np.histogram2d(x, y, bins=bins)#, range=range)
    xc, yc = (binsx[:-1]+ binsx[1:])/2, (binsy[:-1]+ binsy[1:])/2
    XC, YC = np.meshgrid(xc, yc)
    im = ax.pcolormesh(XC, YC, h2d.T, norm=LogNorm(vmin=1, vmax=np.max(h2d)))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, extend='max')
    ax.set_xlabel('tan($\\phi$)')
    ax.set_ylabel('tan($\\theta_y$)')
    ax.set_title(f'{label}', fontsize='xx-large', color = 'red', fontweight='bold')
    ax.invert_yaxis()
    foutname = f"tthetax_tthetay_{ec.detector_ensemble.name}_"+ basename + '.png'
    fout = dout / foutname
    fig.savefig(fout, dpi=300)
    print(f'Save figure {fout}')
    plt.close()

    x0, y0 = final.content['x0'][mask_track], final.content['y0'][mask_track]
    fig, ax = plt.subplots(figsize=(10,10))
    bins, range = 150, [[-170, 170], [-170, 170]]
    # print(f"x0, y0 = [{np.nanmin(x0)}, {np.nanmax(x0)}], [{np.nanmin(y0)}, {np.nanmax(y0)}] mm ")
    hxy0, binsx, binsy = np.histogram2d(x0, y0, bins=bins, range=range)
    xc, yc = (binsx[:-1]+ binsx[1:])/2, (binsy[:-1]+ binsy[1:])/2
    # print(xc.shape, yc.shape, hxy0.shape)
    XC, YC = np.meshgrid(xc, yc)
    # print(f"XC, YC = [{np.nanmin(XC)}, {np.nanmax(XC)}], [{np.nanmin(YC)}, {np.nanmax(YC)}] mm ")
    im = ax.pcolormesh(XC, YC, hxy0.T, norm=LogNorm(vmin=1, vmax=np.max(hxy0)))
    ax.invert_yaxis()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, extend='max')
    ax.set_xlabel('X0 [mm]')
    ax.set_ylabel('Y0 [mm]')
    ax.set_title(f'{label}', fontsize='xx-large', color = 'red', fontweight='bold')
    ax.invert_yaxis()
    foutname = f"xy0_{ec.detector_ensemble.name}_"+ basename + '.png'
    fout = dout / foutname
    fig.savefig(fout, dpi=300)
    print(f'Save figure {fout}')
    plt.close()


    #RESIDUALS (track-cluster) in each panel
    
    npan = len(tel.detectors)

    for ipan in np.arange(0, npan):
        
        resx, resy = final.content['resx'][mask_track][:,ipan], final.content['resy'][mask_track][:,ipan]
        # fig, ax = plt.subplots(figsize=(10,10))
        fig = plt.figure(figsize=(12, 10))
        gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                            left=0.1, right=0.9, bottom=0.1, top=0.9,
                            wspace=0.05, hspace=0.05)
        ax = fig.add_subplot(gs[1, 0])
        ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
        ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
        bins, range = 150, [[-5, 5], [-5, 5]]

        hxy, binsx, binsy = np.histogram2d(resx, resy, bins=bins, range=range)
        xc, yc = (binsx[:-1]+ binsx[1:])/2, (binsy[:-1]+ binsy[1:])/2
        wx, wy = (binsx[1:]-binsx[:-1]), (binsy[1:]-binsy[:-1])
        XC, YC = np.meshgrid(xc, yc)
        im = ax.pcolormesh(XC, YC, hxy.T, norm=LogNorm(vmin=1, vmax=np.max(hxy0)))
        ax.invert_yaxis()
        ax.set_xlabel('res x [mm]')
        ax.set_ylabel('res y [mm]')

        hx, _, _ = ax_histx.hist(resx, bins=bins, range=range[0], orientation='vertical', histtype='step', facecolor="None", edgecolor="blue", lw=1)
        fran = (range[0][0]/2 < xc) & (xc < range[0][1]/2)
        popt, pcov = curve_fit(gaus, xc[fran], hx[fran], bounds=(0, [1e6, 5, 25]))
        errx = np.sqrt(np.diagonal(pcov))
        xnew = np.linspace(range[0][0]/2, range[0][1]/2, bins)
        # ax_histx.plot(xnew, gaus(xnew, *popt), color='red', label=f'fit:\nmean={popt[1]:.2f}$\\pm${errx[1]:.2f} mm ,\n sigma={popt[2]:.2f}$\\pm${errx[2]:.2f} mm')
        # ax_histx.legend(loc='upper right')

        hy, biny, _ = ax_histy.hist(resy, bins=bins, range=range[1], orientation='horizontal', histtype='step', facecolor="None", edgecolor="blue", lw=1)
        yc = (binsy[:-1]+ binsy[1:])/2
        fran = (range[1][0]/2 < yc) & (yc < range[1][1]/2)
        popty, pcovy = curve_fit(gaus, yc[fran], hy[fran], bounds=(0, [1e6, 5, 25]))
        erry = np.sqrt(np.diagonal(pcovy))
        ynew = np.linspace(range[1][0]/2, range[1][1]/2,  bins)
        # ax_histy.plot( gaus(yc, *popty), ynew, color='red', label=f'fit:\nmean={popty[1]:.2f}$\\pm${erry[1]:.2f} mm,\n sigma={popty[2]:.2f}$\\pm${erry[2]:.2f} mm')
        # ax_histy.legend(loc='upper left')

        ax_histx.tick_params(axis="x", labelbottom=False)
        ax_histy.tick_params(axis="y", labelleft=False)
        plt.gcf().text(x=0.1, y=0.95, s=f"{label} : residuals in X and Y in panel {ipan}",  
                            fontsize='xx-large', 
                            rotation='horizontal', 
                            color = 'red',
                            fontweight='bold', ha='left')  

        ax.invert_yaxis()
        foutname = f"res2d_{ec.detector_ensemble.name}_"+ basename + f'_panel{ipan}.png'
        fout = dout / foutname
        fig.savefig(fout, dpi=300)
        print(f'Save figure {fout}')
        plt.close()
# '''