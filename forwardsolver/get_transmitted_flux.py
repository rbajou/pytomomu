#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import matplotlib.lines as mlines  #use for legend settings
import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d
import seaborn as sns
import time
#package module(s)
from stoppingpower import StoppingPower
from fluxmodel import FluxModel
from tomomuLib_Simulation.OpInv import f_op
from material import str2medium

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


if __name__ ==  "__main__":

    parser=argparse.ArgumentParser(
    description='''Plot muon stopping power (-dE/dx) as a function of kinetic energy for standard rock and water''', epilog="""All is well that ends well.""")
    parser.add_argument('--medium', '-m', default="rock", help='Check available medium in material', type=str2medium)
    parser.add_argument('--rho', '-r', default=None, help='Medium density [g/cm3]', type=float)
    args=parser.parse_args()
    

    t0 = time.time()

    sp = StoppingPower(I=126)

    energy_bins= np.logspace(np.log10(1), np.log10(1e5), 100) #
    theta_bins = np.linspace(0, 90, 90)*np.pi/180 #rad
    op_bins= np.logspace(np.log10(1), np.log10(3000), 600) #mwe = hg/cm^2
  
    fm = FluxModel(altitude=0.)

    mu = 0.1056 #muon mass in GeV

      
    medium = args.medium
    rho = medium.rho if not args.rho else args.rho

    out_dir = Path(__file__).parent / "files"
    out_dir_dedx = out_dir / medium.name / str(rho) / "dedx"
    out_dir_emin = out_dir / medium.name / str(rho) / "emin"
    out_dir_emin.mkdir(parents=True, exist_ok=True)
    out_dir_openflux = out_dir / medium.name / str(rho) / "flux" / "opensky"
    out_dir_openflux.mkdir(parents=True, exist_ok=True)
    out_dir_transflux = out_dir / medium.name / str(rho) / "flux" / "transmitted"
    out_dir_transflux.mkdir(parents=True, exist_ok=True)
    
    ftotout = out_dir_dedx / "dEdx_tot.csv"
    ene_range, dEdx_ion, dEdx_rad = np.loadtxt(ftotout, unpack=True, skiprows=1) 
    print(f"Load {ftotout}")
    dEdx_tot = dEdx_ion + dEdx_rad #MeV/(g/cm^2)
    dEdx_tot_interp = interp1d(ene_range*1e-3, dEdx_tot*1e-3, kind='linear') #GeV/(g/cm^2)
    fdEdx_tot = lambda o, e: dEdx_tot_interp(e)

    opmin, opmax= np.min(op_bins), np.max(op_bins)
    f_emin_out = out_dir_emin / f"emin_{int(opmin)}_{int(opmax)}mwe.txt"
    if not f_emin_out.exists():
        print(op_bins[np.newaxis].shape)
        vec_emin = sp.minimum_energy(func=fdEdx_tot, opacity=op_bins*100, E0=mu) #opacity in g/cm^2
        print(f"emin = {vec_emin} GeV")
        np.savetxt(f_emin_out, vec_emin)
        print(f"Save {f_emin_out}")
    else : 
        vec_emin = np.loadtxt(f_emin_out)
        print(f"Load {f_emin_out}")
    
    kwargs = {} #kwargs={'limit':1000}
    f_flux_trans = out_dir_transflux / f"flux.txt"
    emax = 1e5 #GeV
    # if not f_flux_trans.exists():
    mat_flux_trans = fm.ComputeTransmittedFlux(emin=vec_emin.flatten(), emax=emax, theta=theta_bins, model="guan", **kwargs)
    print(mat_flux_trans, mat_flux_trans.shape)
    np.savetxt(f_flux_trans, mat_flux_trans)
    print(f"Save {f_flux_trans}")
    # else : 
    #     mat_flux_trans = np.loadtxt(f_flux_trans)
    #     print(f"Load {f_flux_trans}")
   
    f_flux_open = out_dir_openflux / f"flux.txt"
    if not f_flux_open.exists():
        vec_flux_open = fm.ComputeOpenSkyFlux( emin=mu, emax=emax, theta=theta_bins, model="guan", **kwargs)
        print(vec_flux_open, vec_flux_open.shape)
        np.savetxt(f_flux_open, vec_flux_open)
        print(f"Save {f_flux_open}")
    else : 
        vec_flux_open = np.loadtxt(f_flux_open)
        print(f"Load {f_flux_open}")


    ###PLOTS
        
    fig, ax = plt.subplots(figsize=(12,10))
    x = op_bins  
    leg_handles = []
    y = vec_emin.flatten()
    ax.plot(x,y, color='blue', linestyle="-")
    ax.tick_params(which="both", bottom=True, top=True, left=True, right=True)
    ax.tick_params(which="both", labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    ax.grid(True, which='both',linestyle='dotted', linewidth="0.25", color='grey')
    ax.axis([np.min(x),np.max(x), 0.1, 1e3])
    ax.set_yscale('log')
    ax.set_xlabel("Opacity $\\varrho$ [mwe]")
    ax.set_ylabel("$E_{min}$ [GeV]")
    ax.legend(handles=leg_handles, loc='upper right')
    fout = f_emin_out.parent / "emin.png"
    plt.savefig(fout)
    print(f"Save {fout}")

    fig, ax = plt.subplots(figsize=(12,10))
    l_theta = [0, 70, 75, 80, 85] #simu_ze[:,0][::30]
    palette = sns.color_palette("Spectral_r", len(l_theta)).as_hex()
    x = op_bins  
    leg_handles = []
    for i, t in enumerate(l_theta) : 
        irow = np.nanargmin(np.abs(theta_bins - t*np.pi/180))
        y= mat_flux_trans[irow,:]
        color = palette[i]
        ax.plot(x,y, color=color, linestyle="-")
        #legend handle
        handle = mlines.Line2D([], [], color=color,  marker=None, linestyle='-',
        markersize=10, label=f"{t:.0f}°")
        leg_handles.append(handle)

    ax.tick_params(which="both", bottom=True, top=True, left=True, right=True)
    ax.tick_params(which="both", labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    ax.grid(True, which='both',linestyle='dotted', linewidth="0.25", color='grey')
    ax.axis([np.min(x),np.max(x), 1e-8, 1e-2])
    ax.set_yscale('log')
    ax.set_xlabel("Opacity $\\varrho$ [mwe]")
    ax.set_ylabel("Transmitted muon flux I [cm$^{-2}$.sr$^{-1}$.s$^{-1}$]")
    ax.legend(handles=leg_handles, loc='upper right')
    fout = f_flux_trans.parent / "flux.png"
    plt.savefig(fout)
    print(f"Save {fout}")

    fig, ax = plt.subplots(figsize=(12,10))
    l_theta = [0, 20, 30, 40, 45, 50] 
    palette = sns.color_palette("Spectral_r", len(l_theta)).as_hex()
    x = op_bins  
    leg_handles = []
    for i, t in enumerate(l_theta) : 
        irow = np.nanargmin(np.abs(theta_bins - t*np.pi/180))
        ynum= mat_flux_trans[irow,:]
        yden = np.tile(vec_flux_open[irow], len(ynum)) #repeat along column axis
        y = ynum/yden
        color = palette[i]
        ax.plot(x,y, color=color, linestyle="-")
        #legend handle
        handle = mlines.Line2D([], [], color=color,  marker=None, linestyle='-',
        markersize=10, label=f"{t:.0f}°")
        leg_handles.append(handle)
    
        y_seb = np.logspace(np.log10(0.01), np.log10(1), 100) #transmission
        x_seb = f_op(theta=np.tan(t*np.pi/180), tx=y_seb) #g.cm⁻2->mwe
        ax.plot(x_seb, y_seb, linestyle='--', color=color) 
 
    ax.set_xscale('log')
    ax.tick_params(which="both", bottom=True, top=True, left=True, right=True)
    ax.tick_params(which="both", labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    ax.grid(True, which='both',linestyle='dotted', linewidth="0.25", color='grey')
    ax.axis([np.min(x),np.max(x), -1, 1])
    ax.set_xlabel("Opacity $\\varrho$ [mwe]")
    ax.set_ylabel("Transmission")
    ax.legend(handles=leg_handles, loc='upper right')
    fout = f_flux_trans.parent / "transmission.png"
    plt.savefig(fout)
    print(f"Save {fout}")
    print(f"End -- {(time.time() - t0)/60:.2f}  min")  