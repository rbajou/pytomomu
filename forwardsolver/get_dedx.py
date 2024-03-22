
#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines  #use for legend settings
import matplotlib.patches as mpatch  #use for legend settings
import os
import time 
import argparse
from pathlib import Path

#package module(s)
from forwardsolver.stoppingpower import StoppingPower
from material import DICT_MEDIUM, str2medium

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

if __name__ == "__main__":
    home_dir = Path(os.environ["HOME"])
    work_dir = Path(__file__).parent
    
    t_start = time.time()
    print("Start: ", time.strftime("%H:%M:%S", time.localtime()))#start time

    parser=argparse.ArgumentParser(
    description='''Plot muon stopping power (-dE/dx) as a function of kinetic energy for standard rock and water''', epilog="""All is well that ends well.""")
    parser.add_argument('--out_dir', '-o', default=str(work_dir / "files"), help='Path to processing output', type=str)
    parser.add_argument('--medium', '-m', default="rock", help='Check available medium in material', type=str2medium)
    parser.add_argument('--rho', '-r', default=None, help='Medium density [g/cm3]', type=float)

    args=parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    T = np.logspace(np.log10(1e3),np.log10(1e9),100) ####MeV
    mu = 105.65839 #MeV 
    E_MeV = T+mu 
    E_GeV = E_MeV*1e-3
    
    

    medium = args.medium
    print(args.rho, medium.rho)
  
    try :   
        f_dedx_table = work_dir/ "files" / "dedx" /f"dEdx_{medium.name}_pdg.csv"
        ene_pdg, dEdx_ion_pdg = np.loadtxt(f_dedx_table, skiprows=2, usecols=(0,2), unpack=True)
    except : 
        ene_pdg = np.logspace(np.log10(1), np.log(1e8), 100) #MeV
       
    I = medium.I
    A, Z, Z_A = medium.A, medium.Z, medium.Z_A 
    sp = StoppingPower(I=I)
    
    rho = medium.rho if not args.rho else args.rho

    out_dir_tmp = out_dir / medium.name / str(rho) / "dedx"
    (out_dir_tmp).mkdir(parents=True, exist_ok=True)   
    ftotout = out_dir_tmp / "dEdx_tot.csv"
    faout = out_dir_tmp / "dEdx_ion.csv"
    fbout = out_dir_tmp / "b_rad.csv"
    cut_rad = (ene_pdg  > 20*mu)
    
    if ftotout.exists():
        _, dEdx_ion, dEdx_rad = np.loadtxt(ftotout, unpack=True, skiprows=1) 
        dEdx_tot = dEdx_rad+dEdx_ion
        
    else :
        dEdx_ion = sp.bethe_bloch(T=ene_pdg, Z_A=Z_A, rho=rho, brem_corr=True, rho0=medium.rho, dens_effect=True) #in MeV/(g/cm^2)
        arr_ion = np.vstack((ene_pdg, dEdx_ion)).T
        np.savetxt(faout, arr_ion, fmt="%.5e", delimiter="      ", header="T(MeV)      dEdx")
        print(f"save {faout}")
        dEdx_tot = np.copy(dEdx_ion)
        gamma = ene_pdg/mu 
        beta = np.sqrt(1-1/gamma**2)
        ene_pdg_GeV = ene_pdg*1e-3
        sub_E_GeV =  ene_pdg_GeV[cut_rad][::5]
        b_pair = [sp.b_pair(e, A,Z) for e in sub_E_GeV] #in cm^2/g
        b_pair_interp = np.interp(ene_pdg_GeV[cut_rad], sub_E_GeV, b_pair)  #in cm^2/g
        b_brems = [sp.b_brems(e, A,Z) for e in ene_pdg_GeV[cut_rad]] #in cm^2/g
        b_phnuc = [sp.b_nuc(e, A,Z) for e in ene_pdg_GeV[cut_rad]] #in cm^2/g
        b_tot =   np.sum(np.vstack((b_pair_interp, b_brems, b_phnuc)).T,axis=1)
        dEdx_rad = np.zeros(len(dEdx_tot))
        arr_rad = np.vstack((ene_pdg[cut_rad], b_pair_interp, b_brems, b_phnuc, b_tot)).T
        np.savetxt(fbout, arr_rad,  fmt="%.5e", 
                   delimiter="      ", 
                   header="E(MeV)      b_pair      b_brem     b_nucl     b_total(cm^2/g)")
        print(f"save {fbout}")
        
        dEdx_rad[cut_rad] = b_tot*ene_pdg[cut_rad]
        dEdx_tot[cut_rad] += b_tot*ene_pdg[cut_rad]
        arr_tot = np.vstack((ene_pdg, dEdx_ion, dEdx_rad)).T
        np.savetxt(ftotout, arr_tot, fmt="%.5e",
                    delimiter="      ", 
                   header="E(MeV)      dEdx_ion(MeV.cm^2/g)      dEdx_rad(MeV.cm^2/g)")
        print(f"save {ftotout}")

print(f'End -- {(time.time()-t_start)/60:.2f} min')
