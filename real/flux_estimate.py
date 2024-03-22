#!/usr/bin/python3
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import quad

from forwardsolver import FluxModel

params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'medium',
         'axes.labelpad':10,
         'axes.grid':True,
         'mathtext.fontset': 'stix',
         'font.family': 'STIXGeneral',
         }
plt.rcParams.update(params)


def acceptance_willis(length, x1, x2, y1, y2):
    z = length * 1e-1
    G = (z**2 + 2*(x1 + x2)**2) / (2*np.sqrt(z**2 + (x1 + x2)**2)) * ( (y1+y2)*np.arctan((y1+y2)/np.sqrt(z**2+(x1+x2)**2))  -  (y2-y1)*np.arctan((y2-y1)/np.sqrt(z**2+(x1+x2)**2))  )
    G-= (z**2 + 2*(x2 - x1)**2) / (2*np.sqrt(z**2 + (x2 - x1)**2)) * ( (y1+y2)*np.arctan((y1+y2)/np.sqrt(z**2+(x2-x1)**2))  -  (y2-y1)*np.arctan((y2-y1)/np.sqrt(z**2+(x2-x1)**2))  )
    G+= (z**2 + 2*(y1 + y2)**2) / (2*np.sqrt(z**2 + (y1 + y2)**2)) * ( (x1+x2)*np.arctan((x1+x2)/np.sqrt(z**2+(y1+y2)**2))  -  (x2-x1)*np.arctan((x2-x1)/np.sqrt(z**2+(y1+y2)**2))  )
    G-= (z**2 + 2*(y2 - y1)**2) / (2*np.sqrt(z**2 + (y2 - y1)**2)) * ( (x1+x2)*np.arctan((x1+x2)/np.sqrt(z**2+(y2-y1)**2))  -  (x2-x1)*np.arctan((x2-x1)/np.sqrt(z**2+(y2-y1)**2))  )
    return G



if __name__ == "__main__":

    fm = FluxModel()
    model = 'guan'

    arr_energy = np.logspace(np.log10(1), np.log10(100), 100) #GeV
    arr_theta = np.zeros(1) #rad
    # arr_difflux = np.zeros(shape=(len(arr_energy), len(arr_theta)))
    # for i, ene in enumerate(arr_energy):
    #     for j, theta in enumerate(arr_theta):
    #         arr_difflux[i,j] = fm.ComputeDiffFlux(ene, theta)

    


    # fig, ax = plt.subplots()
    # ax.plot(arr_energy, arr_difflux[:,0])

    # ax.set_xscale('log')
    # ax.set_yscale('log')

    # fig.tight_layout()

    # figname = 'tflux_test'
    # dout = Path(__file__).parent / 'out'
    # dout.mkdir(parents=True, exist_ok=True)
    # fout = dout / f'{figname}.png'
    # plt.savefig(fout)
    # plt.close()


    # theta_c = 0
    theta_c = 0*np.pi/180
    s, L = 170, 180 #mm
    # theta_min, theta_max = 0, np.arctan(s/(2*L))#
    theta_min, theta_max = 0, 90*np.pi/180 # theta_c - np.arctan(s/(2*L)), theta_c + np.arctan(s/(2*L))
    emin, emax = 0.106, 1e3

    # phi_min, phi_max = - np.arctan(s/(2*L)), np.arctan(s/(2*L))
    # phi_min, phi_max = 0, 2*np.pi #
    print(f"theta_min, theta_max = {theta_min*180/np.pi:.1f}, {theta_max*180/np.pi:.1f}")
    # print(f"phi_min, phi_max = {phi_min*180/np.pi:.1f}, {phi_max*180/np.pi:.1f}")
    # arr_theta = np.linspace(theta_min, theta_max, 100)
    
    # arr_emin = np.ones(1)*emin
    # arr_intflux = np.zeros(len(arr_theta))
    # for i, theta in enumerate(arr_theta):
    #     # for j, emin in enumerate(arr_emin):
    #     arr_intflux[i] = quad(fm.ComputeDiffFlux, emin, emax, args=(theta, model))[0] #TransmittedFlux(ene, theta, model=model)


    # fig, ax = plt.subplots()
    # ax.plot(arr_theta*180/np.pi, arr_intflux)
    # # ax.set_xscale('log')
    # ax.set_yscale('log')
    # ax.legend()
    # fig.tight_layout()

    # figname = 'tflux_test'
    # dout = Path(__file__).parent / 'out'
    # dout.mkdir(parents=True, exist_ok=True)
    # fout = dout / f'{figname}.png'
    # plt.savefig(fout)
    # plt.show()

    integrand = lambda theta : quad(fm.ComputeDiffFlux, emin, emax, args=(theta, model))[0] * np.cos(theta)
    r = 20.0
    print(quad(integrand, np.cos(theta_min), np.cos(theta_max) )[0]) 
    S = np.pi*r**2
    print(f"S = {S:.2f} cm2 = {S*1e-4:.2f} m2")
    rate = 2*np.pi * S * abs(quad(integrand, np.cos(theta_min), np.cos(theta_max) )[0]) 
    nmu = 1e5
    time = nmu / rate
    # # mucount = abs(quad(integrand, np.cos(theta_min), np.cos(theta_max) )[0])
    print(f"time = {time/(3600*24):.2f} d")

    # acc = acceptance_willis(L, -s, s, -s, s)










