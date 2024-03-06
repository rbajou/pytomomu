#!/usr/bin/python3
# -*- coding: utf-8 -*-

from typing import List, Union
import numpy as np
from pathlib import Path
import uproot
import time
import re
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
from matplotlib.offsetbox import AnchoredText
import pandas as pd
import glob
from datetime import datetime
import sys
from scipy.integrate import quad
import json
from scipy import stats 

#package module(s)
from hitmap import HitMap
from rootfile import SimuOutput
from detector import INB72
from analyse import Final
from forwardsolver import FluxModel

t0 = time.time()

prefix = "G4TomoDet"

# sim = "INB72" 
sim_name = "INB72_onsite"
main_path = Path(__file__).parents[2] / "sps" / "tomomu" / "rbajou" 
# main_path = Path(__file__).parents[2]  
# main_path = Path.home()


tel = INB72

g4_dirname = f"{prefix}_{tel.name}" 
g4_simpath = main_path / g4_dirname
# g4_outpath = main_path / f"{prefix}_Data" / sim_name
g4_outpath = main_path / "data" / sim_name
# g4_outpath = main_path / "Data" / "G4" / sim_name



# gen, phi, theta, h, irun= "Guan", 39.6, 46.8, 3350, 'all' 
fmd = "metadata_simu.json"
try :
	assert (Path(__file__).parents[1]/ fmd ).exists()
except AssertionError:
	exit(f"{Path(__file__).parent} does not exist")
try :
	gen, phi, theta, h, irun = [None]*5
	with open(fmd) as f: 
		dict_md = json.load(f)
		gen, phi, theta, h, irun = dict_md['gen'], dict_md['phi'], dict_md['theta'], dict_md['h'], dict_md['run']
except : 
	exit(f"Check keys in {fmd}")

# gen, phi, theta, h, irun= "Guan", 0, 46.8, 3350, 'all' 
# gen, phi, theta, irun= "Guan", 39.6, 45, 'all' #"Guan"
# gen, theta, irun= "Flat", 0, 'all'
datestr = sys.argv[1] 
# datestr = "240207*"
print(datestr)
basename = f"{prefix}_{gen}_phi{int(np.around(phi,0))}_theta{int(np.around(theta,0))}_h{h}_run{irun}"
# basename = f"{prefix}_{gen}_phi{iFnt(np.around(phi,0))}_theta{int(np.around(theta,0))}_run{irun}"
# basename = f"{prefix}_{gen}_{theta}deg_run{irun}"
# basename ="G4TomoDet_Flat_0deg_Test10"

# filename = f"{basename}_{datestr}.root"
filename = f"{basename}*{datestr}*root"

print("Input ", str(g4_outpath / filename))
file_path =  glob.glob(str(g4_outpath / filename))[0]
filename = file_path.split("/")[-1]
dout = main_path / "out" / sim_name

final_path =  dout / "final"  
final_path.mkdir(parents=True, exist_ok=True)

is_test = False

ffinal = final_path / f"{filename.split('.')[0]}_final.root"
print(f"Read file {file_path}")
entry_stop = 1e7
if is_test : 
	ffinal = final_path / f"{filename.split('.')[0]}_final_test.root"
	entry_stop = 1e4

kwargs= {'entry_start':0, 'entry_stop':entry_stop}
mu = SimuOutput(file=file_path, **kwargs)

nevtot = len(mu.content['Event'])
if irun == "all" : mu.content['Event'] = np.arange(1, nevtot+1)
print(f"\nnev : {nevtot}\n")

X, Y, Z = mu.content['HitPosX'], mu.content['HitPosY'], mu.content['HitPosZ']

rotx, rotz = -theta*np.pi/180, phi*np.pi/180
mrotx = np.array([[1,0,0], [0,np.cos(rotx),-np.sin(rotx)], [0,np.sin(rotx),np.cos(rotx)]]) 
mrotz = np.array([[np.cos(rotz),-np.sin(rotz), 0], [np.sin(rotz),np.cos(rotz), 0], [0,0,1]]) 

det = mu.content['HitDet']
xp, yp, zp = mu.content["GenPos"].T
xd, yd, zd = mu.content["GenDir"].T

print(f"GenPos : (min, max) = ({np.nanmin(xp):.2f}, {np.nanmin(yp):.2f}, {np.nanmin(zp):.2f}), ({np.nanmax(xp):.2f}, {np.nanmax(yp):.2f}, {np.nanmax(zp):.2f})")
print(f"GenDir : (min, max) = ({np.nanmin(xd):.2f}, {np.nanmin(yd):.2f}, {np.nanmin(zd):.2f}), ({np.nanmax(xd):.2f}, {np.nanmax(yd):.2f}, {np.nanmax(zd):.2f})")


md = SimuOutput(file=file_path, treename="MetaData")
disk_radius =  md.content["Radius"][0]
phi_gen_min, phi_gen_max =  md.content["PhiMin"][0], md.content["PhiMax"][0] #not correctly set in G4TomoDet
theta_gen_min, theta_gen_max =  md.content["ThetaMin"][0], md.content["ThetaMax"][0]
zgenplane =  md.content["GenPlaneZ"][0]

#position on emission disk

print("\n MetaData Generation: ")
print(f"disk_radius = {disk_radius} mm")
print(f"phi_min, phi_max = {phi_gen_min*180/np.pi:.1f}, {phi_gen_max*180/np.pi:.1f} °")
print(f"theta_min, theta_max = {theta_gen_min*180/np.pi:.1f}, {theta_gen_max*180/np.pi:.1f} °")
print(f"zgenplane = {zgenplane:.0f} mm")

fm = FluxModel()
model = 'guan'
pmin, pmax = 0.1056, 1000 #GeV
integrand = lambda theta : quad(fm.ComputeDiffFlux, pmin, pmax, args=(theta, model))[0] * np.cos(theta)
r = disk_radius*1e-1 #mm->cm
phi_prim, theta_prim = np.arctan(yd/xd),  np.pi - np.arccos(zd)
m = np.abs(np.round(phi_prim*180/np.pi, 1)) != 90.0
# dphi = np.max(phi_prim[m]) - np.min(phi_prim[m]) #2*np.pi if vertical telescope 
dphi = np.pi/2
dtheta = theta_gen_max - theta_gen_min
print(f"dphi, dtheta = {dphi*180/np.pi}, {dtheta*180/np.pi}°")
s_gen = np.pi * r**2 #disk surface
rate =  dphi * s_gen * abs(quad(integrand, np.cos(theta_gen_min), np.cos(theta_gen_max) )[0]) 
texp = nevtot / rate
print(f"Exposure time = {texp/(3600*24):.2f} d")

label = f"{prefix}_{tel.name}_{gen}_phi{phi}deg_theta{theta}deg_h{h}mm \n(exposure {texp/(3600*24):.2f} days)"
kwargs = {'label':label}

print(f"phi, theta in [{np.min(phi_prim[m])*180/np.pi:.2f}°, {np.max(phi_prim[m])*180/np.pi:.2f}]°, [{np.min(theta_prim)*180/np.pi:.2f}, {np.max(theta_prim)*180/np.pi:.2f}]°")
print(f"phi_mean, phi_std = {np.mean(phi_prim)*180/np.pi:.2f}, {np.std(phi_prim)*180/np.pi:.2f}°")
print(f"theta_mean, theta_std = {np.mean(theta_prim)*180/np.pi:.2f}, {np.std(theta_prim)*180/np.pi:.2f}°")


theta_off = 0. #-np.pi/4
# m = (90*np.pi/180 > theta_prim) & (theta_prim > 0)
#histo incident angular direction
# phi_prim [(phi_prim < 0)] = phi_prim[(phi_prim < 0)] + np.pi
# fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
# bins=150
# az, ze = phi_prim, theta_prim
# zbins, abins = np.linspace(0,np.max(ze), bins), np.linspace(0, 2*np.pi, bins)
# h2d, _, _ = np.histogram2d(az, ze , bins=(abins, zbins))
# A, Z = np.meshgrid(abins, zbins)
# im = ax.pcolormesh(A, Z, h2d.T, cmap="viridis")
# ax.set_theta_zero_location("S")
x, y =  xp, yp #np.tan(phi_prim), np.tan(theta_prim+theta_off)
fig, ax = plt.subplots(figsize=(10,10))
range, bins = [[-800, 800], [-800, 800]], 150
h2d, binsx, binsy = np.histogram2d(x, y, bins=bins)#, range=range)
xc, yc = (binsx[:-1]+ binsx[1:])/2, (binsy[:-1]+ binsy[1:])/2
XC, YC = np.meshgrid(xc, yc)
im = ax.pcolormesh(XC, YC, h2d.T, norm=LogNorm(vmin=1, vmax=np.max(h2d)))
ax.set_xlabel('X [mm]')
ax.set_ylabel('Y [mm]')
ax.set_title(f'{label} : XY Gen', fontsize='xx-large', color = 'red', fontweight='bold')
dy_tel_storwell = 2165 #mm
spacing_well = 900 #mm
radius_well_top = 305 #mm
radius_well_bot = 200 #mm
kwargs = {"fill": False, "color": "red", "linestyle":"dashed"}
fontdict = {"size":"x-large", "color": "red", "weight":"bold", "ha":"center"}
telwell1 = plt.Circle((0, 0), radius_well_top, **kwargs)
ax.text(0., 0.,s="63", fontdict=fontdict)
storwell1 = plt.Circle((0., dy_tel_storwell), radius_well_top, **kwargs)
ax.text(0., dy_tel_storwell,s="53", fontdict=fontdict)
storwell2 = plt.Circle((-spacing_well, dy_tel_storwell), radius_well_top, **kwargs)
ax.text(-spacing_well, dy_tel_storwell,s="52", fontdict=fontdict)
storwell3 = plt.Circle((-2*spacing_well, dy_tel_storwell), radius_well_top, **kwargs)
ax.text(-2*spacing_well, dy_tel_storwell,s="51", fontdict=fontdict)
ax.add_patch(telwell1)
ax.add_patch(storwell1)
ax.add_patch(storwell2)
ax.add_patch(storwell3)
ax.set_ylim(0, np.max(YC))
anchored_text = AnchoredText(f"entries: {nevtot:.2e}", loc="upper right", frameon=True, prop=dict(fontsize="medium", color="red"))
ax.add_artist(anchored_text)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax, extend='max')
foutname = f"xy_gen_{basename}.png"
fout = dout / foutname
fig.savefig(fout, dpi=300)
print(f"Save figure {fout}")
Nmu0 = h2d


mu.fill_event_collection(ndet=len(tel.detectors), nh=8, mrotx=mrotx, mrotz=mrotz)
ec = mu.event_collection
ec.detector_ensemble = tel

hitmap = HitMap(event_collection=ec)
foutname = f"mosaic_hit_xy_pos_{basename}.png"
fout = dout / foutname
hit_pos = np.array([ev.xyz for i, (_, ev) in enumerate(ec.events.items())])
# print(np.min(hit_pos[:,0]), np.max(hit_pos[:,0]), np.min(hit_pos[:,1]), np.max(hit_pos[:,1]))
##hit_pos shape : (nev, nlay, nclus)
kwargs = {"label":basename}
nhit =  np.array([ev.hit_collection.nhit for i, (_, ev) in enumerate(ec.events.items())])# if mask[i] == True ])
# mask = np.full(len(ec.events), True) 
# mask_hit_mult = (nhit > 1)
hitmap.plot(hit_pos=hit_pos, mask=None, file=fout, **kwargs)


if not ffinal.exists():
	ec.xyz = np.array([ev.xyz for _, ev in ec.events.items() ])
	XYZ = np.expand_dims(ec.xyz[:,:,0,:], axis=2) #expand dim to keep four axis: XYZ.shape = (nev, ndet, nhit, 3)
	ec.get_tracks(XYZ)
	print(f"Get tracks -- {time.time()-t0:.1f} s")
	print(f"ntrack = {len(ec.tracks)}")
	dict_final = {}
	lid, lev = list(ec.events.keys()), list(ec.events.values())#[:10]
	dict_final["evn"] = list(lid)
	dict_final["evttime"] = np.zeros(len(lev))
	zdet = np.array([d.position[-1] for d in tel.detectors])
	ix_bottom = np.argmin(zdet)
	z0 = tel.detectors[ix_bottom].position[-1]
	ltrack = [ev.track for ev in lev]
	a, b = np.array([trk.a for trk in ltrack]), np.array([trk.b for trk in ltrack])
	dict_final['tthetay'], dict_final['tthetax'] = np.array([trk.tthetay for trk in ltrack]), np.array([trk.tthetax for trk in ltrack])
	(ax, ay), (bx, by) = a.T, b.T
	x0, y0 = ax * z0 + bx, ay * z0 + by 
	dict_final["x0"], dict_final["y0"]  = x0, y0
	r = np.array([trk.r for trk in ltrack])
	dict_final["resx"], dict_final["resy"]  = r[:, 0], r[:,1]
	root_out = uproot.recreate(str(ffinal))
	root_out['T'] = dict_final
	root_out.close()
	print(f"Save file {ffinal}")

'''
dout = main_path / "out" / sim 
(dout/"evd").mkdir(parents=True, exist_ok=True)
for i in np.arange(1, 11):
	trk = ec.tracks[i]
	id_evt = trk.id
	ev = ec.events[id_evt]
	xyz0 = ev.xyz
	m = np.all(xyz0 != -999, axis=2)
	edep = ev.hit_collection.edep
	xyz_t = np.array([ev.track.intersection(z) for z in np.linspace(0, 180, 50)])
	fig = plt.figure()
	ax = Axes3D(fig)
	ec.detector_ensemble.plot3D(fig=fig, ax=ax, rotx=theta*np.pi/180)
	ax.scatter(
			xyz0[m, 0],
			xyz0[m, 1],
			xyz0[m, 2],
			s= edep[m]*1e5,
			c= 'limegreen',
			marker='o',
			edgecolor='green',
			label= 'hit' 
	)
	ax.set_xlim(-175, 175)
	ax.set_ylim(-100, 200)
	ax.set_zlim(-50., 175)
	ax.set_title(f"Event {id_evt}",fontsize='xx-large', color = 'red', fontweight='bold')
	# ax.legend(handles=handles, fontsize='large',  loc=(0.5,0.85))#'upper right')#, pad=-10)
	# ax.view_init(elev=10., azim=-60)
	ax.view_init(elev=0., azim=0.)
	ax.plot(xyz_t[:,0],xyz_t[:,1], xyz_t[:,-1],
			c="blue", linewidth=0.75)

	fout = dout/"evd"/f"evt{id_evt}.png"
	fig.savefig(fout, dpi=300)
	print(f"Save evd#{id_evt} in {fout}")
'''

final = Final(ffinal)
final.get_content(n_events=int(1e7))
tthetax, tthetay = final.content['tthetax'], final.content['tthetay']
mask_track = (tthetax != -999.) & (tthetay!= -999.) 
ntrack = len(mask_track) 
x0, y0 = final.content['x0'][mask_track], final.content['y0'][mask_track]

x, y = x0, y0
fig, ax = plt.subplots(figsize=(10,10))
xmin, xmax, ymin, ymax = np.nanmin(x), np.nanmax(x), np.nanmin(y), np.nanmax(y)
bins, range = 150, [[-170, 170], [-170, 170]]#[[xmin, xmax], [ymin, ymax]]
h2d, binsx, binsy = np.histogram2d(x0, y0, bins=bins, range=range)
xc, yc = (binsx[:-1]+ binsx[1:])/2, (binsy[:-1]+ binsy[1:])/2
XC, YC = np.meshgrid(xc, yc)
im = ax.pcolormesh(XC, YC, h2d.T, norm=LogNorm(vmin=1, vmax=np.max(h2d)))
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax, extend='max')
anchored_text = AnchoredText(f"entries: {ntrack:.2e}", loc="upper right", frameon=True, prop=dict(fontsize="medium", color="red"))
ax.add_artist(anchored_text)
ax.set_title(f'{label}', fontsize='xx-large', color = 'red', fontweight='bold')
ax.set_xlabel('X0 [mm]')
ax.set_ylabel('Y0 [mm]')
# ax.invert_yaxis()
foutname = f"xy0_{ec.detector_ensemble.name}_{basename}.png"
fout = dout / foutname
fig.savefig(fout, dpi=300)
print(f'Save figure {fout}')
plt.close()

x, y = tthetax[mask_track], tthetay[mask_track]
fig, ax = plt.subplots(figsize=(10,10))
xmin, xmax, ymin, ymax = np.nanmin(x), np.nanmax(x), np.nanmin(y), np.nanmax(y)
bins, range = 150, [[-1, 1], [-1, 1]]#[[xmin, xmax], [ymin, ymax]]
h2d, binsx, binsy = np.histogram2d(x, y, bins=bins, range=range)
xc, yc = (binsx[:-1]+ binsx[1:])/2, (binsy[:-1]+ binsy[1:])/2
XC, YC = np.meshgrid(xc, yc)
im = ax.pcolormesh(XC, YC, h2d.T, norm=LogNorm(vmin=1, vmax=np.max(h2d)))
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax, extend='max')
anchored_text = AnchoredText(f"entries: {ntrack:.2e}", loc="upper right", frameon=True, prop=dict(fontsize="medium", color="red"))
ax.add_artist(anchored_text)
ax.set_title(f'{label}', fontsize='xx-large', color = 'red', fontweight='bold')
ax.set_xlabel('tan($\\theta_x$)')
ax.set_ylabel('tan($\\theta_y$)')
xlim, ylim = [-1., 1.], [-1.5, 1.5]
ax.set_xlim(xlim)
ax.set_ylim(ylim)

##region GAP @ top of storage well 
hgap = 530 #mm
htel = h
dist = dy_tel_storwell
tthetax_gap_min, tthetax_gap_max = -radius_well_bot/dist, radius_well_bot/dist 
tthetay_gap_min, tthetay_gap_max =  -(htel-hgap)/dist, -htel/dist
kwargs = {"linewidth": 1., "linestyle":"dashed", "color":"red"}
fontdict = {'family': 'serif', 'color':  'red', 'weight': 'normal','size': 'large', 'ha':'center',}
ax.text((tthetax_gap_max+tthetax_gap_min)/2,(tthetay_gap_max+tthetay_gap_min)/2, s="gap",  fontdict=fontdict) #default is data coordinates

nwell, nbar = 3, 10
hbarrel, hpot = 617, 285 #mm
radius_barrel, radius_pot = 180, 110 #mm
h0 = (htel-hgap)
# r = plt.Rectangle(xy=(0, 0), width=0.1, height=0.25, color="red", visible=True)
# ax.add_patch(r)
for l in np.arange(0, nwell): #storage well	
	hb = h0
	hp = hb-(hbarrel-hpot)
	dist1 = np.sqrt(dist**2 + (l*spacing_well)**2) 
	tthetax_well_min, tthetax_well_max = -(radius_well_bot+l*spacing_well)/dist1, (radius_well_bot-l*spacing_well)/dist1 
	ax.axvline(tthetax_well_min, **kwargs)
	ax.axvline(tthetax_well_max, **kwargs)
	for k in np.arange(0, nbar): #barrels
		tthetax_bar_min, tthetax_bar_max = -(radius_barrel+l*spacing_well)/dist1, (radius_barrel-l*spacing_well)/dist1
		tthetay_bar_min, tthetay_bar_max =  -(hb-(k)*hbarrel)/dist1, -(hb-(k+1)*hbarrel)/dist1
		tthetax_pot_min, tthetax_pot_max = -(radius_pot+l*spacing_well)/dist1, (radius_pot-l*spacing_well)/dist1
		tthetay_pot_min, tthetay_pot_max =  -(hp-k*(hbarrel))/dist1, -(hp-hpot-k*hbarrel)/dist1

		xr, yr = tthetax_bar_min, tthetay_bar_min
		w, h = abs(tthetax_bar_max-tthetax_bar_min), abs(tthetay_bar_max-tthetay_bar_min)
		# print(k, tthetax_bar_min, tthetax_bar_max, tthetay_bar_min, tthetay_bar_max)
		if (ylim[0] <tthetay_bar_max) & (tthetay_bar_max < ylim[1]): 
			rb = plt.Rectangle(xy=(xr,yr), width=w, height=h,  facecolor="blue", edgecolor='darkgrey', linewidth=2., alpha=0.1)
			ax.add_patch(rb)
			xt, yt = (tthetax_bar_min+tthetax_bar_max)/2, (tthetay_bar_min+tthetay_bar_max)/2
			ax.text(xt, yt, s=f"B{k+1}", fontdict=fontdict)
	
		# if k < 2:
			xr, yr = tthetax_pot_min, tthetay_pot_min
			w, h = abs(tthetax_pot_max-tthetax_pot_min), abs(tthetay_pot_max-tthetay_pot_min)
			# print(k, tthetax_bar_min, tthetax_bar_max, tthetay_bar_min, tthetay_bar_max)
			if (ylim[0] <tthetay_pot_max) & (tthetay_pot_max < ylim[1]): 
				rp = plt.Rectangle(xy=(xr,yr), width=w, height=h,  facecolor="red", edgecolor='darkred', linewidth=1., alpha=0.1, hatch='//')
				ax.add_patch(rp)


ax.invert_yaxis()


foutname = f"thetaxy_{ec.detector_ensemble.name}_{basename}.png"
fout = dout / foutname
fig.savefig(fout, dpi=300)
print(f'Save figure {fout}')
plt.close()
Nmuf = h2d


if sim_name=='INB72' :
	acc_geo = dphi * s_gen * Nmuf/Nmu0  #cm2.sr

	print(f"acc_geo \n\t (mean, std) = {np.mean(acc_geo):.2f}, {np.std(acc_geo):.2f} cm2.sr\n\t (min, max) = {np.min(acc_geo):.2f}, {np.max(acc_geo):.2f} cm2.sr")

	fig, ax = plt.subplots(figsize=(10,10))
	im = ax.pcolormesh(XC, YC, acc_geo, norm=LogNorm(vmin=1, vmax=np.max(h2d)))
	ax.invert_yaxis()
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	cbar = fig.colorbar(im, cax=cax, extend='max')
	anchored_text = AnchoredText(f"entries: {ntrack:.2e}", loc="upper right", frameon=True, prop=dict(fontsize="medium", color="red"))
	ax.add_artist(anchored_text)
	ax.set_title(f'{label}', fontsize='xx-large', color = 'red', fontweight='bold')
	ax.set_xlabel('tan($\\theta_x$)')
	ax.set_ylabel('tan($\\theta_y$)')
	ax.set_xlim(-1,1)
	ax.set_ylim(-1,1)
	foutname = f"acceptance_geo_{ec.detector_ensemble.name}_{basename}.png"
	fout = dout / foutname
	fig.savefig(fout, dpi=300)
	print(f'Save figure {fout}')
	plt.close()

	fig = plt.figure(figsize=(12, 7))
	ax = fig.add_subplot(projection='3d')
	im = ax.plot_surface(
				XC, YC, acc_geo,
				# **kwargs
			)
	cbar = plt.colorbar(im,  shrink=0.5, orientation="vertical")
	cbar.ax.tick_params(labelsize=12)
	cbar.set_label(label='Acceptance [cm².sr]', size=14)
	ax.set_title(f'{label}', fontsize='xx-large', color = 'red', fontweight='bold')
	ax.set_xlabel('tan($\\theta_x$)')
	ax.set_ylabel('tan($\\theta_y$)')
	ax.set_xlim(-1,1)
	ax.set_ylim(-1,1)
	ax.invert_yaxis()
	foutname = f"acceptance_geo_3D_{ec.detector_ensemble.name}_{basename}.png"
	fout = dout / foutname
	fig.savefig(fout, dpi=300)
	print(f'Save figure {fout}')
	plt.close()