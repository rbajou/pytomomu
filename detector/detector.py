#!/usr/bin/python3
# -*- coding: utf-8 -*-

from dataclasses import dataclass, field
from abc import abstractmethod
from typing import List, Tuple, Union
from enum import Enum, auto
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from math import copysign
import argparse

#package module(s)
from utils import tools


@dataclass
class StripLayer:
    index : int
    nstrips : int
    length : float #mm
    pitch : float #mm
    
    def __str__(self):
        coord = 'X'
        if self.index%2 != 0 : coord = 'Y' 
        sout = f"{coord}{self.index}:\n\t- nstrips = {self.nstrips} \n\t- length = {self.length} mm \n\t- pitch = {self.pitch} microns"
        return sout

class Parent(Enum):
    veto = auto()
    tel = auto()

@dataclass
class Detector:
    version : int
    nfam : int
    nch : int = field(default=61)
    side : float = field(default=546.)#mm 
    thickness : float = field(default=10.) #mm
    color : str = field(default='blue')

    def __post_init__(self):
        self._name = None
        self._position = None
        self._layers = None
        self._index = None
        self._parent = None

    def __str__(self):
        sout = f"MGv{self.version}:\n\t- layers = {self._layers} \n\t- n multiplexing families = {self.nfam}"
        return sout

    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, value:str):
        self._name = value

    @property
    def position(self):
        """xyz position."""
        return self._position

    @position.setter
    def position(self, value:np.ndarray):
        if value.shape != (3,) : raise ValueError("Incorrect value array shape, must be (3,).")
        self._position = value

    @property
    def strip_layers(self):
        return self._layers
    
    @strip_layers.setter
    def strip_layers(self, value:tuple):
        """Tuple of 'StripLayer' objects (X, Y)"""
        self._layers = value

    @property
    def index(self):
        return self._index
    
    @index.setter
    def index(self, value:int):
        self._index = value

    @property
    def parent(self):
        return self._parent
    
    @parent.setter
    def parent(self, value:str):
        self._parent = value

class DetectorEnsemble: 

    def __init__(self, name:str, detectors:List[Detector]) -> None:
        self.name = name
        self.detectors = detectors
        self._position = None
        self.npanels = len(self.detectors)
        self._xy_offset = None

    @property
    def position(self):
        return self._position
    
    @position.setter
    def position(self, value:np.ndarray):
        if value.shape != (3,) : raise ValueError("Incorrect value array shape, must be (3,).")
        self._position = value

    @property
    def xy_offset(self, value):
        self._xy_offset = value
        
        return self._xy_offset

    @xy_offset.setter
    def xy_offset(self, value:float):
        """Offset to detector central position in mm."""
        for d in self.detectors: 
            p = d.position.copy()
            for c in [0,1]: 
                d.position[c] = p[c] + copysign(1,p[c])*value

    def plotXY(self,  ax, offset:np.ndarray=np.zeros(3), rotx:float=0., roty:float=0., rotz:float=0., npt:int = 100, **kwargs):
        '''
        Input:
        - 'ax' (plt.Axes) 
        - 'offset' xyz (np.ndarray)
        '''
        for i, p in enumerate(self.detectors):
            xlay, ylay = p.strip_layers
            xp, yp = xlay.pitch, ylay.pitch
            xn, yn = xlay.nstrips, ylay.nstrips
            sx, sy = xp*xn, yp*yn
            pos = p.position + offset
            x, y = np.linspace(pos[0]-sx/2, offset[0]+sx/2, npt), np.linspace(pos[1]-sy/2, offset[1]+sy/2, npt)
            ax.plot(x, y, **kwargs)

    def plotXZ(self,  ax, offset:np.ndarray=np.zeros(3), rotx:float=0., roty:float=0., rotz:float=0., npt:int = 100, **kwargs):
        '''
        Input:
        - 'ax' (plt.Axes) 
        - 'offset' xyz (np.ndarray)
        '''
        zticks=[]
        for i, p in enumerate(self.detectors):
            xlay, _ = p.strip_layers
            xp = xlay.pitch
            xn = xlay.nstrips
            sx = xp*xn
            pos = p.position + offset
            x = np.linspace(pos[0]-sx/2, pos[0]+sx/2, npt)
            z = np.ones(x.shape)*p.position[2]
            ax.plot(x, z, **kwargs)

            side_box, thickness_box = p.side, p.thickness
            rb = plt.Rectangle(xy=(pos[0]-side_box/2,p.position[2]-thickness_box/2), width=side_box, height=thickness_box,  facecolor="none", edgecolor='darkgrey', linewidth=2., alpha=0.1)
            ax.add_patch(rb)

            zticks.append(z[-1])
    
    def plotYZ(self,  ax, offset:np.ndarray=np.zeros(3), rotx:float=0., roty:float=0., rotz:float=0., npt:int = 100,**kwargs):
        '''
        Input:
        - 'ax' (plt.Axes) 
        - 'offset' xyz (np.ndarray)
        '''
        zticks=[]
        for i, p in enumerate(self.detectors):
            _, ylay = p.strip_layers
            yp = ylay.pitch
            yn = ylay.nstrips
            sy = yp*yn
            pos = p.position + offset
            y = np.linspace(pos[1]-sy/2, pos[1]+sy/2, npt)
            z = np.ones(y.shape)*p.position[2]
            ax.plot(y, z, **kwargs)

            side_box, thickness_box = p.side, p.thickness
            rb = plt.Rectangle(xy=(pos[0]-side_box/2,p.position[2]-thickness_box/2), width=side_box, height=thickness_box,  facecolor="none", edgecolor='darkgrey', linewidth=2., alpha=0.1)
            ax.add_patch(rb)

            zticks.append(z[-1])

        

    def plot3D(self, ax, offset:np.ndarray=np.zeros(3), rotx:float=0., roty:float=0., rotz:float=0., npt:int = 100, **kwargs):
        '''
        Input:
        - 'ax' (plt.Axes3D) : e.g 'ax = fig.add_subplot(111, projection='3d')'
        - 'offset' xyz (np.ndarray)
        '''
        zticks=[]
        mrotx = np.array([[1,0,0], [0,np.cos(rotx),-np.sin(rotx)], [0,np.sin(rotx),np.cos(rotx)]]) 
        mroty = np.array([[np.cos(roty),0,np.sin(roty)], [0,1,0],  [-np.sin(roty),0,np.cos(roty)]]) 
        mrotz = np.array([[np.cos(rotz),-np.sin(rotz), 0], [np.sin(rotz),np.cos(rotz), 0], [0,0,1]]) 

        for i, p in enumerate(self.detectors):

            xlay, ylay = p.strip_layers
            xp, yp = xlay.pitch, ylay.pitch
            xn, yn = xlay.nstrips, ylay.nstrips
            sx, sy = xp*xn, yp*yn
            x, y = np.linspace(offset[0]-sx/2, offset[0]+sx/2, npt), np.linspace(offset[1]-sy/2, offset[1]+sy/2, npt)
            z = np.ones(x.shape)*p.position[2]
            # xyz = np.vstack((x,y,z)).T
            # xyz = (mrotz @  xyz.T)
            # xyz = (mrotx @  xyz)
            # xyz = (np.dot(np.dot(xyz, mrotx), mrotz ))
            # xyz = np.dot( xyz, mrotx )
            # xyz = np.matmul( xyz, mrotz )
            # xyz = (np.dot(xyz, np.dot(mrotx, mrotz)))
            # xyz = xyz @ mrotz
            # xyz = np.dot(xyz, mrotz)
            (X, Y), Z = np.meshgrid(x, y), np.tile(z, (len(z), 1))  
            Y, Z = np.einsum('ij,klj->ikl', mrotx[1:,1:],  np.dstack([Y, Z])) 
            X, Y = np.einsum('ij,klj->ikl', mrotz[:-1,:-1], np.dstack([X, Y])) 
            ax.plot_surface(X,Y,Z, **kwargs)
            zticks.append(z[-1])
          
        ax.set_xlabel("X [mm]", labelpad=20)
        ax.set_ylabel("Y [mm]", labelpad=20)
        ax.tick_params(axis='both', which='major', pad=15)  
        # ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        # ax.invert_zaxis()
        ax.set_zticks(zticks)
        # ax = fig.add_axes(tools.MyAxes3D(ax, 'l'))
        # ax.annotate("Z", xy=(0.5, .5), fontsize='x-large', xycoords='axes fraction', xytext=(0.04, .78),)
        
        return ax



class VetoPlane(DetectorEnsemble):

    def __init__(self, name, detectors):
        DetectorEnsemble.__init__(self, name, detectors)


class Telescope(DetectorEnsemble):
    
    def __init__(self, name, detectors):
        DetectorEnsemble.__init__(self, name, detectors)


class Bench:

    def __init__(self, name:str, detector_ensemble:List[DetectorEnsemble]):
        self.name = name
        self.detector_ensemble = detector_ensemble
        self.veto_plane = [ de for de in detector_ensemble if isinstance(de, VetoPlane)]
        self.telescope = [ de for de in detector_ensemble if isinstance(de, Telescope)]
        self.detectors = [ vp.detectors for vp in self.veto_plane ]
        self.detectors.extend([ te.detectors for te in self.telescope ])
        self.detectors = [item for row in self.detectors for item in row] #flatten list of lists

    def plot3D(self, fig, ax, offset:np.ndarray=np.zeros(3), rotx:float=0., roty:float=0., rotz:float=0., **kwargs):
        '''
        Input:
        - 'ax' (plt.Axes3D) : e.g 'ax = fig.add_subplot(111, projection='3d')'
        - 'offset' xyz (np.ndarray)
        '''
        zticks=[]
        mrotx = np.array([[1,0,0], [0,np.cos(rotx),-np.sin(rotx)], [0,np.sin(rotx),np.cos(rotx)]]) 
        mroty = np.array([[np.cos(roty),0,np.sin(roty)], [0,1,0],  [-np.sin(roty),0,np.cos(roty)]]) 
        mrotz = np.array([[np.cos(rotz),-np.sin(rotz), 0], [np.sin(rotz),np.cos(rotz), 0], [0,0,1]]) 
        
        fontdict = {'family': 'serif', 'color':  'red', 'weight': 'normal','size': 'large', 'ha':'center',}
     
        npt = 100

        for i, p in enumerate(self.detectors):

            xlay, ylay = p.strip_layers
            xp, yp = xlay.pitch, ylay.pitch
            xn, yn = xlay.nstrips, ylay.nstrips
            sx, sy = xp*xn, yp*yn
            pos = p.position + offset
            x, y = np.linspace(pos[0]-sx/2, pos[0]+sx/2, npt), np.linspace(pos[1]-sy/2, pos[1]+sy/2, npt)
            z = np.ones(x.shape)*p.position[2]
            (X, Y), Z = np.meshgrid(x, y), np.tile(z, (len(z), 1))  
            Y, Z = np.einsum('ij,klj->ikl', mrotx[1:,1:],  np.dstack([Y, Z])) 
            X, Y = np.einsum('ij,klj->ikl', mrotz[:-1,:-1], np.dstack([X, Y])) 
            ax.plot_surface(X,Y,Z, **kwargs)
            zticks.append(z[-1])
            # ax.text(x=pos[0], y=pos[1], z=pos[2], s=f"{i}",  fontdict=fontdict) #default is data coordinates


        ax.set_xlabel("X [mm]", labelpad=20)
        ax.set_ylabel("Y [mm]", labelpad=20)
        ax.tick_params(axis='both', which='major', pad=15)  
        # ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        # ax.invert_zaxis()
        ax.set_zticks(zticks)
        # ax = fig.add_axes(tools.MyAxes3D(ax, 'l'))
        # ax.annotate("Z", xy=(0.5, .5), fontsize='x-large', xycoords='axes fraction', xytext=(0.04, .78),)
        
        return ax


#IZEN (Bench): VETO BOTTOM, VETO TOP, TELESCOPE
z_vb, z_vt, z0_tel = 120, 342, 3974  #mm
x_c = y_c = 250 #mm
xy_off = 23. #mm
det_xyz_vb = np.array([[-x_c,-y_c,z_vb], [x_c,-y_c,z_vb], [-x_c,y_c,z_vb], [x_c,y_c,z_vb]])
det_xyz_vt = np.array([[-x_c,-y_c,z_vt], [x_c,-y_c,z_vt], [-x_c,y_c,z_vt], [x_c,y_c,z_vt]])
det_xyz_tel = np.array([[0,0,3974],[0.151,0.023,4044-0.0894],[0,0,4124],[0.050,0.0594,4194-0.11]])

nstrips_mgv4 = 732
strip_length_mgv4 = 500. #mm
pitch_mgv4 = 0.683 #mm
kwargs_mgv4 = {"nstrips" : nstrips_mgv4, "length": strip_length_mgv4, "pitch" : pitch_mgv4}
mgv4 = Detector(version=4, nfam=12)
detectors_vb, detectors_vt, detectors_tel = [], [], [] #[deepcopy(mgv4) for _ in range(4)], [deepcopy(mgv4) for _ in range(4)], [deepcopy(mgv4) for _ in range(4)]

for ilay in range(0, 24, 2): 
    strip_layer_x = StripLayer(index=ilay, **kwargs_mgv4)
    strip_layer_y = StripLayer(index=ilay+1, **kwargs_mgv4)
    det = deepcopy(mgv4)
    det.strip_layers = (strip_layer_x, strip_layer_y)
    if ilay in range(0,7):
        det.parent = Parent.veto.name
        detectors_vb.append(det)
    elif ilay in range(8,15):
        det.parent = Parent.veto.name
        detectors_vt.append(det)
    elif ilay in range(16,23):
        det.parent = Parent.tel.name
        detectors_tel.append(det)
    else : pass

for i, (dvb,dvt,dtel)  in enumerate(zip(detectors_vb, detectors_vt, detectors_tel)):  
    dvb.name, dvt.name, dtel.name = f'dvb_{i}', f'dvt_{i}', f'dtel_{i}'
    dvb.position, dvt.position, dtel.position = det_xyz_vb[i], det_xyz_vt[i], det_xyz_tel[i]
    dvb.index, dvt.index, dtel.index = i, i+len(detectors_vb), i+len(det_xyz_vb)+len(detectors_vt)

VB = VetoPlane(name="VB", detectors = detectors_vb)
VB.xy_offset = xy_off
VT = VetoPlane(name="VT", detectors = detectors_vt)
VT.xy_offset = xy_off
TEL = Telescope(name="TEL", detectors = detectors_tel)
IZEN = Bench(name="IZEN", detector_ensemble = [VB,VT,TEL])


#INB72 (Telescope)
inb72 = Detector(version=5, nfam=5, side=200.)
npan = 4 
nlay = 2*npan
detectors_inb72 = []
nstrips_inb72 = 305
strip_length_inb72 = 170. #mm
pitch_inb72 = 0.557 #mm
kwargs_mgv4 = {"nstrips" : nstrips_inb72, "length": strip_length_inb72, "pitch" : pitch_inb72}
for ilay in range(0, nlay, 2): 
    strip_layer_x = StripLayer(index=ilay, **kwargs_mgv4)
    strip_layer_y = StripLayer(index=ilay+1, **kwargs_mgv4)
    det = deepcopy(inb72)
    det.parent = Parent.tel.name
    det.strip_layers = (strip_layer_x, strip_layer_y)
    detectors_inb72.append(det)
  
det_xyz_inb72 = np.array([[0,0,0],[0,0,50],[0,0,130],[0,0,180]])

for i, dtel  in enumerate(detectors_inb72):  
    dtel.name = f'inb72_pan{i}'
    dtel.position = det_xyz_inb72[i]
    dtel.index = i

INB72 = Telescope(name="INB72", detectors = detectors_inb72)


DICT_DET = {'INB72': INB72, 'IZEN': IZEN}


def str2detector(v):
    '''
    Convert 'str' to 'Telescope' or 'Bench' object type
    '''
   
    if isinstance(v, Telescope) or isinstance(v, Bench):
       return v

    if v in list(DICT_DET.keys()):
        return DICT_DET[v]
    elif v in [ k.lower() for k in list(DICT_DET.keys())]:
        return DICT_DET[v.upper()]
    elif v in [f"tel_{k}" for k in list(DICT_DET.keys()) ]:
        return DICT_DET[v[4:]]
    elif v in [f"tel_{k.lower()}" for k in list(DICT_DET.keys()) ]:
        return DICT_DET[v[4:].upper()]
    else:
        raise argparse.ArgumentTypeError('Input detector does not exist.')



if __name__ == "__main__":
    
    #Draw tel
    tel = IZEN.telescope[0]
    vetobot= IZEN.veto_plane[0]
    vetotop= IZEN.veto_plane[1]

    fig, ax = plt.subplots(figsize=(12,8))
   
    kwargs = dict(alpha=1., color='greenyellow')
    tel.plotXZ(ax=ax, **kwargs) 
    vetobot.plotXZ(ax=ax, **kwargs) 
    vetotop.plotXZ(ax=ax, **kwargs) 
    ax.set_xlim(-500, 500)
    ax.set_ylim(0, 4500)
    fout = f"{tel.name}_xz.png"
    fig.savefig(fout)
    # plt.show(block=True)
    print(f'Save {fout}')

    # tel = IZEN
    # fig = plt.figure(figsize=(12,8))
    # ax = fig.add_subplot(111, projection='3d' )
    # kwargs = dict(alpha=0.2, color='greenyellow', edgecolor='none')
    # tel.plot3D( ax=ax, **kwargs) 
    # ax.set_xlim(-500, 500)
    # ax.set_ylim(-500, 500)
    # ax.set_zlim(0, 4500)
    # fout = f"{tel.name}.png"
    # fig.savefig(fout)
    # plt.show(block=True)
    # print(f'Save {fout}')



