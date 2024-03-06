'''
Adapted from Paul's script
'''

from tomomuLib_Telemodel import Detecteur
from tomomuLib_Telemodel import PlanDetection
from tomomuLib_Telemodel import Hodoscope


def def_INB72(arr_pos, arr_size, points_tel:bool=True) -> Hodoscope:

    #####   Definition of Tel ################################################
    #Each of the detector is part of its own plane of detection

    colors_tel = ['orange', 'green', 'brown', 'purple']
    positions_tel = [[0,0,0],[0,0,50],[0,0,130],[0,0,180]]

    tel_det = [None] * 4
    tel_planes = [None] * 4
    for i in range(4):
        tel_det[i] = Detecteur(f'Tel_det{i}', colors_tel[i])
        if points_tel==True :
            tel_det[i].add_events(arr_pos[:, 2*i], arr_pos[:, 2*i+1], arr_size[:, 2*i], arr_size[:, 2*i+1])
        tel_planes[i] = PlanDetection(f'Tel_{i}', [tel_det[i]], [positions_tel[i]])

    tel = Hodoscope('INB72', tel_planes)

    return tel

