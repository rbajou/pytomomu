import numpy as np
import matplotlib.pyplot as plt
import uproot
from pathlib import Path

#tomomuLib libraries
from pyTomomu_Tomography import Tomography_py


def display_tomo(title, opacity, muonflux, os_muonflux, distmat_red, distmat_col, distmat_ysm):
    fig, (ax1, ax2, ax2b, ax3, ax4, ax5) = plt.subplots(1, 6, figsize=(20,5))

    ax1.matshow(opacity, vmax = 60)
    ax1.set_title('opacity')

    ax2.matshow(muonflux)
    ax2.set_title('muonflux')

    ax2b.matshow(os_muonflux)
    ax2b.set_title('o.sky muonflux')

    ax3.hist(distmat_red, 50)
    ax3.set_title('distmat_red')

    ax4.hist(distmat_col, 50)
    ax4.set_title('distmat_col')

    ax5.hist(distmat_ysm, 50)
    ax5.set_title('distmat_ysm')

    fig.tight_layout()
    plt.savefig(title+'.png')
    # plt.show()
    plt.close()


def display_2(data, title):
    plt.figure(figsize=(16,5))

    plt.subplot(131)
    plt.subplots_adjust(bottom=0.1, left=0.01, top=0.9)
    plt.imshow(data.sum(0).T[-1::-1], cmap='jet')#, vmin=0, vmax=vmax)
    plt.title('integral(X)')

    plt.subplot(132)
    plt.imshow(data.sum(1).T[-1::-1], cmap='jet')#, vmin=0, vmax=vmax)
    plt.title('integral(Y)')

    plt.subplot(133)
    plt.imshow(data.sum(2), cmap='jet')#, vmin=0, vmax=vmax)
    plt.title('integral(Z)')
 
    plt.subplots_adjust(bottom=0.1, right=0.83, top=0.9)

    #plt.show()
    plt.savefig(title+'_projections.png')
    plt.close()



if __name__ == "__main__":
   
   
    home = Path.home()
    data_path = home / "Data" / "tomomu" / "inb72" 
    ana_path = data_path / "analyse" 
    final_path = data_path / "final"  
    final_path.mkdir(parents=True, exist_ok=True)


    # Load the positions and orientations
    file = open('/local/home/pl275970/IZEN_analysis/telpos_summary.txt')
    #file = open('/local/home/bl268739/Documents/izen_analysis/telpos_summary.txt')  # TODO Paul

    config = "P1_Q1"

    telescopes, positions, quart, x, y, z, alpha, beta, gamma, days = [],[],[],[],[],[],[],[],[],[]

    head = True
    for line in file:
        if not head:
            name,x_,y_,z_,alpha_, beta_, gamma_, days_ = line.split(',')
            tel, pos, quartier = name.split('_')
            telescopes.append(tel)
            positions.append(pos[1:])
            quart.append(quartier[1:])
            x.append(float(x_))
            y.append(float(y_))
            z.append(float(z_))
            alpha.append(float(alpha_))
            beta.append(float(beta_))
            gamma.append(float(gamma_))
            days.append(float(days_))
        head = False

    print(len(telescopes))

    opensky_sim_time = 100


    ## Do the tomographies

    #nvx, nvy, nvz = 16,12,14

    nvx, nvy, nvz = 20,20,20 #80,60,70 #

    Xv0, Xv1 = -0.8, 0.8
    Yv0, Yv1 = -0.8, 0.8
    Zv0, Zv1 = 0.4, 3.0

    nttr = 100
    ntpr = 100

    #TTr0, TTr1 = -0.06, +0.19
    #TPr0, TPr1 = -0.06, +0.19

    TTr0, TTr1 = -0.2, 0.2
    TPr0, TPr1 = -0.2, 0.2

    nvx, nvy, nvz = int(nvx), int(nvy), int(nvz)
    nv = 'nv_' + str(nvx)+'_'+str(nvy)+'_'+str(nvz)


    tomographies_processed = []

    qx = [-0.125, -0.125, 0.125, 0.125]
    qy = [-0.125, 0.125, -0.125, 0.125]

    basename = "INB72_2024_run1"
    dout = data_path / "tomography"
    dout.mkdir(parents=True, exist_ok=True) 

    for i in range(len(telescopes)):
        tel_pos = telescopes[i]+'_P'+positions[i]+'_Q'+quart[i]
        config = f"P{positions[i]}_Q{quart[i]}"
        duration_days = 1 #days[i]

        # path = '/local/home/pl275970/data/final3_files/MID/'
        #path = '/local/home/bl268739/Documents/izen_analysis/'  # TODO Paul
        filename = f"{basename}_final.root"
        fin = final_path / filename #'G4TomoDet_IZEN_MID_OsBall_002_' + config +'_final2.root'
        
        #filename = 'sim_'+telescopes[i]+'_P'+positions[i]+'_d'+str(days[i])+'.root'
        try :
            data_tr = uproot.open(str(file))['Tmeta']
            max_ts = np.array(data_tr['max_time'])
            duration_days = max_ts.sum()/60/60/24
        except:
            print('! Cannot find file or open Tmeta tree for ', tel_pos)

        if duration_days > 0:
            assert beta[i]==0, 'Not prepared for that yet'

            tomography = Tomography_py()
            tomography.input_data = str(fin)
            tomography.input_data_is_simulated = False
            tomography.input_opensky_is_simulated = False

            tomography.telescope_position_x = x[i]
            tomography.telescope_position_y = y[i]
            tomography.telescope_position_z = z[i]
            tomography.telescope_position_alpha_rad = alpha[i]
            tomography.telescope_position_beta_rad = beta[i]
            tomography.telescope_position_gamma_rad = gamma[i]

            tomography.nvx = nvx
            tomography.nvy = nvy
            tomography.nvz = nvz


            tomography.nttr = nttr
            tomography.ntpr = ntpr

            tomography.Xv0 = Xv0
            tomography.Xv1 = Xv1

            tomography.Yv0 = Yv0
            tomography.Yv1 = Yv1

            tomography.Zv0 = Zv0
            tomography.Zv1 = Zv1

            tomography.TTr0 = TTr0
            tomography.TTr1 = TTr1
            tomography.TPr0 = TPr0
            tomography.TPr1 = TPr1

            tomography.geometry_simulation_days = duration_days
            tomography.geometry_tomography_days = duration_days

            tomography.max_coeff_nb_divisor = 1


            # tomography.output_file = '/local/home/pl275970/data/tomography_files/MID/OsBall/G4TomoDet_IZEN_MID_OsBall_002_'+config+'_tomography.root'
            #tomography.output_file = path + 'tomography.root'  # TODO Paul
            tomography.output_file = str(dout /  f"{basename}_tomography.root")

            # tomography.input_opensky = str(flux_path / 'G4TomoDet_IZEN_MID_OsBall_002_'+config+'_OS_final2.root' ) 

            tomography.opensky_simulation_days = duration_days # days[i]
            tomography.opensky_tomography_days = duration_days # days[i]

            tomography.compute_opacities = 1  # TODO Paul
            
            tomography.temp_izen_behavior = 1
            tomography.temp_izen_x_tracker = qx[int(quart[i])-1] # TODO Paul
            tomography.temp_izen_y_tracker = qy[int(quart[i])-1] # TODO Paul
            tomography.temp_izen_z_tracker = 3.974
            tomography.temp_izen_halfsize_tracker = 0.125
            tomography.compute_distances = 1
            r=tomography.run()

            if r==0:
                tomographies_processed.append(str(tomography.output_file)[2:-1])
                root_file = uproot.open(tomographies_processed[-1])
                opacity = np.array(root_file['T']['opacity'].array()).reshape((nttr,ntpr))
                muonflux = np.array(root_file['T']['muonflux'].array()).reshape((nttr,ntpr))
                os_muonflux = np.array(root_file['T']['os_muonflux'].array()).reshape((nttr,ntpr))
                try:
                    distmat_red = np.array(root_file['T']['distmat_red'].array()).flatten()
                    distmat_col = np.array(root_file['T']['distmat_col'].array()).flatten()
                    distmat_ysm = np.array(root_file['T']['distmat_ysm'].array()).flatten()
                except:
                    distmat_red = np.zeros(1)
                    distmat_col = np.zeros(1)
                    distmat_ysm = np.zeros(1)
                display_tomo(tomographies_processed[-1], opacity, muonflux,os_muonflux, distmat_red, distmat_col, distmat_ysm)
                
                
                            
                ray_start=0
                ray_stop=10000
                
                total_distances = np.zeros(20*20*20)
                max_distances = np.zeros(20*20*20)
                for i in range(distmat_ysm[ray_start], distmat_ysm[ray_stop],1):#distmat_col.size):
                    voxel_i = distmat_col[i]
                    distance = distmat_red[i]
                    
                    # print(voxel_i, distance)
                    total_distances[voxel_i] += distance
                    max_distances[voxel_i] = max(max_distances[voxel_i], distance)


                display_2(total_distances.reshape((20,20,20), order='F'), tomographies_processed[-1]+'_rays')  # TODO Paul

