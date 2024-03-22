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
import pandas as pd
import glob
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")
import subprocess
#package module(s)
from detector import INB72
from survey import set_survey
from real.rootfile import RootFile



params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large',
         'axes.labelpad':10,
         'mathtext.fontset': 'stix',
         'font.family': 'STIXGeneral',
         'axes.grid' : True 
         }
plt.rcParams.update(params)

# class Monitoring(RootFile):

#     def __init__(self, file:Union[Path, str], **kwargs) -> None:
      
#         RootFile.__init__(self, file=file)
#         self.open(treename = "HV_Mon")
#         print('rootfile branches : ', self._branch_names)



class Monitoring :
    '''
    Class to fecth monitoring files recorded and saved on nas within a given run period (ti_mon, tf_mon)
    '''

    def __init__(self, **kwargs) -> None:

        self.files_glob = []
        self.files_read = []
        self.files_corrupt = []
        self.nevt_valid = 0
        self.dict_nevt = {}

    def fetch_files(self, path:Path, str_wildcard:str="*HV_mon.root", regex:str=r'(\d{8})_(\d{2}H\d{2})'):
        '''
        
        '''
        files_glob = glob.glob(str(path / str_wildcard))
        # print(f"files_glob ({len(files_glob)}) = {files_glob}")
        self.ti, self.tf = np.zeros((len(files_glob),)), np.zeros((len(files_glob),))
        self.t_range_all = np.zeros((len(files_glob), 2))  
        self.t_gap =       np.zeros(len(files_glob)-1,)
        self.files_glob = files_glob

        for i, f in enumerate(files_glob):
            if i > 10 : continue
            # Define the regex pattern to extract date and time
            pattern = re.compile(regex)
            # Match the pattern in the filename
            match = pattern.search(f)

            if match:
                extracted_date = match.group(1)
                extracted_time = match.group(2)
                # print("Extracted date:", extracted_date)
                # print("Extracted time:", extracted_time)
                formatted_time = extracted_time.replace('H', ':')
                datetime_str = f"{extracted_date} {formatted_time}"

                try:
                    # Convert the combined date and time to a Unix timestamp
                    ti_mon_mon = int(datetime.strptime(datetime_str, "%Y%m%d %H:%M").timestamp())
                    # unix_timestamp=$(stat -c "%Y" ${indir}/${prefix}_run${irun}_analyse.root)
                    tf_mon_mon = Path(f).lstat().st_mtime # subprocess.run(['bash', f'stat -c "%Y" ${f}'], check=True, stdout=subprocess.PIPE).stdout.decode('utf-8')
                    self.ti[i], self.tf[i] =  ti_mon_mon, tf_mon_mon
                    self.t_range_all[i, :] = (ti_mon_mon, tf_mon_mon)
                    self.files_read.append(f)
                except ValueError:
                    print("Error converting date and time to Unix timestamp.")
            else:
                print("Date and/or time not found in the filename.")

        if len(self.files_read) > 0:             
            arg_nnull_sort = np.argsort(self.ti[self.ti != 0])
            self.ti, self.tf, self.t_range_all = self.ti[arg_nnull_sort], self.tf[arg_nnull_sort], self.t_range_all[arg_nnull_sort]
            self.files_read = np.asarray(self.files_read)[arg_nnull_sort]
            self.t_gap = self.ti[1:] - self.tf[:-1] #time gap between i-th run and (i-1)-th run

    def read_rootfile(self, file:Union[str, Path], treename:str="HV_Mon", n_events:int=int(1e9))->dict:
        
        rootfile = RootFile(file)
        content = {}
        try : 
            rootfile.open(treename=treename)
            rootfile.get_content(n_events=n_events)
            content = rootfile.content
    
        except:
            print(f"Failed to open rootfile {file}")
            print(f"Fail to fetch content of {file}")
    
        finally:
            return content
    
    def get_t_range_from_tree(self, files:list, treename:str="HV_Mon", n_events:int=int(1e9)):
        '''
        
        '''
        t_range = np.zeros(shape=(len(files), 2))  
       
        for i, f in enumerate(files): 
            content = self.read_rootfile(f, treename, n_events)
            try : 
                t = content['time']
                t_range[i] = np.array([np.nanmin(t), np.nanmax(t)])
            except:
                print(f"No time branch found in {f}")    
                continue
        
        return t_range

    def to_df(self, files_in:list=None, file_out:Union[Path,str]=None):
        '''
        Save list monitoring files as dateframe with time periods and datetime start index.
        '''
       
        files_read = self.files_read
        ti, tf, t_gap = self.ti, self.tf, self.t_gap
     
        if len(files_in) >0 :
            files_read = files_in
            ti, tf = self.get_t_range_from_tree(files_in).T
            arg_nnull_sort = np.argsort(ti[ti != 0])
            files_read = np.asarray(files_read)[arg_nnull_sort]
            ti = ti[arg_nnull_sort]
            tf = tf[arg_nnull_sort]
            t_gap = ti[1:] - tf[:-1]

        self.df = pd.DataFrame(data={"name": [Path(f).name for f in files_read ], "date_init":pd.to_datetime(ti, unit='s'), "date_fin":pd.to_datetime(tf, unit='s'), "ti": ti ,"tf": tf})
        self.df['gap'] = np.zeros_like(ti)
        self.df['gap'].iloc[1:] = t_gap
        self.df['gap'][(self.df['ti'] == 0)&(self.df['tf'] == 0)] = 0
        if file_out : 
            self.df.to_csv(file_out)
            print(f"Save df {file_out}")
    
    def get_run_files(self, t_range:tuple[float, float]=(0, 1e10)):
        '''
        
        '''
        ti_run, tf_run = t_range
        self.files_run = []
        for _, (tr, f) in enumerate(zip(self.t_range_all, self.files_read)): 
            ti_mon, tf_mon = tr
            if (ti_mon < ti_run) & (tf_mon <= tf_run): self.files_run.append(f)

    def hadd_files(self, files_in:list, file_out:Union[str, Path], opt='') -> None:
        '''
        
        '''
        args = ['hadd', opt, f' {file_out}']
        for f in files_in : args.append(f"{f}")
        stdout = subprocess.run(args=args, check=True, stdout=subprocess.PIPE).stdout.decode('utf-8')
        print(f"hadd stdout: {stdout}")

    def is_rootfile_valid(self, file:Union[str, Path], treename:str="HV_Mon", n_events:int=int(1e9), **kwargs):
        valid = True
        try : 
            content = self.read_rootfile(file=file, treename=treename, n_events=n_events)
            nev = len(content['time'])
            self.nevt_valid += nev
            self.dict_nevt[Path(file).name] = nev
        except ValueError as e:
            print(f"Root file {file} not valid")
            valid = False
        return valid
    
    def plot_timeserie(self, time:np.ndarray, dict_var:dict, file_out:Union[str, Path], step:int=1, title="Monitoring", **kwargs ):
        
        nvar = len(dict_var)
        nrow, ncol = nvar, 1 
        fig = plt.figure(figsize=(16,11), )
        left, right, top, bottom = 0.2/(ncol+1), 0.95-0.1/(ncol+1), 0.94-0.1/(nrow+1), 0.4/(nrow+1)
        kwargs_size = dict(wspace=0.0, hspace=0.1, 
                    top=top, bottom=bottom, 
                    left=left, right=right)
        gs = fig.add_gridspec(nrow, ncol, **kwargs_size )

        for ivar, (key, var) in enumerate(dict_var.items()):
            ax = fig.add_subplot(gs[ivar, 0])
            ax.plot(time[::step], var[::step], **kwargs)
            ax.set_ylabel(f'{key}')
            if ivar != nvar-1 : ax.set_xticklabels([])
            else : 
                datetime_ticks = [datetime.fromtimestamp(int(ts)).strftime('%d/%m/%Y %H:%M') for ts in ax.get_xticks()]
                ax.set_xticklabels(datetime_ticks)
                ax.set_xlabel('time')   
                for label in ax.get_xticklabels(which='major'):
                    label.set(rotation=30, horizontalalignment='right')

        plt.gcf().text(x=left, y=0.985-0.1/(nrow+1), s=title,  
                    fontsize='xx-large', 
                    rotation='horizontal', 
                    color = 'red',
                    fontweight='bold')  
        
        fig.savefig(file_out, dpi=300)
        print(f'Save figure {file_out}')
        plt.close()
    

    def main_hadd(self, files_in:list=None, file_out:Union[str,Path]=None):
        
        if files_in is None:
            files_in = self.files_read

        try : 
            print('\nhadd - 1st try\n')
            self.hadd_files(files_in=files_in, file_out=file_out)

        except subprocess.CalledProcessError as e:  
            print('\nhadd - 2nd try\n')
            new_files_in = []
            for f in files_in: 
                fout = f[:-5] + '_new' + '.root'
                print(f"\n{f} \n\t--> {fout}\n")
                try : 
                    self.hadd_files(files_in=[f,], file_out=fout, opt='-f',)
                    valid = self.is_rootfile_valid(fout)
                    if valid == True: 
                        new_files_in.append(fout)
                    else: 
                        # print(f"\nNot valid : {f}\n")
                        self.files_corrupt.append(fout)
                except: 
                    # print(f"\nError during hadd of {f}\n")
                    self.files_corrupt.append(f)
                    continue
            
            # print("new_files_in = ", new_files_in)
            self.hadd_files(files_in=new_files_in, file_out=file_out, opt='-f',)
            files_in = new_files_in

        self.to_df(files_in)
        self.df['name'] = [Path(f).name for f in files_in]
        self.df.set_index('name', inplace=True)
        self.df['nevt'] = self.df.index.map(self.dict_nevt)
    
        print(f"files_corrupt ({len(self.files_corrupt)}) = {self.files_corrupt}")
        print(f"nevt_valid = {self.nevt_valid} ")
            

    def main_timeseries(self, file_in:Union[str, Path], path_out:Union[str, Path], t_range:None=tuple, treename:str="HV_Mon", n_events:int=int(1e9), prefix:str=None, **kwargs):

        content = self.read_rootfile(file_in, treename, n_events, **kwargs)
        t = content['time']
        
        order = np.argsort(t) 
        ran = np.full(shape=t.shape, fill_value=True) if t_range is None else (min(t_range) < t[order]) & (t[order] < max(t_range))  
        nnull = t[order] != 0 
        for k, v in content.items():
            content[k] = v[order][ran & nnull]
        
        flow = content['flow']
        dict_monflow = {'flowIn [L/h]':flow[:,0,0], 'flowOut [L/h]': flow[:,0,1]}
        foutname = f"mosaic_monitoring_flow.png"
        fout = path_out / foutname

        self.plot_timeserie(time=content['time'], dict_var=dict_monflow, file_out=fout, title=f"{prefix} Monitoring gas flow")

        HVPS_5V  = content['HVPS_5V']
        HVPS_15V  = content['HVPS_15V']
        HVPS_3V3  = content['HVPS_3V3']
        Vmon  = content['V'] 
        Iset = content['I']
        feedback_type = content['feedback_type']
        # print("voltage var : ", HVPS_5V.shape, HVPS_15V.shape, HVPS_3V3.shape, Vmon.shape, Iset.shape, feedback_type.shape)
        dict_montens = {'HVPS_5V [V]':HVPS_5V, 'HVPS_15V [V]': HVPS_15V, 'HVPS_3V3 [V]' : HVPS_3V3, 'feedback_type': feedback_type}
        foutname = f"mosaic_Vmon.png"
        fout = path_out / foutname
        self.plot_timeserie(time=content['time'], dict_var=dict_montens, file_out=fout, title=f"{prefix} Monitoring HV")

        T_PC = content['T_PC']
        T_HV = content['T_HV']
        dict_montemp = {'T_PC [°C]':T_PC, 'T_HV [°C]': T_HV}
        foutname = f"mosaic_monitoring_temp.png"
        fout = path_out / foutname
        self.plot_timeserie(time=content['time'], dict_var=dict_montemp, file_out=fout, title=f"{prefix} Monitoring Temp")





if __name__ == "__main__":

    t0 = time.time()
    main_path = Path.home()
    CURRENT_SURVEY = set_survey("izen")
    tel = CURRENT_SURVEY.detector
    tel_name = tel.name.lower()
    data_path = main_path / "Projects" / "tomomu" / "Data" / CURRENT_SURVEY.name.lower()  #e.g INB72_monitoring_20240122_12H21_HV_mon.root


    prefix = CURRENT_SURVEY.run_prefix

    mon_path = data_path.parents[1] / tel_name / "mnt"
    mon_path = data_path / "3_monitoring"
    
    mon = Monitoring()
    print("mon_path = ", mon_path)
    mon.fetch_files(mon_path)
    mon.to_df( file = mon_path / "df_mon.csv")
    files_in = [str(mon_path / n) for n in mon.df['name'].to_list()]
    file_hadd_out = str(mon_path / "hadd_attempt.root")     
    lnevt = []
    ##test
    f = files_in[3]
    mon.read_rootfile(f)
    ###
    valid = mon.is_rootfile_valid(f)
    print(valid)
    mon.main_hadd(files_in=[files_in[3],], file_out=file_hadd_out)
    mon.main_timeseries(file=file_hadd_out)

    exit()

