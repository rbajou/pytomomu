#!/usr/bin/python3
# -*- coding: utf-8 -*-

from dataclasses import dataclass, field
import numpy as np
from pathlib import Path
from typing import List, Union
import glob
import yaml
import json
import glob
import re
#package module(s)
from detector import Telescope, Bench, DICT_DET, str2detector
from .run import Run

@dataclass
class Survey: 
   
    name : str
    data_path : Union[Path, str]
    detector : Union[Telescope, Bench]

    def __post_init__(self):
        self.runs = {}

    def __setitem__(self, name:str, run:Run): 
        self.runs[name] = run

    def __getitem__(self,name:str): 
        run = self.runs[name]
        return run

    def __str__(self): 
        sout = f"Survey: {self.name}\n\n - "+ f"\n - ".join(v.__str__() for _,v in self.runs.items())
        return sout



class RealSurvey(Survey):
    
    def __init__(self, name:str, data_path:Union[Path, str], detector:Union[Telescope, Bench], run_prefix, run_num:list, clusvar_prefix:str):
        Survey.__init__(self, name, data_path, detector)
        self.run_prefix = run_prefix 
        self.run_num = run_num
        self.clusvar_prefix = clusvar_prefix

    def __str__(self): 
        sout = f"Real Survey: {self.name}\n\n - "+ f"\n - ".join(v.__str__() for _,v in self.runs.items())
        return sout



class SimuSurvey(Survey):

    def __init__(self, name:str, data_path:Union[Path, str], detector:Union[Telescope, Bench]):
        Survey.__init__(self, name, data_path, detector)

    def get_metadata(self, file:Union[Path, str] = Path(__file__).parents[1] / "simu" / "metadata_simu.json"):
        try :
            assert file.exists()
        except AssertionError:
            exit(f"{file} does not exist")

        self.metadata = None 
        try :
            with open(file) as f: 
                self.metadata = json.load(f)
        except : 
            exit(f"Issue when reading {file}")


    def __str__(self): 
        sout = f"Simu Survey: {self.name}\n\n - "+ f"\n - ".join(v.__str__() for _,v in self.runs.items())
        return sout


def set_real_survey(name:str) -> RealSurvey:

    with open( Path(__file__).parent / "survey.yaml") as f: #same
            try:
                survey_config = yaml.load(f, Loader=yaml.SafeLoader) # The FullLoader parameter handles the conversion from YAML scalar values to Python the dictionary format
            except yaml.YAMLError as exc:
                print(exc) 

     #Check survey name
    try :
        assert name in survey_config.keys()
    except AssertionError:
        print( f"'{name}' survey was not found (available: {list(survey_config.keys())}). If needed, edit 'survey.yaml' file.")

    DET_NAME = survey_config[name]["detector"]
    print(f"DETECTOR = {DET_NAME}")
    #Check detector
    try : 
        assert DET_NAME in DICT_DET.keys()
    except AssertionError: 
        print(f"'{DET_NAME} is not defined (available: {list(DICT_DET.keys())}). If needed, edit 'detector.py' file.")
    DETECTOR = str2detector(DET_NAME)

    dict_real_survey = survey_config[name]['real']
    #Check fields
    lkeys = ["data_path", "run_num", "run_prefix", "clusvar_prefix"]
    for k in lkeys : assert k in list(dict_real_survey.keys()), f"No field '{k}' found, check 'survey.yaml' file."
        

    DATA_PATH = dict_real_survey["data_path"]
    print(f"DATA_PATH = {DATA_PATH}")
    try : 
        assert Path(DATA_PATH).exists()
    except AssertionError: 
        print(f"'{DATA_PATH}' does not exist.")

    l_avail_run = glob.glob(str(Path(DATA_PATH)/"run*"))
    l_avail_run_num = sorted([int(re.findall(r'\d+', r.split('/')[-1])[0]) for r in l_avail_run])
    input_run_num = dict_real_survey['run_num']
    RUN_NUM = []
    if input_run_num is None: 
        print('PROCESS RUN_NUM (default) = ', RUN_NUM)    
        try : 
            RUN_NUM.append(l_avail_run_num[0])
        except IndexError:
            print(f"No run found in {DATA_PATH}")
    
    else : 
        print('INPUT RUN_NUM = ', input_run_num)
        for irun in RUN_NUM:
            try : 
                assert irun in l_avail_run_num
                RUN_NUM.append(irun)
            except AssertionError: 
                print(f"run{irun} does not exist, will not be processed.")
                print(f"Available run numbers are: {l_avail_run_num}")

    RUN_PREFIX = dict_real_survey["run_prefix"]
    print(f"RUN_PREFIX = {RUN_PREFIX}")
    CLUSVAR_PREFIX = dict_real_survey["clusvar_prefix"]
    print(f"CLUSVAR_PREFIX = {CLUSVAR_PREFIX}")

    
    CURRENT_SURVEY = RealSurvey(name=name, 
                        data_path=DATA_PATH, 
                        detector=DETECTOR,
                        run_prefix=RUN_PREFIX,
                        run_num=RUN_NUM,
                        clusvar_prefix=CLUSVAR_PREFIX)
    
    return CURRENT_SURVEY



def set_simu_survey(name:str) -> RealSurvey:

    with open( Path(__file__).parent / "survey.yaml") as f: #same
            try:
                survey_config = yaml.load(f, Loader=yaml.SafeLoader) # The FullLoader parameter handles the conversion from YAML scalar values to Python the dictionary format
            except yaml.YAMLError as exc:
                print(exc) 

     #Check survey name
    try :
        assert name in survey_config.keys()
    except AssertionError:
        print( f"'{name}' survey was not found (available: {list(survey_config.keys())}). If needed, edit 'survey.yaml' file.")

    DET_NAME = survey_config[name]["detector"]
    print(f"DETECTOR = {DET_NAME}")
    #Check detector
    try : 
        assert DET_NAME in DICT_DET.keys()
    except AssertionError: 
        print(f"'{DET_NAME} is not defined (available: {list(DICT_DET.keys())}). If needed, edit 'detector.py' file.")
    DETECTOR = str2detector(DET_NAME)

    dict_simu_survey = survey_config[name]['simu']
    #Check fields
    lkeys = ["data_path", "name"]
    for k in lkeys : assert k in list(dict_simu_survey.keys()), f"No field '{k}' found, check 'survey.yaml' file."
        
    DATA_PATH = dict_simu_survey["data_path"]
    print(f"DATA_PATH = {DATA_PATH}")
    try : 
        assert Path(DATA_PATH).exists()
    except AssertionError: 
        print(f"'{DATA_PATH}' does not exist.")


    CURRENT_SURVEY = SimuSurvey(name=name, 
                        data_path=DATA_PATH, 
                        detector=DETECTOR)
    
    return CURRENT_SURVEY


def set_survey(name:str) :

    set_real_survey(name)
    set_simu_survey(name)
