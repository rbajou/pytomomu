from pathlib import Path
import sys
import yaml
import json
import glob
import re

#package module(s)
# from .config import set_config
from detector import DICT_DET, str2detector

with open("current_survey.json", 'r') as f: 
    metadata_survey = json.load(f)
survey = metadata_survey['current_survey']

def set_config(name:str):

    global SURVEY_NAME
    global DATA_PATH
    global DETECTOR
    global RUN_PREFIX
    global RUN_NUM
    global CLUSVAR_PREFIX 

    SURVEY_NAME = name

    with open( "survey.yaml" ) as f:
        try:
            survey_config = yaml.load(f, Loader=yaml.SafeLoader) # The FullLoader parameter handles the conversion from YAML scalar values to Python the dictionary format
        except yaml.YAMLError as exc:
            print(exc) 
    
    #Check survey name
    try :
        assert name in survey_config.keys()
    except AssertionError:
        print( f"'{name}' survey was not found (available: {list(survey_config.keys())}). If needed, edit 'survey.yaml' file.")

    dict_args_survey = survey_config[name]
    #Check fields
    lkeys = ["detector", "data_path", "run_num", "run_prefix", "clusvar_prefix"]
    for k in lkeys : assert k in list(dict_args_survey.keys()), f"No field '{k}' found, check 'survey.yaml' file."
        
    DET_NAME = dict_args_survey["detector"]
    print(f"DETECTOR = {DET_NAME}")
    #Check detector
    try : 
        assert DET_NAME in DICT_DET.keys()
    except AssertionError: 
        print(f"'{DET_NAME} is not defined (available: {list(DICT_DET.keys())}). If needed, edit 'detector.py' file.")
    DETECTOR = str2detector(DET_NAME)


    DATA_PATH = dict_args_survey["data_path"]
    print(f"DATA_PATH = {DATA_PATH}")
    try : 
        assert Path(DATA_PATH).exists()
    except AssertionError: 
        print(f"'{DATA_PATH}' does not exist.")

    l_avail_run = glob.glob(str(Path(DATA_PATH)/"run*"))
    l_avail_run_num = sorted([int(re.findall(r'\d+', r.split('/')[-1])[0]) for r in l_avail_run])
    input_run_num = dict_args_survey['run_num']
    print('INPUT RUN_NUM = ', input_run_num)
    RUN_NUM = []
    if input_run_num is None: 
        try : 
            RUN_NUM.append(l_avail_run_num[0])
        except IndexError:
            print(f"No run found in {DATA_PATH}")
    else : 
        for irun in RUN_NUM:
            try : 
                assert irun in l_avail_run_num
                RUN_NUM.append(irun)
            except AssertionError: 
                print(f"run{irun} does not exist, will not be processed.")
        print(f"Available run numbers are: {l_avail_run_num}")
    print('PROCESS RUN_NUM (default) = ', RUN_NUM)

    RUN_PREFIX = dict_args_survey["run_prefix"]
    print(f"RUN_PREFIX = {RUN_PREFIX}")
    CLUSVAR_PREFIX = dict_args_survey["clusvar_prefix"]
    print(f"CLUSVAR_PREFIX = {CLUSVAR_PREFIX}")


set_config(survey)