#!/usr/bin/python3
# -*- coding: utf-8 -*-

from dataclasses import dataclass, field
import numpy as np
from pathlib import Path
from typing import List, Union
import glob
#package module(s)
from detector import Telescope, Bench
from .run import Run
from .__init__ import SURVEY_NAME, DATA_PATH, DETECTOR, RUN_PREFIX, RUN_NUM, CLUSVAR_PREFIX

@dataclass
class Survey: 
   
    name : str
    data_path : Union[Path, str]
    detector : Union[Telescope, Bench]
    run_prefix : str
    run_num : list
    clusvar_prefix : str

    def __post_init__(self):
        self.runs = {}

    def __setitem__(self, name:str, run:Run): 
        self.runs[name] = run

    def __getitem__(self,name:str): 
        run = self.runs[name]
        return run

    def __str__(self): 
        sout = f"\nSurvey: {self.name}\n\n - "+ f"\n - ".join(v.__str__() for _,v in self.runs.items())
        return sout
    
CURRENT_SURVEY = Survey(name=SURVEY_NAME, 
                        data_path=DATA_PATH, 
                        detector=DETECTOR,
                        run_prefix=RUN_PREFIX,
                        run_num=RUN_NUM,
                        clusvar_prefix=CLUSVAR_PREFIX)

if __name__ == "__main__": 
    print(CURRENT_SURVEY)

