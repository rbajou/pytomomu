#!/usr/bin/python3

from pathlib import Path
import sys
import yaml
import json
#package module(s)
from __init__ import set_config, SURVEY_NAME

if __name__ == "__main__":
    
    survey = sys.argv[1] if len(sys.argv) > 1  else 'inb72'
    
    with open("current_survey.json", 'w+') as f:
        json.dump({"current_survey": survey}, f, indent=4)

    if survey != SURVEY_NAME: 
        print("\nNEW SURVEY CONFIG:\n")
        set_config(survey)
        print()
    


