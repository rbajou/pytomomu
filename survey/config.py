#!/usr/bin/python3

from pathlib import Path
import sys
import yaml
import json
#package module(s)
from survey import set_survey


if __name__ == "__main__":

    survey_in = sys.argv[1] if len(sys.argv) > 1  else 'inb72'

    with open( Path(__file__).parent / "current_survey.json", 'r') as f:
        current_survey = json.load(f)
   
    if survey_in != current_survey: 
        print("\nNEW SURVEY CONFIG:\n")
        set_survey(survey_in)

        with open( "survey/current_survey.json", 'w+') as f:
                json.dump({"current_survey": survey_in}, f, indent=4)