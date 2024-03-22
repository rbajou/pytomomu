#!/usr/bin/python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import glob
from pathlib import Path
import sys
import subprocess
import json


# survey = "inb72"
# opt = "--survey"
# if opt in sys.argv:
#     # Find the index of the custom option
#     index = sys.argv.index(opt)
#     # Check if the option has a value next to it
#     if index + 1 < len(sys.argv):
#         survey = sys.argv[index + 1]
#         ###Load default script arguments stored in .yaml file
#         # Remove the option and its value from sys.argv
#         sys.argv.pop(index) #first time for '--opt'
#         sys.argv.pop(index) #second time for 'value'

REQUIREMENTS = [
    # 'bs4', #to fetch online files
    # 'cx_Freeze', #to make executable files
    'iminuit', #used for advanced fit in 'utils/functions.py'
    'matplotlib',
    'numpy',
    'palettable', # for nice color maps
    'pandas',
    'pillow', #for animation
    'pyjson',
    'pylandau',  #used in 'utils/functions.py'
    'pyyaml', #for .yaml config files
    'scikit-learn', #ml library
    'scikit-image', #for fit with ransac  
    'scipy',
    'seaborn',#for nice plot templates
    'uproot', #to read .root file
]


setup(
    name='pytomomu',
    version='0.1.0',
    description="This package is for processing and analyzing muography data recorded in Micromegas detectors developed at CEA-IRFU, Saclay, France.",
    author="Raphaël Bajou",
    author_email='raphael.bajou@cea.fr',
    url='https://github.com/rbajou/pytomomu.git',
    packages=find_packages(),
    package_dir={
        'pytomomu': [
                     'analyse', 'cluster', 'detector', 'event', 'forwardsolver', 'real', 'rootfile', 'simu', 'survey', 'utils', 'tracking',
                     ]
    },  
    include_package_data=True,
    package_data={
          #'simu':  ['*.json', '*.yaml'],
          'survey':  [ '*.yaml'],
    },
    install_requires=REQUIREMENTS,
    keywords='Muography, Micromegas, Tracking',
    classifiers=[
        'Programming Language :: Python :: 3.10',
    ], 
)

# basename = "current_survey"
# path = Path(__file__).parent / "survey" / basename
# print(f"\n{basename.upper()}: {survey}")
# #subprocess.run(["python3", "config/config.py", f"{survey}" ], stdout=subprocess.PIPE).stdout.decode('utf-8')
# with open(f"{path}.json", 'w') as f:
#     json.dump({basename: survey}, f, indent=4)
