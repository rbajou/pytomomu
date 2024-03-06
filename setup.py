#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import glob
from pathlib import Path
import sys
import subprocess
import json


survey = "inb72"
opt = "--survey"
if opt in sys.argv:
    # Find the index of the custom option
    index = sys.argv.index(opt)
    # Check if the option has a value next to it
    if index + 1 < len(sys.argv):
        survey = sys.argv[index + 1]
        ###Load default script arguments stored in .yaml file
        # Remove the option and its value from sys.argv
        sys.argv.pop(index) #one time for '--opt'
        sys.argv.pop(index) #two time for 'value'

print(sys.argv)

REQUIREMENTS = [
    'pandas',
    'numpy',
    'scipy',
    'matplotlib',
    'scikit-learn',
    'scikit-image',
    'pyjson',
    'pylandau',  #used in 'utils/functions.py'
    'iminuit', #used in 'utils/functions.py'
    'pyyaml',
    'palettable', 
    'requests', 
    'bs4', 
    'cx_Freeze', #make executable files
    'mat73', #read v7.3 mat files
    'pillow', #for animation
    'uproot', #to read .root file
]


setup(
    name='pytomomu',
    version='0.1.0',
    description="",
    author="Raphaël Bajou",
    author_email='r.bajou2@gmail.com',
    url='https://github.com/rbajou/pytomomu.git',
    packages=find_packages(),
    package_dir={
        'pytomomu': [
                     'analyse', 'cluster', 'config', 'detector', 'event', 'real', 'simu', 'utils', 'tracking',
                     ]
    },  
    package_data={
          'files': [],
      },
    include_package_data=True,
    install_requires=REQUIREMENTS,
    keywords='Muography Nuclear Package',
    classifiers=[
        'Programming Language :: Python :: 3.10',
    ], 
)
basename = "current_survey"
print(f"\n{basename.upper()}: {survey}")
#stdout = subprocess.run(["python3", "config/config.py", f"{survey}" ], stdout=subprocess.PIPE).stdout.decode('utf-8')
#print(stdout)
with open(f"{basename}.json", 'w') as f:
    json.dump({basename: survey}, f, indent=4)
