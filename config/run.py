#!/usr/bin/python3
# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
#!/usr/bin/env python3

import numpy as np
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Union
from enum import Enum, auto
import yaml
#package module(s)


@dataclass
class Run:
   
    name : str
    path : Path #List[RawData] = field(default_factory=lambda : list())
    
    def __str__(self): 
        sout = f"Run: {self.name}"
        return sout