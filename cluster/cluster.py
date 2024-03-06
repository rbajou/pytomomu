#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dataclasses import dataclass, field
import numpy as np


@dataclass
class Strip:
    index : int
    amplitude : float
    _position : float = field(init=False)

    @property
    def position(self):
        return self._position
    
    @position.setter
    def position(self, interstep:float):
        self._position = self.index * interstep



class StripCollection:

    @property
    def amplitude(self):
        return self._amplitude
    
    @amplitude.setter
    def amplitude(self, value:np.ndarray):
        """
        value (np.ndarray) : shape (nchip, nchannel, nsample) 
        """
        self._amplitude = value
        self._nchip, self._nchannel, self._nsample = value.shape

    @property
    def nchip(self):
        return self._nchip
    
    @property
    def nchannel(self):
        return self._nchannel
    
    @property
    def nsample(self):
        return self._nsample


# class Cluster:
    
#     @abstractmethod
#     def cluster(self, model):
#         pass

class ClusterCollection:


    @property
    def nclus(self):
        """Number of clusters"""
        return self._nclus
    
    @nclus.setter
    def nclus(self, value:int):
        self._nclus = value

    @property
    def position(self):
        """Barycenter cluster strips position strips"""
        return self._position
    
    @position.setter
    def position(self, value:float):
        self._position = value


    @property
    def amplitude(self):
        """Sum amplitude of cluster strips [ADC]"""
        return self._max_ampl
    
    @amplitude.setter
    def amplitude(self, value:float):
        self._amplitude = value

    @property
    def max_amplitude(self):
        """Maximal amplitude in cluster strips [ADC]"""
        return self._max_ampl
    
    @max_amplitude.setter
    def max_amplitude(self, value:float):
        self._max_ampl = value

    @property
    def max_amplitude_strip(self):
        """Strip with maximal amplitude in cluster strips [ADC]"""
        return self._max_ampl_strip
    
    @max_amplitude_strip.setter
    def max_amplitude_strip(self, value:float):
        self._max_ampl_strip = value
  
    @property
    def size(self):
        """Size cluster [nstrips]"""
        return self._size
    
    @size.setter
    def size(self, value:float):
        self._size = value

    @property
    def tot(self):
        """Time-over-threshold [nbins]"""
        return self._tot
    
    @tot.setter
    def tot(self, value:int):
        self._tot = value

    @property
    def max_sample(self):
        """Maximum samples [nbins]"""
        return self._max_sample
    
    @max_sample.setter
    def max_sample (self, value:int):
        self._max_sample = value

    @property
    def time(self):
        """Time-over-threshold [float]"""
        return self._time
    
    @time.setter
    def time(self, value:float):
        self._time = value

