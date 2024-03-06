#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from abc import abstractmethod

class Hit:
    
    @abstractmethod
    def hit(self):
        pass

class HitCollection: 

    @property
    def nhit(self):
        """Number of hits"""
        return self._nhit
    
    @nhit.setter
    def nhit(self, value:int):
        self._nhit = value

    @property
    def pdgcode(self):
        """Pdg code (13 for muons, 11 for electrons)"""
        return self._pdgcode
    
    @pdgcode.setter
    def pdgcode(self, value:int):
        self._pdgcode = value

    @property
    def det_ix(self):
        """Index detector panel"""
        return self._det_ix
    
    @det_ix.setter
    def det_ix(self, value:int):
        self._det_ix = value


    @property
    def position(self):
        """Hit positions [mm]"""
        return self._position
    
    @position.setter
    def position(self, value:float):
        self._position = value

    @property
    def edep(self):
        """Energy deposits [GeV]"""
        return self._edep
    
    @edep.setter
    def edep(self, value:float):
        self._edep = value
