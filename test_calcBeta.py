# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:17:18 2019

@author: adming
"""
import unittest
import numpy as np
from pbm_MagnitudeAnalysis import setMatPlotLib, setParams, drawBeta, calcBeta

class Test_calcBeta(unittest.TestCase):
    def test_init(self):
        global params
        params = setParams(unittestcode=0)
        print params # 
        self.assertEqual(params['flu_T'], 328.)
#        self.assertEqual(params['tur_eps'], 0.59)
        self.assertEqual(params['par_rho'], 2250.)
        self.assertEqual(params['rel_vel'], 0.)
#        tur_eps = 1000. # turbulent flow epsilon
#        kB = 1.381e-23
#        flu_T = 140 + 273 # K
#        flu_niu = 2.5e-5 # fluid niu
#        flu_rho = 0.88
#        flu_P = 101325.
#        flu_miu = flu_niu * flu_rho
#        flu_moleculediam = 364. * 1e-12 # using kinetic diameter, https://en.wikipedia.org/wiki/Kinetic_diameter
#        flu_lamda = kB * flu_T / (np.sqrt(2.) * np.pi * flu_moleculediam**2 * flu_P) # https://en.wikipedia.org/wiki/Mean_free_path
#        par_rho = 2250. # kg/m3
        
    def test_beta(self):
        tar_dp_ = 0.0015e-6
#        dp_logRng = np.linspace(np.log10(tar_dp_ * 1.0e6), 2, num=int(4/0.02))
        dp_logRng = np.linspace(np.log10(0.0015), np.log10(100.), num=300)
        dp = np.power(np.ones_like(dp_logRng)*10, dp_logRng) * 1e-6 # in m
        dp_tar = dp[np.argmax(dp>=tar_dp_)] # dp[0]
        beta_B, beta_TS, beta_TI, beta_IC, beta_IM, beta_Total = calcBeta(params, dp=dp, dp_tars=[dp_tar], \
                                          dp_icAndim_ignore={'dp_tar_upbound':2.5e-6,'dp_src_lowbound':10e-6}, ldraw=True)
        print beta_B[0][0], beta_B[0][-1]
        print beta_TS[0][0] + beta_TI[0][0], beta_TS[0][-1] + beta_TI[0][-1]
#        print flu_T
#        print beta_B
        
if __name__ == '__main__':
    unittest.main()        