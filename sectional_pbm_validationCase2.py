# -*- coding: utf-8 -*-
"""
Created on Mon Jul 02 10:56:07 2018

@author: gliu

reference: Direct simulation Monte Carlo method for particle coagulation and aggregation.Kruis 2000.pdf
comparison with the Case 2 in that paper, i.e. beta(i,j) = A*(i+j)
"""
import numpy as npy
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
matplotlib.rcParams['font.family'] = 'Times New Roman'

#%%
# keep in mind that v is volume (not the diameter and of course NOT the velocity) of the particle (aggregate).
v_min = 1.0
v_max = 100.0
v_ratio = 1.08
Nv = int(npy.log10(v_max/v_min) / npy.log10(v_ratio))
ns = npy.arange(Nv+3)
Nv = len(ns) -1
vs = v_ratio ** ns
vs_tmp = (vs[:-1] + vs[1:])/2
vs_h = (vs[:-1] + vs[1:])/2
vs_l = vs_h[:-1].copy()
vs_l = npy.concatenate(([vs[0]],vs_l))
vs = vs[:-1]
vs_bounds = npy.column_stack((vs, vs_l, vs_h))
dps = npy.power((vs * 6 / npy.pi), 1./3.)
#%%
N0 = 1.e+10
A = 1./N0
dt = 0.001
max_iters = 1000 # 必须足够长，长到计算终点时刻的v_max尺寸的颗粒团数目仍然非常少，否则tot_count的统计值会不准，因为超出了粒径范围
#
#end_t = dt * max_iters
#P_v_max_estimate = (N0*A*end_t/2)**(v_max/v_min-1) / (1 + N0*A*end_t/2)**(v_max/v_min)
#if P_v_max_estimate >= 0.001:
#    print 'Warning, end_t is too long, please make max_iters smaller or make v_max larger!'
#    exit(-1)
#
par_counts = npy.zeros_like(vs)
par_counts[0] = N0
tot_counts = npy.zeros((max_iters+1,1))
tot_counts[0] = npy.sum(par_counts)
parN_t = npy.zeros((max_iters+1, Nv+1)) # first column is time
parN_t[0, :] = npy.insert(par_counts, 0, 0.)
#
progress_report_times = 20
iter_per_report = max_iters / progress_report_times
i_report = 0
i_reportlast = 0
#
for i_t in range(max_iters):
    t = (i_t+1) * dt
    dndt = npy.zeros_like(vs)
    for i in range(Nv):        
        ni = par_counts[i]
        vi = vs[i]
        for j in range(i, Nv):
            nj = par_counts[j]
            vj = vs[j]
            beta_ij = A * (vi/v_min + vj/v_min) #((i+1) + (j+1))
            vij = vi + vj
            k = npy.argmax(vs_h > vij)
            if i != j:
                dndt[k] = dndt[k] + vij / vs[k] * ni * nj * beta_ij
                dndt[i] = dndt[i] - ni * nj * beta_ij
                dndt[j] = dndt[j] - ni * nj * beta_ij
            else:
                dndt[k] = dndt[k] + vij / vs[k] * (0.5 * ni * nj * beta_ij)
                dndt[i] = dndt[i] - (0.5 * ni * nj * beta_ij) * 2
    d_counts = dndt * dt
    par_counts = par_counts + d_counts
    parN_t[i_t+1, :] = npy.insert(par_counts, 0, t) 
    tot_counts[i_t+1] = npy.sum(par_counts)            
    # report (print) the progress of calculation
    i_report = (i_t + 1) / iter_per_report
    if i_report > i_reportlast:
        print '   {}% finished'.format(i_report*100./progress_report_times)
        i_reportlast = i_report
#            
primary_par_counts = vs / v_min
count_percent = par_counts / tot_counts[-1]
#
NtoN0ratio_calc = tot_counts[-1]/tot_counts[0]  
NtoN0ratio_expect = npy.exp(-N0*A*(dt*max_iters))     
print '   N/N0 at end of simulation given by PBM: ', NtoN0ratio_calc
print '   Expected (theoretical) value: ', NtoN0ratio_expect 
pass_relative_threshld = 0.02
if abs(NtoN0ratio_calc/(NtoN0ratio_expect + 1.0e-13) - 1.) > pass_relative_threshld:
    print '[Warning] Relative error is too big!'
else:
    print '[Good] Relative error is acceptable'  