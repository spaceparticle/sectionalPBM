# -*- coding: utf-8 -*-
"""
Created on Mon Jul 02 10:56:07 2018

@author: gliu

reference: Direct simulation Monte Carlo method for particle coagulation and aggregation.Kruis 2000.pdf
"""
import numpy as npy
import matplotlib.pyplot as plt
#from matplotlib.pyplot import *
from pbm_MagnitudeAnalysis import setMatPlotLib, setParams, calcBeta
import copy
from scipy import interpolate
import pandas as pd
from pandas import DataFrame, Series
#matplotlib.rcParams['font.family'] = 'Times New Roman'
import time
from datetime import datetime

figsize_dflt, total_markers = setMatPlotLib()
params = setParams()

def interpOnParRelVel(relvel_src, dp_tar, ldbg=False):
    f_tar = interpolate.interp1d(relvel_src.index.values, relvel_src.values, kind='slinear')
    if ldbg: print type(f_tar(dp_tar[0])), f_tar(dp_tar[0]).shape
    relvel_tar = [float(f_tar(dp_)) for dp_ in dp_tar]
    if ldbg: print relvel_tar
    relvel_tar = Series(relvel_tar, index=dp_tar) 
    return relvel_tar

def readParRelVel(filename='2D_xt4D_Standard.csv', ldbg=False):
    if ldbg: print 'in <readParRelVel> '
    if not filename:
        return
    df = pd.read_csv(filename)
    cols = df.columns
    dps = []
    relvel_mean = []
    for col in cols:
        col_split = col.split('um_rel_vel')
        if len(col_split) == 2:
            dp = float(col_split[0]) 
            dps.append(dp)
            relvel_mean.append(df[col].mean())
    if ldbg: print 'dps: ', dps
    relvel_series = Series(relvel_mean, index=dps) 
    return relvel_series

def drawPSD(dps, np, figttl=''):
    plt.figure(figsize=figsize_dflt)
    plt.plot(dps, np, 'b-', markersize = 5)
#    plt.xticks(fontsize=tick_size) #, rotation=90) # plot.tick_params(axis='both', which='major', labelsize=10)
#    plt.yticks(fontsize=tick_size)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('particle diamter, micron')
#    plt.ylabel('particle number density(?)',fontsize=label_size)
#    legend(['Brownian', 'Turbulent shearing', 'Turbulent inertia', 'Total'], fontsize=tick_size)
#    foo = 1e7
#    plt.title(figttl, fontsize=label_size)
    
#%%
def getnp(dps=[1.0], linmiu=True):
    if not linmiu:
        dps = dps * 1e6 # convert to in miu m
    # 粒径分布采用赵海波公式6.38设定（红皮专著p244页）
    Np = npy.array([5e14, 1e11, 1e9])
    dpg = npy.array([0.08, 2.0, 10.0]) # miu m
    sgmpg = npy.array([1.5, 2.0, 1.5])
    ln_sgmpg = npy.log(sgmpg)
    np = [npy.sum(Np / npy.sqrt(2*npy.pi)/ln_sgmpg * npy.exp(-(npy.log(dp/dpg))**2/2/ln_sgmpg**2) \
                  / dp) for dp in dps]
#    delta_dp = dps[1:] - dps[:-1]
    delta_dp = npy.diff(dps)
    np_bin = delta_dp * np[:-1]
    return np, np_bin

def getInitialPSD(dp=None):
    if dp is None:
        #dp=10**(npy.arange(-1, 2, 0.05))
        dp = 10**npy.linspace(npy.log10(0.08), npy.log10(20), 200) # miu m，粒径分布采用赵海波公式6.38设定（红皮专著p244页）
    np, np_bin = getnp(dps=dp)
    #drawPSD(dp, np, 'np ~ dp')
    log_dp = npy.log10(dp*1e-6)
    par_m = npy.pi / 6 * (dp*1e-6)**3 * params['par_rho'] * 1e6 # mg
#    delta_dp=dp[1:]-dp[:-1]
#    np_bin = delta_dp*np[:-1]
    Np0 = npy.sum(np_bin)
    par_m_bin = par_m[:-1] * np_bin
    print 'Np0: ', Np0
    print 'total mass: ', npy.sum(par_m_bin)
    #
    #drawPSD(dp[:-1], np_bin/Np0, 'np_bin/Np0 ~ np')
    #delta_m = par_m[1:]- par_m[:-1]
    delta_log_dp = log_dp[1:]- log_dp[:-1]
    drawPSD(dp[:-1]*1e-6, par_m_bin/delta_log_dp, 'dM/dlog(dp) ~ dp')
    return

#%%

relvel_mean = readParRelVel(ldbg=True)
dptar = npy.array([10, 12, 15, 21], dtype=npy.float)
relvel = interpOnParRelVel(relvel_mean, dptar)
#time.sleep(5)    
#getInitialPSD()

laddOneBig = False # add One group of big particles etc 50um
if laddOneBig:
    mass_OB = 5e3 # mg
    dpOB = 35.e-6 # m
#
dp_norm0 = 1e-6 # in m
v_norm0 = npy.pi/6 * dp_norm0**3
#    
v_min = 0.08**3 #0.01 miu
v_max = 20.0**3 #100 miu
v_ratio = 1.08
Nv = int(npy.log10(v_max/v_min) / npy.log10(v_ratio))
ns = npy.arange(Nv+3)
Nv = len(ns) -1
vs = v_min * v_ratio ** ns
vs_tmp = (vs[:-1] + vs[1:])/2
vs_h = (vs[:-1] + vs[1:])/2
vs_l = vs_h[:-1].copy()
vs_l = npy.concatenate(([vs[0]],vs_l))
vs = vs[:-1]
vs_bounds = npy.column_stack((vs, vs_l, vs_h))
dps = npy.power((vs * v_norm0 * 6 / npy.pi), 1./3.)
#
beta_B, beta_TS, beta_TI, beta_IC, beta_IM, B = calcBeta(params, dp=dps, dp_tars=dps, \
             dp_icAndim_ignore={'dp_tar_upbound':1.67e-6,'dp_src_lowbound':10e-6}, ldraw=False) #calcBeta(dpin = dps, linmiu=False)
#
if npy.max(npy.abs(beta_B-beta_B.T)) > 1e-23 or npy.max(npy.abs(beta_TS-beta_TS.T)) > 1e-23 or \
    npy.max(npy.abs(beta_TI-beta_TI.T)) > 1e-23: # B must be symmetric (B.T must = B)
    print 'Potential ERROR, certain betas are not symmetric, pausing for couple of seconds, better STOP the program and check !'
    time.sleep(20)
log_dp = npy.log10(dps)
#
np_, Np_bin = getnp(dps=dps, linmiu=False)

print dps[-1], len(dps)
#
#if laddOneBig: npy.append(dps, dpOB) # BUG!
if laddOneBig: 
    dps = npy.append(dps, dpOB)
print dps[-1], len(dps), len(Np_bin)
par_m = npy.pi / 6 * dps**3 * params['par_rho'] * 1e6 # mass of each single particle, mg
if laddOneBig:
    Np_bin = npy.append(Np_bin, mass_OB/par_m[-1])
    par_m_bin = par_m * Np_bin   
else:
    par_m_bin = par_m[:-1] * Np_bin
#    
print '>>total mass of particles, mg/Nm3 (not counting the mass of the "addOneBig"): ', \
npy.sum(par_m_bin[:-1]) if laddOneBig else npy.sum(par_m_bin)
mass_conc_inNm3 = npy.sum(par_m_bin)
conv_div_ = params['flu_T']/273.15 * 101325./params['flu_P'] * 1./(1 - 0.07) # H2O vfrac of flue gas at ESP is ~7%
print '...(raw Np): ', npy.sum(Np_bin)
print '...converting Np_bin and par_m_bin from units in Nm3 to real m3 with conv_div_: {}'.format(conv_div_)
Np_bin = Np_bin / conv_div_
par_m_bin = par_m_bin / conv_div_
mass_conc_inm3 = npy.sum(par_m_bin)
assert npy.abs(mass_conc_inm3 * conv_div_ / mass_conc_inNm3 - 1.) <= 1e-8, 'Error in converting'
#
delta_log_dp = log_dp[1:]- log_dp[:-1]
#drawPSD(dps[:-1], par_m_bin/delta_log_dp, 'dM/dlog(dp) ~ dp')
#%%
lEvolv = True
if lEvolv:
    print 'evolving ...'
    wt0 = datetime.now()
    #
    if not laddOneBig: Np_bin = npy.append(Np_bin, 0)
    dt = 0.001
    max_iters = 50 # 不能过长，使得计算终点时刻的v_max尺寸的颗粒团数目仍然非常少，否则tot_count的统计值会不准，因为超出了粒径范围
    #    
    progress_report_times = 50
    iter_per_report = max(max_iters / progress_report_times, 1)
    i_report = 0
    i_reportlast = 0
    #    
    # par_counts = npy.zeros_like(vs)
    par_counts = Np_bin
    tot_counts = npy.zeros((max_iters+1,1))
    tot_counts[0] = npy.sum(par_counts)
    parN_t = npy.zeros((max_iters+1, Nv+1)) # first column is time
    parN_t[0, :] = npy.insert(par_counts, 0, 0.)
    for i_t in range(max_iters):
        t = (i_t+1) * dt
        dndt = npy.zeros_like(vs)
        for i in range(Nv):        
            ni = par_counts[i]
            vi = vs[i]
            for j in range(i, Nv):
                nj = par_counts[j]
                vj = vs[j]
                beta_ij = B[i,j] #A * (vi/v_min + vj/v_min) #((i+1) + (j+1))
                vij = vi + vj
                k = npy.argmax(vs_h > vij)
#                if i == Nv -1 and j == Nv -1 :
#                    print 'k:',k
#                    print vs_h
#                    print 'vij:',vij
                if k == 0:
#                    print i_t, i,j
                    k = Nv-1
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
            wt1 = datetime.now()
            percent_finished = i_report*100./progress_report_times
            print '   {:.1f}% finished, ETA: {:.1f} seconds'.format(percent_finished, (wt1-wt0).total_seconds() * (100. - percent_finished)/percent_finished)
            i_reportlast = i_report
    ts = npy.arange(max_iters + 1) * dt
    primary_par_counts = vs / v_min
    count_percent = par_counts / tot_counts[-1]
    parM_t = copy.deepcopy(parN_t)
    parM_t[:,1:] = parN_t[:,1:] * par_m
    tot_mass = npy.sum(parM_t[:,1:], 1)
    yita = 100*(1 - parM_t[-1,1:-2]/parM_t[0,1:-2]) # percentage of removal
    rmvl = npy.column_stack((dps[:-2], yita))
    print 'N/N0 at end of simulation given by PBM: ', tot_counts[-1]/tot_counts[0]     
    print 'Ratio of second biggest particle counts between end and start of simulation: ', \
    parN_t[-1,-2]/parN_t[0,-2]
    #%%
    plt.figure(figsize=figsize_dflt)
    plt.plot(ts, tot_counts, 'b-', markersize = 5)
#    plt.xticks(fontsize=tick_size) #, rotation=90) # plot.tick_params(axis='both', which='major', labelsize=10)
#    plt.yticks(fontsize=tick_size)
#    plt.xscale('')
#    plt.yscale('log')
    plt.xlabel('time (s)')#,fontsize=label_size)
    plt.ylabel('particle total number')#, fontsize=label_size)
    plt.title('Evolution of particle total counts with time')#, fontsize=label_size)
    #%%
    plt.figure(figsize=figsize_dflt)
    plt.plot(dps[0:-1]*1e6, parM_t[0,1:-1]/npy.diff(npy.log10(dps)), 'b-', markersize = 5)
    t1_ = ts[-1] * 0.5
    it1_ = npy.argmax(ts>=t1_)
    plt.plot(dps[0:-1]*1e6, parM_t[it1_, 1:-1]/npy.diff(npy.log10(dps)), 'm-', markersize = 5)
    plt.plot(dps[0:-1]*1e6, parM_t[-1, 1:-1]/npy.diff(npy.log10(dps)), 'r-', markersize = 5)
    plt.legend(['t=0 s','t=' + str(ts[it1_]) + 's', 't=' + str(ts[-1]) + 's'])#, fontsize=tick_size)
#    plt.xticks(fontsize=tick_size) #, rotation=90) # plot.tick_params(axis='both', which='major', labelsize=10)
#    plt.yticks(fontsize=tick_size)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\it{d}_p$ ($\mu$m)')#,fontsize=label_size)
    plt.ylabel(r'$d\it{M}$/dlog($\it{d}_p$) (mg/m3)')#, fontsize=label_size)
    plt.title('dM/dlog(dp) ~ dp')#, fontsize=label_size)

    #%%
    plt.figure(figsize=figsize_dflt)
    plt.plot(dps[:-2]*1e6, yita, 'xr-', markersize = 5)
    plt.xscale('log')
    plt.xlabel(r'$\it{d}_p$ ($\mu$m)')#,fontsize=label_size)
    plt.ylabel(r'$\it{\eta}$ (%)')#, fontsize=label_size)
    plt.title("Percentage of 'removal' ~ dp")#, fontsize=label_size)
    #%% interpolate to certain dp points
    dp_tar = 10**(npy.arange(-1, 1, 0.05) - 6)
    f_tar = interpolate.interp1d(dps[:-2], yita, kind='slinear')
    yita_tar = f_tar(dp_tar)
    rmvl_tar = npy.column_stack((dp_tar, yita_tar))
    t2_ = ts[-1] * 0.5
    it2_ = npy.argmax(ts>=t2_)
    yita2 = 100*(1 - parM_t[it2_,1:-2]/parM_t[0,1:-2]) 
    f_tar2 = interpolate.interp1d(dps[:-2], yita2, kind='slinear')
    yita_tar2 = f_tar2(dp_tar)
    rmvl_tar2 = npy.column_stack((dp_tar, yita_tar2))
    plt.figure(figsize=figsize_dflt)
    plt.plot(dp_tar*1e6, yita_tar2, 'm-', markersize = 5)
    plt.plot(dp_tar*1e6, yita_tar, 'r-', markersize = 5)
    plt.legend(['t=' + str(ts[it2_]) + 's', 't=' + str(ts[-1]) + 's'])#, fontsize=tick_size)
    plt.xscale('log')
    plt.xlabel(r'$\it{d}_p$ ($\mu$m)')#,fontsize=label_size)
    plt.ylabel(r'$\it{\eta}$ (%)')#, fontsize=label_size)
    plt.title('Agglomeration efficiency')#, fontsize=label_size)
    
#    plt.ylabel(r'$\it{\beta}$ ($m^3$/s)')