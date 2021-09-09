import numpy as np
from numpy.core.fromnumeric import squeeze
from numpy.core.numeric import NaN
import pandas as pd
import seaborn as sns
from scipy import io as sio
from scipy.stats import zscore
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join as pjoin

# ------------------------------------------------------------------
# functions written as utilities for the report_figs.ipynb notebook.
# ------------------------------------------------------------------

# BVP related functions (cfr report_figs.ipynb)

def get_bvp(path):
    matfile = sio.loadmat(path,squeeze_me=True)
    sig  = matfile['ppg']['sig_preproc'].item()
    time = matfile['ppg']['t'].item()
    fs   = matfile['ppg']['fs'].item()
    return time,sig,fs

def bvp_plot(time,sig,size):
    fig,ax = plt.subplots(figsize=size)
    plt.plot(time,sig,lw=1,color='black')
    plt.xlim(0,20)
    plt.xlabel('time [s]')
    plt.ylabel('amplitude [a.u.]')
    plt.tight_layout()
    return fig,ax

def bvp_overview(mid,path_from,path_to,size):
    path_mid = pjoin(path_from,mid)
    time,sig,_ = get_bvp(path_mid)

    fig,ax = bvp_plot(time,sig,size)
    path_fig = pjoin(path_to,mid)
    plt.savefig(path_fig,dpi=80)
    plt.close()
    return

def bvp_subplot(ax,time,sig,fs,T,apply_zscore=True):
    time = time[0:T*fs]
    bvp = sig[0:T*fs]
    if apply_zscore:
        bvp  = zscore(bvp)
    ax.plot(time,bvp,lw=.1,color='black')
    ax.set_ylabel('amplitude [a.u.]')
    ax.set_xlim(0,T)
    return ax

def fid_path(path,mid,fid):
    fullname = mid + '_' + fid
    return pjoin(path,fid,fullname)

def bvp_fids(path,mid,fid_str):
    fid_val = dict()
    for fid in fid_str:
        fpath = fid_path(path,mid,fid)
        matfile = sio.loadmat(fpath,squeeze_me=True)
        fid_val[fid] = matfile[fid]
    return fid_val

def fids_scatter(ax,path_sig,path_fid,fid_str,mid,T,apply_zscore=True):
    
    t,sig,fs = get_bvp(pjoin(path_sig,mid))
    fid_val = bvp_fids(path_fid,mid,fid_str)

    time = t[0:T*fs]
    bvp = sig[0:T*fs]
    if apply_zscore:
        bvp  = zscore(bvp)

    ax.plot(time,bvp,lw=1,color='black')
    ax.set_xlim(0,T)

    for fid in fid_str:
        fid_x = fid_val[fid][fid_val[fid] <= T*fs] -1 # -1 adjusts for indexing shift (matlab vs python)
        fid_y = bvp[fid_x]
        fid_t = fid_x/fs
        ax.scatter(fid_t,fid_y,label=fid[1:3],marker='o')
    return ax


# EDA related functions

def eda_fixlen(eda,T):
    condition_t = np.less(eda['time'],T)    
    for i in eda.keys():

        # indexed series (remove events past T)
        if i=='onset':
            condition_k = np.less(eda['onset'],T)
            eda['onset_idx'] = np.extract(condition_k,eda['onset_idx'])
            eda['impulse']   = np.extract(condition_k,eda['impulse'])
            eda['overshoot'] = np.extract(condition_k,eda['overshoot'])

        # time series (remove signal past T)
        elif i not in ['onset_idx','impulse','overshoot']:
            eda[i] = np.extract(condition_t,eda[i])
    return eda

def get_eda(dir_cda,dir_dda,mid,T):

    cda = dict()
    select_file = pjoin(dir_cda,mid)
    matfile = sio.loadmat(select_file,squeeze_me=True)
    cda['org']          = matfile['data']['conductance'].item()
    cda['time']         = matfile['data']['time'].item()
    cda['tonic']        = matfile['analysis']['tonicData'].item()
    cda['phasic']       = matfile['analysis']['phasicData'].item()
    cda['driver']       = matfile['analysis']['driver'].item()
    cda['tonicDriver']  = matfile['analysis']['phasicDriverRaw'].item()
    cda['phasicDriver'] = matfile['analysis']['tonicDriver'].item()
    cda['residual']     = matfile['analysis']['remainder'].item()

    dda = dict()
    select_file = pjoin(dir_dda,mid)
    matfile = sio.loadmat(select_file,squeeze_me=True)
    dda['org']          = matfile['data']['conductance'].item()
    dda['time']         = matfile['data']['time'].item()
    dda['tonic']        = matfile['analysis']['tonicData'].item()
    dda['phasic']       = matfile['analysis']['phasicData'].item()
    dda['driver']       = matfile['analysis']['driver'].item()
    dda['residual']     = matfile['analysis']['remainder'].item()
    dda['onset']        = matfile['analysis']['onset'].item()
    dda['onset_idx']    = matfile['analysis']['onset_idx'].item()
    dda['impulse']      = matfile['analysis']['impulse'].item()
    dda['overshoot']    = matfile['analysis']['overshoot'].item()

    if T==NaN:
        return cda,dda
    elif T<=cda['time'][-1] and T>0:
        cda = eda_fixlen(cda,T)
        dda = eda_fixlen(dda,T)
        return cda,dda 
    else:
        print('incorrect arguments')
        return NaN

def eda_subplot(ax,time,sig,**kwargs):
    ax.plot(time,sig,lw=.1,**kwargs)
    ax.set_ylabel('amplitude [a.u.]')
    ax.set_xlim(0,round(time[-1]))
    return ax

def plot_impulses(dda,ax):
    '''
    still requires some work, doesn't properly reconstruct the signal yet..
    But first, try running 'overview',1 in Ledalab batch to generate .jpg files!
    '''
    for i in range(len(dda['impulse'])):

        # get start and stop values for impulse i
        start = dda['onset_idx'][i]
        stop  = start + len(dda['impulse'][i])

        # clip corresponding signals
        segment_time  = dda['time'][start:stop]
        segment_tonic = dda['tonic'][start:stop]

        # reconstruct segments for plotting
        impulse   = dda['impulse'][i] + segment_tonic
        overshoot = dda['overshoot'][i] + impulse

        ax.plot(segment_time,impulse,lw=.1,color='green',label='impulse estimations')
        ax.plot(segment_time,overshoot,lw=.1,color='red',label='overshoot pulses')
        return ax


# RSP related functions

def get_rsp_raw(path,mid):
    select_file = pjoin(path,mid)
    matfile = sio.loadmat(select_file,squeeze_me=True)
    rsp = matfile['rsp']
    fs = matfile['fs']
    time = np.arange(0,len(rsp)/fs,1/fs)
    return time,rsp,fs

def rsp_raw_subplot(ax,time,rsp,T,fs,apply_zscore=True):
    time = time[0:T*fs]
    rsp = rsp[0:T*fs]
    if apply_zscore:
        rsp  = zscore(rsp)
    ax.plot(time,rsp,lw=1,color='black')
    ax.set_ylabel('amplitude [a.u.]')
    ax.set_xlim(0,T)
    return ax

def get_rsp(path,mid):
    select_file = pjoin(path,mid)
    matfile = sio.loadmat(select_file,squeeze_me=True)
    rsp = dict()
    rsp['sig'] = matfile['respiration']['respiration'].item()
    rsp['time'] = matfile['respiration']['time'].item()
    rsp['inh_x'] = matfile['respiration']['inh_onsets'].item()
    rsp['exh_x'] = matfile['respiration']['exh_onsets'].item()
    rsp['inh_t'] = rsp['time'][rsp['inh_x']]
    rsp['exh_t'] = rsp['time'][rsp['exh_x']]
    rsp['inh_y'] = rsp['sig'][rsp['inh_x']]
    rsp['exh_y'] = rsp['sig'][rsp['exh_x']]
    return rsp

def rsp_subplot(ax,rsp,fixlength,apply_zscore=True):

# work in progress
#    if apply_zscore:
#        rsp['sig']  = zscore(rsp['sig'])

    ax.plot(rsp['time'],rsp['sig'],lw=1,color='black')
    ax.scatter(rsp['inh_t'],rsp['inh_y'],marker='.',label='inhale onsets',color='green')
    ax.scatter(rsp['exh_t'],rsp['exh_y'],marker='.',label='exhale onsets',color='red')
    ax.set_xlim(0,fixlength)
    ax.set_ylabel('amplitude [a.u.]')
    # ax.legend()
    return ax