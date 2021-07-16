
"""
Run suite2p registration for all sessions
With default registration parameters (nonrigid)
Combine multiple files from the same session
2021/07/15 JK

"""
import h5py
import numpy as np
from suite2p.run_s2p import run_s2p, default_ops
import os, glob, shutil
import matplotlib.pyplot as plt
from skimage import exposure
import napari
from suite2p.io.binary import BinaryFile
import time

# CLAHE each mean images
def clahe_each(img):
    newimg = (img - np.amin(img)) / (np.amax(img) - np.amin(img)) * (2**16-1)
    newimg = exposure.equalize_adapthist(newimg.astype(np.uint16))
    return newimg

def phase_corr(a,b):
    if a.shape != b.shape:
        raise('Dimensions must match')
    R = np.fft.fft2(a) * np.fft.fft2(b).conj()
    R /= np.absolute(R)
    r = np.absolute(np.fft.ifft2(R))
    ymax, xmax = np.unravel_index(np.argmax(r), r.shape)
    cmax = np.amax(r)
    center = r[0, 0]    
    return ymax, xmax, cmax, center, r

h5Dir = 'D:/TPM/JK/h5/'
fastDir = 'C:/JK/' # This better be in SSD
mice =          [25,    27,   30,   36,     37,     38,     39,     41,     52,     53,     54,     56]
refSessions =   [4,     3,    3,    1,      7,      2,      1,      3,      3,      3,      3,      3]
zoom =          [2,     2,    2,    1.7,    1.7,    1.7,    1.7,    1.7,    1.7,    1.7,    1.7,    1.7]
freq =          [7.7,   7.7,  7.7,  7.7,    6.1,    6.1,    6.1,    6.1,    7.7,    7.7,    7.7,    7.7]

ops = default_ops()
ops['tau'] = 1.5
ops['look_one_level_down'] = False
ops['do_bidiphase'] = True
ops['nimg_init'] = 100
ops['batch_size'] = 5000
ops['two_step_registration'] = True
ops['keep_movie_raw'] = True
ops['smooth_sigma_time'] = 2
ops['move_bin'] = True

for mi in [0,3,8]:
    mouse = mice[mi]
    for pn in range(1,9):
        mouseDir = f'{h5Dir}{mouse:03}/'
        planeDir = f'{mouseDir}plane_{pn}/'
        
        # Make a list of session names and corresponding files
        tempFnList = glob.glob(f'{planeDir}{mouse:03}_*_plane_{pn}.h5')    
        fnames = [fn.split('\\')[1].split('.h5')[0] for fn in tempFnList]
        midNum = np.array([int(fn.split('\\')[1].split('_')[1]) for fn in tempFnList])
        trialNum = np.array([int(fn.split('\\')[1].split('_')[2][0]) for fn in tempFnList])
        regularSi = np.where(midNum<1000)[0]
        spontSi = np.where( (midNum>5000) & (midNum<6000) )[0]
        piezoSi = np.where(midNum>9000)[0]
        
        if np.any(spontSi): 
            spontTrialNum = np.unique(trialNum[spontSi]) # used only for mouse > 50
        
        if np.any(piezoSi):
            piezoTrialNum = np.unique(trialNum[piezoSi])
        
        sessionNum = np.unique(midNum)
        regularSni = np.where(sessionNum < 1000)[0]
        
        sessionNames = []
        sessionFiles = []
        
        for sni in regularSni:
            sn = sessionNum[sni]
            sname = f'{mouse:03}_{sn:03}_'
            sessionNames.append(sname)
            sessionFiles.append([fn for fn in tempFnList if sname in fn])
        if mouse < 50:
            for si in spontSi:
                sessionNames.append(tempFnList[si].split('\\')[1].split('.h5')[0][:-8])
                sessionFiles.append([tempFnList[si]])
        else:
            for stn in spontTrialNum:
                sn = midNum[spontSi[0]]
                sname = f'{mouse:03}_{sn}_{stn}'
                sessionNames.append(sname)
                sessionFiles.append([fn for fn in tempFnList if sname in fn])
        for ptn in piezoTrialNum:
            sn = midNum[piezoSi[0]]
            sname = f'{mouse:03}_{sn}_{ptn}'
            sessionNames.append(sname)
            sessionFiles.append([fn for fn in tempFnList if sname in fn])
        
        # Run suite2p nonrigid registration for each session
        for sn, sf in zip(sessionNames, sessionFiles):
            midNum = sn.split('_')[1]
            if int(midNum) < 1000:
                sname = midNum
            else: # spont and piezo sessions
                trialNum = sn.split('_')[2]
                sname = f'{midNum}_{trialNum}'
            opsFn = f'{planeDir}{sname}/plane0/ops.npy'
            binFn = f'{planeDir}{sname}/plane0/data.bin'
            if (not os.path.isfile(opsFn)) & (not os.path.isfile(binFn)): # if suite2p was not run in this session
                db = {'h5py': sf,
                    'h5py_key': ['data'],
                    'data_path': [],
                    'save_path0': planeDir,
                    'save_folder': f'{sname}',
                    'fast_disk': f'{fastDir}',
                    'roidetect': False,
                }
                run_s2p(ops,db)
                rawbinFn = f'{planeDir}{sname}/plane0/data_raw.bin'
                os.remove(rawbinFn)
