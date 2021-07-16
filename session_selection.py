'''
2021/07/10 JK

First, search the best parameters for registration.
Register the reference session (naive 7-angle) with the default parameter.
Register 50 evenly spaced frames from each session to the reference session.
Test 4 combinations of registration parameters
Parameter sets: 
    ops['nonrigid'] = True
    ops['maxregshift'] = 0.3
    ops['snr_threshold'] = [1.2, 1.4]
    ops['maxregshiftNR'] = [5,10]
    ops['block_size'] = [[48,48], [32,32]]
Select the best parameter set from mean phase correlation values
Select FOV-matched sessions using phase correlation and pixel correlation values, as well as using visual inspection
    *Depth estimation using the registered images
    *Upper limiat and lower limit (Between planes)
Save the result for each plane
'regParam.npy'
regParam['allSessionName']
regParam['allSessionFilename']
regParam['selectSessionName']
regParam['selectSessionFilename']
regParam['nonrigid']
regParam['maxregshift']
regParam['snr_threshold']
regParam['maxregshiftNR']
regParam['block_size']
regParam['phaseCorr']
regParam['pixelCorr']
regParam['meanImgSession']
    Only save that from the best parameter
* session_QC.py

Results in /plane_x/regParam00 ~ /regParam07 folders

Based on 210621_auto_reg_param_test_fromAll.py and 210624_check_auto_reg_param.py
'''

import h5py
import numpy as np
from suite2p.run_s2p import run_s2p, default_ops
import os, glob, shutil
import matplotlib.pyplot as plt
from skimage import exposure
import napari
from suite2p.io.binary import BinaryFile
import time
from functools import reduce

import gc
gc.enable()

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

framesPerSession = 50
#%% This part can run in multiple mice
for mi in [0,1,2,3,5,6,7,8,9,10,11]:
    mouse = mice[mi]
    refS = refSessions[mi]
    for planeNum in range(1,9):
        mouseDir = f'{h5Dir}{mouse:03}/'
        planeDir = f'{mouseDir}plane_{planeNum}/'
        
        # #%% First, make ref session registration
        # Check if there's one. If not, then run suite2p registration with the default ops (two-step registration)
        refOpsFn= f'{planeDir}{refS:03}/plane0/ops.npy'
        if not os.path.isfile(refOpsFn):
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
            
            refFnList = glob.glob(f'{planeDir}{mouse:03}_{refS:03}_*_plane_{planeNum}.h5')
            db = {'h5py': refFnList,
                'h5py_key': ['data'],
                'data_path': [],
                'save_path0': planeDir,
                'save_folder': f'{refS:03}',
                # 'fast_disk': f'{planeDir}/{refS:03}',
                'fast_disk': f'{fastDir}', # This better be in SSD
                'roidetect': False,
            }
            run_s2p(ops,db)
        ops = np.load(refOpsFn, allow_pickle=True).item()
        refImg = ops['meanImg']
        
        # #%% pick linearly separated 'framesPerSession' frames from each session, make a new h5 file
        # multiple sbx files (h5 files) from a single session should be merged
        # Piezo sessions should be merged first
        # For mouse > 50, Spont sessions should also be merged
        tempFnList = glob.glob(f'{planeDir}{mouse:03}_*_plane_{planeNum}.h5')    
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
            
        # #%% Select linearly separated frames from each session.
        # Combine if multiple files in a session
        f = h5py.File(tempFnList[0], 'r')
        data = f['data']
        _, height, width = data.shape
        wfn = f'{planeDir}selected.h5'
        if not os.path.isfile(wfn):
            newdata = np.zeros((len(sessionNames)*framesPerSession, height, width), dtype = np.uint16)
            for i, fnlist in enumerate(sessionFiles):
                for j, fn in enumerate(fnlist):
                    f = h5py.File(fn, 'r')
                    if j == 0:
                        data = f['data']
                    else:
                        data = np.concatenate((data, f['data']), axis=0)
                numFrames, height, width = data.shape
                frames = np.linspace(0,numFrames-1, num=framesPerSession, dtype=int)
                for j in range(len(frames)):
                    newdata[i*framesPerSession+j, :, :] = data[frames[j],:,:]    
            with h5py.File(wfn, 'w') as wf:
                wf.create_dataset('data', data=newdata, dtype='uint16')
    
        # #%% Parameter search
        ops = default_ops()
        ops['tau'] = 1.5
        ops['look_one_level_down'] = False
        ops['do_bidiphase'] = True
        ops['batch_size'] = 5000
        ops['two_step_registration'] = True
        ops['keep_movie_raw'] = True
        ops['smooth_sigma_time'] = 2
        ops['move_bin'] = True
        ops['fs'] = freq[mi]
        ops['zoom'] = zoom[mi]
        ops['umPerPix'] = 1.4/ops['zoom']
        ops['force_refImg'] = True
        ops['refImg'] = refImg
        
        maxregshiftNRList = [5, 10]
        snr_threshList = [1.2, 1.4]
        block_sizeList = [[48, 48], [32, 32]]        
        
        paramSetInd = 0
        
        ops['nonrigid'] = True
        ops['maxregshift'] = 0.3
        for mrsn in maxregshiftNRList:
            for st in snr_threshList:
                for bs in block_sizeList:
                    regParamDn = f'regParam{paramSetInd:02}'
                    
                    ops['maxregshiftNR'] = mrsn
                    ops['block_size'] = bs
                    ops['snr_thresh'] = st
                    
                    db = {'h5py': wfn,
                        'h5py_key': ['data'],
                        'data_path': [],
                        'save_path0': planeDir,
                        'save_folder': regParamDn,
                        # 'fast_disk': f'{planeDir}/{regParamDn}', 
                        'fast_disk': f'{fastDir}', # This better be in SSD
                        'roidetect': False,
                        'testFileList': sessionFiles,
                        'testSessionNames': sessionNames,
                        'framesPerSession': framesPerSession,
                    }
                    run_s2p(ops,db)
                    paramSetInd += 1
        
        
        
        
        
        
        
        
        
        
        
#%% Based on the mean phase correlation values, select the best parameter set
# Perform this per mouse
# Requires manual processing
# Based on 210624_check_auto_reg_param.py
mi = 0
mouse = mice[mi]
mouseDir = f'{h5Dir}{mouse:03}/'
planes = range(1,9)

mimgPlanes = []
viewer = napari.Viewer()
for planeNum in planes:
    planeDir = f'{mouseDir}plane_{planeNum}/'
        
    paramSetInd = 0
    regParamDn = f'regParam{paramSetInd:02}'
    dataDir = f'{planeDir}/{regParamDn}/plane0/'
    binfn = f'{dataDir}data.bin'
    opsfn = f'{dataDir}ops.npy'
    
    ops = np.load(opsfn, allow_pickle = True).item()
    framesPerSession = ops['framesPerSession']
    Ly = ops['Ly']
    Lx = ops['Lx']
    nframes = ops['nframes']
    framesPerSession = ops['framesPerSession']
    numSessions = len(ops['testSessionNames'])
    meanImgs = np.zeros((numSessions,Ly,Lx))
    with BinaryFile(Ly = Ly, Lx = Lx, read_filename = binfn) as f:
        for i in range(numSessions):
            inds = np.arange(i*framesPerSession,(i+1)*framesPerSession)
            frames = f.ix(indices=inds).astype(np.float32)
            tempMimg = clahe_each(frames.mean(axis=0))
            meanImgs[i,:,:]= tempMimg
            # viewer.add_image(tempMimg, name = ops['testSessionNames'][i])
    viewer.add_image(meanImgs, name = f'plane#{planeNum}')
    mimgPlanes.append(meanImgs)
        
#%% Set top bottom left right margin in "pixel number" and test the result
topMargin = 70
bottomMargin = 350
leftMargin = 40
rightMargin = 650
Ly = bottomMargin - topMargin
Lx = rightMargin - leftMargin
viewer = napari.Viewer()
for pi, planeNum in enumerate(planes):
    viewer.add_image(mimgPlanes[pi][:,topMargin:bottomMargin, leftMargin:rightMargin], name = f'plane#{planeNum}')

#%% Compare with the ref session, collect and save all reg params and results
# Do this per mouse

mi = 0
mouse = mice[mi]
mouseDir = f'{h5Dir}{mouse:03}/' 
for planeNum in range(1,9):
    tic = time.time()
    savefn = 'regResult.npy'
    planeDir = f'{mouseDir}plane_{planeNum}/'
    
    allPhaseCorr = []
    allPixCorr = []
    allMeanImg = []
    allMrsNR = []
    allSnrThresh = []
    allBlockSize = []
    for paramNum in range(8):            
        regParamDn = f'regParam{paramNum:02}'
        dataDir = f'{planeDir}/{regParamDn}/plane0/'
        binfn = f'{dataDir}data.bin'
        opsfn = f'{dataDir}ops.npy'
        ops = np.load(opsfn, allow_pickle = True).item()
        
        allMrsNR.append(ops['maxregshiftNR'])
        allSnrThresh.append(ops['snr_thresh'])
        allBlockSize.append(ops['block_size'])
        
        framesPerSession = ops['framesPerSession']
        numSessions = len(ops['testSessionNames'])
        meanImgs = np.zeros((numSessions,Ly,Lx))
        with BinaryFile(Ly = ops['Ly'], Lx = ops['Lx'], read_filename = binfn) as f:
            for i in range(numSessions):
                inds = np.arange(i*framesPerSession,(i+1)*framesPerSession)
                frames = f.ix(indices=inds).astype(np.float32)
                tempMimg = clahe_each(frames.mean(axis=0))
                meanImgs[i,:,:]= tempMimg[topMargin:bottomMargin, leftMargin:rightMargin]
       
        phaseCorr = np.zeros(numSessions)
        pixCorr = np.zeros(numSessions)
        
        refSession = refSessions[mi]
        refSname = f'{mouse:03}_{refSession:03}_'
        refSi = [i for i,fn in enumerate(ops['testSessionNames']) if refSname in fn]
        refImg = meanImgs[refSi[0],:,:]
        for i in range(numSessions):
            img1 = meanImgs[i,:,:]
            _, _, _, phaseCorr[i], _ = phase_corr(img1, refImg)
            pixCorr[i] = np.corrcoef(img1.flatten(), refImg.flatten())[0,1]
        
        allPhaseCorr.append(phaseCorr)
        allPixCorr.append(pixCorr)
        allMeanImg.append(meanImgs)
    maxInd = np.argmax([np.mean(pc) for pc in allPhaseCorr]) # same as max param num
    
    result = {}
    result['maxInd'] = maxInd
    result['phaseCorr'] = allPhaseCorr[maxInd]
    result['pixCorr'] = allPixCorr[maxInd]
    result['testSessionNames'] = ops['testSessionNames'] # same for all parameter numbers
    result['refSession'] = refSname # same for all parameters
    result['meanImg'] = allMeanImg[maxInd]
    result['mrsNR'] = allMrsNR[maxInd]
    result['snr_thresh'] = allSnrThresh[maxInd]
    result['block_size'] = allBlockSize[maxInd]
    result['refImg'] = ops['refImg']
    
    np.save(f'{planeDir}{savefn}', result)
    
    toc= time.time()
    timeInMin = np.round((toc-tic)/60)
    print(f'{timeInMin} minutes for plane#{planeNum}.')

#%% Visual inspection and selection of matched FOV (or nonmatched FOV)
# Do this per volume
# Requires manual processing

mi = 0
planeList = range(1,5)
# planeList = range(5,9)

clist  = ['k', 'b', 'c', 'g']

mouse = mice[mi]
mouseDir = f'{h5Dir}{mouse:03}/' 
tempPlane = planeList[0]
resultFn = f'{mouseDir}plane_{tempPlane}/regResult.npy'
result = np.load(resultFn, allow_pickle=True).item()
numSessions = result['meanImg'].shape[0]
allPhaseCorr = np.zeros((len(planeList),numSessions))
allPixCorr = np.zeros((len(planeList),numSessions))
viewer = napari.Viewer()
f, (ax1, ax2) = plt.subplots(1,2)
for i, planeNum in enumerate(planeList):
    resultFn = f'{mouseDir}plane_{planeNum}/regResult.npy'
    result = np.load(resultFn, allow_pickle=True).item()
    viewer.add_image(result['meanImg'], name=f'plane {planeNum}')
    allPhaseCorr[i,:] = result['phaseCorr']
    allPixCorr[i,:] = result['pixCorr']
    ax1.plot(result['phaseCorr'], 'o-', color=clist[i])
    ax2.plot(result['pixCorr'], 'o-', color=clist[i])
maxPhaseCorr = np.max([pc for pc in allPhaseCorr.flatten() if pc != 1])
minPhaseCorr = np.min([pc for pc in allPhaseCorr.flatten() if pc != 1])
maxPixCorr = np.max([pc for pc in allPixCorr.flatten() if pc != 1])
minPixCorr = np.min([pc for pc in allPixCorr.flatten() if pc != 1])
ax1.legend([f'plane {pn}' for pn in planeList])
ax1.set_title('Phase correlation')
ax2.set_title('Pixel correlation')
f.suptitle(f'JK{mouse:03}')
f.tight_layout()
ax1.set_ylim(minPhaseCorr, maxPhaseCorr)
ax2.set_ylim(minPixCorr, maxPixCorr)



