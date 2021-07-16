# -*- coding: utf-8 -*-
"""
2021/07/09 JK
Overall suite2p process
"""

#%%
# session_selection.py
'''
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

Results in /plane_x/regParam00 ~ /regParam07 folders
           /plane_x/regResult.npy
Take notes of selected sessions in an excel file

'''


# session_QC.py
'''
For *, in another script
Save the result for each mouse (because it contains between-plane comparison)
'regQC.npy'
regQC['depth']
    list of 8 plane results, for each session in each plane
regQC['withinPlane']
    list of 8 plane results, for each session in each plane
regQC['betweenPlane']
    list of 6 results, 3 for upper volume 1-plane difference, 3 for lower volume 1-plane difference
    ** The depth difference can be different bewteen mice
'''


# register_selected.py
'''
Run registration for all the frames from the selected sessions

Results in /selected/plane0\ops.npy and  /../..\data.bin
'''



# roi_param_search.py
'''
Search the best ROI selection parameter.
    ops['threshold']
    
    ops['allow_overlap'] = False
    ops['max_overlap'] = 0.3
    ops['threshold_scaling'] = 0.01-1.0, dynamic search. 9 iterations, 0.05 resolution
        0.2 resolution, 0.1 resolution around the peak, 0.05 resolution around the previous peak
        
Results in /selected/plane0/roiParam00 ~ /roiParam08 folders
/selected/plane0\ops.npy /../..\F.npy, Fneu.npy, iscell.npy, spks.npy, stat.npy all these copied from the best parameter folder
'''



# Use suite2p
'''
Manual curation
'''


#%% Test
# register_all_sessions.py
'''
Run suite2p registration for all sessions
With default registration parameters (nonrigid)
Combine multiple files from the same session
'''
