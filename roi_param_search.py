# roi_param_search.py
'''
2021/07/10 JK
Search the best ROI selection parameter.
    ops['threshold']
    
    ops['allow_overlap'] = False
    ops['max_overlap'] = 0.3
    ops['threshold_scaling'] = 0.01-1.0, dynamic search. 9 iterations, 0.05 resolution
        0.2 resolution, 0.1 resolution around the peak, 0.05 resolution around the previous peak

Results in /selected/plane0/roiParam00 ~ /roiParam08 folders
/selected/plane0\ops.npy /../..\F.npy, Fneu.npy, iscell.npy, spks.npy, stat.npy all these copied from the best parameter folder
'''