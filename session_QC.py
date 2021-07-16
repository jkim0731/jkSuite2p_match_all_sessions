'''
2021/07/10 JK

After session_selection.py
Using the saved file 'regParam.npy'
    regParam['selectSessionName']
    regParam['selectSessionFilename']
    regParam['phaseCorr']
    regParam['pixelCorr']
    regParam['meanImgSession']
Also need z-stack of each mouse
Depth diff info

Depth estimation using the registered images
Upper limiat and lower limit (Between planes)

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