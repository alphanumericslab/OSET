# Python function that returns the mid points relative to the R-peaks. (using this convention that a heartbeat is starting at the middle point between two R peaks and ends up at the middle point between the next two R peaks.)

#Input: peaks - the R-peaks positions
#     : signal_length - length of the signal
#Output: mids - the middle points positions


import numpy as np
def MiddleDetector(peaks,signal_length):
    mids = []
    for i in range(len(peaks)-1):        
        mid = peaks[i] + int(np.floor((peaks[i+1] - peaks[i])/2))
        mids.append(mid)
        
    mid = peaks[0] - int(np.floor((peaks[1] - peaks[0])/2))
    if(mid >= 0):
        mids.insert(0, mid)
    
    mid = peaks[-1] + int(np.floor((peaks[-1] - peaks[-2])/2))
    if(mid <= signal_length):
        mids.append(mid)        
    mids = np.array(mids)    
    return mids    

