"""loads in  the foot forward/foot back indices (FFi, FBi) and  outputs the 
longest chunk of continuous events from the initial index (FFi_mod, FBi_mod). 
"""

__author__ = "Reginaldo K Fukuchi, https://github.com/regifukuchi"
__version__ = "1.0.1"
__license__ = "MIT"

import numpy as np

def largest_block(FFi, FBi):
    '''
    Parameters
    ----------
    FFi : list
        Indices of foot-forward events for entire length of trial.
    FBi : list
        Indices of foot-back events for entire length of trial..

    Returns
    -------
    FFI_MOD/FBI_MOD: array
        Abbreviated index of foot-forward/foot-events containing only the longest 
        set of continuous events. 
    BLOCK_START/BLOCK_END: int
        Frame index of when the largest block of continuous events begins and ends
        
    This code was based on the Matlab code largest_block.m written by Allan Brett.

    '''

    # Combine and sort the FF and FB
    allsort = np.vstack((np.array([FFi,np.zeros(len(FFi))]).T,np.array([FBi,np.ones(len(FBi))]).T))
    inds = np.argsort(allsort[:,0])
    allsort = allsort[inds,:].astype(int)
    
    # Remove trailing FF
    while allsort[-1,1] != 1:
        allsort = allsort[:-1,:]
        
    allsort_bak = allsort
    
    i=-1
    skip = -1
    longest_length =[]
    
    # series of while loops will search for longest continuous chunk of data
    #first while loop adds up the length of a continuous segments adding in the
    #indicies that are skipped (because they contain the dicontinuity)
    while (np.sum(longest_length)+skip < allsort_bak.shape[0]) and allsort.shape[0]>1:
    
        # points, in case of two 1s run below
        idx_skip = -1
        k = 0
    
        # allsort must start with a 0
        while allsort[0,1] == 1:
            allsort = allsort[1:,:]
            #skip is an overall counter for the main while loop in conjuction with longest_length
            skip = skip + 1
    
            # idx_skip keeps track of when values are skipped for the purpose of indexing
            idx_skip = idx_skip + 1
    
            # remove discontinuities occuring at the start of allsort
            while (allsort[k,1]==allsort[k+1,1]) and (k+1<allsort.shape[0]):
                allsort = allsort[k+2:,:]
                skip = skip + 2
                idx_skip = idx_skip + 2
    
        k = 2
        while k<=len(allsort) and allsort[0:k:2,1].mean()==0 and allsort[1:k:2,1].mean()==1:
            k = k + 2
    
        i = i+1
    
        # we don't want to use a possibly erroneous point in the data and we
        # must end the sequence on a 1, so when the dicontinuity occurs with two
        # 1s in a row, we must roll back by 2
    
        if (allsort.shape[0] > k) and (k > 2):
            if (allsort[k-2,1]==1) and (allsort[k-3,1]==1):
                allsort = allsort[:k-5,:]
            else:
                allsort = allsort[:k-3,:]
    
        #for the special case where there are two discontinuities of 0s in a row
        index = np.empty(shape=(2,allsort_bak.shape[0])) * np.NaN
        if k==2 and allsort[0,1]==0:
            allsort = allsort[2:,:]
            longest_length.append(0)
    
            # we want to index one passed the discontinuity
            if i==0:
                # if this occurs for the first index, only includes values skipped
                index[0,0] = idx_skip
                index[1,0] = idx_skip
            else:
                index[0,i] = index[1,i-1] + idx_skip+1
                index[1,i] = index[1,i-1] + idx_skip+1
    
            skip = skip + 2
    
        else:
            # otherwise count as normal
            longest_length.append(allsort.shape[0])
    
            # create ordered index of where continuous chunks occur
            if i==0:
                index[0,0] = 1 + idx_skip
                index[1,0] = longest_length[i] + idx_skip
            else:
                index[1,i] = index[1,i-1]+3+idx_skip
    
                # Longest_length can only be 0 when two discontinuities of 1s happen in a row, 
                #below accounts that the index end needs to still progress by 1 (but longest_length
                #still needs to be 0 for the main counter)
                if longest_length[i] > 0:
                    index[1,i] = index[0,i]+longest_length[i]-1
                else:
                    index[1,i] = index[0,i]+longest_length[i]
    
            #reset allsort for next loop iteration to be passed the discontinuity
            index = index[~np.isnan(index)].reshape((2,-1))
            allsort = allsort_bak[index[1,i].astype(int)+3:,:]
    
            #however we want to skip passed the discontinuity to the next footfall.
            #This entails skipping the discontinuity (for example if the
            #discontinuity is two FF, we skip over these two values
            skip = skip + 2
            
    #determine which index has the largest continuous block
    longest_index = np.argmax(np.diff(index, axis=0))
    
    # reorder allsort to contain only this block
    allsort_longest = allsort_bak[int(index[0,longest_index]):int(index[1,longest_index]+1),:]
    allsort = allsort_longest
    
    # Break back into components
    FFi_mod = allsort[allsort[:,1]==0,0]
    FBi_mod = allsort[allsort[:,1]==1,0]
    
    # want to track frame numbers of when the block on continuous data starts
    #and ends
    block_start = FFi_mod[0]
    block_end   = FBi_mod[-1]
    
    return FFi_mod, FBi_mod, block_start, block_end