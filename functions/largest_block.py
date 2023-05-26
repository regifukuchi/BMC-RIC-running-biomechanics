# Generated with SMOP  0.41
from libsmop import *
# 

    
@function
def largest_block(FFi=None,FBi=None,*args,**kwargs):
    varargin = largest_block.varargin
    nargin = largest_block.nargin

    ###########################################################################
# Function LARGEST_BLOCK loads in  the foot forward/foot back indices 
# (FFi, FBi) and  outputs the longest chunk of continuous events from the
# initial index (FFi_mod, FBi_mod).
# Function also outputs the index of the start/end of the entire continuous
# block (block_start, block_end).
    
    # This function is necessary as there are many reason why either a FF or FB
# event could be missed. An imbalanced index (more FFs or more FBs) causes
# many downstream problems with the code.
    
    
    #  INPUTS
#  --------
#  FFI/FBI (mat):   Frame index of foot-forward/foot-back events for
#                   entire length of trial
    
    #  OUTPUTS
#  -------
#  FFI_MOD/FBI_MOD (mat):   Abbreviated index of foot-forward/foot-events 
#                           containing only the longest set of continuous 
#                           events.
    
    #  BLOCK_START/BLOCK_END (int): Frame index of when the largest block of 
#                               continuous events begins and ends
    
    # Copyright (C) 2016-2023 Allan Brett and The Running Injury Clinic
###########################################################################
##
# Combine and sort the FF and FB
    allsort=concat([[ravel(FFi),zeros(length(FFi),1)],[ravel(FBi),ones(length(FBi),1)]])

    __,inds=sort(allsort(arange(),1),nargout=2)
34
    allsort=allsort(inds,arange())
35
    # Remove trailing FF
    while allsort(end(),2) != 1:

        allsort[end(),arange()]=[]
39

    
    allsort_bak=copy(allsort)
42
    # Confirm alternating status and remove trailing chunks that aren't
# ordered
    
    i=0
47
    skip=0
48
    longest_length=0
49
    #series of while loops will search for longest continuous chunk of data
#first while loop adds up the length of a continuous segments adding in the
#indicies that are skipped (because they contain the dicontinuity)
    while (sum(longest_length) + skip < size(allsort_bak,1)) and size(allsort,1) > 1:

        idx_skip=0
56
        k=1
59
        while allsort(1,2) == 1:

            allsort[1,arange()]=[]
63
            #conjuction with longest_length
            skip=skip + 1
66
            #purpose of indexing
            idx_skip=idx_skip + 1
70
            while allsort(k,2) == allsort(k + 1,2) and k + 1 < size(allsort,1):

                allsort[arange(k,k + 1),arange()]=[]
74
                skip=skip + 2
75
                idx_skip=idx_skip + 2
76


        k=2
80
        while k <= length(allsort) and mean(allsort(arange(1,k,2),2)) == 0 and mean(allsort(arange(2,k,2),2)) == 1:

            k=k + 2
82

        i=i + 1
85
        #must end the sequence on a 1, so when the dicontinuity occurs with two
    #1s in a row, we must roll back by 2
        if size(allsort,1) > k and k > 2:
            if allsort(k - 1,2) == 1 and allsort(k - 2,2) == 1:
                allsort=allsort(arange(1,k - 4),arange())
92
            else:
                allsort=allsort(arange(1,k - 2),arange())
94
        #for the special case where there are two discontinuities of 0s in a row
        if k == 2 and allsort(1,2) == 0:
            allsort[arange(1,2),arange()]=[]
100
            longest_length[i,1]=0
101
            if i == 1:
                #if this occurs for the first index, only includes values
            #skipped
                index[1,1]=idx_skip
107
                index[2,1]=idx_skip
108
            else:
                index[1,i]=index(2,i - 1) + idx_skip + 1
110
                index[2,i]=index(2,i - 1) + idx_skip + 1
111
            skip=skip + 2
113
        else:
            #otherwise count as normal
            longest_length[i,1]=size(allsort,1)
118
            if i == 1:
                index[1,1]=1 + idx_skip
123
                index[2,1]=longest_length(i,1) + idx_skip
124
            else:
                index[1,i]=index(2,i - 1) + 3 + idx_skip
126
                #two discontinuities of 1s happen in a row, below accounts that
            #the index end needs to still progress by 1 (but longest_length
            #still needs to be 0 for the main counter)
                if longest_length(i,1) > 0:
                    index[2,i]=index(1,i) + longest_length(i,1) - 1
133
                else:
                    index[2,i]=index(1,i) + longest_length(i,1)
135
            #reset allsort for next loop iteration to be passed the discontinuity
            allsort=allsort_bak(arange(index(2,i) + 3,end()),arange())
141
            #This entails skipping the discontinuity (for example if the
        #discontinuity is two FF, we skip over these two values
            skip=skip + 2
148

    
    #determine which index has the largest continuous block
    __,longest_index=max(diff(index),nargout=2)
154
    #reorder allsort to contain only this block
    allsort_longest=allsort_bak(arange(index(1,longest_index),index(2,longest_index)),arange())
157
    allsort=copy(allsort_longest)
158
    # Break back into components
    FFi_mod=allsort(allsort(arange(),2) == 0,1)
161
    FBi_mod=allsort(allsort(arange(),2) == 1,1)
162
    #want to track frame numbers of when the block on continuous data starts
#and ends
    block_start=FFi_mod(1,1)
166
    block_end=FBi_mod(end(),1)
167
    return FFi_mod,FBi_mod,block_start,block_end
    
if __name__ == '__main__':
    pass
    