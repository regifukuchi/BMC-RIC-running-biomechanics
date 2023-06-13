%function [ FFi_mod, FBi_mod, block_start, block_end ] = largest_block( FFi, FBi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function LARGEST_BLOCK loads in  the foot forward/foot back indices 
% (FFi, FBi) and  outputs the longest chunk of continuous events from the
% initial index (FFi_mod, FBi_mod).
% Function also outputs the index of the start/end of the entire continuous
% block (block_start, block_end).
%
% This function is necessary as there are many reason why either a FF or FB
% event could be missed. An imbalanced index (more FFs or more FBs) causes
% many downstream problems with the code. 
%
%
%  INPUTS
%  --------
%  FFI/FBI (mat):   Frame index of foot-forward/foot-back events for
%                   entire length of trial
%
%  OUTPUTS
%  -------
%  FFI_MOD/FBI_MOD (mat):   Abbreviated index of foot-forward/foot-events 
%                           containing only the longest set of continuous 
%                           events.       
%
%  BLOCK_START/BLOCK_END (int): Frame index of when the largest block of 
%                               continuous events begins and ends

% Copyright (C) 2016-2023 Allan Brett and The Running Injury Clinic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUTS
FFi = [119, 249, 379, 511, 641, 771, 902, 1031, 1162, 1292, 1422, 1553, 1683,... 
       1812, 1943, 2075, 2205, 2336, 2466, 2597, 2729, 2859, 2990, 3120, 3249,... 
       3379, 3510, 3640, 3770, 3901, 4031, 4162, 4291, 4420, 4549, 4680, 4809, 4938];
FBi = [185, 314, 436, 577, 707, 837, 967, 1097, 1229, 1358, 1489, 1619, 1745,... 
       1877, 2009, 2142, 2263, 2401, 2523, 2664, 2794, 2925, 3056, 3186, 3314,... 
       3445, 3574, 3705, 3836, 3966, 4096, 4226, 4356, 4485, 4616, 4746, 4873];
%%
% Combine and sort the FF and FB
allsort = [FFi(:) zeros(length(FFi),1); FBi(:) ones(length(FBi),1)];
[~,inds] = sort(allsort(:,1));
allsort = allsort(inds,:);

% Remove trailing FF
while allsort(end,2) ~= 1
    allsort(end,:) = [];
end

allsort_bak = allsort;

%% Confirm alternating status and remove trailing chunks that aren't
% ordered

i=0;
skip = 0;
longest_length =0;

%% series of while loops will search for longest continuous chunk of data
%first while loop adds up the length of a continuous segments adding in the
%indicies that are skipped (because they contain the dicontinuity)
while (sum(longest_length)+skip < size(allsort_bak,1)) && size(allsort,1)>1
    
    idx_skip = 0;
    
    %points, in case of two 1s run below
    k = 1;
    
    %allsort must start with a 0
    while allsort(1,2) == 1
        allsort(1,:) = [];
        %skip is an overall counter for the main while loop in
        %conjuction with longest_length
        skip = skip+1;
        
        %idx_skip keeps track of when values are skipped for the
        %purpose of indexing
        idx_skip = idx_skip+ 1;
        
        %remove discontinuities occuring at the start of allsort
        while allsort(k,2)==allsort(k+1,2) && k+1<size(allsort,1)
            allsort(k:k+1,:) = [];
            skip = skip + 2;
            idx_skip = idx_skip + 2;
        end
    end
    
    k = 2;
    while k <= length(allsort) && mean(allsort(1:2:k,2)) == 0 && mean(allsort(2:2:k,2)) == 1
        k = k + 2;
    end
    
    i=i+1;
    
    %we don't want to use a possibly erroneous point in the data and we
    %must end the sequence on a 1, so when the dicontinuity occurs with two
    %1s in a row, we must roll back by 2
    if size(allsort,1) > k && k > 2
        if allsort(k-1,2)==1 && allsort(k-2,2) == 1
            allsort = allsort(1:k-4,:);
        else
            allsort = allsort(1:k-2,:);
        end
    end
    
    % for the special case where there are two discontinuities of 0s in a row
    if k == 2 && allsort(1,2) == 0
        allsort(1:2,:) = [];
        longest_length(i,1) = 0;
        
        %we want to index one passed the discontinuity
        if i == 1
            %if this occurs for the first index, only includes values
            %skipped
            index(1,1) = idx_skip;
            index(2,1) = idx_skip;
        else
            index(1,i) = index(2,i-1)+idx_skip+1;
            index(2,i) = index(2,i-1)+idx_skip+1;
        end
        skip = skip + 2;
        
    else
        
        %otherwise count as normal
        longest_length(i,1) = size(allsort,1);
        
        
        %create ordered index of where continuous chunks occur
        if i == 1
            index(1,1) = 1 + idx_skip;
            index(2,1) = longest_length(i,1)+idx_skip;
        else
            index(1,i) = index(2,i-1)+3+idx_skip;
            
            %Longest_length can only be 0 when
            %two discontinuities of 1s happen in a row, below accounts that
            %the index end needs to still progress by 1 (but longest_length
            %still needs to be 0 for the main counter)
            if longest_length(i,1) >0
                index(2,i) = index(1,i)+longest_length(i,1)-1;
            else
                index(2,i) = index(1,i)+longest_length(i,1);
            end
        end
        
        
        %reset allsort for next loop iteration to be passed the discontinuity
        allsort = allsort_bak(index(2,i)+3:end,:);
        
        
        %however we want to skip passed the discontinuity to the next footfall.
        %This entails skipping the discontinuity (for example if the
        %discontinuity is two FF, we skip over these two values
        
        skip = skip+ 2;
    end
    
end

%% determine which index has the largest continuous block
[~,longest_index]  = max(diff(index));

%reorder allsort to contain only this block
allsort_longest = allsort_bak(index(1,longest_index):index(2,longest_index),:);
allsort = allsort_longest;

% Break back into components
FFi_mod = allsort(allsort(:,2) == 0,1);
FBi_mod = allsort(allsort(:,2) == 1,1);

%want to track frame numbers of when the block on continuous data starts
%and ends
block_start = FFi_mod(1,1);
block_end = FBi_mod(end,1);