function M = updateTracked(M)

M.tracked.t(M.i) = M.t;

M.tracked.NT(M.i) = M.NT;
M.tracked.tumor_types(M.i,:,:) = sum(M.mut_log==[0,1] & M.ha_log==reshape([0,1],1,1,2),1); %length of ha_log and mut_log: # of tumor cells

M.tracked.NI(M.i) = M.NI;
