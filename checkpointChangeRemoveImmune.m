function M = checkpointChangeRemoveImmune(M,j)

region_losses = accumarray(M.immunes(j,M.I.region_checkpoint),1,[M.checkpoint.n_regions,1]);

M.checkpoint.aPD1 = (M.checkpoint.aPD1 .* M.checkpoint.volumes.regions + region_losses .* M.checkpoint.PD1aPD1)./M.checkpoint.volumes.regions;
M.checkpoint.volumes.pd1 = M.checkpoint.volumes.pd1 - region_losses;
