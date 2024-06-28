function M = initializeEnvironment(M)

M = initializeGrid(M);
M = initializeBloodVessels(M);
M = initializeFGFR3(M);

M = initializeTumor(M);
M = initializeImmune(M);

M = initializeCheckpoint(M);

if M.setup.NI0>0
    M = placeImmune(M,M.setup.NI0);
    M.NI = size(M.immunes,1);
end

M = initializeImmuneStimulatoryFactor(M);

M.L = full(sparse(M.tumors(:,M.I.ind),1,M.val.tum(M.tumors(:,M.I.tumor_mut)+1),M.V_tot,1)); % this tracks the occupancy of each lattice site


