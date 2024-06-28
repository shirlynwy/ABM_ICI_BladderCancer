function M = initializeTracked(M)

M.tracked.t = 0;
M.tracked.NT = M.NT;
M.tracked.tumor_types = sum(M.mut_log==[0,1] & M.ha_log==reshape([0,1],1,1,2),1); % dims: [time,mut,antigen]; time increases on time dim, mut status goes [non-mut,mut], ant goes [low,high]
M.tracked.NI = M.NI;

M.tracked.tum_apop = zeros(1,2,2); % tracks apoptosis of all 4 tumor types at each time point; [LA,HA;LA Mut,HA Mut]
M.tracked.tum_prolif = zeros(1,2,2); % tracks proliferation of all 4 tumor types at each time point; [LA,HA;LA Mut,HA Mut]
M.tracked.tum_contact_inhibition = zeros(1,2,2); % tracks contact inhibition of all 4 tumor types at each time point; [LA,HA;LA Mut,HA Mut]

M.tracked.phiD_mean = M.fgfr3.all_receptors_without_inhibitor(3) * any(M.mut_log);
M.tracked.phiD_std = 0;

M.tracked.imm_recruit = 0;
M.tracked.imm_prolif = 0;
M.tracked.imm_contact_inhibition = 0;
M.tracked.deactivations = zeros(1,3); % tracks deactivations of [unengaged, FASL (slow)-killing, perforin (fast)-killing] CTLs
M.tracked.imm_cleared = zeros(1,2,2); % tracks tumor cells killed by ctls of all 4 tumor types at each time point; [LA,HA;LA Mut,HA Mut]
M.tracked.imm_movement_halts = 0;
M.tracked.imm_attack_misses = 0;
M.tracked.imm_attack_hits = 0;
M.tracked.imm_aicd = 0;
