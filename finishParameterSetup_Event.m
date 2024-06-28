function M = finishParameterSetup_Event(M,DT)

if DT < 1e-10
    DT = 0;
end

M.Nsteps = ceil(DT/M.pars.max_dt);
M.dt = DT/M.Nsteps;
M.Nsteps_Imm = ceil(M.dt/M.pars.max_dt_Imm); % number of immune steps per one tumor step; last one is when immune cells can proliferate. immune timescale = dt/ImSpeed (in days)
M.dt_Imm = M.dt/M.Nsteps_Imm; % actual duration of one immune step (in days)