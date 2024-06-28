function M = simPatient(M)

M = finishParameterSetup_Patient(M);

%% initialize inputs
M.t = 0;
M.i = 1; % step index
M = initializeEnvironment(M);
M = initializeTracked(M);

%% initialize figure
if M.plot_pars.plot_fig
    M = initializeFigure(M);
end

%%
if M.save_pars.dt < Inf
    M = finishParameterSetup_Saves(M);
    M = saveInitialModelData(M);
else
    M.next_save_time = Inf;
end

for ei = 1:M.events.n % event index
    M = finishParameterSetup_Event(M,M.events.times(ei)-M.t);
    M = simForward(M);

    extendTimeSeriesPlots(M,ei)

    switch M.events.event_index(ei)
        case 0 % begin immune recruitment
            error('have not got this yet')

        case 1 % add dose of anti-PD1
            M.checkpoint.aPD1_circulation = M.checkpoint.aPD1_circulation + M.checkpoint.dose_vals(M.checkpoint.dose_index);
            M.checkpoint.dose_index = M.checkpoint.dose_index + 1;

        case 2 % add dose of anti-FGFR3
            M.fgfr3.circ = M.fgfr3.circ + M.fgfr3.dose_vals(M.fgfr3.dose_index);
            M.fgfr3.dose_index = M.fgfr3.dose_index + 1;

        case Inf % end of simulation
            break;
    end
end

if M.save_pars.dt < Inf
    M = saveFinalModelData(M);
end
