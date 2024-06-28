function M = simForward(M)

if M.Nsteps == 0 || isempty(M.tumors) % nothing to simulate here
    return;
end

%% prepare new tracked values
names = fieldnames(M.tracked);
for i = 1:length(names)
    M.tracked.(names{i}) = cat(1,M.tracked.(names{i}),...
        zeros([M.Nsteps,size(M.tracked.(names{i}),2:ndims(M.tracked.(names{i})))]));
end

%% iterations
for i = 1:M.Nsteps

    if mod(i,100) == 1
        disp(['Current time:' num2str(M.t)]);
    end

    %%% added by Shirlyn to stop then too many tumor cells have hit the boundary
    temp_tumorloc = M.tumors(:,1:3) - 26;
    max_tumorloc = max(temp_tumorloc, [], 2);
    num_hitbound= length(find(max_tumorloc == 25));
    %%%

    if M.NT > M.pars.max_tumor_size || num_hitbound > M.pars.max_tumor_boundary || M.NT <= M.pars.min_tumor_size
        names = fieldnames(M.tracked);
        for ni = 1:length(names)
            colons = repmat({':'},1,ndims(M.tracked.(names{ni}))-1);
            M.tracked.(names{ni})(M.i+1:end,colons{:}) = NaN;
        end
        
        if num_hitbound > M.pars.max_tumor_boundary || M.NT > M.pars.max_tumor_size
            M.escaped = 1;
        end

        break
    end

    M.t = M.t + M.dt;
    M.i = M.i+1;

    M = updateFGFR3(M);

    M = updateTumor(M);

    M = updateImmune(M);

    %% clean up tumor stuff
    M = removeApoptotic(M);

    %% advance tumor proliferation timers
    M.tumors(:,M.I.proliferation_timer) = max(0,M.tumors(:,M.I.proliferation_timer)-M.dt); % update time until tumor cells can proliferate for all that did not proliferate (1) or were just created (0) because those timers were updated in updateTumor

    %% recruit immune cells
    M = recruitImmune(M);

    %% update these vals
    M.NT = size(M.tumors,1);
    M.mut_log = M.tumors(:,M.I.tumor_mut)==1;
    M.ha_log = M.tumors(:,M.I.type)==1;
    M.NI = size(M.immunes,1);

    %% update tracked values
    M = updateTracked(M);

    %% plot
    if M.plot_pars.plot_fig && mod(M.i,M.plot_pars.plot_every)==M.plot_pars.plot_offset
        plotFunction_EndStep(M)
    end

    %% save any "big" data
    if M.t >= M.next_save_time - 0.5 * M.dt % then it appears that this is the closest time to the desired save time
        M = saveModelData(M);
    end

end %%end of for