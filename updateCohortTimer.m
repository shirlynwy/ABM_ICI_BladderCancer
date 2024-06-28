function cohort = updateCohortTimer(cohort,i)

if mod(i,cohort.num_workers)==0
    v_n = toc(cohort.batch_start);
    cohort.mu_n = (((i/cohort.num_workers)-1)*cohort.mu_n+2*v_n)/((i/cohort.num_workers)+1); % this is computing the average duration of each run, weighting the more recent runs more heavily
    etr = cohort.mu_n*(cohort.total_runs-i)/cohort.num_workers;
    fprintf('Finished %d of %d, or %3.2f%%, after %s. ETR: %s for total run time of %s.\n',...
        i,cohort.total_runs,100*i/cohort.total_runs,duration(0,0,v_n),duration(0,0,etr),duration(0,0,etr+toc(cohort.start)))
    cohort.batch_start = tic;
end
