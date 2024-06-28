function M = saveFinalModelData(M)

load(sprintf("data/%s/%s/output_%08d.mat",M.save_pars.batch_name, M.save_pars.sim_identifier,M.save_index-1),"time");
if M.t>time % make new save
    M = saveModelData(M);
end

tracked = M.tracked;
if isfield(M, 'escaped')
    escaped = M.escaped;
else
    escaped = 0;
end
save(sprintf("data/%s/%s/output_final",M.save_pars.batch_name, M.save_pars.sim_identifier),"tracked", "escaped",'-v7.3')