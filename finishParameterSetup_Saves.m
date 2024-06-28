function M = finishParameterSetup_Saves(M)

M.next_save_time = 0;
M.save_index = 0;

if ~isfield(M.save_pars,"sim_identifier")
    M.save_pars.sim_identifier = string(datetime("now","Format","yyMMddHHmmssSSS")); % default to this for determining an id if none given
end

while exist(sprintf("data/%s",M.save_pars.sim_identifier),"dir") % just in case this directory already exists somehow (not sure how to processes could start at the same time to the millisecond and then one create this folder before the other looks for it)
    M.save_pars.sim_identifier = string(datetime("now","Format","yyMMddHHmmssSSS")); % default to this for determining an id if none given
end

mkdir(sprintf("data/%s/%s",M.save_pars.batch_name, M.save_pars.sim_identifier))

max_grid_sub = max(M.grid.size);
max_grid_sub_log2 = ceil(log2(max_grid_sub+1));
if max_grid_sub_log2 <= 8
    M.save_pars.integrify = @uint8;
elseif max_grid_sub_log2 <= 16
    M.save_pars.integrify = @uint16;
elseif max_grid_sub_log2 <= 32
    M.save_pars.integrify = @uint32;
elseif max_grid_sub_log2 <= 64
    M.save_pars.integrify = @uint64;
else
    error("grid too large to store!")
end
