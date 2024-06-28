function simCohort(M,cohort_pars)

nsamps_per_condition = cohort_pars.nsamps_per_condition;

% lattice sampling of fields of M
fn = string(fieldnames(M));
cohort.lattice_parameters = struct("path",{},"values",{});
for i = 1:numel(fn)
    current_struct_path = fn(i);
    cohort.lattice_parameters = grabFields(M.(fn(i)),cohort.lattice_parameters,current_struct_path);
end

% non-lattice parameter varying (for dosing schedule, rectangle sizes...and others?)
if cohort_pars.last_dose_is_no_dose
    dose_start_inds = [];
    vals = {};
    paths = {};
    for i = 1:numel(cohort.lattice_parameters)
        if cohort.lattice_parameters(i).path(end)=="start_day"
            dose_start_inds(end+1) = i;
            vals{end+1} = cohort.lattice_parameters(i).values;
            paths{end+1} = cohort.lattice_parameters(i).path;
        end
    end
    if length(dose_start_inds)>1 % only make this change if two dose start times are varying
        cohort.lattice_parameters(end+1).values = allCombos(vals{:},'matlab');
        cohort.lattice_parameters(end).values(end,:) = Inf;
        cohort.lattice_parameters(end).path = paths;
        cohort.lattice_parameters(dose_start_inds) = [];
    end
end

% cohort = nameSubCohorts(cohort); saves names of the subcohorts in .mat

cohort_size = arrayfun(@(i) size(cohort.lattice_parameters(i).values,1),1:numel(cohort.lattice_parameters));
total_runs = prod(cohort_size) * nsamps_per_condition;
cohort_pars.total_runs = total_runs;

colons = repmat({':'},[1,length(cohort_size)]);
vp_ind = cell(1,length(cohort_size));

cohort_start_time = string(datetime("now","Format","yyMMddHHmmssSSS"));
used_sims = strings(0,1); % track which previously ran sims are being reused here

fprintf("Beginning simulations now...\n")

ndigits = ceil(log10(total_runs+1)); % get enough digits after the underscore (+1 because I start counting from 1, not 0)
sim_id_format = sprintf("_%%0%dd",ndigits);

if total_runs>=cohort_pars.min_parfor_num
    F(1:total_runs) = parallel.FevalFuture;
    ppool = gcp;
    cohort_pars.num_workers = ppool.NumWorkers;
    for ri = 1:total_runs % run index
        [vp_ind{colons{:}},~] = ind2sub([cohort_size,nsamps_per_condition],ri);
        for vpi = 1:numel(cohort.lattice_parameters)
            M = setField(M,cohort.lattice_parameters(vpi).path,cohort.lattice_parameters(vpi).values(vp_ind{vpi},:));
        end
        [sim_this,used_sims] = findSimilarSims(M,used_sims,cohort_start_time);
        if sim_this
            M.save_pars.sim_identifier = string(datetime("now","Format","yyMMddHHmm")) + sprintf(sim_id_format,ri);
            F(ri) = parfeval(ppool,@simPatient,1,M);
        else
            F(ri) = parfeval(ppool,@grabSimData,1,used_sims(end));
        end
    end
else
    cohort_pars.num_workers = 1;
end

cohort_pars.mu_n = 0;
cohort_pars.start = tic;
cohort_pars.batch_start = tic;

for ri = total_runs:-1:1
    if total_runs>=cohort_pars.min_parfor_num
        [idx,out_temp] = fetchNext(F);
    else
        idx = ri;
        [vp_ind{colons{:}},~] = ind2sub([cohort_size,nsamps_per_condition],ri);
        for vpi = 1:numel(cohort.lattice_parameters)
            M = setField(M,cohort.lattice_parameters(vpi).path, cohort.lattice_parameters(vpi).values(vp_ind{vpi},:));
        end
        [sim_this,used_sims] = findSimilarSims(M,used_sims,cohort_start_time);
        if sim_this
%             while sim_this
%                 [sim_this,used_sims] = findSimilarSims(M,used_sims,cohort_start_time);
%             end
%                         out_temp = grabSimData(used_sims(end));
            M.save_pars.sim_identifier = string(datetime("now","Format","yyMMddHHmm")) + sprintf(sim_id_format,ri);
            out_temp = simPatient(M);
        else
            out_temp = grabSimData(used_sims(end));
        end
    end
    cohort_pars = updateCohortTimer(cohort_pars,total_runs-ri+1);
    cohort.tracked(idx) = out_temp.tracked;
    if isfield(out_temp.save_pars,"sim_identifier")
        cohort.ids(idx) = out_temp.save_pars.sim_identifier;
    end
end

cohort.tracked = reshape(cohort.tracked,[cohort_size,nsamps_per_condition,1]);
if M.save_pars.dt < Inf
    cohort.ids = reshape(cohort.ids,[cohort_size,nsamps_per_condition,1]);
end

if ~isfield(cohort_pars,"cohort_identifier")
    cohort_pars.cohort_identifier = string(datetime("now","Format","yyMMddHHmmssSSS")); % default to this for determining an id if none given
end

while exist(sprintf("data/cohort_%s",cohort_pars.cohort_identifier),"dir") % just in case this directory already exists somehow (not sure how to processes could start at the same time to the millisecond and then one create this folder before the other looks for it)
    cohort_pars.cohort_identifier = string(datetime("now","Format","yyMMddHHmmssSSS")); % default to this for determining an id if none given
end

mkdir(sprintf("data/cohort_%s",cohort_pars.cohort_identifier))

% [subcohort_title_names,subcohort_folder_names] = nameSubCohorts(cohort);

save(sprintf("data/cohort_%s/output",cohort_pars.cohort_identifier),"nsamps_per_condition","total_runs")
save(sprintf("data/cohort_%s/output",cohort_pars.cohort_identifier),'-struct',"cohort","-append")

fprintf("Finished cohort_%s.\n",cohort_pars.cohort_identifier)

end

function lattice_parameters = grabFields(S,lattice_parameters,incoming_struct_path)

fn = string(fieldnames(S));
for i = 1:numel(fn)
    current_struct_path = [incoming_struct_path,fn(i)];
    if isstruct(S.(fn(i)))
        lattice_parameters = grabFields(S.(fn(i)),lattice_parameters,current_struct_path);
    elseif numel(S.(fn(i)))>1 % then vary over these parameters
        lattice_parameters(end+1) = struct("path",current_struct_path,"values",S.(fn(i)));
    end
end

end

function S = setField(S,path,val)

if iscell(path)
    for i = 1:numel(path)
        S = setField(S,path{i},val(i));
    end
elseif length(path)>1
    S.(path(1)) = setField(S.(path(1)),path(2:end),val);
else
    S.(path(1)) = val;
end

end

function [sim_this,used_sims] = findSimilarSims(M,used_sims,cohort_start_time)
sim_this = true;
old_sims = dir("data/sims/*");
old_sims = old_sims([old_sims.isdir]);
for i = numel(old_sims):-1:1
    if ~startsWith(old_sims(i).name,digitsPattern(1)) ... % I only want folders that start with a number 
            || (startsWith(old_sims(i).name,"22" + digitsPattern(13)) && str2double(old_sims(i).name) >= str2double(cohort_start_time)) ... % for those that follow the current time numbering, disregard those that started after this cohort began
            || any(strcmp(old_sims(i).name,used_sims)) % make sure this one has not already been used
        old_sims(i) = [];
    end
end
for i = 1:numel(old_sims) % look for the first one that matches the inputs here
    if exist(sprintf("data/sims/%s/output_constants.mat",old_sims(i).name),"file") && exist(sprintf("data/sims/%s/output_final.mat",old_sims(i).name),"file")
        X = load(sprintf("data/sims/%s/output_constants.mat",old_sims(i).name));
        fn = fieldnames(M);
        these_match = true;
        for j = 1:numel(fn)
            if (~strcmp(fn{j},"plot_pars")) % don't worry about plot_pars being equal
                par_fn = fieldnames(M.(fn{j}));
                for k = 1:numel(par_fn)
                    if strcmp(par_fn{k},'n_regions') % n_regions is actually set for each substrate (i.e. I can remove n_regions from base parameter stuff)
                        continue;
                    end
                    if strcmp(par_fn{k},'deactivation_function') % it seems the isequal cannot really check if two anonymous functions are equal
                        continue;
                    end
                    if ~isfield(X,fn{j}) || ~isfield(X.(fn{j}),par_fn{k}) || ~isequal(M.(fn{j}).(par_fn{k}),X.(fn{j}).(par_fn{k}))
%                         disp(fn{j})
%                         disp(par_fn{k})
%                         if strcmp(old_sims(i).name,"221012070214082")
%                             disp('')
%                         end
                        these_match = false;
                        break;
                    end
                end
                if ~these_match
                    break;
                end
            end
        end
        if these_match
            sim_this = false;
            used_sims(end+1) = string(old_sims(i).name);
            return;
        end
    end
end
end

function out = grabSimData(folder_name)

load(sprintf("data/sims/%s/output_final.mat",folder_name),"tracked")
out.tracked = tracked;
out.save_pars.sim_identifier = folder_name;

end

% function [subcohort_title_names,subcohort_folder_names] = nameSubCohorts(cohort)
% 
%     subcohort_title_names = string(1:numel(cohort.ids)/size(cohort.ids,ndims(cohort.ids)));
%     subcohort_folder_names = string(1:numel(cohort.ids)/size(cohort.ids,ndims(cohort.ids)));
% 
% 
% end
