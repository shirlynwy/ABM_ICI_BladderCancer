function simBatchLHS(M,cohort_pars,I, filename)

nsamps_per_condition = cohort_pars.nsamps_per_condition;
load(filename, 'lhs_parameters');
cohort.lhs_parameters = lhs_parameters;

cohort_size = size(cohort.lhs_parameters(1).values,1);
%%%%%%%
% % lattice sampling of fields of M
% fn = string(fieldnames(M))
% cohort.lattice_parameters = struct("path",{},"values",{});
% 
% 
% for i = 1:numel(fn)
%     current_struct_path = fn(i);
%     cohort.lattice_parameters = grabFields(M.(fn(i)),cohort.lattice_parameters,current_struct_path);
% end
%%%%%%%%%%%%

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
    cohort.lattice_parameters(end+1).values = allCombos(vals{:},'matlab');
    cohort.lattice_parameters(end).values(end,:) = Inf;
    cohort.lattice_parameters(end).path = paths;
    cohort.lattice_parameters(dose_start_inds) = [];
end

%%%%%%
% cohort_size = arrayfun(@(i) size(cohort.lattice_parameters(i).values,1),1:numel(cohort.lattice_parameters));
% 
% colons = repmat({':'},[1,length(cohort_size)]);
% vp_ind = cell(1,length(cohort_size));
% 
% [vp_ind{colons{:}},~] = ind2sub([cohort_size,nsamps_per_condition],I);
% for vpi = 1:numel(cohort.lattice_parameters)
%     M = setField(M,cohort.lattice_parameters(vpi).path,cohort.lattice_parameters(vpi).values(vp_ind{vpi},:));
% end
%%%%%

[vp_ind,~] = ind2sub([cohort_size,nsamps_per_condition],I);
for vpi = 1:numel(cohort.lhs_parameters)
    M = setField(M,cohort.lhs_parameters(vpi).path,cohort.lhs_parameters(vpi).values(vp_ind,:));
end

total_runs = prod(cohort_size) * nsamps_per_condition;
ndigits = ceil(log10(total_runs+1)); % get enough digits after the underscore (+1 because I start counting from 1, not 0)
sim_id_format = sprintf("_%%0%dd",ndigits);
M.save_pars.sim_identifier = string(datetime("now","Format","yyMMddHHmm")) + sprintf(sim_id_format,I);
M = simPatient(M);

disp(['last tumor cell # is ' num2str(M.tracked.NT(end))])
disp(['max tumor cell # is ' num2str(max(M.tracked.NT))])


disp(['last immune cell # is ' num2str(M.tracked.NI(end))])
disp(['max immune cell # is ' num2str(max(M.tracked.NI))])

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
