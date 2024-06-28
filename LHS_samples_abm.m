close all
clear all

runs= 2000;

%% LHS parameter ranges
prop_ha0= LHS_Call(0.05, 0.5, 0.95, 0, runs, 'unif');
immune_recruit_rate2 = LHS_Call(2.5, 10, 25, 0, runs, 'unif');
isf_prolif_max = LHS_Call(0.04, 0.15, 0.6, 0, runs, 'unif');
p1 = LHS_Call(0.01, 0.92, 1, 0, runs, 'unif');
p2 = LHS_Call(0.01, 0.33, 1, 0, runs, 'unif');
fast_kill_rate = LHS_Call(12, 48, 120, 0, runs, 'unif');
move_rate_microns = LHS_Call(1*(24*60), 2 * (24*60), 8*(24*60) , 0, runs, 'unif'); % immune movement rate in microns / day

conjugation_rate = LHS_Call(12, 28.8, 96, 0, runs, 'unif'); %immune_pars
la_isf_factor = LHS_Call(0.1, 0.5, 0.9, 0, runs, 'unif'); %immune_pars; 
isf_reach = round(LHS_Call(3, 5, 10, 0, runs, 'unif')); % immune_stimulatory_factor_pars
% isf_length_scale = LHS_Call(5, 20, 30, 0, runs, 'unif'); 
% imm_steps_per_move = round(LHS_Call(2, 4, 12, 0, runs, 'unif')); %immune_pars; must be integers
% NI0 = round(LHS_Call(0, 1, 6, 0, runs, 'unif')); %must be integers

LHSmatrix = [prop_ha0  immune_recruit_rate2  isf_prolif_max  fast_kill_rate...
             p1  p2  move_rate_microns ...
             conjugation_rate  la_isf_factor  isf_reach];

var_paths ={["setup",'prop_ha0' ], ["immune_pars", "immune_recruit_rate2"],...
           ["immune_pars", "isf_prolif_max"],["immune_pars", "fast_kill_rate"] ...
           ["immune_pars", "p1"], ["immune_pars", "p2"],...
           ["immune_pars", "move_rate_microns"],...
           ["immune_pars", "conjugation_rate"], ["immune_pars", "low_antigen_isf_factor"],...
           ["immune_stimulatory_factor_pars", "reach"]} ;
% var_names = {'prop_ha0', 'immune_recruit_rate2', 'isf_prolif_max',...
%              'p1', 'p2', 'fast_kill_rate'};

lhs_parameters = struct("path",{},"values",{});
for i = 1:numel(var_paths)
    curcell = var_paths{i};
    curval = LHSmatrix(:,i);
    lhs_parameters(end+1) = struct("path", curcell, "values", curval);
end


% figure;
% for i = 1:size(LHSmatrix,2)
%     histogram(LHSmatrix(:,i))
%     keyboard;
% end
tempname = ['LHSsamples_abm_' , datestr(now,'mm-dd-yy'), '_%03d.mat'];
file_num = 1;
while exist(sprintf(tempname,file_num),'file')
    file_num = file_num+1;
end
mat_filename =  sprintf(tempname,file_num);
save(mat_filename, 'lhs_parameters');

% 
% mat_filename = ['LHSsamples_cpblk_' , datestr(now,'mm-dd-yy'), '.mat'];
% save(mat_filename, 'LHSmatrix', 'PRCC_var');

