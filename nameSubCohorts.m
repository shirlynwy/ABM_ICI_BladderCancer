function cohort = nameSubCohorts(cohort)


cohort.subcohort_title_names = string(1:numel(cohort.ids)/size(cohort.ids,ndims(cohort.ids)));
cohort.subcohort_folder_names = string(1:numel(cohort.ids)/size(cohort.ids,ndims(cohort.ids)));

for i = 1:length(cohort.lattice_parameters)
    if iscell(cohort.lattice_parameters(i).path)
        cohort.dim_names(i) = "";
        for j = 1:numel(cohort.lattice_parameters(i).path)
            [dim_name_temp,cohort.dim_labels{i}(:,j)] = nameThis(cohort.lattice_parameters(i).path{j},cohort.lattice_parameters(i).values(:,j));
            if j==1
                cohort.dim_names(i) = dim_name_temp;
            else
                if numel(cohort.lattice_parameters(i).path)>2
                    cohort.dim_names(i) = cohort.dim_names(i) + ", ";
                else 
                    cohort.dim_names(i) = cohort.dim_names(i) + " ";
                end
                if j==numel(cohort.lattice_parameters(i).path)
                    cohort.dim_names(i) = cohort.dim_names(i) + "and ";
                end
                cohort.dim_names(i) = cohort.dim_names(i) + dim_name_temp;
            end
        end
    else
        [cohort.dim_names(i),cohort.dim_labels{i}] = nameThis(cohort.lattice_parameters(i).path,cohort.lattice_parameters(i).values);
    end
end

end

function [dim_name,dim_labels] = nameThis(path,values)
switch path(end)
    case "start_day"
        if path(end-1) == "fgfr3"
            dim_name = "FGFR3 Start Day";
        elseif path(end-1) == "checkpoint"
            dim_name = "ICI Start Day";
        else
            error("Start day of what?")
        end
        dim_labels = "Day " + num2str(values);
        if any(values==Inf)
            dim_labels(values==Inf) = "Control";
        end

    case "prop_ha0"
        dim_name = "Initial High Antigen Proportion";
        dim_labels = num2str(round(100*values)) + "%";

    case "prop_mut0"
        dim_name = "Initial FGFR3 Mutant Proportion";
        dim_labels = num2str(round(100*values)) + "%";

    case "fgfr3_affects_immune_recruit"
        if ~isequal(size(values),[2,1])
            error("This value must be a 2x1")
        end
        dim_name = "FGFR3 Recruitment Effect";
        dim_labels = ["Off";"On"];
        if values(1)
            dim_labels = flip(dim_labels);
        end

    case "fgfr3_affects_cytotoxicity"
        if ~isequal(size(values),[2,1])
            error("This value must be a 2x1")
        end
        dim_name = "FGFR3 Cytotoxicity Effect";
        dim_labels = ["Off";"On"];
        if values(1)
            dim_labels = flip(dim_labels);
        end

    otherwise
        error("Not sure what this parameter is")
end
end