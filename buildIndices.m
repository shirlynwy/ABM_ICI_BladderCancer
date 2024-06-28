function I = buildIndices()

I.subs = 1:3;
I.ind = 4;
I.type = 5; % for tumor cells: 0=low antigen; 1=high antigen...for immune cells: -1=deactivated; 0=unengaged; 1=only FasL (slow); 2=perforin + FasL (fast)
I.event = 6;
I.proliferation_timer = 7;
I.region_fgfr3 = 8;
I.region_checkpoint = 9;

I.tumor_clearance = I.region_checkpoint + (1:2); % [amount this tumor cell has been cleared as a probability, rate at which it is currently being cleared]
I.tumor_mut = I.region_checkpoint + 3;
I.tumor_receptors = I.region_checkpoint + (4:9); % column indices in which the receptor information is stored (always have these be the last)

I.immune_target = I.region_checkpoint + 1;
I.immune_seek_time = I.immune_target + 1; % tracks time since last conjugation for this ctl
I.immune_gradient = I.region_checkpoint + (3:5); % always have this be the last (not sure why I said this anymore...some older versions may have used this in intializing the M.immune array with the right number of columns)