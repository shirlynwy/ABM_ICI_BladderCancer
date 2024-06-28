function [M,event_out] = updateTumor(M)

[M,event_out] = eventSelection_Tumor(M);

M = performEvents(M,event_out);
