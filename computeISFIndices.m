function [xx,yy,zz,xind,yind,zind] = computeISFIndices(M,subs)


xx = subs(1) + M.immune_stimulatory_factor.interval;
xind = xx>=1 & xx<=M.grid.size(1);
xx = xx(xind);

yy = subs(2) + M.immune_stimulatory_factor.interval;
yind = yy>=1 & yy<=M.grid.size(2);
yy = yy(yind);

zz = subs(3) + M.immune_stimulatory_factor.interval;
zind = zz>=1 & zz<=M.grid.size(3);
zz = zz(zind);