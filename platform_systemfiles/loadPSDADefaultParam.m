function[d,Sadef,Ddef,SMLIB,default_reg]=loadPSDADefaultParam(ME)

d     = logsp(1,100,10); % displacement in CM
Sadef = 30; % default number of Sa simulations
Ddef  = 40; % default number of D simulations

% this part needs improvement to account for PSHA with non-Sa GMMs
% default regular displacement models
SMLIB(1:0) = struct('id',[],'str',[],'func',[],'isregular',[],'integrator',[],'Safactor',[],'param',[]);
cont = 0;

[~,ind1] = intersect({ME.str},'psda_BMT2017M'); % default subduction
[~,ind2] = intersect({ME.str},'psda_BT2007M');  % default subductio
mth      = ME([ind1,ind2]);       % default for interface and slab sources
default_reg.T3= {'slope1','DEF1','DEF1','DEF2',1};
for j=1:2
    cont=cont+1;
    SMLIB(cont).id         = sprintf('DEF%g',cont);
    SMLIB(cont).str        = mth(j).str;
    SMLIB(cont).func       = mth(j).func;
    SMLIB(cont).isregular  = mth(j).isregular;
    SMLIB(cont).integrator = mth(j).integrator;
    SMLIB(cont).primaryIM  = mth(j).primaryIM;
    SMLIB(cont).Safactor   = mth(j).Safactor;
    SMLIB(cont).param      = [];
end

