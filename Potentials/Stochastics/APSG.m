function [VBack_S,VAdd_S]=APSG(R,t,optsPhys)

% get potential parameters
Vm=optsPhys.Vm;

tSwitch = optsPhys.tSwitch;

% reshape R to separate species;
% note that for stochastic, the number of species will always be 1 as we
% made x a column vector so this does nothing
nSpecies=size(R,2);
R=reshape(R,[],nSpecies);

VBack=Vm.*R.^2;
DVBack=2*Vm.*R;

tSwitch=tSwitch(1);

if(t==0  || t>tSwitch)
    Z=40;
    strength=10;
else
    Z=8;
    strength=10;
end

VAdd=-strength.*exp(-R.^2/2./Z);
DVAdd=-R./Z.*VAdd;

VBack_S = struct('V',VBack, ...
                    'dy',DVBack, ...
                    'grad', DVBack);

               
VAdd_S = struct('V',VAdd, ...
                    'dy',DVAdd, ...
                    'grad', DVAdd);

end