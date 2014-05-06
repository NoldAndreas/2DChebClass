function [VBack_S,VAdd_S]=APSHS(R,t,optsPhys)


% get potential parameters
Vm=optsPhys.Vm;

Z=optsPhys.Z;

% reshape R to separate species;
% note that for stochastic, the number of species will always be 1 as we
% made x a column vector so this does nothing
nSpecies=size(R,2);
R=reshape(R,[],nSpecies);

b=0;
a=50;
rI=2;

if(t==0)
    rt=4;
else
    rt=3;
end


rIh=4;

h       =  (erf((R+rt)/rIh) - erf((R-rt)/rIh))/2;
Dh      =  1/sqrt(pi)/rIh*( exp(-(R+rt).^2/rIh^2) - exp(-(R-rt).^2/rIh^2) );
VAdd    =  Vm.*( h.*(b  -R.^2) - a*Z.*exp(-(R-rt).^2 /rI^2 ) ) ;
VAdd(abs(R)==inf) = 0;
DVAdd   =  Vm.*( Dh.*(b  -R.^2) + h.*(-2*R) ...
                + a* Z.*exp(- ( R-rt).^2 /rI^2 ).*2.*(R-rt)/rI^2 );
DVAdd(abs(R)==inf) = 0;
VBack=Vm.*R.^2;
DVBack=2*Vm.*R;


VBack_S = struct('V',VBack, ...
                    'dy',DVBack, ...
                    'grad', DVBack);

               
VAdd_S = struct('V',VAdd, ...
                    'dy',DVAdd, ...
                    'grad', DVAdd);

end