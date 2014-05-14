function [VBack_S,VAdd_S]=APSG(y,t,optsPhys)


    Vm      = optsPhys.Vm;
    tSwitch = optsPhys.tSwitch;

    %--------------------------------------------------------------------------
    VBack  = Vm.*y.^2;
    DVBack = 2*Vm.*y;
    %--------------------------------------------------------------------------

    if(t==0 || t>tSwitch(1))
    %if(t>tSwitch)
        Z=40;
        strength = 10;
    else
        Z = 8;
        strength = 10;
    end
    
    VAdd  = -strength.*exp(-y.^2/2/Z);
    DVAdd = -y/Z.*VAdd;
        
    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack,'DV',DVBack);

    VAdd_S = struct('V',VAdd,'DV',DVAdd);

end