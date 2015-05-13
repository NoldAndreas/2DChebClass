function uwall = WallBC(t,Ind,BCWall)    

    if(isempty(BCWall))
        uwall = zeros(sum(Ind.finite2),1);
        return;
    end

    tau   = BCWall.tau;
    bc    = BCWall.bc;
    M     = sum(Ind.finite2);
    
    if(strcmp(bc,'exp'))        
        uwall = (1-exp(-t/tau))*ones(M,1);
        if(isfield(BCWall,'u_max'))
            uwall = BCWall.u_max*uwall;
        end
    elseif(strcmp(bc,'sinHalf'))
        sinHalf          = sin(pi*t/tau);
        sinHalf(t > tau) = 0;
        uwall = sinHalf*ones(M,1);  
    end
end