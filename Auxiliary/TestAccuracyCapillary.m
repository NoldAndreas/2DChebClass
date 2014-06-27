function TestAccuracyCapillary(PhysArea,Conv,Pts,Maps,Interp,Ind,optsPhys,mu_sat,rhoLiq_sat,HS_f)
%TestAccuracyRectangle(PhysArea,Conv,Pts,Maps,Interp)        
        he = PhysArea.y2Max - PhysArea.y2Min;
        
        if( he == inf)
            return;
        end
        if(~strcmp(optsPhys.V2DV2,'Phi2DLongRange'))
            return;
        end
        
        %********************** 1st Test: *******************
        val = Conv*ones(size(Pts.y1_kv));
        %Interp0   = SpectralSpectral_Interpolation(0,PhysArea.y2Min+he/2,Pts,Maps); 
        %res = Interp0.InterPol*val-Capillary(he/2,1);
        c  = min(abs(Pts.y2_kv));
        ma = (abs(Pts.y2_kv) <= c) & (abs(Pts.y1_kv) ~= Inf);
        res = max(abs(val(ma)-Capillary(he/2,1)));
        disp(['1. Error: ',num2str(res),' (for Convolution of Long-Range-Potential at origin).']);        
                
        %********************** 2nd Test: *******************
        optsPhysTest           = optsPhys;
        kBT                    = optsPhys.kBT;
        optsPhysTest.epsilon_w = rhoLiq_sat*[1;1;1;1];
        optsPhysTest.V0        = 0;

        y2Inner = Pts.y2_kv(~Ind.left & ~Ind.right);
        
        Vext   = Vext_Cart_Slit_Static([],Pts.y2_kv,[],optsPhysTest);
        
        rho_ig = rhoLiq_sat*ones(size(y2Inner));
        rhoInf = rhoLiq_sat*ones(size(Pts.y1_kv(Ind.left)));        
        
        opts = optimset('Display','off','Jacobian','on');
        rho  = fsolve(@f,rho_ig,opts);
        
        rho_dev  = GetFullRho(rho)-rhoLiq_sat;
        n      = max(abs(rho_dev));
        disp(['2. Error: ',num2str(n),' (for equilibrium liquid density)']);
        doPlots_SC(Interp,Pts,rho_dev); title('Numerical error for computation of equilibrium liquid density.');                
        
        function [mu_s,J] = f(rho)            
            rho = GetFullRho(rho);
            [muSC,~,dmuSC] = HS_f(rho,kBT);
                        
            mu_s  = kBT*log(rho) + muSC + Conv*rho - mu_sat + Vext; 
            J     = diag(kBT./rho + dmuSC) + Conv;

            mu_s  = mu_s(~Ind.left & ~Ind.right);
            J     = J(~Ind.left & ~Ind.right,~Ind.left & ~Ind.right);
        end

    function rhof = GetFullRho(rho)
        rhof = zeros(size(Pts.y1_kv));
        rhof(~Ind.left & ~Ind.right) = rho;
        rhof(Ind.left)  = rhoInf;
        rhof(Ind.right) = rhoInf;
    end
        
end
