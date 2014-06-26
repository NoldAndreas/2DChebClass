function TestAccuracyRectangle(abox,PhysArea,Conv,optsPhys,mu_sat,rhoLiq_sat,HS_f)%Pts,Maps,Interp
%TestAccuracyRectangle(PhysArea,Conv,Pts,Maps,Interp)
        wd = abox.y1Max - abox.y1Min;
        he = abox.y2Max - abox.y2Min;
        
        if(wd == inf || he == inf)
            return;
        end
        if(~strcmp(optsPhys.V2.V2DV2,'Phi2DLongRange'))
            return;
        end
        
        %********************** 1st Test: *******************
        val = Conv*ones(size(abox.Pts.y1_kv));
        
        Interp0   = abox.InterpolationMatrix_Pointwise(PhysArea.y1Min + wd/2,PhysArea.y2Min+he/2);       %Interp0   = SpectralSpectral_Interpolation(PhysArea.y1Min + wd/2,PhysArea.y2Min+he/2,Pts,Maps); 
        res       = Interp0*val-Rectangle(wd/2,he/2,1);
        disp(['1. Accuracy: ',num2str(res),' (for Convolution of Long-Range-Potential at origin).']);        
        
        if(nargin == 6)
            return;
        end
                
        %********************** 2nd Test: *******************
        optsPhysTest           = optsPhys;
        kBT                    = optsPhys.kBT;
        optsPhysTest.epsilon_w = rhoLiq_sat*[1;1;1;1];
        optsPhysTest.V0        = 0;
        
        Vext   = Vext_Cart_Capillary_Static(abox.Pts.y1_kv,abox.Pts.y2_kv,0,optsPhysTest);
        
        rho_ig = rhoLiq_sat*ones(size(abox.Pts.y1_kv));
        x_ig   = kBT*log(rho_ig);
        
        opts = optimset('Display','off','Jacobian','on');
        x      = fsolve(@f,x_ig,opts);
        
        rho_dev  = exp(x/kBT)-rhoLiq_sat;
        n      = max(abs(rho_dev));
        disp(['2. Accuracy: ',num2str(n),' (for equilibrium liquid density)']);
        %doPlots_SC(Interp,Pts,rho_dev); title('Numerical error for computation of equilibrium liquid density.');                
        
        function [mu_s,J] = f(x)
            rho            = exp(x/kBT);
            [muSC,~,dmuSC] = HS_f(rho,kBT);
                        
            mu_s  = x + Conv*rho - mu_sat + muSC + Vext;                        
            J     = diag(1 + dmuSC.*rho/kBT) + Conv*diag(rho)/kBT;
        end
        
end
