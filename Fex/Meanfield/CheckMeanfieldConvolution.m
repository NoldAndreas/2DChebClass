function res = CheckMeanfieldConvolution(this)
    %Test Convolution at infinity      
    global PersonalUserOutput
    
     M         = this.IDC.M;              
     [h1,h2,a] = getV2(0,this.optsPhys.V2);     
     
     %Convolution profile
     if(isa(this.IDC,'HalfSpace'))
         
        marky2Inf = (this.IDC.Pts.y2_kv == inf);
        y2MinCart = min(this.IDC.GetCartPts.y2_kv);
        PrintErrorPos(this.IntMatrV2.Conv(marky2Inf,:)*ones(M,1)- 2*a,'convolution at y2 = infinity',this.IDC.Pts.y1_kv(marky2Inf));          
         
         if(strcmp(this.optsPhys.V2.V2DV2,'Phi2DLongRange'))
             y0R = this.IDC.GetCartPts.y2_kv - y2MinCart;
             h   = conv_Phi2DLongRange(y0R);             
             [res.error_conv1,ind_conv1] = PrintErrorPos(h-this.IntMatrV2.Conv*ones(M,1),'Phi2DLongRange*1',this.IDC.GetCartPts);
             res.error_conv1_posy2 = this.IDC.GetCartPts.y2_kv(ind_conv1);
         elseif(strcmp(this.optsPhys.V2.V2DV2,'BarkerHenderson_2D'))                 
             conv  = this.IntMatrV2.Conv(this.IDC.Pts.y1_kv==inf,:);
             y2_h  = this.IDC.GetCartPts.y2_kv(this.IDC.Pts.y1_kv==inf) - y2MinCart;
             check = conv_BarkerHenderson2D(y2_h);   

             [res.error_conv1,ind_conv1] = PrintErrorPos(conv*ones(M,1) - check,'convolution at y1 = infinity',y2_h);                 
             res.error_conv1_posy2 = y2_h(ind_conv1);         
         elseif(strcmp(this.optsPhys.V2.V2DV2,'BarkerHendersonHardCutoff_2D'))    
             res.error_conv1 = inf;
             res.error_conv1_posy2 = 1;             
         elseif(strcmp(this.optsPhys.V2.V2DV2,'ExponentialDouble'))    
             conv  = this.IntMatrV2.Conv(this.IDC.Pts.y1_kv==inf,:);
             y2_h  = this.IDC.GetCartPts.y2_kv(this.IDC.Pts.y1_kv==inf) - y2MinCart;
             check = conv_ExpDouble2D(y2_h);

             [res.error_conv1,ind_conv1] = PrintErrorPos(conv*ones(M,1) - check,'convolution at y1 = infinity',y2_h);                 
             res.error_conv1_posy2 = y2_h(ind_conv1);             
             %res.error_conv1 = inf;
             %res.error_conv1_posy2 = 1;
         elseif(strcmp(this.optsPhys.V2.V2DV2,'BarkerHendersonCutoff_2D'))                 
             conv  = this.IntMatrV2.Conv(this.IDC.Pts.y1_kv==inf,:);
             y2_h  = this.IDC.GetCartPts.y2_kv(this.IDC.Pts.y1_kv==inf) - y2MinCart;
             check = conv_BarkerHendersonCutoff2D(y2_h);   

             [res.error_conv1,ind_conv1] = PrintErrorPos(conv*ones(M,1) - check,'convolution at y1 = infinity',y2_h);                 
             res.error_conv1_posy2 = y2_h(ind_conv1);         
         elseif(strcmp(this.optsPhys.V2.V2DV2,'ConstShortRange'))                 
             conv  = this.IntMatrV2.Conv(this.IDC.Pts.y1_kv==inf,:);
             y2_h  = this.IDC.GetCartPts.y2_kv(this.IDC.Pts.y1_kv==inf) - y2MinCart;
             check = conv_ConstShortRange(y2_h);   
             [res.error_conv1,ind_conv1] = PrintErrorPos(conv*ones(M,1) - check,'convolution at y1 = infinity',y2_h);                 
             res.error_conv1_posy2 = y2_h(ind_conv1);                      
         else                          
             cprintf('m','CheckMeanfieldConvolution: Case not yet implemented\n');
         end  
     elseif(isa(this.IDC,'InfCapillary'))
         if(PersonalUserOutput)
             figure('Position',[0 0 800 500]);
             title('Convolution of Barker-Henderson potential with rho=1');
             marky1Inf = (this.IDC.Pts.y1_kv == inf); 
             convOne   = this.IntMatrV2.Conv*ones(M,1);
             this.IDC.plotLine([0 0],[this.IDC.y2Min this.IDC.y2Max],convOne); hold on;
             
             %analytical comparison
             y2_h  = this.IDC.GetCartPts.y2_kv(this.IDC.Pts.y1_kv==inf) - this.IDC.y2Min;
             check  = conv_BarkerHenderson_IC(y2_h);
             
            if(strcmp(this.optsPhys.V2.V2DV2,'BarkerHenderson_2D'))                
                plot(y2_h,check,'g','linewidth',2);                 
                [res.error_conv1,ind_conv1] = PrintErrorPos(convOne - conv_BarkerHenderson_IC(this.IDC.Pts.y2_kv),'convolution',this.IDC.Pts.y2_kv);     
                res.error_conv1_posy2 = this.IDC.Pts.y2_kv(ind_conv1);                
            else(strcmp(this.optsPhys.V2.V2DV2,'Phi2DLongRange'))
                 cprintf('m','CheckMeanfieldConvolution: Case not yet implemented\n');
            end
         end
     else
         cprintf('m','CheckMeanfieldConvolution: Case not yet implemented\n');
     end
          
     
    function z = conv_Phi2DLongRange(y2)
         z = this.optsPhys.V2.epsilon*(- pi^2/2 + ...
                                      + pi*atan(-y2)+...
                                      - pi*(y2)./(1+y2.^2));
        z(y2==inf) =  this.optsPhys.V2.epsilon*(- pi^2);
    end
    function C = conv_BarkerHenderson_IC(y2)        
        C = zeros(size(y2));
        y2Min = this.IDC.y2Min;
        y2Max = this.IDC.y2Max;
        d     = 1;
        
        mark1    = (y2 <= y2Min + d);
        C(mark1) = -6/5*pi*(y2(mark1)+d-y2Min) + ...
                   BH_Int(y2Max - y2(mark1)) - BH_Int(d);
        
        mark2    = ((y2 > y2Min + d) & (y2 < y2Max - d));
        C(mark2) = -6/5*pi*2*d + ...
                    BH_Int(y2Max-y2(mark2)) - BH_Int(d) + ...
                    BH_Int(y2(mark2)-y2Min) - BH_Int(d);
                
        mark3    = (y2 >= y2Max - d);
        C(mark3) = -6/5*pi*(y2Max-y2(mark3)+d) + ...                    
                    + BH_Int(y2(mark3)-y2Min) - BH_Int(d);                

        C = C*this.optsPhys.V2.epsilon;
    end

    function z = conv_BarkerHendersonCutoff2D(y2_h)
        z = 2*a-this.optsPhys.V2.epsilon*BH_Psi(y2_h);
    end
    function z = conv_ConstShortRange(y2_h)
        z = zeros(size(y2_h));
        
        lambda = this.optsPhys.V2.lambda;
        z(y2_h > lambda) = 2*a;
        
        alpha  = -2/3*pi*(lambda^3 - 1);    
        c      = (-16/9*pi)/alpha;    
        
        markG1      = (y2_h > 1) & (y2_h < lambda);
        hh          = lambda - y2_h(markG1);
        z(markG1)   = 2*a + c*pi/3*hh.^2.*(3*lambda-hh);
        
        markL1      = (y2_h <= 1);
        hh          = lambda - y2_h(markL1);
        h1          = 1 - y2_h(markL1);
        z(markL1) = - c*(4/3*pi*lambda^3 - pi/3*hh.^2.*(3*lambda-hh) -...
                        (4/3*pi - pi/3*h1.^2.*(3-h1)));
                                                                 
    end

    function z = conv_ExpDouble2D(y2_h)
        lambda      = 1;
        alpha      = 2*pi*(1/(2*lambda*exp(lambda))-(1/4)*sqrt(pi)*erf(sqrt(lambda))/lambda^(3/2)+(1/4)*sqrt(pi)/lambda^(3/2));        
        c           = (-16/9*pi)/alpha;
                
        Psi(y2_h < 1)  = 0.5*(pi/lambda)^(3/2)*(1-erf(sqrt(lambda))) + pi/lambda*exp(-lambda)*(1-y2_h(y2_h < 1));
        Psi(y2_h >= 1) = 0.5*(pi/lambda)^(3/2)*(1-erf(sqrt(lambda)*y2_h(y2_h >= 1)));
        z              = 2*a - c*this.optsPhys.V2.epsilon*Psi';
    end

    function z = conv_BarkerHenderson2D(y2_h)
        Psi(y2_h < 1)  = -16/9*pi +6/5*pi*y2_h(y2_h < 1);         
        Psi(y2_h >= 1) = 4*pi*(1./(45*y2_h(y2_h >= 1).^9) - 1./(6*y2_h(y2_h >= 1).^3));
        z              = 2*a-this.optsPhys.V2.epsilon*Psi';
    end
    function v = BH_Int(z)
        v = 2*pi/5*(-2./(9*z.^9) + 5./(3*z.^3));
    end
    function z = conv_BarkerHenderson2D_rg1(y2_h)        
        z = 2/3*pi./y2_h.^3 - 4/45*pi./y2_h.^9 - 25/32*pi^2;
    end

    function z = conv_BarkerHenderson2D_rl1(y2_h)        
        z = -32/9*pi+25/32*pi^2;
    end


end