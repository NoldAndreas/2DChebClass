classdef (Abstract) TrefSpectralSpectral < M1SpectralSpectral

    properties
        %these are the trefethen-map-parameters.
        %they have the length N1
        d_x1  = 0
        ep_x1 = 0 
        SafetyFactor_ep = 0.75
        ep_threshold    = 0.5
        
        %note that this.Pts.x1 and this.Pts.x2
        %are the Chebychev-collocation points!!
    end
    %**********************************************
    %************   Initializations    *************
    %**********************************************
    methods (Access = public) 
        function this = TrefSpectralSpectral(N1,N2)
            this@M1SpectralSpectral(N1,N2);                                      
        end
        
        function  [fx_new,Pts,Diff,Int,Ind,Interp] = UpdatePadeValues(this,fx,PlotArea)
            fx_new = ComputePadeValues(this,fx);            
            InitializationPts(this);
            if(nargin >= 3)
                [Pts,Diff,Int,Ind,Interp] = ComputeAll(this,PlotArea);    
            else
                [Pts,Diff,Int,Ind] = ComputeAll(this);    
            end
        end
        
        function fx_new = ComputePadeValues(this,fx)
            
            N1   = this.Pts.N1;     N2 = this.Pts.N2;
            
            ep_new = zeros(N1,1);
            d_new  = zeros(N1,1);
            
            if(this.ep_x1 == 0)                
                for i=1:N1   
                    [d_new(i),ep_new(i)] = GetPolesPadeApproximation(fx(1+(i-1)*N2:i*N2),2);                
                end
            else
                %Search for pole in a subinterval of [-1,1] for each x1
                %(1) Iteration through each x1:
                for i=1:N1 
                    
%                     y1 = PhysSpace1Tref(this,this.Pts.x1(i));
%                     
%                     xInterp = (-1:0.01:1)';
%                     IP      = barychebevalMatrix(this.Pts.x2,xInterp);                    
%                     
%                     d_old  = this.d_x1(i);
%                     ep_old = this.ep_x1(i);
%                     
%                     %Plot Physical grid before: 
%                     x_old = M1Tref(this.Pts.x2,d_old,ep_old);
%                     z_old = PhysSpace2Tref(this,x_old);
%                                                
%                     plot(x_old,Vext17(y1,z_old),'o'); hold on;                                                               
%                     plot([d_old,d_old],[-1 1],'--');
%                     
%                     xIPold  = M1Tref(xInterp,d_old,ep_old);   
%                     plot(xIPold,IP*Vext17(y1,z_old));                    
                    
                    %(a)define the subinterval [subx1Min,subx1Max]
                    D = 10*this.ep_x1(i);
                    
                    subx1Min = max(-1,this.d_x1(i)-D);
                    subx1Max = min(1,this.d_x1(i)+D);
                    
                    %(b) map Chebychebpoints to subinterval
                    xCheb     = this.Pts.x2;
                    x2Sub     = subx1Min + (subx1Max-subx1Min)*(1+xCheb)/2;
                    x2SubCheb = InvM1Tref(x2Sub,this.d_x1(i),this.ep_x1(i));
                    
                    %(c) interpolate onto these points
                    Interp = ComputeInterpolationMatrix(this,this.Pts.x1(i),x2SubCheb,false,false);
                    fxSub  = Interp.InterPol*fx;                    
%                    plot(x2Sub,fxSub,'g');
                    
                    %(d) Compute delta and epsilon on Subinterval
                    [d,ep]         = GetPolesPadeApproximation(fxSub,2);
                    
                    %(e) Translate to full interval
                    d_new(i)   = subx1Min + (subx1Max-subx1Min)*(1+d)/2;
                    ep_new(i)  = ep*(subx1Max-subx1Min)/2;  
                    
                    %TEST:
                    %Plot Physical grid after: 
%                     d_new  = this.d_x1(i);
%                     ep_new = this.ep_x1(i); 
%                     
%                     x_new = M1Tref(this.Pts.x2,d_new,ep_new);
%                     z_new = PhysSpace2Tref(this,x_new);
%                                                
%                     plot(x_new,Vext17(y1,z_new),'or'); hold on;                                                               
%                     plot([d_new,d_new],[-1 1],'--r');
%                     
%                     xIPnew  = M1Tref(xInterp,d_new,ep_new);   
%                     plot(xIPnew,IP*Vext17(y1,z_new),'r');
%                     
%                     xlim([-1 1]);
%                                         
%                     plot(xInterp,Vext17(y1,PhysSpace2Tref(this,xInterp)),'k');
                end
            end
            
            ep_new = ep_new*this.SafetyFactor_ep;
            
            %Intermediate Step: Only set new ep/d if ep<ep_threshold
            if(max(ep_new) < this.ep_threshold)

                %(2) Interpolate onto new set of points
                %  a. get new set of points (x1 points are the same)   
                d_newFull  = kron(d_new,ones(N2,1));
                ep_newFull = kron(ep_new,ones(N2,1));

                y2_new     = PhysSpace2Tref(this,M1Tref(this.Pts.x2_kv,d_newFull,ep_newFull));
                IP_new     = InterpolationMatrix_Pointwise(this,this.Pts.y1_kv,y2_new);
                fx_new     = IP_new*fx;

                %(3) Finally update the position of the poles
                this.d_x1  = d_new;
                this.ep_x1 = ep_new;
            else
                fx_new = fx;
                disp(['Maximal value of imaginary part of singularity position is ',...
                    num2str(max(ep_new)),'. No adaptation was done.']);
            end
        end
        
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)
            %Global Map:
            % y1 = PhysSpace1Tref(x1)
            % y2 = PhysSpace2Tref(M1Tref(x2,d(x1),ep(x1)));
            n                = length(x1);            
            [y1_kv,dy1,ddy1] = PhysSpace1Tref(this,x1);
            
            if(this.ep_x1 == 0)
                [y2_kv,dy2,ddy2] = PhysSpace2Tref(this,x2);
                if(nargout > 2)                    
                    J        = zeros(n,2,2);
                    J(:,1,1) = dy1;
                    J(:,2,2) = dy2;  
                end
                if(nargout > 3)
                    dH1        = zeros(n,2,2);
                    dH1(:,1,1) = ddy1;
                    dH2        = zeros(n,2,2);                    
                    dH2(:,2,2) = ddy2;
                end                
            else
                X1Diff           = barychebdiff(this.Pts.x1);
                X1Interp         = barychebevalMatrix(this.Pts.x1,x1);
                [x2_t,DiffTref]  = M1Tref(x2,X1Interp*this.d_x1,X1Interp*this.ep_x1);             
                [y2_kv,dy2,ddy2] = PhysSpace2Tref(this,x2_t);
                

                dd_dx1           = X1Interp*X1Diff.Dx*this.d_x1;
                dep_dx1          = X1Interp*X1Diff.Dx*this.ep_x1;
                
                dd_ddx1          = X1Interp*X1Diff.DDx*this.d_x1;
                dep_ddx1         = X1Interp*X1Diff.DDx*this.ep_x1;
                
                if(nargout > 2)
                    n        = length(x1);
                    J        = zeros(n,2,2);
                    J(:,1,1) = dy1;
                    J(:,2,1) = dy2.*(DiffTref.dyd_d.*dd_dx1 + DiffTref.dyde.*dep_dx1);
                    J(:,2,2) = dy2.*DiffTref.dydx;
                end
                
                if(nargout > 3)
                    dH1        = zeros(n,2,2);
                    dH1(:,1,1) = ddy1;
                    
                    dH2        = zeros(n,2,2);
                    dH2(:,1,1) = ddy2.*(DiffTref.dyd_d.*dd_dx1 + DiffTref.dyde.*dep_dx1).^2+...
                                  dy2.*(DiffTref.dydd_d.*(dd_dx1).^2 + DiffTref.dydde.*(dep_dx1).^2 + 2*DiffTref.dyd_ddep.*dd_dx1.*dep_dx1)+...
                                  dy2.*(DiffTref.dyd_d.*dd_ddx1 + DiffTref.dyde.*dep_ddx1);
                    dH2(:,1,2) = ddy2.*DiffTref.dydx.*(DiffTref.dyd_d.*dd_dx1 + DiffTref.dyde.*dep_dx1)+...
                                 dy2.*(DiffTref.dydxd_d.*dd_dx1 + DiffTref.dydedx.*dep_dx1);
                    dH2(:,2,1) = dH2(:,1,2);
                    dH2(:,2,2) = ddy2.*(DiffTref.dydx).^2 + dy2.*DiffTref.dyddx;                    
                end                
            end                                    
        end
        
        function [x1,x2] = CompSpace(this,y1,y2)
            x1              = CompSpace1Tref(this,y1);
            if(this.ep_x1 == 0)
                x2         = CompSpace2Tref(this,y2);                
            else
                X1Interp    = barychebevalMatrix(this.Pts.x1,x1);
                x2T         = CompSpace2Tref(this,y2);
                x2          = InvTref(x2T,X1Interp*this.d_x1,X1Interp*this.ep_x1);            
            end
        end        
        
        function PlotLineOfPoles(this,V,PlotArea)
                        
            x1plot      = (-1:0.01:1)';
            X1Interp    = barychebevalMatrix(this.Pts.x1,x1plot);            
            
            
            y1_d        = PhysSpace1Tref(this,x1plot);
            y2_d        = PhysSpace2Tref(this,X1Interp*this.d_x1);                                    
            
            y2_d_epP    = PhysSpace2Tref(this,X1Interp*(this.d_x1-this.ep_x1));            
            y2_d_epM    = PhysSpace2Tref(this,X1Interp*(this.d_x1+this.ep_x1));                        
            
            y2_d_epP    = min(max(min(this.Pts.y2_kv),y2_d_epP),max(this.Pts.y2_kv));
            y2_d_epM    = min(max(min(this.Pts.y2_kv),y2_d_epM),max(this.Pts.y2_kv));
            
            if(nargin > 2)
                mark_d      = ((y1_d <= PlotArea.y1Max) & (y1_d >= PlotArea.y1Min) & ...
                               (y2_d <= PlotArea.y2Max) & (y2_d >= PlotArea.y2Min));

                mark_epP      = ((y1_d <= PlotArea.y1Max) & (y1_d >= PlotArea.y1Min) & ...
                               (y2_d_epP <= PlotArea.y2Max) & (y2_d_epP >= PlotArea.y2Min));

                mark_epM    = ((y1_d <= PlotArea.y1Max) & (y1_d >= PlotArea.y1Min) & ...
                               (y2_d_epM <= PlotArea.y2Max) & (y2_d_epM >= PlotArea.y2Min));
            else
                mark_d      = true(size(y1_d));
                mark_epP    = mark_d;
                mark_epM    = mark_d;
            end
            
            hold on;
            
            if(nargin >= 2)
                IP     = InterpolationMatrix_Pointwise(this,y1_d,y2_d);
                IP_epP = InterpolationMatrix_Pointwise(this,y1_d,y2_d_epP);
                IP_epM = InterpolationMatrix_Pointwise(this,y1_d,y2_d_epM);
                plot3(y1_d(mark_d),y2_d(mark_d),IP(mark_d,:)*V,'k','linewidth',1.8);
                plot3(y1_d(mark_epM),y2_d_epM(mark_epM),IP_epM(mark_epM,:)*V,'k-.','linewidth',2);
                plot3(y1_d(mark_epP),y2_d_epP(mark_epP),IP_epP(mark_epP,:)*V,'k-.','linewidth',2);
            else
                plot(y1_d,y2_d,'k','linewidth',1.5);
                plot(y1_d(mark_epM),y2_d_epM(mark_epM),'k-.','linewidth',2);
                plot(y1_d(mark_epP),y2_d_epP(mark_epP),'k-.','linewidth',2);
            end
            
            
        end        
        function PlotGrid(this)            
            PlotGrid@Shape(this); 
            hold on;
            PlotLineOfPoles(this);
        end
        
        function PlotXGrid(this)
            
            figure('Color','white','Position',[0 0 800 800]);
            x1X = this.Pts.x1_kv;
            
            if(this.d_x1 == 0)
                x2X = this.Pts.x2_kv;                
            else
                d   = kron(this.d_x1,ones(this.N2,1));
                ep  = kron(this.ep_x1,ones(this.N2,1));
                x2X = M1Tref(this.Pts.x2_kv,d,ep);
                
                x1plot      = (-1:0.01:1)';
                X1Interp    = barychebevalMatrix(this.Pts.x1,x1plot);            
                plot(x1plot,X1Interp*this.d_x1,'k','linewidth',2); hold on;
                plot(x1plot,X1Interp*(this.d_x1 + this.ep_x1),'k-.','linewidth',2); hold on;
                plot(x1plot,X1Interp*(this.d_x1 - this.ep_x1),'k-.','linewidth',2); hold on;
            end
            
            h = scatter(x1X,x2X,'ob');              
            set(h,'MarkerEdgeColor','k','MarkerFaceColor','g');
            
            box on
            xlim([-1 1]);   ylim([-1 1]);                
            xlabel('$x_1$','Interpreter','Latex','fontsize',20);
            if(this.d_x1 == 0)
                ylabel('$x_2$','Interpreter','Latex','fontsize',20);
            else
                ylabel('$T(x_2)$','Interpreter','Latex','fontsize',20);
            end

            pbaspect([1 1 1]);
            set(gca,'fontsize',15);
        end    
        
    end
    
    %**********************************************
    %************   Mapping functions *************
    %**********************************************
    methods (Abstract = true,Access = public)         
         [y,dy,ddy] = PhysSpace1Tref(x);
         [y,dy,ddy] = PhysSpace2Tref(x);
         x = CompSpace1Tref(y);
         x = CompSpace2Tref(y);
    end  
    
end


