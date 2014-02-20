function M_conv  = ComputeConvolutionMatrix_PointwiseX(this,f)
%*********************************************************************
%ONLY TESTED FOR HALF SPACE
%DESCRIPTION:
% computes matrix M_conv such that
%           M_conv * g_i = int f(y1-y1t,y2-y2t)*g(y1t,y2t) dy1t dy2t,
%INPUT:
%  f   = function for which convolution is to be computed
%  Pts = data for grid points in computational and physical space:
%        - 'y1_kv' - y1-variable for each point of 2D-grid 
%        - 'y2_kv' - y2-variable for each point of 2D-grid
%        - 'x1'    - grid in computational space only fo 1st variable
%        - 'x2'    - grid in computational space only for 2nd variable
%  Int = Integration weight vector
%  MapSub.PhysSpace1(y10,y20,x1,optsPhys)
%  MapSub.PhysSpace2(y10,y20,x2,optsPhys)
%OUTPUT:
% M_conv = convolution matrix, size (N1*N2,N1*N2) , see above
%*********************************************************************

   %*******************************
    %Initialization
	disp('Computing Convolution matrices...');     

    N1  = this.N1;  
    N2  = this.N2;   
    Pts = this.Pts;
    
    n1  = this.ConvN(1);  
    n2  = this.ConvN(2);
    
    M_conv = zeros(N1*N2,N1*N2);
    
    %1) Compute Subgrid    
    [x1,wInt1]  = ClenCurtFlip(n1-1);  
    [x2,wInt2]  = ClenCurtFlip(n2-1); 	                        
    
    %ptsC_Cart = GetCartPts(this,ptsC.y1_kv,ptsC.y2_kv);
    PtsCartY1_0  = GetCartPts(this,zeros(N2,1),this.Pts.y2);    
    PtsCart      = GetCartPts(this);
    %*******************************
    opts.y2Min  = this.y2Min*sin(this.alpha);
    opts.LConv  = this.LConv;
    opts.L2Conv = this.L2Conv;
    
    hw = waitbar(0,'Computing Convolution Matrix');
    
    [y1DC,dy1DC]  = PhysSpace1_Sub(this,[],[],x1,x1,@fy1);
    y1_kvDC       = kronecker(y1DC,ones(n2,1));
    
    %Go through all points in box
    for i2Pts = 1:N2
        
        y20     = Pts.y2_kv(i2Pts);
        y20Cart = PtsCartY1_0.y2_kv(i2Pts);
                
        %Get x2 Points for subgrid
        if(y20 ~= inf)
            [y2DC,dy2DC]   = PhysSpace2_Sub(this,[],y20Cart,x2,x2,@fy2,opts); 
            y2C            = y2DC + y20Cart;    
            pts            = GetInvCartPts(this,zeros(size(y2C)),y2C);
            xx2            = CompSpace2(this,pts.y2_kv);
            Interp2        = barychebevalMatrix(Pts.x2,xx2);  
        else
            [y2DC,dy2DC]   = PhysSpace2_Sub(this,[],y20,x2,x2,@fy2);                        
            Interp2        = barychebevalMatrix(Pts.x2,ones(size(y2DC)));  
        end
        
        y2_kvDC = kronecker(ones(n1,1),y2DC);
        
        ptsD = GetInvCartPts(this,y1_kvDC,y2_kvDC);
        
        for i1Pts = 1:N1
            
            iPts    = i2Pts + (i1Pts-1)*N2;
                        
            y2_kvDC = kronecker(ones(n1,1),y2DC);
                        
            fP      = f(sqrt(y1_kvDC.^2+y2_kvDC.^2));

            %Interpolate onto the new set of points
            IP = zeros(n1*n2,N1*N2);
            
            if(Pts.y1(i1Pts) == inf)
                y1Interp = inf*ones(n1*n2,1);
            elseif(Pts.y1(i1Pts) == -inf)
                y1Interp = -inf*ones(n1*n2,1);
            else
                y1Interp = Pts.y1(i1Pts) + ptsD.y1_kv;%y1DC;
            end
            %y1DC+Pts.y1(i1Pts)
            
            interp1 = CompSpace1(this,y1Interp);
            for i = 1:n1*n2
                i2_sub  = mod(i,n2);
                if(i2_sub == 0)
                    i2_sub = n2;
                end
                %i1_sub  = ceil((i/n2));
                Interp1 = barychebevalMatrix(this.Pts.x1,interp1(i));
                % kron form of the two interpolations                                                            
                IP(i,:) = kronecker(Interp1,Interp2(i2_sub,:));
            end              
            
            Int_i   = kron(dy1DC'.*wInt1,dy2DC'.*wInt2);
            Int_i(Int_i==Inf) = 0;        

            %Finally get correct integration weight            
            M_conv(iPts,:) = (Int_i.*fP')*IP;

        end
        hw = waitbar(i2Pts/N2);
    end
    
    close(hw);
   
    function y = fy2(y2)
        y = f(y2);
    end
 
    function y = fy1(y1)
        y = f(y1);
    end
end

    %********************************************************
    %********************************************************
    %********************************************************
    %********************************************************
    %********************************************************
% 
% 
%     disp('Computing Convolution matrices...'); 
%     boolPlot = false;
% 
%     N1  = this.N1;  
%     N2  = this.N2;   
%     Pts = this.Pts;
%     
%     n1  = this.ConvN(1);  
%     n2  = this.ConvN(2);
%     
%     M_conv = zeros(N1*N2,N1*N2);
%     
%     %1) Compute Subgrid    
%     [x1,wInt1]  = ClenCurtFlip(n1-1);  
%     [x2,wInt2]  = ClenCurtFlip(n2-1); 	        
%     
%     %2) Loop over all points
%     for i=1:(N1*N2) 
%         y10 = Pts.y1_kv(i); y20 = Pts.y2_kv(i);
%         
%         %Analyze epsilon
%         
%         if(boolPlot)
%             figure('Color','white');
%             subplot(1,2,1);
%             [y2D,dy2D,d_Pade,ep_Pade] = PhysSpace2_Sub(this,y10,y20,x2,x2,@fy2);
%             fP                        = fy2(y2D);            
%             
%             yy = (max(-5,this.y2Min-y20):0.01:5)';
%             plot(y2D,fP,'o','MarkerSize',10,'MarkerFace','b'); hold on;
%             AInterp = barychebevalMatrix(x2,CompSpace2_Sub(this,y10,y20,yy,d_Pade,ep_Pade));
%             plot(yy,AInterp*fP,'b');
%             
%             fAna = fy2(yy);
%             plot(yy,fAna,'--g');
%             if(~isempty(d_Pade))
%                 d_y = PhysSpace2_Sub(this,y10,y20,InvTref(d_Pade,d_Pade,ep_Pade),x2,@fy2);
%                 plot([d_y,d_y],[min(fAna),max(fAna)],'-.k');
%             end
%             xlim([min(yy)  max(yy)]);
%             
%             subplot(1,2,2);            
%             xx = (-1:0.01:1)';
%             AInterp = barychebevalMatrix(x2,xx);
%             plot(x2,fP,'o','MarkerSize',5,'MarkerFace','r'); hold on;
%             plot(xx,AInterp*fP,'r');
% 
%             %Plot the actual function
%             y2Inter  = PhysSpace2_Sub(this,y10,y20,xx,x2,@fy2);
%             fAna     = f(zeros(size(y2Inter)),y2Inter);
%             plot(xx,fAna,'--g','LineWidth',1);            
%         end
%         if(y20 ~= inf)
%             [y2D,dy2D]     = PhysSpace2_Sub(this,y10,y20,x2,x2,@fy2);
%             y2             = y2D + y20;
%             xx2            = CompSpace2(this,y2);
%             Interp2        = barychebevalMatrix(Pts.x2,xx2);  
% %            dy2D    = dy2D.*dx2Tee;
%         else
%             y2         = y20*ones(size(y2D));
%             Interp2    = barychebevalMatrix(Pts.x2,CompSpace2(this,y2));  
%             
%             [y2D,dy2D] = PhysSpace2_Sub(this,y10,y20,x2);
%         end
%         
%         %2a) Compute points in subgrid        
%         [y1D,dy1D]  = PhysSpace1_Sub(this,y10,y20,x1,x1,@fy1);        
%         y1D_kv      = kronecker(y1D,ones(size(y2D)));
%         y2D_kv      = kronecker(ones(size(y1D)),y2D);
%         
%         %2b) Compute f-values        
%         fP          = f(GetDistance(this,y1D_kv,y2D_kv));
%         
%         %2c) Compute Interpolation Matrix onto Subspace
%         if((y10 == inf) || (y10 == -inf))
%             y1 = y10*ones(size(y1D));
%         else
%             y1 = y1D + y10;
%         end
%         
%         Interp1   = barychebevalMatrix(Pts.x1,CompSpace1(this,y1));                                       
%         Interp_i  = kron(Interp1,Interp2);
%         
%         %2d) Compute integration weight vector
%         Int_i       = kron(dy1D'.*wInt1,dy2D'.*wInt2);
%         Int_i(Int_i==Inf) = 0;
%                         
%         %2e) Combine
%         M_conv(i,:) = (Int_i.*fP')*Interp_i;    
%     end
%     M_conv(isnan(M_conv)) = 0;
%    
%     if(nargin < 3)
%         return;
%     end    
%     
% 
% end

% figure('Color','white');
%             [y2D,dy2D]  = PhysSpace2_Sub(this,y10,y20,x2);
%             fP          = f(zeros(size(y2D)),y2D);
%             [d,ep]      = GetPolesPadeApproximation(fP,2);
% 
%             [x2Tee,dx2Tee]    = Tref(x2,d,ep);
%             [y2DTee,dy2DTee]  = PhysSpace2_Sub(this,y10,y20,x2Tee);
%             fPTee             = f(zeros(size(y2DTee)),y2DTee);        
% 
%             xx = (-1:0.01:1)';
%             AInterp = barychebevalMatrix(x2,xx);
%             plot(x2,fP,'o','MarkerSize',10,'MarkerFace','b'); hold on;
%             plot(xx,AInterp*fP,'b');
% 
%             xxTee       = InvTref(xx,d,ep);
%             ATee_Interp = barychebevalMatrix(x2,xxTee);
%             plot(x2Tee,fPTee,'o','MarkerSize',5,'MarkerFace','r');
%             plot(xx,ATee_Interp*fPTee,'r');
% 
%             %Plot the actual function
%             y2Inter  = PhysSpace2_Sub(this,y10,y20,xx);
%             fAna     = f(zeros(size(y2Inter)),y2Inter);
%             plot(xx,fAna,'--g','LineWidth',1);
% 
%             plot([d,d],[min(fAna),max(fAna)],'-.k');