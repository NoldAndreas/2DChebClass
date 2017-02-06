classdef (Abstract) SpectralSpectral < M1SpectralSpectral

    properties
        ConvN
    end
    %**********************************************
    %*** *********   Initializations    *************
    %**********************************************
    methods (Access = public) 
        function this = SpectralSpectral(N1,N2)
            this@M1SpectralSpectral(N1,N2);
        end
                
        function this = InitializationPts(this)
            this = InitializationPts@Shape(this);
            
            this.Pts.y1  = PhysSpace1(this,this.Pts.x1);
            this.Pts.y2  = PhysSpace2(this,this.Pts.x2);           
        end            
        function Diff = ComputeDifferentiationMatrix(this,Sel)
            
            order = 4;

            Diff1 = barychebdiff(this.Pts.x1,order);
            Diff2 = barychebdiff(this.Pts.x2,order);    

            if(nargin == 1)
                Sel = {'Dy1' ;'DDy1' ;'DDDy1';'DDDDy1'; 'Dy2'; 'DDy2' ;'DDDy2';'DDDDy2';...
                   'Dy1Dy2'; 'DDy1Dy2'; 'Dy1DDy2' ;'Lap' ;'grad' ;'div';...
                   'gradLap' ;'gradDiv'; 'LapVec';'Lap2'};
            end
               
            Diff = PhysicalDerivatives(this,this.Pts,Sel,Diff1,Diff2,order);
            this.Diff = Diff;
        end          
        function [Int,Int1,Int2]    = ComputeIntegrationVector(this,t1Odd,t2Odd)
            if((nargin >= 2))
                Int    = ComputeIntegrationVector@M1SpectralSpectral(this,t1Odd,t2Odd);                
            else
                [h1,wInt1]  = ClenCurtFlip(this.N1-1);  
                [h1,wInt2]  = ClenCurtFlip(this.N2-1);                  
                
                [h1,dy1] = PhysSpace1(this,this.Pts.x1);
                [h1,dy2] = PhysSpace2(this,this.Pts.x2);
                
                dy1(dy1 == inf) = 0;
                dy2(dy2 == inf) = 0;
                
                Int1 = dy1'.*wInt1;
                Int2 = dy2'.*wInt2;
                                                
                Int      = kron(Int1,Int2);
                this.Int = Int;
            end
        end                
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)
           
            [y1_kv,dy1] = PhysSpace1(this,x1);
            [y2_kv,dy2] = PhysSpace2(this,x2);
            
            if(nargout > 2)
                n        = length(x1);
                J        = zeros(n,2,2);
                J(:,1,1) = dy1;
                J(:,2,2) = dy2;                
            end
            
            if(nargout > 3)
                dH1 = 0;
                dH2 = 0;
                disp('SpectralFourierClass:Comp_to_Phys: dH1/2 NOT YET IMPLEMENTED');
            end
            
        end
        function [x1,x2] = CompSpace(this,y1,y2)
            x1 = CompSpace1(this,y1);
            x2 = CompSpace2(this,y2);
        end        
        %M_conv  = ComputeConvolutionMatrix(this,f,ptsCheck);  
        %TestConvolutionMatrix(this,ptsCheck,f);
        
%        function doPlotCut1(this,f,iC)            
%             xIP = (-1:0.005:1)';
%             
%             %Const in y1-direction
%             yIP = PhysSpace2(this,xIP);
%             IP  = ComputeInterpolationMatrix(this,this.Pts.x1(iC),xIP,false); 
%                         
%             plot(yIP,IP.InterPol*f,'linewidth',2); hold on;
%             mark = (this.Pts.x1_kv == this.Pts.x1(iC));
%             plot(this.Pts.y2_kv(mark),f(mark),'o','MarkerEdgeColor','k',...
%                                                   'MarkerFaceColor','g');            
%             xlabel('$y_2$','Interpreter','Latex');
%             title(['y1 = ',num2str(this.Pts.y1(iC))]);            
%         end
%         
%         function doPlotCut2(this,f,iC)            
%             xIP = (-1:0.005:1)';
%             
%             %Const in y1-direction
%             yIP = PhysSpace1(this,xIP);
%             IP  = ComputeInterpolationMatrix(this,xIP,this.Pts.x2(iC),false); 
%                         
%             plot(yIP,IP.InterPol*f,'linewidth',2); hold on;
%             mark = (this.Pts.x2_kv == this.Pts.x2(iC));
%             plot(this.Pts.y1_kv(mark),f(mark),'o','MarkerEdgeColor','k',...
%                                                   'MarkerFaceColor','g');            
%             xlabel('$y_1$','Interpreter','Latex');
%             title(['y2 = ',num2str(this.Pts.y2(iC))]);
%         end
    end
    
    %**********************************************
    %************   Mapping functions *************
    %**********************************************
    methods (Abstract = true,Access = public)         
         [y,dy,dx,ddx,dddx,ddddx] = PhysSpace1(x);
         [y,dy,dx,ddx,dddx,ddddx] = PhysSpace2(x);
         x = CompSpace1(y);
         x = CompSpace2(y);
    end  
    
end


