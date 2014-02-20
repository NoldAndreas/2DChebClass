function TestConvolutionMatrix(this,ptsCheck,f,VisualOutput)
    
        [x1,wInt1]  = ClenCurtFlip(this.ConvN(1)-1);  
        [x2,wInt2]  = ClenCurtFlip(this.ConvN(2)-1); 	

        PtsT = struct('x1',x1,'x2',x2);
                  
        s = size(ptsCheck);
        for j=1:s(1)
        
            y10 = ptsCheck(j,1);  y20 = ptsCheck(j,2);
            [y1D,dy1D]     = PhysSpace1_Sub(this,y10,y20,x1,x1,@fy1);
            [y2D,dy2D,d_Pade,ep_Pade]     = PhysSpace2_Sub(this,y10,y20,x2,x2,@fy2);
            PtsT.y1_kv     = kron(y1D,ones(size(y2D)));
            PtsT.y2_kv     = kron(ones(size(y1D)),y2D);  
    
            y1Min  = max(- 4,min(PtsT.y1_kv)); 
            y1Max  = min(  4,max(PtsT.y1_kv));    
            y1Plot = y1Min + (0:100)'/100*(y1Max - y1Min );
            
            y2Min  = max(-4,min(PtsT.y2_kv));
            y2Max  = min(+4,max(PtsT.y2_kv));
            y2Plot = y2Min + (0:100)'/100*(y2Max - y2Min );
                
            fP          = f(GetDistance(this,PtsT.y1_kv,PtsT.y2_kv)).*...
                               TestFunction(PtsT.y1_kv+y10,PtsT.y2_kv+y20);
                                    
            Interp1 = barychebevalMatrix(PtsT.x1,CompSpace1_Sub(this,y10,y20,y1Plot));  
            Interp2 = barychebevalMatrix(PtsT.x2,CompSpace2_Sub(this,y10,y20,y2Plot,d_Pade,ep_Pade));       

            % kron form of the two interpolations                                                            
            Interp_test = struct('InterPol',kron(Interp1,Interp2),...
                        'pts1',kron(y1Plot,ones(size(y2Plot))),...
                        'pts2',kron(ones(size(y1Plot)),y2Plot),...
                        'Nplot1',length(y1Plot),'Nplot2',length(y2Plot));                    
            %Interp_test = SpectralSpectral_Interpolation(y1Plot,y2Plot,PtsT,MapsTest);
            Int_i       = kron(dy1D'.*wInt1,dy2D'.*wInt2);
            
            fPExact = f(GetDistance(this,Interp_test.pts1,Interp_test.pts2)).*...
                        TestFunction(Interp_test.pts1+y10,Interp_test.pts2+y20);
                
            posStr  = ['(y10,y20) = (',num2str(y10),' , ',num2str(y20),')'];            
            if((nargin == 4) && VisualOutput)                
                figure;
                doPlots_SC(Interp_test,PtsT,fP);
                title(posStr);
            end
            IPPts.y1_kv = Interp_test.pts1;
            IPPts.y2_kv = Interp_test.pts2;
            PrintErrorPos(fPExact - Interp_test.InterPol*fP,...
                ['Interpolation for Convolution at ',posStr],IPPts);              
        end
    

    function z = TestFunction(y1,y2)
        ep = 1;
        %z = 0.5*(1-tanh(y1/(2*ep)));
        z = ones(size(y1));
    end

    function doPlots_SC(Interp,Pts,V)

        y1s = Interp.pts1;
        y2s = Interp.pts2;
        
        ptsCart = GetCartPts(this,y1s,y2s);
        y1s = ptsCart.y1_kv;
        y2s = ptsCart.y2_kv;
        
        Pts = GetCartPts(this,Pts.y1_kv,Pts.y2_kv);

        xl = [(min(y1s)-0.5) (max(y1s)+0.5)];
        yl = [(min(y2s)-0.5) (max(y2s)+0.5)];

        y1M     = reshape(y1s,Interp.Nplot2,Interp.Nplot1);
        y2M     = reshape(y2s,Interp.Nplot2,Interp.Nplot1);            

        mask = (Pts.y1_kv <= max(y1s) & Pts.y1_kv >= min(y1s) & ...
                Pts.y2_kv <= max(y2s) & Pts.y2_kv >= min(y2s));
        if(length(V) == length(Interp.pts1))
            z = V;
        else
            z = real(Interp.InterPol*V);        
        end

        mesh(y1M,y2M,reshape(z,Interp.Nplot2,Interp.Nplot1));  hold on;
        h = scatter3(Pts.y1_kv(mask),Pts.y2_kv(mask),V(mask),'o');
        set(h,'MarkerEdgeColor','k','MarkerFaceColor','g');
        hold off;
        xlabel('y_1'); ylabel('y_2');
        pbaspect([(xl(2)-xl(1)) (yl(2)-yl(1)) 1/2*min((xl(2)-xl(1)),(yl(2)-yl(1)))]);

    end

	function y = fy1(y1)
        y = f(GetDistance(this,y1,zeros(size(y1))));
    end
    function y = fy2(y2)
        y = f(GetDistance(this,zeros(size(y2)),y2));
    end
end