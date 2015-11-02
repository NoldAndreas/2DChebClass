function SimulationSplitDisk_BarkerHenderson()

    y2_h = 1;

    shapeHS.N     = [60,60];
    shapeHS.L1    = 2;
    shapeHS.L2    = 2;
    shapeHS.y2Min = 0.5;
    HS            = HalfSpace(shapeHS);
    
    [~,~,alpha]     = BarkerHenderson_2D(0);

    Conv = GetConv(y2_h,[30,30],1);    

	%y2_h           = PtsCart.y2_kv(Pts.y1_kv==inf) - R;
    y2_dw          = y2_h - shapeHS.y2Min;
    Psi(y2_dw < 1)  = -16/9*pi +6/5*pi*y2_dw(y2_dw < 1);         
    Psi(y2_dw >= 1) = 4*pi*(1./(45*y2_dw(y2_dw >= 1).^9) - 1./(6*y2_dw(y2_dw >= 1).^3));
    check          = 2*alpha-Psi;
    PrintErrorPos(sum(Conv) - check','Error for convolution at y2_dw = ',y2_dw);
    
    %Check convergence of integration for e^(-2*x)
    z = test(HS.Pts.y2_kv);        
    
    n = 16:4:40;
    L = 1;    
    for i=1:length(n);
        Conv   = GetConv(y2_h,[n(i),n(i)],L);
        v(i)   = Conv*z;
    end
    
    plot(n,v); hold on;
    
    L = 1.5;    
    for i=1:length(n);
        Conv   = GetConv(y2_h,[n(i),n(i)],L);
        v(i)   = Conv*z;
    end
    plot(n,v,'r'); hold on;
    
    L = 0.5;    
    for i=1:length(n);
        Conv   = GetConv(y2_h,[n(i),n(i)],L);
        v(i)   = Conv*z;
    end
    plot(n,v,'m'); hold on;
    
    
    function z= test(y2)        
        z = 0.5+3.5*exp(-3*(y2-0.5))+cos(pi*0.5+1.2*(y2-0.5)*pi)./(y2);
        z(y2 == inf) = 0.5;
    end


        
    function Conv = GetConv(offset,N,L)
        
        shapeAnn.RMin = 1;
        shapeAnn.L    = L;
        shapeAnn.N    = N;%[30,30];
        shapeAnn.Origin = [0,offset];
        annArea       = InfAnnulus(shapeAnn);

        annPts        = Intersect(HS,annArea);

        shapeDisc.R   = 1;
        shapeDisc.N   = N;%[30,30];
        shapeDisc.Origin = [0,offset];
        diskArea      = Disc(shapeDisc);    

        diskPts       = Intersect(HS,diskArea);
        
%        scatter(annPts.pts.y1_kv,annPts.pts.y2_kv); hold on;
%        scatter(diskPts.pts.y1_kv,diskPts.pts.y2_kv);
%        xlim([-5 5]); ylim([0 5]);

        %[annBH]         = BarkerHenderson_2D(annPts.ptsPolLoc.y1_kv,1);
        %[diskBH]        = BarkerHenderson_2D(diskPts.ptsPolLoc.y1_kv,1);
        [annBH]         = BarkerHenderson_2D(annPts.ptsPolLoc);
        [diskBH]        = BarkerHenderson_2D(diskPts.ptsPolLoc);        
        
        IPAnn  = HS.SubShapePtsCart(annPts.pts);
        IPDisk = HS.SubShapePtsCart(diskPts.pts);

        Conv            = [diskPts.int.*diskBH',annPts.int.*annBH']*...
                                                    [IPDisk;IPAnn];
                                                
        %markAnn = ((annPts.pts.y1_kv < 5) & (annPts.pts.y1_kv > -5) & (annPts.pts.y2_kv < 5));
        %markDisk = ((diskPts.pts.y1_kv < 5) & (diskPts.pts.y1_kv > -5) & (diskPts.pts.y2_kv < 5));

%        figure;
%        scatter3(annPts.pts.y1_kv(markAnn),annPts.pts.y2_kv(markAnn),annBH(markAnn)); hold on;
%        scatter3(diskPts.pts.y1_kv(markDisk),diskPts.pts.y2_kv(markDisk),diskBH(markDisk));       

        %Get function value
        %Check Integration with 1-function:    
        %disp(diskPts.int*diskBH + annPts.int*annBH - 2*alpha);                                                    
    end

end