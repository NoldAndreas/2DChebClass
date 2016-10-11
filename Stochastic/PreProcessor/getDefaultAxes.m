function D=getDefaultAxes(D)

switch D.geom
   
    case 'planar'
        
        D.rMin=-10;
        D.rMax=10;
        
        D.RMMin=D.rMin;
        D.RMMax=D.rMax;
        
    case 'spherical'

        D.rMin=0;
        D.rMax=10;
        
        D.RMMin=D.rMin;
        D.RMMax=D.rMax;
                
    case 'planar2D'
        
        D.rMin={{-10,-10}};
        D.rMax={{10,10}};
        
        D.RMMin=D.rMin;
        D.RMMax=D.rMax;
              
    case 'polar2D'
        
        D.rMin={{-10,-10}};
        D.rMax={{10,10}};
        
        D.RMMin={{0,0}};
        D.RMMax={{10,2*pi}};

    case 'full'
        
        D.rMin={{-10,-10,-10}};
        D.rMax={{10,10,10}};
        
        D.RMMin=D.rMin;
        D.RMMax=D.rMax;

        
end

D.pMin=D.rMin;
D.pMax=D.rMax;

D.RMin=0;
D.RMax=1;

switch D.dim
    
    case 1
        
        D.PMin=-1;
        D.PMax=1;

    case 2
        
        D.PMin={{-1,-1}};
        D.PMax={{1,1}};
        
    case 3
        
        D.PMin={{-1,-1,-1}};
        D.PMax={{1,1,1}};
        
end

D.PMMin=D.PMin;
D.PMMax=D.PMax;

end