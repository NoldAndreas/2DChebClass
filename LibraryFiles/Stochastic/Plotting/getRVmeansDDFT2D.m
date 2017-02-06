function [meanR,meanV]=getRVmeansDDFT2D(ddft,geom)

switch geom
    
    case 'planar2D'
        
        [meanR,meanV]=getRVmeansDDFTPlanar2D(ddft);
        
    case 'polar2D'
        
        [meanR,meanV]=getRVmeansDDFTPolar2D(ddft);
end

end