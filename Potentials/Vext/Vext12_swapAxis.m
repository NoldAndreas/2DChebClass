function [V,VDiff,VInt] = Vext12_swapAxis(y1,y2,intBound)

    if(nargin == 3)
        [V,VDiffS,VInt] = Vext12(y2,y1,intBound);
    else
        [V,VDiffS] = Vext12(y2,y1);
    end
    VDiff          = struct('dy1',VDiffS.dy2,'dy2',VDiffS.dy1,'ddy1',VDiffS.ddy2,'ddy2',VDiffS.ddy1,'dy1dy2',VDiffS.dy1dy2);     

end