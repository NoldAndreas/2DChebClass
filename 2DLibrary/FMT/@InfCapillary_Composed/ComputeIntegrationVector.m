function [Int,Int1,Int2] = ComputeIntegrationVector(this)    
    [IntB,IntB1,IntB2]    = this.Bottom_Strip.ComputeIntegrationVector();
    [IntM,IntM1,IntM2]    = this.Main_Strip.ComputeIntegrationVector();
    [IntT,IntT1,IntT2]    = this.Top_Strip.ComputeIntegrationVector();    

    if(max(abs(IntB1-IntM1)) ~= 0)
        PrintErrorPos(IntB1-IntM1,...
            'InfCapillary_Composed:ComputeIntegrationVector: Error a');
    end
    
    if(max(abs(IntB1-IntT1)) ~= 0)
        PrintErrorPos(IntB1-IntT1,...
            'InfCapillary_Composed:ComputeIntegrationVector: Error a');
    end
    
    Int  = [IntB,IntM,IntT];
    Int1 = IntB1;
    Int2 = [IntB2,IntM2,IntT2];
end