function [Int,Int1,Int2] = ComputeIntegrationVector(this)    
    [Int_SS,Int1_SS,Int2_SS]    = this.Sub_Strip.ComputeIntegrationVector();
    [Int_SHS,Int1_SHS,Int2_SHS] = this.Sub_HalfSpace.ComputeIntegrationVector();

    if(max(abs(Int1_SS-Int1_SHS)) ~= 0)
        PrintErrorPos(Int1_SS-Int1_SHS,...
            'Composed Half Space: Agreement between upper and lower integreation of 1st dimension');
    end
    
    Int  = [Int_SS,Int_SHS];
    Int1 = Int1_SS;
    Int2 = [Int2_SS,Int2_SHS];
end