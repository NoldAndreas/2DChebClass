function InitAnalysisGrid(this,y1Int,y2Int)

    if(~isempty(y1Int))
        [this.y1,this.Int_y1,this.DiffY1,this.DiffYY1] = InitAnalysisGridY(this,y1Int,100,'CHEB');        
    end
    
    if(~isempty(y2Int))
        [this.y2,this.Int_y2,this.DiffY2] = InitAnalysisGridY(this,y2Int,100);                      
    end
        
end