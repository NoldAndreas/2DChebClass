function ComputeHardSphereMatrices(this)

    if(strcmp(this.optsNum.FexNum.Fex,'CarnahanStarling'))
        this.IntMatrFex = [];
        return;
    else

        fprintf(1,'Computing Fex matrices ...\n');   
        params.sigmaS   = this.optsPhys.sigmaS;
        params.FexNum   = this.optsNum.FexNum;
        params.PhysArea = this.optsNum.PhysArea;
        params.Polar    = 'cart';
        params.Comments = this.configName;
        func            = str2func(['FexMatrices_',this.optsNum.FexNum.Fex]);
        [this.IntMatrFex,recFex] = DataStorage(['HalfSpace_FMT' filesep func2str(func)],func,params,this.HS); 

        %Test
        if(recFex)
            CheckAverageDensities_Rosenfeld_3D(this.HS,this.IntMatrFex,true);
        else
            CheckAverageDensities_Rosenfeld_3D(this.HS,this.IntMatrFex,false);
        end
    end
end