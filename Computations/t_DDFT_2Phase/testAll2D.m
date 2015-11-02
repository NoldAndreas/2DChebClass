function testAll2D

    AddPaths();                
  
    testSimulations();
    testDiffusion();
    test2PhaseDDFT();        
    
%    testDDFT();       
%    testNSpeciesDDFT();           
%    testFMT();    
  
    %testDDFTInertia();
    %testCahnHilliard();    
end

function testDDFTInertia()        
    DDFT_InertiaBox_2Phase_1_simple();    
end
function testCahnHilliard()        
    Test_CahnHillardPolar();
    Test_CahnHillardBox_Vel();
    Test_CahnHillardBox();
end

function testNSpeciesDDFT()
        
%    Test_DDFT_DiffusionDisk_NSpecies();    
    Test_DDFT_DiffusionBox_NSpecies();    
    Test0_DDFT_DiffusionBox_2Species();    
    Test0_DDFT_DiffusionDisk_2Species();


    DDFT_DiffusionHalfSpace_NSpecies();
    Test_DDFT_DiffusionPolar_NSpecies();
    
    %Test_DDFT_DiffusionPlanar_NSpecies();
    
    

    
%    TestFMT_DDFT_DiffusionHalfSpace_NSpecies();
%    TestFMT_DDFT_DiffusionPolar_NSpecies();
    
%    Test_DDFT_DiffusionDisk_NSpecies();
%    Test_DDFT_DiffusionBox_NSpecies();
    
    %Test0_DDFT_DiffusionBox_2Species();    
  
end
function testSimulations()
    SimulationBox_M1();
    SimulationInfBox_M1(); 

    SimulationBox_Tref();
   % SimulationHalfInfiniteCapillary();
    SimulationHalfSpace();
    SimulationHalfSpace_Composed();
    SimulationHalfSpace_Tref();
    SimulationInfiniteCapillary();
    
    SimulationTriangle();
%    SimulationWedge_M1();         
    
    SimulationDisk();
    SimulationWedge();  
    
    SimulationBox();
    SimulationBoxFourier(); 
    SimulationEvenBoxFourier();

    SimulationBall();
    SimulationBigSegment();
    SimulationSegment();
    SimulationSegment2();
    SimulationTriangle();

    SimulationInfinityWedge();


 %   SimulationBoxFD();
    
    SimulationPolarInfinity();
    
    SimulationHalfSpaceSkewed();
           
end
function testDiffusion()

  %  DiffusionClosedBoxFD();

    %Diffusion Equations    
    DiffusionClosedBox();
    

    DiffusionWedge();
    DiffusionAdvectionPolarInfinity();    
    	
    DiffusionCapillary();
    DiffusionPeriodicSlit(); 
    DiffusionDisk();           
end

function testDDFT()
    %*********************************************
    %TODO:
    %Default_InfiniteCapillary_DDFT_DiffusionClosedBox();    
    %*********************************************           
    
    Default1_DDFT_DiffusionWedge();    
    
    Default1_DDFT_DiffusionPolarInfinity();
    
    
    Default1_DDFT_DiffusionPolarDisk();        
    
	Default1_DDFT_DiffusionBox();          
    Default1_DDFT_DiffusionSector();   
    Default_InfiniteCapillary_DDFT_DiffusionClosedBox();    
    %Default1_DDFT_DiffusionWedgeInfinity();        
    
    
end

function test2PhaseDDFT()
    DDFT_DiffusionBox_2Phase_1();        
    DDFT_DiffusionBox_2Phase_Sat_1();    
    %Test_DDFT_EqBoxInf_1Phase();        
    DDFT_DiffusionBox_2Phase_2();
    DDFT_DiffusionBox_2Phase_3();
    
    
    %Test1_DFT_1DBox_2Phase();
end

function testFMT()        
    DDFT_DiffusionInfSpace_NSpecies();
    CheckFMTComputation();
    DDFT_DiffusionHalfSpace_NSpecies();    
    
end