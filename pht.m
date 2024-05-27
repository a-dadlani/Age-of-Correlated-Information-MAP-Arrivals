function [ph1,al,mm] = pht(x,y)
%%%%%%% E1
switch x 
case 1  
    mm = 1;
    TT = (1/y)*[-1];
    ala = zeros(1,mm);
    ala(1,1) = 1;
    al =  ala;
    ph1 = TT;
%%%%%%%%%%% E04
case 2 
     mm = 4;
     TT = (4/y)*[-1,1,0,0;0,-1,1,0;0,0,-1,1;0,0,0,-1];
     ala = zeros(1,mm);
     ala(1,1) = 1;
     al =  ala;
     ph1 = TT;
   
%%%%%%%% Hyperexponential
case 3 
    mm = 2;
    TT = (0.19/y)*[-10,0;0,-1];
    ala = zeros(1,mm);
    ala = [0.9,0.1];
    al =  ala;
    ph1 = TT;
end
    
