%% Quiz: Parameter Update Course Sebastian Thrun
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu1 = 10;                            
simgaSquare1 = 8;                
mu2 = 13;           
simgaSquare2 = 2;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiply Gaussian1 and Gaussian2
mu=(1/(simgaSquare1+simgaSquare2))*(simgaSquare2*mu1+simgaSquare1*mu2)
simgaSquare=1/((1/simgaSquare1)+(1/simgaSquare2))