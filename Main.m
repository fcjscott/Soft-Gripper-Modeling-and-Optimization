%% MAE 259 
clear all;
clc;
%-----------------------------------------------
%% Global Variables
global graspRange fingerForce gripperMass runTime E
global consInd Pressure

E = 1e7;
graspRange = [0.015 0.065];  % Meters
fingerForce = [0 11.12];   % Newtons
gripperMass = 0.368;        % kilogram
runTime = 0.3;
Pressure = 40
fingerSpacing = 0.02;        % Meter
grasp = false;
object = "Ball";
weight = 1*9.81;     % Newton
miu = 0.25;
i = 1;
j = 1;
%DER(fingerSpacing, verticalMove, object, runTime);
for angle = 45:5:80
    i = 1;
    for verticalMove = -0.1:-0.002:-0.12
        [nodesNum, validNum, crossNum] = DER(fingerSpacing, verticalMove, object, runTime, angle, E);
        fixedNodes(i,j) = nodesNum - 4;
        if validNum == 0 || crossNum == 0
            fixedNodes(i,j) = 0;
        end
        if object == "Cubic"
            N_G_ration(i,j) = 1/((fixedNodes(i,j)/2)/miu);
%          elseif object == "Ball"
%              continue;
%              
         end
        i = i + 1;
    end
    j = j + 1;
end
angle = 45:5:80;
verticalMove = -0.1:-0.002:-0.12;
[M, I] = max(fixedNodes(:));
[I_row, I_col] = ind2sub(size(fixedNodes), I);
maxNodes = DER(fingerSpacing, verticalMove(I_row), object,... 
runTime, angle(I_col), E);
fprintf("optimal angle is %d and optimal vertical move is %f m \n",...
angle(I_col),verticalMove(I_row));
consInd;

