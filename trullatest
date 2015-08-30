% OpenJones This is a sandbox containing several examples of ball bearings
% found in the relevant literature.
% Copyright (C) 2015 Samuel Sudhof
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.    

%Zenker Example
% clear;
% trulla=BallBearing();
% D=6.35e-3*2;
% trulla.setGeometry(D, 90e-3, 25, 6.80e-3, 6.66e-3, 18, 0, 0);
% trulla.setPhysical(7750.37, 210080000000, 210080000000, 0.3, 0.3);
% init_cond=[];
% load=[4150,0,0,0,0,1.1111111e4];
% [c,init_cond]=trulla.calcDisplacement(load,'outer',init_cond);

%Gupta ADORE Example
% clear;
% trulla=BallBearing();
% D=8e-3;
% trulla.setGeometry(D, 3.1e-2, 24, 0.56*D, 0.52*D, 6, 1.01662e-4,0,0);
% trulla.setPhysical(7750.37, 1.99948e11, 1.99948e11, 0.25, 0.25);
% init_cond=[];
% load=[2000,0,400,0,0,1.2e5];
% [c,init_cond]=trulla.calcDisplacement(load,'outer',init_cond);
% load=[2000,400,0,0,0,1.2e5];
% stiffy=trulla.calcStiffness(load, 1,0.001, init_cond);

%MA Example
% clear;
% trulla=BallBearing();
% trulla.setGeometry(6.35e-3, 31.031e-3, 25, 3.43e-3, 3.3e-3, 12, 0, 0);
% trulla.setPhysical(7750.37, 310e9, 1.99948e11, 0.27, 0.25);
% load=[500,200,0,0,0,30000];
% init_cond=[];
% [c,b]=trulla.calcDisplacement(load,'outer',init_cond);

%Nosaka 99 Tribo-Characteristics of Cryogenic Hybrid Ceramic Ball
%Bearings for Rocket Turbopumps: Bearing Wear and
%Transfer Film

% All-Steel Bearing

clear;
D=7.938e-3;
trulla=BallBearing();
trulla.setGeometry(D,38.5e-3,20,D*0.56,D*0.52,10,0,0);
trulla.setPhysical(7750.37, 1.99948e11, 1.99948e11, 0.25, 0.25);
load=[2670,0,0,0,0,50000];
init_cond=[];
[c,b]=trulla.calcDisplacement(load,'outer',init_cond);
stiffy=trulla.calcStiffness(load, b, 100*sqrt(eps));

% Hybrid Bearing

% clear;
% D=7.938e-3;
% trulla=BallBearing();
% trulla.setGeometry(D,38.5e-3,20,D*0.56,D*0.52,10,0,0);
% trulla.setPhysical(3160, 320e9, 1.99948e11, 0.26, 0.25);
% load=[2700,0,0,0,0,50000];
% init_cond=[];
% [d,e]=trulla.calcDisplacement(load,'outer',init_cond);
