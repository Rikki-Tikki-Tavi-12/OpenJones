# OpenJones
OpenJones implements a simulation of high speed ball bearings as described by A. Jones and T. Harris (race control theory). Assumes rotating inner race and stationary outer race.

Copyright (C) 2015 Samuel Sudhof

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.    

# About the simulation
OpenJones implements the race control theory first described by Jones(doi:10.1115/1.3662587) and later expanded upon by Harris (ISBN 9780849371837 and ISBN 9780849371820). To solve the non-linear systems of this theory it uses solvers integrated with Matlab.
Variable names in the program do not follow standart conventions of capitalization, but instead are meant to represent Harris' names as best as possible in ASCII.

# Usage
OpenJones contains first and foremost a Matlab class called BallBearing, which is the meat of the simulation. To use this in a program you need to generate an instance, set the geometrical and other physical properties of your bearing, define a load case and then run the simulation your want. Currenlty available are calcDisplacement, which calculates the external and internal displacements of a ball bearing given a certain load and calcStiffness, which calculates a numerical derivative of the external displacements in three cartesian directions.
Several example of existing bearings are included in trullatest, which is meant as a sandbox.

One example from a ball bearing of the Japanese LE-7 Rocket:

%set a ball diameter (because the grove radii are given as a fraction of this)
D=7.938e-3;
%crate an instance
trulla=BallBearing();
%set the geometry
trulla.setGeometry(D,38.5e-3,20,D*0.56,D*0.52,10,0,0);
%set the meterial properties
trulla.setPhysical(7750.37, 1.99948e11, 1.99948e11, 0.25, 0.25);
%define the applied load
load=[2670,0,0,0,0,50000];
%supply an initial guess for the solution (in this case: none)
init_cond=[];
% run the simulation
[c,b]=trulla.calcDisplacement(load,'outer',init_cond); 

This yields a struct c, which contains a wealth of output information in a human-readable format. The variable names correspond to those used by Harris.

The inputs (geometry, physical and load) are described in detail in the source file.

# Compatibility
No part of OpenJones was tested with any Matlab-Compatible interpreter, such as Octave. Octave's limited selection of non-linear solvers is likely a problem for convergence in most cases.
