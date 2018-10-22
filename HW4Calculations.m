% -----------------------------------------------------------------
%  HW4 Calculations MIME262 -- Can be run in Matlab or Octave
% -----------------------------------------------------------------


clear all;
close all;

% ----------------------------------------------------------------- 
% Q5 
C_MolarMass = 12.01;
H_MolarMass = 1.008;

M = [7500:5000:47500]'

x = [0.01;
     0.09;
     0.17;
     0.18;
     0.20;
     0.17;
     0.09;
     0.06;
     0.03]

w = x.*M/(sum(x.*M))

Mw = sum(w.*M)

Dimer_MolarMass = 3*(C_MolarMass)+6*(H_MolarMass)

Nw = Mw/Dimer_MolarMass

% ----------------------------------------------------------------- 
% Q6
O_MolarMass = 16;

PETDimer_MolarMass = 4*O_MolarMass + 10*C_MolarMass + 8*H_MolarMass

Nmol = 1000/PETDimer_MolarMass

Ndimer = Nmol*6.0221413e+23

PETDimerWithH2O_MolarMass = PETDimer_MolarMass + 2*O_MolarMass + 4*H_MolarMass

MassReactants = Nmol*PETDimerWithH2O_MolarMass





