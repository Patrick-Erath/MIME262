clear all;
close all;

% ----------------------------------------------------------------
% Q3
clear all;
q = 1.6e-19;  % Electron charge coulombs
kB = 1.380e-23;  % Boltzmann's constant
NE = 500;
Ef = 0;
E = linspace(-5,5,NE);  % Plotting range in eV
NT = 100;
T = linspace(1e-9,300,NT);  % Temperature range in Kelvin

f = zeros(NT,NE);
for ii=1:NT
    f(ii,:) = 1./(1+exp((E-Ef)/(kB*T(ii)/q)));  % keep same units for energy
end

plot(E,f);

%break;

% ----------------------------------------------------------------
% Q4
clear all;
Rb_Resistivity = 11.5e-8;   % Ohm*m
Rb_ConductingElectronsPerAtom = 1; 
Rb_NumberOfAtomsPerUnitCell = 2;   % Bcc unit cell
Rb_LatticeLength = 0.5705e-9;  % m
Rb_FermiVelocity = 8.1e7;  % m/s
Rb_ElectronMass = 9.109e-31;  % Kg
e = 1.6e-19;

Rb_Conductivity = 1/Rb_Resistivity;

Rb_ElectronDensity = (Rb_NumberOfAtomsPerUnitCell*Rb_ConductingElectronsPerAtom)/(Rb_LatticeLength^3)

Rb_RelaxationTime = (Rb_Conductivity*Rb_ElectronMass)/(Rb_ElectronDensity*e^2)
Rb_MeanFreePath = Rb_RelaxationTime*Rb_FermiVelocity   % Units of m
Rb_MeanFreePathInTermsOfLatticeConstants = (Rb_MeanFreePath/Rb_LatticeLength)
