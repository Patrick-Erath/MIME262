% -----------------------------------------------------------------
%  H1 Calculations MIME262 -- Can be run in Matlab or Octave
% -----------------------------------------------------------------


clear all;
close all;

% -----------------------------------------------------------------
% Q3
alpha = 1.75;
n = 9;
r_NaBr = 0.298e-9;
r_KBr = 0.329e-9;
Z = 1;
q = 1.6e-19;
eps0 = 8.85e-12;
NA = 6.02e23;

prefactor = ((NA*alpha*Z*Z*q*q)/(4*pi*eps0)) * (1 - 1/n);

UL_NaBr = prefactor/r_NaBr
UL_KBr = prefactor/r_KBr

% -----------------------------------------------------------------
% Q4
Au_density = 19281;  % Kg/m^3
Au_molar_mass = 0.19697;   % Kg/mol
NA = 6.02e23;   % Atoms/mol
Ne_per_atom = 1;   % Number of electrons per atom

Au_electron_gas_density = Ne_per_atom*(Au_density/Au_molar_mass)*NA

% -----------------------------------------------------------------
% Q5
A = 13.6;  % eV
h = 4.136e-15;  % eV*s
c = 3e8; % m/s
n_start = 5;
n_end = [1 2 3];  % Array

E_transition = -A*(1/n_start^2 - 1./(n_end.*n_end))

v_transition = E_transition/h

wavelength_transition = c./v_transition

% -----------------------------------------------------------------
% Q6
h = 6.626e-34;  % J*s
m = 9.11e-31;  % kg
Ag_molar_mass = 0.1079 ; %kg/mol
Ag_density = 10500;  % kg/m^3
NA = 6.02e23;

N_per_V = (Ag_density/Ag_molar_mass)*NA

Ef_J = (h^2/(8*m))*((3/pi)*N_per_V)^(2/3)
Ef_eV = Ef_J/1.6e-19



