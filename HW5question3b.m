clear all;
close all;

% 
%Q3 
clear all;
q= 1.6e-19; % Electron charge in coulombs
kB= 1.380e-23; % Boltzmann's constant
NE = 500;
Ef = 0;
E = linspace(-1,1,NE); % Plotting range in eV
NT = 100;
T = linspace(1e-9,300,NT); % Temperature range in Kelvin

f=zeros(NT,NE)
for ii=1:NT
    f(ii,:)=1./(1+exp((E-Ef)/(kB*T(ii)/q))); %keep sace units for energy
end 

plot(E,f)

%break;
%
