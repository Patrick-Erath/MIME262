clear all;
clear all;

%-----------------------------
% Define variables

q = 1.6e-19; % Electron charge in coulombs
kB = 1.380e-23; %Boltzmann's constant
hbar = 1.054e-34; % Units J*s
m = 9.11e-31; % Units in kg
Nc = 3.23e19; % SI (1/cm^3)
vf = 1e6/1e2; % (cm/s)
li = 100e-9*1e2; % (cm)
b0 = 5e3;

%-----------------------------
% Solutions 1b)

Ec = -4.0; % eV Si
NT = 100
Nf = 110;
Ef = -4.2; % Ev, n-type
T = linspace(100,400,NT); % Kelvin
n_metal = 6e23; % approx. one e- per atom (cm^(-3))
rho_metal = ((m*vf)/(Nc*q^2))*(1/li + b0*T).*exp((Ec-Ef)./(kB*T/q)); %convert m to cm

%plot(T,A./rho_metal,'b')
%plot(T,rho_metal,'b',T,rho_si,'r');

semilogy(T,rho_metal,'b',T,rho_si,'r');
legend('\rho_{e,metal}','\rho_{e,Si}');
xlabel('T (K)','FontSize',16);
ylabel('log(\rho)','FontSize',16);
grid on,


