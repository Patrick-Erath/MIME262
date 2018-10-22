% Solutions 1a
clear all;
close all;

% Define variables
q = 1.6e-19;  % Charge of an Electron in coulombs
kB = 1.380e-23;  % Constant : Boltzmann's constant
hbar = 1.054e-34;  % Units J*s
m = 9.1e-31;      % Units Kg
Nc = 3.23e19;  % Si [1/cm^3]
Nv = 1.83e19;  % Si [1/cm^3]
vf = 1e6/1e2; % [cm/s], http://hyperphysics.phy-astr.gsu.edu/hbase/tables/fermi.html
li = 100e-9*1e2;  % [cm]
b0 = 5e3;  

% -------------------------------------------------------------------------
Ec = -4.0; % eV Si   
NT = 100;
Nf = 110;
Ef = linspace(-4.5,-4.0,Nf);   % eV, define zero at the midgap
T = linspace(100,300,NT);  % Kelvin
n = zeros(NT,Nf);
for ii = 1:NT
  n(ii,:) = Nc*exp(-1*(Ec-Ef)/(kB*T(ii)/q));   % Convert kT to eV
end

% Linear scale plot
figure;
surf(Ef,T,n);
zlabel('n','FontSize',16);
shading flat;
ylabel('T (K)','FontSize',16);
xlabel('Ef (eV)','FontSize',16);

% Log scale plot
figure;
surf(Ef,T,log(n/Nc));   % Alternate log plot
zlabel('log(n/Nc)','FontSize',16);
shading flat;
ylabel('T (K)','FontSize',16);
xlabel('Ef (eV)','FontSize',16);


%break;

% =========================================================================
% Solutions 1b and 1c
clear all;
close all;

% Define variables
q = 1.6e-19;  % Electron charge coulombs
kB = 1.380e-23;  % Boltzmann's constant
hbar = 1.054e-34;  % Units J*s
m = 9.1e-31;      % Units Kg
Nc = 3.23e19;  % Si [1/cm^3]
Nv = 1.83e19;  % Si [1/cm^3]
vf = 1e6/1e2; % [cm/s], http://hyperphysics.phy-astr.gsu.edu/hbase/tables/fermi.html
li = 100e-9*1e2;  % [cm]
b0 = 5e3;  

% -------------------------------------------------------------------------
Ec = -4.0; % eV Si   
NT = 100;
Nf = 110;
Ef = -4.2;   % eV, n-type,
T = linspace(100,400,NT);  % Kelvin
n_metal = 6e23;  % Around one electron per atom [1/cm^3] for Cu, http://hyperphysics.phy-astr.gsu.edu/hbase/electric/ohmmic.html#c1
rho_metal = ((m*vf)/(n_metal*q^2))*(1/li + b0*T);
rho_si = ((m*vf)/(Nc*q^2))*(1/li + b0*T).*exp((Ec-Ef)./(kB*T/q));   % Convert m to cm

%plot(T,1./rho_metal,'b')
%plot(T,rho_metal,'b',T,rho_si,'r');
semilogy(T,rho_metal,'b',T,rho_si,'r');
legend('\rho_{e,metal}','\rho_{e,Si}');
xlabel('T (K)','FontSize',16);
ylabel('log(\rho)','FontSize',16);
grid on;

Troom = 300;  % Kelvin
rho_metal_300K = ((m*vf)/(n_metal*q^2))*(1/li + b0*Troom)  % 1c ans
rho_si_300K = ((m*vf)/(Nc*q^2))*(1/li + b0*Troom).*exp((Ec-Ef)./(kB*Troom/q))  % 1c ans

%break;

% =========================================================================
% Solutions 2b and 2c
clear all;
close all;

% -------------------------------------------------------------------------
% Define variables
q = 1.6e-19;  % Electron charge coulombs
kB = 1.380e-23;  % Boltzmann's constant
hbar = 1.054e-34;  % Units J*s
m = 9.1e-31;      % Units Kg
Nc = 3.23e19;  % Si [1/cm^3]
Nv = 1.83e19;  % Si [1/cm^3] 
T = 300;   % [K]

% -------------------------------------------------------------------------
% Solutions 2b
Nc_variable = Nc*cat(2,linspace(0.1,1,10),linspace(1,10,10));
Nv_variable = Nv*cat(2,linspace(0.1,1,10),linspace(1,10,10));
Ec = -4.0*ones(1,max(size(Nc_variable))); % eV Si 
Ev = -5.0*ones(1,max(size(Nc_variable))); % eV
Ef = 0.5*(Ec+Ev)+0.5*kB*T*log(Nv_variable/Nc)/q;

h=plot(Nv_variable/Nc,Ec,Nv_variable/Nc,Ev,Nv_variable/Nc,Ef)
axis([min(Nv_variable/Nc) max(Nv_variable/Nc) -5.5 -3.5]);
legend('Ec','Ev','Ef')
set(h,'linewidth',2.0);
set(gca,'Fontsize',16);
xlabel('Nv/Nc','FontSize',16);
ylabel('(eV)','FontSize',16);
grid;

Ef_2c = (Ec(1)+Ev(1))/2 + 0.5*kB*T*log(1/100)/q   % Ans to part 2c

%break;

% =========================================================================
% Solutions 3b and 3c
clear all;
close all;

% -------------------------------------------------------------------------
% Define variables
q = 1.6e-19;  % Electron charge coulombs
kB = 1.380e-23;  % Boltzmann's constant
hbar = 1.054e-34;  % Units J*s
m = 9.1e-31;      % Units Kg
Nc = 3.23e19;  % Si [1/cm^3]
Nv = 1.83e19;  % Si [1/cm^3] 
T = 300;   % [K]

% -------------------------------------------------------------------------
% Solutions 3b
Ec = -4.0; % eV
Ev = -5.0; % eV
Eg = 0.5*(Ec-Ev);
NT = 100;
T = linspace(100,300,NT);  % Kelvin
mu_e = 1000;  % cm^2/(V*s)
sigma_e = Nc*exp(-1*Eg./(kB*T/q))*q*mu_e;


h=plot(1./T,log(sigma_e))
set(h,'linewidth',2.0);
set(gca,'Fontsize',16);
xlabel('(1/K)','FontSize',16);
ylabel('log(\sigma)','FontSize',16);
grid;

rho_e_300K = inv(Nc*exp(-1*Eg./(kB*300/q))*q*mu_e)   % Ans 3c
