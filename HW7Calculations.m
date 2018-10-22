% =========================================================================
% Solutions 2
clear all;
close all;

% -------------------------------------------------------------------------
% Define variables
q = 1.6e-19;  % Electron charge coulombs
kB = 1.380e-23;  % Boltzmann's constant
hbar = 1.054e-34;  % Units J*s
m = 9.1e-31;      % Units Kg
Nc = 2.5e19;  % Si [1/cm^3]
Nv = 2.5e19;  % Si [1/cm^3]
ND = 1e17;    % Donor desnity [1/cm^3] 
Ec = -4.0; % eV Si   
Ev = -5.0; % eV Si
Ed = -4.0-0.05;

% -------------------------------------------------------------------------
% Generate First plot to show zero
NE = 10000;
Ef = linspace(Ev,Ec,NE);  % Ef must lie some where between the conduction and valance band
T = 300; %300;   % K
n = Nc*exp(-1*(Ec-Ef)/(kB*T/q));
NDplus = ND./(1+2*exp(-1*(Ed-Ef)/(kB*T/q)));
p = Nv*exp(-1*(Ef-Ev)/(kB*T/q));
NeutralityFunction = n - NDplus - p;
[val,index] = min(abs(NeutralityFunction))
Ef_neutrality = Ef(index)

% -------------------------
% Linear scale plot
figure;
h=plot(Ef,NeutralityFunction);
%plot(Ef,NeutralityFunction,'b',Ef,abs(NeutralityFunction),'r');
%legend('Neutrality Function','Absolute Value of Neutrality Function')
grid on;
ylabel('Neutrality Function (/cm^3)','FontSize',16);
xlabel('Ef (eV)','FontSize',16);
set(h,'linewidth',2.0);
set(gca,'Fontsize',16)
%axis([Ev Ec -1e19 1e19]);


%break;

% -------------------------------------------------------------------------
% Plot Ef as a function of Temperature
NT = 1000;
T = linspace(1,2000,NT);
Ef_T = 0*T;
for ii=1:NT
    n = Nc*exp(-1*(Ec-Ef)/(kB*T(ii)/q));
    NDplus = ND./(1+2*exp(-1*(Ed-Ef)/(kB*T(ii)/q)));
    p = Nv*exp(-1*(Ef-Ev)/(kB*T(ii)/q));
    NeutralityFunction = n - NDplus - p;
    [val,index] = min(abs(NeutralityFunction));
    Ef_T(ii) = Ef(index);
end
    

% -------------------------
% Linear scale plot
figure;
%subplot(1,2,1);
h=plot(T,Ef_T,'b', T,Ec*ones(1,NT),'r', T,Ev*ones(1,NT),'g');
legend('Ef(T)','Ec','Ev');
grid on;
ylabel('E','FontSize',16);
xlabel('T (K)','FontSize',16);
set(h,'linewidth',2.0);
set(gca,'Fontsize',16);
axis([min(T) max(T) Ev-0.1 Ec+0.1]);

% subplot(1,2,2);
% h=plot(1./T,Ef_T);
% grid on;
% ylabel('E','FontSize',16);
% xlabel('1/T (/K)','FontSize',16);
% set(h,'linewidth',2.0);
% set(gca,'Fontsize',16);

%break;

% =========================================================================
% Solutions 3
clear all;
%close all;

% -------------------------------------------------------------------------
% Define variables
q = 1.6e-19;  % Electron charge coulombs
kB = 1.380e-23;  % Boltzmann's constant
hbar = 1.054e-34;  % Units J*s
m = 9.1e-31;      % Units Kg
Nc = 2.5e19;  % Si [1/cm^3]
Nv = 2.5e19;  % Si [1/cm^3]
ND = [1e17,5e17,1e18];    % Donor desnity [1/cm^3] 
Ec = -4.0; % eV Si   
Ev = -5.0; % eV Si
Ed = -4.0-0.05;

% -------------------------------------------------------------------------
% Plot Ef as a function of Temperature
NE = 10000;
Ef = linspace(Ev,Ec,NE);  % Ef must lie some where between the conduction and valance band
NT = 1000;
T = linspace(100,1250,NT);
Ef_T = zeros(NT,max(size(ND)));
n_T = zeros(NT,max(size(ND)));

for jj=1:max(size(ND));
    ND(jj)
    for ii=1:NT
        
        n = Nc*exp(-1*(Ec-Ef)/(kB*T(ii)/q));
        NDplus = ND(jj)./(1+2*exp(-1*(Ed-Ef)/(kB*T(ii)/q)));
        p = Nv*exp(-1*(Ef-Ev)/(kB*T(ii)/q));
        NeutralityFunction = n - NDplus - p;
        
        % -------------------------
        % Matlab trick for finding a root in a function with only one root
        [val,index] = min(abs(NeutralityFunction));   
        Ef_T(ii,jj) = Ef(index);
        
        n_T(ii,jj) = Nc*exp(-1*(Ec-Ef_T(ii,jj))/(kB*T(ii)/q));
    end
end

% -------------------------
% log scale plot
figure;
%h=semilogy(T,n_T(:,3));
h=plot(T,n_T(:,1),'b', T,n_T(:,2),'r', T,n_T(:,3),'g');
legend('1e17','5e17','1e18')
grid on;
ylabel('ln(n)','FontSize',16);
xlabel('T (K)','FontSize',16);
set(h,'linewidth',2.0);
set(gca,'Fontsize',16);

