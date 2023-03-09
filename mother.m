function dydt = mother(t,y)
%Units: time(day) length(mm)

V=y(1);  % Viral load
E=y(2);  % virus-specific immune effector cells
RI=y(3); % actively-infected cells
RS=y(4); % susceptible cells
RL=y(5); % latently-infected cells

%parameters
rho=5; % Virion induced immune response cells/(virions*day)
n=40; %Productivity of infected cell
delta=0.4; %Rate of viral induced cell death
c=0.3; %Rate of viral clearance
k=1e-4; %Infection rate constant
m=0.1; %Immune induced cell lysis
a0=0.2; %Exit and reactivation rate for monocytes
kappa=2e-3; %rate of latency
lambdaE=4e-2; %Homeostatic replenishment of immune cells
lambda=1e-3; %Cell replenishment rate
effector=9; %HCMV-specific effector cell term
epsilonS=0; %Level of immune suppression
f=1; %Number of infecting virions per cell
rS=4e2; %Equilibrium level of susceptible cells
rL=4e-2; %Equilibrium level of latent cells
tH=2; %Half-life of virions during antiviral treatment
tD=1.5; %Doubling time of the virions
Ee=10; %Equilibrium level of HCMV-specific effector cells
Ve=1e-2; %Equilibrium level of virions


Vderivative=n*delta*RI-c*V-f*k*RS*V;
Ederivative=(1-epsilonS)*(lambdaE*(1-E/effector)*E + rho *V);
RIderivative=k*RS*V-delta*RI-(1-epsilonS)*m*E*RI+a0*RL-kappa*RI;
RSderivative=lambda*(1-RS/rS)*RS-k*RS*V;
RLderivative=lambda*(1-RL/rL)*RL-a0*RL+kappa*RI;

dydt = [Vderivative; Ederivative; RIderivative; RSderivative; RLderivative];

end