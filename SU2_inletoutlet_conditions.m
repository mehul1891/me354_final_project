clc
close
clear

% Boundary Conditions for inlet/oullet
% Inputs
Re = 22000;
M = 0.2;
%P = 101325;
T = 300;
L = 1;
gamma = 1.4;
R = 287.87;

% Reference Values
Tref        = 273.15;               %[K]
S           = 110.4;                %[-]
muref       = 1.716*10^(-5);        %[N*s/m^2]

% Equations
PtP = @(gamma,M) (1+(gamma-1)/2*M^2)^(gamma/(gamma-1));
TtT = @(gamma,M) (1+(gamma-1)/2*M^2);

% Calcualtions
mu      = muref*(T/Tref)^(3/2)*(Tref+S)/(T+S);%[N*s/m^2]
c       = sqrt(gamma*R*T);  %[m/s]
V       = M*c
rho     = Re*mu/(L*V)
P       = R*T*rho

Pt = P*PtP(gamma,M)
Tt = T*TtT(gamma,M)

% Boundary Layer aproximations
delta_laminar = 4.906*L/sqrt(Re)
delta_turbulent = .382*L/Re^(1/5)




