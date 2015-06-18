%===============================================================================
% CellML file:   D:\Desktop\Models\noble_1962.cellml
% CellML model:  noble_1962
% Date and time: 17/06/2015 at 22:26:30
%-------------------------------------------------------------------------------
% Conversion from CellML 1.0 to MATLAB (init) was done using COR (0.9.31.1409)
%    Copyright 2002-2015 Dr Alan Garny
%    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
%-------------------------------------------------------------------------------
% http://www.cellml.org/
%===============================================================================

function dY = noble_1962(time, Y)

%-------------------------------------------------------------------------------
% Initial conditions
%-------------------------------------------------------------------------------

% Y = [-87.0, 0.01, 0.8, 0.01];

% YNames = {'V', 'n', 'h', 'm'};
% YUnits = {'millivolt', 'dimensionless', 'dimensionless', 'dimensionless'};
% YComponents = {'membrane', 'potassium_channel_n_gate', 'sodium_channel_h_gate', 'sodium_channel_m_gate'};

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: V (millivolt) (in membrane)
% 2: n (dimensionless) (in potassium_channel_n_gate)
% 3: h (dimensionless) (in sodium_channel_h_gate)
% 4: m (dimensionless) (in sodium_channel_m_gate)

%-------------------------------------------------------------------------------
% Constants
%-------------------------------------------------------------------------------

E_L = -60.0;   % millivolt (in leakage_current)
g_L = 75.0;   % microS (in leakage_current)
Cm = 12.0;   % microF (in membrane)
E_Na = 40.0;   % millivolt (in sodium_channel)
g_Na_max = 400000.0;   % microS (in sodium_channel)

%-------------------------------------------------------------------------------
% Computed variables
%-------------------------------------------------------------------------------

% i_Leak (nanoA) (in leakage_current)
% alpha_n (per_second) (in potassium_channel_n_gate)
% beta_n (per_second) (in potassium_channel_n_gate)
% g_K1 (microS) (in potassium_channel)
% g_K2 (microS) (in potassium_channel)
% i_K (nanoA) (in potassium_channel)
% alpha_h (per_second) (in sodium_channel_h_gate)
% beta_h (per_second) (in sodium_channel_h_gate)
% alpha_m (per_second) (in sodium_channel_m_gate)
% beta_m (per_second) (in sodium_channel_m_gate)
% g_Na (microS) (in sodium_channel)
% i_Na (nanoA) (in sodium_channel)

%-------------------------------------------------------------------------------
% Computation
%-------------------------------------------------------------------------------

% time (second)

i_Leak = g_L*(Y(1)-E_L);
g_Na = Y(4)^3.0*Y(3)*g_Na_max;
i_Na = (g_Na+140.0)*(Y(1)-E_Na);
g_K1 = 1200.0*exp((-Y(1)-90.0)/50.0)+15.0*exp((Y(1)+90.0)/60.0);
g_K2 = 1200.0*Y(2)^4.0;
i_K = (g_K1+g_K2)*(Y(1)+100.0);
dY(1, 1) = -(i_Na+i_K+i_Leak)/Cm;
alpha_n = 0.1*(-Y(1)-50.0)/(exp((-Y(1)-50.0)/10.0)-1.0);
beta_n = 2.0*exp((-Y(1)-90.0)/80.0);
dY(2, 1) = alpha_n*(1.0-Y(2))-beta_n*Y(2);
alpha_h = 170.0*exp((-Y(1)-90.0)/20.0);
beta_h = 1000.0/(1.0+exp((-Y(1)-42.0)/10.0));
dY(3, 1) = alpha_h*(1.0-Y(3))-beta_h*Y(3);
alpha_m = 100.0*(-Y(1)-48.0)/(exp((-Y(1)-48.0)/15.0)-1.0);
beta_m = 120.0*(Y(1)+8.0)/(exp((Y(1)+8.0)/5.0)-1.0);
dY(4, 1) = alpha_m*(1.0-Y(4))-beta_m*Y(4);

%===============================================================================
% End of file
%===============================================================================
