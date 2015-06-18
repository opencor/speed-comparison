%===============================================================================
% CellML file:   D:\Desktop\Models\garny_kohl_hunter_boyett_noble_2003.cellml
% CellML model:  garny_2003
% Date and time: 17/06/2015 at 21:18:05
%-------------------------------------------------------------------------------
% Conversion from CellML 1.0 to MATLAB (init) was done using COR (0.9.31.1409)
%    Copyright 2002-2015 Dr Alan Garny
%    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
%-------------------------------------------------------------------------------
% http://www.cellml.org/
%===============================================================================

function dY = garny_kohl_hunter_boyett_noble_2003(time, Y)

%-------------------------------------------------------------------------------
% Initial conditions
%-------------------------------------------------------------------------------

% Y = [0.04804900895, 0.48779845203, 0.42074047435, 0.038968420558, 0.29760539675, 0.064402950262, 0.03889291759, -39.013558536, 0.13034201158, 0.46960956028, 0.87993375273, 0.082293827208, 0.015905380261, 0.01445216109, 0.092361701692];

% YNames = {'d_L', 'f_L', 'd_T', 'f_T', 'q', 'r', 'y', 'V', 'P_af', 'P_as', 'P_i', 'xs', 'h1', 'h2', 'm'};
% YUnits = {'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'millivolt', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless'};
% YComponents = {'L_type_Ca_channel_d_gate', 'L_type_Ca_channel_f_gate', 'T_type_Ca_channel_d_gate', 'T_type_Ca_channel_f_gate', 'four_AP_sensitive_currents_q_gate', 'four_AP_sensitive_currents_r_gate', 'hyperpolarisation_activated_current_y_gate', 'membrane', 'rapid_delayed_rectifying_potassium_current_P_af_gate', 'rapid_delayed_rectifying_potassium_current_P_as_gate', 'rapid_delayed_rectifying_potassium_current_P_i_gate', 'slow_delayed_rectifying_potassium_current_xs_gate', 'sodium_current_h_gate', 'sodium_current_h_gate', 'sodium_current_m_gate'};

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: d_L (dimensionless) (in L_type_Ca_channel_d_gate)
% 2: f_L (dimensionless) (in L_type_Ca_channel_f_gate)
% 3: d_T (dimensionless) (in T_type_Ca_channel_d_gate)
% 4: f_T (dimensionless) (in T_type_Ca_channel_f_gate)
% 5: q (dimensionless) (in four_AP_sensitive_currents_q_gate)
% 6: r (dimensionless) (in four_AP_sensitive_currents_r_gate)
% 7: y (dimensionless) (in hyperpolarisation_activated_current_y_gate)
% 8: V (millivolt) (in membrane)
% 9: P_af (dimensionless) (in rapid_delayed_rectifying_potassium_current_P_af_gate)
% 10: P_as (dimensionless) (in rapid_delayed_rectifying_potassium_current_P_as_gate)
% 11: P_i (dimensionless) (in rapid_delayed_rectifying_potassium_current_P_i_gate)
% 12: xs (dimensionless) (in slow_delayed_rectifying_potassium_current_xs_gate)
% 13: h1 (dimensionless) (in sodium_current_h_gate)
% 14: h2 (dimensionless) (in sodium_current_h_gate)
% 15: m (dimensionless) (in sodium_current_m_gate)

%-------------------------------------------------------------------------------
% Constants
%-------------------------------------------------------------------------------

E_Ca_L = 46.4;   % millivolt (in L_type_Ca_channel)
g_Ca_L_Centre_0DCapable = 0.0057938;   % microS (in L_type_Ca_channel)
g_Ca_L_Centre_1DCapable = 0.0082;   % microS (in L_type_Ca_channel)
g_Ca_L_Centre_Published = 0.0058;   % microS (in L_type_Ca_channel)
g_Ca_L_Periphery_0DCapable = 0.06588648;   % microS (in L_type_Ca_channel)
g_Ca_L_Periphery_1DCapable = 0.0659;   % microS (in L_type_Ca_channel)
g_Ca_L_Periphery_Published = 0.0659;   % microS (in L_type_Ca_channel)
E_Ca_T = 45.0;   % millivolt (in T_type_Ca_channel)
g_Ca_T_Centre_0DCapable = 0.00427806;   % microS (in T_type_Ca_channel)
g_Ca_T_Centre_1DCapable = 0.0021;   % microS (in T_type_Ca_channel)
g_Ca_T_Centre_Published = 0.0043;   % microS (in T_type_Ca_channel)
g_Ca_T_Periphery_0DCapable = 0.0138823;   % microS (in T_type_Ca_channel)
g_Ca_T_Periphery_1DCapable = 0.00694;   % microS (in T_type_Ca_channel)
g_Ca_T_Periphery_Published = 0.0139;   % microS (in T_type_Ca_channel)
g_b_Ca_Centre_0DCapable = 1.3236e-5;   % microS (in calcium_background_current)
g_b_Ca_Centre_1DCapable = 1.323e-5;   % microS (in calcium_background_current)
g_b_Ca_Centre_Published = 1.32e-5;   % microS (in calcium_background_current)
g_b_Ca_Periphery_0DCapable = 4.2952e-5;   % microS (in calcium_background_current)
g_b_Ca_Periphery_1DCapable = 4.29e-5;   % microS (in calcium_background_current)
g_b_Ca_Periphery_Published = 4.3e-5;   % microS (in calcium_background_current)
g_sus_Centre_0DCapable = 6.645504e-5;   % microS (in four_AP_sensitive_currents)
g_sus_Centre_1DCapable = 0.000266;   % microS (in four_AP_sensitive_currents)
g_sus_Centre_Published = 6.65e-5;   % microS (in four_AP_sensitive_currents)
g_sus_Periphery_0DCapable = 0.01138376;   % microS (in four_AP_sensitive_currents)
g_sus_Periphery_1DCapable = 0.0114;   % microS (in four_AP_sensitive_currents)
g_sus_Periphery_Published = 0.0114;   % microS (in four_AP_sensitive_currents)
g_to_Centre_0DCapable = 0.004905;   % microS (in four_AP_sensitive_currents)
g_to_Centre_1DCapable = 0.004905;   % microS (in four_AP_sensitive_currents)
g_to_Centre_Published = 0.00491;   % microS (in four_AP_sensitive_currents)
g_to_Periphery_0DCapable = 0.036495;   % microS (in four_AP_sensitive_currents)
g_to_Periphery_1DCapable = 0.0365;   % microS (in four_AP_sensitive_currents)
g_to_Periphery_Published = 0.03649;   % microS (in four_AP_sensitive_currents)
g_f_K_Centre_0DCapable = 0.0005465;   % microS (in hyperpolarisation_activated_current)
g_f_K_Centre_1DCapable = 0.000437;   % microS (in hyperpolarisation_activated_current)
g_f_K_Centre_Published = 0.000548;   % microS (in hyperpolarisation_activated_current)
g_f_K_Periphery_0DCapable = 0.006875;   % microS (in hyperpolarisation_activated_current)
g_f_K_Periphery_1DCapable = 0.0055;   % microS (in hyperpolarisation_activated_current)
g_f_K_Periphery_Published = 0.0069;   % microS (in hyperpolarisation_activated_current)
g_f_Na_Centre_0DCapable = 0.0005465;   % microS (in hyperpolarisation_activated_current)
g_f_Na_Centre_1DCapable = 0.000437;   % microS (in hyperpolarisation_activated_current)
g_f_Na_Centre_Published = 0.000548;   % microS (in hyperpolarisation_activated_current)
g_f_Na_Periphery_0DCapable = 0.006875;   % microS (in hyperpolarisation_activated_current)
g_f_Na_Periphery_1DCapable = 0.0055;   % microS (in hyperpolarisation_activated_current)
g_f_Na_Periphery_Published = 0.0069;   % microS (in hyperpolarisation_activated_current)
Ca_i = 0.0001;   % millimolar (in ionic_concentrations)
Ca_o = 2.0;   % millimolar (in ionic_concentrations)
K_i = 140.0;   % millimolar (in ionic_concentrations)
K_o = 5.4;   % millimolar (in ionic_concentrations)
Na_i = 8.0;   % millimolar (in ionic_concentrations)
Na_o = 140.0;   % millimolar (in ionic_concentrations)
CmCentre = 2.0e-5;   % microF (in membrane)
CmPeriphery = 6.5e-5;   % microF (in membrane)
F = 96845.0;   % coulomb_per_mole (in membrane)
FCellConstant = 1.0309347;   % dimensionless (in membrane)
R = 8314.0;   % millijoule_per_mole_kelvin (in membrane)
T = 310.0;   % kelvin (in membrane)
Version = 1.0;   % dimensionless (in membrane)
dCell = 0.0;   % dimensionless (in membrane)
i_Ca_p_max_Centre_0DCapable = 0.0;   % nanoA (in persistent_calcium_current)
i_Ca_p_max_Centre_1DCapable = 0.0042;   % nanoA (in persistent_calcium_current)
i_Ca_p_max_Centre_Published = 0.0;   % nanoA (in persistent_calcium_current)
i_Ca_p_max_Periphery_0DCapable = 0.0;   % nanoA (in persistent_calcium_current)
i_Ca_p_max_Periphery_1DCapable = 0.03339;   % nanoA (in persistent_calcium_current)
i_Ca_p_max_Periphery_Published = 0.0;   % nanoA (in persistent_calcium_current)
g_b_K_Centre_0DCapable = 2.523636e-5;   % microS (in potassium_background_current)
g_b_K_Centre_1DCapable = 2.52e-5;   % microS (in potassium_background_current)
g_b_K_Centre_Published = 2.52e-5;   % microS (in potassium_background_current)
g_b_K_Periphery_0DCapable = 8.1892e-5;   % microS (in potassium_background_current)
g_b_K_Periphery_1DCapable = 8.19e-5;   % microS (in potassium_background_current)
g_b_K_Periphery_Published = 8.19e-5;   % microS (in potassium_background_current)
g_K_r_Centre_0DCapable = 0.00079704;   % microS (in rapid_delayed_rectifying_potassium_current)
g_K_r_Centre_1DCapable = 0.000738;   % microS (in rapid_delayed_rectifying_potassium_current)
g_K_r_Centre_Published = 0.000797;   % microS (in rapid_delayed_rectifying_potassium_current)
g_K_r_Periphery_0DCapable = 0.016;   % microS (in rapid_delayed_rectifying_potassium_current)
g_K_r_Periphery_1DCapable = 0.0208;   % microS (in rapid_delayed_rectifying_potassium_current)
g_K_r_Periphery_Published = 0.016;   % microS (in rapid_delayed_rectifying_potassium_current)
g_K_s_Centre_0DCapable = 0.0003445;   % microS (in slow_delayed_rectifying_potassium_current)
g_K_s_Centre_1DCapable = 0.000345;   % microS (in slow_delayed_rectifying_potassium_current)
g_K_s_Centre_Published = 0.000518;   % microS (in slow_delayed_rectifying_potassium_current)
g_K_s_Periphery_0DCapable = 0.0104;   % microS (in slow_delayed_rectifying_potassium_current)
g_K_s_Periphery_1DCapable = 0.0104;   % microS (in slow_delayed_rectifying_potassium_current)
g_K_s_Periphery_Published = 0.0104;   % microS (in slow_delayed_rectifying_potassium_current)
g_b_Na_Centre_0DCapable = 5.81818e-5;   % microS (in sodium_background_current)
g_b_Na_Centre_1DCapable = 5.8e-5;   % microS (in sodium_background_current)
g_b_Na_Centre_Published = 5.8e-5;   % microS (in sodium_background_current)
g_b_Na_Periphery_0DCapable = 0.0001888;   % microS (in sodium_background_current)
g_b_Na_Periphery_1DCapable = 0.000189;   % microS (in sodium_background_current)
g_b_Na_Periphery_Published = 0.000189;   % microS (in sodium_background_current)
d_NaCa = 0.0001;   % dimensionless (in sodium_calcium_exchanger)
gamma_NaCa = 0.5;   % dimensionless (in sodium_calcium_exchanger)
k_NaCa_Centre_0DCapable = 2.7229e-6;   % nanoA (in sodium_calcium_exchanger)
k_NaCa_Centre_1DCapable = 2.8e-6;   % nanoA (in sodium_calcium_exchanger)
k_NaCa_Centre_Published = 2.7e-6;   % nanoA (in sodium_calcium_exchanger)
k_NaCa_Periphery_0DCapable = 8.83584e-6;   % nanoA (in sodium_calcium_exchanger)
k_NaCa_Periphery_1DCapable = 8.8e-6;   % nanoA (in sodium_calcium_exchanger)
k_NaCa_Periphery_Published = 8.8e-6;   % nanoA (in sodium_calcium_exchanger)
g_Na_Centre_0DCapable = 0.0;   % microlitre_per_second (in sodium_current)
g_Na_Centre_1DCapable = 0.0;   % microlitre_per_second (in sodium_current)
g_Na_Centre_Published = 0.0;   % microlitre_per_second (in sodium_current)
g_Na_Periphery_0DCapable = 1.204e-6;   % microlitre_per_second (in sodium_current)
g_Na_Periphery_1DCapable = 3.7e-7;   % microlitre_per_second (in sodium_current)
g_Na_Periphery_Published = 1.2e-6;   % microlitre_per_second (in sodium_current)
K_m_K = 0.621;   % millimolar (in sodium_potassium_pump)
K_m_Na = 5.64;   % millimolar (in sodium_potassium_pump)
i_p_max_Centre_0DCapable = 0.04782545;   % nanoA (in sodium_potassium_pump)
i_p_max_Centre_1DCapable = 0.0478;   % nanoA (in sodium_potassium_pump)
i_p_max_Centre_Published = 0.0478;   % nanoA (in sodium_potassium_pump)
i_p_max_Periphery_0DCapable = 0.1551936;   % nanoA (in sodium_potassium_pump)
i_p_max_Periphery_1DCapable = 0.16;   % nanoA (in sodium_potassium_pump)
i_p_max_Periphery_Published = 0.16;   % nanoA (in sodium_potassium_pump)

%-------------------------------------------------------------------------------
% Computed variables
%-------------------------------------------------------------------------------

% alpha_d_L (per_second) (in L_type_Ca_channel_d_gate)
% beta_d_L (per_second) (in L_type_Ca_channel_d_gate)
% d_L_infinity (dimensionless) (in L_type_Ca_channel_d_gate)
% tau_d_L (second) (in L_type_Ca_channel_d_gate)
% alpha_f_L (per_second) (in L_type_Ca_channel_f_gate)
% beta_f_L (per_second) (in L_type_Ca_channel_f_gate)
% f_L_infinity (dimensionless) (in L_type_Ca_channel_f_gate)
% tau_f_L (second) (in L_type_Ca_channel_f_gate)
% g_Ca_L (microS) (in L_type_Ca_channel)
% i_Ca_L (nanoA) (in L_type_Ca_channel)
% alpha_d_T (per_second) (in T_type_Ca_channel_d_gate)
% beta_d_T (per_second) (in T_type_Ca_channel_d_gate)
% d_T_infinity (dimensionless) (in T_type_Ca_channel_d_gate)
% tau_d_T (second) (in T_type_Ca_channel_d_gate)
% alpha_f_T (per_second) (in T_type_Ca_channel_f_gate)
% beta_f_T (per_second) (in T_type_Ca_channel_f_gate)
% f_T_infinity (dimensionless) (in T_type_Ca_channel_f_gate)
% tau_f_T (second) (in T_type_Ca_channel_f_gate)
% g_Ca_T (microS) (in T_type_Ca_channel)
% i_Ca_T (nanoA) (in T_type_Ca_channel)
% g_b_Ca (microS) (in calcium_background_current)
% i_b_Ca (nanoA) (in calcium_background_current)
% q_infinity (dimensionless) (in four_AP_sensitive_currents_q_gate)
% tau_q (second) (in four_AP_sensitive_currents_q_gate)
% r_infinity (dimensionless) (in four_AP_sensitive_currents_r_gate)
% tau_r (second) (in four_AP_sensitive_currents_r_gate)
% g_sus (microS) (in four_AP_sensitive_currents)
% g_to (microS) (in four_AP_sensitive_currents)
% i_sus (nanoA) (in four_AP_sensitive_currents)
% i_to (nanoA) (in four_AP_sensitive_currents)
% alpha_y (per_second) (in hyperpolarisation_activated_current_y_gate)
% beta_y (per_second) (in hyperpolarisation_activated_current_y_gate)
% g_f_K (microS) (in hyperpolarisation_activated_current)
% g_f_Na (microS) (in hyperpolarisation_activated_current)
% i_f_K (nanoA) (in hyperpolarisation_activated_current)
% i_f_Na (nanoA) (in hyperpolarisation_activated_current)
% Cm (microF) (in membrane)
% FCell (dimensionless) (in membrane)
% i_Ca_p (nanoA) (in persistent_calcium_current)
% i_Ca_p_max (nanoA) (in persistent_calcium_current)
% g_b_K (microS) (in potassium_background_current)
% i_b_K (nanoA) (in potassium_background_current)
% P_af_infinity (dimensionless) (in rapid_delayed_rectifying_potassium_current_P_af_gate)
% tau_P_af (second) (in rapid_delayed_rectifying_potassium_current_P_af_gate)
% P_as_infinity (dimensionless) (in rapid_delayed_rectifying_potassium_current_P_as_gate)
% tau_P_as (second) (in rapid_delayed_rectifying_potassium_current_P_as_gate)
% P_i_infinity (dimensionless) (in rapid_delayed_rectifying_potassium_current_P_i_gate)
% tau_P_i (second) (in rapid_delayed_rectifying_potassium_current_P_i_gate)
% P_a (dimensionless) (in rapid_delayed_rectifying_potassium_current)
% g_K_r (microS) (in rapid_delayed_rectifying_potassium_current)
% i_K_r (nanoA) (in rapid_delayed_rectifying_potassium_current)
% E_Ca (millivolt) (in reversal_and_equilibrium_potentials)
% E_K (millivolt) (in reversal_and_equilibrium_potentials)
% E_K_s (millivolt) (in reversal_and_equilibrium_potentials)
% E_Na (millivolt) (in reversal_and_equilibrium_potentials)
% alpha_xs (per_second) (in slow_delayed_rectifying_potassium_current_xs_gate)
% beta_xs (per_second) (in slow_delayed_rectifying_potassium_current_xs_gate)
% g_K_s (microS) (in slow_delayed_rectifying_potassium_current)
% i_K_s (nanoA) (in slow_delayed_rectifying_potassium_current)
% g_b_Na (microS) (in sodium_background_current)
% i_b_Na (nanoA) (in sodium_background_current)
% i_NaCa (nanoA) (in sodium_calcium_exchanger)
% k_NaCa (nanoA) (in sodium_calcium_exchanger)
% F_Na (dimensionless) (in sodium_current_h_gate)
% h (dimensionless) (in sodium_current_h_gate)
% h1_infinity (dimensionless) (in sodium_current_h_gate)
% h2_infinity (dimensionless) (in sodium_current_h_gate)
% tau_h1 (second) (in sodium_current_h_gate)
% tau_h2 (second) (in sodium_current_h_gate)
% m_infinity (dimensionless) (in sodium_current_m_gate)
% tau_m (second) (in sodium_current_m_gate)
% g_Na (microlitre_per_second) (in sodium_current)
% i_Na (nanoA) (in sodium_current)
% i_p (nanoA) (in sodium_potassium_pump)
% i_p_max (nanoA) (in sodium_potassium_pump)

%-------------------------------------------------------------------------------
% Computation
%-------------------------------------------------------------------------------

% time (second)

if (Version == 0.0)
   FCell = 1.07*(3.0*dCell-0.1)/(3.0*(1.0+0.7745*exp(-(3.0*dCell-2.05)/0.295)));
elseif (Version == 1.0)
   FCell = FCellConstant*dCell/(1.0+0.7745*exp(-(3.0*dCell-2.05)/0.295));
else
   FCell = 1.07*29.0*dCell/(30.0*(1.0+0.7745*exp(-(29.0*dCell-24.5)/1.95)));
end;

if (Version == 0.0)
   g_Ca_L = g_Ca_L_Centre_Published+FCell*(g_Ca_L_Periphery_Published-g_Ca_L_Centre_Published);
elseif (Version == 1.0)
   g_Ca_L = g_Ca_L_Centre_0DCapable+FCell*(g_Ca_L_Periphery_0DCapable-g_Ca_L_Centre_0DCapable);
else
   g_Ca_L = g_Ca_L_Centre_1DCapable+FCell*(g_Ca_L_Periphery_1DCapable-g_Ca_L_Centre_1DCapable);
end;

i_Ca_L = g_Ca_L*(Y(2)*Y(1)+0.006/(1.0+exp(-(Y(8)+14.1)/6.0)))*(Y(8)-E_Ca_L);

if (Version == 0.0)
   d_L_infinity = 1.0/(1.0+exp(-(Y(8)+23.1)/6.0));
elseif (Version == 1.0)
   d_L_infinity = 1.0/(1.0+exp(-(Y(8)+22.3+0.8*FCell)/6.0));
else
   d_L_infinity = 1.0/(1.0+exp(-(Y(8)+22.2)/6.0));
end;

if (Version == 0.0)
   alpha_d_L = -28.38*(Y(8)+35.0)/(exp(-(Y(8)+35.0)/2.5)-1.0)-84.9*Y(8)/(exp(-0.208*Y(8))-1.0);
elseif (Version == 1.0)
   alpha_d_L = -28.39*(Y(8)+35.0)/(exp(-(Y(8)+35.0)/2.5)-1.0)-84.9*Y(8)/(exp(-0.208*Y(8))-1.0);
else
   alpha_d_L = -28.4*(Y(8)+35.0)/(exp(-(Y(8)+35.0)/2.5)-1.0)-84.9*Y(8)/(exp(-0.208*Y(8))-1.0);
end;

if (Version == 1.0)
   beta_d_L = 11.43*(Y(8)-5.0)/(exp(0.4*(Y(8)-5.0))-1.0);
else
   beta_d_L = 11.42*(Y(8)-5.0)/(exp(0.4*(Y(8)-5.0))-1.0);
end;

tau_d_L = 2.0/(alpha_d_L+beta_d_L);
dY(1, 1) = (d_L_infinity-Y(1))/tau_d_L;
f_L_infinity = 1.0/(1.0+exp((Y(8)+45.0)/5.0));

if (Version == 1.0)
   alpha_f_L = 3.75*(Y(8)+28.0)/(exp((Y(8)+28.0)/4.0)-1.0);
else
   alpha_f_L = 3.12*(Y(8)+28.0)/(exp((Y(8)+28.0)/4.0)-1.0);
end;

if (Version == 1.0)
   beta_f_L = 30.0/(1.0+exp(-(Y(8)+28.0)/4.0));
else
   beta_f_L = 25.0/(1.0+exp(-(Y(8)+28.0)/4.0));
end;

if (Version == 1.0)
   tau_f_L = (1.2-0.2*FCell)/(alpha_f_L+beta_f_L);
else
   tau_f_L = 1.0/(alpha_f_L+beta_f_L);
end;

dY(2, 1) = (f_L_infinity-Y(2))/tau_f_L;

if (Version == 0.0)
   g_Ca_T = g_Ca_T_Centre_Published+FCell*(g_Ca_T_Periphery_Published-g_Ca_T_Centre_Published);
elseif (Version == 1.0)
   g_Ca_T = g_Ca_T_Centre_0DCapable+FCell*(g_Ca_T_Periphery_0DCapable-g_Ca_T_Centre_0DCapable);
else
   g_Ca_T = g_Ca_T_Centre_1DCapable+FCell*(g_Ca_T_Periphery_1DCapable-g_Ca_T_Centre_1DCapable);
end;

i_Ca_T = g_Ca_T*Y(3)*Y(4)*(Y(8)-E_Ca_T);
d_T_infinity = 1.0/(1.0+exp(-(Y(8)+37.0)/6.8));
alpha_d_T = 1068.0*exp((Y(8)+26.3)/30.0);
beta_d_T = 1068.0*exp(-(Y(8)+26.3)/30.0);
tau_d_T = 1.0/(alpha_d_T+beta_d_T);
dY(3, 1) = (d_T_infinity-Y(3))/tau_d_T;
f_T_infinity = 1.0/(1.0+exp((Y(8)+71.0)/9.0));

if (Version == 1.0)
   alpha_f_T = 15.3*exp(-(Y(8)+71.0+0.7*FCell)/83.3);
else
   alpha_f_T = 15.3*exp(-(Y(8)+71.7)/83.3);
end;

if (Version == 1.0)
   beta_f_T = 15.0*exp((Y(8)+71.0)/15.38);
else
   beta_f_T = 15.0*exp((Y(8)+71.7)/15.38);
end;

tau_f_T = 1.0/(alpha_f_T+beta_f_T);
dY(4, 1) = (f_T_infinity-Y(4))/tau_f_T;

if (Version == 0.0)
   g_b_Ca = g_b_Ca_Centre_Published+FCell*(g_b_Ca_Periphery_Published-g_b_Ca_Centre_Published);
elseif (Version == 1.0)
   g_b_Ca = g_b_Ca_Centre_0DCapable+FCell*(g_b_Ca_Periphery_0DCapable-g_b_Ca_Centre_0DCapable);
else
   g_b_Ca = g_b_Ca_Centre_1DCapable+FCell*(g_b_Ca_Periphery_1DCapable-g_b_Ca_Centre_1DCapable);
end;

E_Ca = R*T/(2.0*F)*log(Ca_o/Ca_i);
i_b_Ca = g_b_Ca*(Y(8)-E_Ca);

if (Version == 0.0)
   g_to = g_to_Centre_Published+FCell*(g_to_Periphery_Published-g_to_Centre_Published);
elseif (Version == 1.0)
   g_to = g_to_Centre_0DCapable+FCell*(g_to_Periphery_0DCapable-g_to_Centre_0DCapable);
else
   g_to = g_to_Centre_1DCapable+FCell*(g_to_Periphery_1DCapable-g_to_Centre_1DCapable);
end;

if (Version == 0.0)
   g_sus = g_sus_Centre_Published+FCell*(g_sus_Periphery_Published-g_sus_Centre_Published);
elseif (Version == 1.0)
   g_sus = g_sus_Centre_0DCapable+FCell*(g_sus_Periphery_0DCapable-g_sus_Centre_0DCapable);
else
   g_sus = g_sus_Centre_1DCapable+FCell*(g_sus_Periphery_1DCapable-g_sus_Centre_1DCapable);
end;

E_K = R*T/F*log(K_o/K_i);
i_to = g_to*Y(5)*Y(6)*(Y(8)-E_K);
i_sus = g_sus*Y(6)*(Y(8)-E_K);
q_infinity = 1.0/(1.0+exp((Y(8)+59.37)/13.1));

if (Version == 0.0)
   tau_q = 0.0101+0.06517/(0.57*exp(-0.08*(Y(8)+49.0)))+2.4e-5*exp(0.1*(Y(8)+50.93));
elseif (Version == 1.0)
   tau_q = 0.001/3.0*(30.31+195.5/(0.5686*exp(-0.08161*(Y(8)+39.0+10.0*FCell))+0.7174*exp((0.2719-0.1719*FCell)*1.0*(Y(8)+40.93+10.0*FCell))));
else
   tau_q = 0.0101+0.06517/(0.5686*exp(-0.08161*(Y(8)+39.0))+0.7174*exp(0.2719*(Y(8)+40.93)));
end;

dY(5, 1) = (q_infinity-Y(5))/tau_q;
r_infinity = 1.0/(1.0+exp(-(Y(8)-10.93)/19.7));

if (Version == 0.0)
   tau_r = 0.001*(2.98+15.59/(1.037*exp(0.09*(Y(8)+30.61))+0.369*exp(-0.12*(Y(8)+23.84))));
elseif (Version == 1.0)
   tau_r = 0.0025*(1.191+7.838/(1.037*exp(0.09012*(Y(8)+30.61))+0.369*exp(-0.119*(Y(8)+23.84))));
else
   tau_r = 0.001*(2.98+19.59/(1.037*exp(0.09012*(Y(8)+30.61))+0.369*exp(-0.119*(Y(8)+23.84))));
end;

dY(6, 1) = (r_infinity-Y(6))/tau_r;

if (Version == 0.0)
   g_f_Na = g_f_Na_Centre_Published+FCell*(g_f_Na_Periphery_Published-g_f_Na_Centre_Published);
elseif (Version == 1.0)
   g_f_Na = g_f_Na_Centre_0DCapable+FCell*(g_f_Na_Periphery_0DCapable-g_f_Na_Centre_0DCapable);
else
   g_f_Na = g_f_Na_Centre_1DCapable+FCell*(g_f_Na_Periphery_1DCapable-g_f_Na_Centre_1DCapable);
end;

E_Na = R*T/F*log(Na_o/Na_i);

if (Version ~= 2.0)
   i_f_Na = g_f_Na*Y(7)*(Y(8)-E_Na);
else
   i_f_Na = g_f_Na*Y(7)*(Y(8)-77.6);
end;

if (Version == 0.0)
   g_f_K = g_f_K_Centre_Published+FCell*(g_f_K_Periphery_Published-g_f_K_Centre_Published);
elseif (Version == 1.0)
   g_f_K = g_f_K_Centre_0DCapable+FCell*(g_f_K_Periphery_0DCapable-g_f_K_Centre_0DCapable);
else
   g_f_K = g_f_K_Centre_1DCapable+FCell*(g_f_K_Periphery_1DCapable-g_f_K_Centre_1DCapable);
end;

if (Version ~= 2.0)
   i_f_K = g_f_K*Y(7)*(Y(8)-E_K);
else
   i_f_K = g_f_K*Y(7)*(Y(8)+102.0);
end;

if (Version == 0.0)
   alpha_y = 1.0*exp(-(Y(8)+78.91)/26.62);
else
   alpha_y = 1.0*exp(-(Y(8)+78.91)/26.63);
end;

beta_y = 1.0*exp((Y(8)+75.13)/21.25);
dY(7, 1) = alpha_y*(1.0-Y(7))-beta_y*Y(7);
Cm = CmCentre+FCell*(CmPeriphery-CmCentre);

if (Version == 0.0)
   g_Na = g_Na_Centre_Published+FCell*(g_Na_Periphery_Published-g_Na_Centre_Published);
elseif (Version == 1.0)
   g_Na = g_Na_Centre_0DCapable+FCell*(g_Na_Periphery_0DCapable-g_Na_Centre_0DCapable);
else
   g_Na = g_Na_Centre_1DCapable+FCell*(g_Na_Periphery_1DCapable-g_Na_Centre_1DCapable);
end;

if (Version == 0.0)
   F_Na = 0.0952*exp(-0.063*(Y(8)+34.4))/(1.0+1.66*exp(-0.225*(Y(8)+63.7)))+0.0869;
else
   F_Na = 0.09518*exp(-0.06306*(Y(8)+34.4))/(1.0+1.662*exp(-0.2251*(Y(8)+63.7)))+0.08693;
end;

h = (1.0-F_Na)*Y(13)+F_Na*Y(14);
i_Na = g_Na*Y(15)^3.0*h*Na_o*F^2.0/(R*T)*(exp((Y(8)-E_Na)*F/(R*T))-1.0)/(exp(Y(8)*F/(R*T))-1.0)*Y(8);

if (Version == 0.0)
   g_K_r = g_K_r_Centre_Published+FCell*(g_K_r_Periphery_Published-g_K_r_Centre_Published);
elseif (Version == 1.0)
   g_K_r = g_K_r_Centre_0DCapable+FCell*(g_K_r_Periphery_0DCapable-g_K_r_Centre_0DCapable);
else
   g_K_r = g_K_r_Centre_1DCapable+FCell*(g_K_r_Periphery_1DCapable-g_K_r_Centre_1DCapable);
end;

P_a = 0.6*Y(9)+0.4*Y(10);
i_K_r = g_K_r*P_a*Y(11)*(Y(8)-E_K);

if (Version == 0.0)
   g_K_s = g_K_s_Centre_Published+FCell*(g_K_s_Periphery_Published-g_K_s_Centre_Published);
elseif (Version == 1.0)
   g_K_s = g_K_s_Centre_0DCapable+FCell*(g_K_s_Periphery_0DCapable-g_K_s_Centre_0DCapable);
else
   g_K_s = g_K_s_Centre_1DCapable+FCell*(g_K_s_Periphery_1DCapable-g_K_s_Centre_1DCapable);
end;

if (Version == 0.0)
   E_K_s = R*T/F*log((K_o+0.12*Na_o)/(K_i+0.12*Na_i));
else
   E_K_s = R*T/F*log((K_o+0.03*Na_o)/(K_i+0.03*Na_i));
end;

i_K_s = g_K_s*Y(12)^2.0*(Y(8)-E_K_s);

if (Version == 0.0)
   g_b_Na = g_b_Na_Centre_Published+FCell*(g_b_Na_Periphery_Published-g_b_Na_Centre_Published);
elseif (Version == 1.0)
   g_b_Na = g_b_Na_Centre_0DCapable+FCell*(g_b_Na_Periphery_0DCapable-g_b_Na_Centre_0DCapable);
else
   g_b_Na = g_b_Na_Centre_1DCapable+FCell*(g_b_Na_Periphery_1DCapable-g_b_Na_Centre_1DCapable);
end;

i_b_Na = g_b_Na*(Y(8)-E_Na);

if (Version == 0.0)
   g_b_K = g_b_K_Centre_Published+FCell*(g_b_K_Periphery_Published-g_b_K_Centre_Published);
elseif (Version == 1.0)
   g_b_K = g_b_K_Centre_0DCapable+FCell*(g_b_K_Periphery_0DCapable-g_b_K_Centre_0DCapable);
else
   g_b_K = g_b_K_Centre_1DCapable+FCell*(g_b_K_Periphery_1DCapable-g_b_K_Centre_1DCapable);
end;

i_b_K = g_b_K*(Y(8)-E_K);

if (Version == 0.0)
   k_NaCa = k_NaCa_Centre_Published+FCell*(k_NaCa_Periphery_Published-k_NaCa_Centre_Published);
elseif (Version == 1.0)
   k_NaCa = k_NaCa_Centre_0DCapable+FCell*(k_NaCa_Periphery_0DCapable-k_NaCa_Centre_0DCapable);
else
   k_NaCa = k_NaCa_Centre_1DCapable+FCell*(k_NaCa_Periphery_1DCapable-k_NaCa_Centre_1DCapable);
end;

if (Version == 0.0)
   i_NaCa = k_NaCa*(Na_i^3.0*Ca_o*exp(0.03743*Y(8)*gamma_NaCa)-Na_o^3.0*Ca_i*exp(0.0374*Y(8)*(gamma_NaCa-1.0)))/(1.0+d_NaCa*(Ca_i*Na_o^3.0+Ca_o*Na_i^3.0));
else
   i_NaCa = k_NaCa*(Na_i^3.0*Ca_o*exp(0.03743*Y(8)*gamma_NaCa)-Na_o^3.0*Ca_i*exp(0.03743*Y(8)*(gamma_NaCa-1.0)))/(1.0+d_NaCa*(Ca_i*Na_o^3.0+Ca_o*Na_i^3.0));
end;

if (Version == 0.0)
   i_p_max = i_p_max_Centre_Published+FCell*(i_p_max_Periphery_Published-i_p_max_Centre_Published);
elseif (Version == 1.0)
   i_p_max = i_p_max_Centre_0DCapable+FCell*(i_p_max_Periphery_0DCapable-i_p_max_Centre_0DCapable);
else
   i_p_max = i_p_max_Centre_1DCapable+FCell*(i_p_max_Periphery_1DCapable-i_p_max_Centre_1DCapable);
end;

i_p = i_p_max*(Na_i/(K_m_Na+Na_i))^3.0*(K_o/(K_m_K+K_o))^2.0*1.6/(1.5+exp(-(Y(8)+60.0)/40.0));

if (Version == 0.0)
   i_Ca_p_max = i_Ca_p_max_Centre_Published+FCell*(i_Ca_p_max_Periphery_Published-i_Ca_p_max_Centre_Published);
elseif (Version == 1.0)
   i_Ca_p_max = i_Ca_p_max_Centre_0DCapable+FCell*(i_Ca_p_max_Periphery_0DCapable-i_Ca_p_max_Centre_0DCapable);
else
   i_Ca_p_max = i_Ca_p_max_Centre_1DCapable+FCell*(i_Ca_p_max_Periphery_1DCapable-i_Ca_p_max_Centre_1DCapable);
end;

i_Ca_p = i_Ca_p_max*Ca_i/(Ca_i+0.0004);
dY(8, 1) = -1.0/Cm*(i_Na+i_Ca_L+i_Ca_T+i_to+i_sus+i_K_r+i_K_s+i_f_Na+i_f_K+i_b_Na+i_b_Ca+i_b_K+i_NaCa+i_p+i_Ca_p);

if (Version ~= 2.0)
   P_af_infinity = 1.0/(1.0+exp(-(Y(8)+14.2)/10.6));
else
   P_af_infinity = 1.0/(1.0+exp(-(Y(8)+13.2)/10.6));
end;

if (Version ~= 2.0)
   tau_P_af = 1.0/(37.2*exp((Y(8)-9.0)/15.9)+0.96*exp(-(Y(8)-9.0)/22.5));
else
   tau_P_af = 1.0/(37.2*exp((Y(8)-10.0)/15.9)+0.96*exp(-(Y(8)-10.0)/22.5));
end;

dY(9, 1) = (P_af_infinity-Y(9))/tau_P_af;
P_as_infinity = P_af_infinity;

if (Version ~= 2.0)
   tau_P_as = 1.0/(4.2*exp((Y(8)-9.0)/17.0)+0.15*exp(-(Y(8)-9.0)/21.6));
else
   tau_P_as = 1.0/(4.2*exp((Y(8)-10.0)/17.0)+0.15*exp(-(Y(8)-10.0)/21.6));
end;

dY(10, 1) = (P_as_infinity-Y(10))/tau_P_as;

if (Version == 0.0)
   tau_P_i = 0.002;
elseif (Version == 1.0)
   tau_P_i = 0.002;
else
   tau_P_i = 0.006;
end;

P_i_infinity = 1.0/(1.0+exp((Y(8)+18.6)/10.1));
dY(11, 1) = (P_i_infinity-Y(11))/tau_P_i;
alpha_xs = 14.0/(1.0+exp(-(Y(8)-40.0)/9.0));
beta_xs = 1.0*exp(-Y(8)/45.0);
dY(12, 1) = alpha_xs*(1.0-Y(12))-beta_xs*Y(12);
h1_infinity = 1.0/(1.0+exp((Y(8)+66.1)/6.4));
tau_h1 = 3.717e-6*exp(-0.2815*(Y(8)+17.11))/(1.0+0.003732*exp(-0.3426*(Y(8)+37.76)))+0.0005977;
dY(13, 1) = (h1_infinity-Y(13))/tau_h1;
h2_infinity = h1_infinity;
tau_h2 = 3.186e-8*exp(-0.6219*(Y(8)+18.8))/(1.0+7.189e-5*exp(-0.6683*(Y(8)+34.07)))+0.003556;
dY(14, 1) = (h2_infinity-Y(14))/tau_h2;

if (Version == 0.0)
   m_infinity = (1.0/(1.0+exp(-Y(8)/5.46)))^(1.0/3.0);
else
   m_infinity = (1.0/(1.0+exp(-(Y(8)+30.32)/5.46)))^(1.0/3.0);
end;

if (Version == 0.0)
   tau_m = 0.0006247/(0.832*exp(-0.335*(Y(8)+56.7))+0.627*exp(0.082*(Y(8)+65.01)))+4.0e-5;
else
   tau_m = 0.0006247/(0.8322166*exp(-0.33566*(Y(8)+56.7062))+0.6274*exp(0.0823*(Y(8)+65.0131)))+4.569e-5;
end;

dY(15, 1) = (m_infinity-Y(15))/tau_m;

%===============================================================================
% End of file
%===============================================================================
