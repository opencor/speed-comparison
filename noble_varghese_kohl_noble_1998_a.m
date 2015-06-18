%===============================================================================
% CellML file:   D:\Desktop\Models\noble_varghese_kohl_noble_1998_a.cellml
% CellML model:  noble_varghese_kohl_noble_1998_basic
% Date and time: 17/06/2015 at 23:07:56
%-------------------------------------------------------------------------------
% Conversion from CellML 1.0 to MATLAB (init) was done using COR (0.9.31.1409)
%    Copyright 2002-2015 Dr Alan Garny
%    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
%-------------------------------------------------------------------------------
% http://www.cellml.org/
%===============================================================================

function dY = noble_varghese_kohl_noble_1998_a(time, Y)

%-------------------------------------------------------------------------------
% Initial conditions
%-------------------------------------------------------------------------------

% Y = [0.0, 0.9349197, 0.9651958, 1.0, 0.0042614, 0.4068154, 0.9944036, 0.0016203, 0.0005555, 0.0003542, 1.88e-5, 1.4e-5, 0.4481927, 0.4531889, 136.5644281, 7.3321223, -92.849333, 1.03e-5, 2.0e-7, 0.001302, 0.0, 0.9948645];

% YNames = {'d', 'f2', 'f2ds', 'f', 'ActFrac', 'ProdFrac', 'h', 'm', 'Ca_Calmod', 'Ca_Trop', 'Ca_ds', 'Ca_i', 'Ca_rel', 'Ca_up', 'K_i', 'Na_i', 'V', 'xr1', 'xr2', 'xs', 'r', 's'};
% YUnits = {'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'millimolar', 'millimolar', 'millimolar', 'millimolar', 'millimolar', 'millimolar', 'millimolar', 'millimolar', 'millivolt', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless'};
% YComponents = {'L_type_Ca_channel_d_gate', 'L_type_Ca_channel_f2_gate', 'L_type_Ca_channel_f2ds_gate', 'L_type_Ca_channel_f_gate', 'calcium_release', 'calcium_release', 'fast_sodium_current_h_gate', 'fast_sodium_current_m_gate', 'intracellular_calcium_concentration', 'intracellular_calcium_concentration', 'intracellular_calcium_concentration', 'intracellular_calcium_concentration', 'intracellular_calcium_concentration', 'intracellular_calcium_concentration', 'intracellular_potassium_concentration', 'intracellular_sodium_concentration', 'membrane', 'rapid_delayed_rectifier_potassium_current_xr1_gate', 'rapid_delayed_rectifier_potassium_current_xr2_gate', 'slow_delayed_rectifier_potassium_current_xs_gate', 'transient_outward_current_r_gate', 'transient_outward_current_s_gate'};

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: d (dimensionless) (in L_type_Ca_channel_d_gate)
% 2: f2 (dimensionless) (in L_type_Ca_channel_f2_gate)
% 3: f2ds (dimensionless) (in L_type_Ca_channel_f2ds_gate)
% 4: f (dimensionless) (in L_type_Ca_channel_f_gate)
% 5: ActFrac (dimensionless) (in calcium_release)
% 6: ProdFrac (dimensionless) (in calcium_release)
% 7: h (dimensionless) (in fast_sodium_current_h_gate)
% 8: m (dimensionless) (in fast_sodium_current_m_gate)
% 9: Ca_Calmod (millimolar) (in intracellular_calcium_concentration)
% 10: Ca_Trop (millimolar) (in intracellular_calcium_concentration)
% 11: Ca_ds (millimolar) (in intracellular_calcium_concentration)
% 12: Ca_i (millimolar) (in intracellular_calcium_concentration)
% 13: Ca_rel (millimolar) (in intracellular_calcium_concentration)
% 14: Ca_up (millimolar) (in intracellular_calcium_concentration)
% 15: K_i (millimolar) (in intracellular_potassium_concentration)
% 16: Na_i (millimolar) (in intracellular_sodium_concentration)
% 17: V (millivolt) (in membrane)
% 18: xr1 (dimensionless) (in rapid_delayed_rectifier_potassium_current_xr1_gate)
% 19: xr2 (dimensionless) (in rapid_delayed_rectifier_potassium_current_xr2_gate)
% 20: xs (dimensionless) (in slow_delayed_rectifier_potassium_current_xs_gate)
% 21: r (dimensionless) (in transient_outward_current_r_gate)
% 22: s (dimensionless) (in transient_outward_current_s_gate)

%-------------------------------------------------------------------------------
% Constants
%-------------------------------------------------------------------------------

speed_d = 3.0;   % dimensionless (in L_type_Ca_channel_d_gate)
delta_f = 0.0001;   % millivolt (in L_type_Ca_channel_f_gate)
speed_f = 0.3;   % dimensionless (in L_type_Ca_channel_f_gate)
FrICa = 1.0;   % dimensionless (in L_type_Ca_channel)
Km_f2 = 100000.0;   % millimolar (in L_type_Ca_channel)
Km_f2ds = 0.001;   % millimolar (in L_type_Ca_channel)
P_CaK = 0.002;   % dimensionless (in L_type_Ca_channel)
P_CaNa = 0.01;   % dimensionless (in L_type_Ca_channel)
P_Ca_L = 0.1;   % nanoA_per_millimolar (in L_type_Ca_channel)
R_decay = 20.0;   % per_second (in L_type_Ca_channel)
g_bca = 0.00025;   % microS (in calcium_background_current)
K_leak_rate = 0.05;   % per_second (in calcium_release)
K_m_Ca_cyt = 0.0005;   % millimolar (in calcium_release)
K_m_Ca_ds = 0.01;   % millimolar (in calcium_release)
K_m_rel = 250.0;   % per_second (in calcium_release)
Ca_o = 2.0;   % millimolar (in extracellular_calcium_concentration)
K_o = 4.0;   % millimolar (in extracellular_potassium_concentration)
Na_o = 140.0;   % millimolar (in extracellular_sodium_concentration)
shift_h = 0.0;   % millivolt (in fast_sodium_current_h_gate)
delta_m = 1.0e-5;   % millivolt (in fast_sodium_current_m_gate)
g_Na = 2.5;   % microS (in fast_sodium_current)
Calmod = 0.02;   % millimolar (in intracellular_calcium_concentration)
Kdecay = 10.0;   % per_second (in intracellular_calcium_concentration)
Trop = 0.05;   % millimolar (in intracellular_calcium_concentration)
V_ds_ratio = 0.1;   % dimensionless (in intracellular_calcium_concentration)
V_e_ratio = 0.4;   % dimensionless (in intracellular_calcium_concentration)
V_rel_ratio = 0.1;   % dimensionless (in intracellular_calcium_concentration)
V_up_ratio = 0.01;   % dimensionless (in intracellular_calcium_concentration)
alpha_Calmod = 100000.0;   % per_millimolar_second (in intracellular_calcium_concentration)
alpha_Trop = 100000.0;   % per_millimolar_second (in intracellular_calcium_concentration)
beta_Calmod = 50.0;   % per_second (in intracellular_calcium_concentration)
beta_Trop = 200.0;   % per_second (in intracellular_calcium_concentration)
length = 0.074;   % micrometre (in intracellular_calcium_concentration)
radius = 0.012;   % micrometre (in intracellular_calcium_concentration)
Cm = 9.5e-5;   % microF (in membrane)
F = 96485.3415;   % coulomb_per_mole (in membrane)
R = 8314.472;   % joule_per_kilomole_kelvin (in membrane)
T = 310.0;   % kelvin (in membrane)
stim_amplitude = -3.0;   % nanoA (in membrane)
stim_duration = 0.003;   % second (in membrane)
stim_end = 100000.0;   % second (in membrane)
stim_period = 1.0;   % second (in membrane)
stim_start = 0.1;   % second (in membrane)
g_pna = 0.004;   % microS (in persistent_sodium_current)
g_Kr1 = 0.0021;   % microS (in rapid_delayed_rectifier_potassium_current)
g_Kr2 = 0.0013;   % microS (in rapid_delayed_rectifier_potassium_current)
P_kna = 0.03;   % dimensionless (in reversal_potentials)
K_cyca = 0.0003;   % millimolar (in sarcoplasmic_reticulum_calcium_pump)
K_srca = 0.5;   % millimolar (in sarcoplasmic_reticulum_calcium_pump)
K_xcs = 0.4;   % dimensionless (in sarcoplasmic_reticulum_calcium_pump)
alpha_up = 0.4;   % millimolar_per_second (in sarcoplasmic_reticulum_calcium_pump)
beta_up = 0.03;   % millimolar_per_second (in sarcoplasmic_reticulum_calcium_pump)
g_Ks = 0.0026;   % microS (in slow_delayed_rectifier_potassium_current)
g_bna = 0.0006;   % microS (in sodium_background_current)
FRiNaCa = 0.001;   % dimensionless (in sodium_calcium_exchanger)
d_NaCa = 0.0;   % dimensionless (in sodium_calcium_exchanger)
gamma = 0.5;   % dimensionless (in sodium_calcium_exchanger)
k_NaCa = 0.0005;   % nanoA (in sodium_calcium_exchanger)
n_NaCa = 3.0;   % dimensionless (in sodium_calcium_exchanger)
K_mK = 1.0;   % millimolar (in sodium_potassium_pump)
K_mNa = 40.0;   % millimolar (in sodium_potassium_pump)
i_NaK_max = 0.7;   % nanoA (in sodium_potassium_pump)
K_mk1 = 10.0;   % millimolar (in time_independent_potassium_current)
g_K1 = 0.5;   % microS (in time_independent_potassium_current)
g_to = 0.005;   % microS (in transient_outward_current)
g_tos = 0.0;   % dimensionless (in transient_outward_current)

%-------------------------------------------------------------------------------
% Computed variables
%-------------------------------------------------------------------------------

% E0_d (millivolt) (in L_type_Ca_channel_d_gate)
% alpha_d (per_second) (in L_type_Ca_channel_d_gate)
% beta_d (per_second) (in L_type_Ca_channel_d_gate)
% E0_f (millivolt) (in L_type_Ca_channel_f_gate)
% alpha_f (per_second) (in L_type_Ca_channel_f_gate)
% beta_f (per_second) (in L_type_Ca_channel_f_gate)
% i_Ca_L (nanoA) (in L_type_Ca_channel)
% i_Ca_L_Ca_cyt (nanoA) (in L_type_Ca_channel)
% i_Ca_L_Ca_ds (nanoA) (in L_type_Ca_channel)
% i_Ca_L_K_cyt (nanoA) (in L_type_Ca_channel)
% i_Ca_L_K_ds (nanoA) (in L_type_Ca_channel)
% i_Ca_L_Na_cyt (nanoA) (in L_type_Ca_channel)
% i_Ca_L_Na_ds (nanoA) (in L_type_Ca_channel)
% i_b_Ca (nanoA) (in calcium_background_current)
% ActRate (per_second) (in calcium_release)
% CadsReg (dimensionless) (in calcium_release)
% CaiReg (dimensionless) (in calcium_release)
% InactRate (per_second) (in calcium_release)
% PrecFrac (dimensionless) (in calcium_release)
% RegBindSite (dimensionless) (in calcium_release)
% SpeedRel (dimensionless) (in calcium_release)
% VoltDep (dimensionless) (in calcium_release)
% i_rel (millimolar_per_second) (in calcium_release)
% i_trans (millimolar_per_second) (in calcium_translocation)
% alpha_h (per_second) (in fast_sodium_current_h_gate)
% beta_h (per_second) (in fast_sodium_current_h_gate)
% E0_m (millivolt) (in fast_sodium_current_m_gate)
% alpha_m (per_second) (in fast_sodium_current_m_gate)
% beta_m (per_second) (in fast_sodium_current_m_gate)
% i_Na (nanoA) (in fast_sodium_current)
% V_Cell (micrometre3) (in intracellular_calcium_concentration)
% V_i (micrometre3) (in intracellular_calcium_concentration)
% V_i_ratio (dimensionless) (in intracellular_calcium_concentration)
% i_Stim (nanoA) (in membrane)
% i_p_Na (nanoA) (in persistent_sodium_current)
% alpha_xr1 (per_second) (in rapid_delayed_rectifier_potassium_current_xr1_gate)
% beta_xr1 (per_second) (in rapid_delayed_rectifier_potassium_current_xr1_gate)
% alpha_xr2 (per_second) (in rapid_delayed_rectifier_potassium_current_xr2_gate)
% beta_xr2 (per_second) (in rapid_delayed_rectifier_potassium_current_xr2_gate)
% i_Kr (nanoA) (in rapid_delayed_rectifier_potassium_current)
% E_Ca (millivolt) (in reversal_potentials)
% E_K (millivolt) (in reversal_potentials)
% E_Ks (millivolt) (in reversal_potentials)
% E_Na (millivolt) (in reversal_potentials)
% E_mh (millivolt) (in reversal_potentials)
% K_1 (dimensionless) (in sarcoplasmic_reticulum_calcium_pump)
% K_2 (millimolar) (in sarcoplasmic_reticulum_calcium_pump)
% i_up (millimolar_per_second) (in sarcoplasmic_reticulum_calcium_pump)
% alpha_xs (per_second) (in slow_delayed_rectifier_potassium_current_xs_gate)
% beta_xs (per_second) (in slow_delayed_rectifier_potassium_current_xs_gate)
% i_Ks (nanoA) (in slow_delayed_rectifier_potassium_current)
% i_b_Na (nanoA) (in sodium_background_current)
% i_NaCa (nanoA) (in sodium_calcium_exchanger)
% i_NaCa_cyt (nanoA) (in sodium_calcium_exchanger)
% i_NaCa_ds (nanoA) (in sodium_calcium_exchanger)
% i_NaK (nanoA) (in sodium_potassium_pump)
% i_K1 (nanoA) (in time_independent_potassium_current)
% alpha_s (per_second) (in transient_outward_current_s_gate)
% beta_s (per_second) (in transient_outward_current_s_gate)
% i_to (nanoA) (in transient_outward_current)

%-------------------------------------------------------------------------------
% Computation
%-------------------------------------------------------------------------------

% time (second)

i_Ca_L_Ca_cyt = (1.0-FrICa)*4.0*P_Ca_L*Y(1)*Y(4)*Y(2)*(Y(17)-50.0)*F/(R*T)/(1.0-exp(-(Y(17)-50.0)*F*2.0/(R*T)))*(Y(12)*exp(100.0*F/(R*T))-Ca_o*exp(-(Y(17)-50.0)*F*2.0/(R*T)));
i_Ca_L_K_cyt = (1.0-FrICa)*P_CaK*P_Ca_L*Y(1)*Y(4)*Y(2)*(Y(17)-50.0)*F/(R*T)/(1.0-exp(-(Y(17)-50.0)*F/(R*T)))*(Y(15)*exp(50.0*F/(R*T))-K_o*exp(-(Y(17)-50.0)*F/(R*T)));
i_Ca_L_Na_cyt = (1.0-FrICa)*P_CaNa*P_Ca_L*Y(1)*Y(4)*Y(2)*(Y(17)-50.0)*F/(R*T)/(1.0-exp(-(Y(17)-50.0)*F/(R*T)))*(Y(16)*exp(50.0*F/(R*T))-Na_o*exp(-(Y(17)-50.0)*F/(R*T)));
i_Ca_L_Ca_ds = FrICa*4.0*P_Ca_L*Y(1)*Y(4)*Y(3)*(Y(17)-50.0)*F/(R*T)/(1.0-exp(-(Y(17)-50.0)*F*2.0/(R*T)))*(Y(12)*exp(100.0*F/(R*T))-Ca_o*exp(-(Y(17)-50.0)*F*2.0/(R*T)));
i_Ca_L_K_ds = FrICa*P_CaK*P_Ca_L*Y(1)*Y(4)*Y(3)*(Y(17)-50.0)*F/(R*T)/(1.0-exp(-(Y(17)-50.0)*F/(R*T)))*(Y(15)*exp(50.0*F/(R*T))-K_o*exp(-(Y(17)-50.0)*F/(R*T)));
i_Ca_L_Na_ds = FrICa*P_CaNa*P_Ca_L*Y(1)*Y(4)*Y(3)*(Y(17)-50.0)*F/(R*T)/(1.0-exp(-(Y(17)-50.0)*F/(R*T)))*(Y(16)*exp(50.0*F/(R*T))-Na_o*exp(-(Y(17)-50.0)*F/(R*T)));
i_Ca_L = i_Ca_L_Ca_cyt+i_Ca_L_K_cyt+i_Ca_L_Na_cyt+i_Ca_L_Ca_ds+i_Ca_L_K_ds+i_Ca_L_Na_ds;
E0_d = Y(17)+24.0-5.0;

if (abs(E0_d) < 0.0001)
   alpha_d = 120.0;
else
   alpha_d = 30.0*E0_d/(1.0-exp(-E0_d/4.0));
end;

if (abs(E0_d) < 0.0001)
   beta_d = 120.0;
else
   beta_d = 12.0*E0_d/(exp(E0_d/10.0)-1.0);
end;

dY(1, 1) = speed_d*(alpha_d*(1.0-Y(1))-beta_d*Y(1));
dY(2, 1) = 1.0-1.0*(Y(12)/(Km_f2+Y(12))+Y(2));
dY(3, 1) = R_decay*(1.0-(Y(11)/(Km_f2ds+Y(11))+Y(3)));
E0_f = Y(17)+34.0;

if (abs(E0_f) < delta_f)
   alpha_f = 25.0;
else
   alpha_f = 6.25*E0_f/(exp(E0_f/4.0)-1.0);
end;

beta_f = 12.0/(1.0+exp(-1.0*(Y(17)+34.0)/4.0));
dY(4, 1) = speed_f*(alpha_f*(1.0-Y(4))-beta_f*Y(4));
E_Ca = 0.5*R*T/F*log(Ca_o/Y(12));
i_b_Ca = g_bca*(Y(17)-E_Ca);
VoltDep = exp(0.08*(Y(17)-40.0));
CaiReg = Y(12)/(Y(12)+K_m_Ca_cyt);
CadsReg = Y(11)/(Y(11)+K_m_Ca_ds);
RegBindSite = CaiReg+(1.0-CaiReg)*CadsReg;
ActRate = 0.0*VoltDep+500.0*RegBindSite^2.0;
InactRate = 60.0+500.0*RegBindSite^2.0;

if (Y(17) < -50.0)
   SpeedRel = 5.0;
else
   SpeedRel = 1.0;
end;

PrecFrac = 1.0-Y(5)-Y(6);
dY(5, 1) = PrecFrac*SpeedRel*ActRate-Y(5)*SpeedRel*InactRate;
dY(6, 1) = Y(5)*SpeedRel*InactRate-SpeedRel*1.0*Y(6);
i_rel = ((Y(5)/(Y(5)+0.25))^2.0*K_m_rel+K_leak_rate)*Y(13);
i_trans = 50.0*(Y(14)-Y(13));
E_mh = R*T/F*log((Na_o+0.12*K_o)/(Y(16)+0.12*Y(15)));
i_Na = g_Na*Y(8)^3.0*Y(7)*(Y(17)-E_mh);
alpha_h = 20.0*exp(-0.125*(Y(17)+75.0-shift_h));
beta_h = 2000.0/(1.0+320.0*exp(-0.1*(Y(17)+75.0-shift_h)));
dY(7, 1) = alpha_h*(1.0-Y(7))-beta_h*Y(7);
E0_m = Y(17)+41.0;

if (abs(E0_m) < delta_m)
   alpha_m = 2000.0;
else
   alpha_m = 200.0*E0_m/(1.0-exp(-0.1*E0_m));
end;

beta_m = 8000.0*exp(-0.056*(Y(17)+66.0));
dY(8, 1) = alpha_m*(1.0-Y(8))-beta_m*Y(8);
V_Cell = 3.141592654*radius^2.0*length;
V_i_ratio = 1.0-V_e_ratio-V_up_ratio-V_rel_ratio;
V_i = V_Cell*V_i_ratio;
i_NaCa_cyt = (1.0-FRiNaCa)*k_NaCa*(exp(gamma*(n_NaCa-2.0)*Y(17)*F/(R*T))*Y(16)^n_NaCa*Ca_o-exp((gamma-1.0)*(n_NaCa-2.0)*Y(17)*F/(R*T))*Na_o^n_NaCa*Y(12))/((1.0+d_NaCa*(Y(12)*Na_o^n_NaCa+Ca_o*Y(16)^n_NaCa))*(1.0+Y(12)/0.0069));
i_NaCa_ds = FRiNaCa*k_NaCa*(exp(gamma*(n_NaCa-2.0)*Y(17)*F/(R*T))*Y(16)^n_NaCa*Ca_o-exp((gamma-1.0)*(n_NaCa-2.0)*Y(17)*F/(R*T))*Na_o^n_NaCa*Y(11))/((1.0+d_NaCa*(Y(11)*Na_o^n_NaCa+Ca_o*Y(16)^n_NaCa))*(1.0+Y(11)/0.0069));
dY(9, 1) = alpha_Calmod*Y(12)*(Calmod-Y(9))-beta_Calmod*Y(9);
dY(10, 1) = alpha_Trop*Y(12)*(Trop-Y(10))-beta_Trop*Y(10);
K_1 = K_cyca*K_xcs/K_srca;
K_2 = Y(12)+Y(14)*K_1+K_cyca*K_xcs+K_cyca;
i_up = Y(12)/K_2*alpha_up-Y(14)*K_1/K_2*beta_up;
dY(12, 1) = -1.0/(2.0*1.0*V_i*F)*(i_Ca_L_Ca_cyt+i_b_Ca-2.0*i_NaCa_cyt-2.0*i_NaCa_ds)+Y(11)*V_ds_ratio*Kdecay+i_rel*V_rel_ratio/V_i_ratio-dY(9, 1)-dY(10, 1)-i_up;
dY(11, 1) = -1.0*i_Ca_L_Ca_ds/(2.0*1.0*V_ds_ratio*V_i*F)-Y(11)*Kdecay;
dY(14, 1) = V_i_ratio/V_up_ratio*i_up-i_trans;
dY(13, 1) = V_up_ratio/V_rel_ratio*i_trans-i_rel;
E_K = R*T/F*log(K_o/Y(15));
i_K1 = g_K1*K_o/(K_o+K_mk1)*(Y(17)-E_K)/(1.0+exp((Y(17)-E_K-10.0)*F*1.25/(R*T)));
i_Kr = (g_Kr1*Y(18)+g_Kr2*Y(19))*1.0/(1.0+exp((Y(17)+9.0)/22.4))*(Y(17)-E_K);
E_Ks = R*T/F*log((K_o+P_kna*Na_o)/(Y(15)+P_kna*Y(16)));
i_Ks = g_Ks*Y(20)^2.0*(Y(17)-E_Ks);
i_to = g_to*(g_tos+Y(22)*(1.0-g_tos))*Y(21)*(Y(17)-E_K);
i_NaK = i_NaK_max*K_o/(K_mK+K_o)*Y(16)/(K_mNa+Y(16));
dY(15, 1) = -1.0/(1.0*V_i*F)*(i_K1+i_Kr+i_Ks+i_Ca_L_K_cyt+i_Ca_L_K_ds+i_to-2.0*i_NaK);
E_Na = R*T/F*log(Na_o/Y(16));
i_p_Na = g_pna*1.0/(1.0+exp(-(Y(17)+52.0)/8.0))*(Y(17)-E_Na);
i_b_Na = g_bna*(Y(17)-E_Na);
dY(16, 1) = -1.0/(1.0*V_i*F)*(i_Na+i_p_Na+i_b_Na+3.0*i_NaK+3.0*i_NaCa_cyt+i_Ca_L_Na_cyt+i_Ca_L_Na_ds);

if ((time >= stim_start) && (time <= stim_end) && (time-stim_start-floor((time-stim_start)/stim_period)*stim_period <= stim_duration))
   i_Stim = stim_amplitude;
else
   i_Stim = 0.0;
end;

dY(17, 1) = -1.0/Cm*(i_Stim+i_K1+i_to+i_Kr+i_Ks+i_NaK+i_Na+i_b_Na+i_p_Na+i_Ca_L_Na_cyt+i_Ca_L_Na_ds+i_NaCa_cyt+i_NaCa_ds+i_Ca_L_Ca_cyt+i_Ca_L_Ca_ds+i_Ca_L_K_cyt+i_Ca_L_K_ds+i_b_Ca);
alpha_xr1 = 50.0/(1.0+exp(-(Y(17)-5.0)/9.0));
beta_xr1 = 0.05*exp(-(Y(17)-20.0)/15.0);
dY(18, 1) = alpha_xr1*(1.0-Y(18))-beta_xr1*Y(18);
alpha_xr2 = 50.0/(1.0+exp(-(Y(17)-5.0)/9.0));
beta_xr2 = 0.4*exp(-((Y(17)+30.0)/30.0)^3.0);
dY(19, 1) = alpha_xr2*(1.0-Y(19))-beta_xr2*Y(19);
alpha_xs = 14.0/(1.0+exp(-(Y(17)-40.0)/9.0));
beta_xs = 1.0*exp(-Y(17)/45.0);
dY(20, 1) = alpha_xs*(1.0-Y(20))-beta_xs*Y(20);
i_NaCa = i_NaCa_cyt+i_NaCa_ds;
dY(21, 1) = 333.0*(1.0/(1.0+exp(-(Y(17)+4.0)/5.0))-Y(21));
alpha_s = 0.033*exp(-Y(17)/17.0);
beta_s = 33.0/(1.0+exp(-0.125*(Y(17)+10.0)));
dY(22, 1) = alpha_s*(1.0-Y(22))-beta_s*Y(22);

%===============================================================================
% End of file
%===============================================================================
