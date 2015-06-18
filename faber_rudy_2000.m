%===============================================================================
% CellML file:   D:\Desktop\Models\faber_rudy_2000.cellml
% CellML model:  faber_2000
% Date and time: 17/06/2015 at 23:12:53
%-------------------------------------------------------------------------------
% Conversion from CellML 1.0 to MATLAB (init) was done using COR (0.9.31.1409)
%    Copyright 2002-2015 Dr Alan Garny
%    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
%-------------------------------------------------------------------------------
% http://www.cellml.org/
%===============================================================================

function dY = faber_rudy_2000(time, Y)

%-------------------------------------------------------------------------------
% Initial conditions
%-------------------------------------------------------------------------------

% Y = [5.82597094505446e-6, 0.997765362821995, 0.00136737866785149, 0.98881442877378, 9.65910542308504e-196, 5.33944967562997e-195, 0.000129515197402902, 1.12791401197882, 1.76731003671612, 0.000117482029668194, 2.69380318286645e-196, 0.0, 0.0, 0.0, -85.2119207874627, 0.985596581239651, 0.990898461370389, 0.00149183115674257, 141.056872392446, 13.3649235394859, 0.000204700363126417, 0.00660746743356887, 0.0303768241233812, 0.999945568566232, 0.0144622472219576];

% YNames = {'d', 'f', 'b', 'g', 'APtrack', 'APtrack2', 'APtrack3', 'Ca_JSR', 'Ca_NSR', 'Cai', 'Cainfluxtrack', 'OVRLDtrack', 'OVRLDtrack2', 'OVRLDtrack3', 'V', 'h', 'j', 'm', 'Ki', 'Nai', 'xr', 'xs1', 'xs2', 'ydv', 'zdv'};
% YUnits = {'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'millimolar', 'millimolar', 'millimolar', 'millimolar', 'dimensionless', 'dimensionless', 'dimensionless', 'millivolt', 'dimensionless', 'dimensionless', 'dimensionless', 'millimolar', 'millimolar', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless'};
% YComponents = {'L_type_Ca_channel_d_gate', 'L_type_Ca_channel_f_gate', 'T_type_Ca_channel_b_gate', 'T_type_Ca_channel_g_gate', 'calcium_dynamics', 'calcium_dynamics', 'calcium_dynamics', 'calcium_dynamics', 'calcium_dynamics', 'calcium_dynamics', 'calcium_dynamics', 'calcium_dynamics', 'calcium_dynamics', 'calcium_dynamics', 'cell', 'fast_sodium_current_h_gate', 'fast_sodium_current_j_gate', 'fast_sodium_current_m_gate', 'ionic_concentrations', 'ionic_concentrations', 'rapid_delayed_rectifier_potassium_current_xr_gate', 'slow_delayed_rectifier_potassium_current_xs1_gate', 'slow_delayed_rectifier_potassium_current_xs2_gate', 'transient_outward_current_ydv_gate', 'transient_outward_current_zdv_gate'};

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: d (dimensionless) (in L_type_Ca_channel_d_gate)
% 2: f (dimensionless) (in L_type_Ca_channel_f_gate)
% 3: b (dimensionless) (in T_type_Ca_channel_b_gate)
% 4: g (dimensionless) (in T_type_Ca_channel_g_gate)
% 5: APtrack (dimensionless) (in calcium_dynamics)
% 6: APtrack2 (dimensionless) (in calcium_dynamics)
% 7: APtrack3 (dimensionless) (in calcium_dynamics)
% 8: Ca_JSR (millimolar) (in calcium_dynamics)
% 9: Ca_NSR (millimolar) (in calcium_dynamics)
% 10: Cai (millimolar) (in calcium_dynamics)
% 11: Cainfluxtrack (millimolar) (in calcium_dynamics)
% 12: OVRLDtrack (dimensionless) (in calcium_dynamics)
% 13: OVRLDtrack2 (dimensionless) (in calcium_dynamics)
% 14: OVRLDtrack3 (dimensionless) (in calcium_dynamics)
% 15: V (millivolt) (in cell)
% 16: h (dimensionless) (in fast_sodium_current_h_gate)
% 17: j (dimensionless) (in fast_sodium_current_j_gate)
% 18: m (dimensionless) (in fast_sodium_current_m_gate)
% 19: Ki (millimolar) (in ionic_concentrations)
% 20: Nai (millimolar) (in ionic_concentrations)
% 21: xr (dimensionless) (in rapid_delayed_rectifier_potassium_current_xr_gate)
% 22: xs1 (dimensionless) (in slow_delayed_rectifier_potassium_current_xs1_gate)
% 23: xs2 (dimensionless) (in slow_delayed_rectifier_potassium_current_xs2_gate)
% 24: ydv (dimensionless) (in transient_outward_current_ydv_gate)
% 25: zdv (dimensionless) (in transient_outward_current_zdv_gate)

%-------------------------------------------------------------------------------
% Constants
%-------------------------------------------------------------------------------

ATPi = 3.0;   % millimolar (in ATP_sensitive_potassium_current)
hATP = 2.0;   % dimensionless (in ATP_sensitive_potassium_current)
i_K_ATP_on = 1.0;   % dimensionless (in ATP_sensitive_potassium_current)
kATP = 0.00025;   % millimolar (in ATP_sensitive_potassium_current)
nATP = 0.24;   % dimensionless (in ATP_sensitive_potassium_current)
nicholsarea = 5.0e-5;   % dimensionless (in ATP_sensitive_potassium_current)
Km_Ca = 0.0006;   % millimolar (in L_type_Ca_channel_f_Ca_gate)
P_Ca = 0.00054;   % litre_per_farad_millisecond (in L_type_Ca_channel)
P_K = 1.93e-7;   % litre_per_farad_millisecond (in L_type_Ca_channel)
P_Na = 6.75e-7;   % litre_per_farad_millisecond (in L_type_Ca_channel)
gamma_Cai = 1.0;   % dimensionless (in L_type_Ca_channel)
gamma_Cao = 0.341;   % dimensionless (in L_type_Ca_channel)
gamma_Ki = 0.75;   % dimensionless (in L_type_Ca_channel)
gamma_Ko = 0.75;   % dimensionless (in L_type_Ca_channel)
gamma_Nai = 0.75;   % dimensionless (in L_type_Ca_channel)
gamma_Nao = 0.75;   % dimensionless (in L_type_Ca_channel)
c1 = 0.00025;   % microA_per_microF (in Na_Ca_exchanger)
c2 = 0.0001;   % dimensionless (in Na_Ca_exchanger)
gamma = 0.15;   % dimensionless (in Na_Ca_exchanger)
g_CaT = 0.05;   % milliS_per_microF (in T_type_Ca_channel)
g_Cab = 0.003016;   % milliS_per_microF (in calcium_background_current)
CMDN_max = 0.05;   % millimolar (in calcium_dynamics)
CSQN_max = 10.0;   % millimolar (in calcium_dynamics)
CSQNthresh = 0.7;   % dimensionless (in calcium_dynamics)
Ca_NSR_max = 15.0;   % millimolar (in calcium_dynamics)
Cao = 1.8;   % millimolar (in calcium_dynamics)
G_rel_max = 60.0;   % per_ms (in calcium_dynamics)
G_rel_overload = 4.0;   % per_ms (in calcium_dynamics)
I_up = 0.00875;   % millimolar_per_ms (in calcium_dynamics)
K_mCMDN = 0.00238;   % millimolar (in calcium_dynamics)
K_mCSQN = 0.8;   % millimolar (in calcium_dynamics)
K_mTn = 0.0005;   % millimolar (in calcium_dynamics)
K_mrel = 0.0008;   % millimolar (in calcium_dynamics)
K_mup = 0.00092;   % millimolar (in calcium_dynamics)
Logicthresh = 0.98;   % dimensionless (in calcium_dynamics)
Tn_max = 0.07;   % millimolar (in calcium_dynamics)
delta_Ca_ith = 0.00018;   % millimolar (in calcium_dynamics)
tau_tr = 180.0;   % ms (in calcium_dynamics)
F = 96485.0;   % coulomb_per_mole (in cell)
R = 8314.0;   % joule_per_kilomole_kelvin (in cell)
T = 310.0;   % kelvin (in cell)
stim_amplitude = -25.5;   % microA_per_microF (in cell)
stim_duration = 2.0;   % ms (in cell)
stim_end = 9000000.0;   % ms (in cell)
stim_period = 1000.0;   % ms (in cell)
stim_start = 100.0;   % ms (in cell)
delta_m = 1.0e-5;   % millivolt (in fast_sodium_current_m_gate)
g_Na = 16.0;   % milliS_per_microF (in fast_sodium_current)
A_cap = 0.0001534;   % cm2 (in geometry)
preplength = 0.1;   % mm (in geometry)
radius = 0.011;   % mm (in geometry)
Ko = 5.4;   % millimolar (in ionic_concentrations)
Nao = 140.0;   % millimolar (in ionic_concentrations)
K_m_ns_Ca = 0.0012;   % millimolar (in non_specific_calcium_activated_current)
g_Kp = 0.00552;   % milliS_per_microF (in plateau_potassium_current)
G_Kr = 0.02614;   % milliS_per_microF (in rapid_delayed_rectifier_potassium_current)
I_pCa = 1.15;   % microA_per_microF (in sarcolemmal_calcium_pump)
K_mpCa = 0.0005;   % millimolar (in sarcolemmal_calcium_pump)
G_Ks = 0.433;   % milliS_per_microF (in slow_delayed_rectifier_potassium_current)
PNaK = 0.01833;   % dimensionless (in slow_delayed_rectifier_potassium_current)
g_K_Na = 0.12848;   % milliS_per_microF (in sodium_activated_potassium_current)
kdKNa = 66.0;   % millimolar (in sodium_activated_potassium_current)
nKNa = 2.8;   % dimensionless (in sodium_activated_potassium_current)
g_Nab = 0.004;   % milliS_per_microF (in sodium_background_current)
I_NaK = 2.25;   % microA_per_microF (in sodium_potassium_pump)
K_mKo = 1.5;   % millimolar (in sodium_potassium_pump)
K_mNai = 10.0;   % millimolar (in sodium_potassium_pump)
G_K1 = 0.75;   % milliS_per_microF (in time_independent_potassium_current)

%-------------------------------------------------------------------------------
% Computed variables
%-------------------------------------------------------------------------------

% GKbaraATP (milliS_per_microF) (in ATP_sensitive_potassium_current)
% g_K_ATP (milliS_per_microF) (in ATP_sensitive_potassium_current)
% i_K_ATP (microA_per_microF) (in ATP_sensitive_potassium_current)
% pATP (dimensionless) (in ATP_sensitive_potassium_current)
% E0_d (millivolt) (in L_type_Ca_channel_d_gate)
% alpha_d (per_ms) (in L_type_Ca_channel_d_gate)
% beta_d (per_ms) (in L_type_Ca_channel_d_gate)
% d_infinity (dimensionless) (in L_type_Ca_channel_d_gate)
% tau_d (ms) (in L_type_Ca_channel_d_gate)
% f_Ca (dimensionless) (in L_type_Ca_channel_f_Ca_gate)
% alpha_f (per_ms) (in L_type_Ca_channel_f_gate)
% beta_f (per_ms) (in L_type_Ca_channel_f_gate)
% f_infinity (dimensionless) (in L_type_Ca_channel_f_gate)
% tau_f (ms) (in L_type_Ca_channel_f_gate)
% I_CaCa (microA_per_microF) (in L_type_Ca_channel)
% I_CaK (microA_per_microF) (in L_type_Ca_channel)
% I_CaNa (microA_per_microF) (in L_type_Ca_channel)
% i_CaCa (microA_per_microF) (in L_type_Ca_channel)
% i_CaK (microA_per_microF) (in L_type_Ca_channel)
% i_CaNa (microA_per_microF) (in L_type_Ca_channel)
% i_Ca_L (microA_per_microF) (in L_type_Ca_channel)
% i_NaCa (microA_per_microF) (in Na_Ca_exchanger)
% b_inf (dimensionless) (in T_type_Ca_channel_b_gate)
% tau_b (ms) (in T_type_Ca_channel_b_gate)
% g_inf (dimensionless) (in T_type_Ca_channel_g_gate)
% tau_g (ms) (in T_type_Ca_channel_g_gate)
% i_Ca_T (microA_per_microF) (in T_type_Ca_channel)
% E_Ca (millivolt) (in calcium_background_current)
% i_Ca_b (microA_per_microF) (in calcium_background_current)
% G_rel (per_ms) (in calcium_dynamics)
% G_rel_Viswanathan (per_ms) (in calcium_dynamics)
% K_leak (per_ms) (in calcium_dynamics)
% RyRclose (dimensionless) (in calcium_dynamics)
% RyRopen (dimensionless) (in calcium_dynamics)
% i_leak (millimolar_per_ms) (in calcium_dynamics)
% i_rel (millimolar_per_ms) (in calcium_dynamics)
% i_tr (millimolar_per_ms) (in calcium_dynamics)
% i_up (millimolar_per_ms) (in calcium_dynamics)
% I_st (microA_per_microF) (in cell)
% dVdt (microA_per_microF) (in cell)
% alpha_h (per_ms) (in fast_sodium_current_h_gate)
% beta_h (per_ms) (in fast_sodium_current_h_gate)
% alpha_j (per_ms) (in fast_sodium_current_j_gate)
% beta_j (per_ms) (in fast_sodium_current_j_gate)
% E0_m (millivolt) (in fast_sodium_current_m_gate)
% alpha_m (per_ms) (in fast_sodium_current_m_gate)
% beta_m (per_ms) (in fast_sodium_current_m_gate)
% E_Na (millivolt) (in fast_sodium_current)
% i_Na (microA_per_microF) (in fast_sodium_current)
% V_JSR (micro_litre) (in geometry)
% V_NSR (micro_litre) (in geometry)
% V_myo (micro_litre) (in geometry)
% volume (micro_litre) (in geometry)
% I_ns_K (microA_per_microF) (in non_specific_calcium_activated_current)
% I_ns_Na (microA_per_microF) (in non_specific_calcium_activated_current)
% P_ns_Ca (litre_per_farad_millisecond) (in non_specific_calcium_activated_current)
% i_ns_Ca (microA_per_microF) (in non_specific_calcium_activated_current)
% i_ns_K (microA_per_microF) (in non_specific_calcium_activated_current)
% i_ns_Na (microA_per_microF) (in non_specific_calcium_activated_current)
% Kp (dimensionless) (in plateau_potassium_current)
% i_Kp (microA_per_microF) (in plateau_potassium_current)
% tau_xr (ms) (in rapid_delayed_rectifier_potassium_current_xr_gate)
% xr_infinity (dimensionless) (in rapid_delayed_rectifier_potassium_current_xr_gate)
% Rect (dimensionless) (in rapid_delayed_rectifier_potassium_current)
% g_Kr (milliS_per_microF) (in rapid_delayed_rectifier_potassium_current)
% i_Kr (microA_per_microF) (in rapid_delayed_rectifier_potassium_current)
% i_p_Ca (microA_per_microF) (in sarcolemmal_calcium_pump)
% tau_xs1 (ms) (in slow_delayed_rectifier_potassium_current_xs1_gate)
% xs1_infinity (dimensionless) (in slow_delayed_rectifier_potassium_current_xs1_gate)
% tau_xs2 (ms) (in slow_delayed_rectifier_potassium_current_xs2_gate)
% xs2_infinity (dimensionless) (in slow_delayed_rectifier_potassium_current_xs2_gate)
% E_Ks (millivolt) (in slow_delayed_rectifier_potassium_current)
% g_Ks (milliS_per_microF) (in slow_delayed_rectifier_potassium_current)
% i_Ks (microA_per_microF) (in slow_delayed_rectifier_potassium_current)
% i_K_Na (microA_per_microF) (in sodium_activated_potassium_current)
% pona (dimensionless) (in sodium_activated_potassium_current)
% pov (dimensionless) (in sodium_activated_potassium_current)
% i_Na_b (microA_per_microF) (in sodium_background_current)
% f_NaK (dimensionless) (in sodium_potassium_pump)
% i_NaK (microA_per_microF) (in sodium_potassium_pump)
% sigma (dimensionless) (in sodium_potassium_pump)
% K1_infinity (dimensionless) (in time_independent_potassium_current_K1_gate)
% alpha_K1 (per_ms) (in time_independent_potassium_current_K1_gate)
% beta_K1 (per_ms) (in time_independent_potassium_current_K1_gate)
% E_K (millivolt) (in time_independent_potassium_current)
% g_K1 (milliS_per_microF) (in time_independent_potassium_current)
% i_K1 (microA_per_microF) (in time_independent_potassium_current)
% alpha_ydv (per_ms) (in transient_outward_current_ydv_gate)
% beta_ydv (per_ms) (in transient_outward_current_ydv_gate)
% tau_ydv (ms) (in transient_outward_current_ydv_gate)
% ydv_ss (dimensionless) (in transient_outward_current_ydv_gate)
% alpha_zdv (per_ms) (in transient_outward_current_zdv_gate)
% beta_zdv (per_ms) (in transient_outward_current_zdv_gate)
% tau_zdv (ms) (in transient_outward_current_zdv_gate)
% zdv_ss (dimensionless) (in transient_outward_current_zdv_gate)
% g_to (milliS_per_microF) (in transient_outward_current)
% i_to (microA_per_microF) (in transient_outward_current)
% rvdv (dimensionless) (in transient_outward_current)

%-------------------------------------------------------------------------------
% Computation
%-------------------------------------------------------------------------------

% time (ms)

g_K_ATP = i_K_ATP_on*0.000193/nicholsarea;
pATP = 1.0/(1.0+(ATPi/kATP)^hATP);
GKbaraATP = g_K_ATP*pATP*(Ko/4.0)^nATP;
E_K = R*T/F*log(Ko/Y(19));
i_K_ATP = GKbaraATP*(Y(15)-E_K);
I_CaCa = P_Ca*2.0^2.0*Y(15)*F^2.0/(R*T)*(gamma_Cai*Y(10)*exp(2.0*Y(15)*F/(R*T))-gamma_Cao*Cao)/(exp(2.0*Y(15)*F/(R*T))-1.0);
I_CaNa = P_Na*1.0^2.0*Y(15)*F^2.0/(R*T)*(gamma_Nai*Y(20)*exp(1.0*Y(15)*F/(R*T))-gamma_Nao*Nao)/(exp(1.0*Y(15)*F/(R*T))-1.0);
I_CaK = P_K*1.0^2.0*Y(15)*F^2.0/(R*T)*(gamma_Ki*Y(19)*exp(1.0*Y(15)*F/(R*T))-gamma_Ko*Ko)/(exp(1.0*Y(15)*F/(R*T))-1.0);
f_Ca = 1.0/(1.0+Y(10)/Km_Ca);
i_CaCa = Y(1)*Y(2)*f_Ca*I_CaCa;
i_CaNa = Y(1)*Y(2)*f_Ca*I_CaNa;
i_CaK = Y(1)*Y(2)*f_Ca*I_CaK;
i_Ca_L = i_CaCa+i_CaK+i_CaNa;
E0_d = Y(15)+10.0;
d_infinity = 1.0/(1.0+exp(-E0_d/6.24));

if (abs(E0_d) < 1.0e-5)
   tau_d = 1.0/(0.035*6.24*2.0);
else
   tau_d = 1.0*d_infinity*(1.0-exp(-E0_d/6.24))/(0.035*E0_d);
end;

alpha_d = d_infinity/tau_d;
beta_d = (1.0-d_infinity)/tau_d;
dY(1, 1) = alpha_d*(1.0-Y(1))-beta_d*Y(1);
f_infinity = 1.0/(1.0+exp((Y(15)+35.06)/8.6))+0.6/(1.0+exp((50.0-Y(15))/20.0));
tau_f = 1.0/(0.0197*exp(-(0.0337*(Y(15)+10.0))^2.0)+0.02);
alpha_f = f_infinity/tau_f;
beta_f = (1.0-f_infinity)/tau_f;
dY(2, 1) = alpha_f*(1.0-Y(2))-beta_f*Y(2);
i_NaCa = c1*exp((gamma-1.0)*Y(15)*F/(R*T))*(exp(Y(15)*F/(R*T))*Y(20)^3.0*Cao-Nao^3.0*Y(10))/(1.0+c2*exp((gamma-1.0)*Y(15)*F/(R*T))*(exp(Y(15)*F/(R*T))*Y(20)^3.0*Cao+Nao^3.0*Y(10)));
E_Ca = R*T/(2.0*F)*log(Cao/Y(10));
i_Ca_T = g_CaT*Y(3)*Y(3)*Y(4)*(Y(15)-E_Ca);
b_inf = 1.0/(1.0+exp(-(Y(15)+14.0)/10.8));
tau_b = 3.7+6.1/(1.0+exp((Y(15)+25.0)/4.5));
dY(3, 1) = (b_inf-Y(3))/tau_b;
g_inf = 1.0/(1.0+exp((Y(15)+60.0)/5.6));

if (Y(15) <= 0.0)
   tau_g = -0.875*Y(15)+12.0;
else
   tau_g = 12.0;
end;

dY(4, 1) = (g_inf-Y(4))/tau_g;
i_Ca_b = g_Cab*(Y(15)-E_Ca);
E_Na = R*T/F*log(Nao/Y(20));
i_Na = g_Na*Y(18)^3.0*Y(16)*Y(17)*(Y(15)-E_Na);
g_Kr = G_Kr*sqrt(Ko/5.4);
Rect = 1.0/(1.0+exp((Y(15)+9.0)/22.4));
i_Kr = g_Kr*Y(21)*Rect*(Y(15)-E_K);
g_Ks = G_Ks*(1.0+0.6/(1.0+(3.8e-5/Y(10))^1.4));
E_Ks = R*T/F*log((Ko+PNaK*Nao)/(Y(19)+PNaK*Y(20)));
i_Ks = g_Ks*Y(22)*Y(23)*(Y(15)-E_Ks);
pona = 0.85/(1.0+(kdKNa/Y(20))^nKNa);
pov = 0.8-0.65/(1.0+exp((Y(15)+125.0)/15.0));
i_K_Na = g_K_Na*pona*pov*(Y(15)-E_K);
g_to = 0.0*0.5;
rvdv = exp(Y(15)/100.0);
i_to = g_to*Y(25)^3.0*Y(24)*rvdv*(Y(15)-E_K);
g_K1 = G_K1*sqrt(Ko/5.4);
alpha_K1 = 1.02/(1.0+exp(0.2385*(Y(15)-E_K-59.215)));
beta_K1 = 1.0*(0.49124*exp(0.08032*(Y(15)-E_K+5.476))+exp(0.06175*(Y(15)-E_K-594.31)))/(1.0+exp(-0.5143*(Y(15)-E_K+4.753)));
K1_infinity = alpha_K1/(alpha_K1+beta_K1);
i_K1 = g_K1*K1_infinity*(Y(15)-E_K);
Kp = 1.0/(1.0+exp((7.488-Y(15))/5.98));
i_Kp = g_Kp*Kp*(Y(15)-E_K);
i_p_Ca = I_pCa*Y(10)/(K_mpCa+Y(10));
i_Na_b = g_Nab*(Y(15)-E_Na);
sigma = 1.0/7.0*(exp(Nao/67.3)-1.0);
f_NaK = 1.0/(1.0+0.1245*exp(-0.1*Y(15)*F/(R*T))+0.0365*sigma*exp(-Y(15)*F/(R*T)));
i_NaK = I_NaK*f_NaK*1.0/(1.0+(K_mNai/Y(20))^2.0)*Ko/(Ko+K_mKo);
P_ns_Ca = 1.75e-7;
I_ns_Na = P_ns_Ca*1.0^2.0*Y(15)*F^2.0/(R*T)*(gamma_Nai*Y(20)*exp(1.0*Y(15)*F/(R*T))-gamma_Nao*Nao)/(exp(1.0*Y(15)*F/(R*T))-1.0);
i_ns_Na = I_ns_Na*1.0/(1.0+(K_m_ns_Ca/Y(10))^3.0);
I_ns_K = P_ns_Ca*1.0^2.0*Y(15)*F^2.0/(R*T)*(gamma_Ki*Y(19)*exp(1.0*Y(15)*F/(R*T))-gamma_Ko*Ko)/(exp(1.0*Y(15)*F/(R*T))-1.0);
i_ns_K = I_ns_K*1.0/(1.0+(K_m_ns_Ca/Y(10))^3.0);
i_ns_Ca = i_ns_Na+i_ns_K;

if ((time >= stim_start) && (time <= stim_end) && (time-stim_start-floor((time-stim_start)/stim_period)*stim_period <= stim_duration))
   I_st = stim_amplitude;
else
   I_st = 0.0;
end;

dVdt = -(i_Na+i_Ca_L+i_Ca_T+i_Kr+i_Ks+i_K_Na+i_K_ATP+i_to+i_K1+i_Kp+i_NaCa+i_p_Ca+i_Na_b+i_Ca_b+i_NaK+i_ns_Ca+I_st);

if (dVdt > 150.0)
   dY(5, 1) = 100.0*(1.0-Y(5))-0.5*Y(5);
else
   dY(5, 1) = -0.5*Y(5);
end;

if ((Y(5) < 0.2) && (Y(5) > 0.18))
   dY(6, 1) = 100.0*(1.0-Y(6))-0.5*Y(6);
else
   dY(6, 1) = -0.5*Y(6);
end;

if ((Y(5) < 0.2) && (Y(5) > 0.18))
   dY(7, 1) = 100.0*(1.0-Y(7))-0.5*Y(7);
else
   dY(7, 1) = -0.01*Y(7);
end;

volume = pi*preplength*radius^2.0;
V_myo = 0.68*volume;

if (Y(5) > 0.2)
   dY(11, 1) = -1.0*A_cap*(i_CaCa+i_Ca_T-2.0*i_NaCa+i_p_Ca+i_Ca_b)/(2.0*V_myo*F);
elseif ((Y(6) > 0.01) && (Y(5) <= 0.2))
   dY(11, 1) = 0.0;
else
   dY(11, 1) = -0.5*Y(11);
end;

if ((1.0/(1.0+K_mCSQN/Y(8)) > CSQNthresh) && (Y(14) < 0.37) && (Y(7) < 0.37))
   dY(12, 1) = 50.0*(1.0-Y(12));
else
   dY(12, 1) = -0.5*Y(12);
end;

if ((Y(12) > Logicthresh) && (Y(13) < Logicthresh))
   dY(13, 1) = 50.0*(1.0-Y(13));
else
   dY(13, 1) = -0.5*Y(13);
end;

if ((Y(12) > Logicthresh) && (Y(14) < Logicthresh))
   dY(14, 1) = 50.0*(1.0-Y(14));
else
   dY(14, 1) = -0.01*Y(14);
end;

if (Y(11) > delta_Ca_ith)
   G_rel_Viswanathan = G_rel_max*(Y(11)-delta_Ca_ith)/(K_mrel+Y(11)-delta_Ca_ith)*(1.0-Y(6))*Y(6);
elseif ((Y(11) <= delta_Ca_ith) && (Y(13) > 0.0))
   G_rel_Viswanathan = G_rel_overload*(1.0-Y(13))*Y(13);
else
   G_rel_Viswanathan = 0.0;
end;

RyRopen = 1.0/(1.0+exp(2.0)*Y(6));
RyRclose = 1.0-RyRopen;
G_rel = RyRopen*RyRclose*150.0/(1.0+exp((i_CaCa+i_Ca_T-2.0*i_NaCa+i_p_Ca+i_Ca_b+5.0)/0.9));
i_rel = G_rel*(Y(8)-Y(10));
i_up = I_up*Y(10)/(Y(10)+K_mup);
K_leak = I_up/Ca_NSR_max;
i_leak = K_leak*Y(9);
i_tr = (Y(9)-Y(8))/tau_tr;
dY(8, 1) = 1.0/(1.0+CSQN_max*K_mCSQN/(K_mCSQN+Y(8))^2.0)*(i_tr-i_rel);
V_JSR = 0.0048*volume;
V_NSR = 0.0552*volume;
dY(9, 1) = -i_tr*V_JSR/V_NSR-i_leak+i_up;
dY(10, 1) = 1.0/(1.0+CMDN_max*K_mCMDN/(K_mCMDN+Y(10))^2.0+Tn_max*K_mTn/(K_mTn+Y(10))^2.0)*(-1.0*A_cap*(i_CaCa+i_Ca_T-2.0*i_NaCa+i_p_Ca+i_Ca_b)/(2.0*V_myo*F)+i_rel*V_JSR/V_myo+(i_leak-i_up)*V_NSR/V_myo);
dY(15, 1) = dVdt;

if (Y(15) < -40.0)
   alpha_h = 0.135*exp((80.0+Y(15))/-6.8);
else
   alpha_h = 0.0;
end;

if (Y(15) < -40.0)
   beta_h = 3.56*exp(0.079*Y(15))+310000.0*exp(0.35*Y(15));
else
   beta_h = 1.0/(0.13*(1.0+exp((Y(15)+10.66)/-11.1)));
end;

dY(16, 1) = alpha_h*(1.0-Y(16))-beta_h*Y(16);

if (Y(15) < -40.0)
   alpha_j = -(127140.0*exp(0.2444*Y(15))+3.474e-5*exp(-0.04391*Y(15)))*(Y(15)+37.78)/(1.0+exp(0.311*(Y(15)+79.23)));
else
   alpha_j = 0.0;
end;

if (Y(15) < -40.0)
   beta_j = 0.1212*exp(-0.01052*Y(15))/(1.0+exp(-0.1378*(Y(15)+40.14)));
else
   beta_j = 0.3*exp(-2.535e-7*Y(15))/(1.0+exp(-0.1*(Y(15)+32.0)));
end;

dY(17, 1) = alpha_j*(1.0-Y(17))-beta_j*Y(17);
E0_m = Y(15)+47.13;

if (abs(E0_m) >= delta_m)
   alpha_m = 0.32*E0_m/(1.0-exp(-0.1*E0_m));
else
   alpha_m = 3.2;
end;

beta_m = 0.08*exp(-Y(15)/11.0);
dY(18, 1) = alpha_m*(1.0-Y(18))-beta_m*Y(18);
dY(20, 1) = -1.0*(i_Na+i_CaNa+i_Na_b+i_ns_Na+i_NaCa*3.0+i_NaK*3.0)*A_cap/(V_myo*F);
dY(19, 1) = -1.0*(I_st+i_CaK+i_Kr+i_Ks+i_K1+i_Kp+i_K_Na+i_K_ATP+i_to+i_ns_K+-i_NaK*2.0)*A_cap/(V_myo*F);
xr_infinity = 1.0/(1.0+exp(-(Y(15)+21.5)/7.5));
tau_xr = 1.0/(0.00138*(Y(15)+14.2)/(1.0-exp(-0.123*(Y(15)+14.2)))+0.00061*(Y(15)+38.9)/(exp(0.145*(Y(15)+38.9))-1.0));
dY(21, 1) = (xr_infinity-Y(21))/tau_xr;
xs1_infinity = 1.0/(1.0+exp(-(Y(15)-1.5)/16.7));
tau_xs1 = 1.0/(7.19e-5*(Y(15)+30.0)/(1.0-exp(-0.148*(Y(15)+30.0)))+0.000131*(Y(15)+30.0)/(exp(0.0687*(Y(15)+30.0))-1.0));
dY(22, 1) = (xs1_infinity-Y(22))/tau_xs1;
xs2_infinity = 1.0/(1.0+exp(-(Y(15)-1.5)/16.7));
tau_xs2 = 4.0/(7.19e-5*(Y(15)+30.0)/(1.0-exp(-0.148*(Y(15)+30.0)))+0.000131*(Y(15)+30.0)/(exp(0.0687*(Y(15)+30.0))-1.0));
dY(23, 1) = (xs2_infinity-Y(23))/tau_xs2;
alpha_ydv = 0.015/(1.0+exp((Y(15)+60.0)/5.0));
beta_ydv = 0.1*exp((Y(15)+25.0)/5.0)/(1.0+exp((Y(15)+25.0)/5.0));
tau_ydv = 1.0/(alpha_ydv+beta_ydv);
ydv_ss = alpha_ydv/(alpha_ydv+beta_ydv);
dY(24, 1) = (ydv_ss-Y(24))/tau_ydv;
alpha_zdv = 10.0*exp((Y(15)-40.0)/25.0)/(1.0+exp((Y(15)-40.0)/25.0));
beta_zdv = 10.0*exp(-(Y(15)+90.0)/25.0)/(1.0+exp(-(Y(15)+90.0)/25.0));
tau_zdv = 1.0/(alpha_zdv+beta_zdv);
zdv_ss = alpha_zdv/(alpha_zdv+beta_zdv);
dY(25, 1) = (zdv_ss-Y(25))/tau_zdv;

%===============================================================================
% End of file
%===============================================================================
