%===============================================================================
% CellML file:   D:\Desktop\Models\bondarenko_szigeti_bett_kim_rasmusson_2004_apical.cellml
% CellML model:  bondarenko_2004_apical
% Date and time: 17/06/2015 at 23:16:51
%-------------------------------------------------------------------------------
% Conversion from CellML 1.0 to MATLAB (init) was done using COR (0.9.31.1409)
%    Copyright 2002-2015 Dr Alan Garny
%    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
%-------------------------------------------------------------------------------
% http://www.cellml.org/
%===============================================================================

function dY = bondarenko_szigeti_bett_kim_rasmusson_2004_apical(time, Y)

%-------------------------------------------------------------------------------
% Initial conditions
%-------------------------------------------------------------------------------

% Y = [0.124216e-3, 0.578679e-8, 0.119816e-12, 0.497923e-18, 0.345847e-13, 0.185106e-13, 0.930308e-18, 125.29, 11.2684, 1299.5, 1299.5, 0.115001, 0.115001, 0.0, 0.279132e-3, 0.020752, 0.673345e-6, 0.155787e-8, 0.0113879, 0.34278, 0.153176e-3, 0.713483e-6, 0.265563e-2, 0.999977, -82.4202, 0.417069e-3, 1.0, 143720.0, 0.992513e-3, 0.641229e-3, 0.319129e-4, 0.175298e-3, 0.16774e-3, 0.149102e-4, 0.951726e-10, 0.262753e-3, 0.417069e-3, 0.998543, 14237.1, 0.417069e-3, 0.998543];

% YNames = {'C2', 'C3', 'C4', 'I1', 'I2', 'I3', 'O', 'HTRPN_Ca', 'LTRPN_Ca', 'CaJSR', 'CaNSR', 'Cai', 'Cass', 'P_RyR', 'C_Na1', 'C_Na2', 'I1_Na', 'I2_Na', 'IC_Na2', 'IC_Na3', 'IF_Na', 'O_Na', 'ato_f', 'ito_f', 'V', 'aKss', 'iKss', 'Ki', 'C_K1', 'C_K2', 'I_K', 'O_K', 'P_C2', 'P_O1', 'P_O2', 'nKs', 'ato_s', 'ito_s', 'Nai', 'aur', 'iur'};
% YUnits = {'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'micromolar', 'micromolar', 'micromolar', 'micromolar', 'micromolar', 'micromolar', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'millivolt', 'dimensionless', 'dimensionless', 'micromolar', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'micromolar', 'dimensionless', 'dimensionless'};
% YComponents = {'L_type_calcium_current', 'L_type_calcium_current', 'L_type_calcium_current', 'L_type_calcium_current', 'L_type_calcium_current', 'L_type_calcium_current', 'L_type_calcium_current', 'calcium_buffering', 'calcium_buffering', 'calcium_concentration', 'calcium_concentration', 'calcium_concentration', 'calcium_concentration', 'calcium_fluxes', 'fast_sodium_current', 'fast_sodium_current', 'fast_sodium_current', 'fast_sodium_current', 'fast_sodium_current', 'fast_sodium_current', 'fast_sodium_current', 'fast_sodium_current', 'fast_transient_outward_potassium_current', 'fast_transient_outward_potassium_current', 'membrane', 'non_inactivating_steady_state_potassium_current', 'non_inactivating_steady_state_potassium_current', 'potassium_concentration', 'rapid_delayed_rectifier_potassium_current', 'rapid_delayed_rectifier_potassium_current', 'rapid_delayed_rectifier_potassium_current', 'rapid_delayed_rectifier_potassium_current', 'ryanodine_receptors', 'ryanodine_receptors', 'ryanodine_receptors', 'slow_delayed_rectifier_potassium_current', 'slow_transient_outward_potassium_current', 'slow_transient_outward_potassium_current', 'sodium_concentration', 'ultra_rapidly_activating_delayed_rectifier_potassium_current', 'ultra_rapidly_activating_delayed_rectifier_potassium_current'};

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: C2 (dimensionless) (in L_type_calcium_current)
% 2: C3 (dimensionless) (in L_type_calcium_current)
% 3: C4 (dimensionless) (in L_type_calcium_current)
% 4: I1 (dimensionless) (in L_type_calcium_current)
% 5: I2 (dimensionless) (in L_type_calcium_current)
% 6: I3 (dimensionless) (in L_type_calcium_current)
% 7: O (dimensionless) (in L_type_calcium_current)
% 8: HTRPN_Ca (micromolar) (in calcium_buffering)
% 9: LTRPN_Ca (micromolar) (in calcium_buffering)
% 10: CaJSR (micromolar) (in calcium_concentration)
% 11: CaNSR (micromolar) (in calcium_concentration)
% 12: Cai (micromolar) (in calcium_concentration)
% 13: Cass (micromolar) (in calcium_concentration)
% 14: P_RyR (dimensionless) (in calcium_fluxes)
% 15: C_Na1 (dimensionless) (in fast_sodium_current)
% 16: C_Na2 (dimensionless) (in fast_sodium_current)
% 17: I1_Na (dimensionless) (in fast_sodium_current)
% 18: I2_Na (dimensionless) (in fast_sodium_current)
% 19: IC_Na2 (dimensionless) (in fast_sodium_current)
% 20: IC_Na3 (dimensionless) (in fast_sodium_current)
% 21: IF_Na (dimensionless) (in fast_sodium_current)
% 22: O_Na (dimensionless) (in fast_sodium_current)
% 23: ato_f (dimensionless) (in fast_transient_outward_potassium_current)
% 24: ito_f (dimensionless) (in fast_transient_outward_potassium_current)
% 25: V (millivolt) (in membrane)
% 26: aKss (dimensionless) (in non_inactivating_steady_state_potassium_current)
% 27: iKss (dimensionless) (in non_inactivating_steady_state_potassium_current)
% 28: Ki (micromolar) (in potassium_concentration)
% 29: C_K1 (dimensionless) (in rapid_delayed_rectifier_potassium_current)
% 30: C_K2 (dimensionless) (in rapid_delayed_rectifier_potassium_current)
% 31: I_K (dimensionless) (in rapid_delayed_rectifier_potassium_current)
% 32: O_K (dimensionless) (in rapid_delayed_rectifier_potassium_current)
% 33: P_C2 (dimensionless) (in ryanodine_receptors)
% 34: P_O1 (dimensionless) (in ryanodine_receptors)
% 35: P_O2 (dimensionless) (in ryanodine_receptors)
% 36: nKs (dimensionless) (in slow_delayed_rectifier_potassium_current)
% 37: ato_s (dimensionless) (in slow_transient_outward_potassium_current)
% 38: ito_s (dimensionless) (in slow_transient_outward_potassium_current)
% 39: Nai (micromolar) (in sodium_concentration)
% 40: aur (dimensionless) (in ultra_rapidly_activating_delayed_rectifier_potassium_current)
% 41: iur (dimensionless) (in ultra_rapidly_activating_delayed_rectifier_potassium_current)

%-------------------------------------------------------------------------------
% Constants
%-------------------------------------------------------------------------------

E_CaL = 63.0;   % millivolt (in L_type_calcium_current)
Kpc_half = 20.0;   % micromolar (in L_type_calcium_current)
Kpc_max = 0.23324;   % per_millisecond (in L_type_calcium_current)
Kpcb = 0.0005;   % per_millisecond (in L_type_calcium_current)
g_CaL = 0.1729;   % milliS_per_microF (in L_type_calcium_current)
i_CaL_max = 7.0;   % picoA_per_picoF (in L_type_calcium_current)
E_Cl = -40.0;   % millivolt (in calcium_activated_chloride_current)
Km_Cl = 10.0;   % micromolar (in calcium_activated_chloride_current)
g_ClCa = 10.0;   % milliS_per_microF (in calcium_activated_chloride_current)
g_Cab = 0.000367;   % milliS_per_microF (in calcium_background_current)
HTRPN_tot = 140.0;   % micromolar (in calcium_buffering)
LTRPN_tot = 70.0;   % micromolar (in calcium_buffering)
CMDN_tot = 50.0;   % micromolar (in calcium_concentration)
CSQN_tot = 15000.0;   % micromolar (in calcium_concentration)
Km_CMDN = 0.238;   % micromolar (in calcium_concentration)
Km_CSQN = 800.0;   % micromolar (in calcium_concentration)
Km_up = 0.5;   % micromolar (in calcium_fluxes)
k_minus_htrpn = 3.2e-5;   % per_millisecond (in calcium_fluxes)
k_minus_ltrpn = 0.0196;   % per_millisecond (in calcium_fluxes)
k_plus_htrpn = 0.00237;   % per_micromolar_millisecond (in calcium_fluxes)
k_plus_ltrpn = 0.0327;   % per_micromolar_millisecond (in calcium_fluxes)
tau_tr = 20.0;   % millisecond (in calcium_fluxes)
tau_xfer = 8.0;   % millisecond (in calcium_fluxes)
v1 = 4.5;   % per_millisecond (in calcium_fluxes)
v2 = 1.74e-5;   % per_millisecond (in calcium_fluxes)
v3 = 0.45;   % micromolar_per_millisecond (in calcium_fluxes)
Km_pCa = 0.5;   % micromolar (in calcium_pump_current)
i_pCa_max = 1.0;   % picoA_per_picoF (in calcium_pump_current)
g_Na = 13.0;   % milliS_per_microF (in fast_sodium_current)
g_Kto_f = 0.4067;   % milliS_per_microF (in fast_transient_outward_potassium_current)
Acap = 1.534e-4;   % cm2 (in membrane)
Cao = 1800.0;   % micromolar (in membrane)
Cm = 1.0;   % microF_per_cm2 (in membrane)
F = 96.5;   % coulomb_per_millimole (in membrane)
Ko = 5400.0;   % micromolar (in membrane)
Nao = 140000.0;   % micromolar (in membrane)
R = 8.314;   % joule_per_mole_kelvin (in membrane)
T = 298.0;   % kelvin (in membrane)
VJSR = 0.12e-6;   % microlitre (in membrane)
VNSR = 2.098e-6;   % microlitre (in membrane)
Vmyo = 25.84e-6;   % microlitre (in membrane)
Vss = 1.485e-9;   % microlitre (in membrane)
stim_amplitude = -80.0;   % picoA_per_picoF (in membrane)
stim_duration = 0.5;   % millisecond (in membrane)
stim_end = 100000.0;   % millisecond (in membrane)
stim_period = 71.43;   % millisecond (in membrane)
stim_start = 20.0;   % millisecond (in membrane)
g_Kss = 0.05;   % milliS_per_microF (in non_inactivating_steady_state_potassium_current)
g_Kr = 0.078;   % milliS_per_microF (in rapid_delayed_rectifier_potassium_current)
kb = 0.036778;   % per_millisecond (in rapid_delayed_rectifier_potassium_current)
kf = 0.023761;   % per_millisecond (in rapid_delayed_rectifier_potassium_current)
k_minus_a = 0.07125;   % per_millisecond (in ryanodine_receptors)
k_minus_b = 0.965;   % per_millisecond (in ryanodine_receptors)
k_minus_c = 0.0008;   % per_millisecond (in ryanodine_receptors)
k_plus_a = 0.006075;   % micromolar4_per_millisecond (in ryanodine_receptors)
k_plus_b = 0.00405;   % micromolar3_per_millisecond (in ryanodine_receptors)
k_plus_c = 0.009;   % per_millisecond (in ryanodine_receptors)
m = 3.0;   % dimensionless (in ryanodine_receptors)
n = 4.0;   % dimensionless (in ryanodine_receptors)
g_Ks = 0.00575;   % milliS_per_microF (in slow_delayed_rectifier_potassium_current)
g_Kto_s = 0.0;   % milliS_per_microF (in slow_transient_outward_potassium_current)
g_Nab = 0.0026;   % milliS_per_microF (in sodium_background_current)
K_mCa = 1380.0;   % micromolar (in sodium_calcium_exchange_current)
K_mNa = 87500.0;   % micromolar (in sodium_calcium_exchange_current)
eta = 0.35;   % dimensionless (in sodium_calcium_exchange_current)
k_NaCa = 292.8;   % picoA_per_picoF (in sodium_calcium_exchange_current)
k_sat = 0.1;   % dimensionless (in sodium_calcium_exchange_current)
Km_Ko = 1500.0;   % micromolar (in sodium_potassium_pump_current)
Km_Nai = 21000.0;   % micromolar (in sodium_potassium_pump_current)
i_NaK_max = 0.88;   % picoA_per_picoF (in sodium_potassium_pump_current)
g_Kur = 0.16;   % milliS_per_microF (in ultra_rapidly_activating_delayed_rectifier_potassium_current)

%-------------------------------------------------------------------------------
% Computed variables
%-------------------------------------------------------------------------------

% C1 (dimensionless) (in L_type_calcium_current)
% Kpcf (per_millisecond) (in L_type_calcium_current)
% alpha (per_millisecond) (in L_type_calcium_current)
% beta (per_millisecond) (in L_type_calcium_current)
% gamma (per_millisecond) (in L_type_calcium_current)
% i_CaL (picoA_per_picoF) (in L_type_calcium_current)
% O_ClCa (dimensionless) (in calcium_activated_chloride_current)
% i_ClCa (picoA_per_picoF) (in calcium_activated_chloride_current)
% E_CaN (millivolt) (in calcium_background_current)
% i_Cab (picoA_per_picoF) (in calcium_background_current)
% BJSR (dimensionless) (in calcium_concentration)
% Bi (dimensionless) (in calcium_concentration)
% Bss (dimensionless) (in calcium_concentration)
% J_leak (micromolar_per_millisecond) (in calcium_fluxes)
% J_rel (micromolar_per_millisecond) (in calcium_fluxes)
% J_tr (micromolar_per_millisecond) (in calcium_fluxes)
% J_trpn (micromolar_per_millisecond) (in calcium_fluxes)
% J_up (micromolar_per_millisecond) (in calcium_fluxes)
% J_xfer (micromolar_per_millisecond) (in calcium_fluxes)
% i_pCa (picoA_per_picoF) (in calcium_pump_current)
% C_Na3 (dimensionless) (in fast_sodium_current)
% E_Na (millivolt) (in fast_sodium_current)
% alpha_Na11 (per_millisecond) (in fast_sodium_current)
% alpha_Na12 (per_millisecond) (in fast_sodium_current)
% alpha_Na13 (per_millisecond) (in fast_sodium_current)
% alpha_Na2 (per_millisecond) (in fast_sodium_current)
% alpha_Na3 (per_millisecond) (in fast_sodium_current)
% alpha_Na4 (per_millisecond) (in fast_sodium_current)
% alpha_Na5 (per_millisecond) (in fast_sodium_current)
% beta_Na11 (per_millisecond) (in fast_sodium_current)
% beta_Na12 (per_millisecond) (in fast_sodium_current)
% beta_Na13 (per_millisecond) (in fast_sodium_current)
% beta_Na2 (per_millisecond) (in fast_sodium_current)
% beta_Na3 (per_millisecond) (in fast_sodium_current)
% beta_Na4 (per_millisecond) (in fast_sodium_current)
% beta_Na5 (per_millisecond) (in fast_sodium_current)
% i_Na (picoA_per_picoF) (in fast_sodium_current)
% E_K (millivolt) (in fast_transient_outward_potassium_current)
% alpha_a (per_millisecond) (in fast_transient_outward_potassium_current)
% alpha_i_1 (per_millisecond) (alpha_i in fast_transient_outward_potassium_current)
% beta_a (per_millisecond) (in fast_transient_outward_potassium_current)
% beta_i_1 (per_millisecond) (beta_i in fast_transient_outward_potassium_current)
% i_Kto_f (picoA_per_picoF) (in fast_transient_outward_potassium_current)
% i_stim (picoA_per_picoF) (in membrane)
% i_Kss (picoA_per_picoF) (in non_inactivating_steady_state_potassium_current)
% tau_Kss (millisecond) (in non_inactivating_steady_state_potassium_current)
% C_K0 (dimensionless) (in rapid_delayed_rectifier_potassium_current)
% alpha_a0 (per_millisecond) (in rapid_delayed_rectifier_potassium_current)
% alpha_a1 (per_millisecond) (in rapid_delayed_rectifier_potassium_current)
% alpha_i_2 (per_millisecond) (alpha_i in rapid_delayed_rectifier_potassium_current)
% beta_a0 (per_millisecond) (in rapid_delayed_rectifier_potassium_current)
% beta_a1 (per_millisecond) (in rapid_delayed_rectifier_potassium_current)
% beta_i_2 (per_millisecond) (beta_i in rapid_delayed_rectifier_potassium_current)
% i_Kr (picoA_per_picoF) (in rapid_delayed_rectifier_potassium_current)
% P_C1 (dimensionless) (in ryanodine_receptors)
% alpha_n (per_millisecond) (in slow_delayed_rectifier_potassium_current)
% beta_n (per_millisecond) (in slow_delayed_rectifier_potassium_current)
% i_Ks (picoA_per_picoF) (in slow_delayed_rectifier_potassium_current)
% ass (dimensionless) (in slow_transient_outward_potassium_current)
% i_Kto_s (picoA_per_picoF) (in slow_transient_outward_potassium_current)
% iss (dimensionless) (in slow_transient_outward_potassium_current)
% tau_ta_s (millisecond) (in slow_transient_outward_potassium_current)
% tau_ti_s (millisecond) (in slow_transient_outward_potassium_current)
% i_Nab (picoA_per_picoF) (in sodium_background_current)
% i_NaCa (picoA_per_picoF) (in sodium_calcium_exchange_current)
% f_NaK (dimensionless) (in sodium_potassium_pump_current)
% i_NaK (picoA_per_picoF) (in sodium_potassium_pump_current)
% sigma (dimensionless) (in sodium_potassium_pump_current)
% i_K1 (picoA_per_picoF) (in time_independent_potassium_current)
% i_Kur (picoA_per_picoF) (in ultra_rapidly_activating_delayed_rectifier_potassium_current)
% tau_aur (millisecond) (in ultra_rapidly_activating_delayed_rectifier_potassium_current)
% tau_iur (millisecond) (in ultra_rapidly_activating_delayed_rectifier_potassium_current)

%-------------------------------------------------------------------------------
% Computation
%-------------------------------------------------------------------------------

% time (millisecond)

i_CaL = g_CaL*Y(7)*(Y(25)-E_CaL);
alpha = 0.4*exp((Y(25)+12.0)/10.0)*(1.0+0.7*exp(-(Y(25)+40.0)^2.0/10.0)-0.75*exp(-(Y(25)+20.0)^2.0/400.0))/(1.0+0.12*exp((Y(25)+12.0)/10.0));
Kpcf = 13.0*(1.0-exp(-(Y(25)+14.5)^2.0/100.0));
beta = 0.05*exp(-(Y(25)+12.0)/13.0);
gamma = Kpc_max*Y(13)/(Kpc_half+Y(13));
dY(7, 1) = alpha*Y(3)+Kpcb*Y(4)+0.001*(alpha*Y(5)-Kpcf*Y(7))-(4.0*beta*Y(7)+gamma*Y(7));
C1 = 1.0-(Y(7)+Y(1)+Y(2)+Y(3)+Y(4)+Y(5)+Y(6));
dY(1, 1) = 4.0*alpha*C1+2.0*beta*Y(2)-(beta*Y(1)+3.0*alpha*Y(1));
dY(2, 1) = 3.0*alpha*Y(1)+3.0*beta*Y(3)-(2.0*beta*Y(2)+2.0*alpha*Y(2));
dY(3, 1) = 2.0*alpha*Y(2)+4.0*beta*Y(7)+0.01*(4.0*Kpcb*beta*Y(4)-alpha*gamma*Y(3))+0.002*(4.0*beta*Y(5)-Kpcf*Y(3))+4.0*beta*Kpcb*Y(6)-(3.0*beta*Y(3)+alpha*Y(3)+1.0*gamma*Kpcf*Y(3));
dY(4, 1) = gamma*Y(7)+0.001*(alpha*Y(6)-Kpcf*Y(4))+0.01*(alpha*gamma*Y(3)-4.0*beta*Kpcf*Y(4))-Kpcb*Y(4);
dY(5, 1) = 0.001*(Kpcf*Y(7)-alpha*Y(5))+Kpcb*Y(6)+0.002*(Kpcf*Y(3)-4.0*beta*Y(5))-gamma*Y(5);
dY(6, 1) = 0.001*(Kpcf*Y(4)-alpha*Y(6))+gamma*Y(5)+1.0*gamma*Kpcf*Y(3)-(4.0*beta*Kpcb*Y(6)+Kpcb*Y(6));
O_ClCa = 0.2/(1.0+exp(-(Y(25)-46.7)/7.8));
i_ClCa = g_ClCa*O_ClCa*Y(12)/(Y(12)+Km_Cl)*(Y(25)-E_Cl);
E_CaN = R*T/(2.0*F)*log(Cao/Y(12));
i_Cab = g_Cab*(Y(25)-E_CaN);
dY(9, 1) = k_plus_ltrpn*Y(12)*(LTRPN_tot-Y(9))-k_minus_ltrpn*Y(9);
dY(8, 1) = k_plus_htrpn*Y(12)*(HTRPN_tot-Y(8))-k_minus_htrpn*Y(8);
Bi = (1.0+CMDN_tot*Km_CMDN/(Km_CMDN+Y(12))^2.0)^-1.0;
J_leak = v2*(Y(11)-Y(12));
J_xfer = (Y(13)-Y(12))/tau_xfer;
J_up = v3*Y(12)^2.0/(Km_up^2.0+Y(12)^2.0);
J_trpn = k_plus_htrpn*Y(12)*(HTRPN_tot-Y(8))+k_plus_ltrpn*Y(12)*(LTRPN_tot-Y(9))-(k_minus_htrpn*Y(8)+k_minus_ltrpn*Y(9));
i_pCa = i_pCa_max*Y(12)^2.0/(Km_pCa^2.0+Y(12)^2.0);
i_NaCa = k_NaCa*1.0/(K_mNa^3.0+Nao^3.0)*1.0/(K_mCa+Cao)*1.0/(1.0+k_sat*exp((eta-1.0)*Y(25)*F/(R*T)))*(exp(eta*Y(25)*F/(R*T))*Y(39)^3.0*Cao-exp((eta-1.0)*Y(25)*F/(R*T))*Nao^3.0*Y(12));
dY(12, 1) = Bi*(J_leak+J_xfer-(J_up+J_trpn+(i_Cab+i_pCa-2.0*i_NaCa)*Acap*Cm/(2.0*Vmyo*F)));
Bss = (1.0+CMDN_tot*Km_CMDN/(Km_CMDN+Y(13))^2.0)^-1.0;
J_rel = v1*(Y(34)+Y(35))*(Y(10)-Y(13))*Y(14);
dY(13, 1) = Bss*(J_rel*VJSR/Vss-(J_xfer*Vmyo/Vss+i_CaL*Acap*Cm/(2.0*Vss*F)));
BJSR = (1.0+CSQN_tot*Km_CSQN/(Km_CSQN+Y(10))^2.0)^-1.0;
J_tr = (Y(11)-Y(10))/tau_tr;
dY(10, 1) = BJSR*(J_tr-J_rel);
dY(11, 1) = (J_up-J_leak)*Vmyo/VNSR-J_tr*VJSR/VNSR;
dY(14, 1) = -0.04*Y(14)-0.1*i_CaL/i_CaL_max*exp(-(Y(25)-5.0)^2.0/648.0);
E_Na = R*T/F*log((0.9*Nao+0.1*Ko)/(0.9*Y(39)+0.1*Y(28)));
i_Na = g_Na*Y(22)*(Y(25)-E_Na);
C_Na3 = 1.0-(Y(22)+Y(15)+Y(16)+Y(21)+Y(17)+Y(18)+Y(19)+Y(20));
alpha_Na11 = 3.802/(0.1027*exp(-(Y(25)+2.5)/17.0)+0.2*exp(-(Y(25)+2.5)/150.0));
beta_Na12 = 0.2*exp(-(Y(25)-2.5)/20.3);
alpha_Na3 = 7.0e-7*exp(-(Y(25)+7.0)/7.7);
beta_Na11 = 0.1917*exp(-(Y(25)+2.5)/20.3);
alpha_Na12 = 3.802/(0.1027*exp(-(Y(25)+2.5)/15.0)+0.23*exp(-(Y(25)+2.5)/150.0));
beta_Na3 = 0.0084+0.00002*(Y(25)+7.0);
dY(16, 1) = alpha_Na11*C_Na3+beta_Na12*Y(15)+alpha_Na3*Y(19)-(beta_Na11*Y(16)+alpha_Na12*Y(16)+beta_Na3*Y(16));
beta_Na13 = 0.22*exp(-(Y(25)-7.5)/20.3);
alpha_Na13 = 3.802/(0.1027*exp(-(Y(25)+2.5)/12.0)+0.25*exp(-(Y(25)+2.5)/150.0));
dY(15, 1) = alpha_Na12*Y(16)+beta_Na13*Y(22)+alpha_Na3*Y(21)-(beta_Na12*Y(15)+alpha_Na13*Y(15)+beta_Na3*Y(15));
alpha_Na2 = 1.0/(0.188495*exp(-(Y(25)+7.0)/16.6)+0.393956);
beta_Na2 = alpha_Na13*alpha_Na2*alpha_Na3/(beta_Na13*beta_Na3);
dY(22, 1) = alpha_Na13*Y(15)+beta_Na2*Y(21)-(beta_Na13*Y(22)+alpha_Na2*Y(22));
beta_Na4 = alpha_Na3;
alpha_Na4 = alpha_Na2/1000.0;
dY(21, 1) = alpha_Na2*Y(22)+beta_Na3*Y(15)+beta_Na4*Y(17)+alpha_Na12*Y(19)-(beta_Na2*Y(21)+alpha_Na3*Y(21)+alpha_Na4*Y(21)+beta_Na12*Y(21));
beta_Na5 = alpha_Na3/50.0;
alpha_Na5 = alpha_Na2/95000.0;
dY(17, 1) = alpha_Na4*Y(21)+beta_Na5*Y(18)-(beta_Na4*Y(17)+alpha_Na5*Y(17));
dY(18, 1) = alpha_Na5*Y(17)-beta_Na5*Y(18);
dY(19, 1) = alpha_Na11*Y(20)+beta_Na12*Y(21)+beta_Na3*Y(16)-(beta_Na11*Y(19)+alpha_Na12*Y(19)+alpha_Na3*Y(19));
dY(20, 1) = beta_Na11*Y(19)+beta_Na3*C_Na3-(alpha_Na11*Y(20)+alpha_Na3*Y(20));
E_K = R*T/F*log(Ko/Y(28));
i_Kto_f = g_Kto_f*Y(23)^3.0*Y(24)*(Y(25)-E_K);
alpha_a = 0.18064*exp(0.03577*(Y(25)+30.0));
beta_a = 0.3956*exp(-0.06237*(Y(25)+30.0));
dY(23, 1) = alpha_a*(1.0-Y(23))-beta_a*Y(23);
alpha_i_1 = 0.000152*exp(-(Y(25)+13.5)/7.0)/(0.0067083*exp(-(Y(25)+33.5)/7.0)+1.0);
beta_i_1 = 0.00095*exp((Y(25)+33.5)/7.0)/(0.051335*exp((Y(25)+33.5)/7.0)+1.0);
dY(24, 1) = alpha_i_1*(1.0-Y(24))-beta_i_1*Y(24);

if ((time >= stim_start) && (time <= stim_end) && (time-stim_start-floor((time-stim_start)/stim_period)*stim_period <= stim_duration))
   i_stim = stim_amplitude;
else
   i_stim = 0.0;
end;

i_Nab = g_Nab*(Y(25)-E_Na);
sigma = 1.0/7.0*(exp(Nao/67300.0)-1.0);
f_NaK = 1.0/(1.0+0.1245*exp(-0.1*Y(25)*F/(R*T))+0.0365*sigma*exp(-Y(25)*F/(R*T)));
i_NaK = i_NaK_max*f_NaK*1.0/(1.0+(Km_Nai/Y(39))^1.5)*Ko/(Ko+Km_Ko);
i_Kto_s = g_Kto_s*Y(37)*Y(38)*(Y(25)-E_K);
i_K1 = 0.2938*Ko/(Ko+210.0)*(Y(25)-E_K)/(1.0+exp(0.0896*(Y(25)-E_K)));
i_Ks = g_Ks*Y(36)^2.0*(Y(25)-E_K);
i_Kur = g_Kur*Y(40)*Y(41)*(Y(25)-E_K);
i_Kss = g_Kss*Y(26)*Y(27)*(Y(25)-E_K);
i_Kr = g_Kr*Y(32)*(Y(25)-R*T/F*log((0.98*Ko+0.02*Nao)/(0.98*Y(28)+0.02*Y(39))));
dY(25, 1) = -(i_CaL+i_pCa+i_NaCa+i_Cab+i_Na+i_Nab+i_NaK+i_Kto_f+i_Kto_s+i_K1+i_Ks+i_Kur+i_Kss+i_Kr+i_ClCa+i_stim);
ass = 1.0/(1.0+exp(-(Y(25)+22.5)/7.7));
tau_Kss = 39.3*exp(-0.0862*Y(25))+13.17;
dY(26, 1) = (ass-Y(26))/tau_Kss;
dY(27, 1) = 0.0;
dY(28, 1) = -(i_Kto_f+i_Kto_s+i_K1+i_Ks+i_Kss+i_Kur+i_Kr-2.0*i_NaK)*Acap*Cm/(Vmyo*F);
C_K0 = 1.0-(Y(29)+Y(30)+Y(32)+Y(31));
beta_a1 = 0.0000689*exp(-0.04178*Y(25));
alpha_a1 = 0.013733*exp(0.038198*Y(25));
dY(30, 1) = kf*Y(29)+beta_a1*Y(32)-(kb*Y(30)+alpha_a1*Y(30));
alpha_a0 = 0.022348*exp(0.01176*Y(25));
beta_a0 = 0.047002*exp(-0.0631*Y(25));
dY(29, 1) = alpha_a0*C_K0+kb*Y(30)-(beta_a0*Y(29)+kf*Y(29));
beta_i_2 = 0.006497*exp(-0.03268*(Y(25)+5.0));
alpha_i_2 = 0.090821*exp(0.023391*(Y(25)+5.0));
dY(32, 1) = alpha_a1*Y(30)+beta_i_2*Y(31)-(beta_a1*Y(32)+alpha_i_2*Y(32));
dY(31, 1) = alpha_i_2*Y(32)-beta_i_2*Y(31);
P_C1 = 1.0-(Y(33)+Y(34)+Y(35));
dY(34, 1) = k_plus_a*Y(13)^n*P_C1+k_minus_b*Y(35)+k_minus_c*Y(33)-(k_minus_a*Y(34)+k_plus_b*Y(13)^m*Y(34)+k_plus_c*Y(34));
dY(35, 1) = k_plus_b*Y(13)^m*Y(34)-k_minus_b*Y(35);
dY(33, 1) = k_plus_c*Y(34)-k_minus_c*Y(33);
alpha_n = 0.00000481333*(Y(25)+26.5)/(1.0-exp(-0.128*(Y(25)+26.5)));
beta_n = 0.0000953333*exp(-0.038*(Y(25)+26.5));
dY(36, 1) = alpha_n*(1.0-Y(36))-beta_n*Y(36);
tau_ta_s = 0.493*exp(-0.0629*Y(25))+2.058;
dY(37, 1) = (ass-Y(37))/tau_ta_s;
iss = 1.0/(1.0+exp((Y(25)+45.2)/5.7));
tau_ti_s = 270.0+1050.0/(1.0+exp((Y(25)+45.2)/5.7));
dY(38, 1) = (iss-Y(38))/tau_ti_s;
dY(39, 1) = -(i_Na+i_Nab+3.0*i_NaK+3.0*i_NaCa)*Acap*Cm/(Vmyo*F);
tau_aur = 0.493*exp(-0.0629*Y(25))+2.058;
dY(40, 1) = (ass-Y(40))/tau_aur;
tau_iur = 1200.0-170.0/(1.0+exp((Y(25)+45.2)/5.7));
dY(41, 1) = (iss-Y(41))/tau_iur;

%===============================================================================
% End of file
%===============================================================================
