%===============================================================================
% CellML file:   D:\Desktop\Models\nygren_fiset_firek_clark_lindblad_clark_giles_1998.cellml
% CellML model:  nygren_fiset_firek_clark_lindblad_clark_giles_1998
% Date and time: 17/06/2015 at 22:21:51
%-------------------------------------------------------------------------------
% Conversion from CellML 1.0 to MATLAB (init) was done using COR (0.9.31.1409)
%    Copyright 2002-2015 Dr Alan Garny
%    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
%-------------------------------------------------------------------------------
% http://www.cellml.org/
%===============================================================================

function dY = nygren_fiset_firek_clark_lindblad_clark_giles_1998(time, Y)

%-------------------------------------------------------------------------------
% Initial conditions
%-------------------------------------------------------------------------------

% Y = [0.6465, 0.6646, 0.4284, 0.0028, 0.4369, 0.0010678, 0.949, 1.3005e-5, 0.9986, 0.9986, 1.8147, 5.3581, 130.011, 0.0048357, 0.0001, 0.0275, 0.0133, 0.1961, 0.7094, 7.2495e-5, 6.729e-5, 129.435, 8.5547, -74.2525, 0.8814, 0.8742, 0.0032017, 0.00015949, 0.9912];

% YNames = {'Ca_rel', 'Ca_up', 'F1', 'F2', 'O_Calse', 'r', 's', 'd_L', 'f_L_1', 'f_L_2', 'Ca_c', 'K_c', 'Na_c', 'n', 'p_a', 'O_C', 'O_TC', 'O_TMgC', 'O_TMgMg', 'Ca_d', 'Ca_i', 'K_i', 'Na_i', 'V', 'h1', 'h2', 'm', 'r_sus', 's_sus'};
% YUnits = {'millimolar', 'millimolar', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'millimolar', 'millimolar', 'millimolar', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'millimolar', 'millimolar', 'millimolar', 'millimolar', 'millivolt', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless'};
% YComponents = {'Ca_handling_by_the_SR', 'Ca_handling_by_the_SR', 'Ca_handling_by_the_SR', 'Ca_handling_by_the_SR', 'Ca_handling_by_the_SR', 'Ca_independent_transient_outward_K_current_r_gate', 'Ca_independent_transient_outward_K_current_s_gate', 'L_type_Ca_channel_d_L_gate', 'L_type_Ca_channel_f_L1_gate', 'L_type_Ca_channel_f_L2_gate', 'cleft_space_ion_concentrations', 'cleft_space_ion_concentrations', 'cleft_space_ion_concentrations', 'delayed_rectifier_K_currents_n_gate', 'delayed_rectifier_K_currents_pa_gate', 'intracellular_Ca_buffering', 'intracellular_Ca_buffering', 'intracellular_Ca_buffering', 'intracellular_Ca_buffering', 'intracellular_ion_concentrations', 'intracellular_ion_concentrations', 'intracellular_ion_concentrations', 'intracellular_ion_concentrations', 'membrane', 'sodium_current_h1_gate', 'sodium_current_h2_gate', 'sodium_current_m_gate', 'sustained_outward_K_current_r_sus_gate', 'sustained_outward_K_current_s_sus_gate'};

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: Ca_rel (millimolar) (in Ca_handling_by_the_SR)
% 2: Ca_up (millimolar) (in Ca_handling_by_the_SR)
% 3: F1 (dimensionless) (in Ca_handling_by_the_SR)
% 4: F2 (dimensionless) (in Ca_handling_by_the_SR)
% 5: O_Calse (dimensionless) (in Ca_handling_by_the_SR)
% 6: r (dimensionless) (in Ca_independent_transient_outward_K_current_r_gate)
% 7: s (dimensionless) (in Ca_independent_transient_outward_K_current_s_gate)
% 8: d_L (dimensionless) (in L_type_Ca_channel_d_L_gate)
% 9: f_L_1 (dimensionless) (in L_type_Ca_channel_f_L1_gate)
% 10: f_L_2 (dimensionless) (in L_type_Ca_channel_f_L2_gate)
% 11: Ca_c (millimolar) (in cleft_space_ion_concentrations)
% 12: K_c (millimolar) (in cleft_space_ion_concentrations)
% 13: Na_c (millimolar) (in cleft_space_ion_concentrations)
% 14: n (dimensionless) (in delayed_rectifier_K_currents_n_gate)
% 15: p_a (dimensionless) (in delayed_rectifier_K_currents_pa_gate)
% 16: O_C (dimensionless) (in intracellular_Ca_buffering)
% 17: O_TC (dimensionless) (in intracellular_Ca_buffering)
% 18: O_TMgC (dimensionless) (in intracellular_Ca_buffering)
% 19: O_TMgMg (dimensionless) (in intracellular_Ca_buffering)
% 20: Ca_d (millimolar) (in intracellular_ion_concentrations)
% 21: Ca_i (millimolar) (in intracellular_ion_concentrations)
% 22: K_i (millimolar) (in intracellular_ion_concentrations)
% 23: Na_i (millimolar) (in intracellular_ion_concentrations)
% 24: V (millivolt) (in membrane)
% 25: h1 (dimensionless) (in sodium_current_h1_gate)
% 26: h2 (dimensionless) (in sodium_current_h2_gate)
% 27: m (dimensionless) (in sodium_current_m_gate)
% 28: r_sus (dimensionless) (in sustained_outward_K_current_r_sus_gate)
% 29: s_sus (dimensionless) (in sustained_outward_K_current_s_sus_gate)

%-------------------------------------------------------------------------------
% Constants
%-------------------------------------------------------------------------------

I_up_max = 2800.0;   % picoA (in Ca_handling_by_the_SR)
Vol_rel = 4.41e-5;   % nanolitre (in Ca_handling_by_the_SR)
Vol_up = 0.0003969;   % nanolitre (in Ca_handling_by_the_SR)
alpha_rel = 200000.0;   % picoA_per_millimolar (in Ca_handling_by_the_SR)
k_cyca = 0.0003;   % millimolar (in Ca_handling_by_the_SR)
k_rel_d = 0.003;   % millimolar (in Ca_handling_by_the_SR)
k_rel_i = 0.0003;   % millimolar (in Ca_handling_by_the_SR)
k_srca = 0.5;   % millimolar (in Ca_handling_by_the_SR)
k_xcs = 0.4;   % dimensionless (in Ca_handling_by_the_SR)
r_recov = 0.815;   % per_second (in Ca_handling_by_the_SR)
tau_tr = 0.01;   % second (in Ca_handling_by_the_SR)
g_t = 7.5;   % nanoS (in Ca_independent_transient_outward_K_current)
E_Ca_app = 60.0;   % millivolt (in L_type_Ca_channel)
g_Ca_L = 6.75;   % nanoS (in L_type_Ca_channel)
k_Ca = 0.025;   % millimolar (in L_type_Ca_channel)
d_NaCa = 0.0003;   % per_millimolar_4 (in Na_Ca_ion_exchanger_current)
gamma = 0.45;   % dimensionless (in Na_Ca_ion_exchanger_current)
k_NaCa = 0.0374842;   % picoA_per_millimolar_4 (in Na_Ca_ion_exchanger_current)
g_B_Ca = 0.078681;   % nanoS (in background_currents)
g_B_Na = 0.060599;   % nanoS (in background_currents)
Ca_b = 1.8;   % millimolar (in cleft_space_ion_concentrations)
K_b = 5.4;   % millimolar (in cleft_space_ion_concentrations)
Na_b = 130.0;   % millimolar (in cleft_space_ion_concentrations)
tau_Ca = 24.7;   % second (in cleft_space_ion_concentrations)
tau_K = 10.0;   % second (in cleft_space_ion_concentrations)
tau_Na = 14.3;   % second (in cleft_space_ion_concentrations)
g_Kr = 0.5;   % nanoS (in delayed_rectifier_K_currents)
g_Ks = 1.0;   % nanoS (in delayed_rectifier_K_currents)
Mg_i = 2.5;   % millimolar (in intracellular_Ca_buffering)
Vol_i = 0.005884;   % nanolitre (in intracellular_ion_concentrations)
phi_Na_en = -1.68;   % picoA (in intracellular_ion_concentrations)
tau_di = 0.01;   % second (in intracellular_ion_concentrations)
g_K1 = 3.0;   % nanoS (in inward_rectifier)
Cm = 0.05;   % nanoF (in membrane)
F = 96487.0;   % coulomb_per_mole (in membrane)
R = 8314.0;   % millijoule_per_mole_kelvin (in membrane)
T = 306.15;   % kelvin (in membrane)
stim_amplitude = -280.0;   % picoA (in membrane)
stim_duration = 0.006;   % second (in membrane)
stim_end = 100000000.0;   % second (in membrane)
stim_period = 1.0;   % second (in membrane)
stim_start = 0.1;   % second (in membrane)
i_CaP_max = 4.0;   % picoA (in sarcolemmal_calcium_pump_current)
k_CaP = 0.0002;   % millimolar (in sarcolemmal_calcium_pump_current)
P_Na = 0.0016;   % nanolitre_per_second (in sodium_current)
i_NaK_max = 70.8253;   % picoA (in sodium_potassium_pump)
k_NaK_K = 1.0;   % millimolar (in sodium_potassium_pump)
k_NaK_Na = 11.0;   % millimolar (in sodium_potassium_pump)
g_sus = 2.75;   % nanoS (in sustained_outward_K_current)

%-------------------------------------------------------------------------------
% Computed variables
%-------------------------------------------------------------------------------

% i_rel (picoA) (in Ca_handling_by_the_SR)
% i_tr (picoA) (in Ca_handling_by_the_SR)
% i_up (picoA) (in Ca_handling_by_the_SR)
% r_act (per_second) (in Ca_handling_by_the_SR)
% r_inact (per_second) (in Ca_handling_by_the_SR)
% r_infinity (dimensionless) (in Ca_independent_transient_outward_K_current_r_gate)
% tau_r (second) (in Ca_independent_transient_outward_K_current_r_gate)
% s_infinity (dimensionless) (in Ca_independent_transient_outward_K_current_s_gate)
% tau_s (second) (in Ca_independent_transient_outward_K_current_s_gate)
% E_K (millivolt) (in Ca_independent_transient_outward_K_current)
% i_t (picoA) (in Ca_independent_transient_outward_K_current)
% d_L_infinity (dimensionless) (in L_type_Ca_channel_d_L_gate)
% tau_d_L (second) (in L_type_Ca_channel_d_L_gate)
% f_L_infinity (dimensionless) (in L_type_Ca_channel_f_L1_gate)
% tau_f_L1 (second) (in L_type_Ca_channel_f_L1_gate)
% tau_f_L2 (second) (in L_type_Ca_channel_f_L2_gate)
% f_Ca (dimensionless) (in L_type_Ca_channel)
% i_Ca_L (picoA) (in L_type_Ca_channel)
% i_NaCa (picoA) (in Na_Ca_ion_exchanger_current)
% E_Ca (millivolt) (in background_currents)
% i_B_Ca (picoA) (in background_currents)
% i_B_Na (picoA) (in background_currents)
% Vol_c (nanolitre) (in cleft_space_ion_concentrations)
% n_infinity (dimensionless) (in delayed_rectifier_K_currents_n_gate)
% tau_n (second) (in delayed_rectifier_K_currents_n_gate)
% p_a_infinity (dimensionless) (in delayed_rectifier_K_currents_pa_gate)
% tau_p_a (second) (in delayed_rectifier_K_currents_pa_gate)
% p_i (dimensionless) (in delayed_rectifier_K_currents_pi_gate)
% i_Kr (picoA) (in delayed_rectifier_K_currents)
% i_Ks (picoA) (in delayed_rectifier_K_currents)
% dOCdt (per_second) (in intracellular_Ca_buffering)
% dOTCdt (per_second) (in intracellular_Ca_buffering)
% dOTMgCdt (per_second) (in intracellular_Ca_buffering)
% Vol_d (nanolitre) (in intracellular_ion_concentrations)
% i_di (picoA) (in intracellular_ion_concentrations)
% i_K1 (picoA) (in inward_rectifier)
% i_Stim (picoA) (in membrane)
% i_CaP (picoA) (in sarcolemmal_calcium_pump_current)
% h_infinity (dimensionless) (in sodium_current_h1_gate)
% tau_h1 (second) (in sodium_current_h1_gate)
% tau_h2 (second) (in sodium_current_h2_gate)
% m_infinity (dimensionless) (in sodium_current_m_gate)
% tau_m (second) (in sodium_current_m_gate)
% E_Na (millivolt) (in sodium_current)
% i_Na (picoA) (in sodium_current)
% i_NaK (picoA) (in sodium_potassium_pump)
% r_sus_infinity (dimensionless) (in sustained_outward_K_current_r_sus_gate)
% tau_r_sus (second) (in sustained_outward_K_current_r_sus_gate)
% s_sus_infinity (dimensionless) (in sustained_outward_K_current_s_sus_gate)
% tau_s_sus (second) (in sustained_outward_K_current_s_sus_gate)
% i_sus (picoA) (in sustained_outward_K_current)

%-------------------------------------------------------------------------------
% Computation
%-------------------------------------------------------------------------------

% time (second)

i_up = I_up_max*(Y(21)/k_cyca-k_xcs^2.0*Y(2)/k_srca)/((Y(21)+k_cyca)/k_cyca+k_xcs*(Y(2)+k_srca)/k_srca);
i_tr = (Y(2)-Y(1))*2.0*F*Vol_rel/tau_tr;
i_rel = alpha_rel*(Y(4)/(Y(4)+0.25))^2.0*(Y(1)-Y(21));
dY(5, 1) = 480.0*Y(1)*(1.0-Y(5))-400.0*Y(5);
dY(1, 1) = (i_tr-i_rel)/(2.0*Vol_rel*F)-31.0*dY(5, 1);
dY(2, 1) = (i_up-i_tr)/(2.0*Vol_up*F);
r_act = 203.8*((Y(21)/(Y(21)+k_rel_i))^4.0+(Y(20)/(Y(20)+k_rel_d))^4.0);
dY(3, 1) = r_recov*(1.0-Y(3)-Y(4))-r_act*Y(3);
r_inact = 33.96+339.6*(Y(21)/(Y(21)+k_rel_i))^4.0;
dY(4, 1) = r_act*Y(3)-r_inact*Y(4);
E_K = R*T/F*log(Y(12)/Y(22));
i_t = g_t*Y(6)*Y(7)*(Y(24)-E_K);
r_infinity = 1.0/(1.0+exp((Y(24)-1.0)/-11.0));
tau_r = 0.0035*exp(-(Y(24)/30.0)^2.0)+0.0015;
dY(6, 1) = (r_infinity-Y(6))/tau_r;
s_infinity = 1.0/(1.0+exp((Y(24)+40.5)/11.5));
tau_s = 0.4812*exp(-((Y(24)+52.45)/14.97)^2.0)+0.01414;
dY(7, 1) = (s_infinity-Y(7))/tau_s;
f_Ca = Y(20)/(Y(20)+k_Ca);
i_Ca_L = g_Ca_L*Y(8)*(f_Ca*Y(9)+(1.0-f_Ca)*Y(10))*(Y(24)-E_Ca_app);
d_L_infinity = 1.0/(1.0+exp((Y(24)+9.0)/-5.8));
tau_d_L = 0.0027*exp(-((Y(24)+35.0)/30.0)^2.0)+0.002;
dY(8, 1) = (d_L_infinity-Y(8))/tau_d_L;
f_L_infinity = 1.0/(1.0+exp((Y(24)+27.4)/7.1));
tau_f_L1 = 0.161*exp(-((Y(24)+40.0)/14.4)^2.0)+0.01;
dY(9, 1) = (f_L_infinity-Y(9))/tau_f_L1;
tau_f_L2 = 1.3323*exp(-((Y(24)+40.0)/14.2)^2.0)+0.0626;
dY(10, 1) = (f_L_infinity-Y(10))/tau_f_L2;
i_NaCa = k_NaCa*(Y(23)^3.0*Y(11)*exp(gamma*F*Y(24)/(R*T))-Y(13)^3.0*Y(21)*exp((gamma-1.0)*Y(24)*F/(R*T)))/(1.0+d_NaCa*(Y(13)^3.0*Y(21)+Y(23)^3.0*Y(11)));
E_Na = R*T/F*log(Y(13)/Y(23));
i_B_Na = g_B_Na*(Y(24)-E_Na);
E_Ca = R*T/(2.0*F)*log(Y(11)/Y(21));
i_B_Ca = g_B_Ca*(Y(24)-E_Ca);
Vol_c = 0.136*Vol_i;
i_Na = P_Na*Y(27)^3.0*(0.9*Y(25)+0.1*Y(26))*Y(13)*Y(24)*F^2.0/(R*T)*(exp((Y(24)-E_Na)*F/(R*T))-1.0)/(exp(Y(24)*F/(R*T))-1.0);
i_NaK = i_NaK_max*Y(12)/(Y(12)+k_NaK_K)*Y(23)^1.5/(Y(23)^1.5+k_NaK_Na^1.5)*(Y(24)+150.0)/(Y(24)+200.0);
dY(13, 1) = (Na_b-Y(13))/tau_Na+(i_Na+i_B_Na+3.0*i_NaK+3.0*i_NaCa+phi_Na_en)/(Vol_c*F);
i_sus = g_sus*Y(28)*Y(29)*(Y(24)-E_K);
i_K1 = g_K1*(Y(12)/1.0)^0.4457*(Y(24)-E_K)/(1.0+exp(1.5*(Y(24)-E_K+3.6)*F/(R*T)));
p_i = 1.0/(1.0+exp((Y(24)+55.0)/24.0));
i_Kr = g_Kr*Y(15)*p_i*(Y(24)-E_K);
i_Ks = g_Ks*Y(14)*(Y(24)-E_K);
dY(12, 1) = (K_b-Y(12))/tau_K+(i_t+i_sus+i_K1+i_Kr+i_Ks-2.0*i_NaK)/(Vol_c*F);
i_CaP = i_CaP_max*Y(21)/(Y(21)+k_CaP);
dY(11, 1) = (Ca_b-Y(11))/tau_Ca+(i_Ca_L+i_B_Ca+i_CaP-2.0*i_NaCa)/(2.0*Vol_c*F);
n_infinity = 1.0/(1.0+exp((Y(24)-19.9)/-12.7));
tau_n = 0.7+0.4*exp(-((Y(24)-20.0)/20.0)^2.0);
dY(14, 1) = (n_infinity-Y(14))/tau_n;
p_a_infinity = 1.0/(1.0+exp((Y(24)+15.0)/-6.0));
tau_p_a = 0.03118+0.21718*exp(-((Y(24)+20.1376)/22.1996)^2.0);
dY(15, 1) = (p_a_infinity-Y(15))/tau_p_a;
dY(16, 1) = 200000.0*Y(21)*(1.0-Y(16))-476.0*Y(16);
dOCdt = dY(16, 1);
dY(17, 1) = 78400.0*Y(21)*(1.0-Y(17))-392.0*Y(17);
dOTCdt = dY(17, 1);
dY(18, 1) = 200000.0*Y(21)*(1.0-Y(18)-Y(19))-6.6*Y(18);
dOTMgCdt = dY(18, 1);
dY(19, 1) = 2000.0*Mg_i*(1.0-Y(18)-Y(19))-666.0*Y(19);
Vol_d = 0.02*Vol_i;
dY(23, 1) = -(i_Na+i_B_Na+3.0*i_NaK+3.0*i_NaCa+phi_Na_en)/(Vol_i*F);
dY(22, 1) = -(i_t+i_sus+i_K1+i_Kr+i_Ks-2.0*i_NaK)/(Vol_i*F);
i_di = (Y(20)-Y(21))*2.0*F*Vol_d/tau_di;
dY(21, 1) = -(-i_di+i_B_Ca+i_CaP-2.0*i_NaCa+i_up-i_rel)/(2.0*Vol_i*F)-(0.08*dOTCdt+0.16*dOTMgCdt+0.045*dOCdt);
dY(20, 1) = -(i_Ca_L+i_di)/(2.0*Vol_d*F);

if ((time >= stim_start) && (time <= stim_end) && (time-stim_start-floor((time-stim_start)/stim_period)*stim_period <= stim_duration))
   i_Stim = stim_amplitude;
else
   i_Stim = 0.0;
end;

dY(24, 1) = -1.0/Cm*(i_Stim+i_Na+i_Ca_L+i_t+i_sus+i_K1+i_Kr+i_Ks+i_B_Na+i_B_Ca+i_NaK+i_CaP+i_NaCa);
h_infinity = 1.0/(1.0+exp((Y(24)+63.6)/5.3));
tau_h1 = 0.03/(1.0+exp((Y(24)+35.1)/3.2))+0.0003;
dY(25, 1) = (h_infinity-Y(25))/tau_h1;
tau_h2 = 0.12/(1.0+exp((Y(24)+35.1)/3.2))+0.003;
dY(26, 1) = (h_infinity-Y(26))/tau_h2;
m_infinity = 1.0/(1.0+exp((Y(24)+27.12)/-8.21));
tau_m = 4.2e-5*exp(-((Y(24)+25.57)/28.8)^2.0)+2.4e-5;
dY(27, 1) = (m_infinity-Y(27))/tau_m;
r_sus_infinity = 1.0/(1.0+exp((Y(24)+4.3)/-8.0));
tau_r_sus = 0.009/(1.0+exp((Y(24)+5.0)/12.0))+0.0005;
dY(28, 1) = (r_sus_infinity-Y(28))/tau_r_sus;
s_sus_infinity = 0.4/(1.0+exp((Y(24)+20.0)/10.0))+0.6;
tau_s_sus = 0.047/(1.0+exp((Y(24)+60.0)/10.0))+0.3;
dY(29, 1) = (s_sus_infinity-Y(29))/tau_s_sus;

%===============================================================================
% End of file
%===============================================================================
