% =========================================================================
% config_Caldeira_Completa.m - Configuração consolidada (Partes 1, 2 e 3)
% =========================================================================
% Define as estruturas:
%   planta  - parâmetros físicos da caldeira e cronograma do ensaio (260 s)
%   comb    - malha escrava de combustão (cascata I-P Prop A + Cross-Limiting)
%   cfg_vap - GPC mestre de vapor (para ControladorPreditivo)
%   cfg_niv - GPC de nível do tubulão (para ControladorPreditivoDin)
%
% Suporta a variável opcional 'fator_erro_valvulas' no workspace antes da
% chamada (padrão = 1.0). Os controladores sempre discretizam o modelo
% nominal; apenas planta.tau_g e planta.tau_r escalam com esse fator.
% =========================================================================

if ~exist('fator_erro_valvulas', 'var') || isempty(fator_erro_valvulas)
    fator_erro_valvulas = 1.0;
end

%% 1. Planta Física e Cronograma do Ensaio Integrado
planta = struct();
planta.tau_g_nom           = 0.10;
planta.tau_r_nom           = 0.20;
planta.fator_erro_valvulas = fator_erro_valvulas;
planta.tau_g               = planta.tau_g_nom * fator_erro_valvulas;
planta.tau_r               = planta.tau_r_nom * fator_erro_valvulas;

planta.K_g   = 2.0;
planta.K_r   = 1.0;
planta.K_cal = 10.0;
planta.tau_1 = 5.0;
planta.tau_2 = 3.0;

planta.K_a    = 5.0;
planta.tau_L  = 50.0;
planta.tau_d  = 5.0;
planta.K_d    = 0.5;
planta.tau_sw = 20.0;

% Ponto de operação inicial (Fv0 = 0, L0 = 2.0, Va0 = 0.40)
planta.L0      = 2.0;
planta.Va0     = 0.40;
planta.t_final = 260.0;

% Eventos do cenário integrado (0 a 260 s)
planta.t_Lref  = 10.0;
planta.L_ref0  = 2.0;
planta.L_ref1  = 2.5;
planta.t_Fv1   = 70.0;
planta.Fv_ref1 = 2.0;
planta.t_dg    = 130.0;
planta.amp_dg  = 0.05;
planta.t_Fv2   = 180.0;
planta.Fv_ref2 = 1.5;
planta.t_dVa   = 215.0;
planta.amp_dVa = -0.10;
planta.t5_meta = 6.48;

%% 2. Malha 2 — Escravas de Combustão (Cascata I-P no Gás + PI LGR no Ar com FF)
% Alocação via LGR sem cancelamento de polos:
%   Gás (I-P, b_g = 0): polos em -8 e -12  -> Kp_g = 0.50, Ki_g = 4.80
%   Ar  (PI,  b_r = 1): polos em -10 e -16 -> Kp_r = 4.20, Ki_r = 32.00 + FF de ag
comb = struct();
comb.R       = 2.1;
comb.polos_g = [-8.0, -12.0];
comb.polos_r = [-10.0, -16.0];
comb.bg      = 0.0;
comb.br      = 1.0;

s_sum_g  = -(comb.polos_g(1) + comb.polos_g(2));
p_prod_g = comb.polos_g(1) * comb.polos_g(2);
comb.Kp_g = (s_sum_g - 10.0) / 20.0;     % 0.50
comb.Ki_g = p_prod_g / 20.0;             % 4.80

s_sum_r  = -(comb.polos_r(1) + comb.polos_r(2));
p_prod_r = comb.polos_r(1) * comb.polos_r(2);
comb.Kp_r = (s_sum_r - 5.0) / 5.0;       % 4.20
comb.Ki_r = p_prod_r / 5.0;              % 32.00

comb.ag_min = 0.0;
comb.ag_max = 1.0;
comb.ar_min = 0.0;
comb.ar_max = 1.0;
comb.Kb_g   = 50.0;   % Back-calculation contínuo (evita chattering em ode45)
comb.Kb_r   = 50.0;

% Feedforward realizável no ar a partir de ag: 2R*(0.2s+1)/(0.1s+1) = 4R - 2R/(0.1s+1)
comb.ff_ar_D   = 4.0 * comb.R;
comb.ff_ar_num = 2.0 * comb.R;
comb.ff_ar_den = [planta.tau_g_nom, 1.0];

%% 3. Malha 1 — GPC Mestre de Vapor (Fv -> Fg_ref)
% Modelo nominal de malha fechada interna Tg(s) em série com a caldeira
s_tf = tf('s');
den_tg_nom = [1.0, 10.0 + 20.0 * comb.Kp_g, 20.0 * comb.Ki_g];
if comb.bg == 0.0
    num_tg_nom = 20.0 * comb.Ki_g;
else
    num_tg_nom = [20.0 * comb.Kp_g * comb.bg, 20.0 * comb.Ki_g];
end
G_tg_nom  = tf(num_tg_nom, den_tg_nom);
G_cal_nom = planta.K_cal / ((planta.tau_1 * s_tf + 1.0) * (planta.tau_2 * s_tf + 1.0));
G_vap_nom = G_tg_nom * G_cal_nom;

cfg_vap = struct();
cfg_vap.Ts     = 0.5;
cfg_vap.d      = 0;
cfg_vap.N1     = 1;
cfg_vap.N2     = 20;
cfg_vap.Nu     = 2;
cfg_vap.lambda = 1.2;

G_vap_z = c2d(G_vap_nom, cfg_vap.Ts, 'zoh');
[num_vap_z, den_vap_z] = tfdata(G_vap_z, 'v');
cfg_vap.A = den_vap_z / den_vap_z(1);
cfg_vap.B = num_vap_z(2:end) / den_vap_z(1);

cfg_vap.u0    = 0.0;
cfg_vap.umin  = 0.0;
cfg_vap.umax  = 0.46;   % Mantem ar abaixo da saturacao com folga para R=2.1
cfg_vap.dumin = -0.30;
cfg_vap.dumax =  0.30;
cfg_vap.ymin  = 0.0;
cfg_vap.ymax  = 200.0;

%% 4. Malha 3 — GPC de Nível do Tubulão (L -> Va_gpc) + FF Dinâmico Exato
cfg_niv = struct();
cfg_niv.Ts     = 1.0;
cfg_niv.d      = 0;
cfg_niv.N1     = 1;
cfg_niv.N2     = 30;
cfg_niv.Nu     = 1;
cfg_niv.lambda = 150.0;

G_niv_nom = planta.K_a / (planta.tau_L * s_tf + 1.0);
G_niv_z   = c2d(G_niv_nom, cfg_niv.Ts, 'zoh');
[num_niv_z, den_niv_z] = tfdata(G_niv_z, 'v');
cfg_niv.A = den_niv_z / den_niv_z(1);
cfg_niv.B = num_niv_z(2:end) / den_niv_z(1);

cfg_niv.u0    = planta.Va0;
cfg_niv.umin  = 0.0;
cfg_niv.umax  = 1.0;
cfg_niv.dumin = -0.05;
cfg_niv.dumax =  0.05;
cfg_niv.ymin  = -10.0;
cfg_niv.ymax  = 50.0;

% Feedforward dinâmico exato: G_ff(s) = -Gd(s)/Gp(s) = -0.4 + 0.5/(5s + 1)
cfg_niv.ff_D   = -0.4;
cfg_niv.ff_num = 0.5;
cfg_niv.ff_den = [5.0, 1.0];

clear s_tf s_sum_g p_prod_g s_sum_r p_prod_r den_tg_nom num_tg_nom G_tg_nom G_cal_nom G_vap_nom G_vap_z num_vap_z den_vap_z G_niv_nom G_niv_z num_niv_z den_niv_z;
