% =========================================================================
% roda_caldeira.m - Simula Caldeira_Controle.slx e executa a animacao
% =========================================================================
% Para testar erros de sintonia ou falhas fisicas (como explosao na fornalha),
% altere os blocos em Caldeira_Controle.slx (ex.: Gain_R = 1.2) ou os ganhos
% em config_caldeira.m e rode este script novamente.
% =========================================================================

clear;
clc;
close all;

cd(fileparts(mfilename('fullpath')));

%% 1. Carrega configuracao dos GPCs e executa o Simulink
config_caldeira;
out = sim('Caldeira_Controle');

t    = out.tout(:);
r_Fv = out.r_Fv(:);
Fv   = out.Fv(:);
Fg   = out.Fg(:);
Fr   = out.Fr(:);
ag   = out.ag(:);
ar   = out.ar(:);
r_L  = out.r_L(:);
L    = out.L(:);
Va   = out.Va(:);

%% 2. Diagnostico fisico real a partir dos sinais simulados
diag = anima_caldeira(t, r_Fv, Fv, Fg, Fr, ag, ar, r_L, L, Va, 10, 'diagnostico');

fprintf('=== RESUMO DA SIMULACAO (Caldeira_Controle.slx) ===\n');
fprintf('Status fisico       : %s (%s)\n', diag.causa, diag.mensagem);
fprintf('Vapor final         : Fv = %.4f (ref = %.4f)\n', Fv(end), r_Fv(end));
fprintf('Nivel final         : L  = %.4f (ref = %.4f) | L in [%.3f, %.3f]\n', ...
    L(end), r_L(end), diag.L_min, diag.L_max);
fprintf('Razao minima (Fg>0) : R_min = %.3f | Risco maximo gas = %.2f%%\n', ...
    diag.R_min_ativo, diag.risco_max_pct);
fprintf('Atuadores [0, 1]    : ag_max=%.3f | ar_max=%.3f | Va in [%.3f, %.3f]\n', ...
    max(ag), max(ar), min(Va), max(Va));

%% 3. Animacao ao vivo (fator_vel = 10 -> 260 s em ~26 s)
fator_vel = 10;
anima_caldeira(t, r_Fv, Fv, Fg, Fr, ag, ar, r_L, L, Va, fator_vel);
