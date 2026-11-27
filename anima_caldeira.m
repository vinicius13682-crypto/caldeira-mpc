function diag = anima_caldeira(t, r_Fv, Fv, Fg, Fr, ag, ar, r_L, L, Va, fator_vel, modo)
% ANIMA_CALDEIRA  Animacao fisica da caldeira sincronizada com os sinais simulados.
%   diag = anima_caldeira(t, r_Fv, Fv, Fg, Fr, ag, ar, r_L, L, Va, fator_vel, modo)
%   - fator_vel : velocidade de reproducao (padrao: 10 -> 260 s em ~26 s).
%   - modo      : 'gui' (padrao) abre a animacao; 'diagnostico' apenas avalia falhas.

    if nargin < 11 || isempty(fator_vel), fator_vel = 10; end
    if nargin < 12 || isempty(modo), modo = 'gui'; end

    % Normaliza para vetores coluna na resolucao original da simulacao
    t    = t(:);
    r_Fv = r_Fv(:);
    Fv   = Fv(:);
    Fg   = Fg(:);
    Fr   = Fr(:);
    ag   = ag(:);
    ar   = ar(:);
    r_L  = r_L(:);
    L    = L(:);
    Va   = Va(:);

    N = numel(t);
    dt_vec = [0; diff(t)];

    % =====================================================================
    % 1. MODELO FISICO DE SEGURANCA E DETECCAO DE FALHAS (100% PELOS SINAIS)
    % =====================================================================
    E_explosao = 0.08;
    tau_purga  = 1.5;

    R         = nan(N, 1);
    E_gas     = zeros(N, 1);
    risco_pct = zeros(N, 1);

    ativo = Fg > 0.01;
    R(ativo) = Fr(ativo) ./ Fg(ativo);

    for k = 2:N
        dt = max(0, dt_vec(k));
        if Fg(k) > 0.01
            excesso = max(0, Fg(k) - Fr(k) / 2.0);
        else
            excesso = 0;
        end
        purga = double((Fg(k) <= 0.01) || (Fr(k) >= 1.96 * Fg(k)));
        E_gas(k) = max(0, E_gas(k - 1) + dt * (excesso - purga * E_gas(k - 1) / tau_purga));
        risco_pct(k) = min(100, 100 * E_gas(k) / E_explosao);
    end

    idx_exp  = find(E_gas >= E_explosao, 1, 'first');
    idx_seco = find(L <= 0.15, 1, 'first');
    idx_over = find(L >= 3.80, 1, 'first');

    idx_falha = N;
    causa = 'OK';
    msg_falha = 'OPERACAO NORMAL';

    if ~isempty(idx_exp) && idx_exp <= idx_falha
        idx_falha = idx_exp;
        causa = 'EXPLOSAO';
        msg_falha = 'EXPLOSAO NA FORNALHA: EXCESSO DE COMBUSTIVEL';
    end
    if ~isempty(idx_seco) && idx_seco < idx_falha
        idx_falha = idx_seco;
        causa = 'TUBULAO_SECO';
        msg_falha = 'FALHA CRITICA: TUBULAO SECO (SUPERAQUECIMENTO)';
    end
    if ~isempty(idx_over) && idx_over < idx_falha
        idx_falha = idx_over;
        causa = 'TRANSBORDAMENTO';
        msg_falha = 'FALHA CRITICA: TRANSBORDAMENTO DE AGUA NA LINHA DE VAPOR';
    end

    R_ativo = R(1:idx_falha);
    R_ativo = R_ativo(~isnan(R_ativo));
    if isempty(R_ativo)
        R_min_ativo = NaN;
    else
        R_min_ativo = min(R_ativo);
    end

    diag = struct( ...
        'causa',         causa, ...
        'mensagem',      msg_falha, ...
        't_falha',       t(idx_falha), ...
        'idx_falha',     idx_falha, ...
        'risco_max_pct', max(risco_pct(1:idx_falha)), ...
        'R_min_ativo',   R_min_ativo, ...
        'L_min',         min(L(1:idx_falha)), ...
        'L_max',         max(L(1:idx_falha)) ...
    );

    if strcmpi(modo, 'diagnostico')
        return;
    end

    % =====================================================================
    % 2. CONSTRUCAO DA INTERFACE GRAFICA E ESQUEMA DA CALDEIRA
    % =====================================================================
    fig = figure( ...
        'Name', 'Simulacao Interativa — Caldeira de Geracao de Vapor', ...
        'Color', [0.973 0.980 0.988], ...
        'NumberTitle', 'off', ...
        'Position', [50 50 1380 760] ...
    );

    % Eixo esquerdo: planta fisica animada
    ax_p = axes('Parent', fig, 'Position', [0.025 0.05 0.435 0.90], ...
        'XLim', [0 10], 'YLim', [0 10.6], 'Visible', 'off');
    hold(ax_p, 'on');
    axis(ax_p, 'manual');

    text(ax_p, 5.0, 10.30, 'Caldeira Aquatubular — Estado Fisico Ao Vivo', ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.06 0.09 0.16]);

    % Camara da fornalha (inferior)
    patch(ax_p, [2.7 7.3 7.3 2.7], [1.9 1.9 4.85 4.85], [0.16 0.19 0.24], ...
        'EdgeColor', [0.25 0.30 0.38], 'LineWidth', 2.2);
    text(ax_p, 5.0, 4.62, 'FORNALHA / QUEIMADOR', ...
        'HorizontalAlignment', 'center', 'FontSize', 8.5, 'FontWeight', 'bold', 'Color', [0.78 0.83 0.90]);

    % Tubos d'agua ligando fornalha ao tubulao
    for xt = [3.5 4.5 5.5 6.5]
        rectangle(ax_p, 'Position', [xt - 0.11, 4.85, 0.22, 0.72], ...
            'FaceColor', [0.58 0.77 0.92], 'EdgeColor', [0.20 0.25 0.33], 'LineWidth', 1.2);
    end

    % Casco do tubulao superior (fundo interno y = 5.55 .. 8.35, escala L in [0, 4.0])
    y_bot = 5.55;
    y_top = 8.35;
    patch(ax_p, [2.5 7.5 7.5 2.5], [y_bot y_bot y_top y_top], [0.92 0.95 0.98], ...
        'EdgeColor', 'none');

    % Preenchimento dinamico de agua L(t)
    h_agua = patch(ax_p, [2.5 7.5 7.5 2.5], [y_bot y_bot 6.95 6.95], [0.22 0.56 0.89], ...
        'EdgeColor', [0.11 0.35 0.65], 'LineWidth', 1.5, 'FaceAlpha', 0.85);

    % Linha de referencia r_L(t) dentro do tubulao
    h_rL_line = plot(ax_p, [2.55 7.45], [6.95 6.95], '--', 'Color', [0.85 0.15 0.15], 'LineWidth', 1.6);

    % Bolhas de vapor subindo proporcionalmente a Fv
    n_bolhas = 14;
    bx_base = linspace(3.0, 7.0, n_bolhas);
    by_phase = linspace(0, 1, n_bolhas);
    h_bolhas = scatter(ax_p, bx_base, y_bot + 0.4 * ones(1, n_bolhas), 26, ...
        [0.95 0.98 1.0], 'filled', 'MarkerEdgeColor', [0.35 0.65 0.92]);

    % Contorno do tubulao por cima da agua
    rectangle(ax_p, 'Position', [2.5, y_bot, 5.0, y_top - y_bot], ...
        'Curvature', [0.06 0.12], 'EdgeColor', [0.15 0.20 0.28], 'LineWidth', 2.8);
    text(ax_p, 2.68, 8.12, 'TUBULAO (L)', 'FontSize', 8.5, 'FontWeight', 'bold', 'Color', [0.15 0.20 0.28]);
    h_txt_L = text(ax_p, 7.32, 8.12, 'L = 2.00 | r_L = 2.00', ...
        'HorizontalAlignment', 'right', 'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.06 0.22 0.45]);

    % Tubulacao e pluma de saida de vapor (topo)
    rectangle(ax_p, 'Position', [4.72, y_top, 0.56, 0.78], ...
        'FaceColor', [0.82 0.86 0.90], 'EdgeColor', [0.20 0.25 0.33], 'LineWidth', 1.8);
    h_pluma = patch(ax_p, [4.65 5.35 5.45 4.55], [9.15 9.15 9.55 9.55], ...
        [0.78 0.88 0.96], 'EdgeColor', 'none', 'FaceAlpha', 0.70);
    h_txt_Fv = text(ax_p, 5.65, 9.28, 'Vapor Fv = 0.00 (ref 0.00)', ...
        'FontSize', 9.5, 'FontWeight', 'bold', 'Color', [0.10 0.25 0.55]);

    % Linha de alimentacao de agua + bomba (esquerda do tubulao)
    plot(ax_p, [0.45 2.5], [6.60 6.60], 'Color', [0.12 0.45 0.78], 'LineWidth', 5.0);
    rectangle(ax_p, 'Position', [1.15, 6.22, 0.76, 0.76], 'Curvature', [1 1], ...
        'FaceColor', [0.90 0.95 1.0], 'EdgeColor', [0.08 0.32 0.62], 'LineWidth', 2.2);
    h_rotor1 = plot(ax_p, [1.26 1.80], [6.60 6.60], 'Color', [0.08 0.32 0.62], 'LineWidth', 2.5);
    h_rotor2 = plot(ax_p, [1.53 1.53], [6.33 6.87], 'Color', [0.08 0.32 0.62], 'LineWidth', 2.5);
    text(ax_p, 1.53, 7.22, 'Bomba Va', 'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.08 0.32 0.62]);
    h_txt_Va = text(ax_p, 1.53, 5.95, 'Va = 40%', 'HorizontalAlignment', 'center', 'FontSize', 8.8, 'FontWeight', 'bold', 'Color', [0.08 0.32 0.62]);

    % Linhas de Gas e Ar + Indicadores de Valvula (esquerda da fornalha)
    plot(ax_p, [0.45 2.7], [3.55 3.55], 'Color', [0.22 0.26 0.32], 'LineWidth', 4.0);
    plot(ax_p, [0.45 2.7], [2.45 2.45], 'Color', [0.13 0.55 0.55], 'LineWidth', 4.0);
    h_valv_ag = rectangle(ax_p, 'Position', [1.32, 3.30, 0.46, 0.50], ...
        'FaceColor', [0.25 0.28 0.34], 'EdgeColor', [0.08 0.10 0.14], 'LineWidth', 1.8);
    h_valv_ar = rectangle(ax_p, 'Position', [1.32, 2.20, 0.46, 0.50], ...
        'FaceColor', [0.10 0.55 0.55], 'EdgeColor', [0.05 0.30 0.32], 'LineWidth', 1.8);
    h_txt_gas = text(ax_p, 1.55, 3.98, 'Gas ag = 0% | Fg = 0.000', ...
        'HorizontalAlignment', 'center', 'FontSize', 8.5, 'FontWeight', 'bold', 'Color', [0.12 0.15 0.20]);
    h_txt_ar  = text(ax_p, 1.55, 1.98, 'Ar ar = 0% | Fr = 0.000', ...
        'HorizontalAlignment', 'center', 'FontSize', 8.5, 'FontWeight', 'bold', 'Color', [0.06 0.42 0.42]);

    % Chama dinamica (2 camadas dentro da fornalha)
    h_chama_ext = patch(ax_p, [4.6 5.4 5.0], [2.15 2.15 2.18], [0.20 0.72 0.98], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.88);
    h_chama_int = patch(ax_p, [4.8 5.2 5.0], [2.15 2.15 2.17], [0.88 0.96 1.00], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.95);

    % Painel HUD inferior (status ao vivo + barra de acumulo de gas)
    h_hud_box = rectangle(ax_p, 'Position', [0.25, 0.20, 9.50, 1.38], 'Curvature', [0.08 0.20], ...
        'FaceColor', [0.92 0.98 0.94], 'EdgeColor', [0.09 0.64 0.29], 'LineWidth', 2.0);
    h_hud_status = text(ax_p, 0.55, 1.24, 'STATUS: OPERACAO NORMAL', ...
        'FontSize', 10.5, 'FontWeight', 'bold', 'Color', [0.05 0.42 0.18]);
    h_hud_info = text(ax_p, 9.50, 1.24, 't = 0.0 s  |  Razao Fr/Fg = --', ...
        'HorizontalAlignment', 'right', 'FontSize', 9.5, 'FontWeight', 'bold', 'Color', [0.12 0.16 0.23]);

    text(ax_p, 0.55, 0.65, 'Acumulo de gas nao queimado:', 'FontSize', 8.8, 'FontWeight', 'bold', 'Color', [0.20 0.25 0.33]);
    rectangle(ax_p, 'Position', [4.25, 0.46, 4.10, 0.34], 'FaceColor', [0.89 0.91 0.94], ...
        'EdgeColor', [0.55 0.60 0.68], 'LineWidth', 1.2);
    h_risk_bar = rectangle(ax_p, 'Position', [4.27, 0.48, 0.01, 0.30], ...
        'FaceColor', [0.09 0.64 0.29], 'EdgeColor', 'none');
    h_txt_risk = text(ax_p, 8.55, 0.64, '0.0%', 'FontSize', 9.2, 'FontWeight', 'bold', 'Color', [0.12 0.16 0.23]);

    % Objetos pre-alocados da animacao de explosao (invisiveis ate falha)
    h_exp_shock = rectangle(ax_p, 'Position', [4.8, 3.2, 0.4, 0.4], 'Curvature', [1 1], ...
        'FaceColor', 'none', 'EdgeColor', [1.0 0.45 0.0], 'LineWidth', 3.0, 'Visible', 'off');
    h_exp_fire  = patch(ax_p, [4.5 5.5 5.0], [3.0 3.0 4.0], [1.0 0.35 0.05], ...
        'EdgeColor', [1.0 0.85 0.15], 'LineWidth', 2.0, 'FaceAlpha', 0.90, 'Visible', 'off');
    h_exp_smoke = patch(ax_p, [4.4 5.6 5.0], [4.0 4.0 5.2], [0.25 0.25 0.28], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.72, 'Visible', 'off');
    h_exp_crack = plot(ax_p, [4.2 4.8 4.6 5.4 5.1 5.8], [5.8 6.4 7.0 7.5 8.0 8.3], ...
        'Color', [0.85 0.08 0.08], 'LineWidth', 3.0, 'Visible', 'off');
    h_banner = text(ax_p, 5.0, 5.15, '', 'HorizontalAlignment', 'center', ...
        'FontSize', 10.5, 'FontWeight', 'bold', 'Color', [1 1 1], ...
        'BackgroundColor', [0.82 0.08 0.08], 'Margin', 6, 'Visible', 'off');

    % =====================================================================
    % 3. QUATRO GRAFICOS SINCRONIZADOS A DIREITA
    % =====================================================================
    t_max = max(t(end), 1);

    % Grafico 1: Vapor Fv
    ax1 = axes('Parent', fig, 'Position', [0.51 0.76 0.46 0.18]);
    hold(ax1, 'on'); grid(ax1, 'on');
    xlim(ax1, [0 t_max]);
    y1_lim = max([2.5; r_Fv(1:idx_falha); Fv(1:idx_falha)]) * 1.12;
    ylim(ax1, [0 y1_lim]);
    h_l_rFv = plot(ax1, t(1), r_Fv(1), '--r', 'LineWidth', 1.3);
    h_l_Fv  = plot(ax1, t(1), Fv(1), 'b', 'LineWidth', 1.9);
    h_m_Fv  = plot(ax1, t(1), Fv(1), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
    title(ax1, 'Malha 1 — Vazao de Vapor F_v(t)', 'FontSize', 10);
    ylabel(ax1, 'F_v');
    legend(ax1, {'r_{Fv}', 'F_v'}, 'Location', 'northeast', 'Orientation', 'horizontal');

    % Grafico 2: Razao Fr/Fg
    ax2 = axes('Parent', fig, 'Position', [0.51 0.53 0.46 0.18]);
    hold(ax2, 'on'); grid(ax2, 'on');
    xlim(ax2, [0 t_max]);
    ylim(ax2, [0.8 2.8]);
    patch(ax2, [0 t_max t_max 0], [0.8 0.8 1.80 1.80], [0.95 0.25 0.25], ...
        'FaceAlpha', 0.12, 'EdgeColor', 'none');
    plot(ax2, [0 t_max], [2.1 2.1], '--k', 'LineWidth', 1.2);
    plot(ax2, [0 t_max], [2.0 2.0], ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
    R_plot = R;
    R_plot(R_plot > 2.8) = 2.8;
    h_l_R = plot(ax2, t(1), R_plot(1), 'Color', [0.55 0.08 0.68], 'LineWidth', 1.8);
    h_m_R = plot(ax2, t(1), R_plot(1), 'o', 'Color', [0.55 0.08 0.68], ...
        'MarkerFaceColor', [0.55 0.08 0.68], 'MarkerSize', 5);
    title(ax2, 'Malha 2 — Razao Ar/Combustivel F_r / F_g (faixa vermelha < 1,80)', 'FontSize', 10);
    ylabel(ax2, 'F_r / F_g');

    % Grafico 3: Nivel L
    ax3 = axes('Parent', fig, 'Position', [0.51 0.30 0.46 0.18]);
    hold(ax3, 'on'); grid(ax3, 'on');
    xlim(ax3, [0 t_max]);
    y3_min = min([-0.1; L(1:idx_falha) - 0.2]);
    y3_max = max([4.0; L(1:idx_falha) + 0.2]);
    ylim(ax3, [y3_min y3_max]);
    plot(ax3, [0 t_max], [0.15 0.15], ':r', 'LineWidth', 1.1);
    plot(ax3, [0 t_max], [3.80 3.80], ':r', 'LineWidth', 1.1);
    h_l_rL = plot(ax3, t(1), r_L(1), '--r', 'LineWidth', 1.3);
    h_l_L  = plot(ax3, t(1), L(1), 'Color', [0.06 0.52 0.24], 'LineWidth', 1.9);
    h_m_L  = plot(ax3, t(1), L(1), 'o', 'Color', [0.06 0.52 0.24], ...
        'MarkerFaceColor', [0.06 0.52 0.24], 'MarkerSize', 5);
    title(ax3, 'Malha 3 — Nivel do Tubulao L(t)', 'FontSize', 10);
    ylabel(ax3, 'L');
    legend(ax3, {'Min/Max', '', 'r_L', 'L'}, 'Location', 'northeast', 'Orientation', 'horizontal');

    % Grafico 4: Atuadores [0, 1]
    ax4 = axes('Parent', fig, 'Position', [0.51 0.07 0.46 0.18]);
    hold(ax4, 'on'); grid(ax4, 'on');
    xlim(ax4, [0 t_max]);
    ylim(ax4, [-0.05 1.08]);
    h_l_ag = plot(ax4, t(1), ag(1), 'k', 'LineWidth', 1.4);
    h_l_ar = plot(ax4, t(1), ar(1), 'Color', [0.10 0.45 0.85], 'LineWidth', 1.4);
    h_l_Va = plot(ax4, t(1), Va(1), 'Color', [0.92 0.45 0.08], 'LineWidth', 1.5);
    title(ax4, 'Atuadores Fisicos [0, 1]', 'FontSize', 10);
    ylabel(ax4, 'Comando');
    xlabel(ax4, 'Tempo (s)');
    legend(ax4, {'a_g', 'a_r', 'V_a'}, 'Location', 'northeast', 'Orientation', 'horizontal');

    % =====================================================================
    % 4. LOOP DE REPRODUCAO EM TEMPO REAL (fator_vel)
    % =====================================================================
    t_alvo_fim = t(idx_falha);
    t_ini = tic;
    ang_rotor = 0;
    frame = 0;

    while ishandle(fig)
        t_anim = toc(t_ini) * fator_vel;
        if t_anim >= t_alvo_fim
            idx = idx_falha;
        else
            idx = find(t <= t_anim, 1, 'last');
            if isempty(idx), idx = 1; end
        end
        frame = frame + 1;

        % Agua no tubulao (escala L in [0, 4.0] mapeada em [y_bot, y_top])
        L_frac = max(0.02, min(1.0, L(idx) / 4.0));
        rL_frac = max(0.02, min(1.0, r_L(idx) / 4.0));
        y_liq = y_bot + L_frac * (y_top - y_bot);
        y_ref = y_bot + rL_frac * (y_top - y_bot);

        xw = linspace(2.5, 7.5, 25);
        onda = 0.03 * min(1, Fv(idx) / 2.0) * sin(2.4 * xw + 0.45 * frame);
        set(h_agua, 'XData', [2.5, xw, 7.5], 'YData', [y_bot, min(y_top, max(y_bot + 0.03, y_liq + onda)), y_bot]);
        set(h_rL_line, 'YData', [y_ref y_ref]);

        % Bolhas proporcionais a Fv e altura da coluna liquida
        if Fv(idx) > 0.05 && (y_liq - y_bot) > 0.15
            prog = mod(by_phase + 0.06 * frame * min(2.0, Fv(idx) / 1.5), 1);
            by_b = y_bot + 0.08 + prog * (y_liq - y_bot - 0.14);
            set(h_bolhas, 'XData', bx_base, 'YData', by_b, 'Visible', 'on');
        else
            set(h_bolhas, 'Visible', 'off');
        end

        % Pluma de vapor proporcional a Fv
        fv_norm = max(0, min(1.8, Fv(idx) / 2.0));
        set(h_pluma, ...
            'XData', [4.70 - 0.18*fv_norm, 5.30 + 0.18*fv_norm, 5.45 + 0.42*fv_norm, 4.55 - 0.42*fv_norm], ...
            'YData', [9.14, 9.14, 9.14 + 0.72*fv_norm + 0.05*sin(0.5*frame), 9.14 + 0.72*fv_norm - 0.05*sin(0.5*frame)], ...
            'FaceAlpha', min(0.85, 0.15 + 0.45 * fv_norm));
        set(h_txt_Fv, 'String', sprintf('Vapor Fv = %.2f (ref %.2f)', Fv(idx), r_Fv(idx)));
        set(h_txt_L,  'String', sprintf('L = %.2f | r_L = %.2f', L(idx), r_L(idx)));

        % Bomba girando proporcionalmente a Va
        ang_rotor = ang_rotor + 0.45 * Va(idx);
        c = 0.27 * cos(ang_rotor);
        s = 0.27 * sin(ang_rotor);
        set(h_rotor1, 'XData', [1.53 - c, 1.53 + c], 'YData', [6.60 - s, 6.60 + s]);
        set(h_rotor2, 'XData', [1.53 + s, 1.53 - s], 'YData', [6.60 - c, 6.60 + c]);
        set(h_txt_Va, 'String', sprintf('Va = %.0f%%', 100 * Va(idx)));

        % Valvulas de gas e ar + chama dinamica
        set(h_valv_ag, 'FaceColor', [0.25 + 0.65*ag(idx), 0.28, 0.34 - 0.15*ag(idx)]);
        set(h_valv_ar, 'FaceColor', [0.10, 0.45 + 0.45*ar(idx), 0.55 + 0.35*ar(idx)]);
        set(h_txt_gas, 'String', sprintf('Gas ag = %.0f%% | Fg = %.3f', 100 * ag(idx), Fg(idx)));
        set(h_txt_ar,  'String', sprintf('Ar ar = %.0f%% | Fr = %.3f', 100 * ar(idx), Fr(idx)));

        if Fg(idx) > 0.005
            fg_n = min(1.3, Fg(idx) / 0.22);
            trem = 0.06 * sin(0.9 * frame);
            w_f = 0.45 + 0.95 * fg_n;
            h_f = 0.40 + 1.55 * fg_n;

            if isnan(R(idx)) || (R(idx) >= 1.95 && R(idx) <= 2.25)
                cor_ext = [0.15 0.68 0.98];
                cor_int = [0.85 0.96 1.00];
            elseif R(idx) > 2.25
                cor_ext = [0.42 0.80 1.00];
                cor_int = [0.94 0.98 1.00];
            elseif R(idx) >= 1.70
                cor_ext = [0.98 0.56 0.10];
                cor_int = [1.00 0.90 0.35];
            else
                cor_ext = [0.92 0.18 0.08];
                cor_int = [1.00 0.65 0.12];
            end

            xf_e = [5.0 - w_f, 5.0 + w_f, 5.0 + 0.55*w_f, 5.0 + trem, 5.0 - 0.55*w_f];
            yf_e = [2.15, 2.15, 2.15 + 0.55*h_f, 2.15 + h_f, 2.15 + 0.55*h_f];
            xf_i = [5.0 - 0.55*w_f, 5.0 + 0.55*w_f, 5.0 + 0.5*trem];
            yf_i = [2.15, 2.15, 2.15 + 0.62*h_f];
            set(h_chama_ext, 'XData', xf_e, 'YData', yf_e, 'FaceColor', cor_ext, 'Visible', 'on');
            set(h_chama_int, 'XData', xf_i, 'YData', yf_i, 'FaceColor', cor_int, 'Visible', 'on');
        else
            set(h_chama_ext, 'Visible', 'off');
            set(h_chama_int, 'Visible', 'off');
        end

        % HUD de status e risco
        if isnan(R(idx))
            r_str = '--';
        else
            r_str = sprintf('%.3f', R(idx));
        end
        set(h_hud_info, 'String', sprintf('t = %.1f s  |  Razao Fr/Fg = %s', t(idx), r_str));

        rp = risco_pct(idx);
        set(h_risk_bar, 'Position', [4.27, 0.48, max(0.02, 4.06 * rp / 100), 0.30]);
        set(h_txt_risk, 'String', sprintf('%.1f%%', rp));

        if rp >= 50 || (~isnan(R(idx)) && R(idx) < 1.75) || L(idx) <= 0.30 || L(idx) >= 3.50
            set(h_hud_box, 'FaceColor', [1.00 0.95 0.88], 'EdgeColor', [0.85 0.45 0.05]);
            set(h_hud_status, 'String', 'STATUS: ALERTA DE OPERACAO', 'Color', [0.72 0.32 0.02]);
            set(h_risk_bar, 'FaceColor', [0.92 0.45 0.08]);
        else
            set(h_hud_box, 'FaceColor', [0.92 0.98 0.94], 'EdgeColor', [0.09 0.64 0.29]);
            set(h_hud_status, 'String', 'STATUS: OPERACAO NORMAL', 'Color', [0.05 0.42 0.18]);
            set(h_risk_bar, 'FaceColor', [0.09 0.64 0.29]);
        end

        % Atualiza graficos sincronizados
        set(h_l_rFv, 'XData', t(1:idx), 'YData', r_Fv(1:idx));
        set(h_l_Fv,  'XData', t(1:idx), 'YData', Fv(1:idx));
        set(h_m_Fv,  'XData', t(idx),   'YData', Fv(idx));

        set(h_l_R, 'XData', t(1:idx), 'YData', R_plot(1:idx));
        if isnan(R_plot(idx))
            set(h_m_R, 'Visible', 'off');
        else
            set(h_m_R, 'XData', t(idx), 'YData', R_plot(idx), 'Visible', 'on');
        end

        set(h_l_rL, 'XData', t(1:idx), 'YData', r_L(1:idx));
        set(h_l_L,  'XData', t(1:idx), 'YData', L(1:idx));
        set(h_m_L,  'XData', t(idx),   'YData', L(idx));

        set(h_l_ag, 'XData', t(1:idx), 'YData', ag(1:idx));
        set(h_l_ar, 'XData', t(1:idx), 'YData', ar(1:idx));
        set(h_l_Va, 'XData', t(1:idx), 'YData', Va(1:idx));

        drawnow limitrate;
        if idx >= idx_falha
            break;
        end
        pause(0.02);
    end

    % =====================================================================
    % 5. SEQUENCIA VISUAL DE EXPLOSAO / FALHA CRITICA (SE OCORREU FALHA)
    % =====================================================================
    if ishandle(fig) && ~strcmp(causa, 'OK')
        set(h_hud_box, 'FaceColor', [1.00 0.89 0.89], 'EdgeColor', [0.80 0.08 0.08]);
        set(h_hud_status, 'String', ['STATUS: ' msg_falha], 'Color', [0.75 0.05 0.05]);
        set(h_risk_bar, 'FaceColor', [0.82 0.08 0.08]);

        % Marca X vermelho nos graficos exatamente no instante da falha
        plot(ax1, t(idx_falha), Fv(idx_falha), 'rx', 'MarkerSize', 12, 'LineWidth', 2.5);
        if ~isnan(R_plot(idx_falha))
            plot(ax2, t(idx_falha), R_plot(idx_falha), 'rx', 'MarkerSize', 12, 'LineWidth', 2.5);
        end
        plot(ax3, t(idx_falha), L(idx_falha), 'rx', 'MarkerSize', 12, 'LineWidth', 2.5);

        ang_e = linspace(0, 2*pi, 24);
        rad_base = 1 + 0.32 * sin(7 * ang_e);

        set(h_banner, 'String', sprintf('%s (t = %.2f s)', msg_falha, t(idx_falha)), 'Visible', 'on');

        for f_e = 1:18
            if ~ishandle(fig), break; end
            esc = 0.35 + 0.20 * f_e;
            if strcmp(causa, 'EXPLOSAO')
                yc = 3.4;
                set(h_exp_fire, 'XData', 5.0 + esc * rad_base .* cos(ang_e), ...
                                'YData', yc + 0.85 * esc * rad_base .* sin(ang_e), ...
                                'FaceAlpha', max(0.20, 0.95 - 0.03 * f_e), 'Visible', 'on');
                set(h_exp_shock, 'Position', [5.0 - 1.25*esc, yc - 1.0*esc, 2.5*esc, 2.0*esc], ...
                                 'Visible', 'on');
                set(h_exp_smoke, 'XData', 5.0 + 1.15*esc * rad_base .* cos(ang_e), ...
                                 'YData', yc + 0.6 + 0.10*f_e + 0.75*esc * sin(ang_e), 'Visible', 'on');
                set(h_exp_crack, 'Visible', 'on');
            elseif strcmp(causa, 'TUBULAO_SECO')
                set(h_agua, 'FaceColor', [0.90 0.25 0.15]);
                set(h_exp_crack, 'Visible', 'on');
                set(h_exp_smoke, 'XData', 5.0 + 0.9*esc * rad_base .* cos(ang_e), ...
                                 'YData', 7.2 + 0.08*f_e + 0.55*esc * sin(ang_e), 'Visible', 'on');
            else
                set(h_pluma, 'FaceColor', [0.18 0.52 0.88], 'FaceAlpha', 0.90, ...
                    'XData', [4.3 5.7 6.3 3.7], 'YData', [9.14 9.14 10.1 10.1]);
                set(h_exp_shock, 'Position', [5.0 - 0.9*esc, 9.2 - 0.4*esc, 1.8*esc, 0.9*esc], ...
                    'EdgeColor', [0.15 0.55 0.90], 'Visible', 'on');
            end
            drawnow;
            pause(0.04);
        end
    end
end
