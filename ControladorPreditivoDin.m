classdef ControladorPreditivoDin < matlab.System
% CONTROLADORPREDITIVODIN  GPC SISO com limites dinamicos via feedforward externo.
%
%   Entradas do bloco:  r    -> referencia r(k)
%                       y    -> saida medida da planta y(k)
%                       v_ff -> compensacao feedforward v_ff(k)
%   Saida do bloco:     u    -> acao de controle u(k) do GPC (somada externamente a v_ff)
%
% Os limites fisicos [umin, umax] definidos em cfg referem-se ao sinal total
% u + v_ff. Em cada amostra, o otimizador desloca a janela para
% [umin - v_ff, umax - v_ff] e aplica protecao de factibilidade frente aos
% limites de taxa [dumin, dumax].

    properties (Nontunable)
        % Struct de configuracao (criada por config_Caldeira_Completa.m)
        cfg = struct()
    end

    properties (Access = private)
        % Matrizes do GPC
        G
        H
        F
        Ny
        Nu

        % Matrizes e restricoes do Quadprog
        Hqp
        A_qp
        lb
        ub
        options
        umin
        umax
        dumin
        dumax
        ymin
        ymax

        % Memorias de execucao
        nMemY
        nMemDu
        yMem
        duMem
        uAnterior
        primeiraExecucao
    end

    methods
        function obj = ControladorPreditivoDin(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)

        function setupImpl(obj)
            polinomio_A = obj.cfg.A;
            polinomio_B = obj.cfg.B;
            atraso_d    = obj.cfg.d;
            N1          = obj.cfg.N1;
            N2          = obj.cfg.N2;
            Nu          = obj.cfg.Nu;
            lambda      = obj.cfg.lambda;
            horizonte_N = N2 - N1 + 1;
            grau_A      = length(polinomio_A) - 1;

            % Polinomio estendido com acao integral: A_til(z^-1) = A(z^-1)(1 - z^-1)
            A_integrado = conv(polinomio_A, [1, -1]);

            % Recursao da equacao Diofantina ate j = N2
            matriz_E = zeros(N2, N2);
            matriz_F = zeros(N2, grau_A + 1);
            matriz_E(1, 1) = 1;
            matriz_F(1, :) = -A_integrado(2:end);
            for passo = 1:(N2 - 1)
                termo_constante = matriz_F(passo, 1);
                matriz_E(passo + 1, 1:passo) = matriz_E(passo, 1:passo);
                matriz_E(passo + 1, passo + 1) = termo_constante;
                temp_F = [matriz_F(passo, :), 0] - termo_constante * A_integrado;
                matriz_F(passo + 1, :) = temp_F(2:end);
            end

            % Matrizes de predicao G, H e F
            [G_full, H_full, ~] = ref_GPC.calculaGH(matriz_E, polinomio_B, atraso_d);
            [G, H, F] = ref_GPC.recortaMatrizes(G_full, H_full, matriz_F, atraso_d, N1, N2, Nu);

            Qu  = lambda * eye(Nu);
            Hqp = 2 * (G' * G + Qu);

            dumin = obj.cfg.dumin;
            dumax = obj.cfg.dumax;
            lb    = dumin * ones(Nu, 1);
            ub    = dumax * ones(Nu, 1);

            T    = tril(ones(Nu));
            A_qp = [ T;
                    -T;
                     G;
                    -G ];

            options = optimoptions('quadprog', 'Display', 'off');

            obj.G       = G;
            obj.H       = H;
            obj.F       = F;
            obj.Hqp     = Hqp;
            obj.A_qp    = A_qp;
            obj.lb      = lb;
            obj.ub      = ub;
            obj.options = options;
            obj.umin    = obj.cfg.umin;
            obj.umax    = obj.cfg.umax;
            obj.dumin   = dumin;
            obj.dumax   = dumax;
            obj.ymin    = obj.cfg.ymin;
            obj.ymax    = obj.cfg.ymax;
            obj.Nu      = Nu;
            obj.Ny      = horizonte_N;

            obj.nMemY  = length(polinomio_A);
            obj.nMemDu = max(length(polinomio_B) - 1 + atraso_d, 1);
        end

        function resetImpl(obj)
            obj.yMem             = zeros(obj.nMemY, 1);
            obj.duMem            = zeros(obj.nMemDu, 1);
            obj.uAnterior        = obj.cfg.u0;
            obj.primeiraExecucao = true;
        end

        function u = stepImpl(obj, r, y, v_ff)
            if obj.primeiraExecucao
                obj.yMem(:) = y;
                obj.primeiraExecucao = false;
            end
            obj.yMem = desloca(obj.yMem, y);

            w = r * ones(obj.Ny, 1);

            if size(obj.H, 2) > 0
                f = obj.F * obj.yMem + obj.H * obj.duMem(1:size(obj.H, 2));
            else
                f = obj.F * obj.yMem;
            end

            fqp   = -2 * obj.G' * (w - f);
            u_ant = obj.uAnterior;

            % Janela admissivel deslocada pelo feedforward externo e protegida contra
            % infactibilidade quando um degrau de v_ff excede [dumin, dumax]
            umin_dyn = obj.umin - v_ff;
            umax_dyn = obj.umax - v_ff;
            umin_ef  = min(umin_dyn, u_ant + obj.dumax);
            umax_ef  = max(umax_dyn, u_ant + obj.dumin);

            b_qp = [ (umax_ef - u_ant) * ones(obj.Nu, 1);
                     (u_ant - umin_ef) * ones(obj.Nu, 1);
                      obj.ymax * ones(obj.Ny, 1) - f;
                     -obj.ymin * ones(obj.Ny, 1) + f ];

            [delta_u_vetor, ~, exitflag] = quadprog(obj.Hqp, fqp, obj.A_qp, b_qp, ...
                                                    [], [], obj.lb, obj.ub, [], obj.options);

            if exitflag <= 0
                A_u = obj.A_qp(1:2*obj.Nu, :);
                b_u = b_qp(1:2*obj.Nu);
                [delta_u_vetor, ~, exitflag2] = quadprog(obj.Hqp, fqp, A_u, b_u, ...
                                                         [], [], obj.lb, obj.ub, [], obj.options);
                if exitflag2 <= 0
                    du = 0;
                else
                    du = delta_u_vetor(1);
                end
            else
                du = delta_u_vetor(1);
            end

            u = obj.uAnterior + du;
            obj.duMem     = desloca(obj.duMem, du);
            obj.uAnterior = u;
        end

        function especAmostragem = getSampleTimeImpl(obj)
            especAmostragem = createSampleTime(obj, 'Type', 'Discrete', ...
                'SampleTime', obj.cfg.Ts, 'OffsetTime', 0);
        end

        function varargout = isInputDirectFeedthroughImpl(~, varargin)
            varargout = num2cell(true(1, numel(varargin)));
        end

        function tamanho = getOutputSizeImpl(~)
            tamanho = [1 1];
        end

        function tipo = getOutputDataTypeImpl(~)
            tipo = 'double';
        end

        function complexo = isOutputComplexImpl(~)
            complexo = false;
        end

        function fixo = isOutputFixedSizeImpl(~)
            fixo = true;
        end

        function [nomeR, nomeY, nomeFF] = getInputNamesImpl(~)
            nomeR  = 'r';
            nomeY  = 'y';
            nomeFF = 'v_ff';
        end

        function nomeU = getOutputNamesImpl(~)
            nomeU = 'u';
        end

        function icone = getIconImpl(~)
            icone = 'Preditivo';
        end
    end

    methods (Static, Access = protected)
        function cabecalho = getHeaderImpl
            cabecalho = matlab.system.display.Header('ControladorPreditivoDin', ...
                'Title', 'Controlador Preditivo (limites dinâmicos)', ...
                'Text', ['Recebe r(k), y(k) e v_ff(k), deslocando [umin, umax] ' ...
                         'para [umin - v_ff, umax - v_ff] com proteção de taxa.']);
        end

        function modo = getSimulateUsingImpl
            modo = 'Interpreted execution';
        end

        function mostrar = showSimulateUsingImpl
            mostrar = false;
        end
    end
end

function mem = desloca(mem, novo)
    if ~isempty(mem)
        mem = [novo; mem(1:end-1)];
    end
end
