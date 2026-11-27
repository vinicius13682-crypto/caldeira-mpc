classdef ControladorPreditivo < matlab.System
% CONTROLADORPREDITIVO  Modelo (template) de controlador preditivo SISO
% para o bloco MATLAB System do Simulink.
%
%   Entradas do bloco:  r -> referencia r(k)
%                       y -> saida medida da planta y(k)
%   Saida do bloco:     u -> acao de controle u(k)
%
%   >>> Os trechos a implementar estao marcados com  "ALUNO:"  <<<
%
% O QUE JA ESTA PRONTO
%   - Integracao com o Simulink: o bloco executa a cada cfg.Ts, com
%     entradas r e y e saida u.
%   - Leitura da configuracao: struct cfg criada por config_Preditivo.m,
%     acessada aqui como obj.cfg.<campo>.
%   - Memorias de valores passados, atualizadas a cada amostra:
%       obj.yMem  = [y(k); y(k-1); ...]         (tamanho obj.nMemY)
%       obj.duMem = [du(k-1); du(k-2); ...]     (tamanho obj.nMemDu)
%       obj.uAnterior = u(k-1)
%   - Lei de controle incremental: u(k) = u(k-1) + du(k).
%
% O QUE O ALUNO IMPLEMENTA
%   1) setupImpl : calculos feitos UMA vez, antes da simulacao (matrizes
%                  de predicao, ganhos...) e o tamanho das memorias que o
%                  seu algoritmo precisa.
%   2) stepImpl  : resposta livre f e calculo do incremento du(k).
%   3) (se precisar) novas propriedades no bloco "Variaveis do algoritmo".
%
% Enquanto nada for implementado o bloco devolve du = 0 (u = cfg.u0) e o
% modelo ja simula sem erros.
%
% CICLO DE VIDA DE UM SYSTEM OBJECT NO SIMULINK
%   setupImpl  -> 1 vez, antes da simulacao
%   resetImpl  -> 1 vez, antes da simulacao: inicializa as memorias
%   stepImpl   -> a cada Ts: le r(k) e y(k), devolve u(k)
%   Os demais metodos apenas informam ao Simulink o tempo de amostragem,
%   o tamanho/tipo das portas, os nomes das portas e o icone do bloco
%   (nao precisam ser alterados).
%
% DICAS
%   - O codigo roda em modo interpretado: pode usar qualquer funcao do
%     MATLAB (inclusive quadprog) e colocar breakpoints no stepImpl.
%   - Para depurar fora do Simulink:
%         config_Preditivo;  ctrl = ControladorPreditivo('cfg', cfg);
%         u = ctrl(r, y);    % uma chamada = uma amostra
%   - Se renomear a classe/arquivo, altere tambem o campo "System object
%     name" do bloco no Simulink.

    % =====================================================================
    % Parametro do bloco (nao alterar)
    % =====================================================================
    % Nontunable: nao pode mudar durante a simulacao. Na janela do bloco
    % digita-se o NOME da variavel do workspace que contem a struct (cfg).
    properties (Nontunable)
        % Struct de configuracao (criada por config_Preditivo.m)
        cfg = struct()
    end

    % =====================================================================
    % Variaveis do algoritmo
    % =====================================================================
    properties (Access = private)
        % ----- ALUNO: declare aqui o que for calculado no setupImpl -----
        % Matrizes do GPC
        G
        H
        F
        Ny
        Nu
        
        % Matrizes estáticas da Programação Quadrática (QP)
        Hqp
        A_qp
        lb
        ub
        options
        umin
        umax
        ymin
        ymax
        % ------------------------------------------------------------------

        % --- Memorias (ja gerenciadas pelo template) ---
        nMemY             % numero de saidas passadas guardadas em yMem
        nMemDu            % numero de incrementos passados guardados em duMem
        yMem              % [y(k); y(k-1); ...; y(k-nMemY+1)]
        duMem             % [du(k-1); du(k-2); ...; du(k-nMemDu)]
        uAnterior         % u(k-1)
        primeiraExecucao  % true ate a primeira execucao do stepImpl
    end

    % =====================================================================
    % Construtor (nao alterar)
    % =====================================================================
    methods
        function obj = ControladorPreditivo(varargin)
            % Permite criar o objeto tambem no MATLAB, fora do Simulink:
            %   ctrl = ControladorPreditivo('cfg', cfg);   u = ctrl(r, y);
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)

        % =================================================================
        % Execucao do algoritmo
        % =================================================================
                function setupImpl(obj)
            % Executado UMA vez, antes da primeira amostra.
            % Calculos que nao dependem das medicoes (feitos "offline").
            c = obj.cfg; %#ok<NASGU>

            %% Minha implementação 
            polinomio_A = obj.cfg.A;
            polinomio_B = obj.cfg.B;
            atraso_d    = obj.cfg.d;
            N1          = obj.cfg.N1;
            N2          = obj.cfg.N2;
            Nu          = obj.cfg.Nu;
            lambda      = obj.cfg.lambda;
            horizonte_N = N2 - N1 + 1;
            grau_A      = length(polinomio_A) - 1;

            % --- Polinômio estendido com integrador ---
            A_integrado = conv(polinomio_A, [1, -1]);
            
            % --- Diofantina ---
            matriz_E = zeros(N2, horizonte_N);
            matriz_F = zeros(N2, grau_A + 1);
            matriz_E(1, 1) = 1;
            matriz_F(1, :) = -A_integrado(2:end);
            for passo = 1:(horizonte_N - 1)
                termo_constante = matriz_F(passo, 1);
                matriz_E(passo + 1, 1:passo) = matriz_E(passo, 1:passo);
                matriz_E(passo + 1, passo + 1) = termo_constante;
                temp_F = [matriz_F(passo, :), 0] - termo_constante * A_integrado;
                matriz_F(passo + 1, :) = temp_F(2:end);
            end
            
            % --- Matrizes G, H, F ---
            [G_full, H_full, M] = ref_GPC.calculaGH(matriz_E, polinomio_B, atraso_d);
            [G, H, F] = ref_GPC.recortaMatrizes(G_full, H_full, matriz_F, atraso_d, N1, N2, Nu);
            
            % --- Estruturas do Quadprog ---
            Qu = lambda * eye(Nu);
            Hqp = 2 * (G' * G + Qu);
            
            % Restricoes lidas do config_Preditivo.m (struct cfg)
            dumin = obj.cfg.dumin; 
            dumax = obj.cfg.dumax; 
            lb = dumin * ones(Nu, 1);
            ub = dumax * ones(Nu, 1);
            
            umin = obj.cfg.umin;
            umax = obj.cfg.umax;
            ymin = obj.cfg.ymin;
            ymax = obj.cfg.ymax;
            
            T = tril(ones(Nu));
            A_qp = [  T;     % Limite superior de u
                     -T;     % Limite inferior de u
                      G;     % Limite superior da saída y
                     -G ];   % Limite inferior da saída y
                 
            options = optimoptions('quadprog', 'Display', 'off');

            % --- Salvando na memória do objeto para usar no stepImpl ---
            obj.G       = G;
            obj.H       = H;
            obj.F       = F;
            obj.Hqp     = Hqp;
            obj.A_qp    = A_qp;
            obj.lb      = lb;
            obj.ub      = ub;
            obj.options = options;
            obj.umin    = umin;
            obj.umax    = umax;
            obj.ymin    = ymin;
            obj.ymax    = ymax;
            obj.Nu      = Nu;
            obj.Ny      = horizonte_N;

            % --- Memórias gerenciadas pelo template ---
            obj.nMemY  = length(polinomio_A);
            obj.nMemDu = max(length(polinomio_B) - 1 + atraso_d, 1);
        end

        function resetImpl(obj)
            % Valores iniciais das memorias (nao precisa alterar).
            % Executado antes da primeira amostra e a cada reinicio.
            obj.yMem             = zeros(obj.nMemY, 1);
            obj.duMem            = zeros(obj.nMemDu, 1);
            obj.uAnterior        = obj.cfg.u0;
            obj.primeiraExecucao = true;
        end

        function u = stepImpl(obj, r, y)
            % Lei de controle, executada uma vez por periodo de amostragem.
            %   r -> referencia r(k)       (entrada 1 do bloco)
            %   y -> saida medida y(k)     (entrada 2 do bloco)
            %   u -> acao de controle u(k) (saida do bloco)
            c = obj.cfg;

            % --- 1) Memoria de y (pronto no template) -------------------
            if obj.primeiraExecucao
                obj.yMem(:) = y;
                obj.primeiraExecucao = false;
            end
            obj.yMem = desloca(obj.yMem, y);   % [y(k); y(k-1); ...]

            % --- 2) Trajetoria de referencia (pronto no template) -------
            Ny = c.N2 - c.N1 + 1;
            w  = r * ones(Ny, 1); 

            % --- 3) ALUNO: resposta livre f (Ny x 1) --------------------
            if size(obj.H, 2) > 0
                f = obj.F * obj.yMem + obj.H * obj.duMem(1:size(obj.H, 2));
            else
                f = obj.F * obj.yMem;
            end

            % --- 4) ALUNO: otimizacao -> du(k) ---------------------------
            fqp = -2 * obj.G' * (w - f);

            u_ant = obj.uAnterior;
            b_qp = [ (obj.umax - u_ant) * ones(obj.Nu, 1);
                     (u_ant - obj.umin) * ones(obj.Nu, 1);
                      obj.ymax * ones(obj.Ny, 1) - f;
                     -obj.ymin * ones(obj.Ny, 1) + f ];

            [delta_u_vetor, ~, exitflag] = quadprog(obj.Hqp, fqp, obj.A_qp, b_qp, ...
                                                    [], [], obj.lb, obj.ub, [], obj.options);

            if exitflag <= 0
                % Fallback: se limites de saida y causarem infactibilidade devido a perturbacao,
                % resolve respeitando os limites fisicos do atuador (u e du) para trazer a planta de volta
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

            % --- 5) Lei de controle incremental (pronto no template) ----
            u = obj.uAnterior + du;

            % --- 6) Atualiza memorias para a proxima amostra ------------
            obj.duMem     = desloca(obj.duMem, du);
            obj.uAnterior = u;
        end

        % =================================================================
        % Informacoes para o Simulink (nao precisa alterar)
        % =================================================================
        function especAmostragem = getSampleTimeImpl(obj)
            % Bloco discreto, executado a cada cfg.Ts a partir de t = 0.
            especAmostragem = createSampleTime(obj, 'Type', 'Discrete', ...
                'SampleTime', obj.cfg.Ts, 'OffsetTime', 0);
        end

        function varargout = isInputDirectFeedthroughImpl(~, varargin)
            % As entradas r(k) e y(k) sao usadas no calculo de u(k) na
            % mesma execucao (alimentacao direta).
            varargout = num2cell(true(1, numel(varargin)));
        end

        function tamanho = getOutputSizeImpl(~)
            tamanho = [1 1];          % u e escalar
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

        function [nomeR, nomeY] = getInputNamesImpl(~)
            % Rotulos das portas de entrada (na ordem das entradas do stepImpl)
            nomeR = 'r';
            nomeY = 'y';
        end

        function nomeU = getOutputNamesImpl(~)
            % Rotulo da porta de saida
            nomeU = 'u';
        end

        function icone = getIconImpl(~)
            % Texto no centro do icone do bloco (pode trocar pelo nome do
            % seu algoritmo, ex.: 'DMC')
            icone = 'Preditivo';
        end
    end

    methods (Static, Access = protected)
        function cabecalho = getHeaderImpl
            % Texto exibido no topo da janela do bloco
            cabecalho = matlab.system.display.Header('ControladorPreditivo', ...
                'Title', 'Controlador Preditivo (template)', ...
                'Text', ['Parametros definidos em config_Preditivo.m. ' ...
                         'No campo abaixo vai o nome da struct (cfg).']);
        end

        function modo = getSimulateUsingImpl
            % Execucao interpretada: o codigo roda no MATLAB
            modo = 'Interpreted execution';
        end

        function mostrar = showSimulateUsingImpl
            % Esconde da janela do bloco a opcao de trocar o modo de execucao
            mostrar = false;
        end
    end
end

% =========================================================================
% Funcoes locais (podem ser usadas pelos metodos da classe)
% =========================================================================
function mem = desloca(mem, novo)
% Insere "novo" no topo de um vetor-coluna de memoria e descarta o valor
% mais antigo (o ultimo). Memoria vazia continua vazia.
    if ~isempty(mem)
        mem = [novo; mem(1:end-1)];
    end
end

% ALUNO: se quiser, crie aqui funcoes auxiliares do seu algoritmo
% (ex.: function [Ej, Fj] = diofantina(A, N2) ... end)
