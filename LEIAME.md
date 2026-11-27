# Caldeira de Geração de Vapor — Controle Preditivo (GPC SISO)

Projeto de controle da Caldeira Aquatubular (UFMT — Controle Preditivo).

## Estrutura de Arquivos

```text
Projeto_Caldeira_MPC/
├── Caldeira_Controle.slx           # Modelo Simulink organizado em 3 áreas (planta editável no diagrama)
├── config_caldeira.m               # Sintonia dos GPCs (cfg_vap, cfg_niv) e parâmetros nominais
├── roda_caldeira.m                 # Script principal: simula o modelo e abre a animação ao vivo
├── anima_caldeira.m                # Animação física da caldeira com detecção real de falhas/explosão
├── ControladorPreditivo.m          # Bloco MATLAB System do GPC SISO (Vapor)
├── ControladorPreditivoDin.m       # Bloco MATLAB System do GPC SISO com limites dinâmicos (Nível)
├── ref_GPC.p                       # Biblioteca auxiliar usada por ControladorPreditivo.m
├── LEIAME.md                       # Este guia
├── colab/
│   └── Caldeira_Validacao_Colab.ipynb  # Protocolo de validação (Testes T1–T7)
└── docs/
    └── arquitetura_caldeira.html   # Diagrama interativo dos blocos Simulink
```

## Como Executar

1. Abra o MATLAB nesta pasta (`Projeto_Caldeira_MPC/`).
2. Execute `roda_caldeira` no Command Window (ou abra `Caldeira_Controle.slx` e clique em **Run**).
3. Para testar uma falha física real na animação (por exemplo, excesso de combustível na fornalha), altere `Gain_R` de `2.1` para `1.2` dentro de `Caldeira_Controle.slx` e execute `roda_caldeira` novamente.
