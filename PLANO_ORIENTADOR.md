# Plano — Estudo da Viga (revisão do orientador)

> Cópia local do plano publicado em `https://claude.ai/artifact/C2ZG2HXnoXc4XqeLV1deHu`.
> Ver esse link para a versão viva/editável.

O Prof. Thales pediu para transformar o exemplo da viga de um caso ilustrativo
único em uma investigação sistemática de seis fatores: passo de tempo,
condição inicial, amortecimento físico, amortecimento numérico, redução modal
e acúmulo de erro. Duas perguntas organizam o resultado final:

1. Como os métodos de integração se comportam quando a frequência relevante,
   o amortecimento físico, a condição inicial e o passo de tempo mudam?
2. Como uma redução modal que elimina frequências pouco relevantes altera a
   escolha do passo de tempo e a resposta calculada?

## Bateria de estudos A–F

| Estudo | Fixo | Variado | Métodos |
|---|---|---|---|
| A — Referência | CI impulso, ξ=0,2, h/T₁=0,01 (implícitos), h_RK=0,9h_cr | nada — reproduz o resultado já no paper | todos os 9 |
| B — Amortecimento | CI impulso; h de cada método (RK compartilhado, baseado no caso ξ=0,2) | ξ ∈ {0; 0,01; 0,2} | todos os 9 |
| C — CI modal, h fixo | h = 0,01·T₁ (absoluto); ξ=0,01 | modo ∈ {1,2,3,5,10} | todos os 9 |
| D — CI modal, h variável | ξ=0,01 | modo ∈ {1,3,10}; h/T_modo ∈ {0,01; 0,1; 0,9} | PE/AD para todos; resposta no tempo só nos 2 extremos |
| E — Rayleigh alternativo | CI impulso | (a₀,a₁) ∈ {(0,1e-4), (0,1e-5)} | todos os 9 |
| F — Redução modal | CI impulso, ξ=0,2 | modelo ∈ {completo, 4 modos, 3 modos} | todos os 9 + comparação item 32 (RK4) |

## Implementação (Setembro 2026)

Todos os estudos B–F foram implementados como scripts MATLAB independentes,
reaproveitando três funções novas (`viga_espectro.m`, `viga_pead_tabela.m`,
`viga_resposta_tempo.m`) para não duplicar a lógica de
`viga_main_simulacao.m`, que permanece intocado como referência do Estudo A.
Veja `README_viga.md` para a lista completa de arquivos novos/modificados e
as decisões tomadas nos pontos que o orientador deixou em aberto.

## Decisões em aberto — resolvidas com uma recomendação (ver README)

1. ξ=0,1 ou 0,2 no conjunto mínimo → mantido 0,2.
2. 3 ou 4 modos na redução → os dois rodados; 4 modos é o caso "seguro".
3. Rayleigh alternativo: os dois casos (a₁=1e-4 e 1e-5) rodados.
4. h fixo no Estudo C: absoluto (0,01·T₁), não relativo ao modo.
5. Estudo D: modos 1, 3, 10 (não os 5), para conter o custo em RK.
6. Referência de Bathe (item 8): ainda não localizada — pendente.
7. Estrutura no paper: Estudo F como possível Seção 4.5 separada — pendente
   até a escrita do texto.

## Riscos de custo computacional

O gargalo é sempre o mesmo: Runge-Kutta no limite de estabilidade custa
~2×10⁵ passos por método por caso. O Estudo D foi deliberadamente limitado
(poucos períodos, só os extremos de h) para não explodir esse custo; o
Estudo F na verdade *reduz* custo (menos GDL). Recomenda-se sempre rodar
`viga_espectro`/`viga_pead_tabela` (baratos, forma fechada) antes de decidir
quais casos merecem a simulação completa no tempo.
