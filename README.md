# Cantilever beam impulse response — Section 4.4 replacement

A uniform steel cantilever, clamped at `x = 0`, receives an impulsive
transverse velocity at the free end. The response is read at the free end,
i.e. as far as possible from the support. All nine integrators of the paper
are applied to the resulting MDOF finite element system.

> **Sept 2026 update.** The advisor reviewed the paper and asked for a much
> larger battery of simulations (damping sweep, modal initial conditions,
> modal reduction — see `PLANO_ORIENTADOR.md` for the full breakdown). This
> file now documents both the **original study (Estudo A, unchanged)** and
> the **new Estudos B–F**, added without touching the numerics of
> `viga_main_simulacao.m` — Estudo A remains the exact reference already in
> the paper. (Two things about its FIGURES did change, 21/09/2026: the
> spectral-radius figures put direct and Runge-Kutta schemes in one figure,
> and the displayed time window now follows the same rule as Estudos B–F.
> Every stored number — `res.time`, `res.timeRK`, `res.pead`, `res.stability`
> — is unchanged; see "Resposta inteira..." below.)

## Onde mexer nos parâmetros

Dois arquivos, com responsabilidades separadas — e **só esses dois**:

| Arquivo | Contém | Mexer aqui muda |
|---|---|---|
| `viga_params.m` | Material (E, ρ), geometria (L, b, h), malha (N), amortecimento de referência do Estudo A (ξ alvo, modos-âncora), CI de impulso (V0) | **o problema físico** — que viga está sendo resolvida |
| `viga_config_estudos.m` | Passo de tempo (`hRat`, `fatorSegurancaRK`), duração e janela de plot (`nPerImpl`, `nPerRK`, `X_ciclos_final_plot`, `tolDecaimento`, `ciclosDecaimentoRapido`, e as versões `_sem_xi`), figuras (`FS`, `FIGPOS`, `YCLIP`), amostragem (`nPontosExato`, `nMaxSalvar`, `nMaxPlot`), e os casos de cada estudo em `cfg.B` … `cfg.F` | **o experimento numérico** feito sobre ela |

Nenhum `viga_estudo_*.m` tem número de simulação chumbado: todos começam com
`cfg = viga_config_estudos();` e leem tudo de lá. Cada campo está comentado no
próprio arquivo explicando o que é, o que muda se for alterado e qual o custo.

Os campos de estudo (`p.xi_list`, `p.modal_modes`, `p.hD_ratios`,
`p.rayleigh_alt`, `p.reduction_modes`) **saíram** de `viga_params.m` em
20/09/2026 e viraram `cfg.B.xi_list`, `cfg.C.modos`, `cfg.D.hRatios`,
`cfg.E.rayleigh` e `cfg.F.reduction_modes`. A tabela de equivalência está no
cabeçalho de `viga_params.m`.

**Estudo A não lê daí, de propósito**: `viga_main_simulacao.m` segue com a
configuração própria embutida e sem mudança numérica, para continuar servindo
de regressão. Os valores equivalentes estão listados na seção 6 de
`viga_config_estudos.m`.

Para baixar o custo num teste rápido, o parâmetro que domina o tempo é
`cfg.nPerRK` (cada período custa ~2·10⁵ passos por método); `cfg.nPerImpl` é
barato.

## Running

Drop every file into the existing `Arquivo/` folder (no name collides) and run:

```matlab
viga_teste_janela.m       % ~10 s — confere a janela de plot ANTES de gastar tempo
viga_main_validacao.m     % checks the model and the integrators
viga_main_simulacao.m     % Estudo A — produces the results already in the paper
viga_main_estudos.m       % Estudos B-F — the advisor's new batteries (see below)
viga_tabelas.m            % consolidates every study into CSV tables
viga_report_parametros.m  % markdown com todos os parametros, para uso em IA
```

`viga_teste_janela.m` roda só a família de métodos diretos (segundos, não
minutos) nos três regimes de amortecimento e mostra onde a janela de plot
caiu em cada um. É o script para iterar em `cfg.X_ciclos_final_plot` sem
pagar o custo de uma bateria com Runge-Kutta. Salva `.eps` + `.png` de
conferência em `ImgViga/TesteJanela/`.

Figures go to `ImgViga/` as `.eps` (same export settings as `ImgCOBEM/`) **and
as a `.png` copy at 300 dpi** next to each one. The new studies write into
their own subfolders, `ImgViga/EstudoB/` … `ImgViga/EstudoF/`, each with its
own `.mat` file, so nothing here collides with Estudo A's output. The output
folder is always `<folder of the code>/ImgViga`, whatever MATLAB's current
folder is (see "Figuras: .eps + .png" below).

## Estudos B–F (advisor's revision)

Each study is a **self-contained script** (`viga_estudo_*.m`) that builds its
own model and calls three new, reusable functions instead of duplicating
Estudo A's code:

| New file | Role |
|---|---|
| `viga_ci_modal.m` | Displacement IC proportional to a single mode shape (Estudos C, D) |
| `viga_reducao_modal.m` | Modal truncation to the first `m` modes — exact reduced `K,M,C`, own `w_max` (Estudo F) |
| `viga_rho_beam.m` | Extracted, reusable version of the spectral-radius helper already used by Estudo A |
| `viga_espectro.m` | Spectral-radius figure for an arbitrary `(w, xi)` case — one figure with all nine schemes (direct solid, Runge-Kutta dashed) |
| `viga_janela_plot.m` | Decides which stretch of the response is plotted: whole response (strong damping), window anchored on the decay (weak damping) or the end of the run (no damping) |
| `viga_pead_tabela.m` | PE/AD table for an arbitrary `(w, xi)` case and method/step list |
| `viga_resposta_tempo.m` | Tip-displacement figures for an arbitrary case, with an optional **zoom inset** on a late oscillation (item 21) and an optional **time-axis cap** (item 22) |

Modified (backwards compatible — every existing call site still works
unchanged):

| File | Change |
|---|---|
| `viga_params.m` | Added `xi_list`, `rayleigh_alt`, `modal_modes`, `modal_modes_D`, `hD_ratios`, `reduction_modes` |
| `viga_amortecimento.m` | Now accepts two optional extra arguments `(a0, a1)` to build Rayleigh damping directly, instead of always solving for two anchor modes |

| Estudo | What it does | Script |
|---|---|---|
| **B** | Sweeps `xi = 0, 0.01, 0.2` on the SAME impulse IC as Estudo A. Runge-Kutta uses one shared step sized from the most-damped case (0.9·h_cr at ξ=0.2), exactly as Estudo A already does across its two cases — so RK1/RK2 are shown diverging at low damping and stable at high damping with the *same* h. Zoom inset produced for ξ=0.01. | `viga_estudo_B_amortecimento.m` |
| **C** | Displaces the beam along a single mode shape (modes 1,2,3,5,10; `p.modal_modes`) at ξ=0.01, with h **fixed in absolute terms** (0.01·T₁, the value Estudo A already uses) across every mode — isolates how h/T_modo changes as the excited mode gets higher. | `viga_estudo_C_modal_h_fixo.m` |
| **D** | Same modal ICs, but now h/T_modo itself is swept (0.01, 0.10, 0.90) for 3 representative modes (1, 3, 10). PE/AD/spectral radius (closed-form) are computed for every combination; full time-domain runs are limited to the two extreme ratios and 3 periods, wrapped in `try/catch`, to keep an intentionally-unstable case from stopping the script. | `viga_estudo_D_modal_h_variavel.m` |
| **E** | Rayleigh damping with `a0=0` and `a1=1e-4` / `1e-5` (not anchored to a target ξ) — reports the resulting ξ per mode instead of imposing it. | `viga_estudo_E_rayleigh_alt.m` |
| **F** | Modal reduction to 4 and 3 modes (`viga_reducao_modal.m`). Part 1 repeats the full stability/PE-AD/time-response pipeline for the full model and both reductions. Part 2 isolates the advisor's item 32 question directly: RK4 on (i) the full model, (ii) the 4-mode model at the SAME h, (iii) the 4-mode model at its OWN, larger, critical step — all three compared against the full model's exact response. | `viga_estudo_F_reducao_modal.m` |

Run all five with `viga_main_estudos.m` (skips nothing, but each script can
also be run on its own). `viga_tabelas.m` then reads whichever `.mat` files
exist and writes `ImgViga/Tabelas/tabela_{frequencias,amortecimento,
passos_criticos,N_erro}.csv` — safe to re-run at any point, studies not yet
executed are skipped with a warning rather than stopping the script.

### Decisions taken where the advisor left the choice open

(See `PLANO_ORIENTADOR.md`, Section 6, for the full reasoning.)

1. Damping set kept at `{0, 0.01, 0.2}` (not 0.1), for direct comparability
   with Estudo A.
2. Both 4-mode and 3-mode reductions are run; 4 modes is the one used for
   the item-32 comparison (3 modes puts ξ₃ already close to 1 when anchored
   at modes 1 and 3, which complicates the comparison).
3. Both Rayleigh alternatives (`a1=1e-4` and `1e-5`) are run — cheap, same
   mesh and methods as everything else.
4. Estudo C's fixed h is **absolute** (0.01·T₁), not relative to each
   mode's own period — the point of the study is to let h/T_modo vary.
5. Estudo D covers modes {1, 3, 10} only (not the full {1,2,3,5,10}), to
   keep the Runge-Kutta sweep bounded.
6. The Bathe *Finite Element Procedures* citation for item 8 (h ≲ T_min/10)
   still needs a page reference — not automated here.
7. The reduction study became its own Estudo F rather than folding into
   Estudo A's Section 4.4 text; how the paper's LaTeX is organised (a new
   4.5, or an expanded 4.4) is left for when the writing starts.

## Results produced

Only three, as requested:

1. **Spectral radius** of every scheme on the FE beam — 2 figures
   (`viga_rho_{xi0,xi02}.eps`), one per damping case, each with the direct
   and the Runge-Kutta schemes together (they used to be 4 figures,
   `viga_rho_{rk,direto}_*.eps`, mirroring Figs. 2 and 3).
2. **PE and AD** for every mode, printed for the first three and stored in
   full.
3. **Tip displacement** against the exact response — 4 figures
   (`viga_desloc_{direto,rk}_{xi0,xi02}.eps`).

All displacement axes are clipped at `±1.6` of the analytical peak, so an
unstable curve leaves the frame instead of flattening everything else.

## Saved data

Everything computed goes into `ImgViga/viga_resultados.mat` (saved with
`-v7.3`), next to the figures:

```matlab
S = load('ImgViga/viga_resultados.mat');
res = S.res;
```

| Field | Content |
|---|---|
| `res.model` | `K`, `M`, `C`, `w`, `Phi`, damping data, initial conditions, `T1`, `wmax` |
| `res.config` | every parameter of the run, plus step counts and a timestamp |
| `res.stability.hSweep` | time steps of the spectral radius sweep |
| `res.stability.rho.<case>.<method>` | spectral radius curves |
| `res.stability.hcr_xi0`, `.hcr_xi02` | critical time steps |
| `res.stability.hcr_table` | the same, as a MATLAB `table` |
| `res.pead.PE.<case>` | `[nMethods × nModes]` period elongation |
| `res.pead.AD.<case>` | `[nMethods × nModes]` amplitude decay |
| `res.pead.rho.<case>` | `[nMethods × nModes]` scalar spectral radii |
| `res.pead.table.<case>` | first three modes, as a MATLAB `table` |
| `res.time.<case>.<method>` | implicit runs: `t`, `h`, `u_tip`, `v_tip`, `a_tip`, `err` |
| `res.time.<case>.exact` | exact reference on the same grid |
| `res.timeRK.<case>.<method>` | Runge-Kutta runs, decimated |

`<case>` is `xi0` or `xi02`; `<method>` is `DC`, `Houbolt`, `Newmark`,
`Wilson`, `Bathe`, `RK1`…`RK4`.

Three notes on the stored histories:

- Only the **tip DOF** is stored by default. Set `saveFullField = true` to
  also store all 20 DOFs of the implicit runs.
- The RK histories are decimated to `nSaveRKMax` points (default `2e5`),
  since a 3-period run is already ~654 000 steps. Set it to `Inf` to keep
  every step. The flag `res.timeRK.<case>.<method>.decimated` records whether
  decimation happened, and `.u_tip_max` keeps the true peak of the
  undecimated run.
- **Central Differences is computed and stored** at `h/T₁ = 0.01` even though
  it is kept out of the figure — at that step it is 53× above its limit and
  would flatten the comparison. Its diverging history is in
  `res.time.<case>.DC`.

The two PE/AD tables can be exported straight to LaTeX or CSV:

```matlab
writetable(res.pead.table.xi02, 'pead_xi02.csv')
```

## File map

| File | Role |
|---|---|
| `viga_params.m` | geometry, material, mesh, damping target, impulse |
| `viga_matrizes_EF.m` | 2-node beam element, consistent mass, assembly, clamping |
| `viga_autovalores_EF.m` | `K φ = ω² M φ`, mass-normalised |
| `viga_amortecimento.m` | Rayleigh `C = a₀M + a₁K` |
| `viga_ci_impulso.m` | impulse initial conditions |
| `viga_solucao_exata.m` | **exact** modal response of the semi-discrete system |
| `viga_modos_analiticos.m` | continuous Euler-Bernoulli modes (validation only) |
| `viga_modelo.m` | one-call assembly |
| `viga_operador.m` | scalar `Â` of all nine schemes (Sec. 2.2 + Appendix A) |
| `viga_pe_ad.m` | PE and AD per Section 3.3 |
| `viga_h_critico.m` | critical step by bisection on `ρ(Â) ≤ 1` |
| `mdof_{dc,houbolt,newmark,wilson,bathe,rk}.m` | matrix integrators |
| `viga_main_validacao.m` | script 1 |
| `viga_main_simulacao.m` | script 2 |

### Why `viga_operador.m` is enough for the beam

Rayleigh damping is proportional, so `M`, `C`, `K` are simultaneously
diagonalised by the mode shapes and every integrator's MDOF operator is
block-diagonal in the modal basis. Hence

```
ρ_MDOF(h) = max over modes k of ρ_scalar(ω_k, ξ_k, h)
```

exactly. Script 1 verifies this against the full `2n` block operator — the
difference is at round-off level (`~1e-15`). This is what lets the spectral
radius, PE and AD of the beam be computed with the same scalar operators
already used in Sections 4.1–4.3.

## Model data

Steel `E = 210 GPa`, `ρ = 7850 kg/m³`; `L = 1 m`, section `50 × 10 mm`
(`L/h = 100`); 10 elements, 20 DOF; `V₀ = 1 m/s` at the tip.

| | value |
|---|---|
| `ω₁` | 52.4971 rad/s (`T₁ = 0.119686 s`) |
| `ω₂` | 328.994 rad/s |
| `ω₃` | 921.206 rad/s |
| `ω_max` | 8.9402e4 rad/s (`ω_max/ω₁ ≈ 1703`) |

Rayleigh `ξ = 0.2` anchored at modes 1 and 3: `a₀ = 19.867`,
`a₁ = 4.1071e-4`, giving `ξ₁ = 0.200`, `ξ₂ = 0.098`, `ξ₃ = 0.200`,
`ξ₂₀ = 18.36`.

## Expected output

### Critical time steps

| method | `h_cr` (ξ=0) | `h_cr/T₁` | `h_cr` (ξ=0.2) | `h_cr/T₁` |
|---|---|---|---|---|
| RK1 | 1.58e-09 * | 1.32e-08 * | 6.0970e-07 | 5.09e-06 |
| RK2 | 1.88e-07 * | 1.57e-06 * | 6.0970e-07 | 5.09e-06 |
| RK3 | 1.9374e-05 | 1.62e-04 | 7.6601e-07 | 6.40e-06 |
| RK4 | 3.1637e-05 | 2.64e-04 | 8.4910e-07 | 7.09e-06 |
| DC | 2.2371e-05 | 1.87e-04 | 2.2371e-05 | 1.87e-04 |

Undamped limits reproduce the closed forms exactly: `DC = 2/ω_max`,
`RK3 = √3/ω_max`, `RK4 = √8/ω_max`.

\* Tolerance artifacts, not real limits: for RK1 and RK2 undamped,
`ρ > 1` for every `h > 0`. Quote them as "no stable step".

### PE and AD, ξ = 0 (percent)

Implicit at `h/T₁ = 0.01`, RK at `h = 4.58e-6 T₁`.

| method | PE m1 | AD m1 | PE m2 | AD m2 | PE m3 | AD m3 |
|---|---|---|---|---|---|---|
| Central Differences | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Houbolt | 0.0000 | 0.0767 | 0.0182 | 11.2851 | 0.6096 | 50.0852 |
| Newmark | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Wilson-θ | 0.0000 | 0.0130 | 0.0009 | 2.5944 | 0.0751 | 21.6136 |
| Bathe | 0.0000 | 0.0005 | 0.0000 | 0.1296 | 0.0007 | 2.3772 |
| RK1 | 0.0000 | −0.0091 | 0.0000 | −0.0567 | 0.0000 | −0.1590 |
| RK2–RK4 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

### PE and AD, ξ = 0.2 (ratios; 1 = physical period / decay reproduced)

| method | PE m1 | AD m1 | PE m2 | AD m2 | PE m3 | AD m3 |
|---|---|---|---|---|---|---|
| Central Differences | 0.9998 | 1.0002 | 0.9932 | 1.0039 | 0.9337 | 1.0673 |
| Houbolt | 1.0014 | 1.0063 | 1.0643 | 0.9645 | 1.3576 | 0.7493 |
| Newmark | 1.0003 | 1.0008 | 1.0123 | 1.0156 | 1.0821 | 1.2312 |
| Wilson-θ | 1.0005 | 1.0026 | 1.0271 | 1.0192 | 1.1764 | 1.1619 |
| Bathe | 1.0001 | 1.0004 | 1.0063 | 1.0066 | 1.0444 | 1.0952 |
| RK1 | 1.0000 | 1.0001 | 1.0000 | 1.0006 | 0.9999 | 1.0016 |
| RK2–RK4 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |

## What the example shows

1. **Physical damping makes RK1 and RK2 usable.** Undamped, neither has a
   stable step; with `ξ = 0.2` a finite `h_cr` appears. The tip-displacement
   figure pair shows this directly: at the same step, RK1 and RK2 leave the
   frame at `ξ = 0` and track the exact response at `ξ = 0.2`.
2. **Physical damping penalises RK3 and RK4** — their stability range shrinks
   by more than an order of magnitude.
3. **Central Differences is insensitive to damping** — `h_cr` unchanged to five
   digits, matching the SDOF `h/Tₙ ≈ 0.3175` observation.
4. **The limit is set by `ω_max`, not `ω₁`.** Conditionally stable schemes need
   `h/T₁ ~ 1e-4` to `1e-6` on the mesh, orders of magnitude below what the SDOF
   analysis alone would suggest.
5. **The PE/AD ordering of Sections 4.2–4.3 is reproduced.** Undamped, amplitude
   decay at mode 3 goes Houbolt (50.1%) > Wilson-θ (21.6%) > Bathe (2.4%) >
   Newmark = DC (0%); period elongation goes Houbolt > Wilson-θ > Bathe, with DC
   the only scheme that shortens the period when damping is present
   (`PE = 0.9337` at mode 3).
6. **`AD > 1` under damping for Newmark, Bathe and Wilson-θ** — amplitude growth
   relative to the physical decay while remaining stable — whereas Houbolt gives
   `AD < 1`. Exactly the method-dependent, counterintuitive behaviour claimed in
   Section 4.3, now on a real structure.
7. **RK dispersion and dissipation errors are negligible** at the steps their
   stability forces on them (`PE = AD = 1.0000` to four digits). Their limitation
   on this problem is purely the step-size restriction, not accuracy.

## Tuning

Everything lives in `viga_params.m`: `p.N` (mesh), `p.mode_i`/`p.mode_j`
(Rayleigh anchors), `p.xi_target`, `p.V0`.

In the configuration block at the top of `viga_main_simulacao.m`:

| Variable | Default | Meaning |
|---|---|---|
| `nPerImpl` | 10 | periods **simulated**, implicit schemes |
| `nPerRK` | 3 | periods **simulated**, Runge-Kutta schemes |
| `n_amostras` | 10 | periods **plotted**, counted from the end of the run |
| `YCLIP` | 1.6 | displacement axis limit, normalised |
| `nSaveRKMax` | 2e5 | cap on stored RK samples (`Inf` = store every step) |
| `saveFullField` | false | also store all DOFs of the implicit runs |
| `FS`, `FIGPOS` | 18, `[0 0 20 14]` | font size and figure size |

`hRat` (0.01) sets the implicit time step further down, in Section 2.

### `n_amostras`

The full run is always carried out; only the displayed window is trimmed.
With `nPerImpl = 100` and `n_amostras = 10`, the implicit schemes integrate
100 periods and the figures show periods 90 to 100. Setting
`n_amostras >= nPerImpl` shows the whole run.

Each figure is normalised by the peak of the exact response **inside the
displayed window**, not by the global peak. Without this, a damped tail
plotted after 90 periods would collapse onto the axis; with it, the window
stays centred on the analytical curve, which is the point of looking at the
tail in the first place.

The same parameter trims the Runge-Kutta figures, relative to `nPerRK`.

**Cost warning.** `nPerImpl` is cheap: at `h/T₁ = 0.01` one period is ~100
steps. `nPerRK` is not: stability forces `h ≈ 5.5e-7 s`, so each period costs
about 2.2e5 steps *per method per damping case* — roughly 1.7 million steps
per period across the 4 RK schemes and 2 cases. The default `nPerRK = 3`
already runs ~5.2 million steps; raising it to 10 triples that.

## Execução completa (20/09/2026)

Todos os Estudos A–F e `viga_tabelas.m` foram executados de ponta a ponta no
MATLAB do usuário (R2026a), com sucesso. Resultados em `ImgViga/` (uma
subpasta por Estudo, mais `ImgViga/Tabelas/` com os 4 CSVs). Três problemas
reais apareceram durante essa primeira execução e foram corrigidos no
código (não são apenas notas — as correções já estão nos arquivos):

1. **`viga_main_estudos.m` só rodava o primeiro Estudo.** Cada
   `viga_estudo_*.m` começa com `clc; clear; close all;` para poder rodar
   sozinho; como `run(script)` executa no MESMO workspace de quem chamou,
   esse `clear` apagava as variáveis do próprio laço do driver (`scripts`,
   `k`). Corrigido isolando a chamada dentro de uma função local
   (`run_estudo_isolado`), que dá a `run()` um workspace próprio para
   limpar. Efeito só aparece rodando `viga_main_estudos.m` inteiro (rodar
   cada `viga_estudo_*.m` avulso, como sempre foi possível, não tinha esse
   problema).

2. **Estouro de memória no Runge-Kutta (Estudos B e F).** O passo do RK é
   limitado pela estabilidade (não pela precisão) e, na malha completa,
   fica da ordem de `1e-6`–`1e-7` s: para 12 períodos isso é ~1–3 milhões
   de passos. `mdof_rk.m` sempre alocava as três matrizes `U`, `V`, `A`
   de `[n x nSteps+1]`, mesmo quando só `U` era pedido — no computador
   usado (memória disponível ~0.4–1.5 GB durante a sessão, ver `memory()`
   no MATLAB), isso estourava. Duas correções, sem mudar nenhum resultado
   numérico:
   - `mdof_rk.m` agora só calcula/aloca `V` e `A` quando o chamador pede
     esses argumentos de saída (`nargout`) — quem já chamava com
     `[t,U,V,A]=mdof_rk(...)` (Estudo A) continua recebendo exatamente o
     mesmo resultado.
   - Novo `viga_rk_stream_tip.m`: faz a mesma recursão do RK (mesmo
     operador `Ahat`, resultado numericamente idêntico) mas nunca guarda o
     histórico completo — só os pontos já decimados (≤ 2·10⁵, os mesmos
     que `viga_resposta_tempo.m`/o item 32 do Estudo F já usavam para
     salvar/plotar) e o `max(|resposta|)` corrido. É o que
     `viga_resposta_tempo.m` e o item 32 (`viga_estudo_F_reducao_modal.m`)
     usam agora no lugar de `mdof_rk.m` direto.

3. **Títulos das figuras novas quebrando na exportação.** `viga_espectro.m`
   e `viga_resposta_tempo.m` montavam o título com `'Interpreter','latex'`,
   mas os rótulos vindos de cada Estudo (`'\xi = 0.2'`, `'a_0=0, a_1=...'`,
   `'modo 3, \xi=0.01'`) misturam texto comum com comandos LaTeX sem
   `$...$` — LaTeX de verdade exige isso, e a exportação (`print -depsc2`)
   falhava silenciosamente nesse texto (aviso "must have valid interpreter
   syntax", título ficava ausente/errado no `.eps`). Trocado para
   `'Interpreter','tex'` nos títulos (os eixos, que já usavam `$...$`
   corretamente, continuam em `'latex'`) — mais tolerante a essa mistura,
   sem mudar a aparência dos eixos.

Também: `viga_estudo_B_amortecimento.m` não estava salvando
`res.damp.(caso)` (só `viga_estudo_E` salvava), o que quebrava
`tabela_amortecimento.csv` em `viga_tabelas.m`. Corrigido — a tabela agora
sai com os 3 casos de B mais os 2 de E.

**Pontos ainda em aberto** (sem mudança nesta rodada): a referência
bibliográfica do método de Bathe (item 8) e a decisão sobre o Estudo F
virar uma seção própria no paper (ambos já listados acima).

## Centralização da configuração (20/09/2026, depois da primeira execução)

Os blocos de configuração que estavam repetidos no topo de cada
`viga_estudo_*.m` (`FS`, `FIGPOS`, `YCLIP`, `nPerImpl`, `nPerRK`,
`n_amostras`, `hRat`, o fator `0.90` do passo crítico) foram unificados em
`viga_config_estudos.m`, junto com os parâmetros de estudo que estavam em
`viga_params.m`. Ver a seção "Onde mexer nos parâmetros" no começo deste
arquivo.

Além de eliminar a duplicação, isso destravou dois parâmetros que antes
estavam invisíveis dentro das funções:

- `cfg.nMaxSalvar` / `cfg.nMaxPlot` / `cfg.nPontosExato` — eram literais
  dentro de `viga_resposta_tempo.m` (2e5, 5000, 4000). Continuam com esses
  valores por padrão, mas agora são configuráveis. Eles controlam só
  quantos pontos são **guardados e desenhados**; todo passo continua sendo
  integrado.
- `cfg.D.ratiosSimuladas` — a escolha "só os extremos de h/T_modo ganham
  simulação no tempo" estava implícita no código do Estudo D. Agora é
  `'extremos'` (padrão) ou `'todas'`.

Cada `.mat` de estudo passou a guardar `res.cfg` com a configuração
completa da corrida, para rastreabilidade.

Validação: Estudos B–F e `viga_tabelas.m` reexecutados após o refactor —
todos os números (passos críticos, PE/AD, ξ por modo, erros do item 32)
idênticos aos de antes.

## Janela de plot condicionada ao decaimento (20/09/2026)

Antes, a figura sempre mostrava os últimos `n_amostras` ciclos da
simulação. Com amortecimento isso é ruim: a resposta morre no meio do
caminho e a figura vira uma linha reta. Agora a janela é **ancorada no
decaimento**.

**Com amortecimento (ξ ≠ 0):** o código procura o último instante em que a
resposta **exata** ainda está acima de `cfg.tolDecaimento` vezes o próprio
pico, e mostra os `cfg.X_ciclos_final_plot` ciclos anteriores a esse ponto.
A detecção usa sempre a solução exata, nunca a resposta de um método — um
esquema instável divergiria e arrastaria a janela para o lugar errado,
justamente no caso em que a figura mais precisa estar certa.

**Sem amortecimento (ξ = 0):** não há decaimento a detectar, então a janela
são os últimos `cfg.X_ciclos_final_plot_sem_xi` ciclos, com contagem de
ciclos calculados própria (`nPerImpl_sem_xi`, `nPerRK_sem_xi`). A escolha
entre os dois conjuntos é automática — `viga_resposta_tempo.m` olha o ξ do
caso que recebeu, nenhum estudo precisa decidir isso.

Cada figura imprime no console onde a janela caiu, e avisa quando ciclos
foram calculados depois do decaimento e jogados fora:

```
[janela direct schemes] resposta decaiu a 1.0e-03 do pico em t/T1 = 5.47;
                        plotando t/T1 em [0.47, 5.47]
[janela direct schemes] AVISO: 24.5 ciclos calculados depois do decaimento
                        foram descartados - da para reduzir nPer para ~7
```

**Quantos ciclos são necessários.** A envoltória cai como exp(−ξ·2π·n)
depois de n ciclos, então para chegar a `tolDecaimento`:

    n = −ln(tolDecaimento) / (2π·ξ)

Com tol = 10⁻³: ξ=0.2 morre em ~5.5 ciclos, ξ=0.05 em ~22, ξ=0.01 só em
~110. Ou seja, **para ξ=0.01 a detecção não dispara** com `nPerImpl = 30` —
o código avisa e mostra o fim da simulação. Se quiser que dispare nesse
caso, ou sobe `nPerImpl` (caro no Runge-Kutta) ou afrouxa `tolDecaimento`
para 10⁻² (~22 ciclos).

### Erro relativo dos Runge-Kutta

Aproveitando, a família RK passou a gravar `err` como a direta já fazia
(antes só guardava `u_tip_max`, e a coluna de erro da tabela N×erro saía
vazia para RK1–RK4). A solução exata é avaliada nos mesmos instantes
salvos, uma vez só, reaproveitada pelos quatro métodos.

### Tudo em inglês

Títulos, eixos e legendas de todas as figuras estão em inglês: `direto` →
`direct schemes` nos títulos de raio espectral, `modo N` → `mode N`,
`modelo: completo` → `model: full (20 DOF)`, e as legendas do item 32
viraram `Case 1: full, h_1` etc. Os identificadores internos (nomes de
campo no `.mat`, sufixos de arquivo) seguem sendo os mesmos onde já eram
ASCII, para não quebrar `viga_tabelas.m`.

## Regra de janela também no Estudo A, e paleta de cores única (21/09/2026)

**Janela do Estudo A.** As figuras de deslocamento de `viga_main_simulacao.m`
agora usam a mesma `viga_janela_plot.m` dos Estudos B–F, em vez da lógica
antiga ("sempre os últimos 10 ciclos da simulação"). Efeito prático: no caso
ξ=0,2 (que decai a 1% do pico em ~3,7 ciclos, dentro do limiar de 6), a
figura passa a mostrar a resposta **inteira** — inclusive o transiente
inicial (overshoot do Wilson-θ, por exemplo), que antes ficava fora da
janela de "últimos 10 ciclos". O caso ξ=0 é undamped, cai no regime "fim da
simulação" e o resultado visual não muda. O skip de Central Differences da
figura também passou a ser decidido dinamicamente (compara `hImpl` contra o
`h_cr` de DC, já calculado na Seção 2), em vez de fixo — sem efeito prático
aqui, já que DC sempre esteve acima do próprio limite nesse `hImpl`, mas
deixa A consistente com B–F. `tolDecaimento`/`ciclosDecaimentoRapido` ficam
**duplicados** no cabeçalho de A (não lidos de `viga_config_estudos.m`, por
quê: ver a nota "Estudo A não lê daí, de propósito" acima) — mantenha-os
em sincronia manualmente se mudar os defaults lá. Nada além da FIGURA
muda: `res.time`, `res.timeRK`, `res.pead`, `res.stability`, os `.err`, os
passos críticos — todos idênticos a antes. Dois campos novos, só
informativos: `res.window.<caso>.direto` / `.rk` (o diagnóstico da janela,
mesmo formato usado em `tresp.config.janelaDireto` nos Estudos B–F).

`viga_print_janela.m` (a linha de diagnóstico no console, "[janela ...]")
foi extraída de dentro de `viga_resposta_tempo.m` para um arquivo próprio,
para que o Estudo A também pudesse chamá-la — mesmo padrão já usado por
`viga_rho_beam.m`.

**Paleta de cores única.** Todo arquivo que desenha mais de um método (A e
B–F, resposta no tempo e raio espectral) agora lê as cores de
`viga_cores.m`, em vez de ter o array `coresDir`/`coresRK` repetido em cada
arquivo:

```matlab
coresDir = [1.0 0.0 0.0;   % DC           - vermelho
            0.0 0.2 1.0;   % Houbolt      - azul
            0.1 0.7 0.1;   % Newmark      - verde
            0.6 0.0 0.6;   % Wilson-Theta - roxo
            1.0 0.5 0.0];  % Bathe        - laranja

coresRK  = [0.6 0.0 0.0;   % RK1 - vermelho escuro (NÃO o 'r' puro - ver abaixo)
            0.0 0.0 1.0;   % RK2 - azul   ('b')
            0.0 1.0 0.0;   % RK3 - verde  ('g')
            1.0 0.0 1.0];  % RK4 - magenta ('m')
```

As cores dos métodos diretos já eram exatamente essas em todo o projeto —
nada mudou aí. As do Runge-Kutta viraram as primárias do MATLAB (`r`, `b`,
`g`, `m`), **exceto RK1**: `[1,0,0]` é literalmente idêntico ao vermelho da
DC, e nas figuras que mostram diretos e RK juntos (raio espectral) as duas
curvas ficariam com a mesma cor — só essa foi trocada, para um vermelho
mais escuro (grena, `[0.6,0,0]`). Os outros pares parecidos (RK2/Houbolt,
RK3/Newmark, RK4/Wilson-Theta) não são literalmente iguais e continuam
distinguíveis pelo estilo de linha (RK tracejado, direto contínuo) nas
figuras onde aparecem juntos; nas figuras de resposta no tempo, direto e RK
nunca aparecem na mesma figura, então não há conflito ali de qualquer jeito.
Isso muda a aparência das figuras de raio-espectral combinado (RK deixa de
ser preto/cinza/ciano/marrom e passa a ser grena/azul/verde/magenta) e das
figuras de resposta-no-tempo do RK (RK2/RK3/RK4 ficam mais saturados que
antes; RK1 vai de vermelho puro para grena).

## Resposta inteira com amortecimento forte + raio espectral combinado (21/09/2026)

**1. Resposta inteira quando o amortecimento é forte.** Se a resposta exata
cai abaixo de `cfg.tolDecaimento` do próprio pico em até
`cfg.ciclosDecaimentoRapido` ciclos de T₁ (padrão 6), a figura mostra a
resposta **inteira**, de t = 0 até o instante do decaimento, em vez de uma
janela de `X_ciclos_final_plot` ciclos. É assim que o transiente inicial —
onde os métodos mais diferem entre si nesse caso (p. ex. o overshoot do
Wilson-θ) — deixa de ficar fora da figura. Com decaimento mais lento vale a
janela ancorada no decaimento de antes; sem amortecimento, o fim da
simulação. `ciclosDecaimentoRapido = 0` desliga a regra.

    n = −ln(tol) / (2π·ξ)      com tol = 1e-2:  ξ=0.2 → 3.7 ciclos (resposta inteira)
                                                ξ=0.1 → 7.3 ciclos (janela)

A lógica mora em `viga_janela_plot.m` e é usada por `viga_resposta_tempo.m`
e pelo item 32 do Estudo F (que antes plotava os `nPerRK` ciclos inteiros,
quase todos achatados com ξ=0.2). O Estudo A **não** usa essa regra (nem a
config central): continua mostrando os últimos 10 ciclos.

**2. Raio espectral: diretos e Runge-Kutta na mesma figura.** Todas as
figuras de raio espectral (Estudo A e B–F) agora saem como **uma figura por
caso** com os nove esquemas: diretos em linha cheia, Runge-Kutta tracejados.
Nove curvas pedem nove cores, então aqui as cores do RK (preto, cinza,
ciano, marrom) não repetem as dos diretos — nas figuras de deslocamento
continuam as cores de antes. Arquivos: `viga_rho_<caso>.eps`
(antes `viga_rho_{rk,direto}_<caso>.eps`). A estrutura
`stab.rho.{direto,rk}.<método>` é a mesma. Só a figura de raio espectral
foi combinada; as de deslocamento seguem separadas por família.

### Correções feitas na revisão

- **Falso "decaimento" na detecção da janela.** O código achava o último
  ponto com `|u| > tol·pico`; se o último ponto da simulação caísse num
  cruzamento por zero da oscilação, o "último ponto acima da tolerância" era
  o penúltimo, e uma resposta ainda longe de decair (ξ=0.01 aos 50 ciclos
  está a ~4% do pico) era dada como decaída, com mensagem errada no console.
  Agora só conta como decaída a resposta que fica abaixo da tolerância por
  pelo menos um período T₁ inteiro antes do fim. Com a regra da resposta
  inteira isso passaria a importar: numa simulação curta (≤ 6 ciclos) uma
  resposta que ainda não decaiu seria tratada como "decaiu rápido".
- **Normalização do `err`.** `err` (coluna de erro da tabela N×erro) era
  dividido pelo pico da resposta exata *dentro da janela plotada*. Com
  amortecimento a janela cai no rabo da resposta, onde o pico é uma fração
  pequena do global, então o erro saía inflado por esse fator. Agora é
  normalizado pelo pico **global** da exata, a mesma definição do Estudo A.
  Não afeta ξ=0. Os `.mat` gerados antes dessa correção têm `err` da
  definição antiga para casos amortecidos.
- **Comentários da configuração.** `viga_config_estudos.m` dizia que o
  Estudo D conta ciclos em T do *modo excitado*; o código conta em T₁ (o
  mesmo eixo `t/T₁` dos outros estudos). Foi o comentário que mudou, não o
  comportamento.

## Figuras: .eps + .png, pasta fixa e conferência (21/09/2026)

**Sintoma investigado:** "o código não salvou as imagens". Não era o
overwrite (`print` sobrescreve normalmente); o que havia de fato:

- o destino era `fullfile(pwd,'ImgViga',...)`: com o `pwd` do MATLAB fora da
  pasta do código, as figuras iam parar em outra pasta sem aviso e o
  `ImgViga` "de sempre" ficava com os arquivos antigos;
- um arquivo aberto em outro programa (visualizador de EPS/PNG, OneDrive
  sincronizando) faz o `print` falhar — e a falha só aparecia como erro no
  meio da rodada;
- `viga_teste_janela.m` reimprimia **todas** as figuras abertas no mesmo
  `preview_<caso>.png`; como não eram fechadas entre os casos, os três
  previews saíam idênticos (a figura do primeiro caso);
- sobras de nomes antigos (`viga_rho_{rk,direto}_*.eps`,
  `..._teste_amortecido.eps`) ficavam misturadas com os arquivos novos.

**O que mudou:**

- `viga_salva_figura.m` é a única porta de saída de figuras. Grava
  `<nome>.eps` (vetorial, como antes) **e** `<nome>.png` (300 dpi), e
  confere que os dois existem, não estão vazios e foram gravados *agora*.
  Se algo falhar (figura inválida, arquivo bloqueado, sem permissão), emite
  `warning('viga_salva_figura:falhou', ...)` com o caminho absoluto e a
  rodada continua — não passa em silêncio.
- `viga_pasta_saida.m`: a saída é sempre `<pasta do viga_params.m>/ImgViga`,
  independente do `pwd`. Cada script imprime no início
  `Figuras (eps + png) e .mat em: <caminho absoluto>`.
- `viga_lista_obsoletos.m`: no fim de cada estudo, lista (sem apagar) os
  `.eps`/`.png` da pasta que **não** foram regravados nessa execução — em
  geral sobras de um nome que mudou. Apagar essas sobras é decisão sua.
- `viga_teste_janela.m`: um PNG por caso (`viga_desloc_direto_teste_<caso>.png`),
  figuras fechadas entre os casos.

## Report de parâmetros para IA

`viga_report_parametros.m` lê os `.mat` de todos os estudos e escreve
`ImgViga/Tabelas/PARAMETROS_ESTUDOS.md`: modelo físico, frequências
naturais, configuração numérica completa, e uma seção por estudo com a
pergunta que ele responde, os casos varridos, os passos efetivamente
usados, os passos críticos e a tabela N×erro por método. Termina com um
glossário dos símbolos. Escrito em inglês e auto-contido — a ideia é poder
colar o arquivo inteiro numa IA como contexto, sem precisar explicar nada
à parte. Estudos ainda não rodados são pulados com aviso, então dá para
gerar a qualquer momento.
