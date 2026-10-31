# Time integration schemes applied to a cantilever beam

**Author:** Murilo Barros Scalari (RA 256830)
**Advisor:** Prof. Thales Peixoto
**Work:** paper *Time Integration Schemes*

This repository holds the MATLAB code for the beam example in the paper. A steel
cantilever, discretized with Euler-Bernoulli finite elements, is solved with nine
time integration schemes, and each one is compared against the exact response of
the semi-discrete system.

Two questions drive the work:

1. How do the schemes behave as the relevant frequency, the physical damping, the
   initial condition and the time step change?
2. How does a modal reduction, which removes less relevant frequencies, change the
   choice of time step and the computed response?

## Project structure

```text
cantilever-beam-study/
├── cantilever-beam-study.prj   MATLAB project (adds src/ to the path when opened)
├── README.md
├── src/
│   ├── viga.m                  entry point, validation and Studies A to F
│   ├── modelo.m                beam, configuration, finite elements, exact solution
│   ├── metodos.m               integrators, spectral radius, PE/AD, critical step
│   ├── graficos.m              figures, plot window, output files
│   └── resultados.m            tables, parameter report and tests
├── ImgViga/                    figures, .mat files and tables (created on run)
└── resources/                  MATLAB project metadata
```

Each library file returns a struct of handles to its public functions, so one file
can call another as `modelo().params()` or `metodos().h_critico(K, M, C, 'DC')`. Functions
that are not in that list stay private to the file.

## How to run

Requirement: MATLAB (run on R2026a).

Open `cantilever-beam-study.prj` in MATLAB, or add the folder to the path:

```matlab
addpath('src')
```

Each stage is called by name. This is the recommended order:

```matlab
viga('validacao')   % checks the model and the integrators
viga('testes')      % quick assert-based tests
viga('janela')      % checks the plot window, direct schemes only
viga('A')           % Study A, the reference case in the paper
viga('estudos')     % Studies B to F
viga('tabelas')     % consolidates the results into CSV files
viga('relatorio')   % report with every parameter
```

A single study also runs on its own, with `viga('B')` to `viga('F')`. `help viga`
prints this list.

Runge-Kutta dominates the cost. Its step is limited by stability, so each period
costs about 2·10⁵ steps per method. For a quick test, lower `cfg.nPerRK` and
`cfg.nPerRK_sem_xi` in `config_estudos` (`modelo.m`).

## The model

Steel beam clamped at `x = 0` and free at `x = L`. The impulse is an initial
transverse velocity at the tip, and the response is read at the tip.

| Property | Value |
|---|---|
| Young's modulus `E` | 210 GPa |
| Density `ρ` | 7850 kg/m³ |
| Length `L` | 1 m |
| Section `b × h` | 50 × 10 mm (`L/h = 100`) |
| Mesh | 10 elements, 20 DOF |
| Initial velocity `V₀` | 1 m/s |

| Frequency | Value |
|---|---|
| `ω₁` | 52.4971 rad/s (`T₁ = 0.119686 s`) |
| `ω₂` | 328.994 rad/s |
| `ω₃` | 921.206 rad/s |
| `ω_max` | 8.9402·10⁴ rad/s (`ω_max/ω₁ ≈ 1703`) |

In the reference case, Rayleigh damping `C = a₀M + a₁K` is anchored at `ξ = 0.2` on
modes 1 and 3. This gives `a₀ = 19.867`, `a₁ = 4.1071·10⁻⁴`, `ξ₂ = 0.098` and
`ξ₂₀ = 18.36`.

### The nine schemes

| Family | Schemes | Stability |
|---|---|---|
| Direct | Newmark, Houbolt, Wilson-θ (θ = 1.4), Bathe | unconditional |
| Direct | Central Differences | conditional |
| Runge-Kutta | RK1 (Euler), RK2 (modified Euler), RK3 (Heun), RK4 | conditional |

### Why the scalar operator is valid for the beam

Rayleigh damping is proportional, so `M`, `C` and `K` are diagonalized by the same
modes and every integrator's operator is block-diagonal in the modal basis:

```text
ρ_MDOF(h) = max over modes k of ρ_scalar(ω_k, ξ_k, h)
```

`viga('validacao')` checks this equality against the full operator of order `2n`,
and the difference stays at round-off level. This is why the spectral radius, PE and
AD of the beam use the same scalar operators as the paper.

## The studies

| Study | Held fixed | Varied | Command |
|---|---|---|---|
| **A**, reference | impulse IC, `h/T₁ = 0.01` for implicit schemes, `h_RK = 0.9·h_cr` | `ξ = 0` and `0.2` | `viga('A')` |
| **B**, physical damping | impulse IC, a single `h_RK` | `ξ ∈ {0; 0.002; 0.005; 0.01; 0.2}` | `viga('B')` |
| **C**, modal IC, fixed `h` | absolute `h = 0.01·T₁`, `ξ = 0.01` | mode ∈ {1, 2, 3, 5, 10} | `viga('C')` |
| **D**, modal IC, swept `h` | `ξ = 0.01` | mode ∈ {1, 3, 10}, `h/T_mode ∈ {0.01; 0.1; 0.5; 0.9}` | `viga('D')` |
| **E**, alternative Rayleigh | impulse IC | `a₀ = 0`, `a₁ ∈ {10⁻⁴; 10⁻⁵}` | `viga('E')` |
| **F**, modal reduction | impulse IC, `ξ = 0.2` | full model, 4 modes, 3 modes | `viga('F')` |

The reasoning behind each choice:

- **B:** the same `h_RK`, sized from the most damped case, is used for every `ξ`. That
  is what shows RK1 and RK2 diverging with little damping and stable with a lot of
  damping at the same step.
- **C:** the step is absolute, not relative to the mode. The point is to let `h/T_mode`
  grow as the mode gets higher.
- **D:** PE, AD and the spectral radius are closed-form for every combination. The
  time-domain cycles are counted in `T₁`, so for mode 10 one cycle holds about 250
  oscillations. High ratios may diverge, which is part of the result, and each run is
  protected by `try/catch`.
- **E:** with `a₀ = 0`, `ξ_k = a₁ω_k/2` grows linearly with frequency. The `ξ` of each
  mode is computed, not imposed.
- **F:** Part 1 repeats the analysis for the three models. Part 2 (item 32) uses RK4 to
  separate the effect of a larger `h` from the effect of the reduction: full model with
  `h₁`, reduced model with the same `h₁`, and reduced model with its own larger `h₃`. All
  three are compared with the exact response of the full model.

### Decisions on open points

1. The damping set keeps `ξ = 0.2`, not `0.1`, for direct comparison with Study A.
2. Both the 4-mode and 3-mode reductions are run. The 4-mode one is used in item 32,
   because with 3 modes `ξ₃` is already close to 1 when anchoring on modes 1 and 3.
3. Both alternative Rayleigh cases are run, since the cost is low.
4. The Study C step is absolute, as explained above.
5. Study D stays on modes {1, 3, 10}, to contain the Runge-Kutta cost.

### Open items

- The Bathe reference for item 8 (`h ≲ T_min/10`) still needs the exact page.
- It is still undecided whether Study F becomes its own section of the paper or joins
  Section 4.4.

## Where to change parameters

| Location | Contains | Effect |
|---|---|---|
| `params` in `modelo.m` | material, geometry, mesh, target `ξ`, anchor modes, `V₀` | changes the physical problem |
| `config_estudos` in `modelo.m` | step, cycles, plot window, figures, sampling, cases `cfg.B` to `cfg.F` | changes the numerical experiment |
| start of `estudo_A` in `viga.m` | Study A cycles and window | changes Study A only |

Study A has its own configuration on purpose: it is the fixed point of comparison for
the others. Every field of `config_estudos` is commented in the code, with the effect
and cost of changing it.

## Figure time window

Period and amplitude errors accumulate, so the simulation runs many cycles and the
figure shows only a stretch of it. `janela_plot` picks that stretch using the exact
solution, since an unstable scheme would drag the window to the wrong place.

| Regime | Condition | What is shown |
|---|---|---|
| Strong damping | decays within `ciclosDecaimentoRapido` cycles | the whole response, from `t = 0` to the decay |
| Weak damping | decays later | `X_ciclos_final_plot` cycles before the decay, or the end of the run if it never decays |
| No damping | `ξ = 0` | the last `X_ciclos_final_plot_sem_xi` cycles |

The number of cycles to decay to `tol` is `n = −ln(tol)/(2πξ)`. With `tol = 10⁻²`,
`ξ = 0.2` takes 3.7 cycles and `ξ = 0.01` takes 73.

Each figure is normalized by the peak of the exact response inside the window, and the
y axis is clipped at `±1.6`, so an unstable curve leaves the frame instead of flattening
the others. The relative error uses the global peak and the whole simulation.

## Output

Everything goes to `ImgViga/`, always at the project root, whatever MATLAB's current
folder is. Each figure is saved as `.eps` (vector, for the paper) and `.png` (300 dpi).

| Location | Content |
|---|---|
| `ImgViga/viga_rho_{xi0,xi02}` | spectral radius of the nine schemes in Study A |
| `ImgViga/viga_desloc_{direto,rk}_{xi0,xi02}` | tip displacement, with a one-period inset |
| `ImgViga/viga_desloc_..._zoom` | the same period as its own figure |
| `ImgViga/viga_resultados.mat` | all Study A data |
| `ImgViga/EstudoB/` to `EstudoF/` | figures and `.mat` of each study |
| `ImgViga/Tabelas/` | four CSV files and `PARAMETROS_ESTUDOS.md` |

### Study A saved data

```matlab
S = load('ImgViga/viga_resultados.mat');
res = S.res;
```

| Field | Content |
|---|---|
| `res.model` | `K`, `M`, `C`, `w`, `Phi`, damping, ICs, `T1`, `wmax` |
| `res.config` | run parameters |
| `res.stability.hSweep` | steps of the spectral radius sweep |
| `res.stability.rho.<case>.<method>` | spectral radius curves |
| `res.stability.hcr_xi0`, `.hcr_xi02`, `.hcr_table` | critical steps |
| `res.pead.PE.<case>`, `.AD.<case>`, `.rho.<case>` | `[nMethods × nModes]` matrices |
| `res.pead.table.<case>` | table of the first three modes |
| `res.time.<case>.<method>` | `t`, `h`, `u_tip`, `v_tip`, `a_tip`, `err` |
| `res.time.<case>.exact` | exact reference on the same grid |
| `res.timeRK.<case>.<method>` | Runge-Kutta runs, decimated |
| `res.window.<case>` | plot window diagnostics |

`<case>` is `xi0` or `xi02`. `<method>` is `DC`, `Houbolt`, `Newmark`, `Wilson`, `Bathe`,
`RK1`, `RK2`, `RK3` or `RK4`.

Only the tip DOF is stored by default (`saveFullField = true` stores every DOF of the
implicit runs). RK histories are decimated to `nSaveRKMax` points, and `.u_tip_max` keeps
the true peak of the run. The PE and AD table exports directly:

```matlab
writetable(res.pead.table.xi02, 'pead_xi02.csv')
```

## Expected results

### Critical steps

| Scheme | `h_cr` (ξ=0) | `h_cr/T₁` | `h_cr` (ξ=0.2) | `h_cr/T₁` |
|---|---|---|---|---|
| RK1 | 1.58e-09 * | 1.32e-08 * | 6.0970e-07 | 5.09e-06 |
| RK2 | 1.88e-07 * | 1.57e-06 * | 6.0970e-07 | 5.09e-06 |
| RK3 | 1.9374e-05 | 1.62e-04 | 7.6601e-07 | 6.40e-06 |
| RK4 | 3.1637e-05 | 2.64e-04 | 8.4910e-07 | 7.09e-06 |
| DC | 2.2371e-05 | 1.87e-04 | 2.2371e-05 | 1.87e-04 |

The undamped limits reproduce the closed forms `DC = 2/ω_max`, `RK3 = √3/ω_max` and
`RK4 = √8/ω_max`.

\* Bisection tolerance artifacts. For undamped RK1 and RK2, `ρ > 1` for every `h > 0`, so
the correct statement is that no stable step exists.

### PE and AD with ξ = 0 (percent)

Implicit schemes at `h/T₁ = 0.01`.

| Scheme | PE m1 | AD m1 | PE m2 | AD m2 | PE m3 | AD m3 |
|---|---|---|---|---|---|---|
| Central Differences | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Houbolt | 0.0000 | 0.0767 | 0.0182 | 11.2851 | 0.6096 | 50.0852 |
| Newmark | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Wilson-θ | 0.0000 | 0.0130 | 0.0009 | 2.5944 | 0.0751 | 21.6136 |
| Bathe | 0.0000 | 0.0005 | 0.0000 | 0.1296 | 0.0007 | 2.3772 |
| RK1 | 0.0000 | −0.0091 | 0.0000 | −0.0567 | 0.0000 | −0.1590 |
| RK2 to RK4 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

### PE and AD with ξ = 0.2 (ratios, 1 reproduces the physical value)

| Scheme | PE m1 | AD m1 | PE m2 | AD m2 | PE m3 | AD m3 |
|---|---|---|---|---|---|---|
| Central Differences | 0.9998 | 1.0002 | 0.9932 | 1.0039 | 0.9337 | 1.0673 |
| Houbolt | 1.0014 | 1.0063 | 1.0643 | 0.9645 | 1.3576 | 0.7493 |
| Newmark | 1.0003 | 1.0008 | 1.0123 | 1.0156 | 1.0821 | 1.2312 |
| Wilson-θ | 1.0005 | 1.0026 | 1.0271 | 1.0192 | 1.1764 | 1.1619 |
| Bathe | 1.0001 | 1.0004 | 1.0063 | 1.0066 | 1.0444 | 1.0952 |
| RK1 | 1.0000 | 1.0001 | 1.0000 | 1.0006 | 0.9999 | 1.0016 |
| RK2 to RK4 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |

## What the example shows

1. **Physical damping makes RK1 and RK2 usable.** Undamped, neither has a stable step.
   With `ξ = 0.2` a finite `h_cr` appears, and at the same step the curves leave the frame
   at `ξ = 0` and follow the exact response at `ξ = 0.2`.
2. **Physical damping penalizes RK3 and RK4.** Their stable range shrinks by more than an
   order of magnitude.
3. **Central Differences does not feel the damping.** `h_cr` is unchanged to five digits.
4. **The limit is set by `ω_max`, not `ω₁`.** Conditionally stable schemes need `h/T₁`
   between `10⁻⁴` and `10⁻⁶` on the mesh, far below what the SDOF analysis suggests.
5. **The PE and AD ordering of Sections 4.2 and 4.3 repeats.** At undamped mode 3, AD goes
   Houbolt (50.1%) > Wilson-θ (21.6%) > Bathe (2.4%) > Newmark = DC (0%).
6. **`AD > 1` appears under damping for Newmark, Bathe and Wilson-θ:** amplitude grows
   relative to the physical decay, even though the scheme is stable. Houbolt gives `AD < 1`.
7. **RK dispersion and dissipation errors are negligible** at the steps stability imposes.
   Their limitation here is step size, not accuracy.

## Code map

| File | Main functions |
|---|---|
| `viga.m` | `viga`, `validacao`, `estudo_A` to `estudo_F`, `estudos_BF` |
| `modelo.m` | `params`, `config_estudos`, `matrizes_EF`, `autovalores_EF`, `amortecimento`, `ci_impulso`, `ci_modal`, `monta`, `modos_analiticos`, `solucao_exata`, `reducao_modal` |
| `metodos.m` | `mdof_dc`, `mdof_houbolt`, `mdof_newmark`, `mdof_wilson`, `mdof_bathe`, `mdof_rk`, `rk_stream_tip`, `integra_metodo`, `operador`, `rho_beam`, `pe_ad`, `pead_tabela`, `h_critico` |
| `graficos.m` | `espectro`, `cores`, `resposta_tempo`, `janela_plot`, `salva_figura`, `pasta_saida` |
| `resultados.m` | `tabelas`, `relatorio`, `teste_janela`, `testes` |

Code comments are in Portuguese. Figure titles and legends are in English, as in the paper.
