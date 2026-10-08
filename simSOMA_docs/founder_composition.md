# Clonal composition of branch and organ founders (definition `local_sector_v1`)

Applies from simSOMA 0.2.0-dev (branch `revision/code-fixes`). Draft text for the SI.

## Definition
At each developmental event the current ASC niche (m ordered positions) is amplified into the
SAM-boundary ring of C cells. Niche position i seeds one contiguous **clonal sector** S_i
(equalized quotas; sizes differ by at most one cell). A founding event samples a contiguous block
of P_eff precursor cells: P_b,eff = min(P_b, m, C) for a branch, P_o,eff = min(P_o, C, O) for an organ.

For one event, with n_i founders from sector S_i and p_i = n_i / P_eff:

| field | symbol | definition |
|---|---|---|
| `founder_sector_count` | K | number of sectors with n_i > 0 |
| `founder_polyclonal` | 1(K > 1) | 1 = polyclonal, 0 = monoclonal founding |
| `founder_effective_sectors` | N_eff | 1 / sum_i p_i^2 |
| `founder_diversity` | d | (N_eff - 1) / (m - 1), 0 if m = 1 |
| `founder_dominant_fraction` | | max_i p_i |
| `founder_sector_counts` | | {sector index: n_i}; `founder_lineage_counts` is the same (legacy key) |

Run-level quantities (fitSOMA targets) are means over events of one type:
pi_B, pi_O = mean `founder_polyclonal`; d_B, d_O = mean `founder_diversity`.

"Sector" is **local**: it refers to the niche positions at the time of the event.

## Expected values
Under uniform focal placement (no phyllotaxy) the expectations depend only on (P_eff, C, m).
For C divisible by m and P_eff <= C/m + 1:

    Pr(polyclonal) = min(1, phi),   E[K] = 1 + Pr(polyclonal),   phi = (P_eff - 1) m / C.

`founder_diversity.polyclonal_founding_expectation(P_eff, C, m)` gives exact values for any
(P_eff, C, m). Configured truth exports `phi_B`, `phi_O`, `pi_B_expected`, `pi_O_expected`,
`sector_count_B_expected`, `sector_count_O_expected`.

Note: fitSOMA 0.3.19 computes `phi_B = P_a_eff * m / C` (no -1). The geometric quantity is
(P_eff - 1) m / C; the two should be aligned.

## Detectability
Polyclonal founding leaves a genetic signature only if the sectors involved differ genetically.
Under high turnover (rho) or shortly after a monoclonal branch founding, neighbouring sectors are
recent relatives and polyclonal founding is nearly invisible in VAF data. Recoverability of
pi_B, pi_O must therefore be demonstrated (fitSOMA recovery benchmark), and is expected to be
regime dependent (weaker for branches than for organs).

## Previous definition (deprecated)
simSOMA <= 0.1.x counted distinct **root** lineage labels among founders. These labels are set
once in the root niche and never reset, so the statistic measured descent from different
root-niche cells. It collapses with turnover and branch order (e.g. m = 4, C = 64, P_o = 8:
pi_O = 0.44 at rho = 0, 0.06 at rho = 0.05, 0 at rho = 0.5). It is retained as
`root_lineage_*` fields for comparison only. No manuscript figure used it.
