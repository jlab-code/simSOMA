# Notation: paper symbols, configuration keys and output fields

From v0.3.0 the documentation uses the symbols of the simSOMA paper (main text and Supplementary Text,
revised version). Configuration keys and output fields keep their names; this table maps them. Fields
that were renamed are listed at the end with their deprecated aliases.

Running indices (i, j, k, o, t) are defined locally in each formula: in the observation model, i is a
mutation or site, o an organ and k a SAM layer.

## Developmental parameters

| Symbol | Meaning | Configuration key |
|---|---|---|
| m | stem-cell (ASC) positions in the apical niche | `m` |
| ρ | probability of one displacement per self-renewal round | `rho` |
| λ_vic | victim locality (0 uniform, 1 nearest neighbour) | `victim_locality` |
| w_B | branch-local competition bias of the favoured clone | `bias_mode`, `branch_bias_value`, `branch_bias_mean`, `branch_bias_kappa` |
| μ_unit | mutations per lineage per topology unit | `mu_unit` (alias `mu_year`) |
| κ_sr | self-renewal rounds per topology unit | `topology.mapping_rate` (run output `kappa_sr`) |
| μ_div | mutations per lineage per division, μ_unit/κ_sr | `mu_div` (optional override) |
| C | cells of the sectorized SAM-boundary population | `sam_boundary_cells` |
| P_b, P_b^eff | branch precursors; effective number min{P_b, m, C} | `branch_precursor_number`; output `P_b_eff` (alias `P_a_eff`) |
| P_o | organ precursors | `organ_precursor_number` |
| O | terminal organ cells | `organ_total_cells` |
| n_seq | sampled organ cells | `sequenced_cells`, `seq_fraction` |
| g | organ mutation-rate multiplier | `organ_mu_multiplier` |
| θ_div | phyllotactic divergence angle | `phyllotaxy.divergence_angle_deg` |
| ε_tie | tolerance for grouping events into one phyllotactic site | `phyllotaxy.tie_tol` |

## Topology, events and cell states

| Symbol | Meaning | Field |
|---|---|---|
| B, pa(B), B_root | branch, its parent, root branch | `branch_id`, `parent_id`, `root_id` |
| L_B, T_B | branch length (topology units); duration in self-renewal rounds | `length`; `T` (observed `T_obs`) |
| ι = (t_ι, τ_ι, z_ι) | event: branch-local coordinate, type, target | `time`, `type`, `target_id` |
| pos_ι | fractional position of an event on its branch | `pos` |
| θ_{B,0} | initial phyllotactic phase of branch B | `phyllo_phase_rad`, `phyllo_phase_deg` |
| j, θ_{B,j} | site rank and site angle | `phyllo_site_rank`, `phyllo_angle_rad`, `phyllo_angle_deg` |
| mut(c) | mutation identifiers carried by cell c | `CellState.genotype` |
| root(c) | root-lineage label (diagnostic) | `CellState.lineage_id`; outputs `root_lineage_*` |
| fav(c) | favoured-clone label | `CellState.comp_label` |
| i_B^* | favoured seed position of branch B | `favored_seed_index` |
| t_B^0 | branch-age offset of a branch segment | `age_offset` |

## Founding

| Symbol | Meaning | Field |
|---|---|---|
| k_foc | focal index of a founding event on the SAM-boundary ring | `focal_index` |
| 𝓘(k_foc, P) | positions of the founder block | `sampled_indices` |
| m_j, O_j | quotas: child-niche positions (branch), terminal cells (organ) per founder | `precursor_quotas` |
| b_i, σ_i | first position and size of sector i | `sector_bounds` (pre-branching output) |
| K | sectors spanned by the founder block | `founder_sector_count` |
| 1(K > 1) | polyclonal founding | `founder_polyclonal` |
| N_eff | effective number of sectors | `founder_effective_sectors` |
| (N_eff − 1)/(m − 1) | normalized founder diversity (per event) | `founder_diversity` |
| π_B, π_O; d_B, d_O | run means of `founder_polyclonal`; of `founder_diversity` | event and truth tables |
| χ_B, χ_O | founder quantity (P_eff − 1)m/C | `chi_B`, `chi_O` (deprecated alias `phi_B`, `phi_O`) |
| π | expected probability of polyclonal founding, min(1, χ) | `pi_B_expected`, `pi_O_expected` |
| E[K] | expected number of sectors | `sector_count_B_expected`, `sector_count_O_expected` |

## Organ readout and summaries

| Symbol | Meaning | Field |
|---|---|---|
| k_io | sampled cells of organ o carrying mutation i | `allele_counts_by_mutation` (`{i: k_io}`), `allele_count` |
| n_o | sampled cells of organ o | `n_sampled_cells` |
| v_io | developmental VAF k_io/n_o | `sampled_vaf`, `source_vaf` |
| deg(i) | number of sampled organs carrying mutation i (sharing degree) | `sharing_degree` |
| N_organs | number of sampled organs | `n_organs` |
| f | carrier fraction of a mutation; organ tail M(≥ f) = 2gμ_div(1/f − P_o) | `carrier_fraction` |
| v_min | low-VAF threshold of the formula statistics (default 0.05) | `--low-threshold` |

## Observation model (plantsoma_obs)

| Symbol | Meaning | Field |
|---|---|---|
| a_k, a_ko | contribution of layer k (in organ o) to the sample | `layer_contributions`, `effective_layer_contribution` |
| η | phasing factor (1 phased, 1/2 unphased) | `phase`, `phase_factor` |
| ṽ_io | assay (expected observed) VAF | `assay_vaf`, `observed_vaf` (exact transform) |
| d̄ | mean read depth | `depth.mean` / `depth.value` |
| d_io | read depth | `depth`, `DP` |
| u_i, u_io, u_o | site, site-by-organ and organ depth factors | `site_sdlog`, `sample_sdlog`, `sample_factor_sdlog` |
| y_io | variant reads | `alt_count`, `alt_depth` |
| e | sequencing error rate; p_e(v) = e + (1 − 2e)v | `sequencing_error` |
| k_c, ψ | beta-binomial concentration; intra-class correlation 1/(k_c + 1) | `concentration` (`read_concentration`) |
| q | background artefact VAF | `background.mean`, `background.distribution` |

## Renamed in v0.3.0

| Old (≤ v0.2.3) | New | Note |
|---|---|---|
| ϕ = (P_eff − 1)m/C | χ (χ_B, χ_O) | output keys `chi_B`, `chi_O`; `phi_B`, `phi_O` kept as aliases with identical values |
| module imports `import pipeline_wrapper` | `from simsoma import pipeline_wrapper` | package layout |
