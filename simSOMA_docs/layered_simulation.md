# Per-layer simulation (`simsoma layers`)

`simsoma layers CONFIG` simulates the meristem layers (e.g. L1, L2, L3) as independent cell-lineage
histories on the same plant topology, mixes them in each sampled organ, and draws sequencing reads
with the shared observation model (`plantsoma_obs`). Start from the template:

```bash
simsoma template layered > my_layers.json
```

## Model

Each layer is one simSOMA run with its own random seed and its own mutation rate (`mu_unit`); the
developmental parameters (m, rho, C, P_b, P_o, O, ...) and the topology are shared. There is no cell
exchange between layers, so a mutation arises in, and can become fixed in, one layer only. For a
mutation of layer k with carrier-cell fraction f in organ o, the assay VAF of a bulk sample is

    v = eta * c_ko * f,     eta = 1/2 (unphased diploid) or 1 (phased),

where c_ko is the contribution of layer k to the sample of organ o (sum over layers = 1). Read counts
are drawn from this VAF with the depth, read, background and caller settings of `observation.model`.

## Layer mixtures: `observation.layer_contributions`

One mixture for all organs:

```json
"layer_contributions": {"L1": 0.3, "L2": 0.6, "L3": 0.1}
```

One mixture per organ, with an optional `default` for organs not listed:

```json
"layer_contributions": {
  "leaf_A":  {"L1": 0.2, "L2": 0.7, "L3": 0.1},
  "fruit_1": {"L1": 0.1, "L2": 0.3, "L3": 0.6},
  "default": {"L1": 0.3, "L2": 0.6, "L3": 0.1}
}
```

Each mixture must be non-negative and sum to 1, and may only name simulated layers (layers that are
not named get contribution 0). Without a `default`, every sampled organ must be listed. Layer-specific
sampling of one layer is the special case with contribution 1 for that layer.

The layers themselves are defined in `simulation.layers`, each with its own `mu_unit` (mutations per
cell lineage per topology unit). Any set of layers can be used, e.g. only L1 and L2.

## Template defaults (leaves) and their sources

The template values are defaults for **leaves**, derived from the literature below. They are a
starting point, not estimates for a particular species or tissue: adapt them to your system.

| Setting | Template value | Basis |
|---|---|---|
| `layer_contributions` | L1 0.13, L2 0.84, L3 0.03 | mean of the three quantitative leaf compositions below, each rescaled to L1 + L2 + L3 = 1 |
| `simulation.layers.*.mu_unit` | L1 0.95, L2 0.55, L3 0.55 | illustrative relative values: L1 accumulates more mutations than the inner layers (Goel et al. 2024; Amundson et al. 2025) |

Quantitative leaf compositions (L1 = epidermis, L2 = mesophyll, L3 = vasculature):

| Reference | System | Reported | Rescaled L1 / L2 / L3 |
|---|---|---|---|
| Pyke, Marrison & Leech 1991 | *Arabidopsis thaliana*, first leaf | tissue volume mesophyll : airspace : epidermis : vasculature = 61 : 26 : 12 : 1 | 0.162 / 0.824 / 0.014 |
| Tolleter et al. 2024 | *A. thaliana*, mature leaf atlas | photosynthetic cells 86% and epidermis 14% of cellular volume; veins 0.7% of leaf volume | 0.139 / 0.854 / 0.007 |
| Goel et al. 2024 | apricot leaves | leaf scRNA-seq: mesophyll 78.2%, epidermis 8.9%, vasculature 7.2% of cells | 0.094 / 0.829 / 0.076 |
| **Mean** | | | **0.132 / 0.836 / 0.032** |

Further observations on layers and VAFs:

| Reference | System | Finding |
|---|---|---|
| Goel et al. 2024 | apricot fruits and leaves | >90% of somatic mutations layer-specific; higher mutation load in L1 |
| Ren et al. 2021 | *Salix suchowensis* leaves | most leaf somatic mutations at VAF < 0.3, consistent with leaves founded by several cell lineages |
| Amundson et al. 2025 | potato (tetraploid), whole leaf and layer-enriched fractions | L1-like mutations near VAF 0.06, L2/L3-like near 0.20; higher L1 mutation rate |

Notes. The compositions are tissue volumes or cell counts, not shares of sequenced DNA; differences
in nuclear DNA content between cell types (e.g. endoreduplication) can shift the DNA share. The three
studies give an L2 share of 0.82-0.85; allowing for species differences, a conservative range is
about 0.60-0.85. A fixed heterozygous L2 mutation in unphased bulk leaf DNA is then expected at VAF
0.5 x c_L2, i.e. about 0.30-0.43 (0.42 with the template value). In woody species (willow), leaf
somatic mutations are mostly found below VAF 0.3.

## Outputs

Per replicate (`<output.dir>/replicate_NNNN/`):

| File | Content |
|---|---|
| `layer_carriers.csv.gz` | true carrier fraction of every mutation in every organ, per layer |
| `read_evidence.csv.gz` | depth, alternative reads, observed VAF, called / ascertained, per site and organ |
| `vafsoma_dp_*_vaf.csv` | read tables in vafSOMA input format, one per depth tier |

`layered_settings.json` in the output folder records the settings used (including the seed) and the software versions.
