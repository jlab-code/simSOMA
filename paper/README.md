# Scripts for the figures of the simSOMA paper

Run from the repository root after `pip install .`. Main-text Fig. 6C–I and Supplementary Figs. S1–S4 of
the paper were produced with v0.2 (simulation code identical to v0.2.3); the scripts here are adapted to the
package layout of v0.3.0, which draws different random numbers for a given seed (seed fix, `CHANGELOG.md`).
The v0.2.3 versions are Supplementary file 3 of the paper.

| Script | Output |
|---|---|
| `fig6/make_fig6_data.py` | data of Fig. 6C–I (topology 01, three independently simulated layers, observation levels): `python3 paper/fig6/make_fig6_data.py --reps 10 --out paper/fig6/data` |
| `fig6/fig6_settings.json` | settings record of the data shown in Fig. 6 |
| `fig6/fig6_stats.py` | numbers quoted in the main text: `python3 paper/fig6/fig6_stats.py paper/fig6/data` |
| `fig6/plot_fig6_CI_styled.py` | panels C–I: `python3 paper/fig6/plot_fig6_CI_styled.py paper/fig6/data Fig6_CI.pdf` |
| `fig6/plot_fig6_draft.py` | draft layout of the same panels |
| `si_schematics/si_schematics.py` | Supplementary Figs. S1, S2, S4: `python3 paper/si_schematics/si_schematics.py OUTDIR` |
| `si_organ_tail/organ_tail.py` | Supplementary Fig. S3 (organ 1/f tail): `python3 paper/si_organ_tail/organ_tail.py FigureS3.pdf` |
| `si_organ_tail/organ_multiplier_counts.py` | numbers for the organ mutation-rate multiplier (Supplementary Text 7.6) |
