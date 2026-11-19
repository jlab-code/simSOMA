"""summaries.py

Dedicated result summarization for the pipeline.

This module is intentionally downstream of the simulator modules: it consumes
raw per-organ mutation allele counts and derives within-organ VAF spectra,
private-versus-shared spectra across organs, and related summary statistics.
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple
import math
import statistics


def _vaf_bin_index(vaf: float, nbins: int) -> int:
    if nbins <= 0:
        raise ValueError("nbins must be positive")
    if vaf >= 1.0:
        return nbins - 1
    if vaf <= 0.0:
        return 0
    i = int(math.floor(vaf * nbins))
    if i < 0:
        return 0
    if i >= nbins:
        return nbins - 1
    return i


def allele_count_spectrum_from_mutation_counts(allele_counts_by_mutation: Dict[int, int], sequenced_cells: int) -> Dict[int, int]:
    counts_by_k: Dict[int, int] = {}
    for _, c in allele_counts_by_mutation.items():
        cc = int(c)
        if cc < 1 or cc > int(sequenced_cells):
            raise ValueError(f"allele count {cc} outside 1..sequenced_cells={sequenced_cells}")
        counts_by_k[cc] = counts_by_k.get(cc, 0) + 1
    return counts_by_k


def _bin_spectrum_counts_by_k(counts_by_k: Dict[int, int], sequenced_cells: int, nbins: int) -> List[int]:
    bins = [0 for _ in range(nbins)]
    for k, c in counts_by_k.items():
        kk = int(k)
        cc = int(c)
        if cc <= 0:
            continue
        vaf = float(kk) / float(sequenced_cells)
        idx = _vaf_bin_index(vaf, nbins)
        bins[idx] += cc
    return bins


def _add_count_spectrum_row(
    rows: List[Dict[str, Any]],
    *,
    level: str,
    organ_id: str,
    variant_class: str,
    sharing_degree: Any,
    sequenced_cells: int,
    allele_count: int,
    n_variants: int,
) -> None:
    if n_variants <= 0:
        return
    sc = int(sequenced_cells)
    ac = int(allele_count)
    if sc <= 0 or ac <= 0:
        return
    rows.append({
        "level": str(level),
        "organ_id": str(organ_id) if organ_id is not None else "",
        "variant_class": str(variant_class),
        "sharing_degree": sharing_degree if sharing_degree is not None else "",
        "n_sampled_cells": sc,
        "allele_count": ac,
        "sampled_vaf": float(ac) / float(sc),
        "n_variants": int(n_variants),
    })


def _aggregate_count_spectrum_rows(summaries: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Aggregate replicate-level allele-count spectra with zero-filled missing keys.

    Each replicate summary stores compact rows, not mutation IDs. When aggregating
    across replicates, absent rows must count as zero so that means and population
    SDs are computed correctly after later re-binning.
    """
    replicate_maps: List[Dict[Tuple[str, str, str, str, int, int], int]] = []
    all_keys = set()
    for s in summaries:
        rows = s.get("vaf_count_spectrum", {}).get("rows", [])
        cur: Dict[Tuple[str, str, str, str, int, int], int] = {}
        for r in rows:
            key = (
                str(r.get("level", "")),
                str(r.get("organ_id", "")),
                str(r.get("variant_class", "")),
                str(r.get("sharing_degree", "")),
                int(r.get("n_sampled_cells", 0)),
                int(r.get("allele_count", 0)),
            )
            cur[key] = cur.get(key, 0) + int(r.get("n_variants", 0))
            all_keys.add(key)
        replicate_maps.append(cur)

    out_rows: List[Dict[str, Any]] = []
    n_reps = len(replicate_maps)
    for key in sorted(all_keys, key=lambda x: (x[0], x[1], x[2], str(x[3]), x[4], x[5])):
        vals = [m.get(key, 0) for m in replicate_maps]
        mean = statistics.mean(vals) if vals else 0.0
        sd = statistics.pstdev(vals) if vals else 0.0
        level, organ_id, variant_class, sharing_degree, n_sampled_cells, allele_count = key
        out_rows.append({
            "level": level,
            "organ_id": organ_id,
            "variant_class": variant_class,
            "sharing_degree": sharing_degree,
            "n_sampled_cells": int(n_sampled_cells),
            "allele_count": int(allele_count),
            "sampled_vaf": float(allele_count) / float(n_sampled_cells) if int(n_sampled_cells) > 0 else None,
            "n_variants_mean": mean,
            "n_variants_sd": sd,
        })
    return {
        "available": True,
        "n_replicates": int(n_reps),
        "rows": out_rows,
    }



class _RunningStats:
    """Population mean/SD accumulator using Welford updates."""

    __slots__ = ("n", "mean", "M2")

    def __init__(self) -> None:
        self.n = 0
        self.mean = 0.0
        self.M2 = 0.0

    def update(self, value: Any) -> None:
        if value is None:
            return
        x = float(value)
        if not math.isfinite(x):
            return
        self.n += 1
        delta = x - self.mean
        self.mean += delta / float(self.n)
        delta2 = x - self.mean
        self.M2 += delta * delta2

    def result(self) -> Dict[str, Any]:
        if self.n <= 0:
            return {"mean": None, "sd": None}
        return {"mean": self.mean, "sd": math.sqrt(self.M2 / float(self.n))}


class _RunningListStats:
    """Fixed-width vector of population mean/SD accumulators."""

    def __init__(self, width: int) -> None:
        self.width = int(width)
        self.stats = [_RunningStats() for _ in range(self.width)]

    def update(self, values: Any) -> None:
        vals = list(values or [])
        if len(vals) < self.width:
            vals = vals + [0] * (self.width - len(vals))
        for i in range(self.width):
            self.stats[i].update(vals[i])

    def result(self) -> Dict[str, Any]:
        return {
            "mean": [st.result()["mean"] for st in self.stats],
            "sd": [st.result()["sd"] for st in self.stats],
        }


class _RunningSparseCountStats:
    """Sparse count accumulator with zero-filled missing keys per replicate.

    The object preserves the semantics used by aggregate_replicate_summaries:
    if a key is absent from a replicate, that replicate contributes zero for the
    key. This is essential for correct population SDs of sparse spectra.
    """

    def __init__(self) -> None:
        self.n_updates = 0
        self.stats: Dict[Any, _RunningStats] = {}

    def update(self, values: Dict[Any, Any]) -> None:
        vals = dict(values or {})
        existing = set(self.stats.keys())
        current = set(vals.keys())
        for key in current.difference(existing):
            st = _RunningStats()
            for _ in range(self.n_updates):
                st.update(0)
            self.stats[key] = st
        for key, st in self.stats.items():
            st.update(vals.get(key, 0))
        self.n_updates += 1

    def result_by_key(self) -> Dict[Any, Dict[str, Any]]:
        return {key: st.result() for key, st in self.stats.items()}


class SummaryAccumulator:
    """Bounded-memory aggregator for replicate summaries.

    This is the streaming equivalent of aggregate_replicate_summaries(). It keeps
    only online statistics and sparse-key accumulators, not the list of all
    replicate summaries. It is intended for large grid runs where n_sim can be
    high.
    """

    def __init__(self) -> None:
        self.n = 0
        self.nbins: int | None = None
        self.organ_ids: List[str] = []
        self.per_organ: Dict[str, Dict[str, Any]] = {}
        self.aggregate: Dict[str, Any] = {
            "total_variants": _RunningStats(),
            "total_unique_variants": _RunningStats(),
            "total_variant_occurrences": _RunningStats(),
            "vaf_bins_organ_occurrences": None,
        }
        self.vaf_count_spectrum = _RunningSparseCountStats()
        self.private_shared_available: bool | None = None
        self.private_shared_reason: Any = None
        self.private_shared: Dict[str, Any] = {}

    def _ensure_initialized(self, summary: Dict[str, Any]) -> None:
        if self.nbins is not None:
            return
        self.nbins = int(summary.get("nbins", 20))
        self.organ_ids = [str(x) for x in summary.get("organ_ids", [])]
        self.aggregate["vaf_bins_organ_occurrences"] = _RunningListStats(self.nbins)
        scalar_diag_fields = [
            "n_lineages_sampled", "dominant_lineage_id", "dominant_fraction", "lineage_entropy",
            "boundary_cells", "organ_total_cells", "organ_precursor_number",
        ]
        for oid in self.organ_ids:
            self.per_organ[oid] = {
                "sequenced_cells": None,
                "counts_by_allele_count": _RunningSparseCountStats(),
                "vaf_bins": _RunningListStats(self.nbins),
                "total_variants": _RunningStats(),
                "singletons": _RunningStats(),
                "fixed": _RunningStats(),
            }
            for fld in scalar_diag_fields:
                self.per_organ[oid][fld] = _RunningStats()

    def update(self, summary: Dict[str, Any]) -> None:
        if not isinstance(summary, dict):
            raise TypeError("SummaryAccumulator.update requires a summary dict.")
        self._ensure_initialized(summary)
        assert self.nbins is not None

        per_organ_summary = summary.get("per_organ", {})
        for oid in self.organ_ids:
            os = per_organ_summary.get(oid, {})
            po = self.per_organ[oid]
            if po["sequenced_cells"] is None and os.get("sequenced_cells") is not None:
                po["sequenced_cells"] = int(os.get("sequenced_cells"))
            counts = {int(k): int(v) for k, v in os.get("counts_by_allele_count", {}).items()}
            po["counts_by_allele_count"].update(counts)
            po["vaf_bins"].update(list(map(int, os.get("vaf_bins", [0 for _ in range(self.nbins)]))))
            for fld in ["total_variants", "singletons", "fixed", "n_lineages_sampled", "dominant_lineage_id", "dominant_fraction", "lineage_entropy", "boundary_cells", "organ_total_cells", "organ_precursor_number"]:
                po[fld].update(os.get(fld))

        agg = summary.get("aggregate", {})
        self.aggregate["total_variants"].update(int(agg.get("total_variants", 0)))
        self.aggregate["total_unique_variants"].update(int(agg.get("total_unique_variants", agg.get("total_variants", 0))))
        self.aggregate["total_variant_occurrences"].update(int(agg.get("total_variant_occurrences", 0)))
        self.aggregate["vaf_bins_organ_occurrences"].update(list(map(int, agg.get("vaf_bins_organ_occurrences", agg.get("vaf_bins", [0 for _ in range(self.nbins)])))))

        spec = summary.get("vaf_count_spectrum", {}) if isinstance(summary.get("vaf_count_spectrum"), dict) else {}
        spec_values: Dict[Tuple[str, str, str, str, int, int], int] = {}
        for r in spec.get("rows", []):
            key = (
                str(r.get("level", "")),
                str(r.get("organ_id", "")),
                str(r.get("variant_class", "")),
                str(r.get("sharing_degree", "")),
                int(r.get("n_sampled_cells", 0)),
                int(r.get("allele_count", 0)),
            )
            spec_values[key] = spec_values.get(key, 0) + int(r.get("n_variants", 0))
        self.vaf_count_spectrum.update(spec_values)

        ps = summary.get("private_shared")
        ps_available = isinstance(ps, dict) and bool(ps.get("available"))
        if self.private_shared_available is None:
            self.private_shared_available = ps_available
            if not ps_available:
                self.private_shared_reason = ps.get("reason") if isinstance(ps, dict) else "requires at least 2 sampled organs"
                self.private_shared = {}
            else:
                self.private_shared = {
                    "n_shared_mutations": _RunningStats(),
                    "n_private_mutations_total": _RunningStats(),
                    "n_private_by_organ": {oid: _RunningStats() for oid in self.organ_ids},
                    "pairwise_shared_counts": _RunningSparseCountStats(),
                    "sharing_degree_counts": _RunningSparseCountStats(),
                    "private_vaf_bins_global": _RunningListStats(self.nbins),
                    "shared_vaf_bins_global": _RunningListStats(self.nbins),
                    "private_vaf_bins_by_organ": {oid: _RunningListStats(self.nbins) for oid in self.organ_ids},
                    "shared_vaf_bins_by_organ": {oid: _RunningListStats(self.nbins) for oid in self.organ_ids},
                }
        if self.private_shared_available:
            if not ps_available:
                raise ValueError("private_shared availability changed across replicates; cannot aggregate safely.")
            assert isinstance(ps, dict)
            self.private_shared["n_shared_mutations"].update(int(ps.get("n_shared_mutations", 0)))
            self.private_shared["n_private_mutations_total"].update(int(ps.get("n_private_mutations_total", 0)))
            for oid in self.organ_ids:
                self.private_shared["n_private_by_organ"][oid].update(int(ps.get("n_private_by_organ", {}).get(oid, 0)))
            self.private_shared["pairwise_shared_counts"].update({str(k): int(v) for k, v in ps.get("pairwise_shared_counts", {}).items()})
            self.private_shared["sharing_degree_counts"].update({str(k): int(v) for k, v in ps.get("sharing_degree_counts", {}).items()})
            self.private_shared["private_vaf_bins_global"].update(list(map(int, ps.get("private_vaf_bins_global", [0 for _ in range(self.nbins)]))))
            self.private_shared["shared_vaf_bins_global"].update(list(map(int, ps.get("shared_vaf_bins_global", [0 for _ in range(self.nbins)]))))
            for oid in self.organ_ids:
                self.private_shared["private_vaf_bins_by_organ"][oid].update(list(map(int, ps.get("private_vaf_bins_by_organ", {}).get(oid, [0 for _ in range(self.nbins)]))))
                self.private_shared["shared_vaf_bins_by_organ"][oid].update(list(map(int, ps.get("shared_vaf_bins_by_organ", {}).get(oid, [0 for _ in range(self.nbins)]))))
        self.n += 1

    def finalize(self) -> Dict[str, Any]:
        if self.n <= 0 or self.nbins is None:
            return {}
        out: Dict[str, Any] = {"nbins": int(self.nbins), "organ_ids": list(self.organ_ids)}
        per_organ_out: Dict[str, Any] = {}
        for oid in self.organ_ids:
            po = self.per_organ[oid]
            counts_stats = po["counts_by_allele_count"].result_by_key()
            keys = sorted(int(k) for k in counts_stats.keys())
            per_organ_out[oid] = {
                "sequenced_cells": int(po["sequenced_cells"] or 0),
                "counts_by_allele_count": {
                    "keys": keys,
                    "mean": [counts_stats[k]["mean"] for k in keys],
                    "sd": [counts_stats[k]["sd"] for k in keys],
                },
                "vaf_bins": po["vaf_bins"].result(),
                "total_variants": po["total_variants"].result(),
                "singletons": po["singletons"].result(),
                "fixed": po["fixed"].result(),
                "n_lineages_sampled": po["n_lineages_sampled"].result(),
                "dominant_lineage_id": po["dominant_lineage_id"].result(),
                "dominant_fraction": po["dominant_fraction"].result(),
                "lineage_entropy": po["lineage_entropy"].result(),
                "boundary_cells": po["boundary_cells"].result(),
                "organ_total_cells": po["organ_total_cells"].result(),
                "organ_precursor_number": po["organ_precursor_number"].result(),
            }
        out["per_organ"] = per_organ_out
        out["aggregate"] = {
            "total_variants": self.aggregate["total_variants"].result(),
            "total_unique_variants": self.aggregate["total_unique_variants"].result(),
            "total_variant_occurrences": self.aggregate["total_variant_occurrences"].result(),
            "vaf_bins_organ_occurrences": self.aggregate["vaf_bins_organ_occurrences"].result(),
        }
        spec_rows: List[Dict[str, Any]] = []
        for key, st in sorted(self.vaf_count_spectrum.result_by_key().items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], str(kv[0][3]), kv[0][4], kv[0][5])):
            level, organ_id, variant_class, sharing_degree, n_sampled_cells, allele_count = key
            spec_rows.append({
                "level": level,
                "organ_id": organ_id,
                "variant_class": variant_class,
                "sharing_degree": sharing_degree,
                "n_sampled_cells": int(n_sampled_cells),
                "allele_count": int(allele_count),
                "sampled_vaf": float(allele_count) / float(n_sampled_cells) if int(n_sampled_cells) > 0 else None,
                "n_variants_mean": st["mean"],
                "n_variants_sd": st["sd"],
            })
        out["vaf_count_spectrum"] = {"available": True, "n_replicates": int(self.n), "rows": spec_rows}
        if not self.private_shared_available:
            out["private_shared"] = {"available": False, "reason": self.private_shared_reason or "requires at least 2 sampled organs"}
            return out
        ps_out: Dict[str, Any] = {
            "available": True,
            "n_shared_mutations": self.private_shared["n_shared_mutations"].result(),
            "n_private_mutations_total": self.private_shared["n_private_mutations_total"].result(),
            "n_private_by_organ": {oid: self.private_shared["n_private_by_organ"][oid].result() for oid in self.organ_ids},
            "pairwise_shared_counts": self.private_shared["pairwise_shared_counts"].result_by_key(),
            "sharing_degree_counts": self.private_shared["sharing_degree_counts"].result_by_key(),
            "private_vaf_bins_global": self.private_shared["private_vaf_bins_global"].result(),
            "shared_vaf_bins_global": self.private_shared["shared_vaf_bins_global"].result(),
            "private_vaf_bins_by_organ": {oid: self.private_shared["private_vaf_bins_by_organ"][oid].result() for oid in self.organ_ids},
            "shared_vaf_bins_by_organ": {oid: self.private_shared["shared_vaf_bins_by_organ"][oid].result() for oid in self.organ_ids},
        }
        # Keep keys as strings, matching aggregate_replicate_summaries output.
        ps_out["pairwise_shared_counts"] = {str(k): v for k, v in ps_out["pairwise_shared_counts"].items()}
        ps_out["sharing_degree_counts"] = {str(k): v for k, v in ps_out["sharing_degree_counts"].items()}
        out["private_shared"] = ps_out
        return out

def organ_vaf_summary(organ_event: Dict[str, Any], *, nbins: int = 20) -> Dict[str, Any]:
    sequenced_cells = int(organ_event["sequenced_cells"])
    allele_counts_by_mutation = {int(k): int(v) for k, v in organ_event["allele_counts_by_mutation"].items()}
    counts_by_k = allele_count_spectrum_from_mutation_counts(allele_counts_by_mutation, sequenced_cells)
    total = int(sum(counts_by_k.values()))
    singletons = int(counts_by_k.get(1, 0))
    fixed = int(counts_by_k.get(sequenced_cells, 0))
    bins = _bin_spectrum_counts_by_k(counts_by_k, sequenced_cells, nbins)
    return {
        "sequenced_cells": sequenced_cells,
        "total_variants": total,
        "singletons": singletons,
        "fixed": fixed,
        "vaf_bins": bins,
        "counts_by_allele_count": counts_by_k,
        "n_lineages_sampled": int(organ_event.get("n_lineages_sampled", len(organ_event.get("lineage_counts", {})))),
        "dominant_lineage_id": int(organ_event.get("dominant_lineage_id", -1)),
        "dominant_fraction": float(organ_event.get("dominant_fraction", 0.0)),
        "lineage_entropy": float(organ_event.get("lineage_entropy", 0.0)),
        "boundary_cells": int(organ_event.get("boundary_cells", 0)),
        "organ_total_cells": int(organ_event.get("organ_total_cells", 0)),
        "organ_precursor_number": int(organ_event.get("precursor_number", 0)),
    }


def compute_summary(
    result: Dict[str, Any],
    *,
    nbins: int = 20,
    private_shared: bool = True,
) -> Dict[str, Any]:
    organ_events: List[Dict[str, Any]] = list(result.get("organ_events", []))
    organ_ids = [str(oe["organ_id"]) for oe in organ_events]

    per_organ: Dict[str, Any] = {}
    agg_bins_occ = [0 for _ in range(nbins)]
    mut_to_orgcounts: Dict[int, Dict[str, int]] = {}
    vaf_count_rows: List[Dict[str, Any]] = []
    global_all_counts: Dict[Tuple[int, int], int] = {}

    for oe in organ_events:
        oid = str(oe["organ_id"])
        organ_summary = organ_vaf_summary(oe, nbins=nbins)
        bins = organ_summary["vaf_bins"]
        for i in range(nbins):
            agg_bins_occ[i] += int(bins[i])
        per_organ[oid] = organ_summary
        sequenced_cells = int(organ_summary["sequenced_cells"])

        for allele_count, n_variants in organ_summary.get("counts_by_allele_count", {}).items():
            _add_count_spectrum_row(
                vaf_count_rows,
                level="organ",
                organ_id=oid,
                variant_class="all",
                sharing_degree="",
                sequenced_cells=sequenced_cells,
                allele_count=int(allele_count),
                n_variants=int(n_variants),
            )
            key = (sequenced_cells, int(allele_count))
            global_all_counts[key] = global_all_counts.get(key, 0) + int(n_variants)

        mc = oe.get("allele_counts_by_mutation")
        if isinstance(mc, dict):
            for mid, cnt in ((int(k), int(v)) for k, v in mc.items()):
                if cnt > 0:
                    mut_to_orgcounts.setdefault(mid, {})[oid] = int(cnt)

    for (sequenced_cells, allele_count), n_variants in sorted(global_all_counts.items()):
        _add_count_spectrum_row(
            vaf_count_rows,
            level="global",
            organ_id="",
            variant_class="all",
            sharing_degree="",
            sequenced_cells=sequenced_cells,
            allele_count=allele_count,
            n_variants=n_variants,
        )

    total_unique = int(len(mut_to_orgcounts))
    total_occ = int(sum(agg_bins_occ))
    summary: Dict[str, Any] = {
        "nbins": int(nbins),
        "organ_ids": organ_ids,
        "per_organ": per_organ,
        "aggregate": {
            "vaf_bins": agg_bins_occ,
            "vaf_bins_organ_occurrences": agg_bins_occ,
            "total_variants": total_unique,
            "total_unique_variants": total_unique,
            "total_variant_occurrences": total_occ,
        },
        "vaf_count_spectrum": {
            "available": True,
            "private_shared_available": False,
            "rows": vaf_count_rows,
        },
    }

    if not private_shared:
        summary["private_shared"] = None
        return summary
    if len(organ_events) <= 1:
        summary["private_shared"] = {
            "available": False,
            "reason": "requires at least 2 sampled organs",
        }
        return summary

    for oe in organ_events:
        oid = str(oe["organ_id"])
        mc = oe.get("allele_counts_by_mutation")
        if mc is None:
            summary["private_shared"] = {
                "available": False,
                "reason": "per-organ allele_counts_by_mutation missing",
            }
            return summary

    private_mutations: Dict[str, List[Tuple[int, int, int]]] = {oid: [] for oid in organ_ids}
    shared_mutations: List[Tuple[int, List[Tuple[str, int]]]] = []
    sharing_degree_counts: Dict[int, int] = {d: 0 for d in range(1, len(organ_ids) + 1)}

    for mid, oc in mut_to_orgcounts.items():
        present = list(oc.items())
        degree = len(present)
        sharing_degree_counts[degree] = sharing_degree_counts.get(degree, 0) + 1
        if degree == 1:
            oid, cnt = present[0]
            sequenced_cells = int(per_organ[oid]["sequenced_cells"])
            private_mutations[oid].append((mid, cnt, sequenced_cells))
        else:
            shared_mutations.append((mid, present))

    pairwise: Dict[str, int] = {}
    for i, a in enumerate(organ_ids):
        for b in organ_ids[i + 1:]:
            pairwise[f"{a}|{b}"] = 0
    for _, present in shared_mutations:
        oset = sorted([oid for oid, _ in present])
        for i, a in enumerate(oset):
            for b in oset[i + 1:]:
                key = f"{a}|{b}" if f"{a}|{b}" in pairwise else f"{b}|{a}"
                pairwise[key] = pairwise.get(key, 0) + 1

    private_bins_by_organ: Dict[str, List[int]] = {oid: [0 for _ in range(nbins)] for oid in organ_ids}
    shared_bins_by_organ: Dict[str, List[int]] = {oid: [0 for _ in range(nbins)] for oid in organ_ids}
    shared_bins_global = [0 for _ in range(nbins)]
    private_bins_global = [0 for _ in range(nbins)]
    class_count_accumulator: Dict[Tuple[str, str, str, str, int, int], int] = {}

    def add_class_count(level: str, organ_id: str, variant_class: str, sharing_degree: Any, sequenced_cells: int, allele_count: int, n_variants: int = 1) -> None:
        key = (str(level), str(organ_id), str(variant_class), str(sharing_degree) if sharing_degree is not None else "", int(sequenced_cells), int(allele_count))
        class_count_accumulator[key] = class_count_accumulator.get(key, 0) + int(n_variants)

    for oid, lst in private_mutations.items():
        for (_, cnt, sequenced_cells) in lst:
            vaf = float(cnt) / float(sequenced_cells)
            idx = _vaf_bin_index(vaf, nbins)
            private_bins_by_organ[oid][idx] += 1
            private_bins_global[idx] += 1
            add_class_count("organ", oid, "private", 1, sequenced_cells, cnt)
            add_class_count("global", "", "private", 1, sequenced_cells, cnt)
            add_class_count("organ", oid, "degree", 1, sequenced_cells, cnt)
            add_class_count("global", "", "degree", 1, sequenced_cells, cnt)

    for _, present in shared_mutations:
        degree = len(present)
        for oid, cnt in present:
            sequenced_cells = int(per_organ[oid]["sequenced_cells"])
            vaf = float(cnt) / float(sequenced_cells)
            idx = _vaf_bin_index(vaf, nbins)
            shared_bins_by_organ[oid][idx] += 1
            shared_bins_global[idx] += 1
            add_class_count("organ", oid, "shared", ">=2", sequenced_cells, cnt)
            add_class_count("global", "", "shared", ">=2", sequenced_cells, cnt)
            add_class_count("organ", oid, "degree", degree, sequenced_cells, cnt)
            add_class_count("global", "", "degree", degree, sequenced_cells, cnt)

    for (level, organ_id, variant_class, sharing_degree, sequenced_cells, allele_count), n_variants in sorted(class_count_accumulator.items()):
        _add_count_spectrum_row(
            vaf_count_rows,
            level=level,
            organ_id=organ_id,
            variant_class=variant_class,
            sharing_degree=sharing_degree,
            sequenced_cells=sequenced_cells,
            allele_count=allele_count,
            n_variants=n_variants,
        )

    summary["vaf_count_spectrum"] = {
        "available": True,
        "private_shared_available": True,
        "rows": vaf_count_rows,
    }

    summary["private_shared"] = {
        "available": True,
        "n_mutations_total": int(len(mut_to_orgcounts)),
        "n_shared_mutations": int(len(shared_mutations)),
        "n_private_mutations_total": int(sum(len(v) for v in private_mutations.values())),
        "n_private_by_organ": {oid: int(len(private_mutations[oid])) for oid in organ_ids},
        "pairwise_shared_counts": pairwise,
        "sharing_degree_counts": {str(k): int(v) for k, v in sorted(sharing_degree_counts.items())},
        "private_vaf_bins_by_organ": private_bins_by_organ,
        "shared_vaf_bins_by_organ": shared_bins_by_organ,
        "private_vaf_bins_global": private_bins_global,
        "shared_vaf_bins_global": shared_bins_global,
    }

    return summary


def aggregate_replicate_summaries(
    summaries: List[Dict[str, Any]],
) -> Dict[str, Any]:
    if not summaries:
        return {}

    nbins = int(summaries[0].get("nbins", 20))
    organ_ids = list(summaries[0].get("organ_ids", []))

    def _mean_sd_lists(vs: List[List[int]]) -> Dict[str, Any]:
        if not vs:
            return {"mean": [], "sd": []}
        width = len(vs[0])
        m = [statistics.mean([v[i] for v in vs]) for i in range(width)]
        sd = [statistics.pstdev([v[i] for v in vs]) for i in range(width)]
        return {"mean": m, "sd": sd}

    def _mean_sd_scalars(vs: List[Any]) -> Dict[str, Any]:
        vs2 = [v for v in vs if v is not None]
        if not vs2:
            return {"mean": None, "sd": None}
        return {"mean": statistics.mean(vs2), "sd": statistics.pstdev(vs2)}

    def _mean_sd_dicts_by_keys(dicts: List[Dict[int, int]], keys: List[int]) -> Dict[str, Any]:
        if not dicts:
            return {"keys": list(keys), "mean": [], "sd": []}
        key_ints = [int(k) for k in keys]
        mean = [statistics.mean([int(d.get(k, 0)) for d in dicts]) for k in key_ints]
        sd = [statistics.pstdev([int(d.get(k, 0)) for d in dicts]) for k in key_ints]
        return {"keys": key_ints, "mean": mean, "sd": sd}

    out: Dict[str, Any] = {
        "nbins": nbins,
        "organ_ids": organ_ids,
    }

    per_organ: Dict[str, Any] = {}
    scalar_diag_fields = [
        "n_lineages_sampled", "dominant_lineage_id", "dominant_fraction", "lineage_entropy",
        "boundary_cells", "organ_total_cells", "organ_precursor_number",
    ]
    bool_diag_fields = []
    text_fields = []
    for oid in organ_ids:
        organ_summaries = [s["per_organ"][oid] for s in summaries]
        allele_keys = sorted({int(k) for os in organ_summaries for k in os.get("counts_by_allele_count", {}).keys()})
        counts_by_k_dicts = [
            {int(k): int(v) for k, v in os.get("counts_by_allele_count", {}).items()}
            for os in organ_summaries
        ]
        vaf_bin_lists = [list(map(int, os.get("vaf_bins", []))) for os in organ_summaries]
        row = {
            "sequenced_cells": int(organ_summaries[0]["sequenced_cells"]),
            "counts_by_allele_count": _mean_sd_dicts_by_keys(counts_by_k_dicts, allele_keys),
            "vaf_bins": _mean_sd_lists(vaf_bin_lists),
            "total_variants": _mean_sd_scalars([int(os["total_variants"]) for os in organ_summaries]),
            "singletons": _mean_sd_scalars([int(os["singletons"]) for os in organ_summaries]),
            "fixed": _mean_sd_scalars([int(os["fixed"]) for os in organ_summaries]),
        }
        for fld in scalar_diag_fields:
            row[fld] = _mean_sd_scalars([os.get(fld) for os in organ_summaries])
        for fld in bool_diag_fields:
            row[fld] = _mean_sd_scalars([1 if bool(os.get(fld)) else 0 for os in organ_summaries])
        for fld in text_fields:
            vals = [os.get(fld) for os in organ_summaries if os.get(fld) is not None]
            row[fld] = vals[0] if vals else None
        per_organ[oid] = row
    out["per_organ"] = per_organ

    agg = [s.get("aggregate", {}) for s in summaries]
    out["aggregate"] = {
        "total_variants": _mean_sd_scalars([int(a.get("total_variants", 0)) for a in agg]),
        "total_unique_variants": _mean_sd_scalars([int(a.get("total_unique_variants", a.get("total_variants", 0))) for a in agg]),
        "total_variant_occurrences": _mean_sd_scalars([int(a.get("total_variant_occurrences", 0)) for a in agg]),
        "vaf_bins_organ_occurrences": _mean_sd_lists([list(map(int, a.get("vaf_bins_organ_occurrences", a.get("vaf_bins", [])))) for a in agg]),
    }
    out["vaf_count_spectrum"] = _aggregate_count_spectrum_rows(summaries)

    ps_all = [s.get("private_shared") for s in summaries]
    ps0 = ps_all[0]
    if not isinstance(ps0, dict) or not ps0.get("available"):
        reason = None
        if isinstance(ps0, dict):
            reason = ps0.get("reason")
        out["private_shared"] = {
            "available": False,
            "reason": reason or "requires at least 2 sampled organs",
        }
        return out

    out_ps: Dict[str, Any] = {
        "available": True,
        "n_shared_mutations": _mean_sd_scalars([int(ps["n_shared_mutations"]) for ps in ps_all]),
        "n_private_mutations_total": _mean_sd_scalars([int(ps["n_private_mutations_total"]) for ps in ps_all]),
    }

    n_private_by_organ: Dict[str, Any] = {}
    for oid in organ_ids:
        n_private_by_organ[oid] = _mean_sd_scalars([
            int(ps.get("n_private_by_organ", {}).get(oid, 0)) for ps in ps_all
        ])
    out_ps["n_private_by_organ"] = n_private_by_organ

    pairwise_keys = sorted({k for ps in ps_all for k in ps.get("pairwise_shared_counts", {}).keys()})
    pairwise_shared_counts: Dict[str, Any] = {}
    for key in pairwise_keys:
        pairwise_shared_counts[key] = _mean_sd_scalars([
            int(ps.get("pairwise_shared_counts", {}).get(key, 0)) for ps in ps_all
        ])
    out_ps["pairwise_shared_counts"] = pairwise_shared_counts

    degree_keys = sorted({str(k) for ps in ps_all for k in ps.get("sharing_degree_counts", {}).keys()}, key=lambda x: int(x))
    sharing_degree_counts: Dict[str, Any] = {}
    for key in degree_keys:
        sharing_degree_counts[key] = _mean_sd_scalars([
            int(ps.get("sharing_degree_counts", {}).get(str(key), ps.get("sharing_degree_counts", {}).get(int(key), 0))) for ps in ps_all
        ])
    out_ps["sharing_degree_counts"] = sharing_degree_counts

    def _aggregate_bins_dict(key: str, *, by_organ: bool) -> Any:
        if by_organ:
            out_by_organ: Dict[str, Any] = {}
            for oid in organ_ids:
                out_by_organ[oid] = _mean_sd_lists([
                    list(map(int, ps.get(key, {}).get(oid, [0 for _ in range(nbins)]))) for ps in ps_all
                ])
            return out_by_organ
        return _mean_sd_lists([list(map(int, ps.get(key, [0 for _ in range(nbins)]))) for ps in ps_all])

    out_ps["private_vaf_bins_global"] = _aggregate_bins_dict("private_vaf_bins_global", by_organ=False)
    out_ps["shared_vaf_bins_global"] = _aggregate_bins_dict("shared_vaf_bins_global", by_organ=False)
    out_ps["private_vaf_bins_by_organ"] = _aggregate_bins_dict("private_vaf_bins_by_organ", by_organ=True)
    out_ps["shared_vaf_bins_by_organ"] = _aggregate_bins_dict("shared_vaf_bins_by_organ", by_organ=True)

    out["private_shared"] = out_ps
    return out
