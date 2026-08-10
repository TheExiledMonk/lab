"""Discrete channel/deposition compatibility analysis for Dev Doc 113."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import combinations
import math
import statistics
from typing import Iterable


CLUSTERS = ("Abell2744", "MACS0416", "MACS1149", "AbellS1063", "Abell370")
SURVIVORS = (
    "hard_bin_half_open", "nearest_center", "bilinear_cic", "tsc_3x3",
    "gaussian_sigma_half_cell",
)
COMPLEXITY = {name: rank for rank, name in enumerate(SURVIVORS, 1)}
FINAL_STATUS = "OBSERVER_STABLE_DEPOSITION_CANDIDATE_ESTABLISHED"


@dataclass(frozen=True)
class ObservableCompatibility:
    observable: str
    deposition: str
    channel_source: str
    reconstruction: str
    stability_score: float
    cross_cluster_score: float
    information_score: float
    morphology_score: float
    complexity_rank: int


def validate_final_audit(audit: dict) -> dict[str, bool]:
    clusters = audit.get("clusters", {})
    methods = [clusters.get(c, {}).get("methods", {}) for c in CLUSTERS]
    checks = {
        "canonical_five_cluster_inventory": tuple(clusters) == CLUSTERS,
        "five_stable_deposition_survivors_loaded": tuple(audit.get("stability_survivors", ())) == SURVIVORS,
        "45_channel_inventory_present": all(all(m.get(s, {}).get("channel_count") == 45 for s in SURVIVORS) for m in methods),
        "68_candidate_inventory_present": all(all(m.get(s, {}).get("candidate_count") == 68 for s in SURVIVORS) for m in methods),
        "backend_matrix_present": all(c in audit.get("backend_matrix", {}) for c in CLUSTERS),
        "cache_stats_present": all(c in audit.get("cache_stats", {}) for c in CLUSTERS),
        "scientific_checks_present": isinstance(audit.get("checks"), dict) and bool(audit["checks"]),
        "final_status_valid": audit.get("status") == FINAL_STATUS,
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError("incomplete final.log evidence: " + ", ".join(failed))
    return checks


def component_availability(audit: dict) -> dict[str, bool]:
    available = {name: True for name in ("kappa", "gamma1", "gamma2")}
    for cluster in CLUSTERS:
        for method in SURVIVORS:
            diagnostics = audit["clusters"][cluster]["methods"][method].get("observational_diagnostics", {})
            for observable in available:
                available[observable] &= bool(diagnostics) and all(
                    observable in row for row in diagnostics.values()
                )
    return available


def channel_source(reconstruction: str) -> str:
    if reconstruction.startswith("family_"):
        return reconstruction[len("family_"):].rsplit("_", 1)[0]
    if reconstruction.startswith("explicit3d_"):
        return "explicit_3d"
    if reconstruction.startswith("established2d_"):
        return "established_2d"
    if reconstruction.startswith("nodepth_"):
        return "all_except_depth_3d"
    return "all_45_channels"


def _robustness(rows: list[dict]) -> dict:
    pearson = [float(r["pearson"]) for r in rows]
    spearman = [float(r["spearman"]) for r in rows]
    finite = all(math.isfinite(x) for x in pearson + spearman)
    return {
        "median_pearson": statistics.median(pearson),
        "minimum_pearson": min(pearson),
        "median_spearman": statistics.median(spearman),
        "minimum_spearman": min(spearman),
        "pearson_stddev": statistics.pstdev(pearson),
        "clusters_positive": sum(x > 0 for x in pearson),
        "clusters_above_0_2": sum(x > .2 for x in pearson),
        "clusters_above_0_4": sum(x > .4 for x in pearson),
        "coverage_fraction_median": statistics.median(float(r["coverage_fraction"]) for r in rows),
        "finite": finite,
        "signed_pearson": pearson,
        "absolute_pearson_diagnostic": [abs(x) for x in pearson],
    }


def _rank_key(row: dict) -> tuple:
    r = row["robustness"]
    # Lexicographic, auditable ordering: consistency before peak performance.
    return (r["clusters_positive"], r["clusters_above_0_2"], r["clusters_above_0_4"],
            r["minimum_pearson"], r["median_pearson"], r["minimum_spearman"],
            row["morphology_score"], row["information_score"],
            -r["pearson_stddev"], -row["complexity_rank"])


def build_matrix(audit: dict, component_result: dict | None = None) -> list[dict]:
    matrix = []
    for observable in ("kappa", "gamma1", "gamma2"):
        source = audit if observable == "kappa" else component_result
        if source is None:
            continue
        candidate_names = None
        for cluster in CLUSTERS:
            for deposition in SURVIVORS:
                diagnostics = (audit["clusters"][cluster]["methods"][deposition]["observational_diagnostics"]
                    if observable == "kappa" else
                    source["clusters"][cluster]["methods"][deposition]["observational_diagnostics"])
                names = set(diagnostics)
                candidate_names = names if candidate_names is None else candidate_names & names
        for deposition in SURVIVORS:
            method_rows = audit["clusters"]
            info = statistics.median(method_rows[c]["methods"][deposition]["information"]["effective_rank"] for c in CLUSTERS)
            morph = statistics.median(method_rows[c]["methods"][deposition]["morphology"]["median_channel_morphology_pearson"] for c in CLUSTERS)
            for reconstruction in sorted(candidate_names or ()):
                rows = []
                for cluster in CLUSTERS:
                    diagnostics = (audit["clusters"][cluster]["methods"][deposition]["observational_diagnostics"]
                        if observable == "kappa" else
                        source["clusters"][cluster]["methods"][deposition]["observational_diagnostics"])
                    rows.append(diagnostics[reconstruction][observable])
                robustness = _robustness(rows)
                record = ObservableCompatibility(observable, deposition,
                    channel_source(reconstruction), reconstruction, 1.0,
                    robustness["median_pearson"], info, morph, COMPLEXITY[deposition])
                matrix.append({**asdict(record), "robustness": robustness,
                    "provenance": "ROOT_FINAL_LOG" if observable == "kappa" else "TARGETED_COMPONENT_COMPLETION"})
    return matrix


def equivalence_classes(matrix: list[dict]) -> dict:
    """Compare complete scalar metric vectors; no array equivalence is claimed."""
    result = {}
    for observable in ("kappa", "gamma1", "gamma2"):
        subset = [r for r in matrix if r["observable"] == observable]
        by_dep = {d: sorted((r["reconstruction"], r["robustness"]) for r in subset if r["deposition"] == d) for d in SURVIVORS}
        pairs = []
        graph = {d: {d} for d in SURVIVORS}
        for a, b in combinations(SURVIVORS, 2):
            same = by_dep[a] == by_dep[b]
            classification = "EXACT_EQUIVALENT" if same else "MEANINGFULLY_DISTINCT"
            pairs.append({"methods": [a, b], "classification": classification,
                          "basis": "logged_scalar_metrics"})
            if same:
                graph[a].add(b); graph[b].add(a)
        classes, seen = [], set()
        for method in SURVIVORS:
            if method not in seen:
                group = tuple(m for m in SURVIVORS if m in graph[method])
                seen.update(group); classes.append(group)
        result[observable] = {"classes": classes, "pairwise": pairs,
                              "array_equivalence_tested": False}
    return result


def select(matrix: list[dict], equivalence: dict) -> dict:
    selections = {}
    for observable in ("kappa", "gamma1", "gamma2"):
        all_rows = [r for r in matrix if r["observable"] == observable and r["robustness"]["finite"]]
        # PCA candidates remain in the descriptive matrix, but selecting one
        # after inspecting these five targets would make the benchmark act as
        # post-hoc training.  Structural L1/L2/signed candidates stay eligible.
        rows = [r for r in all_rows if "_pc" not in r["reconstruction"] and "pca_" not in r["reconstruction"]]
        if not rows:
            selections[observable] = {"status": "UNRESOLVED", "candidates": []}
            continue
        winner = max(rows, key=_rank_key)
        peers = [r for r in rows if r["reconstruction"] == winner["reconstruction"] and
                 r["deposition"] in next(c for c in equivalence[observable]["classes"] if winner["deposition"] in c)]
        representative = min(peers, key=lambda r: r["complexity_rank"])
        omitted_winners = []
        for omitted in range(len(CLUSTERS)):
            def loo_key(row):
                p = [x for i, x in enumerate(row["robustness"]["signed_pearson"]) if i != omitted]
                return (sum(x > 0 for x in p), sum(x > .2 for x in p), sum(x > .4 for x in p),
                        min(p), statistics.median(p), -statistics.pstdev(p), -row["complexity_rank"])
            omitted_winners.append((max(rows, key=loo_key)["deposition"], max(rows, key=loo_key)["reconstruction"]))
        same = sum(x == (winner["deposition"], winner["reconstruction"]) for x in omitted_winners)
        too_weak = (winner["robustness"]["median_pearson"] < .2 and
                    winner["robustness"]["clusters_above_0_2"] == 0)
        confidence = ("UNRESOLVED" if too_weak else
            "ROBUST_EQUIVALENCE_CLASS" if len(peers) > 1 and same == 5 else
            "ROBUST_UNIQUE" if same == 5 else "WEAK_PREFERENCE")
        common = {
            "status": confidence,
            "equivalent_candidates": [{"deposition": r["deposition"], "reconstruction": r["reconstruction"]} for r in peers],
            "reason": "cross-cluster robustness, then morphology, information retention, observational compatibility, and complexity",
            "leave_one_cluster_out": {"winners": omitted_winners, "same_winner_count": same,
                "classification": "same winner 5/5" if same == 5 else "same winner 4/5" if same == 4 else "unstable"},
        }
        if too_weak:
            common["reason"] = "no structurally eligible candidate exceeds Pearson 0.2 in any cluster; a production choice is not justified"
            common["descriptive_leader"] = representative
            common["candidates"] = [
                {"deposition": r["deposition"], "channel_family": r["channel_source"],
                 "reconstruction": r["reconstruction"]}
                for r in sorted(rows, key=_rank_key, reverse=True)[:3]
            ]
        else:
            common["selected"] = representative
        selections[observable] = common
    return selections


def iter_selected_channels(selection: dict) -> Iterable[str]:
    for observable in ("kappa", "gamma1", "gamma2"):
        selected = selection.get(observable, {}).get("selected")
        if selected:
            yield selected["channel_source"]
