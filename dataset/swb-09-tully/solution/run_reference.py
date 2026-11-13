"""Reproduce one selected Fig. 2 temperature from the self-contained package."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from tully_ba_ref import (  # noqa: E402
    AU_VELOCITY_TO_M_PER_S,
    FS_TO_AU,
    Tully2012ReactionModel,
    estimate_mpes_qtst,
    run_transmission_batch,
    sample_constrained_configurations,
    thermal_adiabatic_rate,
    thermal_nonadiabatic_rate,
)
from tully_ba_ref.rate_analysis import (  # noqa: E402
    assemble_eq40_rate,
    combine_independent_batches,
)


def _barrier(model: Tully2012ReactionModel) -> tuple[float, float]:
    result = minimize_scalar(
        lambda x: -float(model.get_energies(float(x))[0, 0]),
        bounds=(-3.0, 3.0),
        method="bounded",
        options={"xatol": 1.0e-13},
    )
    if not result.success:
        raise RuntimeError(result.message)
    return float(result.x), float(-result.fun)


def _exact_level(
    temperature: float,
    xmin: float,
    xmax: float,
    ngrid: int,
    nquadrature: int,
) -> dict:
    model = Tully2012ReactionModel()
    started = time.perf_counter()
    adiabatic = thermal_adiabatic_rate(
        model,
        temperature=temperature,
        xmin=xmin,
        xmax=xmax,
        ngrid=ngrid,
        nquadrature=nquadrature,
    )
    nonadiabatic = thermal_nonadiabatic_rate(
        model,
        temperature=temperature,
        xmin=xmin,
        xmax=xmax,
        ngrid=ngrid,
        nquadrature=nquadrature,
    )
    return {
        "xmin_bohr": xmin,
        "xmax_bohr": xmax,
        "ngrid": ngrid,
        "nquadrature": nquadrature,
        "adiabatic_rate_au": adiabatic.total_rate,
        "nonadiabatic_rate_au": nonadiabatic.total_rate,
        "nonadiabatic_rate_by_incoming_state_au": nonadiabatic.rate_by_incoming_state.tolist(),
        "wall_seconds": time.perf_counter() - started,
    }


def _batch_worker(job: dict) -> dict:
    model = Tully2012ReactionModel()
    started = time.perf_counter()
    configurations = sample_constrained_configurations(
        model,
        centroid=job["dividing_surface"],
        state=job["surface_state"],
        temperature=job["temperature"],
        nbeads=job["nbeads"],
        nconfigurations=job["ntraj"],
        burnin_steps=job["burnin_steps"],
        steps_between=job["steps_between"],
        seed=job["seed"],
    )
    batch = run_transmission_batch(
        model,
        configurations.positions,
        surface_state=job["surface_state"],
        incoming_state=job["incoming_state"],
        temperature=job["temperature"],
        dt=job["dt_fs"] * FS_TO_AU,
        approximation="BA",
        dividing_surface=job["dividing_surface"],
        reactant_boundary=job["reactant_boundary"],
        product_boundary=job["product_boundary"],
        max_steps=job["max_steps"],
        seed=job["seed"] + 1_000_000,
        nblocks=min(5, max(2, job["ntraj"] // 2)),
        event_convention=job["event_convention"],
        frustrated_hop_policy=job["frustrated_hop_policy"],
        electronic_integrator=job["electronic_integrator"],
    )
    return {
        "surface_state": job["surface_state"],
        "incoming_state": job["incoming_state"],
        "chunk": job["chunk"],
        "seed": job["seed"],
        "ntraj": job["ntraj"],
        "value": batch.total.value,
        "standard_error": batch.total.standard_error,
        "configuration_acceptance": configurations.acceptance_rate,
        "configuration_scale": configurations.final_proposal_scale,
        "reactive_count": int(np.sum(batch.reactive)),
        "invalid_history_count": int(np.sum(batch.final_states < 0)),
        "weights": batch.backward_weights.tolist(),
        "crossing_counts": batch.crossing_counts.tolist(),
        "backward_attempts": int(np.sum(batch.backward_attempts)),
        "backward_accepted": int(np.sum(batch.backward_accepted)),
        "forward_attempts": int(np.sum(batch.forward_attempts)),
        "forward_accepted": int(np.sum(batch.forward_accepted)),
        "wall_seconds": time.perf_counter() - started,
    }


def _rate_payload(rate_au: float, standard_error_au: float = 0.0) -> dict:
    converted = rate_au * AU_VELOCITY_TO_M_PER_S
    relative = standard_error_au / rate_au if rate_au > 0.0 else math.inf
    return {
        "rate_au": rate_au,
        "standard_error_au": standard_error_au,
        "rate_m_per_molecule_s": converted,
        "log10_rate_m_per_molecule_s": math.log10(converted) if converted > 0.0 else None,
        "log10_standard_error_dex": relative / math.log(10.0) if math.isfinite(relative) else None,
    }


def _summarize_transmission(results: list[dict], qtst_rate: np.ndarray, qtst_se: np.ndarray) -> dict:
    channels: dict[str, dict] = {}
    values: dict[tuple[int, int], float] = {}
    errors: dict[tuple[int, int], float] = {}
    for surface in range(2):
        for incoming in range(2):
            selected = [
                item
                for item in results
                if item["surface_state"] == surface and item["incoming_state"] == incoming
            ]
            estimate = combine_independent_batches(
                np.array([item["value"] for item in selected]),
                np.array([item["standard_error"] for item in selected]),
                np.array([item["ntraj"] for item in selected]),
            )
            key = (surface, incoming)
            values[key] = estimate.value
            errors[key] = estimate.standard_error
            channels[f"surface_{surface}_incoming_{incoming}"] = {
                "value": estimate.value,
                "standard_error": estimate.standard_error,
                "within_batch_standard_error": estimate.within_batch_standard_error,
                "between_batch_standard_error": estimate.between_batch_standard_error,
                "trajectories": estimate.total_trajectories,
                "chunks": estimate.batch_count,
            }
    total, total_se, missing = assemble_eq40_rate(qtst_rate, qtst_se, values, errors)
    all_weights = np.concatenate([np.asarray(item["weights"], dtype=float) for item in results])
    positive = all_weights[all_weights > 0.0]
    return {
        "channels": channels,
        "total": _rate_payload(total, total_se),
        "missing_channels": [list(item) for item in missing],
        "diagnostics": {
            "trajectory_count": int(sum(item["ntraj"] for item in results)),
            "reactive_count": int(sum(item["reactive_count"] for item in results)),
            "invalid_history_count": int(sum(item["invalid_history_count"] for item in results)),
            "configuration_acceptance_range": [
                float(min(item["configuration_acceptance"] for item in results)),
                float(max(item["configuration_acceptance"] for item in results)),
            ],
            "backward_weight_quantiles": np.quantile(all_weights, [0.0, 0.5, 0.9, 0.99, 1.0]).tolist(),
            "log_backward_weight_quantiles": (
                np.quantile(np.log(positive), [0.0, 0.5, 0.9, 0.99, 1.0]).tolist()
                if positive.size
                else []
            ),
            "backward_attempts": int(sum(item["backward_attempts"] for item in results)),
            "backward_accepted": int(sum(item["backward_accepted"] for item in results)),
            "forward_attempts": int(sum(item["forward_attempts"] for item in results)),
            "forward_accepted": int(sum(item["forward_accepted"] for item in results)),
            "chunk_wall_seconds": [item["wall_seconds"] for item in results],
        },
        "chunks": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("smoke", "reference"), default="reference")
    parser.add_argument("--config", type=Path, default=ROOT / "config.json")
    parser.add_argument("--output", type=Path, default=ROOT / "output" / "reference_250K.json")
    parser.add_argument("--nbeads", type=int, choices=(16, 32), help="override the bead count in config.json")
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    if args.nbeads is not None:
        config["nbeads"] = args.nbeads
    if args.profile == "smoke":
        config["qtst"].update({"windows": 7, "steps_per_window": 300, "seeds": [9101, 9102]})
        config["transmission"].update(
            {"trajectories_per_channel": 8, "chunks_per_channel": 2, "burnin_steps": 300, "steps_between": 10, "workers": 4}
        )
        # The macOS ARM64 smoke run follows the same deterministic seeds as
        # Windows.  A small subset can remain near the barrier longer because
        # of platform-level floating-point divergence, so retain a generous
        # completion margin while keeping the sampling profile tiny.
        config["max_steps"] = 400000

    overall_started = time.perf_counter()
    model = Tully2012ReactionModel()
    dividing_surface, barrier_energy = _barrier(model)

    exact_started = time.perf_counter()
    exact_levels = [
        _exact_level(config["temperature_K"], -12.0, 12.0, 1600, 64),
        _exact_level(config["temperature_K"], -16.0, 16.0, 2400, 96),
    ]
    exact_wall = time.perf_counter() - exact_started
    exact_final = exact_levels[-1]

    qtst_started = time.perf_counter()
    qtst, chains = estimate_mpes_qtst(
        model,
        temperature=config["temperature_K"],
        nbeads=config["nbeads"],
        reactant_surface=config["qtst"]["reactant_centroid_bohr"],
        dividing_surface=dividing_surface,
        nwindows=config["qtst"]["windows"],
        nsteps=config["qtst"]["steps_per_window"],
        seeds=tuple(config["qtst"]["seeds"]),
        burnin_fraction=config["qtst"]["burnin_fraction"],
    )
    qtst_wall = time.perf_counter() - qtst_started

    transmission = config["transmission"]
    total_per_channel = int(transmission["trajectories_per_channel"])
    chunks = int(transmission["chunks_per_channel"])
    if total_per_channel % chunks:
        raise ValueError("trajectories_per_channel must be divisible by chunks_per_channel")
    per_chunk = total_per_channel // chunks
    jobs = []
    for surface in range(2):
        for incoming in range(2):
            for chunk in range(chunks):
                seed = int(transmission["seed"] + 100000 * surface + 10000 * incoming + chunk)
                jobs.append(
                    {
                        "temperature": config["temperature_K"],
                        "nbeads": config["nbeads"],
                        "dt_fs": config["dt_fs"],
                        "reactant_boundary": config["reactant_boundary_bohr"],
                        "product_boundary": config["product_boundary_bohr"],
                        "max_steps": config["max_steps"],
                        "event_convention": config["event_convention"],
                        "frustrated_hop_policy": config["frustrated_hop_policy"],
                        "electronic_integrator": config["electronic_integrator"],
                        "dividing_surface": dividing_surface,
                        "surface_state": surface,
                        "incoming_state": incoming,
                        "chunk": chunk,
                        "seed": seed,
                        "ntraj": per_chunk,
                        "burnin_steps": transmission["burnin_steps"],
                        "steps_between": transmission["steps_between"],
                    }
                )

    transmission_started = time.perf_counter()
    chunk_results = []
    with ProcessPoolExecutor(max_workers=int(transmission["workers"])) as pool:
        futures = [pool.submit(_batch_worker, job) for job in jobs]
        for future in as_completed(futures):
            result = future.result()
            chunk_results.append(result)
            print(
                f"completed s={result['surface_state']} in={result['incoming_state']} "
                f"chunk={result['chunk']} kappa={result['value']:.6g} "
                f"wall={result['wall_seconds']:.1f}s",
                flush=True,
            )
    transmission_wall = time.perf_counter() - transmission_started
    chunk_results.sort(key=lambda item: (item["surface_state"], item["incoming_state"], item["chunk"]))
    ba = _summarize_transmission(chunk_results, qtst.rate_mean, qtst.rate_standard_error)

    fine_a = exact_final["adiabatic_rate_au"]
    fine_na = exact_final["nonadiabatic_rate_au"]
    coarse_a = exact_levels[-2]["adiabatic_rate_au"]
    coarse_na = exact_levels[-2]["nonadiabatic_rate_au"]
    payload = {
        "status": "complete",
        "profile": args.profile,
        "method": "Shushkov-Li-Tully 2012 Eqs. 36, 40, 50-56; standalone implementation",
        "temperature_K": config["temperature_K"],
        "nbeads": config["nbeads"],
        "dividing_surface_bohr": dividing_surface,
        "lower_barrier_hartree": barrier_energy,
        "exact": {
            "levels": exact_levels,
            "adiabatic": _rate_payload(fine_a),
            "nonadiabatic": _rate_payload(fine_na),
            "adiabatic_refinement_change_dex": abs(math.log10(fine_a / coarse_a)),
            "nonadiabatic_refinement_change_dex": abs(math.log10(fine_na / coarse_na)),
        },
        "qtst": {
            "rate_au_by_surface": qtst.rate_mean.tolist(),
            "standard_error_au_by_surface": qtst.rate_standard_error.tolist(),
            "minimum_chain_ess": qtst.minimum_chain_ess,
            "acceptance_range": list(qtst.acceptance_range),
            "windows": int(qtst.grid.size),
            "steps_per_chain": config["qtst"]["steps_per_window"],
            "chain_count": len(chains),
        },
        "rpsh_ba": ba,
        "runtime": {
            "exact_seconds": exact_wall,
            "qtst_seconds": qtst_wall,
            "transmission_seconds": transmission_wall,
            "total_seconds": time.perf_counter() - overall_started,
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "logical_cpu_count": os.cpu_count(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "provenance": {
            "config_sha256": hashlib.sha256(args.config.read_bytes()).hexdigest(),
            "command": " ".join(sys.argv),
            "source_runtime_dependency": "none; only files under this standalone directory plus NumPy/SciPy",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(args.output), "runtime": payload["runtime"], "rates": {"exact_A": payload["exact"]["adiabatic"], "exact_NA": payload["exact"]["nonadiabatic"], "BA": payload["rpsh_ba"]["total"]}}, indent=2), flush=True)


if __name__ == "__main__":
    main()
