# SciWorkBench

SciWorkBench is a scientific agent benchmark. Agents analyze scientific inputs, produce reproducible artifacts and justify conclusions. This release bundles the inputs and Harbor verifiers for ten selected tasks across biology, chemistry and physics.

| ID | Domain | Task | Input coverage |
| --- | --- | --- | --- |
| SWB-01 | Biology | [Soma-free 3D neurite tracing](dataset/swb-01-neurite/instruction.md) | Bundled source inputs |
| SWB-02 | Biology | [AGO2 perturbation and RNA kinetics](dataset/swb-02-ago2/instruction.md) | Bundled source inputs |
| SWB-03 | Biology | [In vivo calcium imaging](dataset/swb-03-calcium/instruction.md) | Bundled source inputs |
| SWB-04 | Biology | [Dual-gate FLIM cardiac electrophysiology](dataset/swb-04-flim/instruction.md) | Sparse source samples |
| SWB-05 | Chemistry | [Organic structure elucidation with NMR NOESY and HRMS](dataset/swb-05-structure/instruction.md) | Bundled source inputs |
| SWB-06 | Chemistry | [Molecular generation claims and docking verification](dataset/swb-06-docking/instruction.md) | Bundled source inputs |
| SWB-07 | Chemistry | [Raman spectroscopy of Fe-doped monolayer MoS₂](dataset/swb-07-raman/instruction.md) | Bundled source inputs |
| SWB-08 | Physics | [Flare-aware transit validation with TESS](dataset/swb-08-tess/instruction.md) | Bundled source inputs |
| SWB-09 | Physics | [Quantum reaction rates with ring-polymer surface hopping](dataset/swb-09-tully/instruction.md) | Bundled source inputs |
| SWB-10 | Physics | [Interface-localized electronic modulation in GdTe₂](dataset/swb-10-gdte2/instruction.md) | Bundled source inputs |

The three refreshed physics tasks occupy SWB-08–10. SWB-04 provides the original sparse FLIM export and evaluates what can be inferred from its observed samples. The supplied export cannot support a full-frame or complete beat-level analysis. See [TASKS.md](TASKS.md) for provenance and scoring scope.

## Run with Harbor

Install uv, start Docker and run from the repository root.

```sh
uvx --from harbor==0.25.0 harbor lint dataset/*
uvx --from harbor==0.25.0 harbor run -p dataset -a "$AGENT" -m "$MODEL"
```

Set `AGENT`, `MODEL` and the model provider credentials before running. See [RUN.md](RUN.md) for single-task and reference runs. The [task index](tasks.csv) gives each task's stable ID and path.

## Reuse and citation

The repository includes grading references and reference solutions for reproducibility. When evaluating agents, make only the task instruction and Harbor agent environment available to the agent. Keep the verifier and reference solution on the evaluation host.

See [LICENSE](LICENSE) for the project code and documentation and [DATA.md](DATA.md) for task materials and upstream terms.

When citing the benchmark, use *SciWorkBench: A Benchmark for Real Scientific Workflows*. For a run on this release, also report the repository commit and task IDs.
