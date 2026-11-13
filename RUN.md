# Run SciWorkBench with Harbor

Use `uv` and a running Docker engine. From the repository root, validate the task configuration and run an agent:

```sh
uvx --from harbor==0.25.0 harbor lint dataset/*
# Set AGENT to a Harbor agent name and MODEL to the provider/model identifier.
# Configure credentials for that provider in your shell.
uvx --from harbor==0.25.0 harbor run -p dataset -a "$AGENT" -m "$MODEL"
```

Use `-p dataset/swb-08-tess` to run just TESS. The [Harbor CLI reference](https://github.com/harbor-framework/docs/blob/main/cli/run.mdx) describes the agent, model and output options. The dataset path selects all ten entries in [tasks.csv](tasks.csv). SWB-04 is labeled `runnable_sparse` because its original photon-count export contains sampled observations.

To check reference workflows without calling a model:

```sh
uvx --from harbor==0.25.0 harbor run -p dataset/swb-09-tully -a oracle
uvx --from harbor==0.25.0 harbor run -p dataset/swb-08-tess -a oracle
uvx --from harbor==0.25.0 harbor run -p dataset/swb-10-gdte2 -a oracle
```

TESS and GdTe₂ recompute outputs from their inputs. Tully rebuilds its deliverables from frozen calculations and includes the underlying calculation code. The other tasks include grading references under `tests/`.

Agent images include the input data and dependencies. Harbor runs each task's verifier separately and passes it the submitted artifacts. Graders write a reward from 0 to 1. Some verifiers automate a subset or approximation of the original scientific rubric. See [TASKS.md](TASKS.md) for task-specific coverage. Report the task IDs, commit, agent, model and per-task results. Compare runs only across matching task and verifier versions.
