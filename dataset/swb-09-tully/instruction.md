# Quantum reaction rates with bead-approximation ring-polymer surface hopping

Reproduce the 143 K and 250 K rate points from Fig. 2 of Shushkov, Li and Tully (2012). Work from the supplied one-dimensional, two-state diabatic Hamiltonian. Perform open-boundary quantum scattering, thermal averaging, ring-polymer sampling and surface-hopping dynamics. Calculate the rates from the model rather than digitizing the paper.

Use `data/problem.json` for constants, model parameters, temperatures and units, `data/scientific_spec.md` for the scientific methods and equations, and `data/output_schema.json` for exact deliverable fields.

- Compute `exact_adiabatic` on the lower adiabatic PES and `exact_nonadiabatic` with the full two-channel Hamiltonian. Use a nuclear mass of 2000 atomic units.
- Compute `rpsh_ba` and its statistical standard error at both temperatures with 32 beads and the bead approximation.
- Place the dividing surface at the lower adiabatic barrier top, with reactants to its left. Compute in atomic units and report rates in both atomic units and m molecule⁻¹ s⁻¹. Round final base-10 log rates and standard errors to three decimals.
- Plot all six rates as log10[k/(m molecule⁻¹ s⁻¹)] versus 1000/T (K⁻¹). Explain whether coupling suppresses the exact rate relative to the adiabatic approximation and whether RPSH-BA reproduces the trend at both temperatures.

Write `results.csv`, `figure2_subset.png` (at least 150 dpi), `method.md`, `conclusions.md`, `diagnostics.json`, `run_all.py` and the source-code directory `src/` to `output/`. Document boundary conditions, discretization, thermal integration, sampling and uncertainties. Cite your computed results in the conclusions. The entrypoint and source must reproduce the reported outputs from the supplied inputs.
