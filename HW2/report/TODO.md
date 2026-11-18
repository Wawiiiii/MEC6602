# HW2 report: to do

State 2026-11-18 midday: `report.tex` compiles with no warnings (22 pages; body pages 1-8).
Body: Section 1 (numerical method, BCs introduced, details in Appendix C), Section 2 (grid
convergence, M 0.5, Hemker-Koren), Section 3 (subsonic / transonic / supersonic, shortened by the
user), Section 4 (fixed lift, Table 5, Mach and Cp at NC 1024). Section 5 (performance) deleted.
Appendices: A Roe, B MUSCL, C boundary conditions, D GMRES, E FMG, F fixed-lift algorithm,
G Richardson, H-J grid convergence (M 0.5 VA, M 0.8 HK, M 0.8 VA), K Van Albada force tables.

## To do

1. [done] Section 4: fixed lift coefficient (user edits 2026-11-18).
2. [done] Section 5 performance deleted (user). The statement asks to "report on the convergence
   and execution time": run times are in the force tables (t columns); convergence plots pending
   (item 7).
3. [done] Table 5 without Iter. and N_alpha; Cp* removed from every Cp plot; Table 4 without
   the Cauchy daggers.
4. **Code cleanup** (solver and `post/` scripts).
5. **Reread the whole report.**
6. [done] Figure and table titles simplified (2026-11-18, user); "Figs." -> "Figures"; no figure
   reference after an appendix reference.
7. [done] **Convergence plots** (residual drop against iteration, one figure per case of the report body,
   grids 129 to 2049): `post/convergence_history.py` -> `report/figures/conv_<case>.png`, 600 DPI.
   Ready: M05_a125, M08_a0, M08_a125, M125_a0, M125_a125, M08_CL04 (129 to 1025).
   No M05_a0 figure (no terminal output kept; user: the alpha 1.25 one is enough).
   Where they go in the report: user to decide.
8. **Correct improper information** (check every number and statement against the results).

## Notes

- Figures: `post/grid_convergence_plots.py` (Section 2, grid-convergence appendices),
  `post/flow_results.py` (Section 3 Mach and Cp, force tables), `post/fixed_cl_results.py`
  (Section 4), `post/convergence_history.py` (convergence plots). Run from the HW2 folder.
- Run times in the tables are the fine-grid solve only (no FMG start), on different machines and
  thread counts.
- The 2049 fixed-CL run was stopped 2026-11-18 (not needed, user); its files are kept.
- Still pending: commit and push of the code, notes and report (only when the user says so).
