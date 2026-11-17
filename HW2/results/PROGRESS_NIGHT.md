# Night of 2026-11-14/15: grid-study progress (one entry per hour)

Machines (shared disk, results in this folder): j5015-01 = queue j01 (here), j5015-18 = Mach_08_alpha_125 (Hemker-Koren),
j5015-17 = Mach_08_alpha_125_vanalbada, j5015-04 = Mach_125_alpha_0 (Hemker-Koren, 9 -> 1025). Lab queues at nice 19, one case each.

## 2026-11-15 00:48 start
- j5015-01: Mach_08_alpha_0 2049 (HK), iteration 13, drop 3.7e-6 (GMRES at its 400-iteration cap, ~19 min per iteration); then M 0.5 alpha 0 VA, M 0.5 alpha 1.25 HK / VA, M 0.8 alpha 0 VA, M 1.25 alpha 0 VA, M 1.25 alpha 0 HK 2049.
- j5015-18 / -17 / -04: started 00:46, first grids converged (M 0.8 alpha 1.25 HK 9x9 = previous campaign).

## 2026-11-15 01:02 settings changed (user: no more cold starts, stronger GMRES)
- Every run now: FMG start (4 levels) from 257 up; GMRES restart 200 / 800 iterations from 1025 up (was 150 / 400).
- Retries: 2 = GMRES 200/800 on every grid + max_iter 400; 3 = same + cfl_max 1e3, max_iter 1000 (low CFL caps removed).
- All four queues restarted 00:59 with these settings; the solvers already running were taken over (one solver per machine, checked).
- 513 attempt 1 (cold start) had stalled on all three lab machines; FMG retries converged: M 0.8 a 1.25 HK CL 0.354756, VA CL 0.355072.

## 2026-11-15 01:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_08_alpha_0/2049x2049_input.txt | it 15, drop 1.288e-07, CFL 7.76e+06, lin its 400
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/1025x1025_input.txt | it 14, drop 8.616e-04, CFL 1.16e+03, lin its 89
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/1025x1025_input.txt | it 13, drop 3.069e-03, CFL 3.26e+02, lin its 23
- j5015-04 (queue j04: running): results/Mach_125_alpha_0/1025x1025_input.txt 
- alerts in the last hour:
    2026-11-15 00:52:56  NOT CONVERGED Mach_125_alpha_0 513x513 attempt 1: stall at iteration 90 (best drop 3.03e-02, no 10x gain in 50 iterations)
    2026-11-15 00:53:25  NOT CONVERGED Mach_08_alpha_125_vanalbada 513x513 attempt 1: stall at iteration 98 (best drop 2.10e-03, no 10x gain in 50 iterations)
    2026-11-15 00:53:26  NOT CONVERGED Mach_08_alpha_125 513x513 attempt 1: stall at iteration 93 (best drop 2.06e-03, no 10x gain in 50 iterations)

## 2026-11-15 01:10 hourly check (Claude)
- All 4 queues alive, 1 solver each, no new alerts since 00:55.
- j5015-01, M 0.8 a0 HK 2049: it 15, drop 1.3e-7 (3-8x per iteration, GMRES at its 400 cap, CFL 7.8e6), CD 8.34736e-3 steady. Converging; ~3-4 iterations (~1 h) left.
- j5015-18, M 0.8 a1.25 HK 1025 (FMG + GMRES 200/800): it 15, drop 4-9e-4 for 4 iterations, CFL bouncing 1.2e3-2.3e3, CL still moving (0.3547 -> 0.3540). Shock settling (B20-like); watching. Finished: 513 CL 0.35475593, CD 2.26277e-2, CM -3.89402e-2 (FMG retry, 23 it).
- j5015-17, M 0.8 a1.25 VA 1025: it 13, drop 2e-4 then back to 3e-3 when CFL passed 4.4e3 (CFL cut to 326), CL moving. Same pattern; watching. Finished: 513 CL 0.35507207, CD 2.26457e-2, CM -3.90358e-2 (23 it).
- j5015-04, M 1.25 a0 HK 1025: in the FMG start (129, 257 done). Finished: 513 CD 9.61440e-2 (FMG retry, 30 it).

## 2026-11-15 02:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_08_alpha_0/2049x2049_input.txt | it 18, drop 4.230e-08, CFL 2.364e+05, lin its linear
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt 
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt 
- j5015-04 (queue j04: running): results/Mach_125_alpha_0/1025x1025_input.txt | it 40, drop 3.367e-02, CFL 2.97e+01, lin its 57

## 2026-11-15 02:10 hourly check (Claude)
- All 4 queues alive, no new alerts.
- j5015-01, M 0.8 a0 HK 2049: still at it 18 (3 steps rejected: GMRES 150/400 failed at CFL 2.4e7, 2.4e6, 2.4e5; ~19 min each). Next try ~CFL 2.4e4. Forces converged (CD 8.347367e-3 since it 15, CL -6e-8); only the 1e-12 drop is missing. Left running (a restart would cost 5+ h).
- j5015-18, M 0.8 a1.25 HK: 1025 CONVERGED (25 it, drop 4.1e-14): CL 0.35347803, CD 2.254034e-2, CM -3.855271e-2. Now 2049, FMG start (513 level done).
- j5015-17, M 0.8 a1.25 VA: 2049, FMG start (513 level done).
- j5015-04, M 1.25 a0 HK 1025: it 41, drop bouncing 5e-3 to 3e-2, CFL collapses each time it reaches ~150-200 (B20 pattern; best 4.4e-3 at it 32). 513 did the same and broke through at it 22. Stall rule would stop it at it ~82; attempt 3 then = same + max_iter 1000, window 150. Watching.

## 2026-11-15 02:50 stall handled (Claude)
- j5015-04, M 1.25 a0 HK 1025: attempt 1 stalled at it 84 (best drop 3.9e-3). B20 pattern: every step taken at CFL > ~130 made the residual jump 4-15x (it 9->10, 16->17, 20->21, 32->33), steps at CFL 40-120 reduced it ~2x.
- The queue's next retry (cfl_max 1e3) could not help (cap above the jump threshold): stopped in its FMG start, files kept as 1025x1025_*_old_cfl1e3_stopped.txt.
- New retry added to run_queue.py (attempt 4: cfl_max 100, max_iter 1500, GMRES 200/800, FMG); j04 queue restarted (nice 19) and running it since 02:49.

## 2026-11-15 03:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_08_alpha_0/2049x2049_input.txt | it 21, drop 2.097e-09, CFL 4.73e+03, lin its 400
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 10, drop 4.575e-04, CFL 2.19e+03, lin its 181
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 12, drop 7.280e-04, CFL 1.37e+03, lin its 102
- j5015-04 (queue j04: running): results/Mach_125_alpha_0/1025x1025_input_a4.txt | it 6, drop 4.992e-02, CFL 2.00e+01, lin its 800
- alerts in the last hour:
    2026-11-15 02:47:53  NOT CONVERGED Mach_125_alpha_0 1025x1025 attempt 1: stall at iteration 84 (best drop 3.89e-03, no 10x gain in 50 iterations)

## 2026-11-15 03:10 hourly check (Claude)
- All 4 queues alive, no new alerts.
- j5015-01, M 0.8 a0 HK 2049: backtracking reached CFL 2.4e3 (4 rejected steps in all), step accepted, it 21 drop 2.1e-9 (20x gain). Converging again.
- j5015-18, M 0.8 a1.25 HK 2049: FMG done, it 10, drop 4.6e-4, CFL 2.2e3, converging.
- j5015-17, M 0.8 a1.25 VA 2049: FMG done, it 12, drop 2.4e-4 -> 7.3e-4 (CFL 4.2e3 -> 1.4e3, shock settling, CL still moving 0.3537 -> 0.3533). Watching.
- j5015-04, M 1.25 a0 HK 1025, attempt 4 (cfl_max 100): FMG done, it 6, drop 5e-2, CFL ramping (20). Too early to judge; watcher set at it 40.

## 2026-11-15 04:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_08_alpha_0/2049x2049_input.txt | it 24, drop 1.122e-09, CFL 3.782e+04, lin its linear
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 13, drop 5.391e-04, CFL 1.85e+03, lin its 163
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 16, drop 6.457e-04, CFL 1.55e+03, lin its 75
- j5015-04 (queue j04: running): results/Mach_125_alpha_0/1025x1025_input_a4.txt | it 48, drop 9.742e-03, CFL 1.00e+02, lin its 130

## 2026-11-15 04:12 hourly check + M 1.25 fix (Claude)
- j5015-04, M 1.25 a0 HK 1025: attempt 4 (cfl_max 100) bounced exactly like attempt 1 (drop 3e-3..3e-2, jumps also at CFL 100). Diagnosis agent: not GMRES (30 of 34 jumps follow a solve converged to 1e-3), not a CFL threshold, not a limit cycle of the scheme (257 / 513 converge to 1e-13): the linearized step is locally wrong at the fish-tail and bow shocks once they move by more than a cell (worse with refinement); 513 and 257 escaped only once the CFL passed ~1000. Forces already carry only ~3-5e-7 iteration error in CD (< 1% of the 257->513 change).
  Action: attempt 4 stopped (files *_old_cfl100_stopped.txt), new retry 5 = FMG + GMRES 200/800 + cfl0 10 (CFL ~1e3 at the bouncing level) + gmres_fail 1e-2 + Cauchy stop on the forces (window 20, tol 2e-7; homework criterion; status "cauchy" accepted) + max_iter 300. Running since 04:11.
- j5015-01, M 0.8 a0 HK 2049: it 24, drop 1.1e-9 (accepted / rejected cycles at CFL 2e4-4e4, GMRES 150/400 too weak there; ~19 min per step). Forces converged (CD 8.347370e-3, CL -5e-9). Left running: stopping it would lose the run (outputs are written only at the end).
- j5015-18, M 0.8 a1.25 HK 2049: it 13, drop 5.4e-4, converging.
- j5015-17, M 0.8 a1.25 VA 2049: it 16, drop 9e-5..6e-4 (shock settling), CL 0.3529. Watching.
- No new alerts since 03:10.

## 2026-11-15 04:14 2049 safety (Claude)
- Queues j01, j18, j17 restarted (solvers untouched, taken over): a taken-over 2049 run now gets a stall window of 30 iterations instead of 15. Reason: 2049 iterations take ~10-20 min and the files are only written at the end, so a stop near convergence (e.g. M 0.8 a0 on j01, forces converged, drop 1.1e-9) would lose ~8 h of work.

## 2026-11-15 05:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_08_alpha_0/2049x2049_input.txt | it 28, drop 7.484e-11, CFL 3.026e+04, lin its linear
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 17, drop 1.222e-03, CFL 8.18e+02, lin its 41
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 18, drop 3.080e-04, CFL 3.25e+03, lin its 114
- j5015-04 (queue j04: running): results/Mach_125_alpha_0/1025x1025_input_a5.txt | it 27, drop 4.124e-02, CFL 2.42e+02, lin its 100

## 2026-11-15 05:12 hourly check (Claude)
- No new alerts; all queues alive.
- j5015-01, M 0.8 a0 HK 2049: it 28, drop 7.5e-11 (from 1.1e-9 at it 24), CD 8.347370e-3, CL -1.8e-9. Converging; a few more steps.
- j5015-18, M 0.8 a1.25 HK 2049: it 17, drop bouncing 1e-4..1.2e-3, CL still drifting (0.35278 -> 0.35260): shock settling. Watching.
- j5015-17, M 0.8 a1.25 VA 2049: it 18, drop 3e-5..6e-4, CL 0.35292 -> 0.35282. Watching.
- j5015-04, M 1.25 a0 HK 1025, retry 5: it 27, still bouncing (drop 1.6e-2..7.5e-2, CFL reached ~600). Changed: the last retry is never stopped by the stall rule (a stop leaves no output files); it now ends at the Cauchy stop or max_iter 300 (~10:30) and writes its files either way. j04 queue restarted, solver taken over.

## 2026-11-15 06:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_08_alpha_0/2049x2049_input.txt | it 31, drop 6.664e-11, CFL 1.210e+04, lin its linear
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 19, drop 8.040e-04, CFL 1.24e+03, lin its 52
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 20, drop 4.748e-05, CFL 2.11e+04, lin its 185
- j5015-04 (queue j04: running): results/Mach_125_alpha_0/1025x1025_input_a5.txt | it 59, drop 1.797e-02, CFL 5.56e+02, lin its 182

## 2026-11-15 06:10 hourly check (Claude)
- No new alerts; all queues alive.
- j5015-01, M 0.8 a0 HK 2049: it 31, drop 4-7e-11 for 4 iterations (accept at CFL 3-6e3, reject at ~1.2e4: GMRES 150/400 too weak this close to machine precision; ~20 min per step). Forces converged to ~1e-10 (CD 8.347370e-3, CL -7e-10). Already below the 1e-10 the previous campaign used on 2049, but this run's target is 1e-12 and the solver writes its files only when it ends, so it is left running (estimate ~11:00-12:00). It holds up the j01 queue (M 0.5 cases): candidates to move to free lab machines this morning.
- j5015-18, M 0.8 a1.25 HK 2049: it 19, drop bouncing 7e-5..1.2e-3, CL still drifting (0.35259 -> 0.35246 over 2 it): shock still settling. Watching.
- j5015-17, M 0.8 a1.25 VA 2049: it 21, drop 3.2e-7 (breaking through, CFL 3e6).
- j5015-04, M 1.25 a0 HK 1025 retry 5: it 59, still bouncing (drop 2e-2..8e-2); CD 9.614412e-2 steady, CL 3-6e-7. Will end at max_iter 300 (~10:30) or the Cauchy stop, and write its files.

## 2026-11-15 06:50 MORNING SUMMARY (Claude)
Finished overnight (all converged to drop <= 1e-12 unless stated):
- j5015-18, M 0.8 a1.25 HK: 9 -> 1025 (513 via FMG retry). 1025: CL 0.35347803, CD 2.254034e-2, CM -3.855271e-2.
- j5015-17, M 0.8 a1.25 VA: 9 -> 1025 (513 via FMG retry). 1025: CL 0.35368382, CD 2.255202e-2, CM -3.861392e-2.
- j5015-04, M 1.25 a0 HK: 9 -> 513 (513 via FMG retry, CD 9.6143957e-2).
- j5015-01: nothing new (M 0.8 a0 HK 2049 still running).
Triplet 256-512-1024 (p / f*): M 0.8 a0 HK: CD 2.62 / 8.342797e-3. M 0.8 a1.25 HK: CL 0.93 / 0.352061, CD 1.56 / 2.249546e-2, CM 0.97 / -0.0381462. M 0.8 a1.25 VA: CL 1.10 / 0.352468, CD 1.57 / 2.250462e-2, CM 1.13 / -0.0382609 (shock: ~first order for CL, CM). No 512-1024-2048 triplet yet.
Running now:
- j5015-01 M 0.8 a0 HK 2049: it 33, drop 6.7e-11 (crawling, forces converged CD 8.347370e-3), end ~11:00-12:00.
- j5015-18 M 0.8 a1.25 HK 2049: it 19, drop 7e-5..1e-3, CL still drifting (0.3525). End uncertain (several hours).
- j5015-17 M 0.8 a1.25 VA 2049: it 21, drop 3.2e-7, breaking through: end ~1-2 h.
- j5015-04 M 1.25 a0 HK 1025 retry 5: it 73, drop ~1e-2 (bouncing), CD 9.614403e-2 steady; ends at max_iter 300 or the Cauchy stop, ~10:30-11:00, files written either way.
Idle lab machines: none yet (-17 first, ~08:00-09:00).
j01 queue if nothing moves: M 0.8 a0 2049 ~11:30, M 0.5 a0 VA ~16:00, M 0.5 a1.25 HK ~20:30, M 0.5 a1.25 VA ~01:00 Mon, M 0.8 a0 VA ~06:00 Mon, M 1.25 a0 VA ~13:00 Mon, M 1.25 a0 HK 2049 ~17:00 Mon. Not queued: M 1.25 a1.25 HK and VA.
Finding: the 2049 solvers use ~1 core (68-94 % CPU), the 1025 one 14 cores: GMRES dot products / Gram-Schmidt and the SGS sweeps are serial (src/linsolve.cpp), and at 2049 with 400-800 GMRES iterations they dominate (20-40 min per iteration).

## 2026-11-15 07:08 j5015-04 finished (Claude)
- M 1.25 a0 HK 1025, retry 5: ended by the Cauchy criterion at it 79 (forces steady within 2e-7 over 20 it; residual drop 1.35e-2): CL -4.2e-10, CD 9.6144021e-2, CM 2.3e-10; surface and field files written. 513 -> 1025 change in CD is +6.4e-8, inside the ~3-5e-7 iteration error (agent estimate), so not resolved. j04 queue done: j5015-04 is idle (nothing new started; waiting for the user).

## 2026-11-15 07:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_08_alpha_0/2049x2049_input.txt | it 34, drop 5.522e-11, CFL 4.84e+03, lin its 400
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 21, drop 2.869e-04, CFL 3.49e+03, lin its 109
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 22, drop 2.430e-08, CFL 4.116e+07, lin its linear
- j5015-04 (queue j04: stopped): no solver running 
  queue j04 finished
- alerts in the last hour:
    2026-11-15 07:07:57  === 2026-11-15 07:07:56  queue done

## 2026-11-15 07:10 hourly check (Claude)
- No new alerts. j01, j18, j17 queues alive; j04 done (idle, waiting for the user).
- j5015-01, M 0.8 a0 HK 2049: it 34, drop 5.5e-11 (crawling ~1.2-2x per accept/reject cycle; forces converged, CD 8.347370e-3). Estimate ~13:00-14:00 for 1e-12.
- j5015-18, M 0.8 a1.25 HK 2049: it 21, drop 3.6e-5..8e-4, CL drift slowing (0.35245 -> 0.35238). Converging slowly.
- j5015-17, M 0.8 a1.25 VA 2049: it 22, drop 2.4e-8; step at CFL 4.1e7 rejected (GMRES 800 not enough there), backtracking. Near the end.

## 2026-11-15 08:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_08_alpha_0/2049x2049_input.txt | it 37, drop 1.537e-11, CFL 1.94e+03, lin its 400
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 22, drop 9.110e-06, CFL 1.10e+05, lin its 800
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 24, drop 2.430e-08, CFL 4.116e+05, lin its linear
- j5015-04 (queue j04: stopped): no solver running 
  queue j04 finished

## 2026-11-15 08:10 hourly check (Claude)
- No new alerts. j01, j18, j17 queues alive; j04 idle (waiting for the user).
- j5015-01, M 0.8 a0 HK 2049: it 37, drop 1.5e-11 (from 5.5e-11), CD 8.347370e-3, CL -1.8e-10. A few more steps.
- j5015-18, M 0.8 a1.25 HK 2049: it 22, drop 9.1e-6, CL 0.352377. Converging.
- j5015-17, M 0.8 a1.25 VA 2049: it 24, drop 2.4e-8; 3 steps rejected (GMRES 800 fails at CFL 4e7, 4e6, 4e5), backtracking to ~4e4 next. Forces converged (CL 0.3527736, CD 2.250205e-2). Same accept/reject cycle as on j01: the cap is cut 10x on a failure but doubles again on each accepted step, so it keeps hitting the CFL where GMRES fails (an improvement for later: keep the cap below the last failed CFL).

## 2026-11-15 09:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_08_alpha_0/2049x2049_input.txt | it 40, drop 4.823e-12, CFL 1.549e+04, lin its linear
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 23, drop 1.805e-05, CFL 5.54e+04, lin its 800
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 25, drop 2.430e-08, CFL 4.116e+04, lin its linear
- j5015-04 (queue j04: stopped): no solver running 
  queue j04 finished

## 2026-11-15 09:10 hourly check (Claude)
- No new alerts. j01, j18, j17 queues alive; j04 idle (waiting for the user).
- j5015-01, M 0.8 a0 HK 2049: it 40, drop 4.8e-12 (1e-12 within ~1-2 accepted steps, ~1 h). CD 8.347370e-3, CL -7e-11.
- j5015-18, M 0.8 a1.25 HK 2049: it 23, drop 1.8e-5 (9e-6 at it 22), CL 0.352359 still drifting slowly; ~30-40 min per iteration (GMRES at its 800 cap).
- j5015-17, M 0.8 a1.25 VA 2049: it 25, drop 2.4e-8 since it 22: 4 steps rejected (GMRES 800 fails at CFL 4e7 down to 4e4); next try ~4e3. Forces converged (CL 0.3527736, CD 2.250205e-2). Likely cause near the end: the residual is so small that the finite-difference Jacobian-vector product of the JFNK GMRES is no longer accurate enough relative to it at high CFL; at lower CFL (diagonal-dominant) it works. Not stopped (the files are written only at the end).

## 2026-11-15 10:05 day split (user: tol 1e-10, use the lab machines, add M 1.25 alpha 1.25, done by Mon 09:00)
- tol 1e-10 on 2049 (forces converged long before 1e-12; the last two orders cost ~6 h per 2049 run). M 1.25 from 1025 up starts with the settings that ended the bouncing (cfl0 10, gmres_fail 1e-2, Cauchy 20 it / 2e-7, max_iter 250).
- A 2049 run uses ~1 core (serial GMRES / SGS), so -18 and -17 run two queues each (threads 26 / 14 for the new runs).
- j01d: M 0.8 a0 HK 2049 (taken over) -> M 0.5 a0 VA (9-2049) -> M 0.5 a1.25 HK, VA (9-1025). ~16:30.
- j18d: M 0.8 a1.25 HK 2049 (taken over) -> M 0.8 a0 VA 2049 -> M 0.5 a1.25 VA 2049. ~21:00.
- j18e: M 1.25 a0 HK 2049 -> M 1.25 a1.25 VA 2049. ~00:00.
- j17d: M 0.8 a1.25 VA 2049 (taken over) -> M 0.5 a1.25 HK 2049. ~17:00.
- j17e: M 1.25 a1.25 HK 2049 -> M 1.25 a0 VA 2049. ~00:00.
- j04d (up to 1025): M 1.25 a0 VA -> M 1.25 a1.25 HK -> M 1.25 a1.25 VA -> M 0.8 a0 VA. ~22:30.

## 2026-11-15 10:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_08_alpha_0/2049x2049_input.txt | it 43, drop 2.322e-12, CFL 6.20e+03, lin its 400
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 24, drop 1.236e-07, CFL 8.091e+06, lin its linear
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 27, drop 1.845e-08, CFL 8.23e+03, lin its 800
- j5015-04 (queue j04: running): results/Mach_125_alpha_0_vanalbada/513x513_input.txt | it 6, drop 4.043e-02, CFL 2.47e+01, lin its 22
  queue j04 finished

## 2026-11-15 10:30 hourly check (Claude)
- No alerts; all 6 queues and solvers alive.
- j01d, M 0.8 a0 HK 2049 (old input, tol 1e-12): it 44, drop 2.0-2.5e-12 for 3 iterations (GMRES at its 400 cap). It is at the noise floor of the JFNK solve; it may not reach 1e-12 soon. Decision: give it until the 11:30 check; if still above 1e-12, stop it and rerun this 2049 with tol 1e-10 (FMG + GMRES 200/800, ~5 h; j01 is not on the critical path). Forces already converged: CL -3e-11, CD 8.347370e-3.
- j18d, M 0.8 a1.25 HK 2049: it 25, drop 1.2e-7, steps rejected at CFL 8e6 / 8e5 (backtracking).
- j17d, M 0.8 a1.25 VA 2049: it 27, drop 1.8e-8, CL 0.3527734.
- j18e, M 1.25 a0 HK 2049: FMG start (513 level done). j17e, M 1.25 a1.25 HK 2049: FMG start (257 level done).
- j04d, M 1.25 a0 VA: 9-513 converged (513 CD 9.6143230e-2, 16 it); 1025 at it 5.

## 2026-11-15 11:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_05_alpha_0_vanalbada/1025x1025_input.txt | it 13, drop 8.610e-08, CFL 1.16e+07, lin its 489
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 25, drop 1.236e-07, CFL 8.091e+05, lin its linear
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 28, drop 1.821e-08, CFL 1.65e+04, lin its 800
- j5015-04 (queue j04: running): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 14, drop 1.619e-03, CFL 6.177e+03, lin its linear
  queue j04 finished

## 2026-11-15 11:35 hourly check (Claude)
- No alerts; all queues alive.
- DONE: M 0.8 a0 HK 2049 converged 10:32 (46 it, drop 2.8e-13): CL -4.6e-12, CD 8.3473703e-3, CM 3.8e-13. No rerun needed. M 0.5 a0 VA 9-1025 done; its 2049 running on j01d since ~11:00.
- j18d, M 0.8 a1.25 HK 2049 (old tol 1e-12): it 26, drop 1.2e-7 for 3 rejected steps (CFL 8e6 -> 8e4); slow. To protect Mon 09:00, M 0.5 a1.25 VA 2049 moved from j18d to the end of j01d (j01d, j18d restarted, solvers taken over).
- j17d, M 0.8 a1.25 VA 2049: it 28, drop 1.8e-8, ~40 min per iteration.
- j18e, M 1.25 a0 HK 2049: FMG start, 1025 level. j17e, M 1.25 a1.25 HK 2049: FMG start, 513 level done after 3430 s (39 it). Watch: on 1025 the M 1.25 cases bounce and end by the Cauchy stop (~80 it), so these FMG levels and then the 2049 grids may be long; they are the critical path.
- j04d, M 1.25 a0 VA 1025: it 23, drop 3.6e-3 at CFL 2.8e3 (higher than the HK run reached).

## 2026-11-15 12:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_05_alpha_0_vanalbada/2049x2049_input.txt | it 9, drop 4.795e-04, CFL 2.09e+03, lin its 265
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 26, drop 1.236e-07, CFL 8.091e+04, lin its linear
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 29, drop 7.771e-09, CFL 3.29e+04, lin its 800
- j5015-04 (queue j04: running): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 37, drop 7.258e-03, CFL 1.38e+03, lin its 183
  queue j04 finished

## 2026-11-15 12:35 hourly check + M 1.25 rescheduling (Claude)
- No alerts. j01d: M 0.5 a0 VA 2049 it 9 (drop 4.8e-4). j18d: M 0.8 a1.25 HK 2049 it 27 (1.2e-7, backtracking). j17d: M 0.8 a1.25 VA 2049 it 30 (2.3e-9). j04d: M 1.25 a0 VA 1025 it 44 (bouncing; CL oscillates ~3e-6).
- Problem: with Cauchy 20 it / 2e-7 and max_iter 250 the M 1.25 runs on 1025 would not stop by Cauchy (CL ~3e-6 oscillation) and a "max_iter" end would have started retries; on 2049, 250 iterations take days.
- New M 1.25 rule (1025 / 2049): Cauchy 10 it / 1e-6, max_iter 150 (1025) / 60 (2049), one attempt, no stall stop, "max_iter" accepted as finished (forces then good to ~1e-6; CD at M 1.25 is grid-converged to ~1e-7 by 513).
- M 1.25 a0 HK 2049 (j18e) and M 1.25 a1.25 HK 2049 (j17e) were still in their FMG start: stopped (files *_old_cauchy2e-7_stopped.txt) and restarted with the new rule (12:31). j04d restarted (its 1025 run taken over under the new rule).
- The second M 1.25 2049 of j18e / j17e now runs in parallel as queues j18f (M 1.25 a1.25 VA 2049, threads 17) and j17f (M 1.25 a0 VA 2049, threads 9): all four M 1.25 2049 runs are going since 12:33 (memory: -18 209 GB free, -17 109 GB before).
- Estimates: M 1.25 2049 runs: FMG ~1-2 h + up to 60 iterations (~10-15 min each) -> ~23:00 to ~04:00 worst case. M 0.8 a1.25 HK / VA 2049 then M 0.8 a0 VA 2049 (j18d) and M 0.5 a1.25 HK 2049 (j17d): ~22:00-02:00. j01d: M 0.5 a0 VA 2049, M 0.5 a1.25 HK / VA to 1025, M 0.5 a1.25 VA 2049: ~20:00. j04d (to 1025): M 1.25 a0 VA, a1.25 HK, a1.25 VA, then M 0.8 a0 VA: ~21:00. Everything before Mon 09:00 with a margin.

## 2026-11-15 13:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_05_alpha_0_vanalbada/2049x2049_input.txt | it 10, drop 5.413e-05, CFL 1.85e+04, lin its 759
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 27, drop 1.236e-07, CFL 8.09e+03, lin its 800
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 30, drop 2.295e-09, CFL 6.58e+04, lin its 800
- j5015-04 (queue j04: running): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 56, drop 6.890e-03, CFL 1.45e+03, lin its 232
  queue j04 finished

## 2026-11-15 13:19 hourly check (Claude)
- No alerts; all 8 queues and solvers alive.
- j01d M 0.5 a0 VA 2049: it 10 (drop 5.4e-5), the next iteration has taken >45 min (GMRES 759 its, serial Gram-Schmidt).
- j18d M 0.8 a1.25 HK 2049: it 27 (1.2e-7), slow. j17d M 0.8 a1.25 VA 2049: it 30 (2.3e-9), slow (~40-50 min per iteration; three 2049 runs share each lab machine's memory bandwidth).
- M 1.25 2049 x4 (j18e, j18f, j17e, j17f): FMG start, 513 level done (44-60 it), now on the 1025 level.
- j04d M 1.25 a0 VA 1025: it 58, drop 4.2e-3, CD 9.61453e-2 (ends by Cauchy or max_iter 250, ~14:45).
- Risk for Mon 09:00: the M 1.25 2049 grids (60 iterations max) at 10-20 min per iteration -> ~00:00-10:00. Re-estimate once they reach the 2049 level; the fix that would help most (parallel Gram-Schmidt / dot products, option A) needs the user's go-ahead.

## 2026-11-15 14:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_05_alpha_0_vanalbada/2049x2049_input.txt | it 11, drop 1.354e-06, CFL 7.39e+05, lin its 800
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 28, drop 2.549e-07, CFL 1.62e+04, lin its 800
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 31, drop 4.897e-10, CFL 1.32e+05, lin its 800
- j5015-04 (queue j04: running): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 74, drop 7.662e-04, CFL 1.31e+04, lin its 345
  queue j04 finished

## 2026-11-15 14:20 hourly check (Claude)
- No alerts; all queues alive.
- j01d M 0.5 a0 VA 2049: it 12, drop 4.8e-8 (tol 1e-10, ~1-2 h left).
- j18d M 0.8 a1.25 HK 2049 (old tol 1e-12): it 29, drop 1.0e-7, crawling (~45 min per iteration). Its next case, M 0.8 a0 VA 2049, now runs in parallel as queue j18g (started 14:20, threads 13; -18 had 172 GB free).
- j17d M 0.8 a1.25 VA 2049 (old tol 1e-12): it 31, drop 4.9e-10. Next on j17d: M 0.5 a1.25 HK 2049.
- M 1.25 2049: a0 HK (j18e) on the 2049 level (it 2); a1.25 HK (j17e), a1.25 VA (j18f), a0 VA (j17f) still on the 1025 FMG level.
- j04d M 1.25 a0 VA 1025: it 75, drop 2.1e-3 (to max_iter 250, ~17:30), then M 1.25 a1.25 HK / VA and M 0.8 a0 VA up to 1025 (~00:30).

## 2026-11-15 15:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_05_alpha_0_vanalbada/2049x2049_input.txt | it 13, drop 1.175e-09, CFL 1.00e+08, lin its 800
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 29, drop 1.032e-07, CFL 3.24e+04, lin its 800
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 32, drop 1.443e-10, CFL 2.634e+05, lin its linear
- j5015-04 (queue j04: running): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 87, drop 2.626e-02, CFL 3.81e+02, lin its 133
  queue j04 finished

## 2026-11-15 15:35 check (Claude)
- j5015-01 j01d: M 0.5 a0 VA 2049 it 13, drop 1.2e-9 (~50 min/it, 1e-10 expected ~17:00). CD 2.72191e-6.
- j5015-18: j18d M 0.8 a1.25 HK 2049 it 30, drop 7.2e-9 (tol 1e-12, ~1 h/it; CL 0.352356 CD 2.248006e-2 CM -3.82226e-2 steady to 6-7 digits); j18e M 1.25 a0 HK 2049 fine level it 4 (drop 9e-2, ~1 h/it); j18f M 1.25 a1.25 VA 2049 on its 1025 FMG level; j18g M 0.8 a0 VA 2049 it 10 drop 4e-4.
- j5015-17: j17d M 0.8 a1.25 VA 2049 it 32, drop 1.4e-10 (tol 1e-12; step rejected at CFL 2.6e5, cap re-doubling as before; M 0.8 a0 HK 2049 reached 2.8e-13 the same way); j17e M 1.25 a1.25 HK 2049 fine level it 2; j17f M 1.25 a0 VA 2049 on its 1025 FMG level.
- j5015-04 j04d: M 1.25 a0 VA 1025 it 92/250 (~150 s/it, old rule 20 it/2e-7 will not trigger: CL varies ~4e-6), ends ~22:00.
- Re-estimate: j04d still had M 1.25 a1.25 HK + VA (to 1025) and M 0.8 a0 VA (to 1025) behind it -> ~Mon 15:00, and M 0.5 a1.25 HK 2049 waited behind j17d (1e-12 run, open end). Action (15:30, backup results/.run_queue_backup_1530.py): new parallel queues
  j18h = M 0.5 a1.25 HK 2049 (threads 13), j18i = M 1.25 a1.25 HK 9-1025 (threads 8),
  j17h = M 1.25 a1.25 VA 9-1025 (threads 8), j17i = M 0.5 a1.25 HK 9-1025 (threads 8);
  removed from j01d / j17d / j04d (restarted, solvers taken over). j04d keeps M 0.8 a0 VA 9-1025 (~01:00). Every grid in exactly one queue (checked). Memory after the FMG starts grow: -18 ~75 GB free, -17 ~30 GB free.
- Expected ends: j01d (then M 0.5 a1.25 VA 9-2049) ~01:00; j18h ~22:00; j17h/j18i (1025 up to 150 it) ~23:00-01:00; j17i ~19:00; M 1.25 2049 runs (60 it, ~1 h/it on fine level) ~Mon 02:00-08:00; j17d/j18d (1e-12) open, forces already steady.

## 2026-11-15 16:09 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_05_alpha_125_vanalbada/1025x1025_input.txt | it 11, drop 2.209e-04, CFL 4.53e+03, lin its 308
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 30, drop 7.199e-09, CFL 6.47e+04, lin its 800
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 33, drop 1.443e-10, CFL 2.63e+04, lin its 800
- j5015-04 (queue j04: running): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 101, drop 9.785e-04, CFL 3.42e+03, lin its 537
  queue j04 finished

## 2026-11-15 17:10 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_05_alpha_125_vanalbada/2049x2049_input.txt | it 7, drop 8.161e-03, CFL 1.23e+02, lin its 27
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 30, drop 7.199e-09, CFL 6.47e+04, lin its 800
- j5015-17 (queue j17: running): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 33, drop 1.443e-10, CFL 2.63e+04, lin its 800
- j5015-04 (queue j04: running): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 118, drop 2.877e-03, CFL 3.48e+03, lin its 800
  queue j04 finished
- alerts in the last hour:
    2026-11-15 17:05:30  === 2026-11-15 17:05:30  queue done

## 2026-11-15 18:10 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_05_alpha_125_vanalbada/2049x2049_input.txt | it 10, drop 6.944e-05, CFL 1.44e+04, lin its 588
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 31, drop 5.194e-09, CFL 1.29e+05, lin its 800
- j5015-17 (queue j17: stopped): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 34, drop 9.155e-11, CFL 5.27e+04, lin its 800
- j5015-04 (queue j04: stopped): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 125, drop 2.350e-03, CFL 4.25e+03, lin its 800
  queue j04 finished

## 2026-11-15 19:10 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_05_alpha_125_vanalbada/2049x2049_input.txt | it 11, drop 4.445e-06, CFL 2.25e+05, lin its 800
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 31, drop 5.194e-09, CFL 1.29e+05, lin its 800
- j5015-17 (queue j17: stopped): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 35, drop 2.447e-11, CFL 1.05e+05, lin its 800
- j5015-04 (queue j04: stopped): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 134, drop 4.722e-03, CFL 2.12e+03, lin its 393
  queue j04 finished

## 2026-11-15 20:10 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_05_alpha_125_vanalbada/2049x2049_input.txt | it 12, drop 6.767e-07, CFL 1.48e+06, lin its 800
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 32, drop 3.094e-10, CFL 2.59e+05, lin its 800
- j5015-17 (queue j17: stopped): results/Mach_08_alpha_125_vanalbada/2049x2049_input.txt | it 35, drop 2.447e-11, CFL 1.05e+05, lin its 800
- j5015-04 (queue j04: stopped): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 151, drop 3.138e-02, CFL 3.19e+02, lin its 110
  queue j04 finished

## 2026-11-15 21:10 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_05_alpha_125_vanalbada/2049x2049_input.txt | it 13, drop 1.404e-09, CFL 1.00e+08, lin its 800
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 32, drop 3.094e-10, CFL 2.59e+05, lin its 800
- j5015-17 (queue j17: stopped): results/Mach_125_alpha_125/2049x2049_input.txt | it 12, drop 1.592e-01, CFL 6.28e+01, lin its 99
- j5015-04 (queue j04: stopped): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 167, drop 3.567e-02, CFL 2.80e+02, lin its 92
  queue j04 finished

## 2026-11-15 22:10 automatic status (post/night_status.sh)
- j5015-01 (queue j01: stopped): no solver running 
  queue j01 finished
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 33, drop 1.403e-10, CFL 5.18e+05, lin its 800
- j5015-17 (queue j17: stopped): results/Mach_125_alpha_125/2049x2049_input.txt | it 14, drop 4.483e-02, CFL 2.23e+02, lin its 108
- j5015-04 (queue j04: stopped): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 178, drop 1.718e-03, CFL 5.822e+03, lin its linear
  queue j04 finished
- alerts in the last hour:
    2026-11-15 21:57:57  === 2026-11-15 21:57:57  queue done

## 2026-11-15 23:10 automatic status (post/night_status.sh)
- j5015-01 (queue j01: stopped): no solver running 
  queue j01 finished
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 33, drop 1.403e-10, CFL 5.18e+05, lin its 800
- j5015-17 (queue j17: stopped): results/Mach_125_alpha_125/2049x2049_input.txt | it 16, drop 1.173e-02, CFL 8.53e+02, lin its 343
- j5015-04 (queue j04: stopped): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 195, drop 6.051e-03, CFL 1.65e+03, lin its 800
  queue j04 finished

## 2026-11-16 00:05 check (Claude)
- Found the queue processes on j5015-17 (j17d/e/f/h) and j5015-04 (j04d) dead since the evening (no log line; the solvers kept running, nice 19). Restarted j04d, j17e, j17f, j17h (solvers taken over); j17d: its 2049 had converged -> queue done. M 0.8 a0 VA 9-1025 moved from j04d to a new queue j01e on j5015-01 (idle, started 00:00).
- Finished since 15:35: M 0.5 a0 VA 2049 (CD 2.7219e-6, 6.1e-11); M 0.5 a1.25 HK 9-1025 (j17i); M 0.5 a1.25 VA 9-2049 (2049 CL 0.17936815 CD 2.35589e-5 CM -2.25301e-3); M 0.8 a1.25 VA 2049 (CL 0.35277251 CD 2.2501972e-2 CM -3.8344877e-2, 9.4e-13, 20:26).
- Running: j18d M 0.8 a1.25 HK 2049 it 33 drop 1.4e-10 (~2 h/it, ~4 it to 1e-12 -> ~08:00); j18g M 0.8 a0 VA 2049 it 14 1.5e-8 (~04:00); j18h M 0.5 a1.25 HK 2049 it 12 3.1e-6 (~06:00-08:00); j18e M 1.25 a0 HK 2049 it 18/60 (~10 min-1 h/it; Cauchy or ~10:00); j17e M 1.25 a1.25 HK 2049 it 20/60 (~07:00); j17f M 1.25 a0 VA 2049 it 13/60 (up to 50 min/it, Mon ~19:00 unless Cauchy); j18f M 1.25 a1.25 VA 2049 still on its 1025 FMG level (60 it at ~10 min) -> fine level 60 it: Mon evening / Tue; j18i M 1.25 a1.25 HK 1025 it 50/150 (~02:30); j17h M 1.25 a1.25 VA 1025 it 54/150 (~8 min/it -> ~13:00); j04d M 1.25 a0 VA 1025 it 210/250 (~01:20); j01e (~03:00).
- Iterations are serial-bound (each solver ~1.2-1.7 cores, GMRES Gram-Schmidt), not CPU contention. The M 1.25 Van Albada runs need many more GMRES iterations per step than Hemker-Koren.

## 2026-11-16 00:11 automatic status (post/night_status.sh)
- j5015-01 (queue j01: running): results/Mach_08_alpha_0_vanalbada/1025x1025_input.txt | it 11, drop 1.719e-03, CFL 5.82e+02, lin its 69
  queue j01 finished
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 34, drop 3.298e-11, CFL 1.04e+06, lin its 800
- j5015-17 (queue j17: running): results/Mach_125_alpha_125/2049x2049_input.txt | it 21, drop 4.580e-02, CFL 2.18e+02, lin its 149
- j5015-04 (queue j04: running): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 212, drop 1.321e-03, CFL 7.57e+03, lin its 332
  queue j04 finished
- alerts in the last hour:
    2026-11-16 00:00:55  === 2026-11-16 00:00:55  queue done

## 2026-11-16 00:17 (Claude)
- User: M 1.25 Van Albada 2049 runs not needed. Stopped M 1.25 a0 VA 2049 (j17f, was it 13/60) and M 1.25 a1.25 VA 2049 (j18f, still on its 1025 FMG level); queues j17f, j18f ended, cases removed from run_queue.py; terminal/input kept as 2049x2049_*_old_killed_by_user.txt. M 1.25 VA grid studies end at 1025 (j04d a0 1025, j17h a1.25 1025 still running).

## 2026-11-16 00:27 check (Claude)
- No new alerts; all queues alive (j17d done; j17f/j18f closed by the user's request).
- j01e M 0.8 a0 VA 1025 it 14 drop 1.1e-5 (~01:00). j18g M 0.8 a0 VA 2049 it 15 drop 3.2e-9 (~02:30). j18d M 0.8 a1.25 HK 2049 it 34 drop 3.3e-11 (1e-12 ~04:00). j18h M 0.5 a1.25 HK 2049 it 13 drop 4.4e-7 (CD 2.36899e-5; ~05:00-06:00).
- M 1.25: j04d a0 VA 1025 it 215/250 (~240 s/it, ~02:50); j18e a0 HK 2049 it 18/60; j17e a1.25 HK 2049 it 23/60 (~07:00-10:00). j18i a1.25 HK 1025 it 52 and j17h a1.25 VA 1025 it 57 of 150: bouncing at drop 1e-2..3e-2 (B20), ~450 s/it; forces within ~3e-6 (CL), 3e-7 (CD), 4e-7 (CM), so the Cauchy stop (1e-6) is unlikely; both end at max_iter 150 ~12:00-13:00 (was 02:30 for j18i).
- 00:31: j18i M 1.25 a1.25 HK 1025 ended by the Cauchy stop at it 55: CL 8.9360688e-02 CD 9.7926946e-02 CM -1.8922425e-02 (queue done). j17h (VA 1025) may end the same way.
- 00:45: j01e M 0.8 a0 VA 9-1025 done (1025 CD 8.3475564e-03, 4.3e-13). j5015-01 idle; no case left to start (every remaining grid is already running).

## 2026-11-16 01:11 automatic status (post/night_status.sh)
- j5015-01 (queue j01: stopped): no solver running 
  queue j01 finished
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 34, drop 3.298e-11, CFL 1.04e+06, lin its 800
- j5015-17 (queue j17: running): results/Mach_125_alpha_125/2049x2049_input.txt | it 27, drop 1.115e-01, CFL 8.97e+01, lin its 65
- j5015-04 (queue j04: running): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 226, drop 2.569e-02, CFL 3.89e+02, lin its 124
  queue j04 finished
- alerts in the last hour:
    2026-11-16 00:45:07  === 2026-11-16 00:45:07  queue done
    2026-11-16 00:31:26  === 2026-11-16 00:31:26  queue done

## 2026-11-16 02:11 automatic status (post/night_status.sh)
- j5015-01 (queue j01: stopped): no solver running 
  queue j01 finished
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 35, drop 3.845e-10, CFL 2.071e+06, lin its linear
- j5015-17 (queue j17: running): results/Mach_125_alpha_125/2049x2049_input.txt | it 30, drop 1.336e-02, CFL 7.49e+02, lin its 800
- j5015-04 (queue j04: running): results/Mach_125_alpha_0_vanalbada/1025x1025_input.txt | it 241, drop 3.047e-03, CFL 2.823e+03, lin its linear
  queue j04 finished
- alerts in the last hour:
    2026-11-16 01:58:52  === 2026-11-16 01:58:52  queue done

## 2026-11-16 03:11 automatic status (post/night_status.sh)
- j5015-01 (queue j01: stopped): no solver running 
  queue j01 finished
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 36, drop 3.845e-10, CFL 2.07e+05, lin its 800
- j5015-17 (queue j17: running): results/Mach_125_alpha_125/2049x2049_input.txt | it 35, drop 1.261e-02, CFL 7.93e+02, lin its 197
- j5015-04 (queue j04: stopped): no solver running 
  queue j04 finished
- alerts in the last hour:
    2026-11-16 02:32:51  === 2026-11-16 02:32:51  queue done
    2026-11-16 03:02:21  === 2026-11-16 03:02:21  queue done

## 2026-11-16 04:11 automatic status (post/night_status.sh)
- j5015-01 (queue j01: stopped): no solver running 
  queue j01 finished
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 37, drop 1.384e-10, CFL 4.14e+05, lin its 800
- j5015-17 (queue j17: running): results/Mach_125_alpha_125/2049x2049_input.txt | it 40, drop 2.940e-02, CFL 3.40e+02, lin its 162
- j5015-04 (queue j04: stopped): no solver running 
  queue j04 finished

## 2026-11-16 05:11 automatic status (post/night_status.sh)
- j5015-01 (queue j01: stopped): no solver running 
  queue j01 finished
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 38, drop 7.441e-11, CFL 8.29e+05, lin its 800
- j5015-17 (queue j17: running): results/Mach_125_alpha_125_vanalbada/1025x1025_input.txt | it 106, drop 2.461e-03, CFL 4.06e+03, lin its 800
- j5015-04 (queue j04: stopped): no solver running 
  queue j04 finished
- alerts in the last hour:
    2026-11-16 04:39:20  === 2026-11-16 04:39:20  queue done

## 2026-11-16 06:11 automatic status (post/night_status.sh)
- j5015-01 (queue j01: stopped): no solver running 
  queue j01 finished
- j5015-18 (queue j18: running): results/Mach_08_alpha_125/2049x2049_input.txt | it 39, drop 7.911e-11, CFL 1.66e+06, lin its 800
- j5015-17 (queue j17: running): results/Mach_125_alpha_125_vanalbada/1025x1025_input.txt | it 116, drop 2.924e-03, CFL 3.421e+03, lin its linear
- j5015-04 (queue j04: stopped): no solver running 
  queue j04 finished

## 2026-11-16 06:55 check (Claude)
- Finished since 00:45: M 0.8 a0 VA 2049 (01:58, CD 8.3470348e-3); M 1.25 a0 VA 1025 (02:32, max_iter 250: CD 9.6143239e-2, CL -2.8e-5 = unconverged noise of the bouncing run, exact 0); M 0.5 a1.25 HK 2049 (03:02, CL 0.17935770 CD 2.3688567e-5 CM -2.2507155e-3); M 1.25 a1.25 HK 2049 (04:38, Cauchy: CL 8.9355985e-2 CD 9.7931002e-2 CM -1.8920155e-2); M 0.8 a1.25 HK 2049 (06:35, 1.2e-13: CL 0.35235550 CD 2.2480045e-2 CM -3.8222514e-2).
- Left: j18e M 1.25 a0 HK 2049 it 38/60 (CL within ~1e-6, CD 3e-7 over the last 8 it: Cauchy may stop it soon, else max_iter ~13:00); j17h M 1.25 a1.25 VA 1025 it 124/150 (~10:15). Every other queue done.

## 2026-11-16 07:11 automatic status (post/night_status.sh)
- j5015-01 (queue j01: stopped): no solver running 
  queue j01 finished
- j5015-18 (queue j18: running): results/Mach_125_alpha_0/2049x2049_input.txt | it 38, drop 4.999e-02, CFL 2.00e+01, lin its 53
- j5015-17 (queue j17: running): results/Mach_125_alpha_125_vanalbada/1025x1025_input.txt | it 129, drop 5.741e-03, CFL 1.74e+03, lin its 386
- j5015-04 (queue j04: stopped): no solver running 
  queue j04 finished
- alerts in the last hour:
    2026-11-16 06:36:15  === 2026-11-16 06:36:15  queue done

## 2026-11-16 08:11 automatic status (post/night_status.sh)
- j5015-01 (queue j01: stopped): no solver running 
  queue j01 finished
- j5015-18 (queue j18: running): results/Mach_125_alpha_0/2049x2049_input.txt | it 39, drop 4.675e-02, CFL 4.00e+01, lin its 800
- j5015-17 (queue j17: running): results/Mach_125_alpha_125_vanalbada/1025x1025_input.txt | it 140, drop 5.900e-03, CFL 1.69e+03, lin its 597
- j5015-04 (queue j04: stopped): no solver running 
  queue j04 finished

## 2026-11-16 09:11 automatic status (post/night_status.sh)
- j5015-01 (queue j01: stopped): no solver running 
  queue j01 finished
- j5015-18 (queue j18: stopped): no solver running 
- j5015-17 (queue j17: stopped): no solver running 
- j5015-04 (queue j04: stopped): no solver running 
  queue j04 finished
- alerts in the last hour:
    2026-11-16 08:41:32  === 2026-11-16 08:41:32  queue done
    2026-11-16 08:23:20  === 2026-11-16 08:23:20  queue done

## 2026-11-16 10:11 automatic status (post/night_status.sh)
- j5015-01 (queue j01: stopped): no solver running 
  queue j01 finished
- j5015-18 (queue j18: stopped): no solver running 
- j5015-17 (queue j17: stopped): no solver running 
- j5015-04 (queue j04: stopped): no solver running 
  queue j04 finished

## 2026-11-16 11:11 automatic status (post/night_status.sh)
- j5015-01 (queue j01: stopped): no solver running 
  queue j01 finished
- j5015-18 (queue j18: stopped): no solver running 
- j5015-17 (queue j17: stopped): no solver running 
- j5015-04 (queue j04: stopped): no solver running 
  queue j04 finished

## 2026-11-16 12:11 automatic status (post/night_status.sh)
- j5015-01 (queue j01: stopped): no solver running 
  queue j01 finished
- j5015-18 (queue j18: stopped): no solver running 
- j5015-17 (queue j17: stopped): no solver running 
- j5015-04 (queue j04: stopped): no solver running 
  queue j04 finished
