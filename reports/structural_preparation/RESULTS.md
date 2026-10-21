# Uniform structural-preparation diagnostic

Twelve previously selected models were processed with one frozen protocol: both long-parent poses, both failed double-mutant poses, six matched controls, and both short-parent poses. None was replaced by a different prediction sample.

Amber14/GBn2 minimization holds the receptor fixed, restrains the binder backbone, and permits at most 500 iterations with a 180-second process limit. The declared movement limits are 0.5 Å binder backbone RMS displacement, 1 Å maximum backbone-atom displacement, and 0.01 Å maximum movement of an original receptor atom. Distances are measured in the original coordinate frame, without a new alignment.

All twelve preparations completed. The long parent's seed-809 Lys205–Glu66 close pair resolves, with binder backbone RMS displacement 0.141 Å and maximum displacement 0.868 Å. This is consistent with a local side-chain geometry problem; it does not establish binding or retrospectively change the original failed screen.

The failed double mutant remains misplaced: its prepared models still have more than 1,000 atom-pair clashes with omitted same-subunit receptor regions. The short parent's seed-607 full-context receptor fit remains 2.0457 Å because receptor coordinates are fixed. Preparation does not rescue either failure.

Eleven models satisfy the displacement limits. The unrelated VHH negative at seed 809 has a maximum backbone-atom movement of 1.318 Å and fails that diagnostic guard. Its raw negative-control result is retained. No preparation limit was relaxed.

For the next round, identical preparation is applied once to every selected model, including controls and failed candidates. Raw prediction scores and decisions remain primary; prepared structures support geometry and pKa diagnostics only. This does not model an experimental pH transition, an unbound ensemble, or full-receptor electrostatics.
