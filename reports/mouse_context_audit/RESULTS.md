# Mouse target correspondence audit

The [public AlphaFold Q62351 model](https://alphafold.ebi.ac.uk/entry/Q62351) was downloaded from the official API, with its URL and SHA-256 recorded. Its 763-residue sequence exactly matches the project's mouse UniProt sequence. This is a predicted monomer, not an experimental mouse TfR1–binder complex.

Using 125 conserved apical Cα anchors, excluding human loop residues 198–218 and low-confidence mouse anchors, the mouse model fits the bound human 6WRW crop with RMSD **0.892 Å**. A provisional mouse crop, native residues **192–385**, has 194 residues and is saved with explicit native/crop mapping. It has not been used for a paid binding prediction.

The local alignment caution is substantive. In the conserved-core-fitted frame, human Tyr211 Cα lies 0.729 Å from mouse Asp213 Cα, and human Leu212 lies 0.586 Å from mouse Pro214. These support the existing sequence-alignment correspondence more closely than the alternative published-position labels in the supplied brief. Human Val210 has no aligned mouse residue in that sequence alignment; its nearest mouse Cα is Leu212 at 1.897 Å. A nearest-atom correspondence is not itself an evolutionary alignment or proof of a bound conformation, so no original map was silently rewritten.

A41H's nearby human receptor residues include 212 and 215, both nonconserved in the existing map. That flags a cross-reactivity risk rather than proving mouse nonbinding. A future mouse prediction should use the actual mouse sequence and coordinates, retain insertion/deletion accounting, and include target-conformation/control checks. Human confidence thresholds and a human-bound template cannot establish a calibrated mouse-binding classifier.

This audit provides target-preparation groundwork. No mouse affinity, mouse binding, or cross-reactivity result is claimed.
