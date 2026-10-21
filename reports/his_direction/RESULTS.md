# A41H gives a consistent release-direction hypothesis, but misses one strict heuristic

The two variants passing both selected raw structural checks were evaluated on their independently predicted interfaces. All four prepared models retain the original pose/contact/context criteria and satisfy the frozen movement limits. No histidines were introduced into these models after prediction: PROPKA scores the actual independently predicted mutant sequences.

| Variant | Prediction seed | Free His pKa(s) | Bound His pKa(s) | Bound minus free | His-only Kd(6.0)/Kd(7.4) proxy | Strict direction rule |
|---|---:|---|---|---|---:|---|
| Short A41H | 1103 | 6.406 | 5.500 | −0.905 | 2.48 | Pass |
| Short A41H | 1201 | 7.302 | 6.183 | −1.118 | 4.91 | Fail: free pKa above 7 |
| Long L7H/L82H | 1103 | 6.343, 6.439 | 6.252, 5.976 | −0.091, −0.463 | 2.04 | Fail: shifts too small |
| Long L7H/L82H | 1201 | 7.372, 6.841 | 7.181, 8.311 | −0.191, +1.469 | 0.35 | Fail: inconsistent direction |

**No variant clears the strict direction rule in both poses.** That result and the empty strict-pass FASTA are preserved.

A41H nevertheless has a consistent predicted pKa decrease of about one unit and a release-direction proxy above two in both poses. Its second-model failure is specifically the free-pKa 6–7 window, not loss of the interface, reversal of the pH direction, a clash, or proximity to artificial termini. That window is a screening heuristic rather than a physical boundary. A41H is therefore worth retaining as an explicitly qualified computational hypothesis; its failed gate is not retroactively changed. It also has the strongest structural consistency among the tested variants, passing 9/10 fresh confirmation samples. Neither unmodified parent passes both selected raw checks in that confirmation, so this does not establish improved affinity or a successful matched-parent/variant qualification.

L7H/L82H has both poorer sample consistency (3/10 raw structural passes) and an inconsistent pKa direction. It should not be prioritized as a pH-switch hypothesis.

Distances from the engineered histidines to artificial receptor termini are 24.5–25.3 Å for A41H and 29.3–29.7 Å for L7H/L82H. An atom-distance audit places the A41H ring near target Leu212, and in one pose Asn215, rather than establishing a direct Arg208-repulsion mechanism. These contacts alone do not identify the energetic cause of the predicted pKa shift.

The numerical ratio is a fixed-conformation, independent-His proxy. It is not a predicted experimental affinity ratio: free binder uses the bound conformation, the receptor is cropped, full proton linkage is omitted, and PROPKA/model uncertainty is not calibrated here. No experimental human binding, mouse cross-reactivity, acidic release, or complete on/off switch is established. The intended challenge compares pH 7.4 with 6.0; these calculations address that interval only as a direction hypothesis.

## Species and tag context

The existing sequence map marks nearby human residues 212 and 215 as different in mouse; position 212 also carries an alignment-gap correspondence caution. This flags a cross-reactivity risk, not a prediction of mouse nonbinding. A mouse comparison requires an audited structural correspondence rather than blindly relabeling the human loop.

The His41 ring is approximately 12.2–12.4 Å from the binder N terminus and 17.1 Å from its C-terminal carboxyl atoms. The binder terminal carbonyl carbon is about 17.8–18.0 Å from the modeled receptor crop. These distances do not model the assay's C-terminal tag, full glycan ensemble, or attachment geometry and therefore do not establish tag compatibility.
