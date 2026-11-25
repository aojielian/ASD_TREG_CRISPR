# Gene-property analysis

The [release workflow](../../docs/REPRODUCIBILITY.md) includes this module. A standalone run uses the extracted v1.2.0 companion data or a verified v1.1.0 input root.

```sh
python3 scripts/gene_property_adjustment/reproduce.py --data-dir DATA --out-dir OUT --verify-r
```

The module reconstructs expression detection and strata, recalculates all 15 Fisher backgrounds and 30 exact conditional comparisons, refits 33 property models, reruns the historical Firth bridge, fixed-caliper matching and candidate checks, and verifies primary CPM ≥1 M2 fits with an independent base-R objective when requested. See the [locked methods](../../config/gene_property_methods.md) and [numerical fixture](../../config/gene_property_expected_results.tsv).
