# 6. Locked software validation

Use this stage when the question is whether the already-trained BiLab model can
score another lab's SEG files without refitting.

```sh
./mbpr inspect --segments cohort.seg --format auto --genome-build GRCh37
./mbpr score --segments cohort.seg --format auto --genome-build GRCh37 \
  --output-dir scores
./mbpr validate --predictions scores/predictions.csv --outcomes outcomes.csv \
  --output-dir validation
```

The model artifact is checksum-locked. `score` emits raw 3-/5-/7-year scores,
derivation-cohort percentiles, provisional probabilities, features, QC, a run
manifest, and a local HTML report. External outcome validation does not modify
coefficients or calibration.

This model is a second, complementary story after hotlist replication: the
hotlist tests a compact biological signal; the locked software tests deployable
prediction from full SEG data.
