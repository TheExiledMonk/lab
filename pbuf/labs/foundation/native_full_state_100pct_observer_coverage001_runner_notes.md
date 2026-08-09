# Runner notes

Run only the committed audit:

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_full_state_100pct_observer_coverage001.py
```

This is expected to be materially heavier than the 25% labs: the 100% lane uses 532×532 = 283,024 rays per cluster, exactly four times the 25% ray count, and runs all five canonical clusters.

Do not alter ray counts, coverage, decoder inventory, coefficients, propagation settings, or output handling if runtime is long. Return the complete stdout/stderr and exit code unchanged.
