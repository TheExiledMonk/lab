# PBUF-FOUNDATION-ENDPOINT-INTERFACE-ROLE-AUDIT-001

**Head:** `af1e48bef95ca74df7fb8ad75dff333c9e85bdcd`

**Outcome:** Outcome A — INTERFACE FIELD IS COORDINATE-SAFE; ENDPOINT FIELD IS ORIENTATION-DEPENDENT

| RC | reversed pairs | E response even | E response odd | E endpoint | E interface | E LOS endpoint | E LOS interface | endpoint closure | interface pairs ok |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| RC0 | 0 | 0.000e+00 | 2.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 1.409e-17 | True |
| RC1 | 0 | 5.681e-16 | 2.000e+00 | 8.578e-16 | 8.993e-16 | 4.069e-16 | 4.715e-16 | 7.000e-18 | True |
| RC2 | 0 | 8.519e-16 | 2.000e+00 | 1.240e-15 | 1.395e-15 | 5.384e-16 | 6.182e-16 | 1.761e-18 | True |
| RC3 | 0 | 7.338e-16 | 2.000e+00 | 1.066e-15 | 1.205e-15 | 4.232e-16 | 4.831e-16 | 2.417e-18 | True |
| RC4 | 32768 | 8.038e-16 | 2.000e+00 | 1.216e+00 | 1.330e-15 | 4.243e-16 | 5.135e-16 | 2.311e-18 | True |
| RC5 | 36288 | 8.386e-16 | 2.000e+00 | 1.147e+00 | 1.374e-15 | 1.459e+00 | 5.799e-16 | 2.459e-18 | True |
| RC6 | 36288 | 5.681e-16 | 2.000e+00 | 1.149e+00 | 8.993e-16 | 1.459e+00 | 4.715e-16 | 1.455e-18 | True |

## Role classification

- pair response orientation-even covariance: `True`
- endpoint orientation dependence pattern: `True`
- interface 3D covariance: `True`
- interface native-z LOS covariance after back-transform: `True`
- endpoint role candidate: `conservation_bookkeeping`
- interface role candidate: `physical_pair_field`

No source change or ray rerun is authorized by this lab.
