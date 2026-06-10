# Headline table: baseline vs. method arms

K9c live-merged population: N=400 designs, 200 unique (kind,id) graphs x 2 sigma rules (full target).

| arm | N | deployable ($\Theta_{max}>0$) | certified ($\varepsilon_{max}>0$) | $\varepsilon_{max}\geq0.1$ rad | median $\varepsilon_{max}$ | max $\varepsilon_{max}$ |
|---|---|---|---|---|---|---|
| Eq.(6)+repairs (K6) | 400 | 0 | 0 | 0 | -- | -- |
| K9 proximity | 400 | 36 | 36 | 33 | 0.253 | 1.997 |
| K9b | 400 | 27 | 27 | 14 | 0.115 | 1.828 |
| K9c range-max | 400 | 307 | 306 | 207 | 0.254 | 3.142 |
| best-of-3 (per-design max over K9/K9b/K9c) | 400 | 307 | 306 | 214 | 0.278 | 3.142 |

