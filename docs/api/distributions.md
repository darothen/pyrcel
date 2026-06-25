# Distributions

Size distribution classes and pre-built aerosol climatology collections.

## Distribution classes

::: pyrcel.distributions.Lognorm

---

::: pyrcel.distributions.MultiModeLognorm

---

!!! note "Gamma distribution"
    `pyrcel.distributions.Gamma` is a placeholder class and is not yet implemented.

## Pre-built collections

Four published aerosol size distribution climatologies are available as module-level
dicts of `Lognorm` or `MultiModeLognorm` objects:

| Name | Source | Keys |
|---|---|---|
| `pyrcel.distributions.FN2005_single_modes` | Fountoukis & Nenes (2005) | `SM1`–`SM5` |
| `pyrcel.distributions.NS2003_single_modes` | Nenes & Seinfeld (2003) | `SM1`–`SM5` |
| `pyrcel.distributions.whitby_distributions` | Whitby (1978) | `marine`, `continental`, `background`, `urban` |
| `pyrcel.distributions.jaenicke_distributions` | Jaenicke (1993) | `Polar`, `Urban`, `Background`, `Maritime`, `Remote Continental`, `Rural` |

`whitby_distributions` values are lists of `Lognorm` objects (nucleation, accumulation,
coarse modes). `jaenicke_distributions` values are `MultiModeLognorm` objects.

```python
from pyrcel.distributions import whitby_distributions, jaenicke_distributions

marine_modes = whitby_distributions["marine"]  # [Lognorm, Lognorm, Lognorm]
urban = jaenicke_distributions["Urban"]        # MultiModeLognorm
```
