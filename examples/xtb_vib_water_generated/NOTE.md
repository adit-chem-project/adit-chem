# About this example (added 2026-09-14)

**The structure in this example has not been geometry-optimized.** It is a water molecule with O–H stretched to 1.27 Å (measured: about 0.96 Å).
So the frequencies come out as 819 / 1100 / 1657 cm⁻¹ and **do not match measured water (1595 / 3657 / 3756 cm⁻¹).**
No vibration corresponding to the O–H stretch appears.

Why it is here: **this example checks whether ADIT notices "a vibrational analysis on a structure with residual forces."**
Running the analysis gives the following note.

```
note: on the starting structure of the vibrational analysis, the largest force on an atom is at least 6.52 eV/Å …
   the convergence threshold for geometry optimization in these settings is 0.00514 eV/Å.
   vibrational analysis is normally done on an optimized structure (forces close to zero)
```

**Order for your own molecule**

1. First run a geometry optimization (`task.type = "geometry_optimization"`)
2. When it converges, run the vibrational analysis (`task.type = "vibrations"`) on that final structure
3. To compare with measurements, use `--freq-scale` (a correction factor such as 0.96; you choose the value) and
   `--spectrum-measured` (overlays a measured table)

This note was added because, in a trial on 2026-09-14, a user with an experimental background took this example as a model and got stuck for 10 minutes.
