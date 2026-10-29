# Controlled pause for larger-sample benchmarks

The scale-study controller was paused while the final `A-invariance/hratt-reml__nonumba` worker was running, to allow the requested larger-sample studies to run without overlapping fits. The worker continued to completion before those studies began. The controller will then resume the remaining scale-study arms.

The affected worker's internal fit time and numerical results do not include the controller pause, but its parent-measured wall time includes the interval before the controller resumes. Do not use that wall time as a resource measurement. A fresh, isolated repeat of this one worker will be recorded separately and compared against its results.

The repeat is in `pause-repeat/hratt-reml__nonumba/`. It took 133.788 s
whole-process wall time and 0.371 GiB peak RSS, and every saved result
array was exactly equal to the paused worker's output. Its own
`measurement.json`, diagnostics, job configuration and log are archived.
The original worker's parent-measured wall time remains in the study's
raw CSV for provenance and must not be interpreted as a runtime.
