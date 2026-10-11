# workflow-900: Sleep

A minimal workflow of CPU steps for exercising job behaviour rather than
doing science. Both run in `debian:bookworm-slim`, with no inputs and no GPU.

`01-sleep` counts up one line per second for a configurable duration, so a
job stays in Processing long enough to exercise long-running job behaviour
such as cancellation.

- Parameter `sleep_seconds` (default 120) sets the duration.

`02-count` prints many lines quickly, for exercising log viewers under volume.

- `lines` (default 1000): prints `1`..`N` on stdout, one per line.
- `lines_per_second` (default 10): prints in bursts of this many lines, each
  followed by `sleep 1` and an `@ <n> <epoch ms>` line saying when line `n`
  was printed.
- `stderr_every` (default 0, off): prints `stderr <k>` on stderr for every
  multiple `k` of this value.

Not for `main`: `main` feeds the production Blueprints list.
