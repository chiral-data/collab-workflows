# workflow-900: Sleep

A minimal workflow with one CPU step (`01-sleep`) that counts up one line per
second for a configurable duration, so a job stays in Processing long enough
to exercise long-running job behaviour such as cancellation.

- Parameter `sleep_seconds` (default 120) sets the duration.
- Runs in `debian:bookworm-slim`; no inputs, no GPU.
- Not for `main`: `main` feeds the production Blueprints list.
