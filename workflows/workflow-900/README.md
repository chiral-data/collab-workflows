# workflow-900: E2E Test - Sleep (Job Cancellation)

Test-only workflow used by `chiral-test-e2e` (issue #74) to cancel a running
job through the UI and verify `Canceled` in the backend and Job History.

- One CPU step (`01-sleep`) that sleeps in `debian:bookworm-slim`; duration is
  configurable via `sleep_seconds` (default 120s).
- No inputs, no GPU, tiny image.
- Lives on the `e2e-epic05` branch only. Do not merge to `main`: `main` feeds
  the Blueprints list in production.
