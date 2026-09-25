# workflow-900: E2E Test - Sleep (Job Cancellation)

Test-only workflow used by `chiral-test-e2e` (issue #74) to cancel a running
job through the UI and verify `Canceled` in the backend and Job History.

- One CPU step (`01-sleep`) that sleeps for 10 minutes in `debian:bookworm-slim`.
- No inputs, no GPU, tiny image.
- Lives on the `e2e-epic05` branch only. Do not merge to `main`: `main` feeds
  the Blueprints list in production.
