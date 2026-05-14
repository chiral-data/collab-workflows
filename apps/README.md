# container-images-silva

## How to use

Application container images are hosted on our public registry at `ghcr.io/chiral-data`.
To pull a specific application image, use the following format:

```bash
docker pull ghcr.io/chiral-data/<app_name>:<date_tag>
```

Where:

- `<app_name>`: The name of the application (e.g., `boltz`, `gromacs`).
- `<date_tag>`: The version or date-based tag for the image. This typically represents a specific build on that day.
  The latest `date_tag` is `2025_09_05`.
  **Example:** To pull the `gromacs` application image using the current tag:

```bash
docker pull ghcr.io/chiral-data/gromacs:2025_09_05
```

## How to build a new image

1. create a new directory as `./a/app_date`. For test builds, append a version suffix, such as `_v1` (e.g., `app_date_v1`).
2. create the `Dockerfile`
3. Execute the `build.sh` script from the project's root directory using the command: `bash build.sh ./a/app_date_v1`.

---

# Workflow Job Images

Container images for workflow job folders are hosted at `ghcr.io/chiral-data` alongside app images, but follow a different naming and build convention.

## Image naming

Image names and tags are defined in the job's `.chiral/job.toml`:

```toml
[container]
image = "sequence-ingest:latest"
```

The CI prepends `ghcr.io/chiral-data/` automatically. The `image` field in `job.toml` should always be a bare name without registry prefix.

## Eligibility

A job folder is built only if it contains **both**:
- `Dockerfile` — the build definition
- `.chiral/job.toml` — with a `[container].image` field

## How to build a new workflow job image

### Automatic (on push to `main`)

Any push to `main` that touches `workflows/*/*/**` triggers `build-changed-workflow-images.yml`, which detects changed eligible job folders and builds them automatically.

### Manual

Trigger `build-workflow-image.yml` via GitHub Actions → Run workflow:

| Field | Description | Example |
|---|---|---|
| `job_dir` | Path relative to `workflows/` | `workflow-019/01-sequence-ingest` |
| `force_rebuild` | Rebuild even if image exists | `false` |

## Build order

If a job's `Dockerfile` uses `FROM <bare-name>` where `<bare-name>` matches another job's image in this repo, it is treated as a **dependent** job and built only after all independent jobs complete. This ensures e.g. `esmfold:latest` is pushed before any job that builds `FROM esmfold:latest`.
