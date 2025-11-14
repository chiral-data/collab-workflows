# container-images-silva

## How to use

Application container images are hosted on our private registry at `chiral.sakuracr.jp`.
To pull a specific application image, use the following format:

```bash
docker pull chiral.sakuracr.jp/<app_name>:<date_tag>
```

Where:

- `<app_name>`: The name of the application (e.g., `boltz`, `gromacs`).
- `<date_tag>`: The version or date-based tag for the image. This typically represents a specific build on that day.
  The latest `date_tag` is `2025_09_05`.
  **Example:** To pull the `gromacs` application image using the current tag:

```bash
docker pull chiral.sakuracr.jp:/gromacs:2025_09_05
```

## How to build a new image

1. create a new directory as `./a/app_date`. For test builds, append a version suffix, such as `_v1` (e.g., `app_date_v1`).
2. create the `Dockerfile`
3. Execute the `build.sh` script from the project's root directory using the command: `bash build.sh ./a/app_date_v1`.
