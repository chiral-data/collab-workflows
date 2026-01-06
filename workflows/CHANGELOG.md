# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2026-01-06]

### Changed

- workflow-007/03_visualize/generate_output.py: Improved HTML output with full-page responsive layout
  - Wrapped 3Dmol viewer in a complete HTML document structure
  - Full-page viewer container (max-width: 1400px)
  - Added zoom-out script to show protein at 70% for better overview
  - Clean header with title, card-style viewer with shadow
  - Added UI controls for protein representation styles (Cartoon, Stick, Sphere, Line, Surface)
  - Fixed: Surface now properly removed when switching to other representations

## [2025-12-08]

### Added

- workflow-007/example.py: Merged 3 visualization scripts with configurable options (pocket_style, render_method, representation)
- workflow-007/example.py: Tested with pocketeer docker image - successfully generated pocket_visualization.html and rotating GIF
- workflow-007/example.py: Added output_format option ("html", "gif", or "both")
- workflow-007: Split example.py into 3 modular scripts (tested):
  - 01_download/download_pdb.py: Download PDB file
  - 02_pocket/calculate_pockets.py: Calculate pockets and save to JSON
  - 03_visualize/generate_output.py: Generate HTML/GIF visualization
- workflow-007/03_visualize: HTML visualization now uses pocket_style and representation options
  - receptor_cartoon/receptor_surface based on representation setting
  - sphere_scale adjusted based on pocket_style (2.0 for single_sphere, 1.0 for filled_surfaces)
- workflow-007: Added Chiral workflow configuration:
- workflow-003: Migrated from old silva config format to new format:
  - Converted workflow.json to workflow.toml
  - Merged node.json into job.toml for all 6 jobs
  - Updated container.docker_image to container.image
  - Added name, description to job.toml files
  - Removed deprecated node.json and workflow.json files
  - Tested successfully with silva
