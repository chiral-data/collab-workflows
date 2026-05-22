# PowerShell Script to Start QSAR Web Application
# Usage: .\start_server.ps1 (from inside 04-prediction folder)

# Configuration
$Port = "5000"
$Image = "ghcr.io/chiral-data/qsar:20260107_v1"

# Determine paths dynamically
# This allows running the script from anywhere (root or inside directory)
$ScriptPath = $MyInvocation.MyCommand.Definition
$ScriptDir = Split-Path -Parent $ScriptPath
$NodeDirName = (Get-Item $ScriptDir).Name
$RootDir = (Get-Item $ScriptDir).Parent.FullName

# Workspace paths for Docker
# We mount the Project Root to /workspace
# And set working directory to /workspace/04-prediction
$ContainerWorkDir = "/workspace/$NodeDirName"

Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "🧪 Starting QSAR Web Application Dashboard" -ForegroundColor Cyan
Write-Host "========================================="
Write-Host ""
Write-Host "DEBUG: Script Directory: $ScriptDir"
Write-Host "DEBUG: Root Directory:   $RootDir"
Write-Host ""

# Check and Copy missing artifacts (Model, Scaler, AD Stats)
$OutputsDir = "$ScriptDir\outputs"
if (-not (Test-Path $OutputsDir)) { New-Item -ItemType Directory -Path $OutputsDir -Force | Out-Null }

# Helper to find and copy file
function Ensure-Artifact ($FileName, $UpstreamRelPath) {
    if (-not (Test-Path "$OutputsDir\$FileName")) {
        $SourcePath = Join-Path $RootDir $UpstreamRelPath
        if (Test-Path $SourcePath) {
            Write-Host "Copying $FileName from upstream..." -ForegroundColor Gray
            Copy-Item $SourcePath -Destination "$OutputsDir\$FileName"
        }
        else {
            Write-Host "Warning: $FileName not found at $SourcePath" -ForegroundColor Yellow
        }
    }
}

# 1. Model (Node 3)
Ensure-Artifact "model.h5" "03-model-training\outputs\model.h5"

# 2. Scaler & AD Stats (Node 2)
Ensure-Artifact "scaler.pkl" "02-feature-engineering\outputs\scaler.pkl"
Ensure-Artifact "ad_stats.json" "02-feature-engineering\outputs\ad_stats.json"

# Final Check
if (-not (Test-Path "$OutputsDir\model.h5")) {
    Write-Host "Error: Model file 'model.h5' missing in $OutputsDir" -ForegroundColor Red
    Write-Host "Please run 'silva .' first to generate the model."
    exit 1
}

Write-Host "• Docker Image: $Image"
Write-Host "• Port: $Port"
Write-Host "• URL: http://localhost:$Port"
Write-Host ""
Write-Host "Starting server... (Press Ctrl+C to stop)" -ForegroundColor Yellow
Write-Host ""

# Run Docker container
# -it: Interactive terminal
# --rm: Remove container after exit
# -p: Map port 5000
# -v: Mount RootDir to /workspace
# -w: Set working directory to the node folder inside container
docker run -it --rm -p ${Port}:${Port} -v "${RootDir}:/workspace" -w $ContainerWorkDir $Image python app.py
