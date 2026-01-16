# Build unified Docker image for Vina workflow
# This script builds a single image containing all tools for all 8 nodes

Write-Host "=== Building Unified Vina Workflow Docker Image ===" -ForegroundColor Green

$IMAGE_NAME = "vina-workflow:latest"

Write-Host "`nBuilding $IMAGE_NAME..." -ForegroundColor Cyan
docker build -t $IMAGE_NAME .

if ($LASTEXITCODE -ne 0) {
    Write-Host "`nBuild failed!" -ForegroundColor Red
    exit 1
}

Write-Host "`n=== Updating job.toml files ===" -ForegroundColor Green

# Update all nodes to use the unified image
$nodes = @(
    "01_ProteinInput",
    "02_LigandInput", 
    "03_ProteinPreparation",
    "04_LigandPreparation",
    "05_PocketPrediction",
    "06_PocketSelection",
    "07_Docking",
    "08_Reporting"
)

foreach ($node in $nodes) {
    $jobToml = "$node/.chiral/job.toml"
    if (Test-Path $jobToml) {
        (Get-Content $jobToml) -replace 'image = ".*"', "image = `"$IMAGE_NAME`"" | Set-Content $jobToml
        Write-Host "  ✓ Updated $node" -ForegroundColor Gray
    }
}

Write-Host "`n✅ Build complete!" -ForegroundColor Green
Write-Host "`nImage built: $IMAGE_NAME" -ForegroundColor Yellow
Write-Host "All nodes configured to use unified image" -ForegroundColor Yellow
Write-Host "`nYou can now run: silva ." -ForegroundColor Cyan
