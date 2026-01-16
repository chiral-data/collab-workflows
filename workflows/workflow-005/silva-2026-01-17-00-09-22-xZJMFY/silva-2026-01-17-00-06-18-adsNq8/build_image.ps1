
$ImageName = "qsar-workflow"

Write-Host "Building Docker Image: $ImageName..." -ForegroundColor Cyan

docker build -t $ImageName .

if ($?) {
    Write-Host "Build Successful!" -ForegroundColor Green
    Write-Host "You can run the workflow with:"
    Write-Host "docker run -it --rm $ImageName"
}
else {
    Write-Host "Build Failed!" -ForegroundColor Red
    exit 1
}
