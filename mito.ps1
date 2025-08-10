# Set parameters
$mitographExe = ".\MitoGraph.exe"
$xy = 0.0645
$z = 0.2
$threshold = 0.15       
$scalesMin = 1.1        
$scalesMax = 1.5        
$scalesCount = 8         
$zAdaptive = $true
$zBlockSize = 8
$rootPath = "Y333 ATP6 ATP2\extracted_cells_dw_30"

Write-Host "Script started"
Write-Host "MitoGraph executable: $mitographExe"
Write-Host "Parameters: xy=$xy, z=$z, threshold=$threshold"
Write-Host "Scales: $scalesMin to $scalesMax with $scalesCount levels"
Write-Host "Z-adaptive: $zAdaptive"
Write-Host "Z-block size: $zBlockSize"
Write-Host "Root path: $rootPath"

# Check if MitoGraph executable exists
if (-not (Test-Path $mitographExe)) {
    Write-Error "MitoGraph executable not found at $mitographExe"
    exit 1
}

# Check if root path exists
if (-not (Test-Path $rootPath)) {
    Write-Error "Root path not found: $rootPath"
    exit 1
}

# Get all subdirectories
Write-Host "Processing $rootPath"
$folders = Get-ChildItem -Path $rootPath -Directory
$folderCount = $folders.Count
Write-Host "Found $folderCount folders to process"

# Process each folder
foreach ($folder in $folders) {
    $folderPath = $folder.FullName
    $folderName = $folder.Name
    Write-Host "Processing $folderPath"
    
    # Build command with comprehensive parameters to improve tubular structure detection
    $cmdArgs = @(
        "-xy", $xy,
        "-z", $z,
        "-path", $folderPath,
        "-threshold", $threshold,
        "-scales", $scalesMin, $scalesMax, $scalesCount
    )
    
    if ($zAdaptive) {
        $cmdArgs += @("-z-adaptive")
        $cmdArgs += @("-z-block-size", $zBlockSize)
    }
    
    # Execute the command
    try {
        & $mitographExe @cmdArgs
        Write-Host "Completed processing $folderName"
    } catch {
        Write-Error "Error processing $folderPath : $_"
    }
}

Write-Host "Script completed"