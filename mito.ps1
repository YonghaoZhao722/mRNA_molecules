# 设置参数
$mitographExe = ".\MitoGraph.exe"
$xy = 0.0645
$z = 0.2
$threshold = 0.3           # 降低阈值以保留更多连接
$scalesMin = 0.8           # 扩展尺度范围下限
$scalesMax = 2.0           # 扩展尺度范围上限
$scalesCount = 8           # 增加尺度数量
$adaptive = 5              # 启用自适应阈值
$rootPath = "Y333 ATP6 ATP2\extracted_cells"

Write-Host "Script started"
Write-Host "MitoGraph executable: $mitographExe"
Write-Host "Parameters: xy=$xy, z=$z, threshold=$threshold"
Write-Host "Scales: $scalesMin to $scalesMax with $scalesCount levels"
Write-Host "Adaptive threshold: $adaptive"
Write-Host "Root path: $rootPath"

# 检查MitoGraph.exe是否存在
if (-not (Test-Path $mitographExe)) {
    Write-Error "MitoGraph.exe not found at $mitographExe"
    exit 1
}

# 检查根路径是否存在
if (-not (Test-Path $rootPath)) {
    Write-Error "Root path not found: $rootPath"
    exit 1
}

# 获取所有子文件夹
$folders = Get-ChildItem -Path $rootPath -Directory
Write-Host "Processing $rootPath"
Write-Host "Found $($folders.Count) folders to process"

foreach ($folder in $folders) {
    $folderPath = $folder.FullName
    Write-Host "Processing $folderPath"
    try {
        # 使用改进的参数来减少骨架断点
        & $mitographExe -xy $xy -z $z -path $folderPath -threshold $threshold -scales $scalesMin $scalesMax $scalesCount -adaptive $adaptive
        Write-Host "Completed processing $($folder.Name)"
    } catch {
        Write-Error "Error processing $folderPath : $_"
    }
}

Write-Host "Script completed"