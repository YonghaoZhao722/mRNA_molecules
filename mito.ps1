# 设置参数
$mitographExe = "F:\atp\MitoGraph.exe"
$xy = 0.0645
$z = 0.2
$rootPath = "F:\atp\Y333 ATP6 ATP3\extracted_cells"

# 获取所有子文件夹
$folders = Get-ChildItem -Path $rootPath -Directory

foreach ($folder in $folders) {
    $folderPath = $folder.FullName
    Write-Host "Processing $folderPath"
    & $mitographExe -xy $xy -z $z -path $folderPath
}