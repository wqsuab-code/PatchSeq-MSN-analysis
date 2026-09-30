param(
    [string]$BundleRoot = (Split-Path -Parent $PSScriptRoot)
)

$ErrorActionPreference = 'Stop'
$manifestDir = Join-Path $BundleRoot 'manifest'
New-Item -ItemType Directory -Path $manifestDir -Force | Out-Null
$manifestPath = Join-Path $manifestDir 'SHA256SUMS.txt'

$files = Get-ChildItem -LiteralPath $BundleRoot -Recurse -File |
    Where-Object { $_.FullName -ne $manifestPath } |
    Sort-Object FullName

$lines = foreach ($file in $files) {
    $relative = $file.FullName.Substring($BundleRoot.Length + 1).Replace('\', '/')
    $hash = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    "$hash  $relative"
}

[System.IO.File]::WriteAllLines($manifestPath, $lines, [System.Text.UTF8Encoding]::new($false))
Write-Output "Wrote $($lines.Count) checksums to $manifestPath"
