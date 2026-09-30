param(
    [string]$BundleRoot = (Split-Path -Parent $PSScriptRoot)
)

$ErrorActionPreference = 'Stop'

function Assert-Equal {
    param([string]$Name, $Observed, $Expected)
    if ($Observed -ne $Expected) {
        throw "${Name}: expected $Expected, observed $Observed"
    }
    Write-Output "PASS $Name = $Observed"
}

$globalDir = Join-Path $BundleRoot 'data/global'
$msnDir = Join-Path $BundleRoot 'data/msn'
$qcDir = Join-Path $BundleRoot 'data/qc'

$scan = Import-Csv -LiteralPath (Join-Path $globalDir 'stage02B_RPCA_scan_summary.csv')
Assert-Equal 'global scan settings' $scan.Count 11
Assert-Equal 'global scan nPC grid' (($scan.nPC -join ',')) '10,12,14,16,18,20,22,24,26,28,30'

$cellSummary = Import-Csv -LiteralPath (Join-Path $globalDir '03A_final_RPCA_CellType_summary.csv')
$rootSummary = Import-Csv -LiteralPath (Join-Path $globalDir '03B_final_RPCA_Root_summary.csv')
Assert-Equal 'stable D1' (($cellSummary | Where-Object Final_CellType_stability -eq 'Stable_D1').n_cells) 272
Assert-Equal 'stable D2' (($cellSummary | Where-Object Final_CellType_stability -eq 'Stable_D2').n_cells) 279
Assert-Equal 'stable IN' (($cellSummary | Where-Object Final_CellType_stability -eq 'Stable_IN').n_cells) 34
Assert-Equal 'stable MSN root' (($rootSummary | Where-Object Final_Root_stability_from_CellType -eq 'Stable_MSN').n_cells) 588
Assert-Equal 'stable IN root' (($rootSummary | Where-Object Final_Root_stability_from_CellType -eq 'Stable_IN').n_cells) 34

$pearsonQc = Import-Csv -LiteralPath (Join-Path $qcDir 'Top15_risk_x040_x050_QC_table.csv')
$npc5Qc = Import-Csv -LiteralPath (Join-Path $qcDir 'Top15CleanMarker_MaxPC_NPC5_QC_table.csv')
Assert-Equal 'Pearson QC rows' $pearsonQc.Count 641
Assert-Equal 'Pearson retained' ($pearsonQc | Where-Object risk_class -eq 'retained').Count 630
Assert-Equal 'five-PC retained' ($npc5Qc | Where-Object risk_class -eq 'retained').Count 623

$bootstrap = Import-Csv -LiteralPath (Join-Path $msnDir 'MSN_top750_PC16_gene_bootstrap_500_per_cell.csv')
Assert-Equal 'MSN bootstrap rows' $bootstrap.Count 588
Assert-Equal 'RPCA-Pearson exact agreement' ($bootstrap | Where-Object RPCA_agrees_bootstrap_top1 -eq 'TRUE').Count 422

$tierCounts = @{
    High = 0
    Moderate = 0
    Low = 0
}
foreach ($row in $bootstrap) {
    $support = [double]$row.bootstrap_top1_fraction
    $margin = [double]$row.bootstrap_margin
    if ($support -ge 0.80 -and $margin -ge 0.50) {
        $tierCounts.High++
    } elseif ($support -ge 0.50 -and $margin -ge 0.20) {
        $tierCounts.Moderate++
    } else {
        $tierCounts.Low++
    }
}
Assert-Equal 'High tier' $tierCounts.High 339
Assert-Equal 'Moderate tier' $tierCounts.Moderate 175
Assert-Equal 'Low tier' $tierCounts.Low 74

$requiredLayouts = @(
    'Figure_1_target.png',
    'Extended_Data_Figure_1_target.png',
    'Extended_Data_Figure_2_target.png'
)
foreach ($name in $requiredLayouts) {
    $path = Join-Path (Join-Path $BundleRoot 'reference_layouts') $name
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Missing reference layout: $name"
    }
    Write-Output "PASS layout present: $name"
}

Write-Output 'ALL FROZEN CHECKS PASSED'
