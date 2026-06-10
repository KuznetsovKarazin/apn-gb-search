# run_batch.ps1
# Robust Magma batch runner for apn-gb-search.
#
# Run this script FROM the batch\ directory:
#   cd <repo>\batch
#   powershell -ExecutionPolicy Bypass -File ..\magma_templates\run_batch.ps1
#
# Optional:
#   powershell -ExecutionPolicy Bypass -File ..\magma_templates\run_batch.ps1 -MaxParallel 6
#   powershell -ExecutionPolicy Bypass -File ..\magma_templates\run_batch.ps1 -Pattern "apn_ds_test_*.m"
#   powershell -ExecutionPolicy Bypass -File ..\magma_templates\run_batch.ps1 -NoClean
#
# It accepts all generated Magma files by default:
#   apn_*.m
#
# Directory layout created inside batch\:
#   results\apn_results_*.txt
#   logs\log_<job>.txt
#   errors\err_<job>.txt

param(
    [string]$MagmaPath = "E:\Magma\magma.exe",
    [int]$MaxParallel = 1,
    [string]$Pattern = "apn_*.m",
    [switch]$NoClean
)

$ErrorActionPreference = "Stop"

function Ensure-Dir($Path) {
    if (-not (Test-Path -LiteralPath $Path)) {
        New-Item -ItemType Directory -Path $Path | Out-Null
    }
}

function Count-Apn-Lines() {
    $count = 0
    Get-ChildItem -Path "results" -Filter "apn_results_*.txt" -File -ErrorAction SilentlyContinue | ForEach-Object {
        $count += (Select-String -Path $_.FullName -Pattern "^APN\|" -ErrorAction SilentlyContinue).Count
    }
    return $count
}

# ----------------------------------------------------------------------
# Basic checks
# ----------------------------------------------------------------------

if (-not (Test-Path -LiteralPath $MagmaPath)) {
    Write-Host "ERROR: Magma executable not found:"
    Write-Host "  $MagmaPath"
    Write-Host ""
    Write-Host "Edit -MagmaPath or run for example:"
    Write-Host '  powershell -ExecutionPolicy Bypass -File ..\magma_templates\run_batch.ps1 -MagmaPath "E:\Magma\magma.exe"'
    exit 1
}

Ensure-Dir "results"
Ensure-Dir "logs"
Ensure-Dir "errors"

if (-not $NoClean) {
    Write-Host "Cleaning old logs, errors, and results..."
    Remove-Item "results\apn_results_*.txt" -ErrorAction SilentlyContinue
    Remove-Item "logs\log_*.txt"           -ErrorAction SilentlyContinue
    Remove-Item "errors\err_*.txt"         -ErrorAction SilentlyContinue
    Write-Host "Done."
    Write-Host ""
}
else {
    Write-Host "NoClean enabled: old logs/results are kept."
    Write-Host ""
}

# ----------------------------------------------------------------------
# Build job list
# ----------------------------------------------------------------------

$allJobs = @(Get-ChildItem -Path "." -Filter $Pattern -File -ErrorAction SilentlyContinue | Sort-Object Name)

if ($allJobs.Count -eq 0) {
    Write-Host "ERROR: no Magma .m files found in current directory."
    Write-Host "Current directory:"
    Write-Host ("  " + (Get-Location).Path)
    Write-Host "Pattern:"
    Write-Host ("  " + $Pattern)
    Write-Host ""
    Write-Host "Examples from repo root:"
    Write-Host "  python src\gen_batch.py --reuse"
    Write-Host "  python src\gen_batch_from_dataset.py --file data\new_apns.txt --count 20 --tag-prefix ds_test"
    Write-Host ""
    Write-Host "Then run:"
    Write-Host "  cd batch"
    Write-Host "  powershell -ExecutionPolicy Bypass -File ..\magma_templates\run_batch.ps1"
    exit 1
}

Write-Host ("Total jobs: " + $allJobs.Count + "  Max parallel: " + $MaxParallel + "  Pattern: " + $Pattern)
Write-Host ("Magma: " + $MagmaPath)
Write-Host ""

# ----------------------------------------------------------------------
# Queue runner
# ----------------------------------------------------------------------

$running = @()
$doneCount = 0
$jobIndex = 0
$failedCount = 0

Write-Host ("{0,-10} {1,-8} {2,-8} {3,-8} {4,-10}" -f "Time","Running","Done","Queued","APN_saved")
Write-Host ("-" * 56)

while ($jobIndex -lt $allJobs.Count -or $running.Count -gt 0) {

    while ($running.Count -lt $MaxParallel -and $jobIndex -lt $allJobs.Count) {
        $jobFile = $allJobs[$jobIndex]
        $tag = [System.IO.Path]::GetFileNameWithoutExtension($jobFile.Name)

        $log = Join-Path "logs"   ("log_" + $tag + ".txt")
        $err = Join-Path "errors" ("err_" + $tag + ".txt")

        $proc = Start-Process `
            -FilePath $MagmaPath `
            -ArgumentList @($jobFile.Name) `
            -RedirectStandardOutput $log `
            -RedirectStandardError  $err `
            -PassThru `
            -WindowStyle Minimized

        $running += ,([PSCustomObject]@{
            Process = $proc
            Tag     = $tag
            Log     = $log
            Err     = $err
            File    = $jobFile.Name
        })

        Write-Host ("  Started: " + $tag + "  PID=" + $proc.Id)
        $jobIndex++
        Start-Sleep -Milliseconds 250
    }

    Start-Sleep -Seconds 30

    $still = @()
    foreach ($item in $running) {
        if ($item.Process.HasExited) {
            $doneCount++

            if ($item.Process.ExitCode -notin @(0,1)) {
                $failedCount++
                Write-Host ("  WARNING: " + $item.Tag + " exited with code " + $item.Process.ExitCode)
                Write-Host ("           log: " + $item.Log)
                Write-Host ("           err: " + $item.Err)
            }
        }
        else {
            $still += ,$item
        }
    }
    $running = $still

    $ts = Get-Date -Format "HH:mm:ss"
    $apnCount = Count-Apn-Lines
    $queued = $allJobs.Count - $jobIndex

    Write-Host ("{0,-10} {1,-8} {2,-8} {3,-8} {4,-10}" -f $ts, $running.Count, $doneCount, $queued, $apnCount)
}

# ----------------------------------------------------------------------
# Final summary
# ----------------------------------------------------------------------

$totalApn = Count-Apn-Lines

Write-Host ""
Write-Host "All done."
Write-Host ("Total jobs             : " + $allJobs.Count)
Write-Host ("Failed jobs            : " + $failedCount)
Write-Host ("Total APN lines saved  : " + $totalApn)
Write-Host ("Results directory      : " + (Join-Path (Get-Location).Path "results"))
Write-Host ("Logs directory         : " + (Join-Path (Get-Location).Path "logs"))
Write-Host ("Errors directory       : " + (Join-Path (Get-Location).Path "errors"))
Write-Host ""

$repoRoot = (Get-Item (Get-Location).Path).Parent.FullName
Write-Host "Next step:"
Write-Host ("  cd " + $repoRoot)
Write-Host "  python src\collect_results.py"
