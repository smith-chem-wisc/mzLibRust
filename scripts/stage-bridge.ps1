<#
.SYNOPSIS
    Stages an mzLib bridge executable so mzLibRust can reach mzLib.

.DESCRIPTION
    The crate compiles and its offline suite passes without a bridge; only calls that actually
    reach mzLib need one. This puts a bridge at _dotnet/<rid>/ where build.rs looks.

    Two sources, in order of what you probably have:

      -FromPyMzLib   copy the payload pyMzLib already staged in its wheel tree. Instant.
      -Build         run pyMzLib's publish-bridge.ps1 to produce a fresh one. Needs .NET.

    Neither is required if you would rather just set MZLIB_BRIDGE to a bridge you already have;
    that always wins over anything staged here.

.PARAMETER PyMzLibRoot
    Path to a pyMzLib checkout (the directory holding pkg/).

.PARAMETER Runtime
    .NET runtime identifier. Defaults to this machine's.

.PARAMETER Build
    Build a fresh bridge rather than copying a staged one.

.EXAMPLE
    .\scripts\stage-bridge.ps1 -PyMzLibRoot E:\CodeReview\pyMzLib\code\pyMzLib

.EXAMPLE
    .\scripts\stage-bridge.ps1 -PyMzLibRoot ..\pyMzLib -Runtime linux-x64 -Build
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$PyMzLibRoot,

    [string]$Runtime,

    [switch]$Build
)

$ErrorActionPreference = 'Stop'

if (-not $Runtime) {
    $arch = if ([System.Runtime.InteropServices.RuntimeInformation]::OSArchitecture -eq 'Arm64') { 'arm64' } else { 'x64' }
    $Runtime = switch ($true) {
        $IsLinux { "linux-$arch"; break }
        $IsMacOS { "osx-$arch"; break }
        default  { "win-$arch" }
    }
}

$exeName = if ($Runtime.StartsWith('win')) { 'mzlib-bridge.exe' } else { 'mzlib-bridge' }
$crateRoot = Split-Path -Parent $PSScriptRoot
$destination = Join-Path $crateRoot "_dotnet/$Runtime"

if (-not (Test-Path $PyMzLibRoot)) {
    throw "pyMzLib checkout not found at '$PyMzLibRoot'."
}

if ($Build) {
    $publish = Join-Path $PyMzLibRoot 'pkg/build/publish-bridge.ps1'
    if (-not (Test-Path $publish)) {
        throw "publish-bridge.ps1 not found at '$publish'. Is -PyMzLibRoot pointing at the pyMzLib repo root?"
    }
    Write-Host "Building the bridge for $Runtime via pyMzLib's publish-bridge.ps1..."
    & $publish -Runtime $Runtime
    if ($LASTEXITCODE -ne 0) { throw "publish-bridge.ps1 failed with exit code $LASTEXITCODE" }
}

$source = Join-Path $PyMzLibRoot "pkg/python/src/pymzlib/_dotnet/$Runtime/$exeName"
if (-not (Test-Path $source)) {
    throw @"
No staged bridge at '$source'.
Either run this script with -Build to produce one, or set MZLIB_BRIDGE to a bridge you already have.
"@
}

New-Item -ItemType Directory -Force -Path $destination | Out-Null
Copy-Item -Path $source -Destination (Join-Path $destination $exeName) -Force

$staged = Join-Path $destination $exeName
$sizeMb = [math]::Round((Get-Item $staged).Length / 1MB, 1)
Write-Host "Staged $exeName ($sizeMb MB) at $staged"

# The whole point is that it runs. A payload that cannot report its own version here will
# certainly fail when the crate calls it.
$isNative = ($Runtime.StartsWith('win') -and $IsWindows -ne $false) -or
            ($Runtime.StartsWith('linux') -and $IsLinux) -or
            ($Runtime.StartsWith('osx') -and $IsMacOS)
if ($isNative) {
    $probe = & $staged version | ConvertFrom-Json
    if (-not $probe.ok) { throw "The staged bridge failed its version probe." }
    Write-Host "Probe OK: bridge $($probe.data.bridge), protocol $($probe.data.protocol), runtime $($probe.data.runtime)"
    Write-Host ""
    Write-Host "Now run:  cargo test --features live"
} else {
    Write-Host "Cross-staged for $Runtime; skipping the version probe (cannot run it here)."
}
