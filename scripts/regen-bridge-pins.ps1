<#
.SYNOPSIS
    Regenerates the pinned bridge version and digests in src/install.rs from a pyMzLib SHA256SUMS
    manifest.

.DESCRIPTION
    pyMzLib publishes `mzlib-bridge-<rid>.tar.gz` for four platforms on every release, alongside a
    `SHA256SUMS` manifest listing their digests. install.rs pins that version and those digests so a
    download can be verified against something that did not arrive down the same connection.

    Transcribing four digests from a browser is the chore the manifest exists to end, so this reads
    them from the same authority that produced the files. It rewrites only the region between the
    BEGIN/END generated bridge pins markers, and refuses to do anything at all if it cannot find
    both — a generator that silently appends when its anchor has moved is worse than one that stops.

    Written in PowerShell to match scripts/stage-bridge.ps1 and because pwsh is present on the CI
    runner and on the maintainer's machine alike.

    Run by .github/workflows/bridge-watch.yml, and runnable by hand.

.PARAMETER Manifest
    Path to a downloaded SHA256SUMS from a pyMzLib release.

.PARAMETER Target
    The file carrying the marked pin block. Defaults to src/install.rs beside this script.

.PARAMETER Version
    The pyMzLib release version, without the leading `v` — for example `0.1.0.dev4`.

.EXAMPLE
    .\scripts\regen-bridge-pins.ps1 -Manifest SHA256SUMS -Version 0.1.0.dev4
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$Manifest,

    [Parameter(Mandatory = $true)]
    [string]$Version,

    [string]$Target
)

$ErrorActionPreference = 'Stop'

if (-not $Target) {
    $Target = Join-Path (Split-Path -Parent $PSScriptRoot) 'src/install.rs'
}

# Read and write through .NET rather than Get-Content/Set-Content, because both of those get this
# file wrong on Windows PowerShell 5.1 and the damage is silent.
#
# Get-Content defaults to the system ANSI codepage, so every em dash in install.rs came back as
# mojibake and was written out double-encoded. Set-Content -Encoding utf8 writes a BOM, which then
# sits in front of the first `//!` of a Rust source file. The first run of this script did both at
# once: a one-line version bump arrived as a 25-line diff, every one of them a comment nobody had
# touched. Reading the script would not have shown it; running it did.
#
# UTF-8 without a BOM, LF endings — what the repository stores and what the Linux runner expects.
$utf8NoBom = [System.Text.UTF8Encoding]::new($false)

# The runtime identifiers this crate pins a bridge for, and the order they appear in. Set by
# pyMzLib's wheels.yml build matrix; if that matrix stops publishing one of these, the "no digest
# for" error below is what says so.
#
# There is no Python packaging vocabulary here: the payload is the raw tarball published for exactly
# this purpose, so the runtime identifier is the whole story.
$rids = @('win-x64', 'osx-arm64', 'osx-x64', 'linux-x64')

# --- read the manifest ------------------------------------------------------------------------

# `sha256sum` writes "<digest>  <name>", and "<digest> *<name>" when it was run in binary mode.
# Accept both rather than assuming which side produced the file.
$digests = @{}
foreach ($line in [System.IO.File]::ReadAllLines($Manifest, $utf8NoBom)) {
    if ($line.Trim() -match '^([0-9a-fA-F]{64})\s+\*?(.+)$') {
        $digests[$Matches[2].Trim()] = $Matches[1].ToLowerInvariant()
    }
}
if ($digests.Count -eq 0) {
    throw "No 'sha256  filename' lines found in $Manifest."
}

# --- build the replacement block --------------------------------------------------------------

# One string per LINE. mzLibR's generator built this by pasting multi-line chunks, which wrote a
# byte-identical file while its no-op check compared a 4-element vector against the 16 lines the
# file actually had — so it reported a change every run and would have opened an empty pull request
# weekly forever. The generated output was correct; only the change DETECTION was wrong, which is
# why reading it would not have caught it. Here the comparison is over the finished text (below), so
# the same mistake cannot be made.
$block = [System.Collections.Generic.List[string]]::new()
$block.Add('// BEGIN generated bridge pins')
$block.Add('/// The pyMzLib release whose published bridge this crate installs by default.')
$block.Add("pub const MZLIB_BRIDGE_VERSION: &str = ""$Version"";")
$block.Add('')
$block.Add('/// The platforms pyMzLib publishes a bridge for, and the SHA-256 of each tarball.')
$block.Add('pub const BRIDGE_ASSETS: &[(&str, &str)] = &[')
foreach ($rid in $rids) {
    $asset = "mzlib-bridge-$rid.tar.gz"
    if (-not $digests.ContainsKey($asset)) {
        throw ("No digest for $asset in $Manifest.`n" +
               "The manifest lists: " + (($digests.Keys | Sort-Object) -join ', '))
    }
    $block.Add('    (')
    $block.Add("        ""$rid"",")
    $block.Add("        ""$($digests[$asset])"",")
    $block.Add('    ),')
}
$block.Add('];')
$block.Add('// END generated bridge pins')

# --- splice it in -----------------------------------------------------------------------------

$source = @([System.IO.File]::ReadAllLines($Target, $utf8NoBom))
$begin = @(0..($source.Count - 1) | Where-Object { $source[$_] -eq '// BEGIN generated bridge pins' })
$end = @(0..($source.Count - 1) | Where-Object { $source[$_] -eq '// END generated bridge pins' })

if ($begin.Count -ne 1 -or $end.Count -ne 1 -or $end[0] -le $begin[0]) {
    throw ("Expected exactly one BEGIN and one END 'generated bridge pins' marker in $Target " +
           "(found $($begin.Count) and $($end.Count)). Nothing was written.")
}

$updated = [System.Collections.Generic.List[string]]::new()
if ($begin[0] -gt 0) { $updated.AddRange([string[]]$source[0..($begin[0] - 1)]) }
$updated.AddRange($block)
if ($end[0] -lt $source.Count - 1) { $updated.AddRange([string[]]$source[($end[0] + 1)..($source.Count - 1)]) }

# Compared as finished text, so "did anything change" cannot disagree with "what will be written".
if (($updated -join "`n") -eq ($source -join "`n")) {
    Write-Host "Already pinned to $Version - nothing to do."
    exit 0
}

[System.IO.File]::WriteAllText($Target, ($updated -join "`n") + "`n", $utf8NoBom)
Write-Host "Repinned $Target to pyMzLib $Version"
