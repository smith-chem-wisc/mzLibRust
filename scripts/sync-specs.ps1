<#
.SYNOPSIS
    Copy the per-verb specs from a bridge checkout into docs/specs/, and record where they came from.

.DESCRIPTION
    The facts about each wire verb - its parameters and their units, its result fields and what a
    null in each one means, its error kinds, caveats, the mzLib code it wraps, citations, and its
    spelling in every binding - are written once, in the bridge repository
    (design/verbs/<module>.<verb>.yaml), so that pyMzLib, mzLibRust and mzLibR cannot drift apart on
    facts. That repository is private and this one is public, so CI here cannot read it. The copy in
    docs/specs/ is what the rendered reference facts (docs/reference/), the doc lint
    (tests/spec_docs.rs) and the doctest replay bridge (tools/replay-bridge) read.

    Never edit docs/specs/*.yaml by hand. A fact that is wrong is wrong in the bridge, and fixing it
    there is what makes every binding's lint flag the pages that repeated it. Fix it there, then run
    this script, then regenerate the fragments:

        $env:MZLIB_RENDER_SPEC_DOCS = 1; cargo test --test spec_docs

    The copy is a mirror: a spec that no longer exists in the source is deleted here. docs/specs/SOURCE
    records the bridge commit, and says so when the source had uncommitted changes, so a reviewer can
    tell a vendored spec that matches a bridge commit from one that matches someone's working tree.

    The pyMzLib counterpart is scripts/sync_specs.py; the two write the same SOURCE file.

.PARAMETER From
    The bridge repository's design/verbs directory.

.PARAMETER WorkingTree
    Vendor the working tree, uncommitted specs and edits included, instead of the committed specs
    at the bridge's HEAD. SOURCE then records uncommitted_changes: yes. The default is HEAD because
    the bridge is shared by several sessions at once, and a spec someone is still writing is not a
    fact yet.

.EXAMPLE
    .\scripts\sync-specs.ps1 -From E:\CodeReview\bridge\design\verbs
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$From,
    [switch]$WorkingTree
)

$ErrorActionPreference = 'Stop'

$source = (Resolve-Path -LiteralPath $From).Path
$dest = Join-Path (Split-Path -Parent $PSScriptRoot) 'docs/specs'

$utf8 = New-Object System.Text.UTF8Encoding($false)
# PowerShell decodes a native command's stdout with the console code page, which would mangle the
# specs' non-ASCII characters (section signs, dashes) on a default Windows console.
[Console]::OutputEncoding = $utf8

# git writes to stderr for a path it does not know, and under 'Stop' PowerShell 5.1 turns that into
# a terminating error; every git call below is checked by its exit code instead.
function Invoke-Git {
    $saved = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    try { & git -C $source @args 2>$null } finally { $ErrorActionPreference = $saved }
}

$listing = @(Invoke-Git ls-tree --name-only HEAD -- .)
$listed = $LASTEXITCODE -eq 0
$committed = @($listing | Where-Object { $_ -like '*.yaml' } | ForEach-Object { Split-Path -Leaf $_ })
$useHead = (-not $WorkingTree) -and $listed -and ($committed.Count -gt 0)
if ($useHead) {
    $specs = @($committed | Sort-Object | ForEach-Object { [pscustomobject]@{ Name = $_ } })
} else {
    $specs = @(Get-ChildItem -LiteralPath $source -Filter '*.yaml' -File | Sort-Object Name)
}
if ($specs.Count -eq 0) {
    Write-Error "No *.yaml specs in $source; is that the bridge's design/verbs?"
}

New-Item -ItemType Directory -Force -Path $dest | Out-Null

$wanted = @{}
foreach ($spec in $specs) { $wanted[$spec.Name] = $true }
foreach ($stale in @(Get-ChildItem -LiteralPath $dest -Filter '*.yaml' -File)) {
    if (-not $wanted.ContainsKey($stale.Name)) {
        Remove-Item -LiteralPath $stale.FullName
        Write-Host "removed  $($stale.Name) (no longer in the bridge)"
    }
}

# Byte for byte from the bridge's git blob: a Windows checkout with core.autocrlf rewrites line
# endings in the working tree, and a vendored copy must not depend on how the maintainer's git is
# configured. With -WorkingTree the files are read from disk and normalised to LF.
function Get-Bytes([string]$name) {
    if ($useHead) {
        $blob = @(Invoke-Git show "HEAD:./$name")
        if ($LASTEXITCODE -ne 0) { Write-Error "git could not show HEAD:./$name" }
        $text = ($blob -join "`n") + "`n"
    } else {
        $text = [System.IO.File]::ReadAllText((Join-Path $source $name)) -replace "`r`n", "`n"
    }
    return $utf8.GetBytes($text)
}

foreach ($spec in $specs) {
    $target = Join-Path $dest $spec.Name
    $before = if (Test-Path -LiteralPath $target) { [System.IO.File]::ReadAllBytes($target) } else { $null }
    $bytes = Get-Bytes $spec.Name
    [System.IO.File]::WriteAllBytes($target, $bytes)
    $same = $before -and ($before.Length -eq $bytes.Length) -and (-not (Compare-Object $before $bytes -SyncWindow 0))
    if (-not $same) {
        Write-Host "$(if ($before) { 'updated' } else { 'added  ' })  $($spec.Name)"
    }
}

$commit = (Invoke-Git rev-parse HEAD) -join ''
if ($LASTEXITCODE -ne 0 -or -not $commit) { $commit = 'unknown (not a git checkout)' }
$dirty = if ($useHead) { $null } else { (Invoke-Git status --porcelain -- .) -join '' }
$skipped = if ($useHead) { @(Invoke-Git ls-files --others --exclude-standard -- . | Where-Object { $_ -like '*.yaml' }) } else { @() }
foreach ($s in $skipped) { Write-Host "skipped  $(Split-Path -Leaf $s) (not committed in the bridge; -WorkingTree takes it)" }

$lines = @(
    '# Written by scripts/sync-specs.ps1. Do not edit; re-run the script.'
    'repository: trishorts/bridge (private)'
    'path: design/verbs'
    "commit: $commit"
    "uncommitted_changes: $(if ($dirty) { 'yes' } else { 'no' })"
    "synced: $(Get-Date -Format 'yyyy-MM-dd')"
    "specs: $($specs.Count)"
)
[System.IO.File]::WriteAllText((Join-Path $dest 'SOURCE'), (($lines -join "`n") + "`n"), $utf8)

Write-Host "$($specs.Count) specs from $($commit.Substring(0, [Math]::Min(12, $commit.Length)))$(if ($dirty) { ' (with uncommitted changes)' })"
Write-Host 'Now regenerate the reference fragments: $env:MZLIB_RENDER_SPEC_DOCS = 1; cargo test --test spec_docs'
