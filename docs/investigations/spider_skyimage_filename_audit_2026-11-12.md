# Original-sky image filename audit — 12 November 2026

Source: investigator-uploaded EcoVis inventory.json, SHA256
`95f920eb5e2e2f4f47e2c320d69e85071d8385f9f94690b9ab6f9e9ef9a4760e`. Filename listing only; no original TIFF pixels
were available in the coding workspace. This is not a full image-integrity check.

43 directories scanned, no warnings, 20 folders containing numbered TIFFs 1–12,
no duplicate numbered files, plus one Waitrose folder with four TIFFs: 244 image
filenames in 21 site-labelled folders. Retained ZIP files are normal; the old
inventory field archives_not_extracted is misleading and is renamed
archive_files_found in code version 2. It never tested whether ZIPs were extracted.

## Observed folder names

| Folder | TIFFs | Filename audit |
|---|---:|---|
| bristol | 12 | 1–12 present |
| churchill | 12 | 1–12 present |
| cottage | 12 | 1–12 present |
| cricketpitch | 12 | 1–12 present |
| east horrington | 12 | 1–12 present |
| gavarnie | 12 | 1–12 present |
| ounsdale | 12 | 1–12 present |
| oxford_mini | 12 | 1–12 present |
| oxford_museum | 12 | 1–12 present |
| oxford_near_mini | 12 | 1–12 present |
| oxford_open_square | 12 | 1–12 present |
| oxford_stannes | 12 | 1–12 present |
| oxford_worcester_lawn | 12 | 1–12 present |
| oxford_worcester_pitch | 12 | 1–12 present |
| oxford_zoology | 12 | 1–12 present |
| sandhammeren | 12 | 1–12 present |
| skillinge | 12 | 1–12 present |
| slab house | 12 | 1–12 present |
| south horrington | 12 | 1–12 present |
| wells | 12 | 1–12 present |
| waitrose_tesco | 4 | Originally DSC_5080.tif–DSC_5083.tif; user subsequently renamed to 1.tif–4.tif |

## Waitrose confirmation

Investigator explicitly confirmed the renamed files 1,2,3,4 correspond to
0,45,90,135 degrees respectively. It is one actual replicate. Numbered 12-file
folders retain three interleaved quartets (1,4,7,10), (2,5,8,11), (3,6,9,12).
No Waitrose data are duplicated to fill absent replicates. Code now accepts
one to three explicit four-angle groups and includes regression tests for the
single-replicate case and rejecting reused files.

## Original site identifiers

Original CSV filename stems are the durable identifiers requested for folder
renaming (case and spelling shown as archived): ChurchHill, Cottage,
EastHorrington, Gavarnie, minifactory, museum, opensquare, Ounsdale,
oxfordbussinesspark, PennCC, Skillinge, slabhouse, southhorrington, stannes,
Waitrose, WClawn, WCpitch, zoology. The spelling oxfordbussinesspark is original.

Names alone do not settle all mappings: churchill likely corresponds to ChurchHill,
cricketpitch likely to PennCC, and oxford_near_mini may be Oxford Business Park,
but these are not automatically assigned as verified correspondences. Bristol,
Sandhammeren and Wells are additional/unresolved folders relative to the 18-site
reference; do not substitute them for missing named calibration sites.
Investigator will rename site folders to resolve identities. Rerun discovery
after any renaming because previous inventories hold the old paths.

## Implemented follow-up

- Single-replicate TIFF input support with explicit 0/45/90/135 grouping.
- Four-file discovery candidates are flagged as requiring layout confirmation,
  rather than declared complete 12-image sites.
- Read-only TIFF-header inspection across an inventory, optional unaligned
  contact-sheet preview for one folder, no required guessed FOV geometry.
- 24 synthetic tests pass locally, including single-replicate image loading,
  no replicate duplication, header/quicklook preservation of original file
  hashes, and reporting unreadable TIFF headers. No EcoVis field PRC results.

Next: after folder naming is settled, rerun discovery to a new output folder,
inspect headers/unaligned inputs and establish the calibrated 60-degree radius
for a selected original site. No source radiance or locked transfer changes.
