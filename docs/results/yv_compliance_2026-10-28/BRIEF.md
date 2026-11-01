# Task brief: YV broadband seafloor compliance with ComPy (2026-10-28)

## Goal
Audit and correct ComPy 2.1, then compute the seafloor compliance function, with
uncertainties and quality flags, for the YV (RHUM-RUM) broadband OBS stations.
The work stops at compliance. There is no inversion.

## Stations
RR28, RR29, RR36, RR38, RR40, RR50, RR52. The data are in /Volumes/DATA_bruitfm/SDS/{2012,2013}/YV and /Volumes/DATA_bruitfm/INVENTORIES/YV.xml.
RR32 was requested, but it is not on the drive, so it is skipped by the owner's decision.

## Decisions (owner, 2026-10-28)
- **DPG gain.** Compliance is computed twice and the results are compared:
  - (a) with the gains in `gain_corrected_2026-10-28/gains.csv`;
  - (b) with ComPy's own teleseismic DPG calibration.
- **Code scope.**
  - Clear defects are fixed, each with a test and a CHANGELOG entry.
  - Methodological alternatives become options that keep the old default, and they are compared in the report.
  - Wayne's conventions apply: explicit Z/P channels, no magic numbers, Google docstrings.
- **Downloads.** Read-only downloads of papers, code, catalogues and metadata are allowed. Nothing is posted.
- **Delivery.**
  - The code goes to ComPy branch `claude/compliance-yv` with a draft PR, which is not merged.
  - The report is published as a private page.
  - An email goes to the owner only.
- **Data safety.** Frozen results are never overwritten. All outputs go in this folder.

## Phases
1. **Study.** Read Aminian et al. 2025 (GJI), the thesis, the Crawford and Zha papers, and TiSKitPy, ATaCR and OBStools.
2. **Audit.** Check ComPy against that theory: units, signs, the infragravity cutoff, coherence, uncertainty, response, tilt, day selection and magic numbers.
3. **Correct.** Fix the defects found in the audit, with tests, and compare RR52 before and after.
4. **Process.** For each station:
   - QC
   - response removal
   - DPG gain
   - transient and earthquake removal
   - tilt and compliance-noise correction
   - compliance with uncertainty
   - flags
5. **Deliver.** Write the report and figures, push the branch, open the draft PR and send the email.
