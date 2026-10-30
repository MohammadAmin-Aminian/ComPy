# From-scratch pipeline: recompute everything with one command

`scripts/run_from_scratch.py` rebuilds every BRUIT-FM product from the raw
archive into one new folder. It covers all four networks (YV, YL, 7D, XO),
PPSDs for all channels with usable responses, the three spectral checks, and
the three calibration methods. Rayleigh runs with both event policies.

The owner chose these settings on 22 October 2026:

* all four networks;
* local saved metadata only, so nothing is downloaded;
* all channels;
* both Rayleigh event policies.

The settings live in [`from_scratch_config.json`](from_scratch_config.json),
and each one reproduces a completed, reviewed run.

## What it never does

* Modify raw waveforms, StationXML, `exports/` masters, frozen results or saved
  products. Those are only read.
* Download metadata or catalogues. The 7D/XO preparations run with
  `--offline`, so a missing local copy stops the run.
* Apply a gain, or set a threshold, score or rating.

## Stages

| # | Stage | What it does | Written to |
|---|---|---|---|
| 1 | `preflight` | Committed code (refuses uncommitted `src/`/`scripts/` changes), Python 3.13, imports, every configured input present, free disk | `run.json` |
| 2 | `workspaces` | Per-network calibration master from local StationXML (YV, YL), or the 7D/XO preparation fed the saved responses and catalogues. Then one standalone response snapshot per usable channel epoch | `calibration/<NET>/` |
| 3 | `identity` | YV Rayleigh with bounded events. **Must return RR52 g = 0.6787278140712512 exactly, or the run stops**, before the long stages | `calibration/YV/rayleigh_bounded/` |
| 4 | `coverage` | Fresh coverage campaign with reuse disabled: PPSDs for all channels, plus the expected-pressure, nearby and infragravity checks | `coverage/campaign/` |
| 5 | `rayleigh` | Every network with `bounded` and `all` events | `calibration/<NET>/rayleigh_<policy>/` |
| 6 | `microseism` | Common-week pressure/vertical check on the fresh masters (YV week 2013-06-10), then the method comparison | `microseism/run_*/` |
| 7 | `local_p` | Reviewed v2 local-P campaign, station table and export | `results/LOCAL_P/campaigns/from_scratch_<run>/` |
| 8 | `registry` | Master-style CSV of every new PPSD, for the plot scripts | `registry/master_inventory.csv` |
| 9 | `summary` | Stage status and Rayleigh comparison with the frozen workflows | `summary.json`, `SUMMARY.md` |

The run folder is `results/FROM_SCRATCH/run_<UTC>_<id>/`. Every command's
output is in `logs/`. If a stage fails, it is recorded and the independent
stages still run; only a failed `workspaces` or `identity` stops the run. The
code state is checked before every stage, so an edit made during the run stops
it.

## How to run (new computer, from `~/Projects/bruit-fm`, env `bruitfm`)

If the archive lives on a separate volume (for example the 20 TB storage),
link it once instead of changing any path. Every tool and saved path then
resolves through the link:

```bash
ln -s /Volumes/<volume>/Data_Bruit_FM ~/Desktop/Data_Bruit_FM
```

```bash
conda activate bruitfm
cd ~/Projects/bruit-fm
python scripts/check_migration.py            # must say READY
python scripts/run_from_scratch.py --check   # preflight only; writes nothing
```

**Small end-to-end test first** (one station, every stage):

```bash
python scripts/run_from_scratch.py --networks YV --stations RR52 --workers 2
```

**Everything:**

```bash
caffeinate -i python scripts/run_from_scratch.py --workers 6
```

`caffeinate -i` keeps the Mac awake. Set `--workers` to about the number of
performance cores minus two.

**After an interruption**, resume the same run. Completed stages are skipped,
and the PPSD and calibration steps resume their own saved progress:

```bash
python scripts/run_from_scratch.py --run-root results/FROM_SCRATCH/run_<id>
```

Other options:

* `--networks YV YL` limits the networks.
* `--only coverage` or `--skip local_p` selects stages.
* `--config FILE` uses another settings file.
* `--allow-dirty` runs uncommitted code. The run records that it did.

## Results

* `SUMMARY.md`: stage table, and for each network the Rayleigh workflow with
  the frozen policy compared station by station with the frozen run.
* Plots of the new PPSDs use the registry, for example:

  ```bash
  python scripts/plot_network_median_ppsd.py --master results/FROM_SCRATCH/run_<id>/registry/master_inventory.csv --network YV --channel BDH
  python -m bruitfm.plot_network_ppsd --master-csv results/FROM_SCRATCH/run_<id>/registry/master_inventory.csv --network YV --channels BDH --probability
  ```

* Coverage tables and heat maps:
  `python -m bruitfm.coverage report --campaign-dir results/FROM_SCRATCH/run_<id>/coverage/campaign --config results/FROM_SCRATCH/run_<id>/coverage/coverage_sources.json`.
* Repeated microseism weeks are optional. They reuse the fresh sources
  automatically:
  `bash scripts/run_microseism_weeks.sh results/FROM_SCRATCH/run_<id>/microseism/run_<…>`.

## Known limits (data, not code)

* **7D horizontals** cannot be computed offline. The only local full 7D
  responses are the saved P/Z snapshots, and `INVENTORIES/7D.xml` has no
  response stages. They are reported as blocked.
* **Blocked epochs** (7D BXH/LDI, XO LA33 HDH) need response metadata that is
  not available locally.
* **RR12's** archived pressure is flat at −3077 dB. YV October–December 2012 is
  missing from the archive.
* The local-P export joins against `exports/master_inventory.csv`, read-only,
  as the frozen campaign did.
* Ratings, thresholds and the noise model remain separate owner decisions.
