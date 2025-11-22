# inmembrane changelog

## v0.96.0-dev (22-Nov-2025)

* Modernized codebase for Python 3.8+ only.
* Quarantined legacy `plugins/` and `tests/` into
  `inmembrane/quarantine_legacy_scripts` (no longer used in the main workflow).
* Replaced legacy predictors (TMHMM, LipoP, SignalP 4, MEMSAT3, web-scrapers)
  with a new, fully local predictor stack:
  - SignalP 6.0 (Sec/Tat signal peptides)
  - TMbed (α-helical TM segments, β-barrels, signal peptides, loop topology)
  - DeepLocPro (Gram-positive / Gram-negative subcellular localization)
  - HMMER 3 + curated Pfam and custom HMM panels (receptor & veto sets)
* Added new Gram-specific protocols:
  - `gram_neg_modern` – OM β-barrels, OM lipoproteins, IM/periplasmic classes,
    and refined phage receptor candidate logic.
  - `gram_pos_modern` – PSE-Membrane / PSE-Cellwall / PSE-Lipoprotein logic
    for thick peptidoglycan cell envelopes, with Gram+ specific Pfam panels.
* Introduced a unified “phage receptor candidate” flag and helper logic shared
  between Gram− and Gram+ protocols.
* Simplified configuration via an updated `inmembrane.config`:
  - explicit paths for external binaries (SignalP6, TMbed, DeepLocPro, hmmscan)
  - shared TMbed exposure thresholds for both protocols
  - Gram-specific Pfam positive/veto panels, custom HMM URL list
* Cleaned up packaging:
  - trimmed legacy web-scraping and SOAP-related Python dependencies
  - aligned `requirements.txt` with the modern codebase
  - refreshed `setup.py`, `MANIFEST.in`, Dockerfile, and README.

## v0.95.0 (24-Jan-2018)
* Resurrected after a period of neglect.
* Removed deprecated CBS-DTU SOAP and TMBETA-NET plugins.
* Fixed TMBETADISC-RBF plugin (pinned twill version).
* Fixed lipop_scrape_web, tmhmm_scrape_web, added new signalp_scrape_web plugins (now used as defaults).
* Added Dockerfile(s).
* Enforce semantic versioning.

## v0.94 (27-Sep-2013)
* Fixed LipoP and TMHMM web scrapers due to server side changes.
* TMB-HUNT plugin deprecated - server appears to be permanently offline
* Added result_poll_retries option to signalp_web, so test fails fast if job is stuck in PENDING state
* Made 'requests' v2.0.0 the lowest version for new installs
* Made default inmembrane.config use SOAP web services (*_web), included commented _scrape_web options in file.

## v0.93.2 (03-Dec-2012)
* Bugfixes to tmhmm_scrape_web and lipop_scrape_web plugins
* Handle Python 2.6 when OrderedDict is absent.

## v0.93 (03-Dec-2012)

* Added --version flag
* Report PSE loop length in gram_pos protocol
* Added tmhmm_scrape_web and lipop_scrape_web plugins (as default since CBS SOAP service seems currently broken)
* Added safe_seqid to be used by programs that munge or fail on seqids with funky punctuation

## v0.92 (22-Oct-2012)

* Started CHANGELOG.
