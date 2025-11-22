# inmembrane (modernized fork)

> **Note**  
> This repository is a **modernized fork** of the original `inmembrane` project, tailored for the SerraPHIM phage–host interaction pipeline.  
>  
> The old plugin-based workflow, web-service wrappers, and unit tests have been retired / quarantined. Only the **modern Gram-negative and Gram-positive protocols** using local command-line tools are supported.

## Overview

*inmembrane* is a pipeline for proteome annotation that prioritizes
bacterial proteins which are likely to be **surface-accessible** and
act as **phage receptor candidates**.

It orchestrates the analysis of protein sequences using several
external predictors and summarizes, for each protein:

- membrane topology / subcellular localization,
- evidence for outer-membrane or cell-wall anchoring,
- a boolean *phage-receptor-candidate* flag, and
- supporting details from all tools.

Current protocols:

- `gram_neg_modern` – modern Gram-negative topology & OM receptor logic;
- `gram_pos_modern` – modern Gram-positive cell wall / surface logic.

Typical usage via the CLI script:

```bash
$ inmembrane_scan --config /path/to/inmembrane.config my_proteome.faa
```

## Quick start

1. **Clone and install** (editable mode recommended during development):

   ```bash
   git clone https://github.com/spyridon-ntokos/inmembrane.git
   cd inmembrane
   python -m pip install -e .
   ```

2. **Install external command-line tools** (not Python packages):

   - [SignalP 6.0](https://services.healthtech.dtu.dk/service.php?SignalP-6.0) (`signalp6` binary)
   - [TMbed](https://github.com/BernhoferM/TMbed) (`tmbed` binary)
   - [DeepLocPro](https://github.com/Jaimomar99/deeplocpro) (`deeplocpro` binary)
   - [HMMER 3.x](http://hmmer.org/) (`hmmscan`, `hmmpress`)

   Make sure they are on your `$PATH` or configure explicit paths in
   `inmembrane.config`.

3. **Create or edit the configuration file**.

   By default, `inmembrane` looks for:

   ```text
   ~/SerraPHIM_v2/tools/inmembrane/inmembrane.config
   ```

   If it does not exist, a **default config** is created there the first time
   you run `inmembrane.get_params`.

   You can override the config path with:

   ```bash
   inmembrane_scan --config /path/to/inmembrane.config my_proteome.faa
   ```

4. **Run a Gram-negative analysis** on a Bakta protein FASTA (example):

   ```bash
   inmembrane_scan \
       --config ~/SerraPHIM_v2/tools/inmembrane/inmembrane.config \
       ~/SerraPHIM_v2/data/bakta_annotations/your_sample/your_sample.faa
   ```

   To switch to Gram-positive logic, apply the following edits in `inmembrane.config`:

   ```python
   'protocol': 'gram_pos_modern',
   ...
   'deeplocpro_group': 'positive',
   ```

## Configuration

The configuration file is a simple **Python dict literal** (not JSON):

```python
{
    # Core IO
    'fasta': '',
    'csv': 'out_surfaceome.csv',
    'json': 'out_surfaceome.json',
    'out_dir': '~/SerraPHIM_v2/data/inmembrane_output/gram_neg_test',

    # Protocol: 'gram_neg_modern' or 'gram_pos_modern'
    'protocol': 'gram_neg_modern',
    'phage_receptor_veto': True,

    # SignalP 6.0
    'signalp_bin': '~/.pyenv/versions/serraphim_inmembrane/bin/signalp6',
    'signalp_organism': 'other',
    'signalp_mode': 'fast',
    'signalp_batch': 10,
    'signalp_write_procs': 8,
    'signalp_threads': 8,
    'signalp_skip_cmd': False,

    # TMbed
    'tmbed_bin': 'tmbed',
    'tmbed_device': 'cuda',
    'tmbed_use_gpu': True,
    'tmbed_cpu_fallback': True,
    'tmbed_threads': 16,
    'tmbed_batch_size': 1500,
    'tmbed_batch_size_long': 200,
    'tmbed_max_gpu_len': 4000,
    'tmbed_chunk_overlap': 500,
    'tmbed_skip_cmd': False,
    'tmbed_debug': False,
    'tmbed_delete_embedding': False,

    # Shared TMbed thresholds
    'tmbed_outside_loop_min_pos': 50,
    'tmbed_outside_loop_min_neg': 14,
    'tmbed_outside_frac_thr': 0.30,
    'tmbed_secreted_outside_frac': 0.80,

    # DeepLocPro
    'deeplocpro_bin': 'deeplocpro',
    'deeplocpro_max_gpu_len': 4000,
    'deeplocpro_truncate_len': 6000,
    'deeplocpro_device_short': 'cuda',
    'deeplocpro_device_long': 'cuda',
    'deeplocpro_group': 'negative',
    'deeplocpro_skip_cmd': False,

    # HMMER / Pfam panels
    'hmmer_bin': 'hmmscan',
    'hmmer_db_root': '~/SerraPHIM_v2/tools/hmmer_db',
    'hmmer_cpu': 16,
    'hmmer_evalue_cutoff': 1e-5,
    'hmmer_skip_cmd': False,

    'hmmer_custom_hmm_urls': [
        "https://gitlab.pasteur.fr/gem/capsuledb/-/archive/master/capsuledb-master.tar.gz?ref_type=heads&path=CapsuleFinder_profiles",
    ],

    # Gram- and Gram+ Pfam panels
    'hmmer_pfam_gram_neg': [...],
    'hmmer_pfam_veto_neg': [...],
    'hmmer_pfam_gram_pos': [...],
    'hmmer_pfam_veto_pos': [...],

    # PhReD keyword behavior
    'phred_appendage_pse': True,

    # Order of predictors
    'predictors': ['signalp6', 'tmbed', 'deeplocpro', 'hmmer'],

    # Metadata
    'citations': 'citations.txt',
}
```

Key points:

- The **FASTA path** can be set in the config (`'fasta'`) or passed via CLI.
- `protocol` selects Gram-negative or Gram-positive logic.
- HMMER panels are curated Pfam lists plus optional custom HMM archives.

## External tools

The modern inmembrane pipeline uses only **local command-line tools** (no web
services). You must install and configure these separately:

- **SignalP 6.0** – signal peptides, cleavage sites
- **TMbed** – topology, transmembrane segments, beta-barrels
- **DeepLocPro** – bacterial subcellular localization
- **HMMER 3.x** – domain annotation against curated Pfam / custom panels

### General recommendations

- Install each tool **into its own environment** (e.g. separate `pyenv` / conda
  env) to avoid dependency conflicts (especially different `torch` versions).
- Put the tool binaries (e.g. `signalp6`, `tmbed`, `deeplocpro`, `hmmscan`) on
  your `PATH` or configure absolute paths in `inmembrane.config`.
- Always refer to the official documentation for platform-specific
  installation details.

### SignalP 6.0

1. Request access and download the package from the [official DTU page](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0i&platform=fast).  
2. After receiving the URL to the software package, download and unpack it in a tools directory, for example:

   ```bash
   cd ~/tools
   wget <received_URL>/signalp-6.0i.fast.tar.gz
   tar -xvzf signalp-6.0i.fast.tar.gz
   cd signalp_fast/
   ```

3. Create a dedicated environment (example with `pyenv`):

   ```bash
   pyenv virtualenv 3.10.12 inmembrane_signalp
   pyenv activate inmembrane_signalp
   ```

4. Install the Python package and dependencies, then copy models and convert
   them for GPU use:

   ```bash
   pip install signalp-6-package/
   pip install "numpy<2"

   SIGNALP_DIR=$(python3 -c "import signalp, os; print(os.path.dirname(signalp.__file__))")
   cp -r signalp-6-package/models/* "$SIGNALP_DIR/model_weights/"

   sudo apt-get install zip  # or equivalent for your OS
   signalp6_convert_models gpu "$SIGNALP_DIR/model_weights/"
   ```

5. Verify:

   ```bash
   signalp6 -h
   ```

Configure `signalp_bin` in `inmembrane.config` to point to the `signalp6`
executable if it is not on your `PATH`, for example (if you installed SignalP 
on a different venv):
```bash
'signalp_bin': '~/.pyenv/versions/inmembrane_signalp/bin/signalp6',
```

### DeepLocPro

1. Clone and install:

   ```bash
   cd ~/tools
   git clone https://github.com/Jaimomar99/deeplocpro
   cd deeplocpro
   pip install .
   ```

2. Verify with a small test FASTA:

   ```bash
   deeplocpro -d cuda \
     -f /path/to/test_proteome.faa \
     -o /path/to/deeplocpro_test_output
   ```

Configure `deeplocpro_bin` in `inmembrane.config` if needed.

### TMbed

1. Clone and install:

   ```bash
   cd ~/tools
   git clone https://github.com/BernhoferM/TMbed.git tmbed
   cd tmbed
   pip install .
   tmbed --help
   ```

2. Minimal verification (embedding + prediction on a small FASTA):

   ```bash
   tmbed embed \
     -f /path/to/test.faa \
     -e /path/to/tmbed_embed_test.h5

   tmbed predict \
     -f /path/to/test.faa \
     -e /path/to/tmbed_embed_test.h5 \
     -p /path/to/tmbed_pred_test.pred
   ```

Configure `tmbed_bin` and GPU/CPU options in `inmembrane.config`.

### HMMER 3.x and custom HMMs

1. Install HMMER:

   ```bash
   sudo apt-get install hmmer
   hmmscan -h
   ```

2. Create a directory for your curated HMM panels:

   ```bash
   mkdir -p ~/tools/hmmer_db
   cd ~/tools/hmmer_db
   ```

3. Download selected Pfam models and build a combined database (example with a
   few OM receptor-related Pfams):

   ```bash
   PFAMS=(
     PF00267  # Porin_1
     PF00593  # TonB_dep_Rec
     PF13505  # OMP_b-brl
   )

   for PF in "${PFAMS[@]}"; do
     echo "Downloading $PF..."
     curl -L "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/${PF}?annotation=hmm" \
          -o "${PF}.hmm.gz"
     gunzip -f "${PF}.hmm.gz"
   done

   cat PF*.hmm > gramneg_receptors.hmm
   hmmpress gramneg_receptors.hmm
   ```

4. Test the panel:

   ```bash
   hmmscan --cpu 8 \
     ~/tools/hmmer_db/gramneg_receptors.hmm \
     /path/to/proteome.faa
   ```

Point `hmmer_bin` and `hmmer_db_root` in `inmembrane.config` to your `hmmscan`
binary and HMM database root, and configure panel-specific lists under
`hmmer_pfam_gram_neg`, `hmmer_pfam_gram_pos`, and corresponding veto panels.

## Protocols and categories

### Modern Gram-negative protocol

Implemented in `inmembrane.protocols.gram_neg_modern`.

Uses:

- TMbed (β-strands + α-helices/loops)
- SignalP 6.0
- DeepLocPro
- HMMER (receptor + veto Pfams, custom HMM archives)  
- Optional PhReD lists  

Predicted categories:

- `OM(barrel)`, `OM`
- `LIPOPROTEIN(OM)`, `LIPOPROTEIN(IM)`
- `IM`, `IM(peri+cyto)`
- `PERIPLASMIC`
- `CYTOPLASMIC`
- `SECRETED`

`is_phage_receptor_candidate` depends on:

- OM-like localization  
- TMbed outside loops  
- HMMER evidence  
- Optional PhReD heuristics  

### Modern Gram-positive protocol

Implemented in `inmembrane.protocols.gram_pos_modern`.

Uses:

- SignalP 6.0  
- TMbed  
- DeepLocPro  
- HMMER Pfam panels (LPXTG, LysM, CW_binding, SLH, GW...)  
- Optional annotation-based veto  

Predicted categories:

- `PSE-Cellwall`
- `PSE-Membrane`
- `PSE-Lipoprotein`
- `MEMBRANE(non-PSE)`
- `SECRETED`
- `CYTOPLASMIC`

## Output

Output files:

- **CSV summary**
- **JSON annotation dictionary**
- **citations.txt**

CSV header:

```text
SeqID,Category,LoopExtent,Details,Name
```

Example row:

```text
SPy_0012,PSE-Cellwall,87,
signalp(Sec/SPI;CS=32;Pr=0.99);
tmbed(H=1;LoopOut=87;OutFrac=0.62);
hmmer(PositiveHits=PF00746+PF01476);
PRC=True,"SPy_0012 hypothetical protein"
```

## Using the library API

```python
import inmembrane

params = inmembrane.get_params("/path/to/inmembrane.config")
params["fasta"] = "/path/to/proteome.faa"

proteins = inmembrane.process(params)
```

Each protein entry contains:

- raw sequence fields  
- predictor outputs  
- protocol outputs  

## Extending the codebase

### Predictors

Predictors live in `inmembrane/predictors` and must define:

```python
def annotate(params, proteins):
    ...
```

Predictors:

- read needed fields  
- add fields in-place  
- return updated dict  

### Protocols

Protocols live in `inmembrane/protocols` and define:

- `get_annotations(params)`
- `post_process_protein(params, protein)`
- utilities (`protein_csv_line`, `summary_table`, …)

## Legacy content

Legacy scripts relying on web services (TMHMM, SignalP 4.1, LipoP, MEMSAT3,
mechanize, twill…) are quarantined and not part of the modern workflow.

Legacy repo: https://github.com/boscoh/inmembrane
