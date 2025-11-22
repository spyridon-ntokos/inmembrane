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

   To switch to Gram-positive logic, edit in `inmembrane.config`:

   ```python
   'protocol': 'gram_pos_modern',
   ```

## Configuration

The configuration file is a simple **Python dict literal** (not JSON):

```python
{
    # Core IO
    'fasta': '',
    'csv': 'out_surfaceome.csv',
    'json': 'out_surfaceome.json',
    'out_dir': '/home/snt/SerraPHIM_v2/data/inmembrane_output/gram_neg_test',

    # Protocol: 'gram_neg_modern' or 'gram_pos_modern'
    'protocol': 'gram_neg_modern',
    'phage_receptor_veto': True,

    # SignalP 6.0
    'signalp_bin': '/home/snt/.pyenv/versions/serraphim_inmembrane/bin/signalp6',
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
    'hmmer_db_root': '/home/snt/SerraPHIM_v2/tools/hmmer_db',
    'hmmer_cpu': 16,
    'hmmer_evalue_cutoff': 1e-5,
    'hmmer_skip_cmd': False,

    'hmmer_custom_hmm_urls': [
        "https://gitlab.pasteur.fr/gem/capsuledb/-/archive/master/"
        "capsuledb-master.tar.gz?ref_type=heads&path=CapsuleFinder_profiles",
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

Tools used:

- **SignalP 6.0**
- **TMbed**
- **DeepLocPro**
- **HMMER 3.x**

Debian/Ubuntu example:

```bash
sudo apt-get install hmmer
hmmscan -h
```

TMbed / DeepLocPro / SignalP 6.0 require manual installation per tool docs.

## Protocols and categories

### Modern Gram-negative protocol

Implemented in `inmembrane.protocols.gram_neg_modern`.

Uses:

- TMbed β-strands + DeepLocPro → OM barrels  
- TMbed helices/loops  
- SignalP 6.0  
- DeepLocPro  
- HMMER (receptor Pfams + veto Pfams)  
- Optional PhReD heuristics  

Example categories:

- `OM(barrel)`, `OM`
- `LIPOPROTEIN(OM)`, `LIPOPROTEIN(IM)`
- `IM`, `IM(peri+cyto)`
- `PERIPLASMIC`, `SECRETED`
- `CYTOPLASMIC`

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

Example categories:

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

Fields:

- **SeqID**
- **Category**
- **LoopExtent**
- **Details**
- **Name**

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
