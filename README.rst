=============================
inmembrane (modernized fork)
=============================

.. note::

   This repository is a **modernized fork** of the original
   :mod:`inmembrane` project, tailored for the SerraPHIM phage–host
   interaction pipeline.

   The old plugin-based workflow, web-service wrappers, and unit tests
   have been retired / quarantined. Only the **modern Gram-negative and
   Gram-positive protocols** using local command-line tools are
   supported.

Overview
========

*inmembrane* is a pipeline for proteome annotation that prioritizes
bacterial proteins which are likely to be **surface-accessible** and
act as **phage receptor candidates**.

It orchestrates the analysis of protein sequences using several
external predictors and summarizes, for each protein:

* membrane topology / subcellular localization,
* evidence for outer-membrane or cell-wall anchoring,
* a boolean *phage-receptor-candidate* flag, and
* supporting details from all tools.

Current protocols:

* ``gram_neg_modern`` – modern Gram-negative topology & OM receptor logic;
* ``gram_pos_modern`` – modern Gram-positive cell wall / surface logic.

Typical usage via the CLI script::

    $ inmembrane_scan --config /path/to/inmembrane.config my_proteome.faa


Quick start
===========

1. **Clone and install** (editable mode recommended during development)::

       git clone https://github.com/spyridon-ntokos/inmembrane.git
       cd inmembrane
       python -m pip install -e .

2. **Install external command-line tools** (not Python packages):

   * `SignalP 6.0`_  (``signalp6`` binary)
   * `TMbed`_        (``tmbed`` binary)
   * `DeepLocPro`_   (``deeplocpro`` binary)
   * `HMMER 3.x`_    (``hmmscan``, ``hmmpress``)

   Make sure they are on your ``$PATH`` or configure explicit paths in
   ``inmembrane.config``.

3. **Create or edit the configuration file**.

   By default, :mod:`inmembrane` looks for::

       ~/SerraPHIM_v2/tools/inmembrane/inmembrane.config

   If it does not exist, a **default config** is created there the first time
   you run :func:`inmembrane.get_params`.

   You can override the config path with::

       inmembrane_scan --config /path/to/inmembrane.config my_proteome.faa

4. **Run a Gram-negative analysis** on a Bakta protein FASTA (example)::

       inmembrane_scan \
           --config ~/SerraPHIM_v2/tools/inmembrane/inmembrane.config \
           ~/SerraPHIM_v2/data/bakta_annotations/your_sample/your_sample.faa

   To switch to Gram-positive logic, edit in ``inmembrane.config``::

       'protocol': 'gram_pos_modern',


Configuration
=============

The configuration file is a simple **Python dict literal** (not JSON):

.. code-block:: python

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
       'deeplocpro_group': 'negative',   # 'negative' for Gram-, 'positive' for Gram+
       'deeplocpro_skip_cmd': False,

       # HMMER / Pfam panels
       'hmmer_bin': 'hmmscan',
       'hmmer_db_root': '/home/snt/SerraPHIM_v2/tools/hmmer_db',
       'hmmer_cpu': 16,
       'hmmer_evalue_cutoff': 1e-5,
       'hmmer_skip_cmd': False,

       'hmmer_custom_hmm_urls': [
           # e.g. CapsuleFinder profiles
           "https://gitlab.pasteur.fr/gem/capsuledb/-/archive/master/"
           "capsuledb-master.tar.gz?ref_type=heads&path=CapsuleFinder_profiles",
       ],

       # Gram- and Gram+ Pfam panels (positive / veto)
       'hmmer_pfam_gram_neg': [...],
       'hmmer_pfam_veto_neg': [...],
       'hmmer_pfam_gram_pos': [...],
       'hmmer_pfam_veto_pos': [...],

       # PhReD keyword behavior
       'phred_appendage_pse': True,

       # Order of predictors to run
       'predictors': ['signalp6', 'tmbed', 'deeplocpro', 'hmmer'],

       # Metadata / bookkeeping
       'citations': 'citations.txt',
   }

Key points:

* The **FASTA path** can be set in the config (``'fasta'``) or passed as
  the positional argument to ``inmembrane_scan``; the CLI argument wins
  if the config value is empty.
* ``protocol`` selects between modern Gram-negative or Gram-positive
  logic.
* HMMER panels are **small curated Pfam lists** plus optional custom HMM
  archives (e.g. CapsuleFinder).


External tools
==============

The modern pipeline uses only **local command-line tools** (no web
services):

* **SignalP 6.0** (Sec/Tat signal peptides, cleavage sites)
* **TMbed** (TM alpha-helices, beta-barrels, signal peptides, topology)
* **DeepLocPro** (bacterial subcellular localization)
* **HMMER 3.x** (hmmscan against curated Pfam & custom panels)

Example installation snippets for a Debian/Ubuntu environment
(adapt as needed):

* HMMER::

    sudo apt-get install hmmer
    hmmscan -h

* TMbed / DeepLocPro / SignalP 6.0:

  See the respective upstream documentation. In a SerraPHIM-style setup,
  each tool lives under ``~/SerraPHIM_v2/tools/<tool>`` and is installed
  into a dedicated Python/conda environment to avoid conflicts (for
  example, SignalP 6.0 on a separate ``pyenv`` than other SerraPHIM tools).


Protocols and categories
========================

Modern Gram-negative protocol
-----------------------------

Implemented in ``inmembrane.protocols.gram_neg_modern``.

Uses:

* TMbed β-strands + DeepLocPro → OM barrel detection
* TMbed helices/loops → IM exposure, periplasmic loops
* SignalP 6.0 → Sec/Tat signal peptides, lipoproteins
* DeepLocPro → outer_membrane, periplasmic, cytoplasmic_membrane, etc.
* HMMER → curated **OM receptor Pfams** and **veto Pfams**
* Optional PhReD keyword heuristics for appendages / flagella (recommended if 
  genomes were annotated with Bakta, so the annotations are included in the 
  "Name" column of the output CSV).

Example high-level categories:

* ``OM(barrel)``, ``OM``
* ``LIPOPROTEIN(OM)``, ``LIPOPROTEIN(IM)``
* ``IM``, ``IM(peri+cyto)``
* ``PERIPLASMIC``, ``SECRETED``
* ``CYTOPLASMIC``

A boolean ``is_phage_receptor_candidate`` flag is computed based on:

* OM-like category / localization,
* **structural exposure** (TMbed outside loops),
* HMMER evidence (receptor Pfams, CapsuleFinder-like models), and
* optional PhReD appendage heuristics.

Modern Gram-positive protocol
-----------------------------

Implemented in ``inmembrane.protocols.gram_pos_modern``.

Uses:

* SignalP 6.0 (including lipoproteins)
* TMbed (membrane segments, exposed loops, outside fraction)
* DeepLocPro (cell_wall, cell_wall_&_surface, extracellular, etc.)
* HMMER Pfam panels enriched for **LPXTG / LysM / CW_binding / SLH / GW**
  and veto Pfams (ribosomal, core metabolism...)
* Optional annotation-based veto to remove obvious phage structural
  proteins and core transporters from receptor candidates.

Example high-level categories:

* ``PSE-Cellwall``
* ``PSE-Membrane``
* ``PSE-Lipoprotein``
* ``MEMBRANE(non-PSE)``
* ``SECRETED``
* ``CYTOPLASMIC``

Here **PSE** stands for *potentially surface-exposed* based on topology,
while the **final receptor-candidate decision** is stored separately as
``is_phage_receptor_candidate``.


Output
======

For each run, :mod:`inmembrane` writes:

* A **CSV summary**, one row per protein
* A **JSON** file with the full per-protein annotation dictionary
* A simple **citations.txt** file summarizing which predictors were used

By default, the CSV header is::

    SeqID,Category,LoopExtent,Details,Name

Example (one line, wrapped for clarity)::

    SPy_0012,PSE-Cellwall,87,
    signalp(Sec/SPI;CS=32;Pr=0.99);
    tmbed(H=1;LoopOut=87;OutFrac=0.62);
    hmmer(PositiveHits=PF00746+PF01476);
    PRC=True,"SPy_0012 hypothetical protein"

Columns:

* **SeqID** – sequence ID (first token from FASTA header)
* **Category** – protocol-specific class (e.g. ``OM(barrel)`` or
  ``PSE-Cellwall``)
* **LoopExtent** – longest TMbed **outside loop** for that protein
* **Details** – semicolon-delimited summary of predictor evidence, including 
  the receptor-candidate flag (e.g. ``PRC=True``, where PRC is "Phage Receptor Candidate")
* **Name** – full FASTA header (minus ``>``), as parsed by the pipeline


Using the library API
=====================

Although :mod:`inmembrane` is primarily used via the CLI, it can also be
used as a Python library::

  import inmembrane

  params = inmembrane.get_params("/path/to/inmembrane.config")
  params["fasta"] = "/path/to/proteome.faa"

  proteins = inmembrane.process(params)

Here ``proteins`` is a dictionary keyed by sequence ID, where each
value is a dict containing:

* raw sequence / metadata (``seq``, ``name``, ``sequence_length``, ...)
* predictor outputs (``signalp_*``, ``tmbed_*``, ``deeplocpro_*``,
  ``hmmer_*``)
* protocol outputs (``category``, ``loop_extent``,
  ``is_phage_receptor_candidate``, ``details``)


Extending the codebase
======================

Modern predictors
-----------------

*Predictors* live in ``inmembrane/predictors`` and must expose::

    def annotate(params, proteins): ...

where:

* ``params`` is the configuration dictionary (from ``inmembrane.config``)
* ``proteins`` is the main dict keyed by sequence ID

Predictors are discovered automatically by
``inmembrane.predictors.get_predictor`` and are listed in
``params["predictors"]`` to control execution order (if it matters).

Each predictor should:

* read only the fields it needs from ``params`` and ``proteins``,
* add new fields to the per-protein dict in place, and
* return the (possibly modified) ``proteins`` dict.


Protocols
---------

*Protocols* live in ``inmembrane/protocols`` and implement the
high-level Gram± logic. A protocol module must define at least:

* ``get_annotations(params)`` – return ordered predictor IDs to run
* ``post_process_protein(params, protein)`` – inspect predictor outputs,
  set ``category``, ``is_phage_receptor_candidate``, ``details``,
  ``loop_extent``, etc.
* utility functions such as ``protein_csv_line`` and ``summary_table``
  used by :func:`inmembrane.process`.

Legacy content
==============

The original plugin-based scripts and test suite that relied on web
services (TMHMM, SignalP 4.1, LipoP, MEMSAT3, mechanize, twill, etc.)
have been removed under a **quarantine** directory and are not part of
the modern workflow.

If you need to reproduce the exact legacy environment, consult the
original upstream repository: https://github.com/boscoh/inmembrane

.. _SignalP 6.0: https://services.healthtech.dtu.dk/service.php?SignalP-6.0
.. _TMbed: https://github.com/BernhoferM/TMbed
.. _DeepLocPro: https://github.com/Jaimomar99/deeplocpro
.. _HMMER 3.x: http://hmmer.org/
