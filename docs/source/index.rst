.. _micaflow-home:

.. raw:: html

    <div class="landing-page">
      <div class="landing-header">
        <img src="_static/images/logo.png" alt="Micaflow Logo" class="landing-logo" />
        <h1>Micaflow 2.0</h1>
        <p class="lead">
          A comprehensive neuroimaging processing pipeline for structural and diffusion MRI analysis
        </p>
        
        <div class="landing-buttons">
          <a href="https://github.com/ian-goodall-halliwell/micaflow2.0" class="landing-button">
            GitHub Repository
          </a>
          <a href="scripts/index.html" class="landing-button">
            API Documentation
          </a>
        </div>
      </div>
    </div>

Overview
========

Micaflow2.0 is a modular neuroimaging pipeline designed for processing structural and diffusion MRI data. Built on Snakemake workflow management system, it provides reproducible and scalable neuroimaging analyses with robust preprocessing, registration, and feature extraction capabilities.

.. raw:: html

    <div class="feature-grid">
      <div class="feature-card">
        <h3>Structural MRI Processing</h3>
        <ul>
          <li>Brain extraction (HD-BET)</li>
          <li>N4 bias field correction</li>
          <li>SynthSeg segmentation</li>
          <li>MNI space registration</li>
        </ul>
      </div>
      
      <div class="feature-card">
        <h3>Diffusion MRI Processing</h3>
        <ul>
          <li>Motion correction</li>
          <li>Denoising with Patch2Self</li>
          <li>FA/MD map calculation</li>
          <li>Distortion correction (HYSCO)</li>
        </ul>
      </div>
      
      <div class="feature-card">
        <h3>Analysis Tools</h3>
        <ul>
          <li>Texture feature extraction</li>
          <li>Registration quality metrics</li>
          <li>ANTs-based transformations</li>
          <li>Automated workflow with Snakemake</li>
        </ul>
      </div>
    </div>

Getting Started
==============

Installation
-----------

.. code-block:: bash

    # Clone the repository
    git clone https://github.com/ian-goodall-halliwell/micaflow2.0.git
    cd micaflow2.0

    # Create and activate a conda environment
    conda env create -f environment.yml
    conda activate micaflow

    # Install additional dependencies
    pip install -r requirements.txt

Basic Usage
----------

.. code-block:: bash

    # Edit the config.yaml file to set your input and output directories
    # Run the pipeline
    snakemake --cores 4

    # Processing a single subject
    export SUBJECT=sub-01
    export SESSION=ses-01
    export DATA_DIRECTORY=/path/to/data
    export OUT_DIR=/path/to/output
    snakemake --cores 4 -C subject=$SUBJECT session=$SESSION

Documentation
============

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   scripts

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`