# =================================================================================================
#     Dependencies
# =================================================================================================

import pandas as pd
import os, sys, pwd
import socket, platform
import subprocess
from datetime import datetime

snakemake.utils.min_version("7.0")

fissioncnv_version = "0.1.0"

# Store sample names into config['global']
if "global" in config:
    raise Exception("Config key 'global' cannot be defined in config file.")
else:
    config["global"] = {}

config["global"]["samples"] = pd.read_csv(
    config["data"]["samples"] if os.path.isabs(config["data"]["samples"])
    else os.path.join(workflow.basedir, config["data"]["samples"]),
    sep='\t', dtype=str
).set_index(["sample"], drop=False)

config["global"]["sample-names"] = list(config["global"]["samples"]["sample"].unique())

# Validate filenames
if not config['params']['bam-input']:
    for index, row in config["global"]["samples"].iterrows():
        if not os.path.isfile(row['fq1']) or not os.path.isfile(row['fq2']):
            raise Exception(
                "Input fastq files not found: " + str(row["fq1"]) + "; " + str(row["fq2"])
            )

# Clean up reference genome path
if config["data"]["genome"].endswith(".gz"):
    config["data"]["genome"] = os.path.splitext(config["data"]["genome"])[0]

# Store the path of FissionCNV
config['params']['absPath'] = workflow.basedir

# =================================================================================================
#     Pipeline User Output
# =================================================================================================

username = pwd.getpwuid(os.getuid())[0]
hostname = socket.gethostname()

try:
    process = subprocess.Popen(['conda', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    conda_ver = out.decode('ascii').strip().replace("conda ", "")
    del process, out, err
    if not conda_ver:
        conda_ver = "n/a"
except:
    conda_ver = "n/a"

logger.info("===========================================================================")
logger.info("    FissionCNV - CNV calling pipeline for haploid fission yeast")
logger.info("")
logger.info("    Date:               " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
logger.info("    User:               " + username)
logger.info("    Host:               " + hostname)
logger.info("    FissionCNV:         " + fissioncnv_version)
logger.info("    Snakemake:          " + str(snakemake.__version__))
logger.info("    Conda:              " + conda_ver)
logger.info("")
logger.info("    Base directory:     " + workflow.basedir)
logger.info("    Working directory:  " + os.getcwd())
logger.info("    Samples:            " + str(len(config['global']['sample-names'])))
logger.info("    Ploidy:             " + str(config['params']['ploidy']))
logger.info("    Bin size:           " + str(config['params']['binSize']))
logger.info("===========================================================================")
logger.info("")

del username, hostname, conda_ver

# =================================================================================================
#     Helper Functions
# =================================================================================================

def get_sample_bam(samples):
    return ["mapped/" + sample + ".bam" for sample in samples]

def get_sample_bai(samples):
    return ["mapped/" + sample + ".bam.bai" for sample in samples]
