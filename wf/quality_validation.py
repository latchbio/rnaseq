"""
Using this paper as a baseline: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8
"""

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Union
import glob
from io import StringIO
import subprocess
import os
import shutil

import pandas as pd

from wf import PairedEndReads, SingleEndReads


class ValidationStage(Enum):
    RAW = 1
    TRIMMED = 2
    ALIGNED = 3


class FastQCReportModules(str, Enum):
    BASIC_STATISTICS = "Basic Statistics"
    PER_BASE_SEQ_QUALITY = "Per base sequence quality"
    PER_SEQ_QUALITY_SCORE = "Per sequence quality scores"
    PER_BASE_SEQ_CONTENT = "Per base sequence content"
    PER_SEQ_GC_CONTENT = "Per sequence GC content"
    PER_BASE_N_CONTENT = "Per base N content"
    SEQ_LENGTH_DIST = "Sequence Length Distribution"
    TOTAL_DEDUPLICATED_PCT = "Total Deduplicated Percentage"
    SEQ_DUP_LEVELS = "Sequence Duplication Levels"
    OVERREP_SEQ = "Overrepresented sequences"
    ADAPTER_CONTENT = "Adapter Content"

    END_MODULE = ">>END_MODULE"


@dataclass()
class FastQCReport:
    total_sequences: int
    poor_quality_sequences: int
    gc_pct: float
    total_deduplicated_pct: float

    per_base_seq_quality_scores: pd.DataFrame
    per_seq_quality_scores: pd.DataFrame
    per_base_seq_content: pd.DataFrame
    per_seq_gc_content: pd.DataFrame
    per_base_n_content: pd.DataFrame

    seq_length_dist: pd.DataFrame
    seq_duplication_levels: pd.DataFrame
    overrepresented_seq: pd.DataFrame
    adapter_content: pd.DataFrame


def parse_fastqc_report(path_to_fastqc_data: str) -> FastQCReport:
    raw_text = Path(path_to_fastqc_data).read_text()
    modules = raw_text.split(FastQCReportModules.END_MODULE)

    report = {}
    for i, mod in enumerate(modules):
        if mod == "\n":
            continue

        if i == 0:
            # parse initial module differently
            x = mod.split("#")
            x = x[2:]  # remove version number
            x[0] = x[0].split("\n")[-2]

        else:
            x = mod.split("#")

        module_title = x[0].split("\t")[0].replace(">>", "").strip("\n ")

        if FastQCReportModules.SEQ_DUP_LEVELS in module_title:
            report[FastQCReportModules.TOTAL_DEDUPLICATED_PCT] = float(
                x[1].split("\t")[1].strip()
            )
            report[module_title] = pd.read_csv(StringIO(x[2]), sep="\t")

        else:
            report[module_title] = pd.read_csv(StringIO(x[1]), sep="\t")

    total_sequences = int(report[FastQCReportModules.BASIC_STATISTICS]["Value"].iloc[3])
    poor_quality_sequences = int(
        report[FastQCReportModules.BASIC_STATISTICS]["Value"].iloc[4]
    )
    gc_pct = float(report[FastQCReportModules.BASIC_STATISTICS]["Value"].iloc[6])

    return FastQCReport(
        total_sequences,
        poor_quality_sequences,
        gc_pct,
        report[FastQCReportModules.TOTAL_DEDUPLICATED_PCT],
        report[FastQCReportModules.PER_BASE_SEQ_QUALITY],
        report[FastQCReportModules.PER_SEQ_QUALITY_SCORE],
        report[FastQCReportModules.PER_BASE_SEQ_CONTENT],
        report[FastQCReportModules.PER_SEQ_GC_CONTENT],
        report[FastQCReportModules.PER_BASE_N_CONTENT],
        report[FastQCReportModules.SEQ_LENGTH_DIST],
        report[FastQCReportModules.SEQ_DUP_LEVELS],
        report[FastQCReportModules.OVERREP_SEQ],
        report[FastQCReportModules.ADAPTER_CONTENT],
    )


def run_fastqc(path_to_fastq: str, output_dir: str) -> FastQCReport:
    os.makedirs(output_dir, exist_ok=True)

    cmd = ["fastqc", path_to_fastq, "-o", output_dir]

    subprocess.run(cmd, check=True)

    # Unzip the raw report and parse
    zip_files = glob.glob(os.path.join(output_dir, "*.zip"))
    shutil.unpack_archive(zip_files[0], output_dir)
    return parse_fastqc_report(os.path.join(output_dir, "fastqc_data.txt"))


def validate_raw_replicate(replicate: Union[PairedEndReads, SingleEndReads]):
    """
    Raw reads are evaluated on:
        - GC Content
        - Duplicated Reads
        - Sequence Quality
    """

    if type(replicate) == PairedEndReads:
        # Each fastq file in a PairedEnd read is analyzed separately

        pass

    else:
        # TODO
        pass
