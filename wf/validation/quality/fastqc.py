from dataclasses import dataclass
from enum import Enum
from pathlib import Path
import glob
from io import StringIO
import subprocess
import os
import shutil
from typing import List

import pandas as pd


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


@dataclass
class FastQCReport:
    filename: str
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

        x = mod.split("#")

        if i == 0:
            # parse initial module differently
            x = x[2:]  # remove version number
            x[0] = x[0].split("\n")[-2]

        module_title = x[0].split("\t")[0].replace(">>", "").strip("\n ")

        if FastQCReportModules.SEQ_DUP_LEVELS in module_title:
            report[FastQCReportModules.TOTAL_DEDUPLICATED_PCT] = float(
                x[1].split("\t")[1].strip()
            )
            report[module_title] = pd.read_csv(StringIO(x[2]), sep="\t")

        else:
            report[module_title] = pd.read_csv(StringIO(x[1]), sep="\t")

    filename = report[FastQCReportModules.BASIC_STATISTICS]["Value"].iloc[0]
    total_sequences = int(report[FastQCReportModules.BASIC_STATISTICS]["Value"].iloc[3])
    poor_quality_sequences = int(
        report[FastQCReportModules.BASIC_STATISTICS]["Value"].iloc[4]
    )
    gc_pct = float(report[FastQCReportModules.BASIC_STATISTICS]["Value"].iloc[6])

    return FastQCReport(
        filename,
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

    cmd = ["fastqc", path_to_fastq, "-o", output_dir]

    subprocess.run(cmd, check=True)

    # Unzip the raw report and parse
    zip_files = glob.glob(os.path.join(output_dir, "*.zip"))
    shutil.unpack_archive(zip_files[0], output_dir)
    dirname = Path(zip_files[0]).name.split(".")[0]

    return parse_fastqc_report(os.path.join(output_dir, dirname, "fastqc_data.txt"))


@dataclass
class FastQCValidationConfig:
    """
    Default values from:
        - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/
    """

    check_poor_quality: bool
    check_overall_gc_content: bool
    check_dedup_pct: bool
    check_adapter_content: bool
    check_sequence_length_dist: bool

    quality_threshold: int = 0
    dedup_pct_threshold: float = 80.0
    adapter_pct_threshold: float = 5.0


def check_poor_quality(report: FastQCReport, threshold: int = 0) -> List[str]:
    if report.poor_quality_sequences > threshold:
        return [
            f"> Found {report.poor_quality_sequences} sequences with bad quality scores."
        ]
    return []


def check_dedup_pct(report: FastQCReport, pct_threshold: float = 80) -> List[str]:
    res = []
    if report.total_deduplicated_pct <= pct_threshold:
        res.append(
            f"> Warning: after deduplicating, only {report.total_deduplicated_pct}% of sequences will remain."
        )

        res.append(
            f"> Top 5 Duplication Levels in the file (by % deduplicated) are: \n"
        )
        report.seq_duplication_levels.sort_values(
            "Percentage of deduplicated", ascending=False, inplace=True
        )

        res.append(report.seq_duplication_levels.head(n=5).to_string(index=False))
        res.append("\n")
    return res


def check_overall_gc_content(report: FastQCReport) -> List[str]:
    return [f"> Overall GC Content in reads: {report.gc_pct}%"]


def check_adapter_content(
    report: FastQCReport, adapter_pct_threshold: float = 5.0
) -> List[str]:
    res = []
    zero_resolution = 1e-10
    non_zero_adapters = report.adapter_content
    cols = list(non_zero_adapters.columns)
    cols.remove("Position")
    cols_with_nonzero_sum = [
        x for x in cols if non_zero_adapters[x].sum() >= zero_resolution
    ]

    non_zero_adapters = non_zero_adapters.loc[
        ~non_zero_adapters.apply(
            lambda row: (row[cols] <= zero_resolution).all(), axis=1
        )
    ]

    if non_zero_adapters.shape[0] > 0:
        res.append(f"> Found adapters in reads: \n")
        res.append(
            non_zero_adapters[["Position", *cols_with_nonzero_sum]].to_string(
                index=False
            )
        )
        res.append("\n")

    gt_pct_adapter = non_zero_adapters.loc[
        (non_zero_adapters[cols] > adapter_pct_threshold).any(1)
    ]

    if gt_pct_adapter.shape[0] > 0:
        res.append(
            f"> Warning: Found adapters that are present in more than {adapter_pct_threshold}% of all reads \n"
        )
        res.append(
            gt_pct_adapter[["Position", *cols_with_nonzero_sum]].to_string(index=False)
        )
        res.append("\n")
    return res


def check_sequence_length_dist(report: FastQCReport) -> List[str]:
    res = []

    if report.seq_length_dist.shape[0] > 1:
        res.append("> Warning: Sequences of multiple lengths detected \n")
        res.append(report.seq_length_dist.to_string(index=False))
        res.append("\n")

    zero_length_seq = report.seq_length_dist[report.seq_length_dist["Length"] == 0]
    if zero_length_seq.shape[0] > 0:
        res.append(
            f"> Warning: Found {zero_length_seq['Count'].iloc[0]} zero-length sequences."
        )

    return res


def analyze_fastqc_report(
    report: FastQCReport, validation_config: FastQCValidationConfig
) -> List[str]:
    """
    Returns a list of errors/warnings based on FastQC Results and Validation Config.

    References:
        - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8
    """

    logging_data: List[str] = []

    for attribute, step, arg in [
        (validation_config.check_poor_quality, check_poor_quality, 0),
        (validation_config.check_overall_gc_content, check_overall_gc_content, None),
        (
            validation_config.check_adapter_content,
            check_adapter_content,
            validation_config.adapter_pct_threshold,
        ),
        (
            validation_config.check_dedup_pct,
            check_dedup_pct,
            validation_config.dedup_pct_threshold,
        ),
        (
            validation_config.check_sequence_length_dist,
            check_sequence_length_dist,
            None,
        ),
    ]:
        if attribute:
            logging_data += step(report, arg) if arg else step(report)

    return logging_data

