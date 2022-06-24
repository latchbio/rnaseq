"""
Using this paper as a baseline: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8
"""

from pathlib import Path
from os import makedirs
from typing import List, Tuple, Union

from wf import PairedEndReads, SingleEndReads
from wf.validation.quality.fastqc import (
    FastQCValidationConfig,
    analyze_fastqc_report,
    run_fastqc,
)


def _fastqc_helper(
    path: str, validation_config: FastQCValidationConfig, output_dir: str
) -> Tuple[str, str]:

    filename = Path(path).name.split(".")[0]
    fastqc_output_dir = Path(output_dir) / filename
    makedirs(fastqc_output_dir, exist_ok=True)

    report = run_fastqc(path, fastqc_output_dir)

    results = analyze_fastqc_report(report, validation_config)

    log_results = "\n".join(
        [
            f"=====> Starting FastQC Validation Report for file {report.filename}\n\n",
            "\n".join(results),
            f"\n\n=====> Completed FastQC Validation Report for file {report.filename}",
        ]
    )

    return (log_results, Path(fastqc_output_dir) / f"{filename}_fastqc.html")


def qc_raw_replicate(
    replicate: Union[PairedEndReads, SingleEndReads], output_dir: str
) -> Tuple[List[str], List[str]]:
    """
    Raw reads are evaluated on:
        - FastQC:
            - GC Content
            - Duplicated Reads
            - Sequence Quality
            - Adapter Content
    """

    raw_validation_config = FastQCValidationConfig(
        check_adapter_content=True,
        check_dedup_pct=True,
        check_poor_quality=True,
        check_overall_gc_content=True,
    )

    logs = []
    report_paths = []

    if type(replicate) == PairedEndReads:
        # Each fastq file in a PairedEnd read is analyzed separately
        for path in [replicate.r1.local_path, replicate.r2.local_path]:
            log_results, report_path = _fastqc_helper(
                path, raw_validation_config, output_dir
            )
            logs.append(log_results)
            report_paths.append(report_path)

    else:
        log_results, report_path = _fastqc_helper(
            replicate.r1.local_path, raw_validation_config, output_dir
        )
        logs.append(log_results)
        report_paths.append(report_path)

    return logs, report_paths

