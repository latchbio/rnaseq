"""latch/rnaseq"""

import gzip
import os
import subprocess
import types
from dataclasses import dataclass
from enum import Enum
from itertools import zip_longest
from pathlib import Path
from typing import List, Optional, Tuple, Union
from urllib.parse import urlparse

import lgenome
from dataclasses_json import dataclass_json
from flytekit import LaunchPlan, task
from flytekitplugins.pod import Pod
from kubernetes.client.models import (
    V1Container,
    V1PodSpec,
    V1ResourceRequirements,
    V1Toleration,
)
from latch import small_task, workflow
from latch.types import LatchDir, LatchFile
from latch.types.glob import file_glob

from wf.validation.input.fastq import (
    LatchValidationError,
    fastq_iterator,
    verify_paired_records,
)
from wf.validation.quality.fastqc import FastQCValidationConfig, fastqc_helper


def run(cmd: List[str]):
    subprocess.run(cmd, check=True)


def _find_locals_in_set(param_set: set) -> List[str]:
    flags = []
    for param, val in locals().items():
        if param not in param_set or val is None:
            continue
        flags.extend(f"--{param}", str(val))
    return flags


# TODO - patch latch with proper def __repr__ -> str
def ___repr__(self):
    return str(self.local_path)


LatchFile.__repr__ = types.MethodType(___repr__, LatchFile)


def _is_valid_url(raw_url: str) -> bool:
    """A valid URL (as a source or destination of a LatchFile) must:
    * contain a latch or s3 scheme
    * contain an absolute path
    """
    try:
        parsed = urlparse(raw_url)
    except ValueError:
        return False
    if parsed.scheme not in ("latch", "s3"):
        return False
    if not parsed.path.startswith("/"):
        return False
    return True


def file_glob(
    pattern: str, remote_directory: str, target_dir: Optional[Path] = None
) -> List[LatchFile]:
    """Constructs a list of LatchFiles from a glob pattern.
    Convenient utility for passing collections of files between tasks. See
    [nextflow's channels](https://www.nextflow.io/docs/latest/channel.html) or
    [snakemake's wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards).
    for similar functionality in other orchestration tools.
    The remote location of each constructed LatchFile will be consructed by
    appending the file name returned by the pattern to the directory
    represented by the `remote_directory`.
    Args:
        pattern: A glob pattern to match a set of files, eg. '*.py'. Will
            resolve paths with respect to the working directory of the caller.
        remote_directory: A valid latch URL pointing to a directory, eg.
            latch:///foo. This _must_ be a directory and not a file.
        target_dir: An optional Path object to define an alternate working
            directory for path resolution
    Returns:
        A list of instantiated LatchFile objects.
    Intended Use: ::
        @small_task
        def task():
            ...
            return file_glob("*.fastq.gz", "latch:///fastqc_outputs")
    """

    if not _is_valid_url(remote_directory):
        return []

    if target_dir is None:
        wd = Path.cwd()
    else:
        wd = target_dir
    matched = sorted(wd.glob(pattern))

    return [LatchFile(str(file), remote_directory + file.name) for file in matched]


def _get_96_spot_pod() -> Pod:
    """[ "c6i.24xlarge", "c5.24xlarge", "c5.metal", "c5d.24xlarge", "c5d.metal" ]"""

    primary_container = V1Container(name="primary")
    resources = V1ResourceRequirements(
        requests={"cpu": "90", "memory": "170Gi"},
        limits={"cpu": "96", "memory": "192Gi"},
    )
    primary_container.resources = resources

    return Pod(
        pod_spec=V1PodSpec(
            containers=[primary_container],
            tolerations=[
                V1Toleration(effect="NoSchedule", key="ng", value="cpu-96-spot")
            ],
        ),
        primary_container_name="primary",
    )


large_spot_task = task(task_config=_get_96_spot_pod())


@dataclass_json
@dataclass
class SingleEndReads:
    r1: LatchFile


@dataclass_json
@dataclass
class PairedEndReads:
    r1: LatchFile
    r2: LatchFile


class ReadType(Enum):
    single = "single"
    paired = "paired"


class Strandedness(Enum):
    auto = "auto"


@dataclass_json
@dataclass
class Sample:
    name: str
    strandedness: Strandedness
    replicates: List[Union[SingleEndReads, PairedEndReads]]


class LatchGenome(Enum):
    RefSeq_hg38_p14 = "Homo sapiens (RefSeq hg38.p14)"
    RefSeq_T2T_CHM13v2_0 = "Homo sapiens (RefSeq T2T-CHM13v2.0)"
    RefSeq_R64 = "Saccharomyces cerevisiae (RefSeq R64)"
    RefSeq_GRCm39 = "Mus musculus (RefSeq GRCm39)"


# TODO - not used
@dataclass_json
@dataclass
class CustomGenome:
    gtf: LatchFile
    ref_genome: LatchFile
    ref_transcript: Optional[LatchFile]
    salmon_index: Optional[LatchFile]
    STAR_index: Optional[LatchFile]


@dataclass_json
@dataclass
class GenomeData:
    gtf: LatchFile


@small_task
def parse_inputs(
    genome: Union[LatchGenome, CustomGenome],
) -> (LatchFile, LatchFile, LatchFile):

    if type(genome) is LatchGenome:
        return LatchFile()

    # retrieves information needed for all other tasks
    # eg. get bed from bam file

    ...


def validate_fastq_reads(
    replicates: List[Union[PairedEndReads, SingleEndReads]]
) -> bool:
    for replicate in replicates:
        try:
            if type(replicate) == SingleEndReads:
                # verify if file exists, is readable and formatted correctly
                it = fastq_iterator(replicate.r1.local_path)
                for _ in it:
                    # iterate over the records to validate
                    continue

            elif type(replicate) == PairedEndReads:
                # verify if the paired end reads are in correct order + match
                r1_records = fastq_iterator(replicate.r1.local_path)
                r2_records = fastq_iterator(replicate.r2.local_path)

                for (record_1, record_2) in zip_longest(r1_records, r2_records):
                    if (record_1 is None and record_2 is not None) or (
                        record_2 is None and record_1 is not None
                    ):
                        longer_file = (
                            replicate.r1 if record_1 is not None else replicate.r2
                        )
                        raise LatchValidationError(
                            f"Error validating paired end reads - {longer_file.local_path} has more reads than the other."
                        )

                    if not verify_paired_records(record_1, record_2):
                        raise LatchValidationError(
                            f"Error validating paired end reads - order mismatch, trying to pair {record_1.identifier} to {record_2.identifier}."
                        )

        except gzip.BadGzipFile as e:
            raise LatchValidationError(
                f"Error opening gzipped file: {e.filename}. Are you sure the gzipped file isn't corrupted?"
            )

        except IOError as e:
            raise LatchValidationError(
                f"Error opening FastQ files: {replicate.r1.local_path}. \nAre you sure if the file exists and the current script has permissions to access it?"
            )

    return True


def qc_raw_replicate(
    replicate: Union[PairedEndReads, SingleEndReads], output_dir: str
) -> Tuple[List[str], List[Path]]:
    """
    Raw reads are evaluated on:
        - FastQC:
            - GC Content
            - Duplicated Reads
            - Sequence Quality
            - Adapter Content
            - Sequence Length != 0
            - num(Sequence Length) == 1
    """

    raw_validation_config = FastQCValidationConfig(
        check_adapter_content=True,
        check_dedup_pct=True,
        check_poor_quality=True,
        check_overall_gc_content=True,
        check_sequence_length_dist=True,
    )

    logs = []
    report_paths = []

    if type(replicate) == PairedEndReads:
        # Each fastq file in a PairedEnd read is analyzed separately
        for path in [replicate.r1.local_path, replicate.r2.local_path]:
            log_results, report_path = fastqc_helper(
                path, raw_validation_config, output_dir
            )
            logs.append(log_results)
            report_paths.append(report_path)

    else:
        log_results, report_path = fastqc_helper(
            replicate.r1.local_path, raw_validation_config, output_dir
        )
        logs.append(log_results)
        report_paths.append(report_path)

    return logs, report_paths


@small_task
def validation(samples: List[Sample], run_name: str) -> List[LatchFile]:
    """Run fastq input + qc validation for all samples."""

    reports = []

    for sample in samples:
        validate_fastq_reads(sample.replicates)
        sample_output_dir = Path(f"fastqc_outputs/{sample.name}").resolve()
        sample_output_dir.mkdir(parents=True, exist_ok=True)

        remote_sample_output_dir = f"latch:///RNA-Seq Outputs/{run_name}/Quality Control Data/FastQC Reports/{sample.name}"

        for replicate in sample.replicates:
            logs, report_paths = qc_raw_replicate(replicate, sample_output_dir)

            # TODO: Do something with the logs
            print("\n".join(logs))

            for report_path in report_paths:
                reports.append(
                    LatchFile(
                        str(report_path),
                        f"{remote_sample_output_dir}/{report_path.name}",
                    )
                )

    return reports


@large_spot_task
def trimgalore(
    samples: List[Sample],
    run_name: str,
    clip_r1: Optional[int] = None,
    clip_r2: Optional[int] = None,
    three_prime_clip_r1: Optional[int] = None,
    three_prime_clip_r2: Optional[int] = None,
    custom_output_dir: Union[None, LatchDir] = None,
) -> (List[Sample], List[LatchFile]):

    sample = samples[0]  # TODO

    trimmed_sample = Sample(
        name=sample.name, strandedness=sample.strandedness, replicates=[]
    )
    trimmed_replicates = []
    trimmed_reports = []
    for reads in sample.replicates:

        single_end_set = ("clip_r1", "three_prime_clip_r1")
        if type(reads) is SingleEndReads:
            flags = _find_locals_in_set(single_end_set)
            run(
                [
                    "trim_galore",
                    "--cores",
                    str(8),
                    *flags,
                    str(reads.r1.local_path),
                ]
            )
        else:
            paired_end_set = single_end_set + ("clip_r2", "three_prime_clip_r2")
            flags = _find_locals_in_set(paired_end_set)
            run(
                [
                    "trim_galore",
                    "--cores",
                    str(8),
                    "--paired",
                    *flags,
                    str(reads.r1.local_path),
                    str(reads.r2.local_path),
                ]
            )

        # Return trimming reports as a side effect.
        report_tail = f"{run_name}/Quality Control Data/Trimming Reports (TrimeGalore)/{sample.name}/"
        read_tail = f"{run_name}/Quality Control Data/Trimming Reads (TrimeGalore)/{sample.name}/"
        if custom_output_dir is None:
            report_literal = "latch:///RNA-Seq Outputs/" + report_tail
            read_literal = "latch:///RNA-Seq Outputs/" + read_tail
        else:
            remote_path = custom_output_dir.remote_path
            if remote_path[-1] != "/":
                remote_path += "/"
            report_literal = remote_path + report_tail
            read_literal = remote_path + read_tail

        trimmed_reports = file_glob("*trimming_report.txt", report_literal)
        trimmed = file_glob("*fq*", read_literal)

        if type(reads) is SingleEndReads:
            trimmed_replicates.append(SingleEndReads(r1=trimmed[0]))
        else:
            # glob results are sorted -  r1 will come first.
            trimmed_replicates.append(PairedEndReads(r1=trimmed[0], r2=trimmed[1]))

    trimmed_sample.replicates = trimmed_replicates

    return [trimmed_sample], trimmed_reports


class InsufficientCustomGenomeResources(Exception):
    pass


class MalformedSalmonIndex(Exception):
    pass


@large_spot_task
def sa_salmon(
    samples: List[Sample],
    run_name: str,
    ref: LatchGenome,
    custom_ref_genome: Optional[LatchFile] = None,
    custom_gtf: Optional[LatchFile] = None,
    custom_ref_trans: Optional[LatchFile] = None,
    custom_salmon_idx: Optional[LatchFile] = None,
    custom_output_dir: Optional[LatchDir] = None,
) -> List[LatchFile]:

    gm = lgenome.GenomeManager(ref.name)

    # lgenome
    # decoy = ref genome + transcript
    # genome + gtf = transcript

    if custom_salmon_idx is not None:
        # TODO: validate provided index...
        run(
            [
                "tar",
                "-xzvf",
                custom_salmon_idx.local_path,
            ]
        )
        if Path("salmon_index").is_dir() is False:
            raise MalformedSalmonIndex(
                "The custom Salmon index provided must be a directory named 'salmon_index'"
            )

    elif custom_ref_genome is not None:

        def _build_gentrome(genome: Path, transcript: Path) -> Path:
            run(
                [
                    "/root/gentrome.sh",
                    str(genome),
                    str(transcript),
                ]
            )
            return Path("/root/gentrome.fa")

        def _build_transcript(genome: Path, gtf: Path) -> Path:
            run(
                [
                    "/root/RSEM-1.3.3/rsem-prepare-reference",
                    "--gtf",
                    str(gtf),
                    "--num-threads",
                    "96",
                    str(genome),
                    "genome",
                ]
            )
            return Path("/root/genome.transcripts.fa")

        def _build_index(gentrome: Path) -> Path:
            run(
                [
                    "salmon",
                    "index",
                    "-t",
                    str(gentrome),
                    "-i",
                    "salmon_index",
                    "--decoys",
                    "decoys.txt",  # Comes from gentrome.sh
                    "-k",
                    "31",
                    "--threads",
                    "96",
                ]
            )
            return Path("/root/salmon_index")

        # Several ways to recover from missing files:
        #   1. Build a gentrome from provided genome + transcriptome
        #   2. First build a transcriptome from a genome + gtf, then 1.
        if custom_ref_trans is not None:
            gentrome = _build_gentrome(
                custom_ref_genome.local_path, custom_ref_trans.local_path
            )
            local_index = _build_index(gentrome)
        elif custom_gtf is not None:
            local_ref_trans = _build_transcript(
                custom_ref_genome.local_path, custom_gtf.local_path
            )
            gentrome = _build_gentrome(custom_ref_genome.local_path, local_ref_trans)
            local_index = _build_index(gentrome)
        else:
            raise InsufficientCustomGenomeResources(
                "Both a custom reference genome + GTF file need to be provided."
            )
    else:
        local_index = gm.download_salmon_index()

    sf_files = []

    sample = samples[0]  # TODO
    for reads in sample.replicates:

        if type(reads) is SingleEndReads:
            reads = ["-r", str(reads.r1.local_path)]
        else:
            reads = ["-1", str(reads.r1.local_path), "-2", str(reads.r2.local_path)]

        for i, read in enumerate(reads):
            if read[-3:] == ".gz":
                run(["gunzip", read])
                reads[i] = read[:-3]

        run(
            [
                "salmon",
                "quant",
                "-i",
                str(local_index),
                "-l",
                "A",
                *reads,
                "--threads",
                str(96),
                "--validateMappings",
                "-o",
                "salmon_quant",
            ]
        )

        path_tail = (
            f"{run_name}/Quantification (salmon)/{sample.name}/{sample.name}_quant.sf"
        )
        if custom_output_dir is None:
            output_literal = "latch:///RNA-Seq Outputs/" + path_tail
        else:
            remote_path = custom_output_dir.remote_path
            if remote_path[-1] != "/":
                remote_path += "/"
            output_literal = remote_path + path_tail

        sf_files.append(
            LatchFile(
                "/root/salmon_quant/quant.sf",
                output_literal,
            )
        )

    return sf_files


class AlignmentTools(Enum):
    star_salmon = "Traditional Alignment + Quantification"
    salmon = "Selective Alignment + Quantification"


@workflow
def rnaseq(
    samples: List[Sample],
    alignment_quantification_tools: AlignmentTools,
    ta_ref_genome_fork: str,
    sa_ref_genome_fork: str,
    output_location_fork: str,
    run_name: str,
    latch_genome: LatchGenome,
    bams: List[List[LatchFile]],
    custom_gtf: Optional[LatchFile] = None,
    custom_ref_genome: Optional[LatchFile] = None,
    custom_ref_trans: Optional[LatchFile] = None,
    star_index: Optional[LatchFile] = None,
    salmon_index: Optional[LatchFile] = None,
    save_indices: bool = False,
    custom_output_dir: Optional[LatchDir] = None,
) -> List[LatchFile]:
    """Produces gene and transcript counts from Bulk RNA-Sequencing reads.

    Bulk RNA-Seq (Alignment & Quantification)
    ----

    This workflow will produce gene and transcript counts from bulk RNA-seq
    sample reads.

    This current iteration of the workflow has two steps:

    1. Automatic read trimming using [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
    2. Selective alignment and quantification of reads using [Salmon](https://salmon.readthedocs.io/en/latest/)

    ### More about Selective Alignment

    This workflow uses Salmon's selective alignment described in this
    [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8),
    which achieves greater accuracy than traditional alignment methods while
    using less computational resources.


    __metadata__:
        display_name: RNAseq
        author:
            name: LatchBio
            email:
            github:
        repository: github.com/latchbio/rnaseq
        license:
            id: MIT
        flow:
        - section: Samples
          flow:
            - text: >-
                  Sample files can be provided and their read type can be
                  inferred from their name or this information can be specified manually.
                  Sample strandedness is inferred automatically (learn more).

            - params:
                - samples
        - section: Alignment & Quantification
          flow:
            - text: >-
                  This workflow uses Salmon's selective alignment described in this
                  [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8),
                  which achieves greater accuracy than traditional alignment methods while
                  using less computational resources.

            - fork: alignment_quantification_tools
              flows:
                traditional:
                    display_name: Selective Alignment
                    flow:
                        - fork: ta_ref_genome_fork
                          flows:
                            database:
                                display_name: Select from Latch Genome Database
                                flow:
                                    - text: >-
                                        We have curated a set of reference
                                        genome data for ease and
                                        reproducibility. More information about
                                        these managed files can be found
                                        [here](https://github.com/latchbio/latch-genomes).
                                    - params:
                                        - latch_genome
                            custom:
                                display_name: Provide Custom Genome
                                _tmp_unwrap_optionals:
                                    - custom_gtf
                                    - custom_ref_genome
                                flow:
                                    - params:
                                        - custom_ref_genome
                                        - custom_gtf
                                    - spoiler: Optional Params
                                      flow:
                                        - text: >-
                                            These files will be generated from the
                                            GTF/Genome files if not provided.
                                        - params:
                                            - salmon_index
                                            - custom_ref_trans
        - section: Output Settings
          flow:
          - params:
              - run_name
          - fork: output_location_fork
            flows:
                default:
                    display_name: Default
                    flow:
                    - text:
                        Output will be at default location in the data
                        viewer - RNA-Seq Outputs/"Run Name"
                custom:
                    display_name: Specify Custom Path
                    _tmp_unwrap_optionals:
                        - custom_output_dir
                    flow:
                    - params:
                        - custom_output_dir
    Args:

        samples:
            Here you can organize your FastQ files by sample and add technical
            replicates for each sample.  Biological replicates should be
            organized as separate samples.

          __metadata__:
            display_name: Sample Sheet
            batch_table_column: true
            _tmp:
                custom_ingestion: auto

        alignment_quantification_tools:

          __metadata__:
            display_name: Alignment & Quantification Method

        latch_genome:
          Curated reference files for specific genome sources and builds.

          __metadata__:
            batch_table_column: true
            display_name: Genome Database Option

        sa_ref_genome_fork:
          Select a reference genome from our curated database or provide your own.

          __metadata__:
            display_name: Reference Genome Source

        ta_ref_genome_fork:
          Select a reference genome from our curated database or provide your own.

          __metadata__:
            display_name: Reference Genome Source

        custom_ref_genome:
          The reference genome you want to align you samples to.

          __metadata__:
            display_name: Reference Genome File
            appearance:
                detail: (.fasta, .fasta.gz, .fa, .fa.gz, .fna, .fna.gz)

        custom_gtf:
          The gene annonation file that corresponds to the reference genome
          provided.

          __metadata__:
            display_name: Annotation File
            appearance:
                detail: (.gtf)

        bams:
          foobar

          __metadata__:
            display_name: bams

        custom_ref_trans:
          If not provided the workflow will generate from the Annotation File
          and Reference Genome File.

          __metadata__:
            display_name: Reference Transcript File (optional)
            appearance:
                detail: (.fasta, .fasta.gz, .fa, .fa.gz, .fna, .fna.gz)

        star_index:
          You are able to provide a zipped prebuilt STAR alignment index for
          your genome. This will speed up run time as the index is generated if
          none is provided. In output settings you are able to save indices
          from a run to be used in future runs.

          __metadata__:
            display_name: Provide Prebuilt STAR Index

        salmon_index:
            You are able to provide a zipped prebuilt Salmon selective alignment
            index for your genome. This will speed up run time as the index is
            generated if none is provided.

          __metadata__:
            display_name: salmon Index
            appearance:
                detail: (.tar.gz is only accepted extension)

        save_indices:
            If you provided a custom genome you can output the alignment
            indexes generated from this run for use in future runs. This will
            speed up runtime since the workflow doesn't have to then regenerate
            the indexes.

          __metadata__:
            display_name: Save Generated Reference Indexes

        run_name:
          A name for this analysis run, this will be used to name outputs from
          this run.

          __metadata__:
            batch_table_column: true
            display_name: Run Name

        output_location_fork:

        custom_output_dir:
          You can provide a custom location where this run's analysis outputs
          will be located.

          __metadata__:
            display_name: Custom Output Location
    """

    validation_reports = validation(samples=samples, run_name=run_name)

    trimmed_samples, reports = trimgalore(
        samples=samples,
        clip_r1=None,
        clip_r2=None,
        three_prime_clip_r1=None,
        three_prime_clip_r2=None,
        run_name=run_name,
        custom_output_dir=custom_output_dir,
    )
    return sa_salmon(
        samples=trimmed_samples,
        ref=latch_genome,
        run_name=run_name,
        custom_gtf=custom_gtf,
        custom_ref_genome=custom_ref_genome,
        custom_ref_trans=custom_ref_trans,
        custom_salmon_idx=salmon_index,
        custom_output_dir=custom_output_dir,
    )


if __name__ == "wf":
    LaunchPlan.create(
        "wf.__init__.rnaseq.Test Data",
        rnaseq,
        default_inputs={
            "samples": [
                Sample(
                    name="test_sample",
                    strandedness=Strandedness.auto,
                    replicates=[
                        PairedEndReads(
                            r1=LatchFile(
                                "s3://latch-public/welcome/rnaseq/r1.fastq",
                            ),
                            r2=LatchFile(
                                "s3://latch-public/welcome/rnaseq/r2.fastq",
                            ),
                        )
                    ],
                )
            ],
            "run_name": "Test Run",
        },
    )
