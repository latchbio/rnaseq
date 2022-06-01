"""latch/rnaseq"""

import os
import subprocess
import types
from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Union

from dataclasses_json import dataclass_json
from flytekit import task
from flytekitplugins.pod import Pod
from kubernetes.client.models import (V1Container, V1PodSpec,
                                      V1ResourceRequirements, V1Toleration)
from latch import small_task, workflow
from latch.types import LatchDir, LatchFile
from latch.types.utils import file_glob


def run(cmd: List[str]):
    subprocess.run(cmd, check=True)


def _find_locals_in_set(param_set: set) -> List[str]:
    flags = []
    for param, val in locals().items():
        if param not in param_set or val is None:
            continue
        flags.extend((f"--{param}"), str(val))
    return flags


# TODO - patch latch with proper def __repr__ -> str
def ___repr__(self):
    return str(self.local_path)


LatchFile.__repr__ = types.MethodType(___repr__, LatchFile)


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


# TODO - not used
@dataclass_json
@dataclass
class CustomGenome:
    gtf: LatchFile
    ref_genome: LatchFile
    ref_transcript: Optional[LatchFile]
    salmon_index: Optional[LatchFile]
    STAR_index: Optional[LatchFile]


@small_task
def parse_inputs(
    genome: Union[LatchGenome, CustomGenome],
) -> (LatchFile, LatchFile, LatchFile):

    if type(genome) is LatchGenome:
        return LatchFile()

    # retrieves information needed for all other tasks
    # eg. get bed from bam file

    ...


@small_task
def trimgalore(
    samples: List[Sample],
    clip_r1: Optional[int] = None,
    clip_r2: Optional[int] = None,
    three_prime_clip_r1: Optional[int] = None,
    three_prime_clip_r2: Optional[int] = None,
) -> List[List[LatchFile]]:

    sample = samples[0]  # TODO

    trimmed_reads = []
    for reads in sample.replicates:

        single_end_set = ("clip_r1", "three_prime_clip_r1")
        if type(reads) is SingleEndReads:
            flags = _find_locals_in_set(single_end_set)
            run(
                [
                    "trim_galore",
                    "--cores",
                    str(96),  # TODO
                    "--gzip",
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
                    str(96),  # TODO
                    "--paired",
                    "--gzip",
                    *flags,
                    str(reads.r1.local_path),
                    str(reads.r2.local_path),
                ]
            )

        # Return trimming reports as a side effect.
        trimmed_reads.append(
            file_glob(
                "*trimming_report.txt",
                f"latch:///Quality Control Data/Trimming Reports (TrimeGalore)/{sample.name}/",
            )
        )

    return trimmed_reads


@large_spot_task
def align_star(
    samples: List[Sample],
    ref: LatchGenome,
) -> List[List[LatchFile]]:

    # handle index

    os.mkdir("STAR_index")

    run(
        [
            "aws",
            "s3",
            "sync",
            "s3://latch-genomes/Homo_sapiens/RefSeq/GRCh38.p14/STAR_index",
            "STAR_index",
        ]
    )

    sample = samples[0]  # TODO

    sample_bams = []
    for reads in sample.replicates:

        if type(reads) is SingleEndReads:
            reads = [str(reads.r1.local_path)]
        else:
            reads = [str(reads.r1.local_path), str(reads.r2.local_path)]

        run(
            [
                "STAR",
                "--genomeDir",
                "STAR_index",
                "--readFilesIn",
                *reads,
                "--runThreadN",
                str(96),  # TODO
                "--outFileNamePrefix",
                sample.name,
                "--outSAMtype",
                "BAM",
                "Unsorted",
                "SortedByCoordinate",
            ]
        )

        sample_bams.append(
            [
                file_glob("*out.bam", f"latch:///Alignment (STAR)/{sample.name}")[0],
                file_glob(
                    "*sortedByCoord.out.bam", f"latch:///Alignment (STAR)/{sample.name}"
                )[0],
            ]
        )
    return sample_bams


@small_task
def quantify_salmon(
    bams: List[List[LatchFile]],
    ref: LatchGenome,
) -> List[LatchDir]:

    run(
        [
            "aws",
            "s3",
            "cp",
            "s3://latch-genomes/Homo_sapiens/RefSeq/GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.transcripts.fna",
            ".",
        ]
    )

    # pseudoaligner
    # ./bin/salmon quant - i transcripts_index - l < LIBTYPE > -1 reads1.fq - 2 reads2.fq - -validateMappings - o transcripts_quant
    # ./bin/salmon quant - t transcripts.fa - l < LIBTYPE > -a aln.bam - o salmon_quant
    # sudo ./salmon-1.5.2_linux_x86_64/bin/salmon quant - t GCF_000001405.40_GRCh38.p14_genomic.transcripts.fna - l A - a foobarAligned.out.bam - o salmon_quant

    quantified_bams = []
    for bam_set in bams:
        bam = bam_set[0]
        run(
            [
                "salmon",
                "quant",
                "-t",
                str("GCF_000001405.40_GRCh38.p14_genomic.transcripts.fna"),
                "-a",
                str(bam.local_path),
                "--threads",
                str(96),  # TODO
                "--libType=A",
                "-o",
                "salmon_quant",  # TODO:
            ]
        )
        quantified_bams.append(
            LatchDir("salmon_quant", "latch:///Quantification (salmon)/foobar/")
        )

    # tx2gene
    # 2gene
    # salmon_tx2gene.py*
    # salmon_tximport.r*
    # summarize_experiments -> .rds

    return quantified_bams


@small_task
def foo(samples: List[Sample]):
    ...


class AlignmentTools(Enum):
    star_salmon = "Traditional Alignment + Quantification"
    salmon = "Selective Alignment + Quantification"


@workflow
def rnaseq(
    samples: List[Sample],
    alignment_quantification_tools: AlignmentTools,
    ta_ref_genome_fork: str,
    output_location_fork: str,
    run_name: str,
    latch_genome: LatchGenome,
    bams: List[List[LatchFile]],
    gtf: Optional[LatchFile] = None,
    ref_genome: Optional[LatchFile] = None,
    ref_transcript: Optional[LatchFile] = None,
    star_index: Optional[LatchFile] = None,
    salmon_index: Optional[LatchFile] = None,
    save_indices: bool = False,
    custom_output_dir: Optional[LatchDir] = None,
):
    """rnaseq

    rnaseq
    ----

    Write some documentation about your workflow in
    markdown here:

    > Regular markdown constructs work as expected.

    # Heading

    * content1
    * content2

    __metadata__:
        display_name: RNAseq
        author:
            name: Kenny Workman
            email:
            github:
        repository: github.com/latchbio/rnaseq
        license:
            id: MIT
        flow:
        - section: Samples
          flow:
            - text:
                  Sample files can be provided and their read type can be
                  inferred from their name (learn more about name formatting used here)
                  or this information can be specified manually. Sample strandedness is
                  inferred automatically (learn more).
            - params:
                - bams
                - samples
        - section: Alignment & Quantification
          flow:
            - text:
                  Two methods are available for the alignment and quantification of your reads.
                  "Traditional alignment" is the more accurate but expensive (in terms of time and
                  computing resources) option. "Selective alignment" is a faster mapping algorithm
                  that is slightly less accurate. Often the differences between accuracy is
                  minimal between these two methods - read more
                  [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8).
            - fork: alignment_quantification_tools
              flows:
                selective:
                    display_name: Selective Alignment
                    flow:
                        - section: Reference Genome
                          flow:
                            - fork: sa_ref_genome_fork
                              flows:
                                database:
                                    display_name: Select from Latch Genome Database
                                    flow:
                                    - text:
                                        We have curated a set of reference
                                        genome data for ease and
                                        reproducibility. More information about
                                        these managed files can be found
                                        [here](https://github.com/latchbio/latch-genomes).
                                    - params:
                                        - latch_genome
                                custom:
                                    display_name: Provide Custom Genome Data
                                    flow:
                                    - text:
                                        When providing custom reference
                                        data to the selective alignment method,
                                        only a reference transcriptome is
                                        needed. This can be provided directly,
                                        generated from a genome + annotation
                                        file or provided pre-built as an index.
                                    - params:
                                        - gtf
                                        - ref_genome
                                    - spoiler: Optional Params
                                      flow:
                                        - text:
                                              These files will be generated from the
                                              GTF/Genome files if not provided.
                                        - params:
                                            - ref_transcript
                                            - salmon_index
                            - spoiler: Advanced Selective Alignment Settings
                              flow:
                                - params:
                                    - star_index
                        - spoiler: General Alignment Settings
                          flow:
                            - params:
                                - star_index
                traditional:
                    display_name: Traditional Alignment
                    flow:
                        - section: Reference Genome
                          flow:
                            - fork: ta_ref_genome_fork
                              flows:
                                database:
                                    display_name: Select from Latch Genome Database
                                    flow:
                                    - text:
                                        We have curated a set of reference
                                        genome data for ease and
                                        reproducibility. More information about
                                        these managed files can be found
                                        [here](https://github.com/latchbio/latch-genomes).
                                    - params:
                                        - latch_genome
                                custom:
                                    display_name: Provide Custom Genome Data
                                    flow:
                                    - params:
                                        - gtf
                                        - ref_genome
                                    - spoiler: Optional Params
                                      flow:
                                        - text:
                                              These files will be generated from the
                                              GTF/Genome files if not provided.
                                        - params:
                                            - ref_transcript
                                            - star_index
                            - spoiler: Advanced Traditional Alignment Settings
                              flow:
                                - params:
                                    - star_index
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
                        viewer - RNA-Seq A&Q Outputs > "Run Name"
                custom:
                    display_name: Specify Custom Path
                    flow:
                    - params:
                        - custom_output_dir
          - spoiler: Advanced Output Settings
            flow:
              - params:
                  - save_indices
    Args:

        bams:
          RNAseq data is generated as a collection of FastQ files. Here you can
          organize your FastQ files by sample and add technical replicates for
          each sample.

          __metadata__:
            display_name: BAMS

        samples:
          RNAseq data is generated as a collection of FastQ files. Here you can
          organize your FastQ files by sample and add technical replicates for
          each sample.

          __metadata__:
            display_name: Sample Sheet
            _tmp:
              custom_ingestion: auto

        alignment_quantification_tools:
          foobar

          __metadata__:
            display_name: Alignment & Quantification Tools

        latch_genome:
          Curated reference files for specific genome sources and builds.
          pull

          __metadata__:
            display_name: Genome Database Option

        ta_ref_genome_fork:
          foobar

          __metadata__:
            display_name: Ref Genome Source

        ref_genome:
          foobar

          __metadata__:
            display_name: Reference Genome File
            detail: (.fasta, .fasta.gz, .fa, .fa.gz, .fna, .fna.gz)

        gtf:
          foobar

          __metadata__:
            display_name: Annotation File

        ref_transcript:
          foobar

          __metadata__:
            display_name: Reference Transcript File

        star_index:
          foobar

          __metadata__:
            display_name: Provide Prebuilt STAR Index

        salmon_index:
          foobar

          __metadata__:
            display_name: salmon Index

        save_indices:
          foobar

          __metadata__:
            display_name: Save Generated Reference Indexes

        run_name:
          foobar

          __metadata__:
            display_name: Run Name

        output_location_fork:
          foobar

        custom_output_dir:
          foobar

          __metadata__:
            display_name: Custom Output Location
    """

    trimmed_samples = trimgalore(
        samples=samples,
        clip_r1=None,
        clip_r2=None,
        three_prime_clip_r1=None,
        three_prime_clip_r2=None,
    )
    bams = align_star(samples=samples, ref=latch_genome)
    quantify_salmon(bams=bams, ref=latch_genome)


if __name__ == "__main__":
    samples = [
        Sample(
            name="foo",
            strandedness=Strandedness.auto,
            replicates=[SingleEndReads(r1=LatchFile("/root/r1.fastq"))],
        )
    ]
    trimgalore(
        samples=samples,
        clip_r1=None,
        clip_r2=None,
        three_prime_clip_r1=None,
        three_prime_clip_r2=None,
    )
