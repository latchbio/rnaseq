"""latch/rnaseq"""

from itertools import groupby
from operator import attrgetter
import os
import subprocess
import types
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, Optional, Tuple, Union, Iterable
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
from latch import small_task, workflow, map_task
from latch.types import LatchDir, LatchFile
from latch.types.glob import file_glob

from wf.gtf_to_gbc import gtf_to_gbc


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


Replicate = Union[SingleEndReads, PairedEndReads]


@dataclass_json
@dataclass
class Sample:
    name: str
    strandedness: Strandedness
    replicates: List[Replicate]


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
) -> Tuple[LatchFile, LatchFile, LatchFile]:

    if type(genome) is LatchGenome:
        return LatchFile()

    # retrieves information needed for all other tasks
    # eg. get bed from bam file

    ...


@dataclass_json
@dataclass
class TrimgaloreInput:
    sample_name: str
    replicate: Replicate
    replicate_index: int
    run_name: str
    clip_r1: Optional[int]
    clip_r2: Optional[int] = None
    three_prime_clip_r1: Optional[int] = None
    three_prime_clip_r2: Optional[int] = None
    custom_output_dir: Optional[LatchDir] = None


@dataclass_json
@dataclass
class TrimgaloreOutput:
    sample_name: str
    trimmed_replicate: Replicate
    reports: List[LatchFile]


@dataclass_json
@dataclass
class TrimgaloreOutputGroup:
    sample_name: str
    outputs: List[TrimgaloreOutput]


@dataclass_json
@dataclass
class MergedSample:
    name: str
    reads: Replicate
    strandedness: Strandedness


@dataclass_json
@dataclass
class SalmonInput:
    merged_sample: MergedSample
    run_name: str
    ref: LatchGenome
    custom_ref_genome: Optional[LatchFile] = None
    custom_gtf: Optional[LatchFile] = None
    custom_ref_trans: Optional[LatchFile] = None
    custom_salmon_idx: Optional[LatchFile] = None
    custom_output_dir: Optional[LatchDir] = None


@small_task
def prepare_trimgalore_inputs(
    samples: List[Sample],
    run_name: str,
    clip_r1: Optional[int],
    clip_r2: Optional[int] = None,
    three_prime_clip_r1: Optional[int] = None,
    three_prime_clip_r2: Optional[int] = None,
    custom_output_dir: Union[None, LatchDir] = None,
) -> List[TrimgaloreInput]:
    _TMP_truncated_samples = [samples[0]]  # TODO
    return [
        TrimgaloreInput(
            sample_name=sample.name,
            replicate=replicate,
            replicate_index=i,
            run_name=run_name,
            clip_r1=clip_r1,
            clip_r2=clip_r2,
            three_prime_clip_r1=three_prime_clip_r1,
            three_prime_clip_r2=three_prime_clip_r2,
            custom_output_dir=custom_output_dir,
        )
        for sample in _TMP_truncated_samples
        for i, replicate in enumerate(sample.replicates)
    ]


@small_task
def prepare_salmon_inputs(
    merged_samples: List[MergedSample],
    run_name: str,
    ref: LatchGenome,
    custom_ref_genome: Optional[LatchFile] = None,
    custom_gtf: Optional[LatchFile] = None,
    custom_ref_trans: Optional[LatchFile] = None,
    custom_salmon_idx: Optional[LatchFile] = None,
    custom_output_dir: Optional[LatchDir] = None,
) -> List[SalmonInput]:
    return [
        SalmonInput(
            merged_sample=x,
            run_name=run_name,
            ref=ref,
            custom_ref_genome=custom_ref_genome,
            custom_gtf=custom_gtf,
            custom_ref_trans=custom_ref_trans,
            custom_salmon_idx=custom_salmon_index,
            custom_output_dir=custom_output_dir,
        )
        for x in merged_samples
    ]


def _merge_replicates(
    replicates: Iterable[Replicate],
    sample_name: str,
    custom_output_dir: Optional[LatchDir],
):
    base = _remote_output_dir(custom_output_dir)
    output_dir = f"{base}{sample_name}/"

    r1_path = os.path.join(output_dir, "r1_merged.fq")
    r1 = _concatenate_files((x.r1 for x in replicates), r1_path)

    if isinstance(replicates[0], SingleEndReads):
        return SingleEndReads(r1=r1)

    r2_path = os.path.join(output_dir, "r2_merged.fq")
    r2 = _concatenate_files((x.r2 for x in replicates), r2_path)
    return PairedEndReads(r1=r1, r2=r2)


def _concatenate_files(files: Iterable[LatchFile], output_file_name: str) -> LatchFile:
    with open(output_file_name, "w") as output_file:
        for file in files:
            with open(file.local_path, "r") as f:
                while chunk := f.read(10000):  # todo(rohankan): better chunk size?
                    output_file.write(chunk)
    return LatchFile(output_file_name, output_file_name)


@small_task
def group_trimgalore_outputs_by_sample_name(
    outputs: List[TrimgaloreOutput],
) -> List[TrimgaloreOutputGroup]:
    return [
        TrimgaloreOutputGroup(sample_name=sample_name, outputs=list(group))
        for sample_name, group in groupby(outputs, attrgetter("sample_name"))
    ]


@large_spot_task
def merge_per_sample_technical_replicates(
    samples: List[Sample],
    groups: List[TrimgaloreOutputGroup],
    custom_output_dir: Optional[LatchDir],
) -> List[MergedSample]:
    sample_name_to_strandedness = {x.name: x.strandedness for x in samples}
    return [
        MergedSample(
            name=x.sample_name,
            reads=_merge_replicates(
                replicates=(y.trimmed_replicate for y in x.outputs),
                sample_name=x.sample_name,
                custom_output_dir=custom_output_dir,
            ),
            strandedness=sample_name_to_strandedness[x.sample_name],
        )
        for x in groups
    ]


def _remote_output_dir(custom_output_dir: Optional[LatchDir]) -> str:
    if custom_output_dir is None:
        return "latch:///RNA-Seq Outputs/"
    remote_path = custom_output_dir.remote_path
    assert remote_path is not None
    if remote_path[-1] != "/":
        remote_path += "/"
    return remote_path


_SINGLE_TRIMGALORE_OUTPUT_FILENAME: str = "r1_trimmed.fq"
_R1_PAIRED_TRIMGALORE_OUTPUT_FILENAME: str = "r1_val_1.fq"
_R2_PAIRED_TRIMGALORE_OUTPUT_FILENAME: str = "r2_val_1.fq"


def _flag_if_not_none(name: str, input: TrimgaloreInput) -> List[str]:
    value = getattr(input, name)
    return [f"--{name}", value] if value is not None else []


@large_spot_task
def trimgalore(input: TrimgaloreInput) -> TrimgaloreOutput:
    reads = input.replicate

    flags = [
        *_flag_if_not_none("clip_r1", input),
        *_flag_if_not_none("three_prime_clip_r1", input),
    ]
    if isinstance(reads, SingleEndReads):
        run(["trim_galore --cores 8", *flags, str(reads.r1.local_path)])
    else:
        flags.extend(_flag_if_not_none("clip_r2", input))
        flags.extend(_flag_if_not_none("three_prime_clip_r2", input))
        run(
            [
                "trim_galore --cores 8 --paired",
                *flags,
                str(reads.r1.local_path),
                str(reads.r2.local_path),
            ]
        )

    base = f"{input.run_name}/Quality Control Data"
    end = f"{input.sample_name}/replicate_{input.replicate_index}/"
    remote_dir = _remote_output_dir(input.custom_output_dir)
    report_dir = remote_dir + f"{base}/Trimming Reports (TrimGalore)/{end}"
    read_dir = remote_dir + f"{base}/Trimming Reads (TrimGalore)/{end}"

    if isinstance(reads, SingleEndReads):
        r1 = LatchFile(
            _SINGLE_TRIMGALORE_OUTPUT_FILENAME,
            os.path.join(read_dir, _SINGLE_TRIMGALORE_OUTPUT_FILENAME),
        )
        trimmed_replicate = SingleEndReads(r1=r1)
    else:
        r1 = LatchFile(
            _R1_PAIRED_TRIMGALORE_OUTPUT_FILENAME,
            os.path.join(read_dir, _R1_PAIRED_TRIMGALORE_OUTPUT_FILENAME),
        )
        r2 = LatchFile(
            _R2_PAIRED_TRIMGALORE_OUTPUT_FILENAME,
            os.path.join(read_dir, _R2_PAIRED_TRIMGALORE_OUTPUT_FILENAME),
        )
        trimmed_replicate = PairedEndReads(r1=r1, r2=r2)

    # Side effect of trimgalore task - will not be used in future tasks
    reports = file_glob("*trimming_report.txt", report_dir)

    return TrimgaloreOutput(
        sample_name=input.sample_name,
        trimmed_replicate=trimmed_replicate,
        reports=reports,
    )


class InsufficientCustomGenomeResources(Exception):
    pass


class MalformedSalmonIndex(Exception):
    pass


@dataclass_json
@dataclass
class SalmonOutput:
    sf_files: List[LatchFile]
    aux_dir: LatchDir


@large_spot_task
def selective_alignment_salmon(input: SalmonInput) -> SalmonOutput:
    gm = lgenome.GenomeManager(input.ref.name)

    # lgenome
    # decoy = ref genome + transcript
    # genome + gtf = transcript

    custom_salmon_idx = input.custom_salmon_idx
    custom_ref_genome = input.custom_ref_genome
    custom_ref_trans = input.custom_ref_trans
    custom_gtf = input.custom_gtf
    sample = input.merged_sample
    run_name = input.run_name
    custom_output_dir = input.custom_output_dir

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

    reads = merged_sample.reads
    if isinstance(reads, SingleEndReads):
        reads = ["-r", str(reads.r1.local_path)]
    else:
        reads = ["-1", str(reads.r1.local_path), "-2", str(reads.r2.local_path)]

    for i, read in enumerate(reads):
        if read.endswith(".gz"):
            run(["gunzip", read])
            reads[i] = read.removesuffix(".gz")

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

    sf_files.append(LatchFile("/root/salmon_quant/quant.sf", output_literal))

    if custom_gtf is not None:
        gtf_path = custom_gtf.local_path
    else:
        gtf_path = gm.download_gtf()

    try:
        subprocess.run(
            [
                "/root/wf/run_tximport.R",
                "--args",
                "/root/salmon_quant/quant.sf",
                gtf_path,
                "/root/salmon_quant/genome_abundance.sf",
            ],
            check=True,
        )

        path_tail = f"{run_name}/Quantification (salmon)/{sample.name}/{sample.name}_genome_abundance.sf"
        if custom_output_dir is None:
            output_literal = "latch:///RNA-Seq Outputs/" + path_tail
        else:
            remote_path = custom_output_dir.remote_path
            if remote_path[-1] != "/":
                remote_path += "/"
            output_literal = remote_path + path_tail

        sf_files.append(
            LatchFile("/root/salmon_quant/genome_abundance.sf", output_literal)
        )
    except subprocess.CalledProcessError as e:
        print(
            f"Unable to produce gene mapping from tximport. Error surfaced from tximport -> {e}"
        )

    path_tail = f"{run_name}/Quantification (salmon)/{sample.name}/Auxilliary Info"
    if custom_output_dir is None:
        output_literal = "latch:///RNA-Seq Outputs/" + path_tail
    else:
        remote_path = custom_output_dir.remote_path
        if remote_path[-1] != "/":
            remote_path += "/"
        output_literal = remote_path + path_tail

    aux_dir = LatchDir("/root/salmon_quant/aux_info", output_literal)

    return SalmonOutput(sf_files=sf_files, aux_dir=aux_dir)


@small_task
def merge_salmon_outputs(
    outputs: List[SalmonOutput],
) -> Tuple[List[LatchFile], LatchDir]:
    _TMP_truncated_output = [outputs[0]]  # todo(rohankan): handle more samples
    return [(x.sf_files, x.aux_dir) for x in _TMP_truncated_output][0]


# @small_task
# def multiqc_task(
#     samples: List[Sample],
#     run_name: str,
#     sf_files: List[LatchFile],
#     aux_files: List[LatchDir],
#     ref: LatchGenome,
#     custom_ref_genome: Optional[LatchFile] = None,
#     custom_gtf: Optional[LatchFile] = None,
#     custom_ref_trans: Optional[LatchFile] = None,
#     custom_salmon_idx: Optional[LatchFile] = None,
#     custom_output_dir: Optional[LatchDir] = None,
# ) -> List[LatchFile]:
#     paths = []
#     for sample in samples:
#         for repl in sample.replicates:
#             paths.append(Path(repl.r1))
#             if hasattr(repl, "r2"):
#                 paths.append(Path(repl.r2))

#     subprocess.run(["/root/FastQC/fastqc", *paths], check=True)

#     quant_paths = []

#     for sf_file in sf_files:
#         paths.append(Path(sf_file))
#         quant_paths.append(Path(sf_file))
#     for aux_dir in aux_files:
#         paths.append(Path(aux_dir))

#     subprocess.run(["multiqc", *paths], check=True)

#     if custom_gtf is not None:
#         gtf_to_gbc(Path(custom_gtf), quant_paths)

#     path_tail = f"{run_name}/"
#     if custom_output_dir is None:
#         output_literal = "latch:///RNA-Seq Outputs/" + path_tail
#     else:
#         remote_path = custom_output_dir.remote_path
#         if remote_path[-1] != "/":
#             remote_path += "/"
#         output_literal = remote_path + path_tail

#     return [
#         LatchFile(
#             "/root/multiqc_report.html",
#             output_literal + "multiqc_report.html",
#         ),
#         LatchFile(
#             "/root/gene_body_coverage.png",
#             output_literal + "gene_body_coverage.png",
#         ),
#     ]


@small_task
def tximport_task(sf_files):
    ...


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
) -> Tuple[List[LatchFile], LatchDir]:
    """Performs alignment & quantification on Bulk RNA-Sequencing reads.

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
        display_name: Bulk RNAseq
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
    trimgalore_inputs = prepare_trimgalore_inputs(
        samples=samples,
        clip_r1=None,
        clip_r2=None,
        three_prime_clip_r1=None,
        three_prime_clip_r2=None,
        run_name=run_name,
        custom_output_dir=custom_output_dir,
    )
    trimgalore_outputs = map_task(trimgalore)(input=trimgalore_inputs)
    groups = group_trimgalore_outputs_by_sample_name(
        outputs=trimgalore_outputs,
    )
    # todo(rohankan): consider making this a map task
    merged_samples = merge_per_sample_technical_replicates(
        samples=samples,
        groups=groups,
        custom_output_dir=custom_output_dir,
    )
    salmon_inputs = prepare_salmon_inputs(
        merged_samples=merged_samples,
        ref=latch_genome,
        run_name=run_name,
        custom_gtf=custom_gtf,
        custom_ref_genome=custom_ref_genome,
        custom_ref_trans=custom_ref_trans,
        custom_salmon_idx=salmon_index,
        custom_output_dir=custom_output_dir,
    )
    salmon_outputs = map_task(selective_alignment_salmon)(input=salmon_inputs)
    return merge_salmon_outputs(outputs=salmon_outputs)
    # return multiqc_task(puts
    #     run_name=run_name,
    #     sf_files=quant_sfs,
    #     aux_files=aux_files,
    #     ref=latch_genome,
    #     custom_ref_genome=custom_ref_genome,
    #     custom_gtf=custom_gtf,
    #     custom_ref_trans=custom_ref_trans,
    #     custom_salmon_idx=salmon_index,
    #     custom_output_dir=custom_output_dir,
    # )


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
