"""latch/rnaseq"""

from collections import defaultdict
import csv
import re
import shutil
import subprocess
import types
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, Optional, Tuple, Union, Iterable
from urllib.parse import urlparse
import os

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
from latch import small_task, workflow, map_task, message
from latch.types import LatchDir, LatchFile

from wf.gtf_to_gbc import gtf_to_gbc


def _capture_output(command: List[str]) -> Tuple[int, str]:
    captured_stdout = []

    with subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
        universal_newlines=True,
    ) as process:
        assert process.stdout is not None
        for line in process.stdout:
            print(line)
            captured_stdout.append(line)
        process.wait()
        returncode = process.returncode

    return returncode, "\n".join(captured_stdout)


def run(command: List[str], check: bool = True, capture_output: bool = False):
    return subprocess.run(command, check=check, capture_output=capture_output)


# TODO - patch latch with proper def __repr__ -> str
def ___repr__(self):
    return str(self.local_path)


LatchFile.__repr__ = types.MethodType(___repr__, LatchFile)


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


class AlignmentTools(Enum):
    star_salmon = "Traditional Alignment + Quantification"
    salmon = "Selective Alignment + Quantification"


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


@dataclass_json
@dataclass
class TrimgaloreSalmonInput:
    sample_name: str
    replicates: List[Replicate]
    run_name: str
    base_remote_output_dir: str
    latch_genome: LatchGenome
    bams: List[List[LatchFile]]
    custom_names: List[LatchFile]
    custom_files: List[LatchFile]
    clip_r1: Optional[int] = None
    clip_r2: Optional[int] = None
    three_prime_clip_r1: Optional[int] = None
    three_prime_clip_r2: Optional[int] = None
    # custom_gtf: Optional[LatchFile] = None
    # custom_ref_genome: Optional[LatchFile] = None
    # custom_ref_trans: Optional[LatchFile] = None
    # star_index: Optional[LatchFile] = None
    # custom_salmon_index: Optional[LatchFile] = None
    save_indices: bool = False


@dataclass_json
@dataclass
class TrimgaloreSalmonOutput:
    passed_salmon: bool
    passed_tximport: bool
    sample_name: str
    sf_files: List[LatchFile]
    auxiliary_directory: List[LatchDir]
    # trimmed_reads: List[Replicate]
    trimgalore_reports: List[LatchFile]


@dataclass_json
@dataclass
class MergedSample:
    name: str
    reads: Replicate
    strandedness: Strandedness


@dataclass_json
@dataclass
class SalmonInput:
    sample_name: str
    trimmed_replicates: List[Replicate]
    strandedness: Strandedness
    run_name: str
    ref: str
    base_remote_output_dir: str
    # todo(rohankan): replace with actual names once Flyte can handle Optional[LatchFile] outputs
    custom_names: List[str]
    custom_files: List[LatchFile]


@small_task
def prepare_trimgalore_salmon_inputs(
    samples: List[Sample],
    run_name: str,
    latch_genome: LatchGenome,
    bams: List[List[LatchFile]],
    save_indices: bool,
    clip_r1: Optional[int] = None,
    clip_r2: Optional[int] = None,
    three_prime_clip_r1: Optional[int] = None,
    three_prime_clip_r2: Optional[int] = None,
    custom_output_dir: Optional[LatchDir] = None,
    custom_gtf: Optional[LatchFile] = None,
    custom_ref_genome: Optional[LatchFile] = None,
    custom_ref_trans: Optional[LatchFile] = None,
    custom_salmon_index: Optional[LatchFile] = None,
) -> List[TrimgaloreSalmonInput]:
    custom_names = []
    custom_files = []
    if custom_ref_genome is not None:
        custom_names.append("genome")
        custom_files.append(custom_ref_genome)
    if custom_ref_trans is not None:
        custom_names.append("trans")
        custom_files.append(custom_ref_trans)
    if custom_salmon_index is not None:
        custom_names.append("index")
        custom_files.append(custom_salmon_index)
    if custom_gtf is not None:
        custom_names.append("gtf")
        custom_files.append(custom_gtf)

    return [
        TrimgaloreSalmonInput(
            sample_name=sample.name,
            replicates=sample.replicates,
            run_name=run_name,
            clip_r1=clip_r1,
            clip_r2=clip_r2,
            three_prime_clip_r1=three_prime_clip_r1,
            three_prime_clip_r2=three_prime_clip_r2,
            base_remote_output_dir=_remote_output_dir(custom_output_dir),
            latch_genome=latch_genome,
            bams=bams,
            custom_names=custom_names,
            custom_files=custom_files,
            save_indices=save_indices,
        )
        for sample in samples
    ]


def _merge_replicates(
    replicates: List[Replicate], sample_name: str
) -> Union[Tuple[Path], Tuple[Path, Path]]:
    local_r1_path = f"{sample_name}_r1_merged.fq"
    r1 = _concatenate_files((str(x.r1.path) for x in replicates), local_r1_path)

    if isinstance(replicates[0], SingleEndReads):
        return (r1,)

    assert all(isinstance(x, PairedEndReads) for x in replicates)
    local_r2_path = f"{sample_name}_r2_merged.fq"
    r2 = _concatenate_files((str(x.r2.path) for x in replicates), local_r2_path)
    return (r1, r2)


def _concatenate_files(filepaths: Iterable[str], output_path: str) -> Path:
    path = Path(output_path).resolve()
    with path.open("w") as output_file:
        for p in filepaths:
            p = p.removesuffix(".gz")
            with open(p, "r") as f:
                shutil.copyfileobj(f, output_file)
    return path


def _remote_output_dir(custom_output_dir: Optional[LatchDir]) -> str:
    if custom_output_dir is None:
        return "latch:///RNA-Seq Outputs/"
    remote_path = custom_output_dir.remote_path
    assert remote_path is not None
    if remote_path[-1] != "/":
        remote_path += "/"
    return remote_path


class TrimgaloreError(Exception):
    pass


def do_trimgalore(
    ts_input: TrimgaloreSalmonInput,
    replicate_index: int,
    reads: Replicate,
) -> Tuple[List[LatchFile], Replicate]:
    def _flag(name: str) -> List[str]:
        value = getattr(ts_input, name)
        return [f"--{name}", value] if value is not None else []

    flags = [*_flag("clip_r1"), *_flag("three_prime_clip_r1")]
    read_paths = [str(reads.r1.local_path)]
    if isinstance(reads, PairedEndReads):
        flags += ["--paired", *_flag("clip_r2"), *_flag("three_prime_clip_r2")]
        read_paths.append(str(reads.r2.local_path))

    local_output = f"{ts_input.sample_name}_replicate_{replicate_index}"
    returncode, stdout = _capture_output(
        [
            "trim_galore",
            "--cores",
            str(8),
            "--dont_gzip",
            "--output_dir",
            f"./{local_output}",
            *flags,
            *read_paths,
        ]
    )

    # todo(rohankan): examine trimgalore for useful warnings and add them here
    if returncode != 0:
        stdout = stdout.rstrip()
        stdout = stdout[stdout.rindex("\n") + 1 :]
        assert reads.r1.remote_path is not None
        path_name = reads.r1.remote_path.split("/")[-1]
        identifier = f"sample {ts_input.sample_name}, replicate {path_name}"
        message(
            "error",
            {
                "title": f"Trimgalore error for {identifier}",
                "body": stdout,
            },
        )
        raise TrimgaloreError(stdout)

    def _output_path(middle: str) -> str:
        base = f"{ts_input.base_remote_output_dir}{ts_input.run_name}"
        tail = f"{ts_input.sample_name}/replicate_{replicate_index}/"
        return f"{base}/Quality Control Data/Trimming {middle} (TrimGalore)/{tail}"

    reads_directory = _output_path("Reads")
    if isinstance(reads, SingleEndReads):
        (r1,) = file_glob(f"{local_output}/*trimmed.fq*", reads_directory)
        trimmed_replicate = SingleEndReads(r1=r1)
    else:
        # File glob sorts files alphanumerically
        r1, r2 = file_glob(f"{local_output}/*val*.fq*", reads_directory)
        trimmed_replicate = PairedEndReads(r1=r1, r2=r2)

    os.remove(reads.r1.local_path)
    if isinstance(reads, PairedEndReads):
        os.remove(reads.r2.local_path)

    reports_directory = _output_path("Reports")
    reports = file_glob("*trimming_report.txt", reports_directory)

    return reports, trimmed_replicate


@large_spot_task
def trimgalore_salmon(
    samples: List[Sample],
    run_name: str,
    latch_genome: LatchGenome,
    bams: List[List[LatchFile]],
    save_indices: bool,
    clip_r1: Optional[int] = None,
    clip_r2: Optional[int] = None,
    three_prime_clip_r1: Optional[int] = None,
    three_prime_clip_r2: Optional[int] = None,
    custom_output_dir: Optional[LatchDir] = None,
    custom_gtf: Optional[LatchFile] = None,
    custom_ref_genome: Optional[LatchFile] = None,
    custom_ref_trans: Optional[LatchFile] = None,
    custom_salmon_index: Optional[LatchFile] = None,
) -> List[TrimgaloreSalmonOutput]:
    custom_names = []
    custom_files = []
    if custom_ref_genome is not None:
        custom_names.append("genome")
        custom_files.append(custom_ref_genome)
    if custom_ref_trans is not None:
        custom_names.append("trans")
        custom_files.append(custom_ref_trans)
    if custom_salmon_index is not None:
        custom_names.append("index")
        custom_files.append(custom_salmon_index)
    if custom_gtf is not None:
        custom_names.append("gtf")
        custom_files.append(custom_gtf)

    inputs = [
        TrimgaloreSalmonInput(
            sample_name=sample.name,
            replicates=sample.replicates,
            run_name=run_name,
            clip_r1=clip_r1,
            clip_r2=clip_r2,
            three_prime_clip_r1=three_prime_clip_r1,
            three_prime_clip_r2=three_prime_clip_r2,
            base_remote_output_dir=_remote_output_dir(custom_output_dir),
            latch_genome=latch_genome,
            bams=bams,
            custom_names=custom_names,
            custom_files=custom_files,
            save_indices=save_indices,
        )
        for sample in samples
    ]

    gm = lgenome.GenomeManager(latch_genome.name)
    gtf_path = (
        custom_gtf.local_path
        if custom_gtf is not None
        else gm.download_gtf(show_progress=False)
    )

    outputs = []
    for i, inp in enumerate(inputs):
        outputs.append(do_trimgalore_salmon(inp, i, gtf_path))
        # for rep in inp.replicates:
        #     os.remove(rep.r1.local_path)
        #     if isinstance(rep, PairedEndReads):
        #         os.remove(rep.r2.local_path)
    return outputs


def do_trimgalore_salmon(
    input: TrimgaloreSalmonInput,
    _tmp_index: int,
    gtf_path: os.PathLike,
) -> TrimgaloreSalmonOutput:
    outputs = [do_trimgalore(input, i, x) for i, x in enumerate(input.replicates)]
    trimmed_replicates = [x[1] for x in outputs]
    trimgalore_reports = [y for x in outputs for y in x[0]]

    gm = lgenome.GenomeManager(input.latch_genome.name)

    # lgenome
    # decoy = ref genome + transcript
    # genome + gtf = transcript

    custom_salmon_index = None
    custom_ref_genome = None
    custom_ref_trans = None
    custom_gtf = None
    for name, file in zip(input.custom_names, input.custom_files):
        if name == "index":
            custom_salmon_index = file
        elif name == "trans":
            custom_ref_trans = file
        elif name == "genome":
            custom_ref_genome = file
        elif name == "gtf":
            custom_gtf = file

    if custom_salmon_index is not None:
        # TODO: validate provided index...
        run(["tar", "-xzvf", custom_salmon_index.local_path])
        if Path("salmon_index").is_dir() is False:
            raise MalformedSalmonIndex(
                "The custom Salmon index provided must be a directory named 'salmon_index'"
            )

    elif custom_ref_genome is not None:

        def _build_gentrome(genome: Path, transcript: Path) -> Path:
            run(["/root/gentrome.sh", str(genome), str(transcript)])
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
                custom_ref_genome.local_path,
                custom_ref_trans.local_path,
            )
            local_index = _build_index(gentrome)
        elif custom_gtf is not None:
            local_ref_trans = _build_transcript(
                custom_ref_genome.local_path,
                custom_gtf.local_path,
            )
            gentrome = _build_gentrome(custom_ref_genome.local_path, local_ref_trans)
            local_index = _build_index(gentrome)
        else:
            message(
                "error",
                {
                    "title": "Unable to build local index",
                    "body": "Both a custom reference genome + GTF file need to be provided",
                },
            )
            raise InsufficientCustomGenomeResources(
                "Both a custom reference genome + GTF file need to be provided."
            )
    else:
        if _tmp_index == 0:
            local_index = gm.download_salmon_index(show_progress=False)
        else:
            local_index = Path("salmon_index")

    sf_files = []

    merged = _merge_replicates(trimmed_replicates, input.sample_name)
    if len(merged) == 1:
        (r1,) = merged
        reads = ["-r", str(r1)]
    else:
        r1, r2 = merged
        reads = ["-1", str(r1), "-2", str(r2)]

    returncode, stdout = _capture_output(
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

    # todo(rohankan): Add info messages too? Decide which Salmon logs
    # are interesting or important enough to show on console
    identifier = f"sample {input.sample_name}"
    for alert_type, alert_message in re.findall(_SALMON_ALERT_PATTERN, stdout):
        data = {"title": f"Salmon {alert_type} for {identifier}", "body": alert_message}
        if alert_type == "warning":
            message("warning", data)
        else:
            message("error", data)

    if returncode != 0:
        raise RuntimeError("Error occurred while running Salmon")

    sf_files.append(
        LatchFile(
            "/root/salmon_quant/quant.sf",
            _salmon_output_path(input, f"{input.sample_name}_quant.sf"),
        )
    )

    try:
        subprocess.run(
            [
                "/root/wf/run_tximport.R",
                "--args",
                "/root/salmon_quant/quant.sf",
                gtf_path,
                "/root/salmon_quant/genome_abundance.sf",
            ],
            capture_output=True,
        )

        sf_files.append(
            LatchFile(
                "/root/salmon_quant/genome_abundance.sf",
                _salmon_output_path(input, f"{input.sample_name}_genome_abundance.sf"),
            )
        )
    except subprocess.CalledProcessError as e:
        message("error", {"title": f"tximport error for {identifier}", "body": str(e)})
        print(
            f"Unable to produce gene mapping from tximport. Error surfaced from tximport -> {e}"
        )
        raise RuntimeError(f"Tximport error: {e}")

    auxiliary_directory = LatchDir(
        "/root/salmon_quant/aux_info",
        _salmon_output_path(input, "Auxiliary Info"),
    )

    for rep in trimmed_replicates:
        os.remove(rep.r1.path)
        if isinstance(rep, PairedEndReads):
            os.remove(rep.r2.path)

    # In case users want access to all the output files
    # all_salmon_outputs_directory = LatchDir(
    #     "/root/salmon_quant",
    #     _salmon_output_path(input, "full"),
    # )

    return TrimgaloreSalmonOutput(
        passed_salmon=True,
        passed_tximport=True,
        sample_name=input.sample_name,
        sf_files=sf_files,
        auxiliary_directory=[auxiliary_directory],
        # trimmed_reads=trimmed_replicates,
        trimgalore_reports=trimgalore_reports,
    )


class InsufficientCustomGenomeResources(Exception):
    pass


class MalformedSalmonIndex(Exception):
    pass


class SalmonError(Exception):
    pass


# Each Salmon warning or error log starts with a timestamp surrounded in square
# brackets ('\d{4}' represents the first part of the timestamp - the year)
_SALMON_ALERT_PATTERN = re.compile(r"\[(warning|error)\] (.+?)(?:\[\d{4}|$)", re.DOTALL)


def _salmon_output_path(inp: SalmonInput, suffix: str) -> str:
    return f"{inp.base_remote_output_dir}{inp.run_name}/Quantification (salmon)/{inp.sample_name}/{suffix}"


_COUNT_TABLE_GENE_ID_COLUMN = "gene_id"


@small_task
def count_matrix_and_multiqc(
    run_name: str,
    outputs: List[TrimgaloreSalmonOutput],
    output_directory: Optional[LatchDir],
    latch_genome: LatchGenome,
    custom_gtf: Optional[LatchFile] = None,
) -> List[LatchFile]:
    output_files = []

    # Beginning creation of combined count matrix
    tximport_outputs = [x for x in outputs if x.passed_tximport]
    if len(tximport_outputs) == len(outputs):
        message(
            "info",
            {
                "title": "Generating count matrix from all samples",
                "body": "\n".join(f"- {x.sample_name}" for x in tximport_outputs),
            },
        )

        combined_counts = defaultdict(dict)
        for output in tximport_outputs:
            genome_abundance_file = next(
                x for x in output.sf_files if x.remote_path.endswith("dance.sf")
            )
            with open(genome_abundance_file.local_path, "r") as f:
                for row in csv.DictReader(f, dialect=csv.excel_tab):
                    gene_name = row["Name"]
                    combined_counts[gene_name][output.sample_name] = row["NumReads"]

        raw_count_table_path = Path("./counts.tsv").resolve()
        with raw_count_table_path.open("w") as file:
            sample_names = (x.sample_name for x in tximport_outputs)
            writer = csv.DictWriter(
                file,
                [_COUNT_TABLE_GENE_ID_COLUMN, *sample_names],
                delimiter="\t",
            )
            writer.writeheader()
            for gene_id, data in combined_counts.items():
                data[_COUNT_TABLE_GENE_ID_COLUMN] = gene_id
                writer.writerow(data)

        base = _remote_output_dir(output_directory) + run_name
        output_path = f"{base}/Quantification (salmon)/counts.tsv"
        count_matrix_file = LatchFile(str(raw_count_table_path), output_path)
        output_files.append(count_matrix_file)
    else:
        message(
            "warning",
            {
                "title": "Unable to create combined count matrix",
                "body": "Some samples failed in the earlier step",
            },
        )

    # Beginning creation of MultiQC report
    salmon_outputs = [x for x in outputs if x.passed_salmon]
    paths = [Path(x.auxiliary_directory[0].local_path) for x in salmon_outputs]
    for output in salmon_outputs:
        paths += [Path(x.local_path) for x in output.sf_files]

    # todo (rohankan): Perform FastQC analysis on the reads?

    try:
        subprocess.run(["multiqc", *paths], check=True)
        base = _remote_output_dir(output_directory) + run_name
        multiqc_report_file = LatchFile(
            "/root/multiqc_report.html",
            f"{base}/multiqc_report.html",
        )
        output_files.append(multiqc_report_file)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while generating MultiQC report -> {e}")
        message(
            "error",
            {
                "title": "Unable to generate MultiQC report",
                "body": "See logs for more information",
            },
        )

    # try:
    #     gtf_path = (
    #         custom_gtf.local_path
    #         if custom_gtf is not None
    #         else lgenome.GenomeManager(latch_genome).download_gtf(show_progress=False)
    #     )
    #     quant_paths = [x for x in paths if str(x).endswith("quant.sf")]
    #     gtf_to_gbc(Path(gtf_path).resolve(), quant_paths)
    #     gbc_file = LatchFile(
    #         "/root/gene_body_coverage.png",
    #         f"{base}/gene_body_coverage.png",
    #     )
    #     output_files.append(gbc_file)
    # except Exception as e:
    #     print(e)
    #     message(
    #         "error",
    #         {
    #             "title": "Unable to generate gene body coverage (GBC) plots",
    #             "body": str(e),
    #         },
    #     )

    return output_files


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
    """Perform alignment and quantification on Bulk RNA-Sequencing reads

    Bulk RNA-Seq (Alignment and Quantification)
    ----

    This workflow will produce gene and transcript counts from bulk RNA-seq
    sample reads.

    This current iteration of the workflow has three steps:

    1. Automatic read trimming using [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
    2. Selective alignment and quantification of reads using [Salmon](https://salmon.readthedocs.io/en/latest/)
    3. Combination of your samples' reads into a count matrix

    ### More about Selective Alignment

    This workflow uses Salmon's selective alignment described in this
    [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8),
    which achieves greater accuracy than traditional alignment methods while
    using less computational resources.


    __metadata__:
        display_name: Bulk RNA-seq
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
        - section: Alignment and Quantification
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
        - section: Output Location
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
            display_name: Alignment and Quantification Method

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
    # trimgalore_salmon_inputs = prepare_trimgalore_salmon_inputs(
    #     samples=samples,
    #     run_name=run_name,
    #     clip_r1=None,
    #     clip_r2=None,
    #     three_prime_clip_r1=None,
    #     three_prime_clip_r2=None,
    #     custom_output_dir=custom_output_dir,
    #     latch_genome=latch_genome,
    #     bams=bams,
    #     custom_gtf=custom_gtf,
    #     custom_ref_genome=custom_ref_genome,
    #     custom_ref_trans=custom_ref_trans,
    #     custom_salmon_index=salmon_index,
    #     save_indices=save_indices,
    # )
    outputs = trimgalore_salmon(
        # inputs=trimgalore_salmon_inputs,
        samples=samples,
        run_name=run_name,
        clip_r1=None,
        clip_r2=None,
        three_prime_clip_r1=None,
        three_prime_clip_r2=None,
        custom_output_dir=custom_output_dir,
        latch_genome=latch_genome,
        bams=bams,
        custom_gtf=custom_gtf,
        custom_ref_genome=custom_ref_genome,
        custom_ref_trans=custom_ref_trans,
        custom_salmon_index=salmon_index,
        save_indices=save_indices,
    )
    return count_matrix_and_multiqc(
        run_name=run_name,
        outputs=outputs,
        output_directory=custom_output_dir,
        latch_genome=latch_genome,
        custom_gtf=custom_gtf,
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
