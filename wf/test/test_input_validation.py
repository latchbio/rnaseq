from dataclasses import dataclass
from wf import Sample, PairedEndReads, Strandedness, SingleEndReads
from wf.input_validation import validate_fastq_reads

import pytest


@dataclass()
class MockLatchFile:
    local_path: str


@pytest.fixture
def paired_read_sample() -> Sample:
    return Sample(
        "Test Data",
        Strandedness.auto,
        [
            PairedEndReads(
                MockLatchFile(
                    local_path="/Users/roshan/Downloads/SRR6357070_1.fastq.gz"
                ),
                MockLatchFile(
                    local_path="/Users/roshan/Downloads/SRR6357070_2.fastq.gz"
                ),
            )
        ],
    )


@pytest.fixture
def single_read_sample() -> Sample:
    return Sample(
        "Test Data",
        Strandedness.auto,
        [
            SingleEndReads(
                MockLatchFile(
                    local_path="/Users/roshan/Downloads/SRR6357070_1.fastq.gz"
                )
            )
        ],
    )


def test_validate_paired_reads(paired_read_sample: Sample):
    assert validate_fastq_reads(paired_read_sample.replicates) == True


def test_validate_single_reads(single_read_sample: Sample):
    assert validate_fastq_reads(single_read_sample.replicates) == True
