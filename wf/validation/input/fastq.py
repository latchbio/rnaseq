from dataclasses import dataclass
from itertools import zip_longest
from typing import Iterator
import gzip

from dataclasses_json import dataclass_json


class LatchValidationError(Exception):
    pass


@dataclass_json
@dataclass
class FastQRecord:
    identifier: str
    sequence: str
    separator: str
    quality_score_seq: str

    def __post_init__(self):
        # Validation checks to see if a FastQ record is actually formatted correctly

        for condition, error in [
            (len(self.identifier) == 0, "Identifier cannot be empty."),
            (not self.identifier.startswith("@"), "Identifier should start with '@'."),
            (len(self.separator) == 0, "Separator should not be empty."),
            (not self.separator.startswith("+"), "Separator start with '+'."),
            (
                len(self.quality_score_seq) != len(self.sequence),
                "Quality score length should match sequence length",
            ),
        ]:
            if condition:
                record = "\n".join(
                    [
                        self.identifier,
                        self.sequence,
                        self.separator,
                        self.quality_score_seq,
                    ]
                )
                raise LatchValidationError(
                    f"Badly formatted fastq record found: {error} \n\n{record}"
                )


def fastq_iterator(fastq_filepath: str, n_records: int = 5) -> Iterator[FastQRecord]:
    is_gzipped = fastq_filepath.endswith(".gz")
    opener, arg = (gzip.open, "rt") if is_gzipped else (open, "r")

    with opener(fastq_filepath, arg) as f:
        # Create a generator for lines needed to get n_records
        head = (next(f).rstrip() for _ in range(n_records * 4))

        # Create a generator which groups and generates tuples of 4 lines at a time
        for record in zip_longest(*[head] * 4):
            yield FastQRecord(*record)


def verify_paired_records(r1: FastQRecord, r2: FastQRecord) -> bool:
    # TODO (r614): Add a non-stupid way to find pattern using regex and extract the id
    uid_r1, uid_r2 = (
        r1.identifier.split("/")[0].split()[0],
        r2.identifier.split("/")[0].split()[0],
    )
    return uid_r1 == uid_r2
