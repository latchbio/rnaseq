from dataclasses import dataclass
from itertools import zip_longest
from typing import Iterator, Union, List
from tqdm import tqdm
import gzip

from wf import PairedEndReads, SingleEndReads


class LatchValidationError(Exception):
    pass


@dataclass
class FastQRecord:
    identifier: str
    sequence: str
    separator: str
    quality_score_seq: str

    def __post_init__(self):
        # Validation checks to see if a FastQ record is actually formatted correctly

        if any(
            [
                len(self.identifier) == 0,
                "@" != self.identifier[0],  # identifier starts with @
                len(self.separator) == 0,
                "+" != self.separator[0],  # separator starts with +
                len(self.quality_score_seq) != len(self.sequence),
            ]
        ):
            record = "\n".join(
                [self.identifier, self.sequence, self.separator, self.quality_score_seq]
            )
            raise LatchValidationError(
                f"Badly formatted fastq record found:\n\n{record}"
            )


def fastq_iterator(fastq_filepath: str, n_records: int = 5) -> Iterator[FastQRecord]:
    is_gzipped = fastq_filepath.endswith(".gz")
    opener = gzip.open if is_gzipped else open

    with opener(fastq_filepath, "r") as f:
        head = (next(f).rstrip() for _ in range(n_records * 4))
        for record in zip_longest(*[head] * 4):
            if is_gzipped:
                yield FastQRecord(*[x.decode("ascii") for x in record])
            else:
                yield FastQRecord(*record)


def verify_paired_records(r1: FastQRecord, r2: FastQRecord) -> bool:
    uid_r1, uid_r2 = r1.identifier.split("/")[0], r2.identifier.split("/")[0]
    return uid_r1 == uid_r2


def validate_fastq_reads(
    replicates: List[Union[PairedEndReads, SingleEndReads]]
) -> bool:
    for replicate in tqdm(replicates):
        if type(replicate) == SingleEndReads:
            # verify if file exists, is readable and formatted correctly
            try:
                it = fastq_iterator(replicate.r1.local_path)
                for _ in it:
                    # iterate over the records to validate
                    continue

            except gzip.BadGzipFile as e:
                raise LatchValidationError(
                    f"Error opening gzipped file: {e.filename}. Are you sure if the gzipped file isn't corrupted?"
                )

            except IOError as e:
                raise LatchValidationError(
                    f"Error opening FastQ file: {replicate.r1.local_path}. \nAre you sure if the file exists/the current script has permissions to access it?"
                )

        elif type(replicate) == PairedEndReads:
            # verify if the paired end reads are in correct order + match
            try:
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
                    f"Error opening gzipped file: {e.filename}. Are you sure if the gzipped file isn't corrupted?"
                )

            except IOError as e:
                raise LatchValidationError(
                    f"Could not open one of the FastQ files in '{replicate.r1.local_path}' or '{replicate.r2.local_path}' \nAre you sure if the files exist/the current script has permissions to access it?"
                )

    return True
