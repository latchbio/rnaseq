"""
Assemble and sort some COVID reads...
"""

import subprocess
from pathlib import Path

from latch import small_task, workflow
from latch.types import LatchFile


@small_task
def assembly_task(
    f1: LatchFile,
    f2: LatchFile,
    f3: LatchFile,
):
    f1.local_path

    open(f2)

    import time

    time.sleep(1000000)


@workflow
def assemble_and_sort(
    f1: LatchFile,
    f2: LatchFile,
    f3: LatchFile,
    ):
    """Description...

    markdown header
    ----

    Write some documentation about your workflow in
    markdown here:

    > Regular markdown constructs work as expected.

    # Heading

    * content1
    * content2

    __metadata__:
        display_name: Assemble and Sort FastQ Files
        author:
            name:
            email:
            github:
        repository:
        license:
            id: MIT

    Args:

        f1:
          Paired-end read 1 file to be assembled.

          __metadata__:
            display_name: Read1

        f2:
          Paired-end read 2 file to be assembled.

          __metadata__:
            display_name: Read2

        f3:
          Paired-end read 2 file to be assembled.

          __metadata__:
            display_name: Read3
    """
    assembly_task(f1=f1, f2=f2, f3=f3)
