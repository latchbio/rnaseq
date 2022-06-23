from pathlib import Path
from wf.quality_validation import parse_fastqc_report

raw_text = Path("/Users/roshan/dev/rnaseq/wf/test/fastqc_data.txt").read_text()


def test_fastqc_parser():
    testfile_path = Path(__file__).parent / "fastqc_data.txt"

    report = parse_fastqc_report(testfile_path)

    print(report)

    assert report.gc_pct == 40.0
    assert report.total_sequences == 50000
    assert report.poor_quality_sequences == 0

