# unit test

reads = SingleEndReads(r1=LatchFile("/root/r1.fastq"))
trimgalore(
    reads=reads,
    clip_r1=None,
    clip_r2=None,
    three_prime_clip_r1=None,
    three_prime_clip_r2=None,
)

reads = PairedEndReads(r1=LatchFile("/root/r1.fastq"))
trimgalore(
    reads=reads,
    clip_r1=None,
    clip_r2=None,
    three_prime_clip_r1=None,
    three_prime_clip_r2=None,
)
