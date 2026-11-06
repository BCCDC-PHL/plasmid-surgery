"""
Deterministic simulation of a fragmented plasmid reconstruction with known truth.

Reference: random sequence containing 3 copies of an IS-like repeat.
Sample: reference + SNPs, a 2 kb deletion, a 4 bp deletion, a 3 kb novel
insertion, and a rotated origin. The sample is broken into contigs at each
repeat copy (as happens with short-read assembly), with one collapsed repeat
contig, one extra random break, and alternating orientations. One contig spans
the reference origin.
"""
import random

COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def revcomp(seq: str) -> str:
    return seq[::-1].translate(COMPLEMENT)


def random_seq(rng: random.Random, n: int) -> str:
    return ''.join(rng.choice('ACGT') for _ in range(n))


def simulate(seed: int = 1, contig_overlap: int = 0) -> dict:
    """
    :param seed: RNG seed
    :param contig_overlap: Number of bp each contig is extended past its end into
                           the next contig (mimics k-1 overlaps from de Bruijn assemblers)
    :return: dict with 'ref' (str), 'sample' (str), 'contigs' (list of (id, seq)),
             'repeat' (str), 'truth_segments' (list of (id, sample-orientation seq))
    """
    rng = random.Random(seed)
    repeat = random_seq(rng, 1300)
    ref = ''.join([random_seq(rng, 15000), repeat, random_seq(rng, 15000), repeat,
                   random_seq(rng, 15000), repeat, random_seq(rng, 12000)])

    sample = list(ref)
    repeat_positions = set()
    # Avoid SNPs inside repeats so all copies stay identical
    repeat_spans = [(15000, 16300), (31300, 32600), (47600, 48900)]
    for start, end in repeat_spans:
        repeat_positions.update(range(start, end))
    candidates = [i for i in range(len(ref)) if i not in repeat_positions]
    for i in rng.sample(candidates, 300):
        sample[i] = rng.choice('ACGT'.replace(sample[i], ''))
    sample = ''.join(sample)
    sample = sample[:5000] + sample[7000:25000] + sample[25004:40000] + random_seq(rng, 3000) + sample[40000:]
    sample = sample[20000:] + sample[:20000]

    # Locate repeat copies in sample
    repeat_starts = []
    pos = sample.find(repeat)
    while pos >= 0:
        repeat_starts.append(pos)
        pos = sample.find(repeat, pos + 1)
    assert len(repeat_starts) == 3

    # Unique segments between repeat copies (circularly)
    L = len(sample)
    doubled = sample + sample
    segments = []
    for k, start in enumerate(repeat_starts):
        seg_start = start + len(repeat)
        seg_end = repeat_starts[(k + 1) % len(repeat_starts)]
        if seg_end <= seg_start:
            seg_end += L
        segments.append((seg_start, seg_end))
    # Extra break in the longest segment
    segments.sort(key=lambda s: s[1] - s[0], reverse=True)
    big_start, big_end = segments.pop(0)
    segments += [(big_start, big_start + 9000), (big_start + 9000, big_end)]

    truth_segments = []
    for start, end in segments:
        truth_segments.append(doubled[start:min(end + contig_overlap, start + L)])
    truth_segments.append(repeat)

    contigs = []
    for i, seg in enumerate(truth_segments):
        seq = revcomp(seg) if i % 2 else seg
        contigs.append((f"contig_{i}", seq))

    return {
        'ref': ref,
        'sample': sample,
        'repeat': repeat,
        'contigs': contigs,
        'truth_segments': [(f"contig_{i}", s) for i, s in enumerate(truth_segments)],
    }


def write_fasta(path, records):
    with open(path, 'w') as f:
        for seq_id, seq in records:
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), 70):
                f.write(seq[i:i + 70] + "\n")
