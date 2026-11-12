# Slicing after a burn: time and memory

An investigation of two differences between cb-compare's two Setonix runs
(Slurm jobs 52637376 and 52790194, controlledburn-rs `56ba259` and
`e1696fc`). Neither is a regression in this crate: the slice code and
the burn's output are the same, and both differences come from the
memory allocator and the kernel. No code changed.

## 1. The serial slice after an approx burn was slower

On Setonix the serial `slice_tiles` was 1.2-2.5x slower in the second
run in every approx-mode case: for example small-1e5 at 256-cell tiles
went from 31 to 67 ms, and rock-32m at 2048-cell tiles from 19 to 35 ms
(the latter alternating 35, 28, 35, 28, 35 ms over the repeats). Coverage
mode was unchanged or faster, the approx burn itself took the same time
(184 and 181 ms on rock-32m), and burnscars' slice in the same jobs
moved by under 1%.

What it is not:

- **Not the burn's output.** Both versions give the same runs in the
  same order (1.70 M on rock-32m, 2.72 M on small-1e5).
- **Not the slice code.** Slicing either version's output with either
  version's `slice_tiles`, warm, takes the same time (rock-32m 24-25 ms,
  small-1e5 85-86 ms on the dev VM). Alone and in a loop the new slice
  is 1-4% slower, not 2x. The serial path of `slice_tiles_par` at one
  thread is `slice_tiles`.
- **Not the CPU or userspace.** The dev VM is also an AMD EPYC Milan
  (Zen 3), and runs the same Ubuntu 24.04 and glibc 2.39 as the
  container.

What it is: **first-touch page faults on the slice's output.** The
slice writes every record into a fresh table, and grows thousands of
tables. Measured on the dev VM, burning then slicing once, as the
harness does:

| input (approx) | tile | 4 KB pages: time, faults | glibc huge pages: time, faults |
|---|---|---|---|
| small-1e5 | 256 | 82 ms, 14,900 | 48 ms, 32 |
| rock-32m | 2048 | 33 ms, 3,150 | 25 ms, 0 |

(Huge pages through `GLIBC_TUNABLES=glibc.malloc.hugetlb=1`, which has
malloc ask for transparent huge pages.) Old and new are the same in
both columns. An approx result is runs only, cheap to route, so the
faults are a large share of its slice; in coverage mode the edges hide
them. Setonix's two runs are consistent with the two columns, the
first with the right and the second with the left (inferred: the runs
did not record faults or huge-page use). Which one a run gets depends on what the burn leaves in
the heap and on the kernel's huge-page policy (SLES, as on Setonix,
defaults to `always`; this VM uses `madvise`), not on the slice.

Tried and dropped: counting each tile's records first and allocating
every table at its exact size. It touches about 30% fewer pages, but
the extra pass made the serial slice 13-25% slower.

## 2. `slice_tiles_par` raised peak memory after a parallel burn

At 128 threads on Setonix, the harness's burn + slice peaked higher
with `slice_tiles_par` than with the serial slice before it: small-1e5
coverage at 256-cell tiles 0.54 to 1.04 GB, rock-32m at 2048-cell tiles
0.87 to 1.32 GB. The parallel burn alone did not change (median ratio
1.00), nor did the serial slice.

Reproduced on the dev VM with the harness at 32 threads (small-1e5
coverage, 256-cell tiles, peak RSS over three repeats):

| slice | default glibc | `MALLOC_ARENA_MAX=1` |
|---|---|---|
| serial | 677 MB | 571 MB |
| parallel | 1,076 MB | 618 MB |

With one malloc arena the parallel slice costs 47 MB more than the
serial one, about its band index (two bytes a record, plus counts),
which is the designed extra. The other ~400 MB is glibc's per-thread
arenas: each slicing thread grows its tiles' tables in its own arena,
and the memory the burn freed sits in other arenas, so little is reused
and every arena keeps its own slack. rock-32m showed the same (about 720 MB serial against
920-970 MB parallel).

Tried and dropped: exact-size tables on the parallel path only (the
counting pass runs in parallel there). It removed a third of the excess
at 32 threads, but at 4 and 16 threads it ranged from 15% less memory
to 2% more, and from 15% faster to 27% slower.

## What a caller can do

- Memory: cap glibc's arenas (`MALLOC_ARENA_MAX=2` to `4`; at 4 the
  parallel slice above peaked at 664 MB, below the default serial
  slice), or use another global allocator. The slice's own extra memory
  is one band index entry per record.
- Time: let malloc use huge pages (`GLIBC_TUNABLES=glibc.malloc.hugetlb=1`,
  or THP `always`), or slice into memory a previous slice freed.

## For the next benchmark run

The harness should record what decided these: the THP settings
(`/sys/kernel/mm/transparent_hugepage/{enabled,defrag}`), the glibc
version and malloc tunables, and the minor page faults of each timed
phase (`getrusage`), so a difference in the first-touch cost shows up as
faults rather than as a slower slice.
