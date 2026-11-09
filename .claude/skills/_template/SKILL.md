---
name: <skill-name>
description: <One sentence: what it does + trigger phrases an agent should match, e.g. "Use when ...".>
---
# <Title>

## Command
```bash
.venv/bin/fibernhat <subcommand> --bam a.bam ... -o outputs/<x>
```

## Outputs
| file | content |
|---|---|

## Notes / pitfalls
- m6A calls: use `fibernhat.io.bam` (`iter_reads`, `read_to_fiber`, `Thresholds`). Do not re-parse MM/ML.
- If a modkit engine exists, add a parity test.

## Provenance
<original script paths this was ported from>
