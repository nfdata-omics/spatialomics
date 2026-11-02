# Reference resolution

## Objective and scope

Reference inputs used to prepare a Space Ranger index are resolved consistently
from explicit pipeline parameters or from the selected genome dictionary.

## Agreed behaviour

- `fasta`, `gtf`, `gff`, and `spaceranger_index` are declared optional
  parameters with `null` defaults.
- An explicit parameter has precedence over the corresponding entry in
  `params.genomes[params.genome]`.
- The resolved `reference_spec` map remains local to
  `NFDATAOMICS_SPATIALOMICS`. It is used for FASTQ validation and passed
  explicitly to `PREPARE_REF` together with the FASTQ trigger.
- Reference channels are enabled by the first validated FASTQ sample. Reference
  preparation is therefore skipped for runs containing only pre-computed Space
  Ranger outputs.
- `PREPARE_REF` converts the selected paths into channels and uses explicit
  conditional branches for a pre-built index, GTF/GFF selection, decompression,
  and index construction. It does not depend directly on pipeline parameters.
- `SPATIALOMICS` receives only the prepared Space Ranger index and the reference
  preparation version channel; it does not handle FASTA, GTF, GFF, archives, or
  reference names.
- A GFF or GFF3 annotation, optionally gzip-compressed, is accepted as an
  alternative to GTF and is converted to GTF before Space Ranger reference
  construction.

## Status

- Implementation: complete.
- Verification: not run after the latest workflow reorganization, per user
  request.
- Functional validation: pending a FASTQ run using direct reference parameters
  or a genome dictionary entry.
