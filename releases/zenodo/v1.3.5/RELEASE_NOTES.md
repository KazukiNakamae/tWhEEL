# tWhEEL v1.3.5

Source/tag date: 2026-09-27. Commit:
`77d7a7c9aa14ef1d5fc8e709293ab86999226d8a`.

GO-KFI reports now display the exact F, S, M, P, N, and total counts from
analysis statistics. Earlier reports reconstructed S and M from percentages
rounded to one decimal place, which could disagree with the underlying tables.
Stacked-bar widths now use unrounded count ratios, preserving small nonzero
categories and a correct total width. Percentage labels retain one decimal
place. The row label is now `GO terms`.

Missing, invalid, or inconsistent exact counts produce an explicit chart error.
Older payloads that retain exact top-level counts remain supported. The tagged
release includes regression tests for eleven retained publication count sets,
small categories, empty data, and invalid or legacy payloads.

This release changes report rendering, not GO-KFI calculation, guide filtering,
or pathway analysis. Existing analyses require report regeneration only for this
fix, using the original configuration and complete prerequisites. Preserve the
original report bundle and record the original analysis version separately
from the v1.3.5 renderer version.

## Archive scope

The source ZIP is generated from the immutable commit above, after checking
the local v1.3.5 tag. Every archived file is checked against its Git blob.
The ZIP retains all 96 tracked files, including previously tracked bytecode and
test-report static assets. Untracked files, `.git`, `LOCAL_DATA`, and other
untracked analysis outputs are absent. Metadata in this preparation directory
was added after the tag and is supplied separately; it is not part of the
original v1.3.5 source tree.

The archive includes configuration templates and Docker build recipes, but
does not include container images, reference databases, or complete publication
benchmark inputs/results. External image tags, dependencies, and online
annotation services are not an immutable execution environment. Source
preservation does not establish end-to-end reproducibility or biological
validation. KFI measures guide-design feasibility, not observed editing
efficiency or proven loss of function.

## License

MIT, as specified in the tagged `LICENSE`. External tools and reference
resources retain their own licenses.
