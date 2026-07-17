"""Canonical release and CLI metadata for tWhEEL."""

from dataclasses import dataclass


# Update this value once when preparing a release. Runtime version reporting,
# documentation, and release-tag checks must refer back to this definition.
__version__ = "1.3.1"


@dataclass(frozen=True)
class StepDefinition:
    """Describe one step accepted by the public CLI dispatcher."""

    name: str
    workflow: str
    summary: str


STEP_DEFINITIONS = (
    StepDefinition(
        "gene_prediction",
        "GO/Pathway",
        "Prepare transcript/protein inputs and run SeqKit statistics and BUSCO.",
    ),
    StepDefinition(
        "annotation",
        "GO",
        "Annotate proteins against the GO protein database and classify groups.",
    ),
    StepDefinition(
        "sgd_annotation",
        "Pathway",
        "Annotate proteins against the SGD database and classify groups.",
    ),
    StepDefinition(
        "guide_design",
        "GO/Pathway",
        "Run CrisFlash and publish a normalized gene/guide generation.",
    ),
    StepDefinition(
        "go_analysis",
        "GO",
        "Add UniProt GO terms and calculate GO-KFI and GO Slim summaries.",
    ),
    StepDefinition(
        "pathway_analysis",
        "Pathway",
        "Match GPML genes and calculate PW-KFI with identifier audits.",
    ),
    StepDefinition(
        "report_generation",
        "GO",
        "Generate the standard GO report, summary, and report data.",
    ),
    StepDefinition(
        "pathway_report_generation",
        "Pathway",
        "Generate the pathway report and pathway summary data.",
    ),
)

STEP_NAMES = tuple(step.name for step in STEP_DEFINITIONS)


def format_step_listing(indent="  "):
    """Return an aligned, human-readable list of implemented steps."""
    name_width = max(len(step.name) for step in STEP_DEFINITIONS)
    return "\n".join(
        f"{indent}{step.name:<{name_width}}  "
        f"[{step.workflow}] {step.summary}"
        for step in STEP_DEFINITIONS
    )
