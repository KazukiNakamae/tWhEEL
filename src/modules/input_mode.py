"""Shared input-mode resolution for the tWhEEL pipeline."""

import os


EXISTING_PROTEIN_MODE = 'existing_protein'
TRANSDECODER_MODE = 'transdecoder'


def get_configured_protein_path(config):
    """Return the normalized existing-protein path from the configuration."""
    protein_path = (config.get('input') or {}).get('protein_sequences', '')
    if protein_path is None:
        return ''
    return os.fspath(protein_path).strip()


def determine_protein_input_mode(config):
    """Resolve the effective protein input mode from one shared contract."""
    protein_path = get_configured_protein_path(config)
    transdecoder_config = config.get('transdecoder') or {}
    skip_transdecoder = transdecoder_config.get(
        'skip_if_protein_provided',
        True,
    )

    if protein_path and os.path.exists(protein_path) and skip_transdecoder:
        return EXISTING_PROTEIN_MODE
    return TRANSDECODER_MODE
