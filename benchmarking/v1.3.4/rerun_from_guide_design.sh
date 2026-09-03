#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MANIFEST="${SCRIPT_DIR}/publication_runs.tsv"
DATA="${ROOT}/LOCAL_DATA/20260715_benchmarking"
OUTPUT_ROOT="${ROOT}/outputs/benchmarking"
PYTHON="${TWHEEL_V134_PYTHON:-${ROOT}/.handoff_runtime/twheel_env-osx-64/bin/python}"
RUN_LOG_DIR="${DATA}/run_logs/v1.3.4_publication"
STATUS_FILE="${DATA}/v1.3.4_publication_rerun_status.tsv"
LOCK_DIR="${DATA}/v1.3.4_publication_rerun.lock"
MODE="preflight"
SELECTED_RUN=""
SKIP_DOCKER=0

SEED_ITEMS=(
    annotation_classified.txt
    annotation_results.txt
    annotation_table.txt
    gene_prediction_results.txt
    query.pep
    transcriptome.fa
    multi_genes.fa
    single_genes.fa
    multi_genes_aligned
)

usage() {
    printf 'Usage: %s [--list | --preflight-only | --prepare | --execute] [--run RUN_KEY] [--skip-docker-check]\n' "$0"
}

while (( $# > 0 )); do
    case "$1" in
        --list)
            MODE="list"
            shift
            ;;
        --preflight-only)
            MODE="preflight"
            shift
            ;;
        --prepare)
            MODE="prepare"
            shift
            ;;
        --execute)
            MODE="execute"
            shift
            ;;
        --run)
            [[ $# -ge 2 && -n "$2" ]] || {
                printf 'Option --run requires a run key.\n' >&2
                exit 2
            }
            SELECTED_RUN="$2"
            shift 2
            ;;
        --skip-docker-check)
            SKIP_DOCKER=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            usage >&2
            exit 2
            ;;
    esac
done

if [[ "${MODE}" == "execute" && "${SKIP_DOCKER}" -eq 1 ]]; then
    printf 'Docker checks cannot be skipped in execute mode.\n' >&2
    exit 2
fi

[[ -f "${MANIFEST}" ]] || {
    printf 'Run manifest not found: %s\n' "${MANIFEST}" >&2
    exit 1
}

selected() {
    local run_key="$1"
    [[ -z "${SELECTED_RUN}" || "${run_key}" == "${SELECTED_RUN}" ]]
}

list_runs() {
    column -s $'\t' -t "${MANIFEST}"
}

item_digest() {
    local item_path="$1"
    if [[ -f "${item_path}" ]]; then
        shasum -a 256 "${item_path}" | awk '{print $1}'
        return
    fi
    (
        cd "${item_path}"
        find . -type f -print0 \
            | LC_ALL=C sort -z \
            | xargs -0 shasum -a 256
    ) | shasum -a 256 | awk '{print $1}'
}

marker_value() {
    local marker="$1"
    local field="$2"
    awk -F '\t' -v field="${field}" \
        '$1 == field { print $2; exit }' "${marker}"
}

validate_existing_seed() {
    local run_key="$1"
    local config_rel="$2"
    local source_run="$3"
    local target_run="$4"
    local source_intermediate="${OUTPUT_ROOT}/${source_run}/intermediate"
    local target_dir="${OUTPUT_ROOT}/${target_run}"
    local target_intermediate="${target_dir}/intermediate"
    local marker="${target_dir}/.v1.3.4-rerun-seed.tsv"
    local config_path="${ROOT}/${config_rel}"
    local item recorded_digest source_digest target_digest

    [[ "$(marker_value "${marker}" schema)" == "1" ]]
    [[ "$(marker_value "${marker}" run_key)" == "${run_key}" ]]
    [[ "$(marker_value "${marker}" tWhEEL_version)" == "1.3.4" ]]
    [[ "$(marker_value "${marker}" config)" == "${config_rel}" ]]
    [[ "$(marker_value "${marker}" config_sha256)" \
        == "$(shasum -a 256 "${config_path}" | awk '{print $1}')" ]]
    [[ "$(marker_value "${marker}" source_run)" == "${source_run}" ]]
    [[ "$(marker_value "${marker}" target_run)" == "${target_run}" ]]
    [[ "$(marker_value "${marker}" restart_step)" == "guide_design" ]]
    [[ "$(marker_value "${marker}" guide_design_sha256)" \
        == "$(shasum -a 256 "${ROOT}/src/modules/guide_design.py" \
            | awk '{print $1}')" ]]

    for item in "${SEED_ITEMS[@]}"; do
        [[ -s "${source_intermediate}/${item}" ]]
        [[ -s "${target_intermediate}/${item}" ]]
        recorded_digest="$(awk -F '\t' -v item="${item}" \
            '$1 == "artifact_sha256" && $2 == item { print $3; exit }' \
            "${marker}")"
        source_digest="$(item_digest "${source_intermediate}/${item}")"
        target_digest="$(item_digest "${target_intermediate}/${item}")"
        [[ -n "${recorded_digest}" ]]
        [[ "${recorded_digest}" == "${source_digest}" ]]
        [[ "${recorded_digest}" == "${target_digest}" ]]
    done
}

check_source_artifacts() {
    local source_run="$1"
    local source_intermediate="${OUTPUT_ROOT}/${source_run}/intermediate"
    local item
    [[ -d "${source_intermediate}" ]] || {
        printf 'Source intermediate directory is missing: %s\n' \
            "${source_intermediate}" >&2
        return 1
    }
    for item in "${SEED_ITEMS[@]}"; do
        [[ -s "${source_intermediate}/${item}" ]] || {
            printf 'Required source artifact is missing or empty: %s\n' \
                "${source_intermediate}/${item}" >&2
            return 1
        }
    done
}

check_common_preflight() {
    local version_output
    [[ -x "${PYTHON}" ]] || {
        printf 'Python runtime is not executable: %s\n' "${PYTHON}" >&2
        return 1
    }
    version_output="$("${PYTHON}" "${ROOT}/src/twheel.py" --version)"
    [[ "${version_output}" == "tWhEEL version 1.3.4" ]] || {
        printf 'Expected tWhEEL version 1.3.4, observed: %s\n' \
            "${version_output}" >&2
        return 1
    }
    command -v shasum >/dev/null
    command -v cp >/dev/null
    if (( SKIP_DOCKER == 0 )); then
        docker info >/dev/null
        docker image inspect kazukinakamae/crisflash:1.2.0 >/dev/null
        docker image inspect kazukinakamae/add-go-columns:v1.2 >/dev/null
        docker image inspect kazukinakamae/compute_go_kfi:v1.2 >/dev/null
        docker image inspect kazukinakamae/count_go_parent_terms:v1.1 >/dev/null
    fi
}

preflight_one() {
    local run_key="$1"
    local config_rel="$2"
    local source_run="$3"
    local target_run="$4"
    local config_path="${ROOT}/${config_rel}"
    local target_dir="${OUTPUT_ROOT}/${target_run}"
    [[ -f "${config_path}" ]] || {
        printf 'Configuration is missing: %s\n' "${config_path}" >&2
        return 1
    }
    check_source_artifacts "${source_run}"
    if [[ -e "${target_dir}" \
        && ! -f "${target_dir}/.v1.3.4-rerun-seed.tsv" ]]; then
        printf 'Target exists without a v1.3.4 seed marker: %s\n' \
            "${target_dir}" >&2
        return 1
    fi
    if [[ -f "${target_dir}/.v1.3.4-rerun-seed.tsv" ]]; then
        validate_existing_seed "${run_key}" "${config_rel}" \
            "${source_run}" "${target_run}" || {
            printf 'Existing v1.3.4 seed failed validation: %s\n' \
                "${target_dir}" >&2
            return 1
        }
    fi
    printf 'READY\t%s\t%s\n' "${run_key}" "${target_run}"
}

prepare_one() {
    local run_key="$1"
    local workflow="$2"
    local config_rel="$3"
    local source_run="$4"
    local target_run="$5"
    local source_intermediate="${OUTPUT_ROOT}/${source_run}/intermediate"
    local target_dir="${OUTPUT_ROOT}/${target_run}"
    local target_intermediate="${target_dir}/intermediate"
    local marker="${target_dir}/.v1.3.4-rerun-seed.tsv"
    local config_path="${ROOT}/${config_rel}"
    local temp_dir item source_digest target_digest

    if [[ -f "${marker}" ]]; then
        validate_existing_seed "${run_key}" "${config_rel}" \
            "${source_run}" "${target_run}"
        printf 'REUSED_VALIDATED_SEED\t%s\t%s\n' \
            "${run_key}" "${target_run}"
        return
    fi
    [[ ! -e "${target_dir}" ]] || {
        printf 'Refusing to replace existing target: %s\n' "${target_dir}" >&2
        return 1
    }

    temp_dir="$(mktemp -d "${OUTPUT_ROOT}/.v1.3.4-seed.XXXXXX")"
    mkdir -p "${temp_dir}/intermediate" "${temp_dir}/web_data"
    {
        printf 'field\tvalue\n'
        printf 'schema\t1\n'
        printf 'run_key\t%s\n' "${run_key}"
        printf 'workflow\t%s\n' "${workflow}"
        printf 'tWhEEL_version\t1.3.4\n'
        printf 'prepared_at_utc\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        printf 'code_commit\t%s\n' "$(git -C "${ROOT}" rev-parse HEAD)"
        printf 'guide_design_sha256\t%s\n' \
            "$(shasum -a 256 "${ROOT}/src/modules/guide_design.py" | awk '{print $1}')"
        printf 'config\t%s\n' "${config_rel}"
        printf 'config_sha256\t%s\n' \
            "$(shasum -a 256 "${config_path}" | awk '{print $1}')"
        printf 'source_run\t%s\n' "${source_run}"
        printf 'target_run\t%s\n' "${target_run}"
        printf 'restart_step\tguide_design\n'
    } > "${temp_dir}/.v1.3.4-rerun-seed.tsv"

    for item in "${SEED_ITEMS[@]}"; do
        source_digest="$(item_digest "${source_intermediate}/${item}")"
        cp -cR "${source_intermediate}/${item}" "${temp_dir}/intermediate/"
        target_digest="$(item_digest "${temp_dir}/intermediate/${item}")"
        [[ "${source_digest}" == "${target_digest}" ]] || {
            printf 'Clone checksum mismatch for %s; partial seed retained at %s\n' \
                "${item}" "${temp_dir}" >&2
            return 1
        }
        printf 'artifact_sha256\t%s\t%s\n' \
            "${item}" "${source_digest}" \
            >> "${temp_dir}/.v1.3.4-rerun-seed.tsv"
    done

    mv "${temp_dir}" "${target_dir}"
    printf 'PREPARED\t%s\t%s\n' "${run_key}" "${target_run}"
}

verify_step() {
    local workflow="$1"
    local step="$2"
    local target_dir="$3"
    local intermediate="${target_dir}/intermediate"
    case "${step}" in
        guide_design)
            [[ -s "${intermediate}/guide_design_results.txt" ]]
            [[ -e "${intermediate}/multi_genes_targets_specific_clean.bed" ]]
            [[ -e "${intermediate}/single_genes_targets_specific_clean.bed" ]]
            [[ -s "${intermediate}/normalized_gene_guide_current.json" ]]
            grep -q 'CrisFlash completion check passed for multi_genes' \
                "${target_dir}/twheel.log"
            grep -q 'CrisFlash completion check passed for single_genes' \
                "${target_dir}/twheel.log"
            PYTHONPATH="${ROOT}/src" "${PYTHON}" - \
                "${intermediate}" <<'PY'
import os
import sys
from modules.gene_guide_normalization import (
    resolve_current_normalized_gene_guide_paths,
)

paths = resolve_current_normalized_gene_guide_paths(sys.argv[1])
if not paths or not all(os.path.isfile(path) for path in paths.values()):
    raise SystemExit(1)
PY
            ;;
        go_analysis)
            for file_name in go_analysis_results.txt go_kfi_table.csv \
                go_kfi_stat.csv go_reverse_lookup.tsv go_parent_counts.csv \
                go_parent_counts.png; do
                [[ -s "${intermediate}/${file_name}" ]]
            done
            ;;
        pathway_analysis)
            for file_name in pathway_analysis_results.txt \
                pathway_kfi_table.csv pathway_gene_mapping.tsv \
                pathway_gene_match_audit.tsv \
                pathway_annotation_identifier_audit.tsv; do
                [[ -s "${intermediate}/${file_name}" ]]
            done
            ;;
        report_generation)
            [[ "${workflow}" == "GO" ]]
            [[ -s "${target_dir}/report.html" ]]
            [[ -s "${target_dir}/web_data/summary.json" ]]
            "${PYTHON}" -m json.tool \
                "${target_dir}/web_data/summary.json" >/dev/null
            ;;
        pathway_report_generation)
            [[ "${workflow}" == "PW" ]]
            [[ -s "${target_dir}/pathway_report.html" ]]
            [[ -s "${target_dir}/web_data/pathway_summary.json" ]]
            "${PYTHON}" -m json.tool \
                "${target_dir}/web_data/pathway_summary.json" >/dev/null
            ;;
        *)
            printf 'Unknown verification step: %s\n' "${step}" >&2
            return 1
            ;;
    esac
}

append_status() {
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$1" "$2" "$3" "$4" "$5" "$6" "$7" >> "${STATUS_FILE}"
}

run_step() {
    local run_key="$1"
    local workflow="$2"
    local config_path="$3"
    local target_run="$4"
    local step="$5"
    local target_dir="${OUTPUT_ROOT}/${target_run}"
    local started_at finished_at stamp log_file exit_code
    if verify_step "${workflow}" "${step}" "${target_dir}" 2>/dev/null; then
        append_status "${run_key}" "${step}" "reused_complete" \
            "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "0" "${target_run}" ""
        return
    fi

    started_at="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    stamp="$(date -u '+%Y%m%dT%H%M%SZ')"
    log_file="${RUN_LOG_DIR}/${run_key}_${step}_${stamp}.log"
    {
        printf 'Repository commit: %s\n' "$(git -C "${ROOT}" rev-parse HEAD)"
        printf 'tWhEEL version: 1.3.4\n'
        printf 'Guide-design SHA-256: '
        shasum -a 256 "${ROOT}/src/modules/guide_design.py"
        printf 'Config SHA-256: '
        shasum -a 256 "${config_path}"
        printf 'Command: %q %q --config %q --step %q --name %q\n' \
            "${PYTHON}" "${ROOT}/src/twheel.py" "${config_path}" \
            "${step}" "${target_run}"
    } > "${log_file}"

    set +e
    "${PYTHON}" "${ROOT}/src/twheel.py" \
        --config "${config_path}" --step "${step}" --name "${target_run}" \
        >> "${log_file}" 2>&1
    exit_code=$?
    set -e
    finished_at="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    if (( exit_code != 0 )) \
        || ! verify_step "${workflow}" "${step}" "${target_dir}"; then
        append_status "${run_key}" "${step}" "failed" "${finished_at}" \
            "${exit_code}" "${target_run}" "${log_file}"
        printf 'Step failed validation: %s %s (log: %s)\n' \
            "${run_key}" "${step}" "${log_file}" >&2
        return 1
    fi
    append_status "${run_key}" "${step}" "completed" "${finished_at}" \
        "0" "${target_run}" "${log_file}"
    printf 'COMPLETED\t%s\t%s\t%s\t%s\n' \
        "${run_key}" "${step}" "${started_at}" "${finished_at}"
}

execute_one() {
    local run_key="$1"
    local workflow="$2"
    local config_rel="$3"
    local target_run="$4"
    local downstream_steps="$5"
    local config_path="${ROOT}/${config_rel}"
    local step
    run_step "${run_key}" "${workflow}" "${config_path}" \
        "${target_run}" guide_design
    IFS=',' read -r -a steps <<< "${downstream_steps}"
    for step in "${steps[@]}"; do
        run_step "${run_key}" "${workflow}" "${config_path}" \
            "${target_run}" "${step}"
    done
}

if [[ "${MODE}" == "list" ]]; then
    list_runs
    exit 0
fi

check_common_preflight
matched=0
while IFS=$'\t' read -r run_key workflow config_rel source_run \
    target_run downstream_steps; do
    [[ "${run_key}" == "run_key" ]] && continue
    selected "${run_key}" || continue
    matched=$((matched + 1))
    preflight_one "${run_key}" "${config_rel}" "${source_run}" "${target_run}"
done < "${MANIFEST}"

(( matched > 0 )) || {
    printf 'No manifest entry matched run key: %s\n' "${SELECTED_RUN}" >&2
    exit 2
}

if [[ "${MODE}" == "preflight" ]]; then
    printf 'Preflight passed for %d publication run(s). No files were changed.\n' \
        "${matched}"
    exit 0
fi

while IFS=$'\t' read -r run_key workflow config_rel source_run \
    target_run downstream_steps; do
    [[ "${run_key}" == "run_key" ]] && continue
    selected "${run_key}" || continue
    prepare_one "${run_key}" "${workflow}" "${config_rel}" \
        "${source_run}" "${target_run}"
done < "${MANIFEST}"

if [[ "${MODE}" == "prepare" ]]; then
    printf 'Prepared %d publication run(s). No analysis step was executed.\n' \
        "${matched}"
    exit 0
fi

mkdir -p "${RUN_LOG_DIR}"
if ! mkdir "${LOCK_DIR}" 2>/dev/null; then
    printf 'Another rerun or stale lock requires diagnosis: %s\n' \
        "${LOCK_DIR}" >&2
    exit 73
fi
trap 'rmdir "${LOCK_DIR}" 2>/dev/null || true' EXIT INT TERM
if [[ ! -f "${STATUS_FILE}" ]]; then
    printf 'run_key\tstep\tstatus\tfinished_at_utc\texit_code\ttarget_run\tlog_file\n' \
        > "${STATUS_FILE}"
fi

while IFS=$'\t' read -r run_key workflow config_rel source_run \
    target_run downstream_steps; do
    [[ "${run_key}" == "run_key" ]] && continue
    selected "${run_key}" || continue
    execute_one "${run_key}" "${workflow}" "${config_rel}" \
        "${target_run}" "${downstream_steps}"
done < "${MANIFEST}"

printf 'All selected v1.3.4 publication reruns completed and passed validation.\n'
