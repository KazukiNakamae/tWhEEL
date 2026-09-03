#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MANIFEST="${SCRIPT_DIR}/publication_runs.tsv"
OUTPUT_ROOT="${ROOT}/outputs/benchmarking"
PYTHON="${TWHEEL_V134_PYTHON:-${ROOT}/.handoff_runtime/twheel_env-osx-64/bin/python}"
failures=0

printf 'run_key\tworkflow\tstatus\ttarget_run\n'
while IFS=$'\t' read -r run_key workflow _config_rel _source_run \
    target_run _downstream_steps; do
    [[ "${run_key}" == "run_key" ]] && continue
    target_dir="${OUTPUT_ROOT}/${target_run}"
    intermediate="${target_dir}/intermediate"
    status="complete"

    for path in \
        "${target_dir}/.v1.3.4-rerun-seed.tsv" \
        "${intermediate}/guide_design_results.txt" \
        "${intermediate}/normalized_gene_guide_current.json"; do
        [[ -s "${path}" ]] || status="incomplete"
    done
    if [[ "${workflow}" == "GO" ]]; then
        for path in \
            "${intermediate}/go_kfi_table.csv" \
            "${intermediate}/go_kfi_stat.csv" \
            "${target_dir}/report.html" \
            "${target_dir}/web_data/summary.json"; do
            [[ -s "${path}" ]] || status="incomplete"
        done
        if [[ -s "${target_dir}/web_data/summary.json" ]]; then
            "${PYTHON}" -m json.tool \
                "${target_dir}/web_data/summary.json" >/dev/null \
                || status="invalid_json"
        fi
    else
        for path in \
            "${intermediate}/pathway_kfi_table.csv" \
            "${intermediate}/pathway_gene_mapping.tsv" \
            "${target_dir}/pathway_report.html" \
            "${target_dir}/web_data/pathway_summary.json"; do
            [[ -s "${path}" ]] || status="incomplete"
        done
        if [[ -s "${target_dir}/web_data/pathway_summary.json" ]]; then
            "${PYTHON}" -m json.tool \
                "${target_dir}/web_data/pathway_summary.json" >/dev/null \
                || status="invalid_json"
        fi
    fi
    if [[ -f "${target_dir}/twheel.log" ]] \
        && grep -Eq ' - ERROR - |Traceback' "${target_dir}/twheel.log"; then
        status="error_in_log"
    fi
    if [[ "${status}" != "complete" ]]; then
        failures=$((failures + 1))
    fi
    printf '%s\t%s\t%s\t%s\n' \
        "${run_key}" "${workflow}" "${status}" "${target_run}"
done < "${MANIFEST}"

if (( failures > 0 )); then
    printf 'Verification failed for %d publication run(s).\n' \
        "${failures}" >&2
    exit 1
fi
printf 'Verification passed for all publication reruns.\n'
