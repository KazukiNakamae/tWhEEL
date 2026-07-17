#!/usr/bin/env python3
"""
Pathway Report Generation Module
HTML report generation and source staging for PW-KFI analysis results
"""

import os
import json
import logging
import csv
import ast
import copy
import hashlib
import shutil
import html
import re
import sys
import math
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime
from pathlib import Path

from .guide_design_contract import CRISFLASH_MAX_MISMATCH
from .report_source_store import ReportSourceStore
try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2 ** 31 - 1)


class PathwayReportGenerator:
    """Pathway report generation processing class"""

    REPORT_SOURCE_FILES = {
        'pathway_kfi_table': 'pathway_kfi_table.csv',
        'pathway_gene_mapping': 'pathway_gene_mapping.tsv',
    }
    
    def __init__(self, config):
        """
        Initialization
        
        Args:
            config (dict): Configuration Dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.intermediate_dir = self.config['output']['intermediate']
        self.output_dir = self.config['output']['base_dir']
        self.web_data_dir = os.path.join(self.output_dir, 'web_data')
        self.pathway_report_config = self.config.get('pathway_report', {})
        self.report_source_store = ReportSourceStore(self.web_data_dir)
        
        # 出力ディレクトリを作成
        Path(self.web_data_dir).mkdir(parents=True, exist_ok=True)

    def load_pathway_analysis_results(self):
        """Load current artifacts or a complete verified pathway snapshot."""
        intermediate_dir = Path(self.intermediate_dir)
        current_paths = {
            key: intermediate_dir / filename
            for key, filename in self.REPORT_SOURCE_FILES.items()
        }

        if intermediate_dir.exists():
            if not intermediate_dir.is_dir():
                raise NotADirectoryError(
                    f"Intermediate path is not a directory: {intermediate_dir}"
                )
            current_present = {
                key: path.is_file() for key, path in current_paths.items()
            }
            if not all(current_present.values()):
                missing = sorted(
                    key for key, present in current_present.items()
                    if not present
                )
                prefix = (
                    "Incomplete pathway report source snapshot in intermediate"
                    if any(current_present.values())
                    else "Pathway report sources are absent from intermediate"
                )
                raise FileNotFoundError(prefix + ': ' + ', '.join(missing))

            pathways, pathway_genes = self._load_pathway_analysis_results(
                current_paths['pathway_kfi_table'],
                current_paths['pathway_gene_mapping'],
            )
            sources = {
                key: (path, f'pathway/{self.REPORT_SOURCE_FILES[key]}')
                for key, path in current_paths.items()
            }
            for pathway in pathways:
                pathway_id = str(
                    pathway.get('pathway_id', '') or ''
                ).strip()
                source_path = Path(
                    str(pathway.get('gpml_file_path', '') or '')
                )
                if source_path.is_file():
                    logical_key, relative_path = self._gpml_source_entry(
                        pathway_id
                    )
                    sources[logical_key] = (source_path, relative_path)

            published = self.report_source_store.publish_namespace(
                'pathway', sources
            )
            self._connect_pathway_gpml(pathways, published)
            return pathways, pathway_genes

        snapshot = self.report_source_store.load_namespace('pathway')
        missing = sorted(set(self.REPORT_SOURCE_FILES) - set(snapshot))
        if missing:
            raise FileNotFoundError(
                "Pathway report source snapshot is incomplete: "
                + ', '.join(missing)
            )
        pathways, pathway_genes = self._load_pathway_analysis_results(
            snapshot['pathway_kfi_table'],
            snapshot['pathway_gene_mapping'],
        )
        self._connect_pathway_gpml(pathways, snapshot)
        return pathways, pathway_genes

    @staticmethod
    def _gpml_source_entry(pathway_id):
        """Return a collision-resistant logical key and GPML destination."""
        safe_id = re.sub(
            r'[^A-Za-z0-9_.-]+', '_', pathway_id
        ).strip('._') or 'pathway'
        digest = hashlib.sha256(pathway_id.encode('utf-8')).hexdigest()[:12]
        return (
            f'gpml.{pathway_id}',
            f'pathway/gpml/{safe_id}-{digest}.gpml',
        )

    @classmethod
    def _connect_pathway_gpml(cls, pathways, source_paths):
        """Point pathway rows only at GPML files from the selected snapshot."""
        for pathway in pathways:
            pathway_id = str(pathway.get('pathway_id', '') or '').strip()
            logical_key, _ = cls._gpml_source_entry(pathway_id)
            pathway['gpml_file_path'] = source_paths.get(logical_key, '')

    def _load_pathway_analysis_results(self, kfi_table_path, mapping_path):
        """Read and validate one canonical pathway table pair."""

        required_pathway_columns = {
            'pathway_id', 'pathway_name', 'total_gene_nodes',
            'evaluable_gene_node_count',
            'unique_xref_count', 'total_genes', 'gene_count_basis',
            'canonicalization_complete', 'canonicalization_status',
            'matched_genes', 'unmatched_genes', 'unmatchable_genes',
            'unmatched_gene_count', 'matched_gene_count',
            'genes_with_guides', 'pw_kfi_score', 'coverage_rate',
            'unmatchable_gene_count',
        }
        pathways = []
        expected_mapping_by_pathway = {}
        expected_pathway_names = {}
        with open(kfi_table_path, 'r', encoding='utf-8', newline='') as f:
            reader = csv.DictReader(f)
            missing_columns = required_pathway_columns - set(
                reader.fieldnames or []
            )
            if missing_columns:
                raise ValueError(
                    "Legacy or incomplete PW-KFI table cannot be reported as "
                    "canonical-gene results. Regenerate pathway analysis; "
                    "missing columns: " + ', '.join(sorted(missing_columns))
                )
            for row in reader:
                pathway_id = row.get('pathway_id', '')
                if not pathway_id:
                    raise ValueError('PW-KFI table contains an empty pathway_id')
                if pathway_id in expected_mapping_by_pathway:
                    raise ValueError(
                        f"PW-KFI table contains duplicate pathway_id: {pathway_id}"
                    )
                if row.get('gene_count_basis') != 'exact_canonical_gene_units':
                    raise ValueError(
                        f"Unsupported gene_count_basis for {pathway_id}: "
                        f"{row.get('gene_count_basis')!r}. Regenerate pathway analysis."
                    )
                if (
                    str(row.get('canonicalization_complete', '')).lower()
                    != 'true'
                    or row.get('canonicalization_status') != 'complete'
                ):
                    raise ValueError(
                        f"PW-KFI is not canonically resolved for {pathway_id}"
                    )

                counts = {
                    field: self._nonnegative_int(
                        row.get(field, ''), field, pathway_id
                    )
                    for field in (
                        'total_gene_nodes', 'evaluable_gene_node_count',
                        'unique_xref_count', 'total_genes',
                        'unmatchable_gene_count', 'matched_gene_count',
                        'unmatched_gene_count', 'genes_with_guides',
                    )
                }
                if not (
                    counts['total_genes'] <= counts['unique_xref_count']
                    <= counts['evaluable_gene_node_count']
                    <= counts['total_gene_nodes']
                ):
                    raise ValueError(
                        f"Invalid pathway count ordering for {pathway_id}: {counts}"
                    )
                if (
                    counts['evaluable_gene_node_count']
                    + counts['unmatchable_gene_count']
                    != counts['total_gene_nodes']
                ):
                    raise ValueError(
                        f"Raw node counts do not reconcile for {pathway_id}"
                    )
                if (
                    counts['matched_gene_count']
                    + counts['unmatched_gene_count']
                    != counts['total_genes']
                ):
                    raise ValueError(
                        f"Canonical gene counts do not reconcile for {pathway_id}"
                    )
                if counts['genes_with_guides'] > counts['matched_gene_count']:
                    raise ValueError(
                        f"Guide-positive genes exceed matched genes for {pathway_id}"
                    )

                serialized_gene_fields = {}
                for field in (
                    'matched_genes', 'unmatched_genes', 'unmatchable_genes'
                ):
                    try:
                        value = self._parse_matched_genes(row.get(field, ''))
                    except (
                        TypeError, ValueError, SyntaxError,
                        json.JSONDecodeError,
                    ) as exc:
                        raise ValueError(
                            f"Cannot deserialize {field} for {pathway_id}"
                        ) from exc
                    if not isinstance(value, list):
                        raise ValueError(
                            f"{field} must be a list for {pathway_id}"
                        )
                    serialized_gene_fields[field] = value

                matched_genes = serialized_gene_fields['matched_genes']
                unmatched_genes = serialized_gene_fields['unmatched_genes']
                unmatchable_genes = serialized_gene_fields['unmatchable_genes']
                if len(matched_genes) != counts['matched_gene_count']:
                    raise ValueError(
                        f"matched_genes does not reconcile for {pathway_id}"
                    )
                if len(unmatched_genes) != counts['unmatched_gene_count']:
                    raise ValueError(
                        f"unmatched_genes does not reconcile for {pathway_id}"
                    )
                if len(unmatchable_genes) != counts['unmatchable_gene_count']:
                    raise ValueError(
                        f"unmatchable_genes does not reconcile for {pathway_id}"
                    )
                if any(not isinstance(gene, dict) for gene in unmatchable_genes):
                    raise ValueError(
                        f"unmatchable_genes contains an invalid record for "
                        f"{pathway_id}"
                    )

                canonical_keys = []
                source_xref_keys = []
                evaluable_node_count = 0
                guide_positive_count = 0
                pathway_source_databases = set()
                unit_canonicalization_statuses = []
                mapping_signatures = {}
                for collection_name, genes, matched in (
                    ('matched_genes', matched_genes, True),
                    ('unmatched_genes', unmatched_genes, False),
                ):
                    for index, gene in enumerate(genes):
                        if not isinstance(gene, dict):
                            raise ValueError(
                                f"{collection_name}[{index}] is invalid for "
                                f"{pathway_id}"
                            )
                        validation = self._validate_canonical_gene_unit(
                            gene, pathway_id, matched=matched
                        )
                        canonical_keys.append(validation['gene_key'])
                        source_xref_keys.extend(validation['source_xref_keys'])
                        pathway_source_databases.update(
                            validation['source_databases']
                        )
                        unit_canonicalization_statuses.append(
                            validation['canonicalization_complete']
                        )
                        evaluable_node_count += validation[
                            'node_occurrence_count'
                        ]
                        guide_positive_count += int(
                            matched and validation['has_guides']
                        )
                        if matched:
                            mapping_signatures[validation['gene_key']] = (
                                self._mapping_signature(gene, pathway_id)
                            )

                if (
                    len(canonical_keys) != counts['total_genes']
                    or len(set(canonical_keys)) != len(canonical_keys)
                ):
                    raise ValueError(
                        f"Canonical gene keys do not reconcile for {pathway_id}"
                    )
                if (
                    len(source_xref_keys) != counts['unique_xref_count']
                    or len(set(source_xref_keys)) != len(source_xref_keys)
                ):
                    raise ValueError(
                        f"Source xrefs do not reconcile for {pathway_id}"
                    )
                if evaluable_node_count != counts['evaluable_gene_node_count']:
                    raise ValueError(
                        f"Evaluable node occurrences do not reconcile for "
                        f"{pathway_id}"
                    )
                if guide_positive_count != counts['genes_with_guides']:
                    raise ValueError(
                        f"Guide-positive genes do not reconcile for {pathway_id}"
                    )
                if (
                    len(pathway_source_databases) > 1
                    and not all(unit_canonicalization_statuses)
                ):
                    raise ValueError(
                        f"Cross-database canonicalization does not reconcile "
                        f"for {pathway_id}"
                    )
                expected_mapping_by_pathway[pathway_id] = mapping_signatures
                expected_pathway_names[pathway_id] = row.get(
                    'pathway_name', ''
                )

                try:
                    pw_kfi_score = float(row.get('pw_kfi_score', ''))
                    coverage_rate = float(row.get('coverage_rate', ''))
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"Invalid pathway rate for {pathway_id}"
                    ) from exc
                if not math.isfinite(pw_kfi_score) or not math.isfinite(
                    coverage_rate
                ):
                    raise ValueError(f"Non-finite pathway rate for {pathway_id}")
                expected_kfi = (
                    counts['genes_with_guides'] / counts['total_genes']
                    if counts['total_genes'] else 0.0
                )
                expected_coverage = (
                    counts['matched_gene_count'] / counts['total_genes']
                    if counts['total_genes'] else 0.0
                )
                if abs(pw_kfi_score - expected_kfi) > 1e-12:
                    raise ValueError(
                        f"PW-KFI does not reconcile for {pathway_id}"
                    )
                if abs(coverage_rate - expected_coverage) > 1e-12:
                    raise ValueError(
                        f"Annotation match coverage does not reconcile for {pathway_id}"
                    )
                pathways.append(row)

        pathway_genes = {}
        if not os.path.exists(mapping_path):
            raise FileNotFoundError(
                f"Canonical pathway gene mapping not found.: {mapping_path}"
            )
        required_mapping_columns = {
            'pathway_id', 'pathway_name', 'pathway_gene_database',
            'pathway_gene_id', 'pathway_gene_name', 'gene_key',
            'pw_kfi_gene_key', 'pathway_xrefs', 'source_databases',
            'node_occurrence_count', 'unique_xref_count',
            'canonical_subject_id', 'canonicalization_methods',
            'canonicalization_complete', 'annotation_source_database',
            'matched_annotation_key', 'matched_annotation_id',
            'resolved_annotation_subject_id', 'match_method',
            'mapping_source', 'sgd_id', 'ensembl_id', 'gene_sgd_id',
            'unique_id', 'gene_name', 'has_guides', 'guide_count',
        }
        with open(mapping_path, 'r', encoding='utf-8', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            missing_columns = required_mapping_columns - set(
                reader.fieldnames or []
            )
            if missing_columns:
                raise ValueError(
                    "Legacy or incomplete pathway mapping cannot be reported. "
                    "Regenerate pathway analysis; missing columns: "
                    + ', '.join(sorted(missing_columns))
                )
            for row in reader:
                pathway_id = row['pathway_id']
                if pathway_id not in expected_mapping_by_pathway:
                    raise ValueError(
                        f"Mapping references an unknown pathway: {pathway_id}"
                    )
                if row.get('pathway_name', '') != expected_pathway_names[pathway_id]:
                    raise ValueError(
                        f"Pathway mapping name does not reconcile for {pathway_id}"
                    )
                gene_key = row.get('pw_kfi_gene_key', '')
                expected_signature = expected_mapping_by_pathway[pathway_id].get(
                    gene_key
                )
                if expected_signature is None:
                    raise ValueError(
                        f"Pathway mapping key does not reconcile for "
                        f"{pathway_id}: {gene_key!r}"
                    )
                if self._mapping_signature(row, pathway_id) != expected_signature:
                    raise ValueError(
                        f"Pathway mapping payload does not reconcile for "
                        f"{pathway_id}: {gene_key}"
                    )
                pathway_genes.setdefault(pathway_id, []).append(row)

        for pathway_id, expected_mapping in expected_mapping_by_pathway.items():
            mapping_keys = [
                row.get('pw_kfi_gene_key', '')
                for row in pathway_genes.get(pathway_id, [])
            ]
            if (
                len(mapping_keys) != len(set(mapping_keys))
                or set(mapping_keys) != set(expected_mapping)
            ):
                raise ValueError(
                    f"Pathway mapping does not reconcile for {pathway_id}"
                )
        
        self.logger.info(f"Pathway data loading complete.: {len(pathways)} Pathway, {len(pathway_genes)} associated genes")
        return pathways, pathway_genes

    
    
    def parse_gpml_file(self, gpml_file_path):
        """
        Analyze GPML files to extract pathway diagram data.
        
        Args:
            gpml_file_path (str): Path to the GPML file
            
        Returns:
            dict: Pathway diagram data (nodes, edges, metadata)
        """
        try:
            if not os.path.exists(gpml_file_path):
                self.logger.warning(f"GPML file not found.: {gpml_file_path}")
                return None
            
            # XMLファイルを解析
            tree = ET.parse(gpml_file_path)
            root = tree.getroot()
            
            # 名前空間を取得
            ns = {'gpml': 'http://pathvisio.org/GPML/2013a'} if root.tag.startswith('{') else {}
            
            # パスウェイメタデータを抽出
            pathway_name = root.get('Name', 'Unknown Pathway')
            organism = root.get('Organism', 'Unknown Organism')
            
            # Graphics要素からキャンバスサイズを取得
            graphics_elem = root.find('.//gpml:Graphics', ns) if ns else root.find('.//Graphics')
            board_width = float(graphics_elem.get('BoardWidth', 800)) if graphics_elem is not None else 800
            board_height = float(graphics_elem.get('BoardHeight', 600)) if graphics_elem is not None else 600
            
            # DataNodeを抽出（遺伝子・代謝物ノード）
            nodes = []
            data_nodes = root.findall('.//gpml:DataNode', ns) if ns else root.findall('.//DataNode')
            
            for node in data_nodes:
                node_type = node.get('Type', 'Unknown')
                node_id = node.get('GraphId', '')
                text_label = node.get('TextLabel', '')
                
                # Graphics情報を取得
                graphics = node.find('gpml:Graphics', ns) if ns else node.find('Graphics')
                center_x = float(graphics.get('CenterX', 0)) if graphics is not None else 0
                center_y = float(graphics.get('CenterY', 0)) if graphics is not None else 0
                width = float(graphics.get('Width', 80)) if graphics is not None else 80
                height = float(graphics.get('Height', 20)) if graphics is not None else 20
                
                # Xref情報を取得（SGD IDなど）
                xref = node.find('gpml:Xref', ns) if ns else node.find('Xref')
                database = (
                    xref.get('Database', '').strip() if xref is not None else ''
                )
                xref_id = xref.get('ID', '').strip() if xref is not None else ''
                
                node_data = {
                    'id': node_id,
                    'type': node_type,
                    'label': text_label,
                    'x': center_x,
                    'y': center_y,
                    'width': width,
                    'height': height,
                    'database': database,
                    'xref_id': xref_id
                }
                
                nodes.append(node_data)
            
            # Interactionを抽出（相互作用エッジ）
            edges = []
            interactions = root.findall('.//gpml:Interaction', ns) if ns else root.findall('.//Interaction')
            
            for interaction in interactions:
                interaction_id = interaction.get('GraphId', '')
                
                # Graphics情報を取得
                graphics = interaction.find('gpml:Graphics', ns) if ns else interaction.find('Graphics')
                connector_type = graphics.get('ConnectorType', 'Straight') if graphics is not None else 'Straight'
                
                # Point情報を取得
                points = []
                if graphics is not None:
                    point_elements = graphics.findall('gpml:Point', ns) if ns else graphics.findall('Point')
                    for point in point_elements:
                        x = float(point.get('X', 0))
                        y = float(point.get('Y', 0))
                        graph_ref = point.get('GraphRef', '')
                        points.append({'x': x, 'y': y, 'ref': graph_ref})
                
                edge_data = {
                    'id': interaction_id,
                    'type': 'interaction',
                    'connector_type': connector_type,
                    'points': points
                }
                
                edges.append(edge_data)
            
            pathway_data = {
                'pathway_name': pathway_name,
                'organism': organism,
                'board_width': board_width,
                'board_height': board_height,
                'nodes': nodes,
                'edges': edges,
                'gpml_file': gpml_file_path
            }
            
            self.logger.debug(f"GPML analysis complete: {pathway_name} - {len(nodes)} node, {len(edges)} edge")
            return pathway_data
            
        except Exception as e:
            self.logger.error(f"GPML parsing error ({gpml_file_path}): {e}")
            return None
    
    def generate_pathway_svg(self, pathway_data, guide_status_map=None):
        """
        Generate SVG strings from pathway data.
        
        Args:
            pathway_data (dict): Data returned from parse_gpml_file
            guide_status_map (dict): Exact ``Database::ID`` -> guide status mapping.
            
        Returns:
            str: SVG string
        """
        if not pathway_data:
            return None
        
        board_width = pathway_data['board_width']
        board_height = pathway_data['board_height']
        nodes = pathway_data['nodes']
        edges = pathway_data['edges']
        
        # SVGヘッダー
        svg_parts = [
            f'<svg width="{board_width}" height="{board_height}" viewBox="0 0 {board_width} {board_height}" xmlns="http://www.w3.org/2000/svg">',
            '<defs>',
            '<style>',
            '.gene-node { stroke-width: 2; cursor: pointer; }',
            '.gene-with-guides { fill: #28a745; stroke: #1e7e34; }',
            '.gene-without-guides { fill: #dc3545; stroke: #bd2130; }',
            '.gene-unknown { fill: #6c757d; stroke: #545b62; }',
            '.metabolite-node { fill: #17a2b8; stroke: #117a8b; stroke-width: 2; }',
            '.interaction-line { stroke: #333; stroke-width: 1.5; fill: none; }',
            '.node-label { font-family: Arial, sans-serif; font-size: 10px; text-anchor: middle; pointer-events: none; }',
            '</style>',
            '</defs>'
        ]
        
        # エッジ（相互作用線）を先に描画
        for edge in edges:
            points = edge['points']
            if len(points) >= 2:
                path_data = f"M {points[0]['x']} {points[0]['y']}"
                for point in points[1:]:
                    path_data += f" L {point['x']} {point['y']}"
                
                svg_parts.append(f'<path d="{path_data}" class="interaction-line" />')
        
        # ノード（遺伝子・代謝物）を描画
        for node in nodes:
            x, y = node['x'], node['y']
            width, height = node['width'], node['height']
            node_type = node['type']
            label = node['label']
            xref_id = node['xref_id']
            database = node['database']
            gene_key = f"{database}::{xref_id}"
            escaped_label = html.escape(str(label), quote=True)
            escaped_database = html.escape(str(database), quote=True)
            escaped_xref_id = html.escape(str(xref_id), quote=True)
            escaped_gene_key = html.escape(gene_key, quote=True)
            
            # Color only an exact database-and-identifier match.
            if node_type == 'GeneProduct':
                gene_status = (
                    guide_status_map.get(gene_key)
                    if guide_status_map else None
                )
                if (
                    gene_status is not None
                    and gene_status.get('match_status', 'matched') == 'matched'
                ):
                    if gene_status['has_guides']:
                        node_class = 'gene-with-guides'
                    else:
                        node_class = 'gene-without-guides'
                else:
                    node_class = 'gene-unknown'
                
                # 遺伝子ノード（楕円形）
                svg_parts.append(
                    f'<ellipse cx="{x}" cy="{y}" rx="{width/2}" ry="{height/2}" '
                    f'class="gene-node {node_class}" '
                    f'data-gene-database="{escaped_database}" '
                    f'data-gene-id="{escaped_xref_id}" '
                    f'data-gene-key="{escaped_gene_key}" '
                    f'data-sgd-id="{html.escape(str((gene_status or {}).get("sgd_id", "")), quote=True)}" '
                    f'data-gene-name="{escaped_label}" '
                    f'title="{escaped_label} ({escaped_database}:{escaped_xref_id})" />'
                )
            
            elif node_type == 'Metabolite':
                # 代謝物ノード（四角形）
                svg_parts.append(
                    f'<rect x="{x - width/2}" y="{y - height/2}" width="{width}" height="{height}" '
                    f'class="metabolite-node" title="{escaped_label}" />'
                )
            
            else:
                # その他のノード（四角形、グレー）
                svg_parts.append(
                    f'<rect x="{x - width/2}" y="{y - height/2}" width="{width}" height="{height}" '
                    f'fill="#e9ecef" stroke="#adb5bd" stroke-width="1" title="{escaped_label}" />'
                )
            
            # ラベルテキスト
            # 長いラベルは切り詰める
            display_label = label[:12] + '...' if len(label) > 15 else label
            svg_parts.append(
                f'<text x="{x}" y="{y + 3}" class="node-label" fill="#000">'
                f'{html.escape(str(display_label))}</text>'
            )
        
        # 凡例を追加
        legend_y = board_height - 80
        svg_parts.extend([
            f'<g id="legend">',
            f'<rect x="10" y="{legend_y}" width="200" height="70" fill="white" stroke="#ccc" stroke-width="1" opacity="0.9" />',
            f'<text x="20" y="{legend_y + 15}" class="node-label" font-weight="bold">Legend</text>',
            f'<ellipse cx="30" cy="{legend_y + 30}" rx="8" ry="6" class="gene-node gene-with-guides" />',
            f'<text x="45" y="{legend_y + 33}" class="node-label" font-size="9px">With Guides</text>',
            f'<ellipse cx="30" cy="{legend_y + 45}" rx="8" ry="6" class="gene-node gene-without-guides" />',
            f'<text x="45" y="{legend_y + 48}" class="node-label" font-size="9px">Without Guides</text>',
            f'<ellipse cx="125" cy="{legend_y + 30}" rx="8" ry="6" class="gene-node gene-unknown" />',
            f'<text x="140" y="{legend_y + 33}" class="node-label" font-size="9px">Unknown</text>',
            f'<rect x="120" y="{legend_y + 39}" width="16" height="12" class="metabolite-node" />',
            f'<text x="140" y="{legend_y + 48}" class="node-label" font-size="9px">Metabolite</text>',
            f'</g>'
        ])
        
        svg_parts.append('</svg>')
        
        return '\n'.join(svg_parts)
    
    @staticmethod
    def _parse_matched_genes(matched_genes_value):
        """Parse the serialized matched-gene list without interpreting free text."""
        if isinstance(matched_genes_value, list):
            return matched_genes_value
        if not isinstance(matched_genes_value, str):
            raise TypeError('Serialized pathway gene data must be a string or list')
        if not matched_genes_value.strip():
            return []
        if matched_genes_value.startswith('[') and matched_genes_value.endswith(']'):
            return ast.literal_eval(matched_genes_value)
        return json.loads(matched_genes_value)

    @staticmethod
    def _exact_pathway_identifier(gene_info):
        """Return an explicit pathway database and ID, with SGD-only legacy support."""
        database = (gene_info.get('pathway_gene_database') or '').strip()
        gene_id = (gene_info.get('pathway_gene_id') or '').strip()
        if database and gene_id:
            return database, gene_id

        # Older pathway output was produced only from an exact SGD GPML ID.
        legacy_id = (gene_info.get('original_gene_id') or '').strip()
        if re.fullmatch(r'S\d{9}', legacy_id):
            return 'SGD', legacy_id
        return '', ''

    @classmethod
    def _pathway_xrefs(cls, gene_info):
        """Return every exact source xref represented by one PW-KFI gene."""
        value = gene_info.get('pathway_xrefs', [])
        if isinstance(value, str) and value.strip():
            value = cls._parse_matched_genes(value)
        if value:
            if not isinstance(value, list):
                raise ValueError('pathway_xrefs must be a list')
            xrefs = []
            for index, xref in enumerate(value):
                if not isinstance(xref, dict):
                    raise ValueError(f'pathway_xrefs[{index}] must be an object')
                database = str(xref.get('database', '') or '').strip()
                gene_id = str(xref.get('gene_id', '') or '').strip()
                if not database or not gene_id:
                    raise ValueError(
                        f'pathway_xrefs[{index}] lacks database or gene_id'
                    )
                expected_gene_key = f'{database}::{gene_id}'
                supplied_gene_key = str(
                    xref.get('gene_key', '') or ''
                ).strip()
                if supplied_gene_key and supplied_gene_key != expected_gene_key:
                    raise ValueError(
                        f'pathway_xrefs[{index}] has an inconsistent gene_key'
                    )
                normalized = xref.copy()
                normalized.update({
                    'database': database,
                    'gene_id': gene_id,
                    'gene_key': expected_gene_key,
                })
                xrefs.append(normalized)
            return sorted(
                xrefs, key=lambda xref: (xref['database'], xref['gene_id'])
            )

        database = str(
            gene_info.get('pathway_gene_database')
            or gene_info.get('database')
            or ''
        ).strip()
        gene_id = str(
            gene_info.get('pathway_gene_id')
            or gene_info.get('gene_id')
            or ''
        ).strip()
        if not database or not gene_id:
            database, gene_id = cls._exact_pathway_identifier(gene_info)
        if not database or not gene_id:
            return []
        return [{
            'database': database,
            'gene_id': gene_id,
            'gene_key': f'{database}::{gene_id}',
            'pathway_gene_name': gene_info.get('pathway_gene_name', ''),
            'match_status': gene_info.get('match_status', 'matched'),
            'reason': gene_info.get('reason', ''),
        }]

    @staticmethod
    def _nonnegative_int(value, field_name, pathway_id, default=None):
        """Parse a report count without silently coercing invalid values."""
        if value is None or not str(value).strip():
            if default is not None:
                return default
            raise ValueError(f'Missing {field_name} for {pathway_id}')
        try:
            parsed = int(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f'Invalid {field_name} for {pathway_id}: {value!r}'
            ) from exc
        if parsed < 0:
            raise ValueError(f'Negative {field_name} for {pathway_id}')
        return parsed

    @staticmethod
    def _strict_bool(value, field_name, pathway_id):
        """Parse an artifact boolean without treating arbitrary text as false."""
        if isinstance(value, bool):
            return value
        normalized = str(value).strip().lower()
        if normalized == 'true':
            return True
        if normalized == 'false':
            return False
        raise ValueError(
            f'Invalid {field_name} for {pathway_id}: {value!r}'
        )

    @staticmethod
    def _normalized_string_list(value):
        """Normalize list-valued artifact fields written as CSV text."""
        if value is None:
            return ()
        if isinstance(value, str):
            items = value.split(',')
        elif isinstance(value, (list, tuple, set)):
            items = value
        else:
            raise ValueError(
                f'Expected a string or list, got {type(value).__name__}'
            )
        return tuple(sorted({
            str(item).strip() for item in items if str(item).strip()
        }))

    @classmethod
    def _freeze_artifact_value(cls, value):
        """Convert nested provenance data to an order-independent value."""
        if isinstance(value, dict):
            return tuple(sorted(
                (str(key), cls._freeze_artifact_value(item))
                for key, item in value.items()
            ))
        if isinstance(value, (list, tuple, set)):
            frozen = [cls._freeze_artifact_value(item) for item in value]
            return tuple(sorted(frozen, key=repr))
        return value

    @staticmethod
    def _guide_entries(value, pathway_id):
        """Normalize serialized guide entries for artifact reconciliation."""
        if value is None:
            return ()
        if isinstance(value, str):
            items = value.split(',')
        elif isinstance(value, (list, tuple, set)):
            items = value
        else:
            raise ValueError(
                f'Invalid guide_targets for {pathway_id}: {value!r}'
            )
        return tuple(sorted({
            str(item).strip() for item in items if str(item).strip()
        }))

    @classmethod
    def _validate_canonical_gene_unit(cls, gene, pathway_id, matched):
        """Validate one canonical PW-KFI denominator unit and its source nodes."""
        gene_key = str(gene.get('pw_kfi_gene_key', '') or '').strip()
        if not gene_key:
            raise ValueError(
                f'Canonical gene provenance is missing for {pathway_id}'
            )
        xrefs = cls._pathway_xrefs(gene)
        if not xrefs:
            raise ValueError(
                f'Canonical gene source xrefs are missing for {pathway_id}'
            )
        source_xref_keys = [xref['gene_key'] for xref in xrefs]
        if len(source_xref_keys) != len(set(source_xref_keys)):
            raise ValueError(
                f'Duplicate source xrefs in canonical gene {gene_key}'
            )

        ensembl_xref_ids = {
            xref['gene_id']
            for xref in xrefs
            if xref['database'] == 'Ensembl'
        }
        annotation_ensembl_id = str(
            gene.get('ensembl_id', '') or ''
        ).strip()
        if (
            ensembl_xref_ids
            and annotation_ensembl_id
            and annotation_ensembl_id not in ensembl_xref_ids
        ):
            raise ValueError(
                'Annotation Ensembl identifier does not reconcile with '
                f'pathway xrefs for {pathway_id}: {gene_key}; '
                f'annotation={annotation_ensembl_id}, '
                f'pathway={sorted(ensembl_xref_ids)}'
            )

        node_occurrence_count = sum(
            cls._nonnegative_int(
                xref.get('node_occurrence_count', ''),
                'pathway_xrefs.node_occurrence_count', pathway_id,
            )
            for xref in xrefs
        )
        if node_occurrence_count != cls._nonnegative_int(
            gene.get('node_occurrence_count', ''),
            'node_occurrence_count', pathway_id,
        ):
            raise ValueError(
                f'Node occurrences do not reconcile for {pathway_id}: {gene_key}'
            )
        if len(xrefs) != cls._nonnegative_int(
            gene.get('unique_xref_count', ''),
            'unique_xref_count', pathway_id,
        ):
            raise ValueError(
                f'Unique xrefs do not reconcile for {pathway_id}: {gene_key}'
            )

        primary_xref = xrefs[0]
        expected_primary = (
            primary_xref['database'], primary_xref['gene_id'],
            primary_xref['gene_key'],
        )
        supplied_primary = (
            str(gene.get('pathway_gene_database', '') or '').strip(),
            str(gene.get('pathway_gene_id', '') or '').strip(),
            str(gene.get('gene_key', '') or '').strip(),
        )
        if supplied_primary != expected_primary:
            raise ValueError(
                f'Primary pathway xref does not reconcile for {pathway_id}: '
                f'{gene_key}'
            )
        expected_databases = tuple(sorted({
            xref['database'] for xref in xrefs
        }))
        if cls._normalized_string_list(
            gene.get('source_databases', ())
        ) != expected_databases:
            raise ValueError(
                f'Source databases do not reconcile for {pathway_id}: {gene_key}'
            )
        canonicalization_complete = cls._strict_bool(
            gene.get('canonicalization_complete', ''),
            'canonicalization_complete', pathway_id,
        )

        guides = cls._guide_entries(gene.get('guide_targets', ''), pathway_id)
        has_guides = cls._strict_bool(
            gene.get('has_guides', ''), 'has_guides', pathway_id
        )
        guide_count = cls._nonnegative_int(
            gene.get('guide_count', ''), 'guide_count', pathway_id
        )
        if has_guides != bool(guides) or guide_count != len(guides):
            raise ValueError(
                f'Guide counts do not reconcile for {pathway_id}: {gene_key}'
            )
        if not matched and (has_guides or guide_count or guides):
            raise ValueError(
                f'Unmatched canonical gene has guides for {pathway_id}: {gene_key}'
            )
        return {
            'gene_key': gene_key,
            'source_xref_keys': source_xref_keys,
            'source_databases': expected_databases,
            'node_occurrence_count': node_occurrence_count,
            'canonicalization_complete': canonicalization_complete,
            'has_guides': has_guides,
        }

    @classmethod
    def _mapping_signature(cls, gene, pathway_id):
        """Return the exact mapping payload represented by a matched gene."""
        xrefs = cls._pathway_xrefs(gene)
        return (
            str(gene.get('pathway_gene_database', '') or ''),
            str(gene.get('pathway_gene_id', '') or ''),
            str(gene.get('pathway_gene_name', '') or ''),
            str(gene.get('gene_key', '') or ''),
            str(gene.get('pw_kfi_gene_key', '') or ''),
            cls._freeze_artifact_value(xrefs),
            cls._normalized_string_list(gene.get('source_databases', ())),
            cls._nonnegative_int(
                gene.get('node_occurrence_count', ''),
                'node_occurrence_count', pathway_id,
            ),
            cls._nonnegative_int(
                gene.get('unique_xref_count', ''),
                'unique_xref_count', pathway_id,
            ),
            str(gene.get('canonical_subject_id', '') or ''),
            cls._normalized_string_list(
                gene.get('canonicalization_methods', ())
            ),
            cls._strict_bool(
                gene.get('canonicalization_complete', ''),
                'canonicalization_complete', pathway_id,
            ),
            str(gene.get('annotation_source_database', '') or ''),
            str(gene.get('matched_annotation_key', '') or ''),
            str(gene.get('matched_annotation_id', '') or ''),
            str(gene.get('resolved_annotation_subject_id', '') or ''),
            str(gene.get('match_method', '') or ''),
            str(gene.get('mapping_source', '') or ''),
            str(gene.get('sgd_id', '') or ''),
            str(gene.get('ensembl_id', '') or ''),
            str(
                gene.get('gene_sgd_id', gene.get('sgd_id', '')) or ''
            ),
            str(gene.get('unique_id', '') or ''),
            str(gene.get('gene_name', '') or ''),
            cls._strict_bool(
                gene.get('has_guides', ''), 'has_guides', pathway_id
            ),
            cls._nonnegative_int(
                gene.get('guide_count', ''), 'guide_count', pathway_id
            ),
        )

    def get_pathway_diagram_data(self, pathway_id, pathways, unique_id_to_sgd=None):
        """
        Obtain diagrammatic data for a specific pathway.
        
        Args:
            pathway_id (str): Pathway ID
            pathways (list): Pathway data list
            unique_id_to_sgd (dict): Deprecated legacy mapping; not used for matching.
            
        Returns:
            dict: Pathway diagram data (SVG string, including metadata)
        """
        # 指定されたパスウェイを検索
        target_pathway = None
        for pathway in pathways:
            if pathway['pathway_id'] == pathway_id:
                target_pathway = pathway
                break
        
        if not target_pathway:
            self.logger.warning(f"Pathway not found.: {pathway_id}")
            return None
        
        # GPMLファイルパスを取得
        gpml_file_path = target_pathway.get('gpml_file_path', '')
        if not gpml_file_path or not os.path.exists(gpml_file_path):
            self.logger.warning(f"GPML file not found.: {gpml_file_path}")
            return None
        
        # GPMLファイルを解析
        pathway_data = self.parse_gpml_file(gpml_file_path)
        if not pathway_data:
            return None
        
        # ガイドRNA有無のマッピングを作成
        guide_status_map = {}
        unrenderable_gene_records = []
        matched_genes_value = target_pathway.get('matched_genes', '[]')
        
        try:
            matched_genes = self._parse_matched_genes(matched_genes_value)
            
            if not isinstance(matched_genes, list):
                raise ValueError('matched_genes must be a list')
            for index, gene_info in enumerate(matched_genes):
                if not isinstance(gene_info, dict):
                    raise ValueError(
                        f"matched_genes[{index}] must be an object"
                    )
                xrefs = self._pathway_xrefs(gene_info)
                if not xrefs:
                    raise ValueError(
                        f"matched_genes[{index}] lacks an exact pathway "
                        "database and identifier"
                    )

                unique_id = gene_info.get('unique_id', '')
                guide_targets = str(gene_info.get('guide_targets', '') or '')
                guides = [
                    item for item in guide_targets.split(',') if item.strip()
                ]
                sgd_id = str(gene_info.get('sgd_id', '') or '')
                if sgd_id.startswith('SGDID:'):
                    sgd_id = sgd_id[len('SGDID:'):]
                ensembl_id = str(gene_info.get('ensembl_id', '') or '')
                pw_kfi_gene_key = gene_info.get(
                    'pw_kfi_gene_key', gene_info.get('gene_key', '')
                )
                source_xref_keys = [xref['gene_key'] for xref in xrefs]
                for xref in xrefs:
                    database = xref['database']
                    gene_id = xref['gene_id']
                    gene_key = xref['gene_key']
                    xref_sgd_id = sgd_id or (
                        gene_id if database == 'SGD' else ''
                    )
                    xref_ensembl_id = ensembl_id or (
                        gene_id if database == 'Ensembl' else ''
                    )
                    previous = guide_status_map.get(gene_key)
                    if previous and previous.get(
                        'pw_kfi_gene_key'
                    ) != pw_kfi_gene_key:
                        raise ValueError(
                            f"Source xref {gene_key} maps to multiple PW-KFI genes"
                        )
                    guide_status_map[gene_key] = {
                            'gene_key': gene_key,
                            'pw_kfi_gene_key': pw_kfi_gene_key,
                            'source_xref_keys': source_xref_keys,
                            'pathway_gene_database': database,
                            'pathway_gene_id': gene_id,
                            'annotation_source_database': gene_info.get(
                                'annotation_source_database', 'Unspecified'
                            ),
                            'matched_annotation_key': gene_info.get(
                                'matched_annotation_key', ''
                            ),
                            'matched_annotation_id': gene_info.get(
                                'matched_annotation_id', ''
                            ),
                            'resolved_annotation_subject_id': gene_info.get(
                                'resolved_annotation_subject_id', ''
                            ),
                            'match_method': xref.get(
                                'match_method', gene_info.get('match_method', '')
                            ),
                            'mapping_source': xref.get(
                                'mapping_source', gene_info.get('mapping_source', '')
                            ),
                            'match_status': 'matched',
                            'failure_reason': '',
                            'sgd_id': xref_sgd_id,
                            'ensembl_id': xref_ensembl_id,
                            'has_guides': bool(guides),
                            'guide_count': len(guides),
                            'gene_name': gene_info.get('gene_name', ''),
                            'description': gene_info.get('description', ''),
                            'unique_id': unique_id,
                    }

            for field_name, match_status in (
                ('unmatched_genes', 'unmatched'),
                ('unmatchable_genes', 'not_evaluable'),
            ):
                unavailable_genes = self._parse_matched_genes(
                    target_pathway.get(field_name, '[]')
                )
                if not isinstance(unavailable_genes, list):
                    raise ValueError(f"{field_name} must be a list")
                for index, gene_info in enumerate(unavailable_genes):
                    if not isinstance(gene_info, dict):
                        raise ValueError(
                            f"{field_name}[{index}] must be an object"
                        )
                    xrefs = self._pathway_xrefs(gene_info)
                    if not xrefs:
                        unrenderable_gene_records.append({
                            'source_field': field_name,
                            'index': index,
                            'pathway_gene_database': '',
                            'pathway_gene_id': '',
                            'pathway_gene_name': gene_info.get('gene_name', ''),
                            'match_status': match_status,
                            'failure_reason': gene_info.get(
                                'reason', 'missing_exact_database_or_identifier'
                            ),
                        })
                        continue
                    for xref in xrefs:
                        database = xref['database']
                        gene_id = xref['gene_id']
                        gene_key = xref['gene_key']
                        if gene_key in guide_status_map:
                            continue
                        guide_status_map[gene_key] = {
                        'gene_key': gene_key,
                        'pw_kfi_gene_key': gene_info.get(
                            'pw_kfi_gene_key', gene_key
                        ),
                        'source_xref_keys': [
                            source_xref['gene_key'] for source_xref in xrefs
                        ],
                        'pathway_gene_database': database,
                        'pathway_gene_id': gene_id,
                        'annotation_source_database': '',
                        'matched_annotation_key': '',
                        'matched_annotation_id': '',
                        'resolved_annotation_subject_id': '',
                        'match_method': '',
                        'mapping_source': '',
                        'match_status': match_status,
                        'failure_reason': xref.get(
                            'reason', gene_info.get('reason', match_status)
                        ),
                        'sgd_id': '',
                        'ensembl_id': gene_id if database == 'Ensembl' else '',
                        'has_guides': False,
                        'guide_count': 0,
                        'gene_name': gene_info.get('gene_name', ''),
                        'description': '',
                        'unique_id': '',
                        }

            node_occurrence_counts = Counter(
                f"{node['database']}::{node['xref_id']}"
                for node in pathway_data['nodes']
                if node['type'] == 'GeneProduct'
                and node['database']
                and node['xref_id']
            )
            for gene_key, gene_status in guide_status_map.items():
                gene_status['node_occurrence_count'] = node_occurrence_counts.get(
                    gene_key, 0
                )
        
        except Exception as exc:
            raise ValueError(
                f"Cannot deserialize exact pathway gene data for {pathway_id}: {exc}"
            ) from exc
        
        # SVGを生成
        svg_content = self.generate_pathway_svg(pathway_data, guide_status_map)
        
        if svg_content:
            total_gene_nodes = len([
                node for node in pathway_data['nodes']
                if node['type'] == 'GeneProduct'
            ])
            evaluable_gene_nodes = self._nonnegative_int(
                target_pathway.get('evaluable_gene_node_count', ''),
                'evaluable_gene_node_count', pathway_id,
                default=sum(
                    value.get('node_occurrence_count', 0)
                    for value in guide_status_map.values()
                    if value.get('match_status') in {'matched', 'unmatched'}
                ),
            )
            total_genes = self._nonnegative_int(
                target_pathway.get('total_genes', ''),
                'total_genes', pathway_id,
                default=len({
                    value.get('pw_kfi_gene_key', key)
                    for key, value in guide_status_map.items()
                    if value.get('match_status') in {'matched', 'unmatched'}
                }),
            )
            unique_xref_count = self._nonnegative_int(
                target_pathway.get('unique_xref_count', ''),
                'unique_xref_count', pathway_id,
                default=sum(
                    value.get('match_status') in {'matched', 'unmatched'}
                    for value in guide_status_map.values()
                ),
            )
            genes_with_guides = self._nonnegative_int(
                target_pathway.get('genes_with_guides', ''),
                'genes_with_guides', pathway_id,
                default=len({
                    value.get('pw_kfi_gene_key', key)
                    for key, value in guide_status_map.items()
                    if value.get('match_status') == 'matched'
                    and value.get('has_guides')
                }),
            )
            gene_nodes_with_guides = sum(
                value.get('node_occurrence_count', 0)
                for value in guide_status_map.values()
                if value.get('match_status') == 'matched'
                and value.get('has_guides')
            )

            return {
                'pathway_id': pathway_id,
                'pathway_name': pathway_data['pathway_name'],
                'organism': pathway_data['organism'],
                'svg_content': svg_content,
                'board_width': pathway_data['board_width'],
                'board_height': pathway_data['board_height'],
                'guide_status_map': guide_status_map,
                'unrenderable_gene_records': unrenderable_gene_records,
                'total_genes': total_genes,
                'evaluable_gene_nodes': evaluable_gene_nodes,
                'unique_xref_count': unique_xref_count,
                'total_gene_nodes': total_gene_nodes,
                'genes_with_guides': genes_with_guides,
                'gene_nodes_with_guides': gene_nodes_with_guides,
            }
        
        return None

    def generate_pathway_guide_csv_data(
        self, pathways, unique_id_to_sgd=None, annotation_data=None, guide_data=None
    ):
        """Generate CSV data for pathway-specific guide designs."""
        pathway_guide_data = {}
        
        for pathway in pathways:
            pathway_id = pathway['pathway_id']
            pathway_name = pathway['pathway_name']
            matched_genes_value = pathway.get('matched_genes', '[]')
            
            # matched_genes列をJSONとして解析
            try:
                matched_genes = self._parse_matched_genes(matched_genes_value)
            except Exception as exc:
                raise ValueError(
                    f"Cannot deserialize matched genes for pathway {pathway_id}: "
                    f"{exc}"
                ) from exc
            if not isinstance(matched_genes, list):
                raise ValueError(
                    f"Matched genes must be a list for pathway {pathway_id}"
                )
            
            guides_data = []
            seen_guides = set()
            
            for gene_info in matched_genes:
                if not isinstance(gene_info, dict):
                    raise ValueError(
                        f"Matched gene entries must be objects for pathway "
                        f"{pathway_id}"
                    )
                unique_id = gene_info.get('unique_id', '')
                gene_name = gene_info.get('gene_name', '')
                description = gene_info.get('description', '')
                xrefs = self._pathway_xrefs(gene_info)
                if not xrefs:
                    raise ValueError(
                        f"Pathway {pathway_id} matched gene lacks an exact "
                        "pathway database and identifier"
                    )
                database = xrefs[0]['database']
                gene_id = xrefs[0]['gene_id']
                gene_key = xrefs[0]['gene_key']
                pw_kfi_gene_key = gene_info.get(
                    'pw_kfi_gene_key', gene_key
                )
                source_xref_keys = [xref['gene_key'] for xref in xrefs]
                sgd_id = str(gene_info.get('sgd_id', '') or '')
                if sgd_id.startswith('SGDID:'):
                    sgd_id = sgd_id[len('SGDID:'):]
                if not sgd_id:
                    sgd_id = next((
                        xref['gene_id'] for xref in xrefs
                        if xref['database'] == 'SGD'
                    ), '')
                ensembl_id = str(gene_info.get('ensembl_id', '') or '')
                if not ensembl_id:
                    ensembl_id = next((
                        xref['gene_id'] for xref in xrefs
                        if xref['database'] == 'Ensembl'
                    ), '')

                guide_entries = sorted({
                    item.strip()
                    for item in str(
                        gene_info.get('guide_targets', '') or ''
                    ).split(',')
                    if item.strip()
                })
                for guide_entry in guide_entries:
                    guide_identity = (pw_kfi_gene_key, guide_entry)
                    if guide_identity in seen_guides:
                        continue
                    seen_guides.add(guide_identity)
                    parts = guide_entry.split('/')
                    if len(parts) >= 6:
                        guide_info = {
                                            'pathway_id': pathway_id,
                                            'pathway_name': pathway_name,
                                            'pathway_gene_database': database,
                                            'pathway_gene_id': gene_id,
                                            'gene_key': gene_key,
                                            'pw_kfi_gene_key': pw_kfi_gene_key,
                                            'source_xref_keys': source_xref_keys,
                                            'source_databases': sorted({
                                                xref['database'] for xref in xrefs
                                            }),
                                            'annotation_source_database': gene_info.get(
                                                'annotation_source_database',
                                                'Unspecified',
                                            ),
                                            'matched_annotation_key': gene_info.get(
                                                'matched_annotation_key', ''
                                            ),
                                            'matched_annotation_id': gene_info.get(
                                                'matched_annotation_id', ''
                                            ),
                                            'resolved_annotation_subject_id': gene_info.get(
                                                'resolved_annotation_subject_id', ''
                                            ),
                                            'match_method': gene_info.get(
                                                'match_method', ''
                                            ),
                                            'sgd_id': sgd_id,
                                            'ensembl_id': ensembl_id,
                                            'unique_id': unique_id,
                                            'gene_name': gene_name,
                                            'description': description,
                                            'guide_unique_id': parts[0],
                                            'start_pos': parts[1],
                                            'end_pos': parts[2],
                                            'guide_sequence': parts[3],
                                            'strand': parts[7] if len(parts) > 7 else '',
                                            'MM0': parts[4],
                                            'MM1-2': parts[5],
                                            'off_target_match_score': parts[6] if len(parts) > 6 else ''
                        }
                        guides_data.append(guide_info)
            
            if guides_data:
                pathway_guide_data[pathway_id] = {
                    'pathway_name': pathway_name,
                    'guides': guides_data
                }
        
        self.logger.info(f"Pathway-specific guide CSV data generation complete.: {len(pathway_guide_data)} Pathway")
        return pathway_guide_data

    def generate_pathway_json_summary(self, pathways, pathway_guide_data=None, unique_id_to_sgd=None):
        """Generate a pathway summary JSON (including SVG data)."""
        # 統計情報計算
        total_pathways = len(pathways)
        fully_targetable = len([p for p in pathways if p['targetable_type'] == 'Fully targetable'])
        partially_targetable = len([p for p in pathways if p['targetable_type'] == 'Partially targetable'])
        not_targetable = len([p for p in pathways if p['targetable_type'] == 'Not targetable'])
        
        # パーセンテージ計算
        if total_pathways > 0:
            fully_targetable_percent = (fully_targetable / total_pathways) * 100
            partially_targetable_percent = (partially_targetable / total_pathways) * 100
            not_targetable_percent = (not_targetable / total_pathways) * 100
        else:
            fully_targetable_percent = partially_targetable_percent = not_targetable_percent = 0
        
        # チャートデータ生成
        chart_data = {
            'targetability': {
                'labels': ['Fully Targetable', 'Partially Targetable', 'Not Targetable'],
                'data': [fully_targetable, partially_targetable, not_targetable],
                'percentages': [fully_targetable_percent, partially_targetable_percent, not_targetable_percent],
                'colors': ['#28a745', '#ffc107', '#dc3545']
            }
        }
        
        # パスウェイSVGデータを事前生成
        pathway_svgs = {}
        svg_generation_stats = {
            'success': 0,
            'failed': 0,
            'no_gpml': 0,
            'failure_reasons': {},
        }
        
        self.logger.info(f"Started generating pathway SVG data.: {total_pathways} pathway")
        
        for i, pathway in enumerate(pathways):
            pathway_id = pathway['pathway_id']
            
            try:
                # 進捗表示
                if (i + 1) % 10 == 0 or i == 0:
                    self.logger.info(f"Pathway SVG generation progress: {i + 1}/{total_pathways}")
                
                # パスウェイ図データを生成
                diagram_data = self.get_pathway_diagram_data(pathway_id, pathways)

                if diagram_data:
                    pathway_svgs[pathway_id] = {
                        'pathway_name': diagram_data['pathway_name'],
                        'organism': diagram_data['organism'],
                        'svg_content': diagram_data['svg_content'],
                        'board_width': diagram_data['board_width'],
                        'board_height': diagram_data['board_height'],
                        'total_genes': diagram_data['total_genes'],
                        'evaluable_gene_nodes': diagram_data[
                            'evaluable_gene_nodes'
                        ],
                        'unique_xref_count': diagram_data[
                            'unique_xref_count'
                        ],
                        'total_gene_nodes': diagram_data['total_gene_nodes'],
                        'genes_with_guides': diagram_data['genes_with_guides'],
                        'gene_nodes_with_guides': diagram_data[
                            'gene_nodes_with_guides'
                        ],
                        'guide_status_map': diagram_data['guide_status_map'],
                        'unrenderable_gene_records': diagram_data[
                            'unrenderable_gene_records'
                        ],
                    }
                    svg_generation_stats['success'] += 1
                else:
                    # GPML file is missing or exact identifier data is unavailable.
                    gpml_file_path = pathway.get('gpml_file_path', '')
                    if not gpml_file_path or not os.path.exists(gpml_file_path):
                        svg_generation_stats['no_gpml'] += 1
                    else:
                        svg_generation_stats['failed'] += 1
                        svg_generation_stats['failure_reasons'][pathway_id] = (
                            'gpml_parse_or_exact_identifier_data_unavailable'
                        )

                    self.logger.debug(
                        f"Pathway {pathway_id} skipping SVG generation: "
                        "GPML or exact identifier issue"
                    )
                    
            except Exception as e:
                svg_generation_stats['failed'] += 1
                svg_generation_stats['failure_reasons'][pathway_id] = str(e)
                self.logger.warning(f"Pathway {pathway_id} SVG generation error: {e}")
        
        self.logger.info(f"Pathway SVG generation completed: Success={svg_generation_stats['success']}, "
                        f"Failure={svg_generation_stats['failed']}, No GPML={svg_generation_stats['no_gpml']}")
        
        summary_data = {
            'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'statistics': {
                'pathway_analysis': {
                    'total_pathways': total_pathways,
                    'fully_targetable_pathways': fully_targetable,
                    'partially_targetable_pathways': partially_targetable,
                    'not_targetable_pathways': not_targetable,
                    'fully_targetable_percent': fully_targetable_percent,
                    'partially_targetable_percent': partially_targetable_percent,
                    'not_targetable_percent': not_targetable_percent,
                    'chart_data': chart_data
                },
                'svg_generation': svg_generation_stats
            },
            'pathway_guide_data': pathway_guide_data if pathway_guide_data else {},
            'pathway_svgs': pathway_svgs
        }
        
        # Keep the machine-readable summary independently atomic.
        json_path = os.path.join(self.web_data_dir, 'pathway_summary.json')
        from .report_generator import ReportGenerator

        ReportGenerator._atomic_write_text(
            json_path,
            json.dumps(summary_data, indent=2, ensure_ascii=False) + '\n',
        )
        
        self.logger.info(f"Pathway summary JSON generation: {json_path}")
        return summary_data
    
    @staticmethod
    def _build_pathway_browser_data(pathways, summary_data):
        """Build the pathway-specific projection published for the browser."""
        pathway_table_data = []
        for pathway in pathways:
            pathway_table_data.append({
                'pathway_id': pathway['pathway_id'],
                'pathway_name': pathway['pathway_name'],
                'total_genes': int(pathway.get('total_genes', 0)),
                'total_gene_nodes': int(
                    pathway.get('total_gene_nodes', 0) or 0
                ),
                'evaluable_gene_nodes': int(
                    pathway.get('evaluable_gene_node_count', 0) or 0
                ),
                'unique_xref_count': int(
                    pathway.get('unique_xref_count', 0) or 0
                ),
                'genes_with_guides': int(pathway.get('genes_with_guides', 0)),
                'pw_kfi_score': float(pathway.get('pw_kfi_score', 0)),
                'targetable_type': pathway.get('targetable_type', 'Not targetable'),
                'coverage_rate': float(pathway.get('coverage_rate', 0))
            })

        # Raw mapping rows remain report sources rather than browser data.
        return {
            'summary': {
                'timestamp': summary_data.get('timestamp', ''),
                'statistics': summary_data.get('statistics', {}),
            },
            'pathways': pathway_table_data,
            'pathway_guide_data': summary_data.get('pathway_guide_data', {}),
            'pathway_svgs': summary_data.get('pathway_svgs', {})
        }

    def _create_legacy_pathway_html_template(
        self, pathways, summary_data, pathway_genes=None
    ):
        """Build pathway-only markup without embedding analysis payloads."""
        html_template = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>tWhEEL - Pathway Analysis Report</title>
    <style>
        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}
        
        html {{
            scroll-behavior: smooth;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            line-height: 1.5;
            background-color: #f8f9fa;
            color: #212529;
        }}
        
        /* Layout */
        .report-layout {{
            display: flex;
            min-height: 100vh;
        }}
        
        /* Sidebar */
        .sidebar {{
            width: 280px;
            background-color: #fff;
            border-right: 1px solid #dee2e6;
            padding: 0;
            position: fixed;
            height: 100vh;
            overflow-y: auto;
            z-index: 1000;
        }}
        
        .sidebar-header {{
            padding: 1.5rem;
            border-bottom: 1px solid #dee2e6;
            text-align: center;
        }}
        
        .sidebar-header img {{
            max-width: 150px;
            height: auto;
            margin-bottom: 0.5rem;
        }}
        
        .sidebar-header h1 {{
            font-size: 1.75rem;
            font-weight: 600;
            color: #5cb85c;
            margin: 0;
        }}
        
        .sidebar-header .subtitle {{
            color: #6c757d;
            font-size: 0.85rem;
            margin: 0.5rem 0;
            font-weight: 400;
        }}
        
        .sidebar-header .meta {{
            color: #6c757d;
            font-size: 0.8rem;
            margin: 0.25rem 0 0 0;
        }}
        
        .sidebar-nav {{
            padding: 0;
        }}
        
        .sidebar-nav ul {{
            list-style: none;
            padding: 0;
            margin: 0;
        }}
        
        .sidebar-nav li {{
            margin: 0;
        }}
        
        .sidebar-nav a {{
            display: block;
            padding: 0.75rem 1.5rem;
            color: #495057;
            text-decoration: none;
            font-size: 0.9rem;
            border-left: 3px solid transparent;
            transition: all 0.2s ease;
        }}
        
        .sidebar-nav a:hover {{
            background-color: #f8f9fa;
            color: #5cb85c;
            border-left-color: #dee2e6;
        }}
        
        .sidebar-nav a.active {{
            background-color: #e8f5e9;
            color: #5cb85c;
            border-left-color: #5cb85c;
            font-weight: 500;
        }}
        
        /* Main content */
        .main-content {{
            margin-left: 280px;
            flex: 1;
            padding: 0;
        }}
        
        .content-header {{
            background-color: #fff;
            border-bottom: 1px solid #dee2e6;
            padding: 1.5rem 2rem;
        }}
        
        .content-header h1 {{
            font-size: 2rem;
            font-weight: 600;
            color: #212529;
            margin: 0;
        }}
        
        .content-header .subtitle {{
            color: #6c757d;
            font-size: 1rem;
            margin: 0.25rem 0 0 0;
        }}
        
        .content-body {{
            padding: 2rem;
        }}
        
        /* Sections */
        .content-section {{
            background-color: #fff;
            border: 1px solid #dee2e6;
            border-radius: 0.375rem;
            margin-bottom: 2rem;
            overflow: hidden;
        }}
        
        .section-header {{
            background-color: #f8f9fa;
            border-bottom: 1px solid #dee2e6;
            padding: 1rem 1.5rem;
        }}
        
        .section-header h2 {{
            font-size: 1.25rem;
            font-weight: 600;
            color: #212529;
            margin: 0;
        }}
        
        .section-body {{
            padding: 1.5rem;
        }}
        
        /* Stats Grid */
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 1.5rem;
            margin-bottom: 1.5rem;
        }}
        
        .stat-card {{
            background: white;
            padding: 1.5rem;
            border-radius: 8px;
            border: 1px solid #dee2e6;
            text-align: center;
        }}
        
        .stat-card.total {{ border-left: 4px solid #007bff; }}
        .stat-card.fully {{ border-left: 4px solid #28a745; }}
        .stat-card.partially {{ border-left: 4px solid #ffc107; }}
        .stat-card.not {{ border-left: 4px solid #dc3545; }}
        
        .stat-number {{
            font-size: 2rem;
            font-weight: bold;
            margin-bottom: 0.5rem;
        }}
        
        .stat-label {{
            font-size: 0.9rem;
            color: #666;
        }}
        
        .chart-container {{
            background: white;
            border-radius: 8px;
            padding: 1.5rem;
            margin-bottom: 1rem;
            border: 1px solid #dee2e6;
        }}
        
        .table-container {{
            background: white;
            overflow: hidden;
        }}
        
        .table-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 1rem 1.5rem;
            flex-wrap: wrap;
            gap: 1rem;
            border-bottom: 1px solid #dee2e6;
        }}
        
        .table-header h3 {{
            margin: 0;
            font-size: 1rem;
            font-weight: 600;
            color: #212529;
        }}
        
        .search-filter {{
            display: flex;
            gap: 1rem;
            align-items: center;
            flex-wrap: wrap;
        }}
        
        input[type="text"], select {{
            padding: 0.5rem;
            border: 1px solid #ced4da;
            border-radius: 4px;
            font-size: 0.9rem;
        }}
        
        .btn {{
            background: #5cb85c;
            color: white;
            border: none;
            padding: 0.5rem 1rem;
            border-radius: 4px;
            cursor: pointer;
            text-decoration: none;
            display: inline-flex;
            align-items: center;
            gap: 0.5rem;
            font-size: 0.9rem;
        }}
        
        .btn:hover {{
            background: #4cae4c;
        }}
        
        .btn.secondary {{
            background: #6c757d;
        }}
        
        .btn.secondary:hover {{
            background: #545b62;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
        }}
        
        th, td {{
            text-align: left;
            padding: 0.75rem 1rem;
            border-bottom: 1px solid #dee2e6;
        }}
        
        th {{
            background-color: #f8f9fa;
            font-weight: 600;
            color: #495057;
            position: sticky;
            top: 0;
            z-index: 10;
        }}
        
        tr:hover {{
            background-color: #f8f9fa;
        }}
        
        .targetable-fully {{ color: #28a745; font-weight: 600; }}
        .targetable-partially {{ color: #856404; font-weight: 600; }}
        .targetable-not {{ color: #dc3545; font-weight: 600; }}
        
        .modal {{
            display: none;
            position: fixed;
            z-index: 1000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0,0,0,0.5);
        }}
        
        .modal-content {{
            background-color: white;
            margin: 5% auto;
            padding: 2rem;
            border-radius: 8px;
            width: 90%;
            max-width: 1000px;
            max-height: 80vh;
            overflow-y: auto;
        }}
        
        .close {{
            color: #aaa;
            float: right;
            font-size: 28px;
            font-weight: bold;
            cursor: pointer;
        }}
        
        .close:hover {{
            color: black;
        }}
        
        .pagination {{
            display: flex;
            justify-content: center;
            align-items: center;
            padding: 1rem;
            gap: 0.5rem;
        }}
        
        .page-btn {{
            padding: 0.5rem 0.75rem;
            border: 1px solid #dee2e6;
            background: white;
            cursor: pointer;
            border-radius: 4px;
        }}
        
        .page-btn:hover, .page-btn.active {{
            background: #5cb85c;
            color: white;
            border-color: #5cb85c;
        }}
        
        .loading {{
            text-align: center;
            padding: 2rem;
            color: #666;
        }}
        
        @media (max-width: 768px) {{
            .sidebar {{
                display: none;
            }}
            
            .main-content {{
                margin-left: 0;
            }}
            
            .stats-grid {{
                grid-template-columns: 1fr;
            }}
            
            .search-filter {{
                flex-direction: column;
                align-items: stretch;
            }}
            
            table {{
                font-size: 0.8rem;
            }}
            
            th, td {{
                padding: 0.5rem;
            }}
        }}
    </style>
</head>
<body>
    <div class="report-layout">
        <!-- Sidebar -->
        <aside class="sidebar">
            <div class="sidebar-header">
                <img src="static/img/logo.png" alt="tWhEEL Logo">
                <p class="subtitle">A pipeline to aggregate results from gene functional analysis, guide RNA design, and pathway-based knock-out feasibility assessment into a single report.</p>
                <p class="meta">Report generated on {html.escape(str(summary_data.get('timestamp', '')), quote=True)}</p>
            </div>
            <nav class="sidebar-nav">
                <ul>
                    <li><a href="#summary-stats" class="nav-link active">Summary Statistics</a></li>
                    <li><a href="#pathway-analysis" class="nav-link">Pathway Analysis</a></li>
                </ul>
            </nav>
        </aside>
        
        <!-- Main Content -->
        <main class="main-content">
            <header class="content-header">
                <h1>tWhEEL</h1>
                <p class="subtitle">transcript-aware Whole-genomic Evaluation for Editing Locus (tWhEEL)</p>
            </header>
            <div class="content-body">
                <section id="summary-stats" class="content-section">
                    <div class="section-header">
                        <h2>Summary Statistics</h2>
                    </div>
                    <div class="section-body">
        <div class="stats-grid">
            <div class="stat-card total">
                <div class="stat-number" id="total-pathways">-</div>
                <div class="stat-label">Total Pathways</div>
            </div>
            <div class="stat-card fully">
                <div class="stat-number" id="fully-targetable">-</div>
                <div class="stat-label">Fully Targetable</div>
            </div>
            <div class="stat-card partially">
                <div class="stat-number" id="partially-targetable">-</div>
                <div class="stat-label">Partially Targetable</div>
            </div>
            <div class="stat-card not">
                <div class="stat-number" id="not-targetable">-</div>
                <div class="stat-label">Not Targetable</div>
            </div>
        </div>
        
        <div class="chart-container">
            <h3>Pathway Targetability Distribution</h3>
            <div id="targetabilityChart" style="display: flex; justify-content: center; align-items: center;"></div>
        </div>
                    </div>
                </section>
                
                <section id="pathway-analysis" class="content-section">
                    <div class="section-header">
                        <h2>Pathway Analysis</h2>
                    </div>
                    <div class="section-body" style="padding: 0;">
        <div class="table-container">
            <div class="table-header">
                <h3>Pathway Analysis Results</h3>
                <div class="search-filter">
                    <input type="text" id="search-input" placeholder="Search pathways...">
                    <select id="targetable-filter">
                        <option value="">All Types</option>
                        <option value="Fully targetable">Fully Targetable</option>
                        <option value="Partially targetable">Partially Targetable</option>
                        <option value="Not targetable">Not Targetable</option>
                    </select>
                    <button class="btn" onclick="downloadPathwayCSV()">
                        📥 Download CSV
                    </button>
                </div>
            </div>
            
            <div id="table-loading" class="loading">Loading pathway data...</div>
            
            <table id="pathway-table" style="display: none;">
                <thead>
                    <tr>
                        <th onclick="sortTable(0)">Pathway ID</th>
                        <th onclick="sortTable(1)">Pathway Name</th>
                        <th onclick="sortTable(2)">Raw GPML Nodes</th>
                        <th onclick="sortTable(3)">Evaluable Nodes</th>
                        <th onclick="sortTable(4)">Unique Xrefs</th>
                        <th onclick="sortTable(5)">Unique PW-KFI Genes</th>
                        <th onclick="sortTable(6)">PW-KFI Genes with Guides</th>
                        <th onclick="sortTable(7)">PW-KFI Score</th>
                        <th onclick="sortTable(8)">Targetability</th>
                        <th onclick="sortTable(9)">Annotation Match Coverage</th>
                        <th>Actions</th>
                    </tr>
                </thead>
                <tbody id="pathway-table-body">
                </tbody>
            </table>
            
            <div class="pagination" id="pagination" style="display: none;">
            </div>
        </div>
                    </div>
                </section>
            </div>
        </main>
    </div>
    
    <!-- パスウェイ遺伝子詳細モーダル -->
    <div id="pathway-genes-modal" class="modal">
        <div class="modal-content" style="max-width: 90%; width: 90%;">
            <span class="close" onclick="closeModal('pathway-genes-modal')">&times;</span>
            <h3 id="modal-pathway-name">Pathway Genes</h3>
            <div class="table-container" style="max-height: 70vh; overflow-y: auto;">
                <table id="pathway-genes-table">
                    <thead>
                        <tr>
                            <th>Gene Name</th>
                            <th>Pathway Xref Database</th>
                            <th>Pathway Xref ID</th>
                            <th>Node Occurrences</th>
                            <th>Annotation Source</th>
                            <th>Description</th>
                            <th>Match Status</th>
                            <th>Failure Reason</th>
                            <th>Guide Status</th>
                            <th>Verified Link</th>
                        </tr>
                    </thead>
                    <tbody id="pathway-genes-table-body">
                    </tbody>
                </table>
                <div id="pathway-genes-empty" style="display: none; text-align: center; padding: 2rem; color: #666;">
                    This pathway does not contain genetic information.
                </div>
            </div>
        </div>
    </div>
    
    <!-- パスウェイ図表示モーダル -->
    <div id="pathway-diagram-modal" class="modal">
        <div class="modal-content" style="max-width: 95%; width: 95%;">
            <span class="close" onclick="closeModal('pathway-diagram-modal')">&times;</span>
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 1rem;">
                <h3 id="diagram-pathway-name">Pathway Diagram</h3>
                <div>
                    <button class="btn secondary" onclick="zoomDiagram(0.8)">🔍-</button>
                    <button class="btn secondary" onclick="zoomDiagram(1.25)">🔍+</button>
                    <button class="btn" onclick="resetZoom()">Reset</button>
                </div>
            </div>
            <div id="pathway-diagram-loading" class="loading">Loading pathway diagram......</div>
            <div id="pathway-diagram-container" style="display: none; border: 1px solid #ddd; overflow: auto; max-height: 70vh; background: white;">
                <div id="pathway-diagram-content"></div>
            </div>
            <div id="pathway-diagram-info" style="margin-top: 1rem; font-size: 0.9rem; color: #666;"></div>
        </div>
    </div>
    
    <!-- 遺伝子詳細情報モーダル -->
    <div id="gene-detail-modal" class="modal">
        <div class="modal-content" style="max-width: 90%; width: 90%;">
            <span class="close" onclick="closeModal('gene-detail-modal')">&times;</span>
            <h3 id="gene-detail-title">Gene Details</h3>
            <div id="gene-basic-info" style="background: #f8f9fa; padding: 1rem; margin-bottom: 1rem; border-radius: 4px;">
                <!-- Basic genetic information is displayed here. -->
            </div>
            <div class="table-container" style="max-height: 60vh; overflow-y: auto;">
                <h4>Guide RNA Details</h4>
                <table id="gene-guides-table">
                    <thead>
                        <tr>
                            <th>Guide Sequence</th>
                            <th>Start Position</th>
                            <th>End Position</th>
                            <th>MM0</th>
                            <th>MM1-2</th>
                            <th>Off-target Match Score</th>
                            <th>Strand</th>
                        </tr>
                    </thead>
                    <tbody id="gene-guides-table-body">
                    </tbody>
                </table>
                <div id="gene-guides-empty" style="display: none; text-align: center; padding: 2rem; color: #666;">
                    This gene does not contain guide RNA design information.
                </div>
            </div>
        </div>
    </div>
    
    <script>
        let pathwayData = [];
        let filteredData = [];
        let pathwayReportData = null;
        let currentPage = 1;
        const itemsPerPage = 25;
        let sortColumn = 0;
        let sortDirection = 'asc';

        function appendTextCell(row, value) {{
            const cell = row.insertCell();
            cell.textContent = String(value ?? '');
            return cell;
        }}

        function getPathwayReportData() {{
            if (pathwayReportData !== null) {{
                return pathwayReportData;
            }}
            const payload = window.__TWHEEL_PATHWAY_PAYLOAD__;
            if (!payload || typeof payload !== 'object' || Array.isArray(payload)) {{
                throw new Error('Pathway browser payload is unavailable.');
            }}
            if (
                !payload.summary
                || typeof payload.summary !== 'object'
                || Array.isArray(payload.summary)
                || !Array.isArray(payload.pathways)
                || !payload.pathway_guide_data
                || typeof payload.pathway_guide_data !== 'object'
                || Array.isArray(payload.pathway_guide_data)
                || !payload.pathway_svgs
                || typeof payload.pathway_svgs !== 'object'
                || Array.isArray(payload.pathway_svgs)
            ) {{
                throw new Error('Pathway browser payload has an invalid shape.');
            }}
            pathwayReportData = payload;
            return pathwayReportData;
        }}

        function validatedSVGMarkup(markup) {{
            const documentObject = new DOMParser().parseFromString(
                String(markup ?? ''), 'image/svg+xml'
            );
            const root = documentObject.documentElement;
            const allowedElements = new Set([
                'svg', 'defs', 'style', 'path', 'ellipse', 'rect', 'text', 'g'
            ]);
            const allowedAttributes = new Set([
                'width', 'height', 'viewbox', 'xmlns', 'd', 'class',
                'cx', 'cy', 'rx', 'ry', 'x', 'y', 'fill', 'stroke',
                'stroke-width', 'opacity', 'id', 'font-weight', 'font-size',
                'data-gene-database', 'data-gene-id', 'data-gene-key',
                'data-sgd-id', 'data-gene-name', 'title'
            ]);
            if (
                !root
                || root.localName !== 'svg'
                || documentObject.querySelector('parsererror')
            ) {{
                throw new Error('Invalid embedded pathway SVG');
            }}
            for (const element of documentObject.querySelectorAll('*')) {{
                if (!allowedElements.has(element.localName)) {{
                    throw new Error(
                        `Forbidden embedded SVG element: ${{element.localName}}`
                    );
                }}
                for (const attribute of [...element.attributes]) {{
                    const name = attribute.name.toLowerCase();
                    if (!allowedAttributes.has(name)) {{
                        throw new Error(
                            `Forbidden embedded SVG attribute: ${{name}}`
                        );
                    }}
                    if (
                        (name === 'fill' || name === 'stroke')
                        && !/^(?:#[0-9a-f]{{3,8}}|none|white)$/i.test(
                            attribute.value
                        )
                    ) {{
                        throw new Error(
                            `Forbidden embedded SVG paint value: ${{attribute.value}}`
                        );
                    }}
                }}
                if (
                    element.localName === 'style'
                    && /(?:url\\s*\\(|@import|expression\\s*\\(|javascript:)/i.test(
                        element.textContent || ''
                    )
                ) {{
                    throw new Error('Forbidden embedded SVG stylesheet');
                }}
            }}
            return new XMLSerializer().serializeToString(root);
        }}
        
        // 初期化
        document.addEventListener('DOMContentLoaded', function() {{
            loadPathwayData();
            setupEventListeners();
        }});
        
        function setupEventListeners() {{
            document.getElementById('search-input').addEventListener('input', filterTable);
            document.getElementById('targetable-filter').addEventListener('change', filterTable);
        }}
        
        function loadPathwayData() {{
            try {{
                const data = getPathwayReportData();
                
                // サマリーデータを表示
                loadPathwaySummary(data.summary);
                
                // パスウェイテーブルデータをロード
                pathwayData = data.pathways;
                filteredData = [...pathwayData];
                
                document.getElementById('table-loading').style.display = 'none';
                document.getElementById('pathway-table').style.display = 'table';
                document.getElementById('pagination').style.display = 'flex';
                
                renderTable();
                
            }} catch (error) {{
                console.error('Error loading data:', error);
                document.getElementById('table-loading').textContent = 'Failed to load data.: ' + error.message;
            }}
        }}
        
        function loadPathwaySummary(summaryData) {{
            try {{
                const stats = summaryData.statistics.pathway_analysis;
                document.getElementById('total-pathways').textContent = stats.total_pathways;
                document.getElementById('fully-targetable').textContent = stats.fully_targetable_pathways;
                document.getElementById('partially-targetable').textContent = stats.partially_targetable_pathways;
                document.getElementById('not-targetable').textContent = stats.not_targetable_pathways;
                
                // チャート生成
                createTargetabilityChart(stats.chart_data.targetability);
                
            }} catch (error) {{
                console.error('Summary data loading error:', error);
            }}
        }}
        
        function createTargetabilityChart(chartData) {{
            const chartContainer = document.getElementById('targetabilityChart');
            if (!chartContainer) return;
            
            // Create horizontal bar chart using SVG
            const chartSVG = createPathwayTargetabilityBarChart(chartData);
            
            chartContainer.innerHTML = '';
            chartContainer.appendChild(chartSVG);
            
            console.log('Pathway Targetability bar chart rendered');
        }}
        
        function createPathwayTargetabilityBarChart(chartData) {{
            const width = 900;
            const height = 300;
            const margin = {{ top: 80, right: 50, bottom: 70, left: 80 }};
            
            // Create SVG element
            const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
            svg.setAttribute('width', width);
            svg.setAttribute('height', height);
            svg.style.fontFamily = 'Arial, sans-serif';
            svg.style.fontSize = '12px';
            
            // Chart dimensions
            const chartWidth = width - margin.left - margin.right;
            const chartHeight = height - margin.top - margin.bottom;
            
            // Bar position and height
            const barY = margin.top + (chartHeight - 60) / 2;
            const barHeight = 60;
            
            // Calculate segment widths based on percentages
            const fullyPercent = chartData.percentages[0] || 0;
            const partiallyPercent = chartData.percentages[1] || 0;
            const notPercent = chartData.percentages[2] || 0;
            
            // Total count for statistics
            const totalPathways = chartData.data.reduce((sum, count) => sum + count, 0);
            
            // Segment data with proper positioning
            const segments = [
                {{
                    label: 'Fully Targetable',
                    percent: fullyPercent,
                    count: chartData.data[0],
                    color: '#28a745',
                    x: 0
                }},
                {{
                    label: 'Partially Targetable',
                    percent: partiallyPercent,
                    count: chartData.data[1],
                    color: '#ffc107',
                    x: (fullyPercent / 100) * chartWidth
                }},
                {{
                    label: 'Not Targetable',
                    percent: notPercent,
                    count: chartData.data[2],
                    color: '#dc3545',
                    x: ((fullyPercent + partiallyPercent) / 100) * chartWidth
                }}
            ];
            
            // Title
            const title = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            title.setAttribute('x', width / 2);
            title.setAttribute('y', 20);
            title.setAttribute('text-anchor', 'middle');
            title.setAttribute('fill', '#333');
            title.setAttribute('font-weight', 'bold');
            title.setAttribute('font-size', '16px');
            title.textContent = 'Pathway Targetability Distribution';
            svg.appendChild(title);
            
            // Legend - arranged in single row
            const legendStartX = 150;
            const legendY = 45;
            const legendSquareSize = 12;
            const legendSpacingX = 200;
            
            segments.forEach((segment, i) => {{
                const legendX = legendStartX + (i * legendSpacingX);
                
                // Legend square
                const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
                rect.setAttribute('x', legendX);
                rect.setAttribute('y', legendY - legendSquareSize / 2);
                rect.setAttribute('width', legendSquareSize);
                rect.setAttribute('height', legendSquareSize);
                rect.setAttribute('fill', segment.color);
                svg.appendChild(rect);
                
                // Legend text
                const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                text.setAttribute('x', legendX + legendSquareSize + 5);
                text.setAttribute('y', legendY + 3);
                text.setAttribute('fill', '#333');
                text.setAttribute('font-size', '11px');
                text.textContent = segment.label;
                svg.appendChild(text);
            }});
            
            // Y-axis label
            const yLabel = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            yLabel.setAttribute('x', 25);
            yLabel.setAttribute('y', barY + barHeight / 2);
            yLabel.setAttribute('text-anchor', 'middle');
            yLabel.setAttribute('fill', '#333');
            yLabel.setAttribute('font-weight', 'bold');
            yLabel.textContent = 'Pathways';
            svg.appendChild(yLabel);
            
            // Draw stacked bar segments
            segments.forEach(segment => {{
                if (segment.percent > 0) {{
                    const segmentWidth = (segment.percent / 100) * chartWidth;
                    
                    // Segment rectangle
                    const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
                    rect.setAttribute('x', margin.left + segment.x);
                    rect.setAttribute('y', barY);
                    rect.setAttribute('width', segmentWidth);
                    rect.setAttribute('height', barHeight);
                    rect.setAttribute('fill', segment.color);
                    rect.setAttribute('stroke', '#fff');
                    rect.setAttribute('stroke-width', '1');
                    svg.appendChild(rect);
                }}
            }});
            
            // Combined statistics text
            const statsText = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            statsText.setAttribute('x', margin.left + chartWidth / 2);
            statsText.setAttribute('y', barY + barHeight / 2);
            statsText.setAttribute('text-anchor', 'middle');
            statsText.setAttribute('fill', '#333');
            statsText.setAttribute('font-weight', 'bold');
            statsText.setAttribute('font-size', '12px');
            statsText.textContent = `F:${{segments[0].count}} P:${{segments[1].count}} N:${{segments[2].count}} n:${{totalPathways}}`;
            svg.appendChild(statsText);
            
            // X-axis scale
            const scaleY = barY + barHeight + 20;
            const scaleMarks = [0, 20, 40, 60, 80, 100];
            
            // X-axis line
            const xAxisLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
            xAxisLine.setAttribute('x1', margin.left);
            xAxisLine.setAttribute('y1', scaleY);
            xAxisLine.setAttribute('x2', margin.left + chartWidth);
            xAxisLine.setAttribute('y2', scaleY);
            xAxisLine.setAttribute('stroke', '#333');
            xAxisLine.setAttribute('stroke-width', '1');
            svg.appendChild(xAxisLine);
            
            // X-axis scale marks and labels
            scaleMarks.forEach(mark => {{
                const x = margin.left + (mark / 100) * chartWidth;
                
                // Scale mark
                const tick = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                tick.setAttribute('x1', x);
                tick.setAttribute('y1', scaleY);
                tick.setAttribute('x2', x);
                tick.setAttribute('y2', scaleY + 5);
                tick.setAttribute('stroke', '#333');
                tick.setAttribute('stroke-width', '1');
                svg.appendChild(tick);
                
                // Scale label
                const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                label.setAttribute('x', x);
                label.setAttribute('y', scaleY + 18);
                label.setAttribute('text-anchor', 'middle');
                label.setAttribute('fill', '#333');
                label.setAttribute('font-size', '10px');
                label.textContent = mark;
                svg.appendChild(label);
            }});
            
            // X-axis label
            const xLabel = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            xLabel.setAttribute('x', margin.left + chartWidth / 2);
            xLabel.setAttribute('y', height - 10);
            xLabel.setAttribute('text-anchor', 'middle');
            xLabel.setAttribute('fill', '#333');
            xLabel.setAttribute('font-weight', 'bold');
            xLabel.textContent = '% PW-KFI';
            svg.appendChild(xLabel);
            
            return svg;
        }}
        
        function renderTable() {{
            const tbody = document.getElementById('pathway-table-body');
            const startIndex = (currentPage - 1) * itemsPerPage;
            const endIndex = startIndex + itemsPerPage;
            const paginatedData = filteredData.slice(startIndex, endIndex);
            
            tbody.innerHTML = '';
            
            paginatedData.forEach(pathway => {{
                const row = tbody.insertRow();
                appendTextCell(row, pathway.pathway_id);
                appendTextCell(row, pathway.pathway_name);
                appendTextCell(row, pathway.total_gene_nodes);
                appendTextCell(row, pathway.evaluable_gene_nodes);
                appendTextCell(row, pathway.unique_xref_count);
                appendTextCell(row, pathway.total_genes);
                appendTextCell(row, pathway.genes_with_guides);
                appendTextCell(row, pathway.pw_kfi_score.toFixed(3));
                const targetableCell = appendTextCell(
                    row, pathway.targetable_type
                );
                const targetableClasses = {{
                    'Fully targetable': 'targetable-fully',
                    'Partially targetable': 'targetable-partially',
                    'Not targetable': 'targetable-not'
                }};
                targetableCell.className = (
                    targetableClasses[pathway.targetable_type]
                    || 'targetable-not'
                );
                appendTextCell(
                    row, `${{(pathway.coverage_rate * 100).toFixed(1)}}%`
                );
                const actionsCell = row.insertCell();
                const actions = [
                    ['btn secondary', 'View Genes', showPathwayGenes],
                    ['btn', '🔍 View Pathway', showPathwayDiagram],
                    ['btn', '📥 Guides', downloadPathwayGuides]
                ];
                actions.forEach(([className, label, handler]) => {{
                    const button = document.createElement('button');
                    button.className = className;
                    button.textContent = label;
                    button.addEventListener(
                        'click', () => handler(pathway.pathway_id)
                    );
                    actionsCell.appendChild(button);
                }});
            }});
            
            renderPagination();
        }}
        
        function renderPagination() {{
            const totalPages = Math.ceil(filteredData.length / itemsPerPage);
            const pagination = document.getElementById('pagination');
            pagination.innerHTML = '';
            
            // Previous button
            if (currentPage > 1) {{
                const prevBtn = document.createElement('button');
                prevBtn.className = 'page-btn';
                prevBtn.textContent = 'Previous';
                prevBtn.onclick = () => {{ currentPage--; renderTable(); }};
                pagination.appendChild(prevBtn);
            }}
            
            // Page numbers
            for (let i = 1; i <= totalPages; i++) {{
                if (i === currentPage || i === 1 || i === totalPages || Math.abs(i - currentPage) <= 2) {{
                    const pageBtn = document.createElement('button');
                    pageBtn.className = 'page-btn' + (i === currentPage ? ' active' : '');
                    pageBtn.textContent = i;
                    pageBtn.onclick = () => {{ currentPage = i; renderTable(); }};
                    pagination.appendChild(pageBtn);
                }}
            }}
            
            // Next button
            if (currentPage < totalPages) {{
                const nextBtn = document.createElement('button');
                nextBtn.className = 'page-btn';
                nextBtn.textContent = 'Next';
                nextBtn.onclick = () => {{ currentPage++; renderTable(); }};
                pagination.appendChild(nextBtn);
            }}
        }}
        
        function filterTable() {{
            const searchTerm = document.getElementById('search-input').value.toLowerCase();
            const targetableFilter = document.getElementById('targetable-filter').value;
            
            filteredData = pathwayData.filter(pathway => {{
                const matchesSearch = !searchTerm || 
                    pathway.pathway_id.toLowerCase().includes(searchTerm) ||
                    pathway.pathway_name.toLowerCase().includes(searchTerm);
                
                const matchesTargetable = !targetableFilter || pathway.targetable_type === targetableFilter;
                
                return matchesSearch && matchesTargetable;
            }});
            
            currentPage = 1;
            renderTable();
        }}
        
        function sortTable(column) {{
            if (sortColumn === column) {{
                sortDirection = sortDirection === 'asc' ? 'desc' : 'asc';
            }} else {{
                sortColumn = column;
                sortDirection = 'asc';
            }}
            
            const columns = [
                'pathway_id', 'pathway_name', 'total_gene_nodes',
                'evaluable_gene_nodes', 'unique_xref_count', 'total_genes',
                'genes_with_guides', 'pw_kfi_score', 'targetable_type',
                'coverage_rate'
            ];
            const sortKey = columns[column];
            
            filteredData.sort((a, b) => {{
                let valueA = a[sortKey];
                let valueB = b[sortKey];
                
                if (typeof valueA === 'string') {{
                    valueA = valueA.toLowerCase();
                    valueB = valueB.toLowerCase();
                }}
                
                if (sortDirection === 'asc') {{
                    return valueA < valueB ? -1 : valueA > valueB ? 1 : 0;
                }} else {{
                    return valueA > valueB ? -1 : valueA < valueB ? 1 : 0;
                }}
            }});
            
            renderTable();
        }}
        
        function showPathwayGenes(pathwayId) {{
            try {{
                // Retrieve pathway SVG information from the browser payload
                const data = getPathwayReportData();
                const pathwaySvgData = data.pathway_svgs[pathwayId];
                
                if (!pathwaySvgData || !pathwaySvgData.guide_status_map) {{
                    alert('Genetic information for this pathway could not be found.: ' + pathwayId);
                    return;
                }}
                
                // パスウェイ名を取得
                const pathway = data.pathways.find(p => p.pathway_id === pathwayId);
                const pathwayName = pathway ? pathway.pathway_name : pathwayId;
                
                // モーダルタイトルを設定
                document.getElementById('modal-pathway-name').textContent = pathwayName + ' (' + pathwayId + ') - Genes';
                
                // 遺伝子基本情報テーブルを構築
                const tbody = document.getElementById('pathway-genes-table-body');
                const emptyDiv = document.getElementById('pathway-genes-empty');
                tbody.innerHTML = '';
                
                // guide_status_mapから遺伝子情報を取得
                const guideStatusMap = pathwaySvgData.guide_status_map;
                const geneList = Object.keys(guideStatusMap);
                
                if (geneList.length === 0) {{
                    // 遺伝子情報がない場合
                    emptyDiv.style.display = 'block';
                    tbody.parentElement.parentElement.style.display = 'none';
                }} else {{
                    // 遺伝子基本情報を表示
                    emptyDiv.style.display = 'none';
                    tbody.parentElement.parentElement.style.display = 'block';
                    
                    geneList.forEach(geneKey => {{
                        const geneInfo = guideStatusMap[geneKey];
                        const hasGuides = geneInfo.has_guides;
                        const geneName = geneInfo.gene_name || 'Unknown';
                        const pathwayDatabase = geneInfo.pathway_gene_database || 'Unavailable';
                        const pathwayGeneId = geneInfo.pathway_gene_id || 'Unavailable';
                        const nodeOccurrenceCount = geneInfo.node_occurrence_count || 0;
                        const annotationSource = geneInfo.annotation_source_database || 'Unspecified';
                        const matchStatus = geneInfo.match_status || 'unavailable';
                        const failureReason = geneInfo.failure_reason || '';
                        
                        // Descriptionを取得（pathway_guide_dataから）
                        let description = geneInfo.description || 'N/A';
                        const pathwayGuideData = data.pathway_guide_data[pathwayId];
                        if (pathwayGuideData && pathwayGuideData.guides) {{
                            const guideInfo = pathwayGuideData.guides.find(g => g.gene_key === geneKey);
                            if (guideInfo && guideInfo.description) {{
                                description = guideInfo.description;
                            }}
                        }}
                        
                        const row = tbody.insertRow();
                        const nameCell = row.insertCell();
                        const nameElement = document.createElement('strong');
                        nameElement.textContent = String(geneName);
                        nameCell.appendChild(nameElement);
                        appendTextCell(row, pathwayDatabase);
                        const idCell = appendTextCell(row, pathwayGeneId);
                        idCell.style.fontFamily = 'monospace';
                        idCell.style.background = '#f8f9fa';
                        appendTextCell(row, nodeOccurrenceCount);
                        appendTextCell(row, annotationSource);
                        const descriptionCell = appendTextCell(
                            row,
                            description && description.length > 40
                                ? description.substring(0, 40) + '...'
                                : description
                        );
                        descriptionCell.title = String(description ?? '');
                        appendTextCell(row, matchStatus);
                        const reasonCell = appendTextCell(
                            row, failureReason || '—'
                        );
                        reasonCell.style.fontFamily = 'monospace';

                        const statusCell = row.insertCell();
                        const statusElement = document.createElement('span');
                        statusElement.style.fontWeight = 'bold';
                        if (matchStatus === 'matched') {{
                            statusElement.style.color = hasGuides
                                ? '#28a745' : '#dc3545';
                            statusElement.textContent = hasGuides
                                ? '✓ Available' : '✗ Not Available';
                            statusCell.appendChild(statusElement);
                            if (hasGuides) {{
                                statusCell.appendChild(document.createTextNode(
                                    ` (${{geneInfo.guide_count}} guides)`
                                ));
                            }}
                        }} else {{
                            statusElement.style.color = '#6c757d';
                            statusElement.textContent = 'Unavailable';
                            statusCell.appendChild(statusElement);
                        }}

                        const linkCell = row.insertCell();
                        if (geneInfo.sgd_id) {{
                            const linkButton = document.createElement('button');
                            linkButton.className = 'btn secondary';
                            linkButton.title = 'Open verified SGD annotation';
                            linkButton.textContent = '🔗 View SGD annotation';
                            linkButton.addEventListener(
                                'click', () => openSgdLink(geneInfo.sgd_id)
                            );
                            linkCell.appendChild(linkButton);
                        }} else {{
                            const unavailable = document.createElement('span');
                            unavailable.title = (
                                'No verified external URL is configured'
                            );
                            unavailable.textContent = 'Unavailable';
                            linkCell.appendChild(unavailable);
                        }}
                    }});
                    
                    // 統計情報を追加
                    const statsRow = tbody.insertRow();
                    statsRow.style.backgroundColor = '#e9ecef';
                    statsRow.style.fontWeight = 'bold';
                    const geneNodesWithGuides = geneList.reduce((total, geneKey) =>
                        total + (guideStatusMap[geneKey].has_guides
                            ? (guideStatusMap[geneKey].node_occurrence_count || 0)
                            : 0), 0
                    );
                    const totalGeneNodes = pathwaySvgData.total_gene_nodes || 0;
                    const evaluableGeneNodes = pathwaySvgData.evaluable_gene_nodes || 0;
                    const unrenderableRecords = pathwaySvgData.unrenderable_gene_records || [];
                    const unrenderableReasons = [...new Set(
                        unrenderableRecords.map(record => record.failure_reason || 'unavailable')
                    )].join(', ');
                    const totalGenes = pathwaySvgData.total_genes || 0;
                    const genesWithGuides = pathwaySvgData.genes_with_guides || 0;
                    const uniqueXrefCount = pathwaySvgData.unique_xref_count || 0;
                    const coverage = totalGenes > 0
                        ? ((genesWithGuides / totalGenes) * 100).toFixed(1)
                        : '0.0';
                    
                    const statsCell = statsRow.insertCell();
                    statsCell.colSpan = 10;
                    statsCell.style.textAlign = 'center';
                    statsCell.style.padding = '1rem';
                    const statsLabel = document.createElement('strong');
                    statsLabel.textContent = '統計情報:';
                    statsCell.append(
                        statsLabel,
                        document.createTextNode(
                            ` ${{totalGeneNodes}} total GPML gene nodes, `
                            + `${{evaluableGeneNodes}} evaluable nodes, `
                            + `${{uniqueXrefCount}} unique Database::ID xrefs, `
                            + `${{totalGenes}} unique PW-KFI genes, `
                            + `${{genesWithGuides}} PW-KFI genes with guide RNA targets, `
                            + `${{geneNodesWithGuides}} colored GPML nodes `
                            + `(${{coverage}}% PW-KFI gene coverage)`
                        )
                    );
                    if (unrenderableRecords.length > 0) {{
                        statsCell.append(
                            document.createElement('br'),
                            document.createTextNode(
                                `${{unrenderableRecords.length}} nodes cannot `
                                + `be keyed for diagram lookup: ${{unrenderableReasons}}`
                            )
                        );
                    }}
                }}
                
                // モーダルを表示
                document.getElementById('pathway-genes-modal').style.display = 'block';
                
            }} catch (error) {{
                console.error('Genetic data display error:', error);
                alert('Failed to display genetic data.: ' + error.message);
            }}
        }}
        
        function openSgdLink(sgdId) {{
            if (sgdId && sgdId !== 'N/A') {{
                const sgdUrl = (
                    'https://www.yeastgenome.org/locus/'
                    + encodeURIComponent(sgdId)
                );
                window.open(sgdUrl, '_blank', 'noopener,noreferrer');
            }} else {{
                alert('No valid SGD ID was found.');
            }}
        }}
        
        function copyToClipboard(text) {{
            // クリップボードにテキストをコピー
            if (navigator.clipboard && window.isSecureContext) {{
                navigator.clipboard.writeText(text).then(() => {{
                    // 成功時のフィードバック（簡易的）
                    const originalTitle = event.target.title;
                    event.target.title = 'Copied!';
                    setTimeout(() => {{
                        event.target.title = originalTitle;
                    }}, 1000);
                }}).catch(err => {{
                    console.error('Clipboard copy error:', err);
                }});
            }} else {{
                // フォールバック: 古いブラウザ対応
                const textArea = document.createElement('textarea');
                textArea.value = text;
                textArea.style.position = 'fixed';
                textArea.style.opacity = '0';
                document.body.appendChild(textArea);
                textArea.focus();
                textArea.select();
                try {{
                    document.execCommand('copy');
                    const originalTitle = event.target.title;
                    event.target.title = 'Copied!';
                    setTimeout(() => {{
                        event.target.title = originalTitle;
                    }}, 1000);
                }} catch (err) {{
                    console.error('Clipboard copy error:', err);
                }}
                document.body.removeChild(textArea);
            }}
        }}
        
        function closeModal(modalId) {{
            if (modalId) {{
                document.getElementById(modalId).style.display = 'none';
            }} else {{
                // 後方互換性のため
                document.getElementById('pathway-genes-modal').style.display = 'none';
            }}
        }}
        
        // パスウェイ図表示機能の変数
        let currentZoom = 1.0;
        let currentPathwayId = null;
        
        function showPathwayDiagram(pathwayId) {{
            try {{
                currentPathwayId = pathwayId;
                
                // パスウェイ名を取得
                const pathway = pathwayData.find(p => p.pathway_id === pathwayId);
                if (!pathway) {{
                    alert('Pathway not found.: ' + pathwayId);
                    return;
                }}
                
                // モーダルタイトルを設定
                document.getElementById('diagram-pathway-name').textContent = 
                    pathway.pathway_name + ' (' + pathwayId + ')';
                
                // ローディング表示
                document.getElementById('pathway-diagram-loading').style.display = 'block';
                document.getElementById('pathway-diagram-container').style.display = 'none';
                document.getElementById('pathway-diagram-info').textContent = '';
                
                // モーダルを表示
                document.getElementById('pathway-diagram-modal').style.display = 'block';
                
                // パスウェイ図を非同期で読み込み
                loadPathwayDiagram(pathwayId, pathway);
                
            }} catch (error) {{
                console.error('Pathway diagram display error:', error);
                alert('Failed to display the pathway diagram.: ' + error.message);
            }}
        }}
        
        function loadPathwayDiagram(pathwayId, pathway) {{
            // Read the requested pathway SVG from the browser payload.
            setTimeout(() => {{
                try {{
                    const data = getPathwayReportData();
                    const pathwaySvgData = data.pathway_svgs[pathwayId];
                    
                    if (!pathwaySvgData) {{
                        // SVGデータが見つからない場合のエラー表示
                        const loading = document.getElementById(
                            'pathway-diagram-loading'
                        );
                        loading.style.color = 'red';
                        loading.textContent = (
                            `The diagram data for this pathway (${{pathwayId}}) `
                            + 'could not be found. The GPML file may not exist, '
                            + 'or the parsing may have failed.'
                        );
                        return;
                    }}
                    
                    // 実際のSVGデータを表示
                    document.getElementById('pathway-diagram-content').innerHTML = (
                        validatedSVGMarkup(pathwaySvgData.svg_content)
                    );
                    document.getElementById('pathway-diagram-loading').style.display = 'none';
                    document.getElementById('pathway-diagram-container').style.display = 'block';
                    
                    // 実際の統計情報を表示
                    const totalGeneNodes = pathwaySvgData.total_gene_nodes || 0;
                    const evaluableGeneNodes = pathwaySvgData.evaluable_gene_nodes || 0;
                    const uniqueXrefCount = pathwaySvgData.unique_xref_count || 0;
                    const totalGenes = pathwaySvgData.total_genes || 0;
                    const genesWithGuides = pathwaySvgData.genes_with_guides || 0;
                    const geneNodesWithGuides = pathwaySvgData.gene_nodes_with_guides || 0;
                    const genesWithoutGuides = totalGenes - genesWithGuides;
                    const coverage = totalGenes > 0
                        ? ((genesWithGuides / totalGenes) * 100).toFixed(1)
                        : '0.0';
                    
                    document.getElementById('pathway-diagram-info').textContent = (
                        `Pathway: ${{pathwaySvgData.pathway_name}} | `
                        + `Total GPML gene nodes: ${{totalGeneNodes}} | `
                        + `Evaluable nodes: ${{evaluableGeneNodes}} | `
                        + `Unique Database::ID xrefs: ${{uniqueXrefCount}} | `
                        + `Unique PW-KFI genes: ${{totalGenes}} | `
                        + `PW-KFI genes with guides: ${{genesWithGuides}} | `
                        + `PW-KFI genes without guides: ${{genesWithoutGuides}} | `
                        + `Colored GPML nodes: ${{geneNodesWithGuides}} | `
                        + `PW-KFI gene coverage: ${{coverage}}%`
                    );
                    
                    // SVG要素にクリックイベントを追加
                    addGeneClickEvents();
                    
                    // ズームをリセット
                    currentZoom = 1.0;
                    updateDiagramZoom();
                    
                }} catch (error) {{
                    console.error('パスウェイ図表示エラー:', error);
                    const loading = document.getElementById(
                        'pathway-diagram-loading'
                    );
                    loading.style.color = 'red';
                    loading.textContent = (
                        'Failed to display the pathway diagram.: '
                        + error.message
                    );
                }}
            }}, 500); // 短い遅延で実際のデータを表示
        }}
        
        function addGeneClickEvents() {{
            const geneNodes = document.querySelectorAll('.gene-node');
            geneNodes.forEach(node => {{
                node.addEventListener('click', function() {{
                    const geneName = this.getAttribute('data-gene-name');
                    const geneDatabase = this.getAttribute('data-gene-database');
                    const geneId = this.getAttribute('data-gene-id');
                    const geneKey = this.getAttribute('data-gene-key');
                    const sgdId = this.getAttribute('data-sgd-id');
                    const hasGuides = this.classList.contains('gene-with-guides');
                    
                    // 詳細なモーダルを表示
                    showGeneDetails(geneName, geneDatabase, geneId, geneKey, sgdId, hasGuides);
                }});
            }});
        }}
        
        function showGeneDetails(geneName, geneDatabase, geneId, geneKey, sgdId, hasGuides) {{
            try {{
                // Read pathway guide data from the browser payload.
                const data = getPathwayReportData();
                
                // 現在表示中のパスウェイのガイドデータを取得
                const pathwayGuideData = data.pathway_guide_data[currentPathwayId];
                const pathwaySvgData = data.pathway_svgs[currentPathwayId];
                const geneStatus = pathwaySvgData && pathwaySvgData.guide_status_map
                    ? (pathwaySvgData.guide_status_map[geneKey] || {{}})
                    : {{}};
                
                // モーダルタイトルを設定
                document.getElementById('gene-detail-title').textContent = `${{geneName}} - Gene Details`;
                
                // 遺伝子基本情報を表示
                const basicInfoDiv = document.getElementById('gene-basic-info');
                const infoGrid = document.createElement('div');
                infoGrid.style.display = 'grid';
                infoGrid.style.gridTemplateColumns = (
                    'repeat(auto-fit, minmax(200px, 1fr))'
                );
                infoGrid.style.gap = '1rem';
                const addInfo = (label, value) => {{
                    const item = document.createElement('div');
                    const labelElement = document.createElement('strong');
                    labelElement.textContent = `${{label}}:`;
                    item.append(
                        labelElement,
                        document.createTextNode(` ${{String(value ?? '')}}`)
                    );
                    infoGrid.appendChild(item);
                }};
                addInfo('Gene Name', geneName || 'N/A');
                addInfo(
                    'Pathway Xref Database', geneDatabase || 'Unavailable'
                );
                addInfo('Pathway Xref ID', geneId || 'Unavailable');
                addInfo(
                    'PW-KFI Gene Key',
                    geneStatus.pw_kfi_gene_key || 'Unavailable'
                );
                addInfo(
                    'Source Xrefs',
                    (geneStatus.source_xref_keys || []).join(', ')
                    || 'Unavailable'
                );
                addInfo(
                    'Resolved SGD annotation ID', sgdId || 'Unavailable'
                );
                addInfo(
                    'Match Status', geneStatus.match_status || 'unavailable'
                );
                addInfo(
                    'Failure Reason', geneStatus.failure_reason || '—'
                );
                const guideStatusItem = document.createElement('div');
                const guideStatusLabel = document.createElement('strong');
                guideStatusLabel.textContent = 'Guide Status:';
                const guideStatusValue = document.createElement('span');
                guideStatusValue.style.color = hasGuides
                    ? '#28a745' : '#dc3545';
                guideStatusValue.style.fontWeight = 'bold';
                guideStatusValue.textContent = hasGuides
                    ? ' ✓ Available' : ' ✗ Not Available';
                guideStatusItem.append(guideStatusLabel, guideStatusValue);
                infoGrid.appendChild(guideStatusItem);
                basicInfoDiv.replaceChildren(infoGrid);
                
                // ガイドRNA詳細テーブルを構築
                const tbody = document.getElementById('gene-guides-table-body');
                const emptyDiv = document.getElementById('gene-guides-empty');
                const tableContainer = document.querySelector('#gene-guides-table').parentElement;
                
                tbody.innerHTML = '';
                
                // 該当遺伝子のガイドRNAを検索
                let geneGuides = [];
                if (pathwayGuideData && pathwayGuideData.guides) {{
                    geneGuides = pathwayGuideData.guides.filter(
                        guide => guide.gene_key === geneKey ||
                            (Array.isArray(guide.source_xref_keys) &&
                                guide.source_xref_keys.includes(geneKey))
                    );
                }}
                
                if (geneGuides.length === 0) {{
                    // ガイドRNAがない場合
                    tableContainer.style.display = 'none';
                    emptyDiv.style.display = 'block';
                }} else {{
                    // ガイドRNA詳細を表示
                    tableContainer.style.display = 'block';
                    emptyDiv.style.display = 'none';
                    
                    geneGuides.forEach(guide => {{
                        const row = tbody.insertRow();
                        const sequence = guide.guide_sequence || '';
                        const sequenceCell = appendTextCell(
                            row, sequence || 'N/A'
                        );
                        sequenceCell.style.fontFamily = 'monospace';
                        sequenceCell.style.fontSize = '0.9em';
                        sequenceCell.style.background = '#f8f9fa';
                        sequenceCell.style.padding = '0.3rem';
                        sequenceCell.style.cursor = 'pointer';
                        sequenceCell.title = `Click to copy: ${{sequence}}`;
                        sequenceCell.addEventListener(
                            'click', () => copySequenceToClipboard(
                                sequenceCell, sequence
                            )
                        );
                        appendTextCell(row, guide.start_pos ?? 'N/A');
                        appendTextCell(row, guide.end_pos ?? 'N/A');
                        appendTextCell(row, guide['MM0'] ?? 'N/A');
                        appendTextCell(row, guide['MM1-2'] ?? 'N/A');
                        appendTextCell(
                            row, guide['off_target_match_score'] ?? 'N/A'
                        );
                        appendTextCell(row, guide['strand'] || 'N/A');
                    }});
                    
                    // 統計情報を追加
                    const statsRow = tbody.insertRow();
                    statsRow.style.backgroundColor = '#e9ecef';
                    statsRow.style.fontWeight = 'bold';
                    const statsCell = statsRow.insertCell();
                    statsCell.colSpan = 7;
                    statsCell.style.textAlign = 'center';
                    statsCell.style.padding = '1rem';
                    const statsLabel = document.createElement('strong');
                    statsLabel.textContent = '統計:';
                    statsCell.append(
                        statsLabel,
                        document.createTextNode(
                            ` ${{geneGuides.length}} Guide RNA design completed`
                        )
                    );
                }}
                
                // モーダルを表示
                document.getElementById('gene-detail-modal').style.display = 'block';
                
            }} catch (error) {{
                console.error('Genetic detail display error:', error);
                alert('Failed to display detailed gene information.: ' + error.message);
            }}
        }}
        
        function copySequenceToClipboard(element, sequence) {{
            if (sequence && sequence !== 'N/A') {{
                navigator.clipboard.writeText(sequence).then(() => {{
                    // 視覚的フィードバック
                    const originalBg = element.style.backgroundColor;
                    const originalColor = element.style.color;
                    element.style.backgroundColor = '#28a745';
                    element.style.color = 'white';
                    
                    setTimeout(() => {{
                        element.style.backgroundColor = originalBg;
                        element.style.color = originalColor;
                    }}, 1000);
                    
                    console.log('Sequence copied:', sequence);
                }}).catch(err => {{
                    console.error('Clipboard copy failed:', err);
                    // フォールバック処理
                    const textArea = document.createElement('textarea');
                    textArea.value = sequence;
                    document.body.appendChild(textArea);
                    textArea.select();
                    document.execCommand('copy');
                    document.body.removeChild(textArea);
                    
                    // 視覚的フィードバック
                    const originalBg = element.style.backgroundColor;
                    element.style.backgroundColor = '#28a745';
                    element.style.color = 'white';
                    setTimeout(() => {{
                        element.style.backgroundColor = originalBg;
                        element.style.color = '';
                    }}, 1000);
                }});
            }}
        }}
        
        function zoomDiagram(factor) {{
            if (!currentPathwayId) return;
            
            currentZoom *= factor;
            // ズーム範囲を制限
            currentZoom = Math.max(0.2, Math.min(5.0, currentZoom));
            
            updateDiagramZoom();
        }}
        
        function resetZoom() {{
            if (!currentPathwayId) return;
            
            currentZoom = 1.0;
            updateDiagramZoom();
        }}
        
        function updateDiagramZoom() {{
            const svgElement = document.querySelector('#pathway-diagram-content svg');
            if (svgElement) {{
                svgElement.style.transform = `scale(${{currentZoom}})`;
                svgElement.style.transformOrigin = '0 0';
                
                // コンテナのサイズを調整
                const container = document.getElementById('pathway-diagram-content');
                const originalWidth = parseFloat(svgElement.getAttribute('width')) || 600;
                const originalHeight = parseFloat(svgElement.getAttribute('height')) || 500;
                
                container.style.width = (originalWidth * currentZoom) + 'px';
                container.style.height = (originalHeight * currentZoom) + 'px';
            }}
        }}
        
        function downloadPathwayCSV() {{
            // CSV ダウンロード機能
            const csvContent = generateCSVContent();
            const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
            const link = document.createElement('a');
            const url = URL.createObjectURL(blob);
            link.setAttribute('href', url);
            link.setAttribute('download', 'pathway_analysis.csv');
            link.style.visibility = 'hidden';
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}

        function escapeCSVField(value) {{
            if (value === null || value === undefined) return '';
            const normalized = String(value).replace(/\\r\\n|\\r|\\n/g, ' ');
            if (!/[",]/.test(normalized)) return normalized;
            return '"' + normalized.replace(/"/g, '""') + '"';
        }}

        function serializeCSVRow(values) {{
            return values.map(escapeCSVField).join(',');
        }}
        
        function generateCSVContent() {{
            const headers = [
                'Pathway ID', 'Pathway Name', 'Raw GPML Nodes',
                'Evaluable Nodes', 'Unique Xrefs', 'Unique PW-KFI Genes',
                'PW-KFI Genes with Guides', 'PW-KFI Score', 'Targetability',
                'Annotation Match Coverage'
            ];
            const rows = filteredData.map(pathway => [
                pathway.pathway_id,
                pathway.pathway_name,
                pathway.total_gene_nodes,
                pathway.evaluable_gene_nodes,
                pathway.unique_xref_count,
                pathway.total_genes,
                pathway.genes_with_guides,
                pathway.pw_kfi_score.toFixed(4),
                pathway.targetable_type,
                (pathway.coverage_rate * 100).toFixed(2) + '%'
            ]);
            
            return [
                headers.join(','),
                ...rows.map(row => serializeCSVRow(row))
            ].join('\\n');
        }}
        
        function downloadPathwayGuides(pathwayId) {{
            try {{
                // Read pathway guide data from the browser payload.
                const data = getPathwayReportData();
                const pathwayGuideData = data.pathway_guide_data[pathwayId];
                
                if (!pathwayGuideData || !pathwayGuideData.guides || pathwayGuideData.guides.length === 0) {{
                    alert('Guide RNA information for this pathway could not be found.: ' + pathwayId);
                    return;
                }}
                
                // パスウェイガイドCSVコンテンツを生成
                const csvContent = generatePathwayGuidesCSVContent(pathwayGuideData);
                
                // CSVファイルダウンロード
                const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
                const link = document.createElement('a');
                const url = URL.createObjectURL(blob);
                
                // ファイル名を生成（パスウェイIDを含む）
                const fileName = `pathway_${{pathwayId}}_guides.csv`;
                
                link.setAttribute('href', url);
                link.setAttribute('download', fileName);
                link.style.visibility = 'hidden';
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
                
                // 成功通知
                console.log(`Guide RNA information for pathway ${{pathwayId}} has been downloaded.: ${{fileName}}`);
                
            }} catch (error) {{
                console.error('Guide RNA download error:', error);
                alert('Guide RNA download failed.: ' + error.message);
            }}
        }}
        
        function generatePathwayGuidesCSVContent(pathwayGuideData) {{
            // CSVヘッダー
            const headers = [
                'Pathway ID',
                'Pathway Name',
                'Pathway Xref Database',
                'Pathway Xref ID',
                'Gene Key',
                'PW-KFI Gene Key',
                'Source Xref Keys',
                'Source Databases',
                'Annotation Source Database',
                'Matched Annotation Key',
                'Matched Annotation ID',
                'Resolved Annotation Subject ID',
                'Match Method',
                'SGD ID',
                'Ensembl ID',
                'Gene Name',
                'Unique ID',
                'Description',
                'Guide Unique ID',
                'Start Position',
                'End Position',
                'Guide Sequence',
                'MM0',
                'MM1-2',
                'off_target_match_score',
                'strand',
            ];
            
            // データ行を生成
            const rows = pathwayGuideData.guides.map(guide => [
                guide.pathway_id || '',
                guide.pathway_name || '',
                guide.pathway_gene_database || '',
                guide.pathway_gene_id || '',
                guide.gene_key || '',
                guide.pw_kfi_gene_key || '',
                (guide.source_xref_keys || []).join(';'),
                (guide.source_databases || []).join(';'),
                guide.annotation_source_database || '',
                guide.matched_annotation_key || '',
                guide.matched_annotation_id || '',
                guide.resolved_annotation_subject_id || '',
                guide.match_method || '',
                guide.sgd_id || '',
                guide.ensembl_id || '',
                guide.gene_name || '',
                guide.unique_id || '',
                guide.description || '',
                guide.guide_unique_id || '',
                guide.start_pos || '',
                guide.end_pos || '',
                guide.guide_sequence || '',
                guide['MM0'] == null ? '' : guide['MM0'],
                guide['MM1-2'] == null ? '' : guide['MM1-2'],
                guide.off_target_match_score || '',
                guide.strand || '',
            ]);
            
            // CSV形式に変換
            return [
                headers.join(','),
                ...rows.map(row => serializeCSVRow(row))
            ].join('\\n');
        }}
    </script>
</body>
</html>'''
        
        return html_template

    def _normalized_shared_report_config(self):
        """Return a report-safe config without mutating the pipeline config."""
        config = copy.deepcopy(self.config)
        config.setdefault('input', {})
        config.setdefault('crispr', {})
        config.setdefault('busco', {})
        config.setdefault('pathway_analysis', {})
        return config

    def _collect_shared_report_context(self):
        """Collect the data used by the common GO/PW report sections."""
        from .report_generator import ReportGenerator

        generator = ReportGenerator(self._normalized_shared_report_config())
        return generator, generator.collect_all_results()

    @staticmethod
    def _pathway_shared_report_results(common_results):
        """Remove GO-only inputs from the report data shared with PW output."""
        base_results = dict(common_results)
        base_results['data'] = dict(common_results.get('data', {}))
        base_results['report_data'] = dict(
            common_results.get('report_data', {})
        )
        base_results['statistics'] = dict(
            common_results.get('statistics', {})
        )
        for key in (
            'go_kfi_stat', 'go_kfi_table', 'go_reverse_lookup',
            'go_parent_counts', 'go_parent_counts_image',
            'go_parent_counts_treemap_image',
        ):
            base_results['data'].pop(key, None)
            base_results['report_data'].pop(key, None)
        base_results['statistics'].pop('go_analysis', None)
        return base_results

    @staticmethod
    def _replace_report_block(content, start_marker, end_marker, replacement):
        """Replace one marker-delimited report block and fail on drift."""
        start = content.find(start_marker)
        if start < 0:
            raise ValueError(f'Report start marker not found: {start_marker}')
        end = content.find(end_marker, start)
        if end < 0:
            raise ValueError(f'Report end marker not found: {end_marker}')
        return content[:start] + replacement + content[end:]

    @staticmethod
    def _pathway_assessment_html(summary_data):
        """Return the PW-KFI assessment section in the shared report style."""
        stats = summary_data.get('statistics', {}).get(
            'pathway_analysis', {}
        )

        def integer(key):
            try:
                return int(stats.get(key, 0) or 0)
            except (TypeError, ValueError):
                return 0

        return f'''                <!-- PW-KFI Assessment Section -->
                <section id="pw-kfi-assessment" class="content-section">
                    <div class="section-header">
                        <h2>PW-KFI Assessment</h2>
                    </div>
                    <div class="section-body">
                        <p>Pathway-based knock-out feasibility assessment (PW-KFI assessment).</p>

                        <div class="stats-grid">
                            <div class="stat-card">
                                <div class="stat-value stat-value-total-go" id="total-pathways">{integer('total_pathways'):,}</div>
                                <div class="stat-label">Total Pathways</div>
                            </div>
                            <div class="stat-card">
                                <div class="stat-value stat-value-fully-targetable" id="fully-targetable">{integer('fully_targetable_pathways'):,}</div>
                                <div class="stat-label">Fully Targetable Pathways</div>
                            </div>
                            <div class="stat-card">
                                <div class="stat-value stat-value-partially-targetable" id="partially-targetable">{integer('partially_targetable_pathways'):,}</div>
                                <div class="stat-label">Partially Targetable Pathways</div>
                            </div>
                            <div class="stat-card">
                                <div class="stat-value stat-value-not-targetable" id="not-targetable">{integer('not_targetable_pathways'):,}</div>
                                <div class="stat-label">Not Targetable Pathways</div>
                            </div>
                        </div>

                        <div class="pw-chart-block">
                            <h3>Pathway Targetability Distribution</h3>
                            <div class="chart-container" id="targetabilityChart"></div>
                        </div>

                        <div class="pw-results-block">
                            <div class="pw-results-header">
                                <h3>Pathway Analysis Results</h3>
                                <div class="filter-group">
                                    <input type="text" class="search-box" id="search-input" placeholder="Search pathways...">
                                    <select class="filter-select" id="targetable-filter">
                                        <option value="">All PW-KFI Classes</option>
                                        <option value="Fully targetable">Fully Targetable</option>
                                        <option value="Partially targetable">Partially Targetable</option>
                                        <option value="Not targetable">Not Targetable</option>
                                    </select>
                                    <button class="btn btn-primary" onclick="downloadPathwayCSV()" title="Download Pathway Analysis CSV">
                                        📥 Download CSV
                                    </button>
                                </div>
                            </div>

                            <div id="table-loading" class="loading">Loading pathway data...</div>
                            <div class="table-responsive">
                                <table class="data-table" id="pathway-table" style="display: none;">
                                    <thead>
                                        <tr>
                                            <th onclick="sortTable(0)">Pathway ID</th>
                                            <th onclick="sortTable(1)">Pathway Name</th>
                                            <th onclick="sortTable(2)">Raw GPML Nodes</th>
                                            <th onclick="sortTable(3)">Evaluable Nodes</th>
                                            <th onclick="sortTable(4)">Unique Xrefs</th>
                                            <th onclick="sortTable(5)">Unique PW-KFI Genes</th>
                                            <th onclick="sortTable(6)">PW-KFI Genes with Guides</th>
                                            <th onclick="sortTable(7)">PW-KFI Score</th>
                                            <th onclick="sortTable(8)">Targetability</th>
                                            <th onclick="sortTable(9)">Annotation Match Coverage</th>
                                            <th>Actions</th>
                                        </tr>
                                    </thead>
                                    <tbody id="pathway-table-body"></tbody>
                                </table>
                            </div>
                            <div class="pagination-controls" id="pagination" style="display: none;"></div>
                        </div>
                    </div>
                </section>

'''

    def _general_statistics_html(self, common_results, summary_data):
        """Return pathway-aware general statistics in the shared layout."""
        stats = common_results.get('statistics', {})
        pathway_stats = (
            summary_data.get('statistics', {}).get('pathway_analysis', {})
        )

        def integer(group, key):
            try:
                return int(stats.get(group, {}).get(key, 0) or 0)
            except (TypeError, ValueError):
                return 0

        total_pathways = int(pathway_stats.get('total_pathways', 0) or 0)
        fully_percent = float(
            pathway_stats.get('fully_targetable_percent', 0) or 0
        )
        return f'''                <!-- General Statistics Section -->
                <section id="general-stats" class="content-section">
                    <div class="section-header">
                        <h2>General Statistics</h2>
                    </div>
                    <div class="section-body">
                        <div class="alert alert-info">
                            <strong>Welcome!</strong> Not sure where to start?
                            <a href="#pw-kfi-assessment">Jump to PW-KFI Assessment</a> section to see detailed results.
                        </div>
                        <div class="table-responsive">
                            <table class="data-table">
                                <thead>
                                    <tr>
                                        <th></th>
                                        <th>Total transcripts</th>
                                        <th>Protein Annotated</th>
                                        <th>Total Annotated Genes</th>
                                        <th>Multi (Redundant)</th>
                                        <th>Single (Non-Redundant)</th>
                                        <th>Pathways</th>
                                        <th>Fully Targetable Pathways [%]</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    <tr>
                                        <td>Counts</td>
                                        <td>{integer('gene_prediction', 'transcript_count'):,}</td>
                                        <td>{integer('gene_prediction', 'protein_count'):,}</td>
                                        <td>{integer('annotation', 'total_genes'):,}</td>
                                        <td>{integer('annotation', 'multi_genes'):,}</td>
                                        <td>{integer('annotation', 'single_genes'):,}</td>
                                        <td>{total_pathways:,}</td>
                                        <td>{fully_percent:.1f}%</td>
                                    </tr>
                                </tbody>
                            </table>
                        </div>
                    </div>
                </section>

'''

    def _analysis_configuration_html(self, common_results):
        """Return PW-specific configuration cards in the shared layout."""
        config = common_results.get('config', {})
        inputs = config.get('input', {})
        crispr = config.get('crispr', {})
        pathway = config.get('pathway_analysis', {})

        def escaped(mapping, key, default='N/A'):
            value = mapping.get(key, default)
            if value in (None, ''):
                value = default
            return html.escape(str(value))

        return f'''                <!-- Configuration Section -->
                <section id="configuration" class="content-section">
                    <div class="section-header">
                        <h2>Analysis Configuration</h2>
                    </div>
                    <div class="section-body">
                        <div class="stats-grid">
                            <div class="stat-card">
                                <h3>Input Files</h3>
                                <ul class="configuration-list">
                                    <li><strong>Genome Assembly:</strong><br>{escaped(inputs, 'genome_assembly')}</li>
                                    <li><strong>Transcriptome:</strong><br>{escaped(inputs, 'transcriptome')}</li>
                                    <li><strong>Protein Sequences:</strong><br>{escaped(inputs, 'protein_sequences')}</li>
                                    <li><strong>Annotation Database:</strong><br>{escaped(inputs, 'sgd_db')}</li>
                                </ul>
                            </div>
                            <div class="stat-card">
                                <h3>CRISPR Settings</h3>
                                <ul class="configuration-list">
                                    <li><strong>PAM Pattern:</strong> {escaped(crispr, 'pam_pattern')}</li>
                                    <li><strong>Maximum Mismatch (fixed by tWhEEL):</strong> {CRISFLASH_MAX_MISMATCH}</li>
                                    <li><strong>Threads:</strong> {escaped(crispr, 'threads')}</li>
                                </ul>
                            </div>
                            <div class="stat-card">
                                <h3>Pathway Analysis Settings</h3>
                                <ul class="configuration-list">
                                    <li><strong>WikiPathways GPML Directory:</strong><br>{escaped(pathway, 'wikipathways_gpml_dir')}</li>
                                    <li><strong>Minimum Pathway Coverage:</strong> {escaped(pathway, 'minimum_pathway_coverage')}</li>
                                    <li><strong>Skip Low-Coverage Pathways:</strong> {escaped(pathway, 'skip_low_coverage_pathways')}</li>
                                </ul>
                            </div>
                        </div>
                    </div>
                </section>
'''

    @staticmethod
    def _pathway_report_css():
        """Return styles limited to PW-only tables, charts, and modals."""
        return '''
        .main-content { min-width: 0; }
        .content-section, .table-responsive { max-width: 100%; }
        #pw-kfi-assessment h3 { margin: 0 0 1rem; color: #333; }
        .pw-chart-block { margin-bottom: 2rem; }
        .pw-chart-block .chart-container { height: 350px; overflow-x: auto; }
        .pw-chart-block svg { max-width: 100%; height: auto; }
        .pw-results-block { border-top: 1px solid #dee2e6; padding-top: 1.5rem; }
        .pw-results-header { display: flex; justify-content: space-between; align-items: flex-start; gap: 1rem; flex-wrap: wrap; }
        .pw-results-header .filter-group { flex: 1; justify-content: flex-end; }
        #pathway-table th { cursor: pointer; white-space: nowrap; }
        #pathway-table th:last-child { cursor: default; }
        #pagination .page-btn { padding: 0.375rem 0.75rem; border: 1px solid #ced4da; border-radius: 0.375rem; background: #fff; color: #495057; cursor: pointer; }
        #pagination .page-btn:hover, #pagination .page-btn.active { background: #007bff; border-color: #007bff; color: #fff; }
        .modal { display: none; position: fixed; z-index: 2000; inset: 0; width: 100%; height: 100%; overflow: auto; background-color: rgba(0, 0, 0, 0.5); }
        .modal-content { position: relative; background-color: #fff; margin: 5vh auto; padding: 2rem; border: 1px solid #dee2e6; border-radius: 0.375rem; max-height: 90vh; overflow: auto; }
        .modal-content .close { position: absolute; top: 0.75rem; right: 1rem; color: #6c757d; font-size: 1.75rem; font-weight: 700; cursor: pointer; }
        .modal-content .close:hover { color: #212529; }
        .modal-content table { width: 100%; border-collapse: collapse; }
        .modal-content th, .modal-content td { padding: 0.75rem; border-bottom: 1px solid #dee2e6; text-align: left; }
        .modal-content th { background-color: #f8f9fa; position: sticky; top: 0; }
        .btn.secondary { background-color: #6c757d; border-color: #6c757d; color: #fff; }
        .configuration-list { text-align: left; font-size: 0.85rem; overflow-wrap: anywhere; }
        .configuration-list li { margin-bottom: 0.75rem; }
        @media (max-width: 768px) {
            .pw-results-header { display: block; }
            .pw-results-header .filter-group { justify-content: stretch; }
            .modal-content { width: calc(100% - 2rem) !important; margin: 1rem auto; padding: 1.25rem; }
        }
'''

    def create_pathway_html_template(
        self, pathways, summary_data, pathway_genes=None,
        common_results=None, common_report_generator=None,
        payload_script_src=None,
    ):
        """Compose PW-KFI content into a standard browser-rendering shell."""
        legacy_html = self._create_legacy_pathway_html_template(
            pathways, summary_data, pathway_genes
        )

        if common_results is None or common_report_generator is None:
            common_report_generator, common_results = (
                self._collect_shared_report_context()
            )

        base_results = self._pathway_shared_report_results(common_results)
        if payload_script_src is None:
            raise ValueError(
                'payload_script_src is required for the pathway HTML shell'
            )

        report_html = common_report_generator._generate_html_content(
            base_results,
            payload_script_src=payload_script_src,
        )

        modal_start = legacy_html.index('<div id="pathway-genes-modal"')
        script_start = legacy_html.index('<script>', modal_start)
        script_end = legacy_html.index('</script>', script_start) + len('</script>')
        pathway_modals = legacy_html[modal_start:script_start]
        pathway_script = legacy_html[script_start:script_end]

        report_html = self._replace_report_block(
            report_html,
            '                <!-- General Statistics Section -->',
            '                <!-- Gene Prediction Section -->',
            self._general_statistics_html(common_results, summary_data),
        )
        report_html = self._replace_report_block(
            report_html,
            '                <!-- GO-KFI Assessment Section -->',
            '                <!-- Configuration Section -->',
            self._pathway_assessment_html(summary_data),
        )
        configuration_end = '            </div>\n        </main>'
        report_html = self._replace_report_block(
            report_html,
            '                <!-- Configuration Section -->',
            configuration_end,
            self._analysis_configuration_html(common_results),
        )

        replacements = {
            '<title>tWhEEL - Transcriptome Functional Analysis Report</title>':
                '<title>tWhEEL - Pathway Analysis Report</title>',
            'gene ontology (GO)-based knock-out feasibility assessment':
                'pathway-based knock-out feasibility assessment',
            '<li><a href="#go-analysis" class="nav-link">GO-KFI Assessment</a></li>':
                '<li><a href="#pw-kfi-assessment" class="nav-link">PW-KFI Assessment</a></li>',
            '<li><a href="#configuration" class="nav-link">Configuration</a></li>':
                '<li><a href="#configuration" class="nav-link">Analysis Configuration</a></li>',
            'TransDecoder was used to predict protein-coding genes from the transcriptome assembly.':
                'Protein-coding sequences used in this analysis are summarized below.',
            'Protein sequences were searched against UniProt Knowledgebase (UniProtKB) database using FASTA for functional annotation.':
                'Functional annotations used for PW-KFI and guide RNA analysis are summarized below.',
            '<th>UniProt Name</th>': '<th>Annotation ID</th>',
        }
        for old, new in replacements.items():
            if report_html.count(old) != 1:
                raise ValueError(f'Unexpected shared report marker count: {old}')
            report_html = report_html.replace(old, new, 1)

        report_html = report_html.replace(
            '</head>',
            f'<style id="pw-kfi-styles">{self._pathway_report_css()}</style>\n</head>',
            1,
        )

        main_script = '    <script src="static/js/main.js"></script>'
        if report_html.count(main_script) != 1:
            raise ValueError('Unexpected shared main script marker count')
        report_html = report_html.replace(
            main_script,
            f'    {pathway_modals}\n\n{main_script}\n    {pathway_script}',
            1,
        )
        return report_html
    
    def generate_pathway_html_report(self, pathways, pathway_genes, summary_data):
        """Publish pathway browser data, then atomically publish its shell."""
        common_generator, common_results = self._collect_shared_report_context()

        base_results = self._pathway_shared_report_results(common_results)
        pathway_browser_data = self._build_pathway_browser_data(
            pathways, summary_data
        )
        report_name = self.pathway_report_config.get(
            'output_name', 'pathway_report.html'
        )
        report_path = os.path.join(self.output_dir, report_name)
        with common_generator.report_payload_store.publication_lock():
            common_generator._create_css_files()
            common_generator._create_js_files()
            common_generator._copy_static_assets(common_results)
            payload = common_generator.publish_browser_payload(
                base_results,
                namespace='pathway',
                extra_globals={
                    '__TWHEEL_PATHWAY_PAYLOAD__': pathway_browser_data,
                },
            )
            html_content = self.create_pathway_html_template(
                pathways,
                summary_data,
                pathway_genes,
                common_results=common_results,
                common_report_generator=common_generator,
                payload_script_src=payload['core_relative_path'],
            )

            # Publish the shell only after its immutable payload is complete.
            common_generator._atomic_write_text(report_path, html_content)

            summary_path = os.path.join(
                self.web_data_dir, 'pathway_summary.json'
            )
            if os.path.isfile(summary_path):
                try:
                    with open(summary_path, 'r', encoding='utf-8') as handle:
                        summary = json.load(handle)
                    summary['browser_payload'] = {
                        'schema_version': (
                            common_generator.report_payload_store.SCHEMA_VERSION
                        ),
                        'namespace': 'pathway',
                        'generation_id': payload['generation_id'],
                        'manifest': os.path.relpath(
                            payload['manifest_path'], self.output_dir
                        ).replace(os.sep, '/'),
                    }
                    common_generator._atomic_write_text(
                        summary_path,
                        json.dumps(
                            summary, ensure_ascii=False, indent=2
                        ) + '\n',
                    )
                except (OSError, TypeError, ValueError) as error:
                    self.logger.warning(
                        "Pathway HTML is complete, but its summary payload "
                        f"pointer could not be updated: {error}"
                    )
            common_generator.report_payload_store.mark_committed(
                'pathway', payload['generation_id']
            )
            common_generator.report_payload_store.cleanup_uncommitted_locked(
                'pathway', payload['generation_id']
            )

        self.logger.info(f"Pathway HTML Report Generation: {report_path}")
        return report_path
    
    def run(self):
        """Full execution of pathway report generation"""
        try:
            self.logger.info("Started generating pathway report....")
            
            # データロード
            pathways, pathway_genes = self.load_pathway_analysis_results()
            manifest_path = str(self.report_source_store.manifest_path)
            
            # パスウェイ別ガイドCSVデータ生成
            pathway_guide_data = self.generate_pathway_guide_csv_data(pathways)
            
            # JSONサマリー生成（ガイドデータとSVGデータを含む）
            summary_data = self.generate_pathway_json_summary(
                pathways, pathway_guide_data
            )
            
            # HTMLレポート生成
            report_path = self.generate_pathway_html_report(pathways, pathway_genes, summary_data)
            
            # 結果まとめ
            results = {
                'summary_json_path': os.path.join(self.web_data_dir, 'pathway_summary.json'),
                'html_report_path': report_path,
                'report_sources_manifest': manifest_path,
                'statistics': summary_data['statistics']['pathway_analysis'],
                'pathway_guide_data': pathway_guide_data
            }
            
            self.logger.info(f"Pathway report generation complete.: {report_path}")
            return results
            
        except Exception as e:
            self.logger.error(f"An error occurred during pathway report generation.: {e}")
            raise
