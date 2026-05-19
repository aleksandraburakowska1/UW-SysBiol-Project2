#!/usr/bin/env python3
"""Batch analyze spatiotemporal signal propagation near field-of-view boundaries.

Algorithm Description:
1. Load single-cell trajectory data pass-by-pass for specified Experiment IDs.
2. For each site within an experiment, reconstruct the spatiotemporal graph 
   (spatial exposure and near-future target cell jumps).
3. Classify cells dynamically at each frame into 'edge' or 'central' categories
   based on coordinate percentiles to account for boundary constraints.
4. Separate the graph nodes by category and calculate isolated propagation 
   metrics: Relative Risk (RR) and Risk Difference (RD) for both zones.
5. Apply a Mann-Whitney U test to evaluate if edge boundary effects lead to a 
   statistically significant drop in local neighbor connectivity counts.
6. Aggregate results for all processed blocks and export a structured CSV summary.
"""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

#Import core processing functions with try/except fallback structure
try:
    from spatiotemporal_signal_propagation import (
        add_track_deltas,
        annotate_spatial_exposure,
        assign_jump_events,
        build_spatial_edges,
        build_temporal_edges,
        compute_future_jump_flags,
        load_metadata,
        resolve_path,
    )
except ModuleNotFoundError:
    from scripts.spatiotemporal_signal_propagation import (
        add_track_deltas,
        annotate_spatial_exposure,
        assign_jump_events,
        build_spatial_edges,
        build_temporal_edges,
        compute_future_jump_flags,
        load_metadata,
        resolve_path,
    )

# Strict column requirements matching the primary schema
REQUIRED_COLUMNS = [
    'Exp_ID',
    'Image_Metadata_Site',
    'track_id',
    'Image_Metadata_T',
    'objNuclei_Location_Center_X',
    'objNuclei_Location_Center_Y',
]


def parse_args() -> argparse.Namespace:
    """Read command-line options in the official project architecture style."""
    parser = argparse.ArgumentParser(
        description='Quantify edge vs central propagation behaviors across multiple experiments.'
    )
    parser.add_argument(
        '--data-path',
        type=str,
        default='single-cell-tracks_exp1-6_noErbB2.csv.gz',
        help='Path to the primary compressed single-cell tracking CSV table.',
    )
    parser.add_argument(
        '--metadata-path',
        type=str,
        default='01-readme-experiment-description_2022-04-05.csv',
        help='Path to the experimental metadata reference document.',
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='analysis_outputs/spatial_heterogeneity',
        help='Directory path where final CSV results table will be stored.',
    )
    parser.add_argument(
        '--signal-col',
        type=str,
        default='ERKKTR_ratio',
        help='Target biosensor column name to compute propagation jumps for.',
    )
    parser.add_argument(
        '--chosen-exp-ids',
        type=int,
        nargs='+',
        default=[1],
        help='Space-separated list of Experiment IDs to run batch processing for.',
    )
    parser.add_argument(
        '--edge-margin-percentile',
        type=float,
        default=10.0,
        help='Percentile cut-off from each field-of-view border to define edge cells.',
    )
    parser.add_argument(
        '--spatial-radius',
        type=float,
        default=60.0,
        help='Spatial network reach radius in microns to map local neighbors.',
    )
    parser.add_argument(
        '--future-window-frames',
        type=int,
        default=3,
        help='Number of forward time-frames monitored to detect target self-jumps.',
    )
    parser.add_argument(
        '--jump-quantile',
        type=float,
        default=0.90,
        help='Quantile threshold on delta values to classify a dynamic signal jump.',
    )
    parser.add_argument(
        '--chunksize',
        type=int,
        default=1000000,
        help='Row size threshold for chunk-by-chunk stream loading optimizations.',
    )
    return parser.parse_args()


def load_single_experiment_data(
    data_path: Path, exp_id: int, target_columns: list[str], chunksize: int
) -> pd.DataFrame:
    """Efficiently streams data chunks from disk to load a single experiment block."""
    parts = []
    for chunk in pd.read_csv(data_path, usecols=target_columns, chunksize=chunksize):
        mask = chunk['Exp_ID'] == exp_id
        if mask.any():
            parts.append(chunk.loc[mask].copy())
    
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def classify_frame_edges(group: pd.DataFrame, percentile: float) -> np.ndarray:
    """Classifies cells into edge or central categories dynamically per frame."""
    x_low, x_high = np.percentile(group['objNuclei_Location_Center_X'], [percentile, 100 - percentile])
    y_low, y_high = np.percentile(group['objNuclei_Location_Center_Y'], [percentile, 100 - percentile])
    
    is_edge = (
        (group['objNuclei_Location_Center_X'] <= x_low) |
        (group['objNuclei_Location_Center_X'] >= x_high) |
        (group['objNuclei_Location_Center_Y'] <= y_low) |
        (group['objNuclei_Location_Center_Y'] >= y_high)
    )
    return np.where(is_edge, 'edge', 'central')


def calculate_propagation_metrics(df: pd.DataFrame) -> dict[str, float]:
    """Computes basic exposed/unexposed probabilities, Relative Risk and Risk Difference."""
    exposed = df[df['neighbor_jump_now'] == True]
    unexposed = df[df['neighbor_jump_now'] == False]
    
    p_exposed = exposed['future_self_jump'].mean() if len(exposed) > 0 else np.nan
    p_unexposed = unexposed['future_self_jump'].mean() if len(unexposed) > 0 else np.nan
    
    if pd.isna(p_exposed) or pd.isna(p_unexposed) or p_unexposed == 0:
        rr = np.nan
    else:
        rr = p_exposed / p_unexposed
        
    rd = p_exposed - p_unexposed if not (pd.isna(p_exposed) or pd.isna(p_unexposed)) else np.nan
    
    return {'RR': rr, 'RD': rd}


def run_neighbor_count_test(edge_df: pd.DataFrame, central_df: pd.DataFrame) -> float:
    """Applies a Mann-Whitney U test to compare physical neighbor count distributions."""
    if edge_df.empty or central_df.empty:
        return np.nan
    stat, p_value = mannwhitneyu(
        edge_df['neighbor_count'], 
        central_df['neighbor_count'], 
        alternative='two-sided'
    )
    return p_value


def main() -> None:
    """Executes the complete spatial stratification batch pipeline orchestration loop."""
    args = parse_args()

    data_path = resolve_path(Path(args.data_path))
    metadata_path = resolve_path(Path(args.metadata_path))
    output_dir = resolve_path(Path(args.output_dir))
    output_dir.mkdir(parents=True, exist_ok=True)
    
    meta = load_metadata(metadata_path)
    frame_to_minutes = float(meta['Acquisition_frequency_min'].iloc[0])
    
    # Track complete feature columns subset for safe stream loading
    load_cols = list(set(REQUIRED_COLUMNS + [args.signal_col]))
    
    all_compiled_site_rows = []
    
    # Iterate sequentially over specified experiments
    for exp_id in args.chosen_exp_ids:
        print(f"LOADING EXPERIMENT {exp_id} INTO MEMORY")
        
        exp_data = load_single_experiment_data(data_path, exp_id, load_cols, args.chunksize)
        if exp_data.empty:
            print(f"Warning: No tracking data entries found for Exp_ID {exp_id}. Skipping.")
            continue
            
        unique_sites = sorted(exp_data['Image_Metadata_Site'].unique())
        
        # Internal iteration block across individual fields of view (sites)
        for site_id in unique_sites:
            # Output minimal clear runtime logging status telemetry
            print(f" -> [Exp {exp_id}, Site {site_id}] Building graph & boundaries...")
            
            # Isolate the current target site block from cached experiment matrices
            block = exp_data[exp_data['Image_Metadata_Site'] == site_id].copy()
            block = block.sort_values(['track_id', 'Image_Metadata_T']).reset_index(drop=True)
            
            if block.empty:
                continue
                
            # Reconstruct spatial and temporal graph linkages natively
            block, threshold = add_track_deltas(block, args.signal_col, frame_to_minutes, args.jump_quantile)
            block = assign_jump_events(block, threshold)
            spatial_edges = build_spatial_edges(block, args.spatial_radius)
            temporal_edges = build_temporal_edges(block, frame_to_minutes)
            block = compute_future_jump_flags(block, args.future_window_frames)
            block = annotate_spatial_exposure(block, spatial_edges)
            
            # Apply dynamic frame-by-frame edge boundary calculations
            block['spatial_category'] = block.groupby('Image_Metadata_T').apply(
                lambda g: pd.Series(classify_frame_edges(g, args.edge_margin_percentile), index=g.index)
            ).reset_index(level=0, drop=True)
            
            # Filter rows into strict spatial category subsets
            edge_df = block[block['spatial_category'] == 'edge']
            central_df = block[block['spatial_category'] == 'central']
            
            # Extract independent propagation summaries
            edge_metrics = calculate_propagation_metrics(edge_df)
            central_metrics = calculate_propagation_metrics(central_df)
            
            # Perform Mann-Whitney U test on network architecture node indices
            p_mw = run_neighbor_count_test(edge_df, central_df)
            
            # Record the compiled metrics dictionary row elements
            all_compiled_site_rows.append({
                'Exp_ID': exp_id,
                'Site_ID': site_id,
                'rr_edge': edge_metrics['RR'],
                'rr_central': central_metrics['RR'],
                'rd_edge': edge_metrics['RD'],
                'rd_central': central_metrics['RD'],
                'mean_neighbors_edge': edge_df['neighbor_count'].mean() if len(edge_df) else np.nan,
                'mean_neighbors_central': central_df['neighbor_count'].mean() if len(central_df) else np.nan,
                'mw_p_value': f"{p_mw:.4e}" if pd.notna(p_mw) else "NaN"
            })
            
        # Aggressively purge completed experiment variables from host RAM
        del exp_data
        gc.collect()
        print(f"Memory cleared for Experiment {exp_id}.")
        
    # Build complete execution dataframe matrices and save outputs to disk
    if all_compiled_site_rows:
        summary_df = pd.DataFrame(all_compiled_site_rows)
        output_csv_file = output_dir / 'spatial_edge_comparison_summary.csv'
        summary_df.to_csv(output_csv_file, index=False)
        
        print(f"\nBatch analysis complete for all selected experiments.")
        print(f"Results table successfully written to: {output_csv_file}")
        print(summary_df.head(10).to_string(index=False))
    else:
        print("\nPipeline finished: No valid processing blocks generated records.")


if __name__ == '__main__':
    main()