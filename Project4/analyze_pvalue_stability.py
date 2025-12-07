"""
P-value Stability Analysis

Analyzes the stability of bootstrap p-values across different random seeds
and iteration counts.

Author: Nicholas Allen
Course: CS274 - Algorithms in Molecular Biology
"""

import subprocess
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def generate_random_seeds(n_seeds=100, seed=42):
    """
    Generate random seeds for reproducibility testing.
    
    Inputs:
        n_seeds: int, number of random seeds to generate (default 100)
        seed: int, seed for generating seeds (default 42)
    
    Returns:
        list of random seed integers
    """
    random.seed(seed)
    return [random.randint(1, 10000) for _ in range(n_seeds)]


def run_pvalue_analysis(protein_a, protein_b, iteration_counts, random_seeds, 
                        drugs_csv='data/drugs.csv', targets_csv='data/targets.csv'):
    """
    Run p-value calculations across multiple seeds and iteration counts.
    
    Inputs:
        protein_a: string, first protein uniprot_accession ID
        protein_b: string, second protein uniprot_accession ID
        iteration_counts: list of ints, bootstrap iteration counts to test
        random_seeds: list of ints, random seeds to use
        drugs_csv: string, path to drugs CSV file
        targets_csv: string, path to targets CSV file
    
    Returns:
        dict mapping iteration_count -> list of p-values
    """
    results = {n: [] for n in iteration_counts}
    total_runs = len(iteration_counts) * len(random_seeds)
    current_run = 0
    
    print(f"Running p-value analysis...")
    print(f"Proteins: {protein_a} vs {protein_b}")
    print(f"Iteration counts: {iteration_counts}")
    print(f"Number of seeds: {len(random_seeds)}")
    print(f"Total runs: {total_runs}\n")
    
    for iteration_count in iteration_counts:
        print(f"Testing with n={iteration_count} iterations:")
        
        for seed in random_seeds:
            # Run pvalue.py with current parameters
            cmd = [
                'python3', 'pvalue.py',
                '-n', str(iteration_count),
                '-r', str(seed),
                drugs_csv, targets_csv,
                protein_a, protein_b
            ]
            
            try:
                # Run command and capture output
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                p_value = float(result.stdout.strip())
                results[iteration_count].append(p_value)
                
                current_run += 1
                if current_run % 10 == 0:
                    print(f"  Progress: {current_run}/{total_runs} runs completed")
                    
            except subprocess.CalledProcessError as e:
                print(f"  Error running seed {seed}: {e}")
            except ValueError as e:
                print(f"  Error parsing p-value for seed {seed}: {e}")
        
        print(f"  Completed {len(results[iteration_count])} runs for n={iteration_count}\n")
    
    return results


def calculate_statistics(results):
    """
    Calculate mean and standard deviation for each iteration count.
    
    Inputs:
        results: dict mapping iteration_count -> list of p-values
    
    Returns:
        pandas DataFrame with statistics
    """
    stats = []
    
    for iteration_count, p_values in sorted(results.items()):
        stats.append({
            'iteration_count': iteration_count,
            'n_samples': len(p_values),
            'mean': np.mean(p_values),
            'std': np.std(p_values, ddof=1),  # Sample standard deviation
            'min': np.min(p_values),
            'max': np.max(p_values),
            'median': np.median(p_values)
        })
    
    return pd.DataFrame(stats)


def save_results(results, stats_df, output_dir='results'):
    """
    Save results and statistics to files.
    
    Inputs:
        results: dict mapping iteration_count -> list of p-values
        stats_df: pandas DataFrame with statistics
        output_dir: string, directory to save results
    
    Returns:
        None (saves files to disk)
    """
    import os
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Save raw p-values for each iteration count
    for iteration_count, p_values in results.items():
        filename = f"{output_dir}/pvalues_n{iteration_count}.txt"
        with open(filename, 'w') as f:
            for p_val in p_values:
                f.write(f"{p_val}\n")
        print(f"Saved p-values to: {filename}")
    
    # Save statistics
    stats_filename = f"{output_dir}/pvalue_statistics.csv"
    stats_df.to_csv(stats_filename, index=False)
    print(f"Saved statistics to: {stats_filename}")


def plot_histograms(results, output_dir='results'):
    """
    Generate histograms for p-value distributions.
    
    Inputs:
        results: dict mapping iteration_count -> list of p-values
        output_dir: string, directory to save plots
    
    Returns:
        None (saves plots to disk)
    """
    import os
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Create subplot for each iteration count
    fig, axes = plt.subplots(1, len(results), figsize=(15, 4))
    
    if len(results) == 1:
        axes = [axes]
    
    for idx, (iteration_count, p_values) in enumerate(sorted(results.items())):
        ax = axes[idx]
        
        # Plot histogram
        ax.hist(p_values, bins=20, edgecolor='black', alpha=0.7)
        ax.axvline(np.mean(p_values), color='red', linestyle='--', 
                   label=f'Mean: {np.mean(p_values):.4f}')
        ax.set_xlabel('P-value')
        ax.set_ylabel('Frequency')
        ax.set_title(f'n={iteration_count} iterations\n(μ={np.mean(p_values):.4f}, σ={np.std(p_values, ddof=1):.4f})')
        ax.legend()
        ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    hist_filename = f"{output_dir}/pvalue_distributions.png"
    plt.savefig(hist_filename, dpi=150)
    plt.close()
    
    print(f"Saved histogram to: {hist_filename}")


if __name__ == "__main__":
    # Configuration
    PROTEIN_A = "P54577"
    PROTEIN_B = "Q7RTX0"
    ITERATION_COUNTS = [100, 500, 1000]
    N_SEEDS = 100
    
    # Generate random seeds
    print("Generating random seeds...")
    random_seeds = generate_random_seeds(N_SEEDS)
    print(f"Generated {len(random_seeds)} random seeds\n")
    
    # Run p-value analysis
    results = run_pvalue_analysis(
        protein_a=PROTEIN_A,
        protein_b=PROTEIN_B,
        iteration_counts=ITERATION_COUNTS,
        random_seeds=random_seeds
    )
    
    # Calculate statistics
    print("\nCalculating statistics...")
    stats_df = calculate_statistics(results)
    
    # Display statistics
    print("\n" + "="*70)
    print("P-VALUE STABILITY ANALYSIS RESULTS")
    print("="*70)
    print(stats_df.to_string(index=False))
    print("="*70)
    
    # Save results
    print("\nSaving results...")
    save_results(results, stats_df)
    
    # Generate histograms
    print("\nGenerating histograms...")
    plot_histograms(results)
    
    print("\n✓ Analysis complete!")
