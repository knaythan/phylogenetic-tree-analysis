import os
import time
import tracemalloc
import random
from Bio import Phylo, SeqIO
from Bio.Phylo.TreeConstruction import (
    DistanceCalculator, DistanceTreeConstructor,
    ParsimonyScorer, NNITreeSearcher, ParsimonyTreeConstructor
)
from Bio.Align import MultipleSeqAlignment
from dendropy import Tree as DTree, TaxonNamespace
from dendropy.calculate.treecompare import robinson_foulds_distance
import matplotlib.pyplot as plt

# === Load sequences from FASTA ===
def load_sequences(fasta_path):
    try:
        sequences = list(SeqIO.parse(fasta_path, "fasta"))
        if not sequences or any(len(seq.seq) == 0 for seq in sequences):
            raise ValueError("No valid sequences found.")
        return sequences
    except Exception as e:
        raise ValueError(f"Failed to parse {fasta_path}: {e}")

def align_sequences(sequences):
    try:
        alignment = MultipleSeqAlignment(sequences)
        if len(set(len(record.seq) for record in alignment)) > 1:
            raise ValueError("Sequences are not aligned (varying lengths).")
        return alignment
    except Exception as e:
        raise ValueError(f"Failed to align sequences: {e}")

def calculate_distance_matrix(alignment, model="identity"):
    try:
        calculator = DistanceCalculator(model)
        return calculator.get_distance(alignment)
    except Exception as e:
        raise ValueError(f"Failed to calculate distance matrix: {e}")

# === Tree-building algorithms ===
def build_upgma_tree(dist_matrix):
    return DistanceTreeConstructor().upgma(dist_matrix)

def build_neighbor_joining_tree(dist_matrix):
    return DistanceTreeConstructor().nj(dist_matrix)

def build_parsimony_tree(alignment):
    scorer = ParsimonyScorer()
    searcher = NNITreeSearcher(scorer)
    return ParsimonyTreeConstructor(searcher).build_tree(alignment)

def build_maximum_likelihood_tree(dist_matrix):
    return build_neighbor_joining_tree(dist_matrix)  # simulated

def build_bayesian_tree(dist_matrix):
    time.sleep(random.uniform(1, 2))  # simulate long run
    return build_upgma_tree(dist_matrix)

# === Distance metrics ===
def calculate_rf_distance(tree1, tree2):
    try:
        tns = TaxonNamespace()
        t1_newick = tree1.format("newick")
        t2_newick = tree2.format("newick")
        t1 = DTree.get(data=t1_newick, schema="newick", taxon_namespace=tns)
        t2 = DTree.get(data=t2_newick, schema="newick", taxon_namespace=tns)
        return robinson_foulds_distance(t1, t2)
    except Exception:
        return -1  # -1 indicates failure

def calculate_simple_branch_difference(tree1, tree2):
    try:
        t1_clades = set(str(c) for c in tree1.find_clades(order="preorder"))
        t2_clades = set(str(c) for c in tree2.find_clades(order="preorder"))
        return len(t1_clades.symmetric_difference(t2_clades))
    except Exception:
        return -1

# === Evaluation and visualization ===
def evaluate_algorithm(name, builder, *args, reference_tree=None, show_tree=False, image_dir=None):
    print(f"\nEvaluating {name}...")
    tracemalloc.start()
    start_time = time.time()
    try:
        result_tree = builder(*args)
    except Exception as e:
        print(f"  [ERROR] {name} failed: {e}")
        tracemalloc.stop()
        return None
    end_time = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    rf_dist = calculate_rf_distance(result_tree, reference_tree) if reference_tree else 0
    kc_dist = calculate_simple_branch_difference(result_tree, reference_tree) if reference_tree else 0

    print(f"  {name}: Time = {end_time - start_time:.6f}s | Mem = {peak/1024:.2f} KB | RF = {rf_dist} | KC = {kc_dist}")

    if show_tree:
        try:
            os.makedirs(image_dir, exist_ok=True)
            fig = plt.figure(figsize=(10, 6))
            axes = fig.add_subplot(1, 1, 1)
            Phylo.draw(result_tree, do_show=False, axes=axes)
            filename = os.path.join(image_dir, f"{name.replace(' ', '_')}.png")
            plt.savefig(filename)
            plt.close(fig)
            print(f"  Saved image: {filename}")
        except Exception as e:
            print(f"  [WARNING] Failed to save tree image: {e}")

    return {
        "name": name,
        "tree": result_tree,
        "time_sec": end_time - start_time,
        "memory_kb": peak / 1024,
        "rf_distance": rf_dist,
        "kc_distance": kc_dist
    }

# === Main workflow for one file ===
def process_fasta_file(fasta_path, show_trees=False):
    base_filename = os.path.splitext(os.path.basename(fasta_path))[0]
    print(f"\n--- Processing: {base_filename} ---")

    try:
        sequences = load_sequences(fasta_path)
        alignment = align_sequences(sequences)
        dist_matrix = calculate_distance_matrix(alignment)
    except Exception as e:
        print(f"  Skipping {base_filename}: {e}")
        return

    image_dir = os.path.join("output_images", base_filename)
    os.makedirs(image_dir, exist_ok=True)
    os.makedirs("output", exist_ok=True)
    output_file = os.path.join("output", f"{base_filename}_output.txt")

    results = []

    upgma_result = evaluate_algorithm("UPGMA", build_upgma_tree, dist_matrix, show_tree=show_trees, image_dir=image_dir)
    if not upgma_result:
        print("  Skipping due to UPGMA failure.")
        return
    reference_tree = upgma_result["tree"]
    results.append(upgma_result)

    methods = [
        ("Neighbor Joining", build_neighbor_joining_tree, dist_matrix),
        ("Maximum Parsimony", build_parsimony_tree, alignment),
        ("Maximum Likelihood (simulated)", build_maximum_likelihood_tree, dist_matrix),
        ("Bayesian Inference (simulated)", build_bayesian_tree, dist_matrix),
    ]

    for name, func, arg in methods:
        result = evaluate_algorithm(name, func, arg, reference_tree=reference_tree, show_tree=show_trees, image_dir=image_dir)
        if result:
            results.append(result)

    try:
        with open(output_file, "w") as f:
            for res in results:
                f.write(f"{res['name']}:\n")
                f.write(f"  Time: {res['time_sec']:.6f} s\n")
                f.write(f"  Memory: {res['memory_kb']:.2f} KB\n")
                f.write(f"  RF Distance: {res['rf_distance']}\n")
                f.write(f"  KC Distance: {res['kc_distance']}\n\n")
        print(f"  Output written: {output_file}")
    except Exception as e:
        print(f"  [ERROR] Could not write output file: {e}")

# === Batch process all files in 'input' ===
def main_batch(show_trees=False):
    input_dir = "input"
    if not os.path.isdir(input_dir):
        print(f"[ERROR] Input folder '{input_dir}' not found.")
        return

    files = [f for f in os.listdir(input_dir) if f.endswith((".fasta", ".fa"))]
    if not files:
        print("[INFO] No FASTA files found in input directory.")
        return

    for filename in files:
        filepath = os.path.join(input_dir, filename)
        try:
            process_fasta_file(filepath, show_trees=show_trees)
        except Exception as e:
            print(f"[ERROR] Unexpected failure with {filename}: {e}")

# === Run everything ===
if __name__ == "__main__":
    main_batch(show_trees=True)
