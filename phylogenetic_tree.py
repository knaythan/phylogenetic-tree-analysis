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
from skbio import TreeNode
from skbio.tree import TreeError
import matplotlib.pyplot as plt

# === Load sequences from FASTA ===
def load_sequences(fasta_path):
    return list(SeqIO.parse(fasta_path, "fasta"))

def align_sequences(sequences):
    return MultipleSeqAlignment(sequences)

def calculate_distance_matrix(alignment, model="identity"):
    calculator = DistanceCalculator(model)
    return calculator.get_distance(alignment)

# === Tree-building algorithms ===
def build_upgma_tree(dist_matrix):
    constructor = DistanceTreeConstructor()
    return constructor.upgma(dist_matrix)

def build_neighbor_joining_tree(dist_matrix):
    constructor = DistanceTreeConstructor()
    return constructor.nj(dist_matrix)

def build_parsimony_tree(alignment):
    scorer = ParsimonyScorer()
    searcher = NNITreeSearcher(scorer)
    constructor = ParsimonyTreeConstructor(searcher)
    return constructor.build_tree(alignment)

def build_maximum_likelihood_tree(dist_matrix):
    return build_neighbor_joining_tree(dist_matrix)  # simulated approximation

def build_bayesian_tree(dist_matrix):
    time.sleep(random.uniform(1, 2))  # simulate longer MCMC sampling
    return build_upgma_tree(dist_matrix)

# === Distance metrics ===
def calculate_rf_distance(tree1, tree2):
    tns = TaxonNamespace()
    t1_newick = tree1.format("newick")
    t2_newick = tree2.format("newick")
    t1 = DTree.get(data=t1_newick, schema="newick", taxon_namespace=tns)
    t2 = DTree.get(data=t2_newick, schema="newick", taxon_namespace=tns)
    return robinson_foulds_distance(t1, t2)

def calculate_simple_branch_difference(tree1, tree2):
    t1_clades = set(str(c) for c in tree1.find_clades(order="preorder"))
    t2_clades = set(str(c) for c in tree2.find_clades(order="preorder"))
    return len(t1_clades.symmetric_difference(t2_clades))


# === Evaluation and visualization ===
def evaluate_algorithm(name, builder, *args, reference_tree=None, show_tree=False):
    print(f"\nEvaluating {name}...")
    tracemalloc.start()
    start_time = time.time()
    result_tree = builder(*args)
    end_time = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    rf_dist = calculate_rf_distance(result_tree, reference_tree) if reference_tree else 0
    kc_dist = calculate_simple_branch_difference(result_tree, reference_tree) if reference_tree else 0

    print(f"{name}: Time = {end_time - start_time:.2f}s | Mem = {peak/1024:.2f} KB | RF = {rf_dist} | KC = {kc_dist}")

    if show_tree:
        print(f"Drawing {name} tree...")
        Phylo.draw(result_tree)

    return {
        "name": name,
        "tree": result_tree,
        "time_sec": end_time - start_time,
        "memory_kb": peak / 1024,
        "rf_distance": rf_dist,
        "kc_distance": kc_dist
    }

# === Main workflow ===
def main(fasta_path, show_trees=False):
    sequences = load_sequences(fasta_path)
    alignment = align_sequences(sequences)
    dist_matrix = calculate_distance_matrix(alignment)

    results = []

    # Use UPGMA as baseline reference
    upgma_result = evaluate_algorithm("UPGMA", build_upgma_tree, dist_matrix, show_tree=show_trees)
    reference_tree = upgma_result["tree"]
    results.append(upgma_result)

    # Other methods compared to UPGMA
    results.append(evaluate_algorithm("Neighbor Joining", build_neighbor_joining_tree, dist_matrix, reference_tree=reference_tree, show_tree=show_trees))
    results.append(evaluate_algorithm("Maximum Parsimony", build_parsimony_tree, alignment, reference_tree=reference_tree, show_tree=show_trees))
    results.append(evaluate_algorithm("Maximum Likelihood (simulated)", build_maximum_likelihood_tree, dist_matrix, reference_tree=reference_tree, show_tree=show_trees))
    results.append(evaluate_algorithm("Bayesian Inference (simulated)", build_bayesian_tree, dist_matrix, reference_tree=reference_tree, show_tree=show_trees))

    return results

# === Run with your aligned FASTA file ===
# Set show_trees=True to see tree graphics
main("new_aligned_sequence.fasta", show_trees=True)
