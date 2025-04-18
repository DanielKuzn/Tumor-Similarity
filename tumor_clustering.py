import argparse
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from scipy.spatial.distance import squareform

# Function to calculate similarity between two samples based on proximity
def calculate_similarity(sample1_variants, sample2_variants, max_distance):
    shared_variants = 0
    proximity_score = 0
    for v1 in sample1_variants:
        for v2 in sample2_variants:
            if v1[0] == v2[0]:
                distance = abs(v1[1] - v2[1])
                if distance <= max_distance:
                    proximity_score += np.exp(-distance / max_distance)
                if v1[1] == v2[1]:
                    shared_variants += 1
    return shared_variants + proximity_score

# Main function to process the data and create plots
def main(input_file, max_distance):
    # Step 1: Load the input TSV data into a DataFrame
    df = pd.read_csv(input_file, delimiter="\t", header=0, names=['chromosome', 'start', 'end', 'tumor_type', 'sample'])

    # Step 2: Create a dictionary of variant positions per sample
    samples = df['sample'].unique()
    sample_variants = {sample: [] for sample in samples}
    
    for _, row in df.iterrows():
        sample_variants[row['sample']].append((row['chromosome'], row['start'], row['end']))

    # Step 3: Create a distance matrix for the samples
    distance_matrix = np.zeros((len(samples), len(samples)))
    for i, sample1 in enumerate(samples):
        for j, sample2 in enumerate(samples):
            if i < j:
                print(f"Comparison of sample {i + 1} with {j + 1}. Total is {len(samples)} samples.")
                sim = calculate_similarity(sample_variants[sample1], sample_variants[sample2], max_distance)
                distance_matrix[i, j] = 1 / (1 + sim)
                print(f"Distance is: {distance_matrix[i, j]}.")
            elif i == j:
                distance_matrix[i, j] = 0
            else:
                distance_matrix[i, j] = distance_matrix[j, i]

    # Save the distance matrix to a file
    np.save('distance_matrix.npy', distance_matrix)

    # Step 4: Perform hierarchical clustering (dendrogram)
    condensed_distance_matrix = squareform(distance_matrix)
    linked = linkage(condensed_distance_matrix, 'ward')

    # Get tumor types for each sample
    tumor_types = df[['sample', 'tumor_type']].drop_duplicates().set_index('sample')['tumor_type']
    colors = tumor_types[samples].values  # Use the tumor types to color the samples

    # Color map for tumor types
    unique_tumor_types = tumor_types.unique()
    color_map = {tumor: idx for idx, tumor in enumerate(unique_tumor_types)}

    # Create color labels for the dendrogram
    dendro_labels = [color_map[tumor] for tumor in colors]

    plt.figure(figsize=(10, 7))
    dendrogram(
        linked,
        labels=samples,
        leaf_rotation=0, orientation="left",
        #color_threshold=0,  # Disable default color threshold
        #above_threshold_color='grey',  # Default color for all labels
        #below_threshold_color='grey',  # Default color for all labels
    )

    # Apply colors for tumor types
    ax = plt.gca()
    for i, label in enumerate(ax.get_xticklabels()):
        label.set_color(plt.cm.viridis(dendro_labels[i] / len(unique_tumor_types)))  # Color labels based on tumor type
    plt.title('Dendrogram of Tumor Sample Similarity')
    plt.xlabel('Distance')
    #plt.xscale('log')
    plt.ylabel('Sample')
    plt.xticks(rotation=90)
    plt.tight_layout()

    # Save the dendrogram plot
    plt.savefig('dendrogram.png')
    plt.show()

    # Step 5: Perform t-SNE for dimensionality reduction
    tsne = TSNE(n_components=2, metric='precomputed', perplexity=np.min([len(samples) - 1, 30]), init='random')
    tsne_results = tsne.fit_transform(distance_matrix)

    # Plot the t-SNE result with tumor type colors
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(tsne_results[:, 0], tsne_results[:, 1], c=[color_map[tumor] for tumor in colors], cmap='viridis')
    plt.title('t-SNE of Tumor Samples')
    plt.xlabel('t-SNE Component 1')
    plt.ylabel('t-SNE Component 2')
    plt.colorbar(scatter, label='Tumor Type')

    # Save the t-SNE plot
    plt.savefig('tsne_plot.png')
    plt.show()

# Command-line interface
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster tumor samples based on proximity of genomic variants")
    parser.add_argument("input_file", help="Path to the input TSV file with sample data", type=str)
    parser.add_argument("max_distance", help="Max distance between variants to assume as similar. Default is 30 (average size of TF binding size)", type=int, default=30, nargs="?")
    args = parser.parse_args()

    # Run the main function with the input file
    main(args.input_file, args.max_distance)
