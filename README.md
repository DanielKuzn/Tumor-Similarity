# Tumor Sample Similarity Analysis

This Python script performs an analysis of tumor sample similarities based on variant positions across chromosomes. It calculates pairwise sample similarity based on both the presence of shared variants and their proximity on the same chromosome. The output includes a hierarchical clustering dendrogram and a t-SNE visualization, with samples colored by tumor type.

## Features

- **Variant Proximity Consideration**  
  Similarity is calculated based on the proximity of variants within the same chromosome, allowing samples to cluster even when variants are not exactly shared but are nearby.

- **Hierarchical Clustering**  
  Generates a dendrogram to visualize relationships between tumor samples.

- **t-SNE Visualization**  
  Reduces dimensionality and plots samples in 2D, colored by tumor type for visual clarity.

- **Customizable Analysis**  
  Parameters like the proximity decay function and clustering method can be tuned within the script.

---

## Prerequisites

- Python 3.x  
- Required Python libraries:
  - `pandas`
  - `numpy`
  - `scipy`
  - `matplotlib`
  - `scikit-learn`

**Install dependencies:**

```bash
pip install pandas numpy scipy matplotlib scikit-learn
```

---

## Input Format

The input file should be in `.csv` or `.tsv` format with the following columns:

| Column Name | Description                                  |
|-------------|----------------------------------------------|
| Chromosome  | Chromosome number                            |
| Start       | Variant start position                       |
| End         | Variant end position (can be same as start)  |
| Tumor Type  | Type of tumor                                |
| Sample      | Tumor sample identifier                      |

**Example input:**

```csv
Chromosome,Start,End,Tumor Type,Sample
1,1230448,1230448,Lymph-BNHL,DO46416
1,1609723,1609723,Lymph-BNHL,DO46416
2,1903276,1903276,Lymph-BNHL,DO46416
...
```

---

## Usage

1. **Prepare the Input File**  
   Ensure your file is properly formatted and includes all required columns.

2. **Run the Script**

```bash
python tumor_sample_similarity_analysis.py input_data.csv
```

3. **View the Outputs**
   - **Dendrogram**: A hierarchical clustering dendrogram will be displayed, showing the similarity between tumor samples.
   - **t-SNE Plot**: A 2D scatter plot will be displayed, with each point representing a sample colored by tumor type.

**Example:**

```bash
python tumor_sample_similarity_analysis.py input_data.csv
```

---

## Customization

You can adjust the following parameters within the script to customize the analysis:

- **Proximity Decay Factor**  
  The similarity score is affected by how close variants are on the chromosome. This is controlled by an exponential decay function:

  ```python
  np.exp(-distance / 30)
  ```

  - `30` is the default decay factor (approximate length of TF binding motifs).  
  - Increase it (e.g., `500`) to calculate similarity over larger regions such as enhancers or silencers.

- **Clustering Method**  
  Modify the method used in hierarchical clustering via `scipy.cluster.hierarchy.linkage()`:

  ```python
  linkage(matrix, method='ward')
  ```

  Supported methods include: `'ward'`, `'single'`, `'complete'`, `'average'`.

---

## Example Output

### Dendrogram
- Displays how tumor samples are grouped based on variant similarity and proximity.
- Samples that are more similar will be clustered closer together.

### t-SNE Plot
- 2D scatter plot of tumor samples.
- Each point represents a sample.
- Color-coded by tumor type.
- Spatial positioning reflects overall similarity.

---

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.
