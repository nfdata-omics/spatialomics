# %% [markdown]
# # Setup: import modules & read data

# %% [markdown]
# # From image segmentation to Bin2cell

# %%
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import os
import seaborn as sns
import matplotlib as mpl

#Segmentation
import geopandas as gpd
import bin2cell as b2c

# Image
from skimage.segmentation import find_boundaries
import scipy.sparse
import tifffile
import czifile
import cv2

import igraph

# %%
# random seed for reproducibility
seed = 42

# %% [markdown]
# # Save mask as sparse NPZ

# %%
# Load original mask
# Mask is the result of segmentation, it should be converted into sparse npz file before being used for bin2cell
mask = tifffile.imread('/home/jovyan/work/custom/utils/mask/NFG021_002_masks.tif')
print(f"Mask: {mask.shape}, dtype: {mask.dtype}")
print(f"Cells: {len(np.unique(mask)) - 1}")

# Convert to sparse matrix
mask_sparse = scipy.sparse.csr_matrix(mask)

# Save as NPZ
scipy.sparse.save_npz('/home/jovyan/work/custom/utils/mask/NFG021_002_ROI_20X_num2_masks_sparse.npz', mask_sparse)

print("Mask saved as sparse NPZ!")

# %% [markdown]
# # Load results of Space Ranger

# %%
# Path to Space Ranger output
path = "/home/jovyan/work/results/alignment_resize/NFG021_005_CyLIB/outs/binned_outputs/square_002um"
source_image_path = "/home/jovyan/work/results/bin2cell/NFG021_005_halfres_for_spaceranger_uint8.tif"
spaceranger_image_path = "/home/jovyan/work/results/alignment_resize/NFG021_005_CyLIB/outs/binned_outputs/square_002um/spatial"

# %%
# Load the Visium data using Bin2Cell
adata = b2c.read_visium(
    path,
    source_image_path=source_image_path,
    spaceranger_image_path=spaceranger_image_path
)

# %% [markdown]
# # Insert labels

# %%
# Insert labels
b2c.insert_labels(
    adata=adata,
    labels_npz_path='/home/jovyan/work/custom/utils/mask/NFG021_005_ROI_20X_num2_masks_sparse.npz',
    basis="spatial",
    spatial_key="spatial",
    mpp=None,
    labels_key="labels"
)

# %%
# Statistics
print(f"\Assignment statistics:")
print(f"   Total bins: {adata.n_obs}")
print(f"   Bins assigned to cells: {np.sum(adata.obs['labels'] > 0)}")
print(f"   Unassigned bins: {np.sum(adata.obs['labels'] == 0)}")
print(f"   Cells detected: {len(np.unique(adata.obs['labels'])) - 1}")

# %%
adata.obs

# %%
# Expand labels to include also the cytoplasm, decide to use or not, decide if use the volume ration or a costant expansion
b2c.expand_labels(adata,
                  labels_key='labels',
                  expanded_labels_key="labels_expanded",
                  volume_ratio = 4
                 )

# %%
# Aggregate bin-level data to cell-level using the provided labels and spatial coordinates
adata_cells = b2c.bin_to_cell(
    adata=adata,
    labels_key="labels_expanded",
    spatial_keys=["spatial"],
    diameter_scale_factor=None
)

# %%
adata_cells.obs

# %%
# Calculate total_counts manually
adata_cells.obs['total_counts'] = np.array(adata_cells.X.sum(axis=1)).flatten()

# Calculate n_genes manually
adata_cells.obs['n_genes'] = np.array((adata_cells.X > 0).sum(axis=1)).flatten()


print(f"   Before: {adata.shape[0]} bins x {adata.shape[1]} genes")
print(f"   After:  {adata_cells.shape[0]} cells x {adata_cells.shape[1]} genes")
print(f"\n   Average bins per cell: {adata_cells.obs['bin_count'].mean():.1f}")
print(f"   Average counts per cell: {adata_cells.obs['total_counts'].mean():.0f}")
print(f"   Average genes per cell: {adata_cells.obs['n_genes'].mean():.0f}")

# %%
mean_bin_count = adata_cells.obs['bin_count'].mean()
median_bin_count = adata_cells.obs['bin_count'].median()

print(f"Average bin_count per cell: {mean_bin_count:.1f}")
print(f"Median bin_count per cell: {median_bin_count:.1f}")

# Plot distribuzione
plt.figure(figsize=(8,5))
plt.hist(adata_cells.obs['bin_count'], bins=50, color='skyblue', edgecolor='black')
plt.axvline(mean_bin_count, color='red', linestyle='--', label=f'Mean: {mean_bin_count:.1f}')
plt.axvline(median_bin_count, color='green', linestyle='--', label=f'Median: {median_bin_count:.1f}')
plt.xlabel('Bin count per cell')
plt.ylabel('Number of cells')
plt.title('Distribution of bin_count per cell')
plt.legend()
plt.show()

# %%
# Save the final result
adata_cells.write_h5ad('/home/jovyan/work/results/object/object_bin2cell/NFG021_005_cells.h5ad')

print("File saved!")
print("Cells x genes matrix!")
