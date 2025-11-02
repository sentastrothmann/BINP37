# **Deciphering myeloid-derived cell states of reprogrammed cells through single cell RNA-seq classification**
## Senta Strothmann
### BINP37: Bioinformatics - Research Project

# **Brief Introduction**
This research project aimed to perform a comparative analysis of various classification methods and their performance on datasets from (Guimarães et al., 2024) and (Miller et al., 2025), consisting of myeloid cells in cancer, which will be subclustered into more finite subsets than only their population, such as dendritic cells, monocytes and macrophages, based on the active gene programs, genes and their respective subsets. 
The following README provides an overview of the scripts and their functionality.

# **Setup**
## Folder Structure
+ Create the folder structure used in the scripts
    + Manually add the right scripts into the right folder, as well as the input data
```bash
# The folder structure can be created with create_folder_structure.sh

# Make sure that the script is executable
chmod +x create_folder_structure.sh

# Run the script
./create_folder_structure.sh
```

## Conda Environment: scseq
+ Set up the scseq environment using conda
+ Automatically installs all tools in their right version
    + The full list of versions can be found at the end of the README file
```bash
# Create the environment
conda env create -f scseq_environment.yml

# Activate the environment
conda activate scseq
```

# **Data**
**NOTE:**
Due to the large size of the data input files (3 GB and 16GB), they are not uploaded to this GitHub Repository
+ scAtlas.h5ad
    + https://datasets.cellxgene.cziscience.com/07405240-4f64-4dd9-83c3-b3db3405b05c.h5ad
        + a multi-tissue single-cell tumor microenvironment atlas
        + Guimarães, G. R., Maklouf, G. R., Teixeira, C. E., de Oliveira Santos, L., Tessarollo, N. G., de Toledo, N. E., ... &  Boroni, M. (2024). Single-cell resolution characterization of myeloid-derived cell states with implication in cancer outcome. Nature Communications, 15(1), 5694.
        + This atlas is an integration of 13 single-cell studies of cells obtained from 8 tumor types and normal tissues, including breast, colorectal, ovary, lung, liver, skin, uvea and PBMC. Four types of samples compose this dataset (normal, tumor, lymph node, and blood) obtained from 3 different types of single-cell technologies (10x, Smart-seq2, and inDrop). About 400 thousand cells were carefully integrated and most cell types annotated in detail, specially myeloid-derived cells and states.
+ adata.h5ad
    + https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003756.v1.p1
        + Miller, T. E., El Farran, C. A., Couturier, C. P., Chen, Z., D’Antonio, J. P., Verga, J., Villanueva, M. A., Gonzalez Castro, L. N., Tong, Y. E., Saadi, T. A., Chiocca, A. N., Zhang, Y., Fischer, D. S., Heiland, D. H., Guerriero, J. L., Petrecca, K., Suva, M. L., Shalek, A. K., & Bernstein, B. E. (2025). Programs, origins and immunomodulatory functions of myeloid cells in glioma. Nature, 1–11. https://doi.org/10.1038/s41586-025-08633-8

# **Brief Description of Each Step**
## Quality Control
Scripts:
+ quality_control.py

Description:
+ The quality control script is used for initial filtering and quality assessment of the scRNA-seq dataset. Metrics such as gene detecion, mitochondrial/ribosomal/hemoglobin content and doublet prediction are computed. Low-quality cells are removed, highly variable genes are selected and PCA and UMAP embeddings are generated for visualization of the cell types. 

## Dimensionality Reduction
Scripts:
+ dimensionality_reduction.py
+ dimensionality_reduction.ipynb

Description:
+ The dimensionality reduction step reduces the dimensionality of the dataset for an improved downstream analysis. PCA, t-SNE and UMAP embeddings are computed for visualization of cellular heterogeneity. Doublets are re-checked, the data is normalized and known immune marker genes are validated.

## Data Integration
Scripts:
+ data_integration.ipynb

Description:
+ The step of data integration correct batch effects and integrates the dataset using Harmony, BBKNN, Scanorama and ComBat, further comparing them. This generated integrated embeddings and batch effects are visualized to show consistency across cell types.

## Clustering
Scripts:
+ clustering.ipynb

Description:
+ Clusters of cell types are built using graph-based Leiden, k-means and hierarchical clustering on the integrated PCA embeddings from Harmony. The cluster relationships and quality are visualized with dendrograms, UMAPs and marker gene expression plots.

## Differential Gene Expression
Scripts:
+ differential_gene_expression.ipynb

Description:
+ Highyl variable genes are extracted and cluster-specific marker genes are identified by multiple statistical tests (t-test. t-test with overestimate variance, Wilcoxon rank-sum). The results are shown by using heatmaps, dot plots and violin plots to evaluate intra- and inter-cluster gene expression differences. 

## Cell Type Prediction
Scripts:
+ cell_type_prediction.ipynb

Description:
+ Cell type labels are assigned by integrating reference datasets and applying Scanorama and CellTypist. The cell type distributions and their cluster correspondence are shown through UMAPs and bar plots. 

## Trajectory Inference
Scripts:
+ trajectory_inference.ipynb

Description:
+ Trajectory inference with PAGA and diffusion pseudtome are computed and visualized for differentiation paths of myeloid, lympoid and stromal cells. Regulateroy programs are mapped, as well as the gene expression levels along trajectories.

## Classifiers
Scripts:
+ Initial Training
    + scAtlas_lr_classifier.py
    + scAtlas_rf_classifier.py
    + scAtlas_svm_classifier.py
    + scAtlas_xgb_classifier.py
+ Preprocessing
    + adata_clustering.ipynb
    + scAtlas_match.ipynb
+ Training
    + lr_filtered.ipynb
    + rf_filtered.ipynb
    + svm_filtered.ipynb
    + xgb_filtered.ipynb
+ Validation
    + combine_pred_adata.ipynb
    + lr_validation.ipynb
    + rf_validation.ipynb
    + svm_validation.ipynb
    + xgb_validation.ipynb

Description:
+ In the final step of building classifiers, four machine learning-based classifiers (logistic regression, random forest, support vector machines and XGBoost) are built and validated for cell type prediction. Initial classifers are trained on the full Guimarães-dataset and in a second attempt, a filtered Miller-dataset using makrer genes for myeloid cells are used for improved performance. The training involves hyperparamter optimization, cross-validation and stratified splitting of the data. The validation compares the predicted labels with the true labels, and performance metrics (accuracy, F1-scores, confusion matrices, heatmap) are computed. 

# **Full list of all versions**
```bash
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       2_gnu    conda-forge
anndata                   0.11.3             pyhd8ed1ab_0    conda-forge
anyio                     4.11.0             pyhcf101f3_0    conda-forge
aom                       3.9.1                hac33072_0    conda-forge
argon2-cffi               25.1.0             pyhd8ed1ab_0    conda-forge
argon2-cffi-bindings      25.1.0          py310h7c4b9e2_1    conda-forge
arpack                    3.9.1           nompi_hf03ea27_102    conda-forge
array-api-compat          1.11.2             pyh29332c3_0    conda-forge
arrow                     1.4.0              pyhcf101f3_0    conda-forge
asttokens                 3.0.0              pyhd8ed1ab_1    conda-forge
async-lru                 2.0.5              pyh29332c3_0    conda-forge
attrs                     25.4.0             pyh71513ae_0    conda-forge
babel                     2.17.0             pyhd8ed1ab_0    conda-forge
bbknn                     1.6.0           py310h1fe012e_4    bioconda
beautifulsoup4            4.14.2             pyha770c72_0    conda-forge
bleach                    6.2.0              pyh29332c3_4    conda-forge
bleach-with-css           6.2.0                h82add2a_4    conda-forge
blosc                     1.21.6               he440d0b_1    conda-forge
brotli                    1.1.0                hb9d3cd8_2    conda-forge
brotli-bin                1.1.0                hb9d3cd8_2    conda-forge
brotli-python             1.1.0           py310hf71b8c6_2    conda-forge
brunsli                   0.1                  h9c3ff4c_0    conda-forge
bzip2                     1.0.8                h4bc722e_7    conda-forge
c-ares                    1.34.4               hb9d3cd8_0    conda-forge
c-blosc2                  2.15.2               h3122c55_1    conda-forge
ca-certificates           2025.10.5            hbd8a1cb_0    conda-forge
cached-property           1.5.2                hd8ed1ab_1    conda-forge
cached_property           1.5.2              pyha770c72_1    conda-forge
celltypist                1.6.3              pyhdfd78af_0    bioconda
certifi                   2025.10.5          pyhd8ed1ab_0    conda-forge
cffi                      1.17.1          py310h8deb56e_0    conda-forge
charls                    2.4.2                h59595ed_0    conda-forge
charset-normalizer        3.4.1              pyhd8ed1ab_0    conda-forge
click                     8.1.8              pyh707e725_0    conda-forge
colorama                  0.4.6              pyhd8ed1ab_1    conda-forge
comm                      0.2.2              pyhd8ed1ab_1    conda-forge
contourpy                 1.3.1           py310h3788b33_0    conda-forge
cycler                    0.12.1             pyhd8ed1ab_1    conda-forge
dav1d                     1.2.1                hd590300_0    conda-forge
debugpy                   1.8.13          py310hf71b8c6_0    conda-forge
decorator                 5.2.1              pyhd8ed1ab_0    conda-forge
defusedxml                0.7.1              pyhd8ed1ab_0    conda-forge
entrypoints               0.4                pyhd8ed1ab_1    conda-forge
et_xmlfile                2.0.0              pyhd8ed1ab_1    conda-forge
exceptiongroup            1.2.2              pyhd8ed1ab_1    conda-forge
executing                 2.1.0              pyhd8ed1ab_1    conda-forge
fa2-modified              0.3.10                   pypi_0    pypi
fbpca                     1.0                pyhd8ed1ab_1    conda-forge
fonttools                 4.56.0          py310h89163eb_0    conda-forge
fqdn                      1.5.1              pyhd8ed1ab_1    conda-forge
freetype                  2.13.3               h48d6fc4_0    conda-forge
geosketch                 1.3                pyhdfd78af_0    bioconda
giflib                    5.2.2                hd590300_0    conda-forge
glpk                      5.0                  h445213a_0    conda-forge
gmp                       6.3.0                hac33072_2    conda-forge
gseapy                    1.1.8           py310hec43fc7_0    bioconda
h11                       0.16.0             pyhd8ed1ab_0    conda-forge
h2                        4.2.0              pyhd8ed1ab_0    conda-forge
h5py                      3.13.0          nompi_py310h60e0fe6_100    conda-forge
harmonypy                 0.0.10             pyhdfd78af_0    bioconda
hdf5                      1.14.3          nompi_h2d575fe_109    conda-forge
hpack                     4.1.0              pyhd8ed1ab_0    conda-forge
httpcore                  1.0.9              pyh29332c3_0    conda-forge
httpx                     0.28.1             pyhd8ed1ab_0    conda-forge
hyperframe                6.1.0              pyhd8ed1ab_0    conda-forge
icu                       75.1                 he02047a_0    conda-forge
idna                      3.10               pyhd8ed1ab_1    conda-forge
igraph                    0.10.15              he44f51b_1    conda-forge
imagecodecs               2024.12.30      py310h78a9a29_0    conda-forge
imageio                   2.37.0             pyhfb79c49_0    conda-forge
importlib-metadata        8.6.1              pyha770c72_0    conda-forge
intervaltree              3.1.0              pyhd8ed1ab_1    conda-forge
ipykernel                 6.29.5             pyh3099207_0    conda-forge
ipython                   8.34.0             pyh907856f_0    conda-forge
isoduration               20.11.0            pyhd8ed1ab_1    conda-forge
jedi                      0.19.2             pyhd8ed1ab_1    conda-forge
jinja2                    3.1.6              pyhd8ed1ab_0    conda-forge
joblib                    1.4.2              pyhd8ed1ab_1    conda-forge
json5                     0.12.1             pyhd8ed1ab_0    conda-forge
jsonpointer               3.0.0           py310hff52083_2    conda-forge
jsonschema                4.25.1             pyhe01879c_0    conda-forge
jsonschema-specifications 2025.9.1           pyhcf101f3_0    conda-forge
jsonschema-with-format-nongpl 4.25.1               he01879c_0    conda-forge
jupyter-lsp               2.3.0              pyhcf101f3_0    conda-forge
jupyter_client            7.4.9              pyhd8ed1ab_0    conda-forge
jupyter_core              5.7.2              pyh31011fe_1    conda-forge
jupyter_events            0.12.0             pyh29332c3_0    conda-forge
jupyter_server            2.17.0             pyhcf101f3_0    conda-forge
jupyter_server_terminals  0.5.3              pyhd8ed1ab_1    conda-forge
jupyterlab                4.4.10             pyhd8ed1ab_0    conda-forge
jupyterlab_pygments       0.3.0              pyhd8ed1ab_2    conda-forge
jupyterlab_server         2.28.0             pyhcf101f3_0    conda-forge
jxrlib                    1.1                  hd590300_3    conda-forge
keyutils                  1.6.1                h166bdaf_0    conda-forge
kiwisolver                1.4.7           py310h3788b33_0    conda-forge
krb5                      1.21.3               h659f571_0    conda-forge
lark                      1.3.1              pyhd8ed1ab_0    conda-forge
lazy-loader               0.4                pyhd8ed1ab_2    conda-forge
lazy_loader               0.4                pyhd8ed1ab_2    conda-forge
lcms2                     2.17                 h717163a_0    conda-forge
ld_impl_linux-64          2.43                 h712a8e2_4    conda-forge
legacy-api-wrap           1.4.1              pyhd8ed1ab_0    conda-forge
leidenalg                 0.10.2          py310hc6cd4ac_0    conda-forge
lerc                      4.0.0                h27087fc_0    conda-forge
libaec                    1.1.3                h59595ed_0    conda-forge
libavif16                 1.2.1                h63b8bd6_0    conda-forge
libblas                   3.9.0           20_linux64_openblas    conda-forge
libbrotlicommon           1.1.0                hb9d3cd8_2    conda-forge
libbrotlidec              1.1.0                hb9d3cd8_2    conda-forge
libbrotlienc              1.1.0                hb9d3cd8_2    conda-forge
libcblas                  3.9.0           20_linux64_openblas    conda-forge
libcurl                   8.12.1               h332b0f4_0    conda-forge
libdeflate                1.23                 h4ddbbb0_0    conda-forge
libedit                   3.1.20250104    pl5321h7949ede_0    conda-forge
libev                     4.33                 hd590300_2    conda-forge
libffi                    3.4.6                h2dba641_0    conda-forge
libgcc                    15.1.0               h767d61c_4    conda-forge
libgcc-ng                 15.1.0               h69a702a_4    conda-forge
libgfortran               15.1.0               h69a702a_4    conda-forge
libgfortran-ng            15.1.0               h69a702a_4    conda-forge
libgfortran5              15.1.0               hcea5267_4    conda-forge
libgomp                   15.1.0               h767d61c_4    conda-forge
libhwloc                  2.11.2          default_h0d58e46_1001    conda-forge
libhwy                    1.1.0                h00ab1b0_0    conda-forge
libiconv                  1.18                 h4ce23a2_1    conda-forge
libjpeg-turbo             3.0.0                hd590300_1    conda-forge
libjxl                    0.11.1               hdb8da77_0    conda-forge
liblapack                 3.9.0           20_linux64_openblas    conda-forge
libleidenalg              0.11.1               h00ab1b0_0    conda-forge
liblzma                   5.6.4                hb9d3cd8_0    conda-forge
libnghttp2                1.64.0               h161d5f1_0    conda-forge
libnsl                    2.0.1                hd590300_0    conda-forge
libopenblas               0.3.25          pthreads_h413a1c8_0    conda-forge
libpng                    1.6.47               h943b412_0    conda-forge
libsodium                 1.0.18               h36c2ea0_1    conda-forge
libsqlite                 3.49.1               hee588c1_2    conda-forge
libssh2                   1.11.1               hf672d98_0    conda-forge
libstdcxx                 15.1.0               h8f9b012_4    conda-forge
libstdcxx-ng              15.1.0               h4852527_4    conda-forge
libtiff                   4.7.0                hd9ff511_3    conda-forge
libuuid                   2.38.1               h0b41bf4_0    conda-forge
libwebp-base              1.5.0                h851e524_0    conda-forge
libxcb                    1.17.0               h8a09558_0    conda-forge
libxcrypt                 4.4.36               hd590300_1    conda-forge
libxgboost                3.0.0            cpu_h97e353d_1    conda-forge
libxml2                   2.13.6               h8d12d68_0    conda-forge
libzlib                   1.3.1                hb9d3cd8_2    conda-forge
libzopfli                 1.0.3                h9c3ff4c_0    conda-forge
llvmlite                  0.44.0          py310h1a6248f_1    conda-forge
louvain                   0.8.2           py310hf71b8c6_1    conda-forge
lz4-c                     1.10.0               h5888daf_1    conda-forge
markupsafe                3.0.3           py310h3406613_0    conda-forge
matplotlib-base           3.10.1          py310h68603db_0    conda-forge
matplotlib-inline         0.1.7              pyhd8ed1ab_1    conda-forge
matplotlib-venn           1.1.2              pyhd8ed1ab_0    conda-forge
mistune                   3.1.4              pyhcf101f3_0    conda-forge
munkres                   1.1.4              pyh9f0ad1d_0    conda-forge
natsort                   8.4.0              pyh29332c3_1    conda-forge
nbclient                  0.10.2             pyhd8ed1ab_0    conda-forge
nbconvert-core            7.16.6             pyhcf101f3_1    conda-forge
nbformat                  5.10.4             pyhd8ed1ab_1    conda-forge
ncurses                   6.5                  h2d0b736_3    conda-forge
nest-asyncio              1.6.0              pyhd8ed1ab_1    conda-forge
networkx                  3.4.2              pyh267e887_2    conda-forge
notebook                  7.4.7              pyhd8ed1ab_0    conda-forge
notebook-shim             0.2.4              pyhd8ed1ab_1    conda-forge
numba                     0.61.0          py310h699fe88_1    conda-forge
numpy                     1.26.4          py310hb13e2d6_0    conda-forge
openjpeg                  2.5.3                h5fbd93e_0    conda-forge
openpyxl                  3.1.5           py310h0999ad4_1    conda-forge
openssl                   3.5.4                h26f9b46_0    conda-forge
overrides                 7.7.0              pyhd8ed1ab_1    conda-forge
packaging                 24.2               pyhd8ed1ab_2    conda-forge
pandas                    2.2.3           py310h5eaa309_1    conda-forge
pandocfilters             1.5.0              pyhd8ed1ab_0    conda-forge
parso                     0.8.4              pyhd8ed1ab_1    conda-forge
patsy                     1.0.1              pyhd8ed1ab_1    conda-forge
pexpect                   4.9.0              pyhd8ed1ab_1    conda-forge
pickleshare               0.7.5           pyhd8ed1ab_1004    conda-forge
pillow                    11.1.0          py310h7e6dc6c_0    conda-forge
pip                       25.0.1             pyh8b19718_0    conda-forge
platformdirs              4.3.7              pyh29332c3_0    conda-forge
prometheus_client         0.23.1             pyhd8ed1ab_0    conda-forge
prompt-toolkit            3.0.50             pyha770c72_0    conda-forge
psutil                    7.0.0           py310ha75aee5_0    conda-forge
pthread-stubs             0.4               hb9d3cd8_1002    conda-forge
ptyprocess                0.7.0              pyhd8ed1ab_1    conda-forge
pure_eval                 0.2.3              pyhd8ed1ab_1    conda-forge
py-xgboost                3.0.0           cpu_pyh1ce2f49_1    conda-forge
pycparser                 2.22               pyh29332c3_1    conda-forge
pygments                  2.19.1             pyhd8ed1ab_0    conda-forge
pynndescent               0.5.13             pyhd8ed1ab_1    conda-forge
pyopengl                  3.1.7              pyhd8ed1ab_0    conda-forge
pyparsing                 3.2.2              pyhd8ed1ab_0    conda-forge
pysocks                   1.7.1              pyha55dd90_7    conda-forge
python                    3.10.16         he725a3c_1_cpython    conda-forge
python-annoy              1.17.3          py310hf71b8c6_2    conda-forge
python-dateutil           2.9.0.post0        pyhff2d567_1    conda-forge
python-fastjsonschema     2.21.2             pyhe01879c_0    conda-forge
python-igraph             0.11.8          py310h5996f03_0    conda-forge
python-json-logger        2.0.7              pyhd8ed1ab_0    conda-forge
python-tzdata             2025.2             pyhd8ed1ab_0    conda-forge
python_abi                3.10                    5_cp310    conda-forge
pytz                      2024.1             pyhd8ed1ab_0    conda-forge
pywavelets                1.8.0           py310hf462985_0    conda-forge
pyyaml                    6.0.3           py310h3406613_0    conda-forge
pyzmq                     24.0.1          py310h330234f_1    conda-forge
qhull                     2020.2               h434a139_5    conda-forge
rav1e                     0.6.6                he8a937b_2    conda-forge
readline                  8.2                  h8c095d6_2    conda-forge
referencing               0.37.0             pyhcf101f3_0    conda-forge
requests                  2.32.3             pyhd8ed1ab_1    conda-forge
rfc3339-validator         0.1.4              pyhd8ed1ab_1    conda-forge
rfc3986-validator         0.1.1              pyh9f0ad1d_0    conda-forge
rfc3987-syntax            1.1.0              pyhe01879c_1    conda-forge
rpds-py                   0.28.0          py310hd8f68c5_1    conda-forge
scanorama                 1.7.4              pyhdfd78af_0    bioconda
scanpy                    1.11.0             pyhd8ed1ab_0    conda-forge
scikit-image              0.25.2          py310h5eaa309_0    conda-forge
scikit-learn              1.5.2           py310h27f47ee_1    conda-forge
scikit-misc               0.5.1           py310h18bf74c_2    conda-forge
scipy                     1.11.4          py310hb13e2d6_0    conda-forge
scrublet                  0.2.3              pyh5e36f6f_1    bioconda
seaborn                   0.13.2               hd8ed1ab_3    conda-forge
seaborn-base              0.13.2             pyhd8ed1ab_3    conda-forge
send2trash                1.8.3              pyh0d859eb_1    conda-forge
session-info2             0.1.2              pyhd8ed1ab_0    conda-forge
setuptools                75.8.2             pyhff2d567_0    conda-forge
six                       1.17.0             pyhd8ed1ab_0    conda-forge
snappy                    1.2.1                h8bd8927_1    conda-forge
sniffio                   1.3.1              pyhd8ed1ab_1    conda-forge
sortedcontainers          2.4.0              pyhd8ed1ab_1    conda-forge
soupsieve                 2.8                pyhd8ed1ab_0    conda-forge
stack_data                0.6.3              pyhd8ed1ab_1    conda-forge
statsmodels               0.14.4          py310hf462985_0    conda-forge
svt-av1                   3.0.1                h5888daf_0    conda-forge
tbb                       2022.0.0             hceb3a55_0    conda-forge
terminado                 0.18.1             pyh0d859eb_0    conda-forge
texttable                 1.7.0              pyhd8ed1ab_1    conda-forge
threadpoolctl             3.6.0              pyhecae5ae_0    conda-forge
tifffile                  2025.3.13          pyhd8ed1ab_0    conda-forge
tinycss2                  1.4.0              pyhd8ed1ab_0    conda-forge
tk                        8.6.13          noxft_h4845f30_101    conda-forge
tomli                     2.3.0              pyhcf101f3_0    conda-forge
tornado                   6.4.2           py310ha75aee5_0    conda-forge
tqdm                      4.67.1             pyhd8ed1ab_1    conda-forge
traitlets                 5.14.3             pyhd8ed1ab_1    conda-forge
typing-extensions         4.12.2               hd8ed1ab_1    conda-forge
typing_extensions         4.12.2             pyha770c72_1    conda-forge
typing_utils              0.1.0              pyhd8ed1ab_1    conda-forge
tzdata                    2025b                h78e105d_0    conda-forge
umap-learn                0.5.7           py310hff52083_1    conda-forge
unicodedata2              16.0.0          py310ha75aee5_0    conda-forge
uri-template              1.3.0              pyhd8ed1ab_1    conda-forge
urllib3                   2.3.0              pyhd8ed1ab_0    conda-forge
wcwidth                   0.2.13             pyhd8ed1ab_1    conda-forge
webcolors                 25.10.0            pyhd8ed1ab_0    conda-forge
webencodings              0.5.1              pyhd8ed1ab_3    conda-forge
websocket-client          1.9.0              pyhd8ed1ab_0    conda-forge
wheel                     0.45.1             pyhd8ed1ab_1    conda-forge
xorg-libxau               1.0.12               hb9d3cd8_0    conda-forge
xorg-libxdmcp             1.1.5                hb9d3cd8_0    conda-forge
yaml                      0.2.5                h280c20c_3    conda-forge
zeromq                    4.3.5                h75354e8_4    conda-forge
zfp                       1.0.1                h5888daf_2    conda-forge
zipp                      3.21.0             pyhd8ed1ab_1    conda-forge
zlib-ng                   2.2.4                h7955e40_0    conda-forge
zstandard                 0.23.0          py310ha75aee5_1    conda-forge
zstd                      1.5.7                hb8e6e7a_2    conda-forge
```
