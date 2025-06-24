# BSc RP 2425 - Dimitar Smenovski
The code in this repository is part of the 2025 CSE3000 Research Project course at TU Delft.

The thesis can be found at <a>https://resolver.tudelft.nl/uuid:0e9e95d2-f304-48e8-8b74-29ff7f7227e8</a>.

## Pipeline
The order in which to run the scripts is as follows:
1. Run combine_rosmap_data.py
2. Run combine_sea-ad_data.py

3. Run prepare_training_data_script.py on ROSMAP data
4. Run prepare_training_data_script.py on SEA-AD data
5. Run build_all_models_script.py for 21 epochs.

6. Run prepare_data_for_imputation_script.py on SEA-AD data
7. Run prepare_data_for_imputation_script.py on full ROSMAP data
8. Run prepare_data_for_imputation_script.py on olygodendroglia ROSMAP data
9. Run impute_xenium_data_script.py for each prepared dataset (from steps 6-8)

10. Run marker_genes.py
11. Run imputation_validation.py
12. Run plot_correlations.py (requires manually inserting correlation results from the previous step)

13. Run spatial_1k_script.py
14. Run classify_xenium_script.py
15. Run results_analysis_script.py with model dataset ROSMAP
16. Run results_analysis_script.py with model dataset SEA-AD
17. Run results_analysis_script.py with model dataset SEA-AD and cell type 'oligodendroglia'

## Datasets
ROSMAP [1-2]: https://www.synapse.org/Synapse:syn3219045

SEA-AD [3-4]: https://registry.opendata.aws/allen-sea-ad-atlas/ - 
Single cell profiling (transcriptomics and epigenomics) data files in a public bucket

## References
[1] D. A. Bennett, A. S. Buchman, P. A. Boyle, L. L. Barnes, R. S. Wilson, and J. A. Schneider, “Religious Orders Study and Rush Memory and Aging Project,” Journal of Alzheimer’s Disease, vol. 64, no. s1, pp. S161–S189, Jun. 2018, doi: 10.3233/JAD-179939.

[2] S. Mostafavi et al., “A molecular network of the aging human brain provides insights into the pathology and cognitive decline of Alzheimer’s disease,” Nature Neuroscience, vol. 21, no. 6, pp. 811–819, Jun. 2018, doi: 10.1038/s41593-018-0154-9.

[3] Dataset: Allen Institute for Brain Science, University of Washington Alzheimer's Disease Research Center, and Kaiser Permanente Washington Health Research Institute (2022). Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD) -- 10x single nucleus RNAseq [Dataset]. AD Knowledge Portal. Available from: https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn26223298 

[4] M. I. Gabitto et al., “Integrated multimodal cell atlas of Alzheimer’s disease,” Nature Neuroscience, vol. 27, no. 12, pp. 2366–2383, Dec. 2024, doi: 10.1038/s41593-024-01774-5.
