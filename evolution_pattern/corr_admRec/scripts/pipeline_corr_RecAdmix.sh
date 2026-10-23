##README:
#this project is for comparing the correlation between admixture and recombination rate
#For each chromosome (e.g., chr1), the original raw data is locate in: ${MEXICANA_INTROGRESSION_ROOT}/0_elai/out_elai_845maize/chr1/mg1000/s10w10/
#The original raw data for each sample (e.g., CIMBL11.amp) is: ${MEXICANA_INTROGRESSION_ROOT}/0_elai/out_elai_845maize/chr1/mg1000/s10w10/CIMBL11.amp.chr1.mg1000.s10w10.ave10runs.admix.txt
#the coordinates is based on B73_v4

## Plotting the correlation scatter between recombination rate (cm/Mb) and admixture dosage
## 1. get the admMap with recRate (cm/Mb)
bsub < addRec_to_admMap_perChr.lsf
## 2. plot the decile box and correlation scatter between rec VS admix
### 2.1 For average admixture across the whole population
bsub < plot_corr_recVSadmix_allChrs.lsf

### 2.2 For each sample (need to update using new rec_map from Jeff)
bsub < correlation_decileScatter_recVSadmix_allChrs.lsf

#### plot the distribution of spearman correlation index for all samples
dir_viz="${MEXICANA_INTROGRESSION_ROOT}/0_elai/out_elai_845maize/visualization/s10w10_mg1000_unmask"
conda activate R4.4_env
Rscript plot_dist_corrRecVSadmix.R $dir_viz/spearman_result_allsamples.txt $dir_viz/spearman_result_allsamples.pdf
