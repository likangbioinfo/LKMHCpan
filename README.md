# LKMHCpan2
A pan-specific peptide-MHC class I binding prediction tool
Antigens/tumor-specific antigens (TAAs/TSAs), which arising from genetic alterations, are crucial in tumor immunotherapy as both key therapeutic targets and vital molecular markers. A critical bottleneck in antigen vaccine development lies in accurately identifying a sufficient number of antigen peptides in each sample while also improving the proportion of peptides that elicit genuine immune responses. To address these limitations, a natural language processing (NLP) based model designed to predict pan-specific peptide-HLA class I binding, named LKMHCpan2, has been developed. Through the application of attention mechanisms, we have identified novel sites that greatly facilitate the representation of HLA molecules. Employing this new representation of HLA, LKMHCpan2 demonstrates significantly improved predictive performance compared to its counterparts, highlighting a substantial advancement in the field.

# Installation
We recommend using the configured Docker:

docker push liqingss/lkmhcpan:v1.0.2

# Usage
# 1) docker
python /workspace/LKMHCpan/dist/LKMHCpan2_predict.py -i peptide_HLA_test.tsv -o peptide_HLA_test_predict.tsv -b 64 -ic 0_1
# 2) singularity
Convert the docker image to a singularity image
singularity exec --nv lkmhcpan.sif python /workspace/LKMHCpan/dist/LKMHCpan2_predict.py -i peptide_HLA_test.tsv -o peptide_HLA_test_predict.tsv -b 64 -ic 0_1
