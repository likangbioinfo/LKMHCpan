
![LKMHCpanlog](https://github.com/user-attachments/assets/0947a35c-7125-4bb3-8c51-82115448f264)


# LKMHCpan
A pan-specific peptide-MHC class I binding prediction tool
Antigens/tumor-specific antigens (TAAs/TSAs), which arising from genetic alterations, are crucial in tumor immunotherapy as both key therapeutic targets and vital molecular markers. A critical bottleneck in antigen vaccine development lies in accurately identifying a sufficient number of antigen peptides in each sample while also improving the proportion of peptides that elicit genuine immune responses. To address these limitations, a natural language processing (NLP) based model designed to predict pan-specific peptide-HLA class I binding, named LKMHCpan, has been developed. Through the application of attention mechanisms, we have identified novel sites that greatly facilitate the representation of HLA molecules. Employing this new representation of HLA, LKMHCpan demonstrates significantly improved predictive performance compared to its counterparts, highlighting a substantial advancement in the field.

# Installation
**We recommend using the zenodo download LKMHCpan singularity**👍:

download LKMHCpan singularity from zenodo：[10.5281/zenodo.18081485](https://doi.org/10.5281/zenodo.18081485)

```LKMHCpan：lkmhcpan1.0.3.sif```
```Train Data：Table_MHCI_LKMHCpandata_Trim_final.tsv```
```Evaluation Data：Evaluation_data_set_from_netMHCpan_trim.csv```

# Usage
The script supports CPU or GPU, Try GPU first; if not, call CPU.

The script is invoked from the command line with various parameters that control its operation. Here is the basic syntax:

```singularity exec --nv lkmhcpan1.0.3.sif python /workspace/LKMHCpan/dist/LKMHCpan_predict.py --help```

## Options
-h, --help show this help message and exit

```-i, --input_file```: [required] Input the file to be predicted, the default format is [peptide HLA]; if the format is not the default, you need to set the [input_header] and [pep_hla_column] parameters.

```-ih, --input_header```: [default] default prediction file is without headers; if there are headers, it needs to be set to "yes".

```-ic, --pep_hla_column```: [default] Specify the column number(s) for peptide_hla in the input file to be predicted, with the default being columns 0 and 1, counting from 0.

```-b, --batch_size```: [default] batch size,default:32

```-l, --hla_list```: [optional] choose whether to input a headerless HLA list file (with each HLA subtype on a separate line). If such a file is provided, the HLA subtypes from this document will be prioritized; otherwise, the HLA subtypes from the `input_file` will be read.

```-o, --output```: [required] output file name

## Examples
### 1）Linux/macOS
```singularity exec --nv lkmhcpan1.0.3.sif python /workspace/LKMHCpan/dist/LKMHCpan_predict.py -i peptide_HLA_test.tsv -o peptide_HLA_test_predict.tsv -b 64 -ic 0_1```
