# NEPTuner

NEPTuner (Nucleosome-Encrypted Pattern Tuner) is a tool for nucleosome profiling and pattern retrieval. First NEPTuner applies a new, simple and fast method for nucleosome profiling; then it detects and reports two types of nucleosome positioning patterns: phased pattern and nucleosome depleted region (NDR-) pattern. 

## Profiling

In simple stacking, each read gives an input that is constant along the read. When the reads are summed up, the over-estimated input of reads ends artificially gives more weight to low-coverage regions. To overcome this, many approaches use kernel functions to model reads. In our approach, we apply the simplest kernel in sake of calculation time. Reads are modeled by triangular kernel with the area of 1 and base that is equal to the read length. Thus, shorter reads produce higher triangles and longer reads give flatter ones; the input of the ends in both cases is much smaller than of the reads centers. The input of all triangles is summed up at each nucleotide to give the final profile curve. Identification of nucleosome positions is accomplished step-wise. In the first step, we describe all local maxima in a window of 40bp; this step ensures that no potential peak will be missed. In the second step, each local maximum is checked for neighboring other peaks in 150bp ( ± 75 bp) window around the peak in question. If there are additional peaks in the window, all peaks are merged in a fuzzy region. Otherwise, the peak is labeled as well positioned. 

## Pattern detection

The first kind of patterns NEPTuner identifies is *phased pattern*. This occurs when nucleosomes lie in regularly spaced arrays compatible with nucleosome width alternating with linker DNA regions.

The other pattern is related to nucleosome depleted regions and we call it NDR-pattern.
This pattern consists of an arrangement of two well-positioned nucleosomes around the transcription start site (TSS). The nucleosome upstream of the TSS is called -1, whereas the one downstream is called +1. These two nucleosomes must be well-defined. We call the space that separates them *nucleosome depleted region* if its length is greater than the usual linker DNA.

## Workflow

```{r algo, echo = F, out.width="50%", fig.show="hold"}
knitr::include_graphics("pictures/algo.png")
```

## File Formats

**Input Files**

*NEPTuner takes two input files:*

| File | Explanation |
| --- | --- |
| Experiment preprocessed file: 	| A .bed file. It comes from an MNase-seq or ChIP-seq experiment and it has already been preprocessed (quality checked, trimmed and mapped). It contains at least three columns: chromosome <string> ; start <int> ; end <int>. The reads should refer to only one chromosome and they should be sorted according to the starting coordinate |
| Annotation file: |	A .gff or a text file, with at least five columns: seqnames <string> ; start <int> ; end <int> ; strand <+ or-> ; ID <string>. The separator character is tab "\t". The contigs should be homogenoeus (e.g. mRNA or CDS). | 

**Additional Parameters**

*If the Discrete modus is chosen*

| File | Explanation |
| --- | --- |
| Threshold value: | A number between 0 and 1. The threshold will be applied to the profile signal in order to retrieve the nucleosome coordinates. |

**Output**

*There are in total four output files:*

| File | Explaination |
| --- | --- |
| Profile: |	A single column text file, each row represents the profile value at the corresponding position. |
| Nucleosome coordinates: |	start <int> ; end <int> ; width <int> |
|Phased pattern: |	start <int> ; end <int> ; width <int> |
|NDR pattern: |	This is an updated version of the annotation file. Original annotation fields ; Pattern <bool> ; NDR length <int> ; Plus1 start <int> ; Plus1 end <int> ; Minus1 start <int> ; Minus1 end <int> |
