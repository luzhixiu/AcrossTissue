found datafile : E-MTAB-2812-atlasExperimentSummary.Rdata
from: https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2812/Result
	- Uses bioconductor S4Vectors and SummarizedExperiment packages
	- Had a hard time navigating and working with the object
	- Saved Rdata object which includes the useful code below.
	- Some useful code
	````r
	> if (!requireNamespace("BiocManager", quietly = TRUE))
	+     install.packages("BiocManager")
	> BiocManager::install(version = "3.10")
	> BiocManager::install("S4Vectors")
	> load("E-MTAB-2812-atlasExperimentSummary.Rdata")
	> experimentSummary
	List of length 1
	names(1): rnaseq> class(experimentSummary)
	[1] "SimpleList"
	attr(,"package")
	[1] "S4Vectors"
	> showClass("SimpleList")
	Class "SimpleList" [package "S4Vectors"]

	Slots:

	Name:           listData       elementType   elementMetadata          metadata
	Class:              list         character DataTable_OR_NULL              list

	Extends: 
	Class "List", directly
	Class "Vector", by class "List", distance 2
	Class "list_OR_List", by class "List", distance 2
	Class "Annotated", by class "List", distance 3
	Class "vector_OR_Vector", by class "List", distance 3

	Known Subclasses: 
	Class "HitsList", directly
	Class "DataFrame", directly
	Class "FilterRules", directly
	Class "SelfHitsList", by class "HitsList", distance 2
	Class "SortedByQueryHitsList", by class "HitsList", distance 2
	Class "DFrame", by class "DataFrame", distance 2
	Class "SortedByQuerySelfHitsList", by class "HitsList", distance 3
	> BiocManager::install("SummarizedExperiment")
	> experimentSummary@listData
	$rnaseq
	Loading required package: SummarizedExperiment
	Loading required package: GenomicRanges
	Loading required package: IRanges
	Loading required package: GenomeInfoDb
	Loading required package: Biobase
	Welcome to Bioconductor

			Vignettes contain introductory material; view with
		'browseVignettes()'. To cite Bioconductor, see
		'citation("Biobase")', and for packages 'citation("pkgname")'.

	Loading required package: DelayedArray
	Loading required package: matrixStats

	Attaching package: ‘matrixStats’

	The following objects are masked from ‘package:Biobase’:

			anyMissing, rowMedians

	Loading required package: BiocParallel

	Attaching package: ‘DelayedArray’

	The following objects are masked from ‘package:matrixStats’:

			colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

	The following objects are masked from ‘package:base’:

			aperm, apply, rowsum

	class: RangedSummarizedExperiment 
	dim: 46772 201 
	metadata(4): pipeline filtering mapping quantification
	assays(1): counts
	rownames(46772): WBGene00000001 WBGene00000002 ... WBGene00269437
WBGene00269438
	rowData names(0):
	colnames(201): SRR006511 SRR006514 ... SRR504339 SRR554453
	colData names(10): AtlasAssayGroup organism ... initial_time_point
technical_replicate_gr
	> experimentSummary$rnaseq
		class: RangedSummarizedExperiment 
		dim: 46772 201 
		metadata(4): pipeline filtering mapping quantification
		assays(1): counts
		rownames(46772): WBGene00000001 WBGene00000002 ... WBGene00269437
	WBGene00269438
		rowData names(0):
		colnames(201): SRR006511 SRR006514 ... SRR504339 SRR554453
		colData names(10): AtlasAssayGroup organism ... initial_time_point
	technical_replicate_group
		> colnames(experimentSummary$rnaseq)
	[1] "SRR006511" "SRR006514" "SRR006515" "SRR006516" "SRR006517" "SRR006518"
	[7] "SRR006519" "SRR006520" "SRR016669" "SRR016670" "SRR016674" "SRR016675"
	[13] "SRR016678" "SRR016690" "SRR016691" "SRR023534" "SRR023535" "SRR023536"
	[19] "SRR023537" "SRR023579" "SRR023580" "SRR023581" "SRR027187" "SRR027188"
	[25] "SRR027189" "SRR027190" "SRR027191" "SRR027192" "SRR027904" "SRR027905"
	[31] "SRR031122" "SRR031123" "SRR089352"
   ````
- Hope this is the pre-averaged data, but not sure.
- If we can get the sd for each gene in each stage, we'd be good to go.
- Added Mikes.Work to AcrossTissue repo

