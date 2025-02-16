# User-specific settings
####################################################################################

# Inputs/Outputs
export GSNAMEBASE='genedex' # Base name of the gene set, which can be custom or one of: genedex, MSigDBv7, MSigDBv6
export DSNAME='gtextoil_iBrain' # name of the background reference dataset: gtextoil, gtextoil_iBrain, HGU133plus2
export OUTDIR_NAME='test_run' # name of directory to store outputs. A folder with this name will be created inside ../output
export COLLECTION='' # Optional, name of collection used to augment the custom geneset.
export NWRK_TYPE='' # Optional, type of network to calculate between custom gene set and collection: unipartite, bipartite, or uni_bipartite

# Resources
export CORES=10 # number of cores (for part 1)
export CORES_P2=25 # number of cores (for part 2)

# PDxN parameters
export MAX_GENES=500 # Largest pathway size allowed
export MIN_GENES=20 # Smallest pathway size allowed 
export MAX_JACQ=85 # Max jacquard index between two pathways
export PCOR_OPTION=0 # choice of partial correlation
export MIN_SAMPLES=10 # min number of samples per tissue

# System settings - DO NOT CHANGE!
####################################################################################
export STD_GSET_TABLE='../input/std_gene_tables/'${GSNAMEBASE}'_pathway_table.csv' # Needed if running 00_prepset_wrapper.sh
export OUTDIR="../output/${OUTDIR_NAME}" # directory to store outputs
export DSETDIR="../input/gene_expression" # directory containing background dataset
export GSETDIR="../input/gene_sets" # directory containing gene set files
export ASETDIR="../input/augment_sets" # directory with standard tables of gene sets used in collections
export GENEUNIV="${DSETDIR}/${DSNAME}/genes_${DSNAME}.csv" # Gene universe file generated automatically when processing background gene expression
export GSNAME=${GSNAMEBASE}__${DSNAME} # Geneset pre-processed for the corresponding dataset

# Create run-specific directories
mkdir -p ${OUTDIR} ${GSETDIR}

# Modify GSNAME when augmenting gene set with collection 

if [[ -n "$COLLECTION" && -n "$NWRK_TYPE" ]]
then
    GSNAME=${GSNAMEBASE}__${COLLECTION}__${NWRK_TYPE}__${DSNAME}
fi
