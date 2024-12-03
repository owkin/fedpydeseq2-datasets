log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(snapcount)
library(recount3)
library(data.table)

download_data = function(project, project_info){
    # Download project data
    exp = recount3::create_rse(
        project_info,
        annotation='gencode_v29',
        type='gene',
    )
    # Save metadata
    meta = as.data.frame(colData(exp))
    fwrite(
        meta,
        file=file.path("results", project, "metadata.tsv.gz"),
        sep="\t",
    )

    # Scale coverage counts into read counts
    read_counts = as.data.frame(recount3::compute_read_counts(exp))
    assays(exp)$counts = recount3::transform_counts(exp)

    read_counts$gene_id = rownames(read_counts)
    read_counts = read_counts[,c('gene_id', colnames(read_counts)[!colnames(read_counts) %in% c('gene_id')])]

    # Save expression data
    fwrite(
        read_counts,
        file=file.path("results", project, paste0("Counts_raw", ".tsv.gz")),
        sep="\t",
    )
    # Save gene names
    fwrite(
        data.frame(
            gene_name=rowData(exp)[['gene_name']],
            gene_id=rowData(exp)[['gene_id']]
        ),
        file=file.path("results", project, paste0("gene_names.tsv.gz")),
        sep="\t",
    )
}

split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))
dataset = split_path(snakemake@output[[1]])[2]

human_projects = recount3::available_projects()
project_info <- subset(
    human_projects,
    project == toupper(dataset),
)
print(project_info)

download_data(dataset, project_info)
