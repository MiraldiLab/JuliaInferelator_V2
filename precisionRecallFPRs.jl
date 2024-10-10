using CSV
using DataFrames
using DelimitedFiles
using Printf
using ColorTypes, ColorSchemes
using Plots


# Inputs and data
output_dir = "/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/MichaelTest"
fig_save_name = " "
gold_standard_path = "/data/miraldiNB/wayman/projects/Tfh10/outs/202112/GRN/inputs/GS/RNA/priors/Log2FC0p5_FDR20_Rank50/prior_RNA_Thelper_Miraldi2019Th17_combine_Log2FC0p5_FDR20_Rank50_Frob_sp.tsv"
common_tfs_path = nothing  # Set to nothing
common_targets_path = nothing  # Set to nothing

inf_grn_paths = [
    "/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs2/Th17_50ss_Combined/combined_50.tsv",
    "/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs2/Th17_100ss_Combined/combined_100.tsv"
]

# Load gold standard data
gold_standard = CSV.read(gold_standard_path, DataFrame; delim = '\t');
# Check if common_tfs_path and common_targets_path exist and not null and read the file
common_tfs = nothing
common_targets = nothing
if !isnothing(common_tfs_path)
    common_tfs = CSV.File(common_tfs_path; delim='\n', header=false)[:, 1]
    common_tfs = collect(common_tfs)  # Convert to Vector
end

if !isnothing(common_targets_path)
    common_targets = CSV.File(common_targets_path; delim='\n', header=false)[:, 1]
    common_targets = collect(common_targets)  # Convert to Vector
end

color_generator = [RGB(rand(), rand(), rand()) for _ in 1:length(inf_grn_paths)]  # Generate random colors for plots

function filter_common_genes(gold_standard, inferred_grn; common_targets=nothing, common_tfs=nothing)

    # 1. Find all genes in common between goldstandard (gs) and inferred_grn
    if isnothing(common_targets)
        # Default behavior: compute common targets
        common_targets = intersect(gold_standard.Target, inferred_grn.Gene)
    end

    # Filter based on common targets
    filtered_gs = filter(row -> row.Target in common_targets, gold_standard)
    filtered_inf = filter(row -> row.Gene in common_targets, inferred_grn)
    
    # Keep only TF, Target/Gene columns
    gs_tfgene = filtered_gs[:, [:TF, :Target]]
    inf_tfgene = filtered_inf[:, [:TF, :Gene]]

    # 2. Keep data points corresponding to shared TFs between GS and inferred GRN
    if isnothing(common_tfs)
        # Default behavior: compute common TFs
        common_tfs = intersect(gs_tfgene.TF, inf_tfgene.TF)
    end

    len_common_tfs = length(common_tfs)

    gs_tfgn = filter(row -> row.TF in common_tfs, gs_tfgene)
    inf_tfgn = filter(row -> row.TF in common_tfs, inf_tfgene)

    return gs_tfgn, inf_tfgn, len_common_tfs
end


# Function to calculate precision, recall, and FPR
function calculate_metrics(gs_tfgn, inf_tfgn, total_possible_interactions)
    precision = zeros(nrow(inf_tfgn))
    recall = zeros(nrow(inf_tfgn))
    fpr = zeros(nrow(inf_tfgn))

    gold_standard_set = Set(zip(gs_tfgn.TF, gs_tfgn.Target))

    # Initialize counters
    true_positives = 0

    for j in 1:nrow(inf_tfgn)
        current_pair = (inf_tfgn[j, :TF], inf_tfgn[j, :Gene])
        
        if current_pair in gold_standard_set
            true_positives += 1
        end
        
        # j is always the total predicted positive at any point j 
        # such that j = TP + FP
        precision[j] = true_positives / j
        recall[j] = true_positives / nrow(gs_tfgn)   # TP + FN = nrow(gs_tfgn)
        
        # obtain FPR
        false_positives = j - true_positives
        true_negatives = total_possible_interactions - nrow(gs_tfgn) - false_positives

        # Breakdown of obtaining true negatives
        # false_negatives_total = nrow(gs_tfgn) - true_positives
        # true_negatives = total_possible_interactions - (true_positives + false_positives + false_negatives_total)

        fpr[j] = false_positives / (false_positives + true_negatives + (true_negatives == 0 ? 1 : 0))
    end

    return precision, recall, fpr
end



# Initialize variables
gs_tfgn = nothing
inf_tfgn = nothing
all_interactions = []
metrics_dict = Dict{String, DataFrame}()

fig = plot(xlims = (0,0.1), ylim= (0,1), dpi = 1000) ; # Initialize the plot
# Loop through inferred GRN paths with an index
for (index, inf_grn_path) in enumerate(inf_grn_paths)
    # Get the base name without extension for the key
    key_name = basename(inf_grn_path)
    key_name = splitext(key_name)[1]  # Remove the extension
    
    println("Processing GRN file: $key_name")
    
    # Read the inferred GRN data
    inferred_grn = CSV.read(inf_grn_path, DataFrame; delim = '\t')
    
    # Filter based on common TFs and targets if provided, else use defaults
    gs_tfgn, inf_tfgn, len_common_tfs = filter_common_genes(gold_standard, inferred_grn; 
                                                              common_targets=common_targets, 
                                                              common_tfs=common_tfs)
    
    total_possible_interactions = len_common_tfs * length(intersect(gs_tfgn.Target, inf_tfgn.Gene))
    precision, recall, fpr = calculate_metrics(gs_tfgn, inf_tfgn, total_possible_interactions)

    # Store the metrics for this GRN
    push!(all_interactions, total_possible_interactions)

    # Create a DataFrame to store metrics
    metrics_df = DataFrame(
        Precision = precision,
        Recall = recall,
        FPR = fpr
    )
    
    # Store the metrics DataFrame in the dictionary
    metrics_dict[key_name] = metrics_df
    
    # Write the metrics to a TSV file
    CSV.write(joinpath(output_dir, "Metrics_GRN_$key_name.tsv"), string.(metrics_df); delim = '\t')
    println("Precision, recall, and FPR saved for GRN $key_name")

    #color_generator = get(ColorSchemes.rainbow,index./length(inf_grn_path))
    # Plot precision-recall curve
    plot!(recall, precision, label=key_name, color=color_generator[index], legend=:topright, lw=2)
    xlabel!("Recall")
    ylabel!("Precision")
    savefig("$(output_dir)/$(fig_save_name)_PRCurve.png")

    # random PR
    # gsRandPR = nrow(gs_tfgn) / total_possible_interactions
    # hline!(y = gsRandPR, linestyle = "dashed", color = "k", label="Random PR")  # Add reference line for random PR
end






































# # Plotting
# fig = plot(xlims = (0,0.1), ylim= (0,1), dpi = 1000) ; # Initialize the plot
# color_generator = [RGB(rand(), rand(), rand()) for _ in 1:length(inf_grn_paths)]  # Generate random colors

# for (index, key_name) in enumerate(keys(metrics_dict))
#     precision = metrics_dict[key_name][!, "Precision"]
#     recall = metrics_dict[key_name][!, "Recall"]
    
#     # Create a new plot for each GRN
#     plot!(recall, precision, label=key_name, linecolor=color_generator[index], legend=:topright, lw =2 );
# end

# # Add a horizontal line for random precision
# # gsRandPR = nrow(gs_tfgn) / total_possible_interactions  # You may need to recalculate this
# # hline!(y = gsRandPR, linestyle = "dashed", color = "k", label="Random")  # Reference line for random PR
# # Add plot details and show/save the plot
# xlabel!("Recall")
# ylabel!("Precision")
# savefig((output_dir * "/Precision_Recall_Curves.png"))


