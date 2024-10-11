using Muon
using  DataFrames, Statistics
using SparseArrays, LinearAlgebra, SparseMatricesCSR
using ArgParse
using Dates
using CSV


"""
    createPseudobulk(adata, ident_var, save_dir)

Create a pseudobulk DataFrame from the given AnnData object.

# Arguments
- `adata::AnnData`: The AnnData object containing the data.
- `ident_var::String`: The identifier variable for grouping cells.
- `save_dir::String`: The directory to save the output pseudobulk DataFrame.

# Returns
- A DataFrame containing the pseudobulk data, where each column corresponds to a unique identifier and the rows represent features (e.g. genes, peaks).
"""

# reference: https://docs.juliahub.com/Muon/QfqCh/0.1.1/objects/

function createPseudobulk(adata, ident_var::String, save_dir::String)
    # Read the AnnData object
    #adata = h5open(dataset, "r")
    # adata = read(adata)

    sp_mat = adata.X  #sparse matrix (cells by feature)
    
    ident_col = adata.obs[!, ident_var]  # element of column corresponding to ident_var. Note: adata.obs -> metadata::dataframe
    unique_idents = unique(ident_col)

    features = adata.var[!, "name"]  # contains rownames
    pseudobulk_df = DataFrame(zeros(length(features), length(unique_idents)), unique_idents)

    for ident in unique_idents
        curr_indices = findall(x -> x == ident, ident_col)
        # Subset data for the current identity using logical indexing
        curr_data = sp_mat[curr_indices, :]
        
        # sum counts across all cells per row/feature
        curr_counts = sum(curr_data, dims=1)
        pseudobulk_df[!, ident] .= vec(curr_counts)
    end

    mat = Matrix(pseudobulk_df);
    row_sums = vec(sum(mat, dims=2));  # Sum along rows;

    pseudobulk_df = insertcols!(pseudobulk_df, 1, :RowNames => features) 
    pseudobulk_df = pseudobulk_df[row_sums .> 0, :]

    # Save the pseudobulk DataFrame
    current_time = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    output_path = joinpath(save_dir, "PseudobulkData_$(ident_var)_$(current_time).tsv")
    # open(output_path, "w") do io
    #     write(io, pseudobulk_df)
    # end
    CSV.write(output_path, string.(pseudobulk_df); delim = '\t')

    println("Pseudobulk DataFrame saved to: ", output_path)
    return pseudobulk_df 
end


# function main(args)
#     # Parse the command line arguments
#     @argparse args begin
#         "--dataset", help="Path to the .h5ad file", required = true, arg_type = String
#         "--ident_var", help="Identifier variable for pseudobulk creation", required=true,
#         arg_type=String
#         "--save_dir", help="Directory to save the output"
#     end
    
#     adata = readh5ad(args.dataset) do file
#         read(file)
#     end 

#     createPseudobulk(adata, args.ident_var, args.save_dir)
# end

# # Allows running the script directly or import it without executing main
# if abspath(PROGRAM_FILE) == abspath(@__FILE__)
#     main(ARGS)
# end



# Testing:

# function main(args)
#     # Load the data from the provided path
#     adata = readh5ad(args["dataset"], "r") do file
#         read(file)
#     end

#     # Call the createPseudobulk function
#     createPseudobulk(adata, args["ident_var"], args["save_dir"])
# end

# args = Dict(
#     "dataset" => "/data/miraldiNB/Katko/Projects/Asthma/Asthma_05_17_2024.h5ad",
#     "ident_var" => "celltype2_offspring_rep_comb",
#     "save_dir" => "/data/miraldiNB/Katko/Projects/Asthma/Analysis_New/Pseudobulk_Analysis/Asthma_05_17"
# )

# main(args)
