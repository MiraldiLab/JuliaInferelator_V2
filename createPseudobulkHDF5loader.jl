using  DataFrames, Statistics
using HDF5, SparseArrays, LinearAlgebra, SparseMatricesCSR
using ArgParse
using Dates
using CSV

include("coordConversionHDF5.jl")

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

function createPseudobulk(adata, ident_var::String, save_dir::String)
    # Read the AnnData object
    #adata = h5open(dataset, "r")
    # adata = read(adata)

    row_index = adata["X"]["indptr"]
    column_index = adata["X"]["indices"]
    coordinates = csrToCoordinates(row_index, column_index)

    rindices =  [coord[1] for coord in coordinates]
    cindices =  [coord[2] for coord in coordinates]

    sp_mat = sparse(rindices, cindices, adata["X"]["data"])

    ident_data = adata["obs"][ident_var]
    unique_idents = unique(ident_data["categories"])  # unique index
    ident_codes = ident_data["codes"]
    unique_codes = sort(unique(ident_codes));  # sorted

    features = adata["var"]["_index"];  # contains rownames
    pseudobulk_df = DataFrame(zeros(length(features), length(unique_idents)), unique_idents);

    for (idx, ident) in enumerate(unique_idents)
        curr_indices = findall(x -> x == unique_codes[idx], ident_codes)
        # Subset data for the current identity using logical indexing
        curr_data = sp_mat[:, curr_indices]
        
        # sum counts across all cells per row/feature
        curr_counts = sum(curr_data, dims=2)
        curr_counts = Array(curr_counts)

        pseudobulk_df[!, ident] .= curr_counts
    end

    mat = Matrix(pseudobulk_df);
    row_sums = vec(sum(mat, dims=2));  # Sum along rows;

    pseudobulk_df = insertcols!(pseudobulk_df, 1, :RowNames => features)   # add features as first column
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


function main(args)
    # Parse the command line arguments
    @argparse args begin
        "--dataset", help="Path to the .h5ad file", required = true, arg_type = String
        "--ident_var", help="Identifier variable for pseudobulk creation", required=true,
        arg_type=String
        "--save_dir", help="Directory to save the output"
    end
    
    adata = h5open(args.dataset) do file
        read(file)
    end 

    createPseudobulk(adata, args.ident_var, args.save_dir)
end

# Allows running the script directly or import it without executing main
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main(ARGS)
end



# Testing:

# function main(args)
#     # Load the data from the provided path
#     adata = h5open(args["dataset"], "r") do file
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


