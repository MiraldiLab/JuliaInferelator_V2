using SparseArrays

function csrToCoordinates(cumsum_nonzero, col_indices)
    coordinates = []  # To store (row, col) pairs
    
    # Number of rows
    n_rows = length(cumsum_nonzero) - 1

    # Adjust for 0-based column indices
    if minimum(col_indices) == 0
        col_indices = col_indices .+ 1  # Julia is 1-based, adjust index
    end
    
    # Iterate through each row
    for row in 1:n_rows
        # Get the start and end index from the cumsum_nonzero array
        start_idx = cumsum_nonzero[row] + 1   # Julia is 1-based, adjust index
        end_idx = cumsum_nonzero[row + 1]
        
        # Extract the column indices for this row
        for col_idx in start_idx:end_idx
            # Julia is 1-based, adjust index
            push!(coordinates, (row, col_indices[col_idx]))  # Store (row, col)
        end
    end
    
    return coordinates
end



