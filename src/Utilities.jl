# Generic utility functions

module Utilities

export adjacent_find


function adjacent_find(binary_predicate, array::AbstractArray)
    # Compare adjacent elements of 'array' using 'binary_predicate'.
    # Returns index of first 'true' from 'binary_predicate'.
    # Returns last index if no match found
    # Following C++ std::adjacent_find(), https://en.cppreference.com/w/cpp/algorithm/adjacent_find
    # No find-like function in Julia that takes a binary predicate?

    if length(array) == 0
        throw(DomainError("invalid argument: empty array"))
    end

    idxs = LinearIndices(array)
    for i in idxs[1:end-1]
        @inbounds if binary_predicate(array[idxs[i]], array[idxs[i+1]])
            return idxs[i]
        end
    end
    return idxs[end] # No match found
end


end # module
