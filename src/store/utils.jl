# convert e.g. numbers [1,2,3,6,7,8] to ranges [1:3, 6:8]
# useful for compressing residue numbering, which often contain unit ranges
function numbers_to_ranges(numbers::AbstractVector{T}) where T <: Integer
    isempty(numbers) && return UnitRange{T}[]

    ranges = UnitRange{T}[]
    start = numbers[begin]
    stop = start
    for (current, next) in zip(numbers, @view(numbers[begin+1:end]))
        if current+1 == next
            stop = next
        else
            push!(ranges, start:current)
            start = next
        end
    end
    push!(ranges, start:last(numbers))

    return ranges
end
