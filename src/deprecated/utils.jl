const InsertionCode = NamedTuple{(:index, :code), Tuple{UInt16, UInt8}}

encode_ins_codes(ins_codes::String) = filter(code -> code.code != 0x20, map(InsertionCode, enumerate(codeunits(ins_codes))))

function decode_ins_codes(ins_codes::Vector{InsertionCode}, len::Integer)
    data = fill(0x20, len)
    for code in ins_codes
        data[code.index] = code.code
    end
    return String(data)
end

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
