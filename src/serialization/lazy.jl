mutable struct Deserializer <: AbstractDict{String,ProteinStructure}
    filename::String
    file::HDF5.File
    keys::Set{String}

    function Deserializer(filename::AbstractString)
        file = HDF5.h5open(filename, "r")
        deserializer = new(filename, file, Set(keys(file)))
        finalizer(close, deserializer)
        return deserializer
    end
end

Base.close(deserializer::Deserializer) = close(deserializer.file)

Base.keys(deserializer::Deserializer) = deserializer.keys
Base.length(deserializer::Deserializer) = length(keys(deserializer))
Base.getindex(deserializer::Deserializer, key::AbstractString) = read(deserializer.file[key], ProteinStructure)

function Base.iterate(deserializer::Deserializer, state...)
    handler(::Nothing) = nothing
    handler((key,state)::Tuple{String,Int}) = (key => deserializer[key], state)
    return handler(iterate(deserializer.keys, state...))
end

Base.show(io::IO, deserializer::Deserializer) = print(io, "$(Deserializer)(\"$(deserializer.filename)\")")
Base.show(io::IO, ::MIME"text/plain", deserializer::Deserializer) = print(io, summary(deserializer))
