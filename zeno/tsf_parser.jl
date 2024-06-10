@Base.kwdef mutable struct ParserState
    header = true
    variable = nothing
    data::Dict{String, Vector{Float64}} = Dict()
    sizes::Dict{String, Vector{Int}} = Dict()
    dtypes::Dict{String, String} = Dict()
    finished = false
    set = "Particles"
end


function parse_tsf(filename::String; kwargs...)
    state = ParserState(;kwargs...)

    open(filename, "r") do file
        for line in eachline(file)
            parse_line(strip(line), state)
        end
    end

    check_consistency(state)
    return state
end


function parse_line(line, state::ParserState)
    if isempty(line)
        return
    end

    if state.header
        if startswith(line, "set " * state.set)
            state.header = false
        end
        return
    end

    if startswith(line, "tes")
        state.variable = nothing
        state.finished = true
        return
    end


    elements = split(line)
    if startswith(line, "float") || startswith(line, "int")
        state.variable = match(r"(\w+)", elements[2]).captures[1]
        var = state.variable
        state.dtypes[var] = elements[1]
        matches = eachmatch(r"\d+", elements[2])
        sizes = [m.match for m in matches]
        state.sizes[var] = parse.(Int, sizes)
        state.data[var] = Float64[]

        append!(state.data[var], parse.(Float64, elements[3:end]))
        return 
    end

    if state.variable != nothing
        var = state.variable
        append!(state.data[var], parse.(Float64, elements))
    end

end


function check_consistency(state::ParserState)
    if !state.finished
        error("incomplete input")
    end
    for (var, data) in state.data
        expected = prod(state.sizes[var])
        N = length(data)
        if N != expected
            error("incomplete input for $var. expected shape $(expected), got $N")
        end
    end

end

function inflate(parser::ParserState, var::String)
    return reshape(parser.data[var], reverse(parser.sizes[var])...)
end
