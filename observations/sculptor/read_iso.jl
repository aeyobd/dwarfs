import DataFrames: DataFrame

"""
A struct representing a isochrone


"""
struct ISOCMD
    "The name of the file"
    filename::String
    "The MIST and MESA version numbers"
    version::Dict{String, String}
    "The photometric system"
    photo_sys::String
    "The input stellar"
    abun::Dict{String, Float64}
    "The extinction value"
    Av_extinction::Float64
    "The rotation value"
    rot::Float64
    "A list of log ages in years"
    log_ages::Vector{Float64}
    "The number of ages"
    num_ages::Int
    "the headers of the isochrones"
    hdr_list::Vector{String}
    "A vector of isochrones for each age"
    isocmds::Vector{DataFrame}

    function ISOCMD(filename::String; verbose::Bool=true)
        if verbose
            println("Reading in: ", filename)
        end
        version, photo_sys, abun, Av_extinction, rot, ages, num_ages, hdr_list, isocmds = read_isocmd_file(filename)
        new(filename, version, photo_sys, abun, Av_extinction, rot, ages, num_ages, hdr_list, isocmds)
    end
end


function Base.getindex(isocmd::ISOCMD, log_age::Real)::DataFrame
    log_age = convert(Float64, log_age)
    idx = age_index(isocmd, log_age)

    return isocmd.isocmds[idx]
end


"""
    read_isocmd_header(filename::String)

Reads the header attributes of the isochrone file
"""
function read_isocmd_header(filename::String)
    header = [split(line) for line in readlines(filename)][1:10]
    version = Dict("MIST" => header[1][end], "MESA" => header[2][end])
    photo_sys = join(header[3][5:end], " ")
    abun = Dict(header[5][i] => parse(Float64, header[6][i]) for i in 2:5)
    rot = parse(Float64, header[6][end])
    num_ages = parse(Int, header[8][end])
    Av_extinction = parse(Float64, header[9][end])

    return version, photo_sys, abun, Av_extinction, rot, num_ages
end



function read_isocmd_file(filename::String)
    version, photo_sys, abun, Av_extinction, rot, num_ages = read_isocmd_header(filename)

    isocmd_set = Vector{Any}()
    log_age = Float64[]
    counter = 1
    data = [split(line) for line in readlines(filename)[11:end]]
    
    hdr_list = data[counter+2][2:end]

    for i_age in 1:num_ages
        num_eeps = parse(Int, data[counter][end-1])
        num_cols = parse(Int, data[counter][end])
        hdr_list = data[counter+2][2:end]
        
        isocmd = Matrix{Float64}(undef, num_eeps, num_cols)
        for eep in 1:num_eeps
            iso_cmd_chunk = data[counter+2+eep]
            iso_cmd_chunk = parse.(Float64, iso_cmd_chunk)
            isocmd[eep, :] = iso_cmd_chunk
        end

        df = DataFrame(isocmd, Symbol.(hdr_list))

        push!(isocmd_set, df)
        push!(log_age, df.log10_isochrone_age_yr[1])
        counter += 3 + num_eeps + 2
    end
    
    return version, photo_sys, abun, Av_extinction, rot, log_age, num_ages, hdr_list, isocmd_set
end


"""
    age_index(isocmd::ISOCMD, log_age::Float64)

Returns the index of the isochrone closest to the requested age
"""
function age_index(isocmd::ISOCMD, log_age::Float64)
    diff_arr = abs.(isocmd.log_ages .- log_age)
    age_index = argmin(diff_arr)
    
    if log_age > maximum(isocmd.log_ages) || log_age < minimum(isocmd.log_ages)
        error("The requested age $log_age is outside the range. Try log age ", minimum(isocmd.log_ages), " and ", maximum(isocmd.log_ages))
    end
    

    log_age_actual = isocmd.log_ages[age_index]

    if diff_arr[age_index] > 0.001
        @info "The requested age is not in the isochrone set. Rounding $log_age to log_age = $log_age_actual"
    end
    return age_index
end

