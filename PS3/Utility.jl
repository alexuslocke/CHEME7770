function contains(string,token)
    return occursin(token,string)
end

function generate_atom_matrix(path_to_atom_file::AbstractString,data_dictionary::Dict{AbstractString,Any})


  # how many metabolite symbols do we have in *the model*?
  list_of_metabolite_symbols_model = data_dictionary["list_of_metabolite_symbols"]
  number_of_metabolites = length(list_of_metabolite_symbols_model)

  # initialize -
  tmp_array::Array{AbstractString} = AbstractString[]
  atom_array = zeros(number_of_metabolites,6)


  # load the atom file -
  try

    open(path_to_atom_file,"r") do model_file
      for line in eachline(model_file)

          if (contains(line,"//") == false && line != "")
            push!(tmp_array,chomp(line))
          end
      end
    end

    # ok, create a local dictionary w/the atom records -
    local_dictionary::Dict{AbstractString,Any} = Dict{AbstractString,Any}()
    for record in tmp_array

      # split -
      split_array = split(record,",")

      # get my key -
      key = split_array[1]  # Metabolite symbol -

      # local array -
      local_atom_array = zeros(6)
      local_atom_array[1] = parse(Float64,split_array[2]) # C
      local_atom_array[2] = parse(Float64,split_array[3]) # H
      local_atom_array[3] = parse(Float64,split_array[4]) # N
      local_atom_array[4] = parse(Float64,split_array[5]) # O
      local_atom_array[5] = parse(Float64,split_array[6]) # P
      local_atom_array[6] = parse(Float64,split_array[7]) # S

      # store -
      local_dictionary[key] = local_atom_array
    end

    # ok, so now we have the local dictionary, we can lookup (in order) the metabolites in the model -
    for (index,model_metabolite_symbol) in enumerate(list_of_metabolite_symbols_model)

      # what is the atom array for *this metabolite*?
      local_atom_array = local_dictionary[model_metabolite_symbol]

      for (atom_index,coefficient) in enumerate(local_atom_array)
        atom_array[index,atom_index] = local_atom_array[atom_index]
      end
    end

  catch err
    showerror(stdout, err, backtrace());println()
  end

  return atom_array
end
