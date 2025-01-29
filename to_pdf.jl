using PlutoPDF

# Stolen from https://github.com/JuliaPluto/PlutoSliderServer.jl/blob/main/src/FileHelpers.jl
import Pluto: is_pluto_notebook

flatmap(args...) = vcat(map(args...)...)

function list_files_recursive(dir=".")
    paths = flatmap(walkdir(dir)) do (root, dirs, files)
        joinpath.([root], files)
    end
    relpath.(paths, [dir])
end

"""
Search recursively for Pluto notebook files.

Return paths relative to the search directory.
"""
function find_notebook_files_recursive(start_dir)
    notebook_files = filter(list_files_recursive(start_dir)) do path
        is_pluto_notebook(joinpath(start_dir, path))
    end

    not_interal_notebook_files = filter(notebook_files) do f
        !occursin(".julia", f) || occursin(".julia", start_dir)
    end

    # reverse alphabetical order so that week5 becomes available before week4 :)
    reverse(not_interal_notebook_files)
end
# End of stolen code

notebooks = find_notebook_files_recursive("src")
for notebook in notebooks
    file = joinpath("src", notebook)
    println("Converting $file to PDF...")
    PlutoPDF.pluto_to_pdf(file; open=false)
end
