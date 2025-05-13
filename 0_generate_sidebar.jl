

function extract_title(plutofile)
    data = read(open(plutofile), String)
    idx  = 1
    while true
        mbeg = match(r"md\"\"\"", data, idx)
        isnothing(mbeg) && error("""Could not find 'md\"\"\"' to begin title search""")
        mend = match(r"md\"\"\"", data, mbeg.offset + 5)
        isnothing(mend) && error("""Could not find terminating '\"\"\"' to delimit search""")
        idx = mend.offset + 3

        region = strip(data[mbeg.offset+5:mend.offset])
        for line in split(region, "\n")
            line = strip(line)
            if startswith(line, "# ")
                return strip(line[3:end])
            end
        end
    end
end

function find_chapters(dir=".")
    res = filter(n -> !isnothing(match(r"^[0-9][0-9]_[a-zA-Z_]+\.jl$", n)), readdir(dir))
    sort(res)
end

function dump_sidebar(file="sidebar.md"; course_title, prefix)
    io = IOBuffer()
    println(io, "**[$(course_title)]($(prefix)/)**", "\n")

    for chapter in find_chapters()
        base, _ = splitext(basename(chapter))
        title = extract_title(chapter)
        println(io, "1. [$(title)]($(prefix)/$(base).html)")
    end
    open(file, "w") do fp
        print(fp, String(take!(io)))
    end
end

function main()
    dump_sidebar(; course_title="Numerical analysis",
                   prefix="https://teaching.matmat.org/numerical-analysis")
end

main()
