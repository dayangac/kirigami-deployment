# diff_results.jl -- compare a C++ results CSV with its Julia (`_julia`) counterpart,
# column by column, and write a Markdown report.
#
#   julia --project=Kirigami/scripts Kirigami/scripts/diff_results.jl CPP.csv JULIA.csv
#         [--rtol 1e-9] [--atol 1e-12] [--config diff_tolerances.toml] [--md OUT.md]
#   julia --project=Kirigami/scripts Kirigami/scripts/diff_results.jl --all
#         [--results DIR] [--md results/JULIA_VS_CPP.md] [--config ...]
#
# Rows are matched by key columns (the first columns of the configured key list that the
# header contains; row order otherwise). Per column the report gives the comparison class
# (exact for strings/integers, numeric within rtol/atol, informational for columns known
# to differ for a documented non-bug reason, skipped for timings), the number of rows
# compared and mismatching, the worst absolute / relative difference and where it occurs,
# plus the rows present on one side only. `--all` walks every `<dir>_julia` under the
# results tree, pairs each of its CSVs with the same-named C++ file (a Julia shard file
# `x_3.csv` falls back to the merged C++ `x.csv`, compared on the shard's keys) and writes
# one section per file to results/JULIA_VS_CPP.md.
using TOML
using Printf

const DEFAULT_CONFIG = joinpath(@__DIR__, "diff_tolerances.toml")

# ---------------------------------------------------------------- csv

struct Table
    path::String
    header::Vector{String}
    rows::Vector{Vector{String}}
end

function read_csv(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && return Table(path, String[], Vector{String}[])
    header = String.(split(lines[1], ','))
    rows = Vector{String}[]
    for l in lines[2:end]
        isempty(strip(l)) && continue
        f = String.(split(l, ','))
        # pad / truncate ragged rows (the C++ K7 main CSV ends its rows with a trailing comma)
        length(f) < length(header) && append!(f, fill("", length(header) - length(f)))
        length(f) > length(header) && (f = f[1:length(header)])
        push!(rows, f)
    end
    return Table(String(path), header, rows)
end

# ---------------------------------------------------------------- config

struct ColRule
    rtol::Float64
    atol::Float64
    info::Bool
    skip::Bool
    note::String
end

# glob -> regex, `*` matches anything
glob_re(p::AbstractString) = Regex("^" * join(map(s -> replace(s, r"([.^\$+?()\[\]{}|\\])" => s"\\\1"), split(p, '*')), ".*") * "\$")
globmatch(p::AbstractString, s::AbstractString) = occursin(glob_re(p), s)

struct Config
    rtol::Float64
    atol::Float64
    keys::Vector{String}
    skip::Vector{String}
    files::Dict{String,Any}   # section name (glob) -> table
end

function load_config(path::AbstractString)
    t = isfile(path) ? TOML.parsefile(path) : Dict{String,Any}()
    g = get(t, "global", Dict{String,Any}())
    files = Dict{String,Any}(k => v for (k, v) in t if k != "global")
    return Config(Float64(get(g, "rtol", 1e-9)), Float64(get(g, "atol", 1e-12)),
                  String.(get(g, "keys", String[])), String.(get(g, "skip", ["secs"])), files)
end

# the sections whose glob matches this file, most specific (longest) last. `base` is
# `<parent dir>/<basename>` (e.g. `k9b_julia/shard_0.csv`); a section name without a `/`
# matches the basename alone, one with a `/` matches the directory-qualified name, so the
# shard files of different apps (all called shard_N.csv) can carry different rules.
function file_sections(cfg::Config, base::AbstractString)
    bn = basename(base)
    secs = [(k, v) for (k, v) in cfg.files if globmatch(k, occursin('/', k) ? base : bn)]
    sort!(secs; by = kv -> length(kv[1]))
    return secs
end

function column_rule(cfg::Config, base::AbstractString, col::AbstractString)
    rtol = cfg.rtol; atol = cfg.atol; info = false; note = ""
    skip = any(p -> globmatch(p, col), cfg.skip)
    for (_, sec) in file_sections(cfg, base)
        # exact column names override globs, globs override the section wildcard
        cands = [(k, v) for (k, v) in sec if k != "rows" && v isa AbstractDict && globmatch(k, col)]
        sort!(cands; by = kv -> (kv[1] == "*" ? 0 : occursin('*', kv[1]) ? 1 : 2))
        for (_, r) in cands
            haskey(r, "rtol") && (rtol = Float64(r["rtol"]))
            haskey(r, "atol") && (atol = Float64(r["atol"]))
            haskey(r, "info") && (info = Bool(r["info"]))
            haskey(r, "note") && (note = String(r["note"]))
        end
    end
    return ColRule(rtol, atol, info, skip, note)
end

# rows marked informational by `column == value` rules
function row_rules(cfg::Config, base::AbstractString)
    out = Tuple{String,String,String}[]
    for (_, sec) in file_sections(cfg, base)
        for r in get(sec, "rows", Any[])
            push!(out, (String(r["column"]), String(r["value"]), String(get(r, "note", ""))))
        end
    end
    return out
end

# ---------------------------------------------------------------- comparison

parse_num(s::AbstractString) = tryparse(Float64, strip(s))
isintstr(s::AbstractString) = occursin(r"^\s*[+-]?\d+\s*$", s)
isnumstr(s::AbstractString) = parse_num(s) !== nothing || occursin(r"^\s*-?(nan|inf)\s*$"i, s)

mutable struct ColReport
    name::String
    class::String          # exact | numeric | info | skipped | missing
    n::Int                 # rows compared
    n_bad::Int             # mismatches (beyond tolerance)
    n_info_bad::Int        # mismatches on informational rows
    max_abs::Float64
    max_rel::Float64
    worst_key::String
    rule::ColRule
    example::String
end

# a file section may name its own `keys = [...]` (e.g. k2c_drift.csv, whose rows are sorted
# by a column with ties); otherwise the global list applies
function key_columns(cfg::Config, header::Vector{String}, base::AbstractString = "")
    keys = cfg.keys
    for (_, sec) in file_sections(cfg, base)
        haskey(sec, "keys") && (keys = String.(sec["keys"]))
    end
    ks = [k for k in keys if k in header]
    return ks
end

row_key(row, idx) = join((row[i] for i in idx), "/")

function compare_tables(cpp::Table, jl::Table, cfg::Config, base::AbstractString; subset::Bool = false,
                        root::AbstractString = "")
    cols = [c for c in jl.header if c in cpp.header]
    only_cpp = [c for c in cpp.header if !(c in jl.header)]
    only_jl = [c for c in jl.header if !(c in cpp.header)]
    ks = key_columns(cfg, cols, base)
    ci = Dict(c => i for (i, c) in enumerate(cpp.header))
    ji = Dict(c => i for (i, c) in enumerate(jl.header))
    kci = [ci[k] for k in ks]; kji = [ji[k] for k in ks]
    # pair rows by key; a repeated key (e.g. the reference cases, all id = -1) is
    # disambiguated by its occurrence number, so a subset run still pairs in order
    function keyed(rows, idx)
        seen_k = Dict{String,Int}()
        out = String[]
        for row in rows
            k = row_key(row, idx)
            n = get(seen_k, k, 0) + 1
            seen_k[k] = n
            push!(out, n == 1 ? k : k * "#" * string(n))
        end
        return out
    end
    ckeys = isempty(ks) ? [string(r) for r in 1:length(cpp.rows)] : keyed(cpp.rows, kci)
    jkeys = isempty(ks) ? [string(r) for r in 1:length(jl.rows)] : keyed(jl.rows, kji)
    cmap = Dict(k => r for (r, k) in enumerate(ckeys))
    pairs = Tuple{Int,Int,String}[]
    missing_cpp = String[]
    seen = Set{Int}()
    for (r, k) in enumerate(jkeys)
        if haskey(cmap, k)
            push!(pairs, (cmap[k], r, k))
            push!(seen, cmap[k])
        else
            push!(missing_cpp, k)
        end
    end
    dup = false
    missing_jl = [ckeys[r] for r in 1:length(cpp.rows) if !(r in seen)]
    rrules = row_rules(cfg, base)
    info_row(row) = any(rr -> haskey(ji, rr[1]) && row[ji[rr[1]]] == rr[2], rrules)

    reports = ColReport[]
    for c in cols
        rule = column_rule(cfg, base, c)
        rep = ColReport(c, "exact", 0, 0, 0, 0.0, 0.0, "", rule, "")
        if rule.skip
            rep.class = "skipped"
            push!(reports, rep)
            continue
        end
        c in ks && (rep.class = "key")
        # class: numeric iff some compared value is a non-integer number on either side
        numeric = false
        for (a, b, _) in pairs
            x = cpp.rows[a][ci[c]]; y = jl.rows[b][ji[c]]
            if (isnumstr(x) && !isintstr(x)) || (isnumstr(y) && !isintstr(y))
                numeric = true
                break
            end
        end
        rep.class == "key" || (rep.class = numeric ? "numeric" : "exact")
        rule.info && (rep.class = "info")
        for (a, b, k) in pairs
            x = cpp.rows[a][ci[c]]; y = jl.rows[b][ji[c]]
            rep.n += 1
            rowinfo = info_row(jl.rows[b])
            bad = false
            if numeric
                fx = parse_num(x); fy = parse_num(y)
                if fx === nothing || fy === nothing
                    bad = strip(x) != strip(y)
                    bad && rep.max_abs < Inf && (rep.max_abs = Inf)
                elseif isnan(fx) || isnan(fy)
                    bad = !(isnan(fx) && isnan(fy))
                elseif isinf(fx) || isinf(fy)
                    bad = fx != fy
                else
                    d = abs(fx - fy)
                    rel = d / max(abs(fx), abs(fy), 1e-300)
                    if d > rep.max_abs
                        rep.max_abs = d
                        rep.max_rel = rel
                        rep.worst_key = k
                        rep.example = "$(x) vs $(y)"
                    end
                    bad = d > rule.atol + rule.rtol * max(abs(fx), abs(fy))
                end
            else
                bad = strip(x) != strip(y)
                if bad && isempty(rep.example)
                    rep.worst_key = k
                    rep.example = "$(x) vs $(y)"
                end
            end
            if bad
                (rule.info || rowinfo) ? (rep.n_info_bad += 1) : (rep.n_bad += 1)
            end
        end
        push!(reports, rep)
    end
    return (cols = reports, only_cpp = only_cpp, only_jl = only_jl, keys = ks,
            n_cpp = length(cpp.rows), n_jl = length(jl.rows), n_pairs = length(pairs),
            missing_cpp = missing_cpp, missing_jl = missing_jl, subset = subset,
            row_rules = rrules, positional = isempty(ks),
            cpp_rel = isempty(root) ? cpp.path : relpath(cpp.path, root),
            jl_rel = isempty(root) ? jl.path : relpath(jl.path, root))
end

# ---------------------------------------------------------------- markdown

fmtg(v::Real) = isfinite(v) ? @sprintf("%.3g", v) : string(v)

function markdown_section(io::IO, title::AbstractString, cpp::Table, jl::Table, res)
    bad_cols = [c for c in res.cols if c.n_bad > 0]
    info_cols = [c for c in res.cols if c.n_info_bad > 0]
    verdict = isempty(bad_cols) ? "PASS" : "DIFF"
    println(io, "## ", title, " -- ", verdict)
    println(io)
    println(io, "- C++: `", res.cpp_rel, "` (", res.n_cpp, " rows); Julia: `", res.jl_rel, "` (", res.n_jl, " rows)")
    println(io, "- rows matched: ", res.n_pairs, " (key: ",
            isempty(res.keys) ? "positional" : join(res.keys, ", "), ")",
            res.subset ? "; Julia shard compared against the merged C++ file (missing-in-Julia not reported)" : "")
    if !isempty(res.missing_cpp)
        println(io, "- rows only in Julia: ", length(res.missing_cpp), " (", join(first(res.missing_cpp, 8), "; "),
                length(res.missing_cpp) > 8 ? "; ..." : "", ")")
    end
    if !isempty(res.missing_jl) && !res.subset
        println(io, "- rows only in C++: ", length(res.missing_jl), " (", join(first(res.missing_jl, 8), "; "),
                length(res.missing_jl) > 8 ? "; ..." : "", ")")
    end
    isempty(res.only_cpp) || println(io, "- columns only in C++: ", join(res.only_cpp, ", "))
    isempty(res.only_jl) || println(io, "- columns only in Julia: ", join(res.only_jl, ", "))
    for (col, val, note) in res.row_rules
        println(io, "- rows with `", col, " == ", val, "` are informational: ", note)
    end
    if isempty(bad_cols)
        println(io, "- verdict: every non-informational column agrees",
                isempty(info_cols) ? "" : "; informational differences in " * join((c.name for c in info_cols), ", "), ".")
    else
        println(io, "- verdict: ", length(bad_cols), " column(s) differ beyond tolerance: ",
                join((c.name for c in bad_cols), ", "), ".")
    end
    println(io)
    println(io, "| column | class | rows | mismatch | info-mismatch | max abs | max rel | worst row | example (C++ vs Julia) | tolerance / note |")
    println(io, "|---|---|---:|---:|---:|---:|---:|---|---|---|")
    for c in res.cols
        c.class == "skipped" && continue
        tol = c.class in ("numeric", "info") ? "rtol " * fmtg(c.rule.rtol) * (c.rule.atol > 0 ? ", atol " * fmtg(c.rule.atol) : "") : ""
        note = isempty(c.rule.note) ? tol : (isempty(tol) ? c.rule.note : tol * "; " * c.rule.note)
        flag = c.n_bad > 0 ? " **" * string(c.n_bad) * "**" : string(c.n_bad)
        println(io, "| ", c.name, " | ", c.class, " | ", c.n, " | ", flag, " | ", c.n_info_bad, " | ",
                c.class in ("numeric", "info") ? fmtg(c.max_abs) : "", " | ",
                c.class in ("numeric", "info") ? fmtg(c.max_rel) : "", " | ",
                (c.n_bad + c.n_info_bad > 0) ? c.worst_key : "", " | ",
                (c.n_bad + c.n_info_bad > 0) ? c.example : "", " | ", note, " |")
    end
    skipped = [c.name for c in res.cols if c.class == "skipped"]
    isempty(skipped) || println(io, "\nskipped (timing): ", join(skipped, ", "))
    println(io)
    return verdict
end

# ---------------------------------------------------------------- pairing under results/

# every (cpp_csv, julia_csv, subset) pair under `root`
function find_pairs(root::AbstractString, cfg::Config)
    out = Tuple{String,String,Bool}[]
    for (dir, subdirs, files) in walkdir(root)
        endswith(dir, "_julia") || continue
        cppdir = dir[1:end - length("_julia")]
        isdir(cppdir) || continue
        for f in sort(files)
            endswith(f, ".csv") || continue
            jl = joinpath(dir, f)
            cpp = joinpath(cppdir, f)
            # a section may redirect to another C++ file (`cpp = "k6_final.csv"`: the
            # archived k6.csv is a 25-row partial, k6_final.csv the merged run)
            for (_, sec) in file_sections(cfg, joinpath(basename(dir), f))
                haskey(sec, "cpp") && (cpp = joinpath(cppdir, String(sec["cpp"])))
            end
            if isfile(cpp)
                push!(out, (cpp, jl, false))
                continue
            end
            # shard file x_N.csv -> merged x.csv
            m = match(r"^(.*)_\d+\.csv$", f)
            if m !== nothing && isfile(joinpath(cppdir, m.captures[1] * ".csv"))
                push!(out, (joinpath(cppdir, m.captures[1] * ".csv"), jl, true))
            else
                push!(out, ("", jl, false))
            end
        end
    end
    return out
end

function diff_pair(io::IO, cpp_path::AbstractString, jl_path::AbstractString, cfg::Config, subset::Bool, root::AbstractString)
    title = relpath(jl_path, root)
    if isempty(cpp_path)
        println(io, "## ", title, " -- NO C++ COUNTERPART\n")
        return "NONE"
    end
    cpp = read_csv(cpp_path); jl = read_csv(jl_path)
    res = compare_tables(cpp, jl, cfg, joinpath(basename(dirname(jl_path)), basename(jl_path));
                         subset = subset, root = root)
    return markdown_section(io, title, cpp, jl, res)
end

# ---------------------------------------------------------------- main

function main(args::Vector{String})
    rtol = nothing; atol = nothing
    config = DEFAULT_CONFIG
    md = ""
    all = false
    results = normpath(joinpath(@__DIR__, "..", "..", "results"))
    files = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--rtol" && i < length(args); rtol = parse(Float64, args[i+1]); i += 2
        elseif a == "--atol" && i < length(args); atol = parse(Float64, args[i+1]); i += 2
        elseif a == "--config" && i < length(args); config = args[i+1]; i += 2
        elseif a == "--md" && i < length(args); md = args[i+1]; i += 2
        elseif a == "--results" && i < length(args); results = args[i+1]; i += 2
        elseif a == "--all"; all = true; i += 1
        else; push!(files, a); i += 1
        end
    end
    cfg = load_config(config)
    (rtol !== nothing || atol !== nothing) &&
        (cfg = Config(something(rtol, cfg.rtol), something(atol, cfg.atol), cfg.keys, cfg.skip, cfg.files))

    pairs = if all
        find_pairs(results, cfg)
    else
        length(files) == 2 || error("usage: diff_results.jl CPP.csv JULIA.csv [--rtol r] [--md out.md] | --all")
        [(files[1], files[2], false)]
    end
    all && isempty(md) && (md = joinpath(results, "JULIA_VS_CPP.md"))

    buf = IOBuffer()
    println(buf, "# Julia vs C++ results\n")
    println(buf, "Generated by `Kirigami/scripts/diff_results.jl` on ", Libc.strftime("%Y-%m-%d %H:%M", time()),
            "; tolerances from `", relpath(config, results), "` (default rtol ", fmtg(cfg.rtol), ", atol ", fmtg(cfg.atol), ").")
    println(buf, "Classes: exact (strings / integers), numeric (within tolerance), info (known, documented non-bug",
            " difference -- reported, not counted), key (row identifier).\n")
    verdicts = Tuple{String,String}[]
    for (cpp, jl, subset) in pairs
        v = diff_pair(buf, cpp, jl, cfg, subset, results)
        push!(verdicts, (relpath(jl, results), v))
    end
    summary = IOBuffer()
    println(summary, "## Summary\n")
    println(summary, "| Julia file | verdict |\n|---|---|")
    for (f, v) in verdicts
        println(summary, "| ", f, " | ", v, " |")
    end
    println(summary)
    text = String(take!(buf))
    # summary first, sections after
    text = replace(text, "Classes: exact" => String(take!(summary)) * "Classes: exact"; count = 1)
    if isempty(md)
        print(text)
    else
        mkpath(dirname(md))
        write(md, text)
        println("wrote ", md)
        for (f, v) in verdicts
            println("  ", rpad(f, 48), v)
        end
    end
    return any(v -> v[2] == "DIFF", verdicts) ? 1 : 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(copy(ARGS)))
end
