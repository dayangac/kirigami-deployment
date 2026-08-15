# gui/app.jl -- desktop explorer: viewer (phase 1) + live design backend (phase 2).
#
#   julia -t auto --project=Kirigami/gui Kirigami/gui/app.jl [frames.json] [--theta rad]
#         [--screenshot out.png] [--design PATTERN|"corpus k1a #ID"] [--sigma file|mc|def]
#         [--method range_max|baseline|constrained]
#
# Three-column layout after web/app.html: inputs and the design table (left), canvas
# (centre) with the θ slider and coloured track below it, actions/readouts/exports
# (right). All numbers come from gui/model.jl and gui/backend.jl; this file only wires
# Observables to widgets. Design jobs run on a worker thread (`-t auto`) and are polled
# from a Timer on the render thread, since GLMakie must not be touched from elsewhere.
using GLMakie
using Observables
import JSON
include(joinpath(@__DIR__, "model.jl"))
include(joinpath(@__DIR__, "backend.jl"))
using .GuiModel
using .GuiBackend

const GUI_DIR = @__DIR__
const DATA_DIR = normpath(joinpath(GUI_DIR, "..", "..", "data"))
const DEFAULT_DATA = joinpath(DATA_DIR, "deploy_frames.json")
const PATTERNS_FILE = joinpath(DATA_DIR, "web_patterns.json")
const CORPUS_FILE = joinpath(DATA_DIR, "corpus", "k1a_200.json")
const EXPORT_DIR = joinpath(GUI_DIR, "exports")
const PLAY_SECONDS = 4.0
const TABLE_ROWS = 16           # relabelable rows; the newest designs are shown first
const CORPUS_IDS = 0:9          # corpus graphs offered in the picker
const AVG_EDGE_MM = 40.0        # export scale, as the web app's defaults
const HINGE_MM = 1.5

col(hex) = Makie.RGBf(parse(Makie.Colorant, hex))
rgba(hex, a) = Makie.RGBAf(col(hex), a)

fmt(x; d = 3) = string(round(x; digits = d))
deg(rad) = rad * 180 / pi

function parse_args(args)
    o = Dict{String,Any}("data" => DEFAULT_DATA, "screenshot" => nothing, "theta" => 0.0,
                         "design" => nothing, "sigma" => "file", "method" => "range_max")
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("--screenshot", "--design", "--sigma", "--method")
            o[a[3:end]] = args[i+1]; i += 2
        elseif a == "--theta"
            o["theta"] = parse(Float64, args[i+1]); i += 2
        else
            o["data"] = a; i += 1
        end
    end
    return o
end

"""Build the figure.  Returns (fig, state) so a driver can screenshot or script it."""
function build_app(frames::GuiModel.AbstractBackend, compute::GuiBackend.ComputeBackend)
    ds, _ = GuiModel.designs(frames)
    ds = collect(ds)

    set_theme!(fontsize = 13, backgroundcolor = col(GuiModel.BG),
               textcolor = col(GuiModel.INK))
    fig = Figure(size = (1400, 960), figure_padding = 0)
    left = fig[1, 1] = GridLayout(tellheight = false, valign = :top)
    centre = fig[1, 2] = GridLayout()
    right = fig[1, 3] = GridLayout(tellheight = false, valign = :top)
    colsize!(fig.layout, 1, Fixed(310))
    colsize!(fig.layout, 3, Fixed(300))
    Box(fig[1, 1], color = col(GuiModel.PANEL), strokecolor = col(GuiModel.LINE), strokewidth = 1)
    Box(fig[1, 3], color = col(GuiModel.PANEL), strokecolor = col(GuiModel.LINE), strokewidth = 1)

    # ---- state
    sel = Observable{Union{Nothing,GuiModel.Design}}(isempty(ds) ? nothing : ds[1])
    theta = Observable(0.0)
    playing = Observable(false)
    job = Ref{Union{Nothing,GuiBackend.DesignJob}}(nothing)

    # ---- left: inputs
    Label(left[1, 1:4], "Kirigami Studio", font = :bold, fontsize = 17, halign = :left, padding = (14, 0, 6, 12))
    Label(left[2, 1:4], "1 · pattern", fontsize = 12, color = col(GuiModel.INK2), halign = :left, padding = (14, 0, 2, 2))
    pattern_names = vcat([p.name for p in compute.patterns], ["corpus k1a #$i" for i in CORPUS_IDS])
    pattern_menu = Menu(left[3, 1:4], options = pattern_names, default = isempty(pattern_names) ? nothing : pattern_names[1],
                        fontsize = 11, width = 280, halign = :left, tellwidth = false)
    Label(left[4, 1], "tile", fontsize = 11, color = col(GuiModel.INK2), halign = :left, padding = (14, 0, 0, 0))
    rep_slider = Slider(left[4, 2:3], range = 1:4, startvalue = 2, color_active = col(GuiModel.ACCENT),
                        color_active_dimmed = rgba(GuiModel.ACCENT, 0.5), color_inactive = col(GuiModel.LINE))
    Label(left[4, 4], lift(r -> "$(r)×$(r)", rep_slider.value), fontsize = 11, halign = :right, padding = (0, 14, 0, 0))
    Label(left[5, 1:4], "2 · orientation σ      3 · method", fontsize = 12, color = col(GuiModel.INK2), halign = :left, padding = (14, 0, 8, 2))
    sigma_menu = Menu(left[6, 1:2], options = ["file", "mc", "def"], default = "file", fontsize = 11, width = 120, halign = :left, tellwidth = false)
    method_menu = Menu(left[6, 3:4], options = [("range-max (paper)", "range_max"), ("baseline Eq. (6)", "baseline"),
                                                ("constrained K9", "constrained")], default = "range-max (paper)",
                       fontsize = 11, width = 150, halign = :right, tellwidth = false)
    design_btn = Button(left[7, 1:2], label = "Design", buttoncolor = col(GuiModel.ACCENT), labelcolor = :white,
                        buttoncolor_hover = col(GuiModel.ACCENT), halign = :left, tellwidth = false, padding = (14, 0, 6, 0))
    cancel_btn = Button(left[7, 3:4], label = "Cancel", halign = :right, tellwidth = false, padding = (0, 14, 6, 0))
    progress = Label(left[8, 1:4], "", fontsize = 10, color = col(GuiModel.INK2), halign = :left, justification = :left,
                     padding = (14, 14, 2, 6), tellwidth = false, word_wrap = true)

    # ---- left: design table (fixed rows, relabelled as designs come and go)
    Label(left[9, 1:4], "Designs", fontsize = 12, color = col(GuiModel.INK2), halign = :left, padding = (14, 0, 8, 2))
    for (c, h) in enumerate(("design", "σ", "faces", "Θmax"))
        Label(left[10, c], h, fontsize = 11, color = col(GuiModel.INK2), halign = c == 1 ? :left : :right,
              padding = (c == 1 ? 14 : 4, c == 4 ? 14 : 4, 2, 2))
    end
    rows = []
    for r in 1:TABLE_ROWS
        b = Button(left[10 + r, 1], label = "", halign = :left, fontsize = 10, buttoncolor = :transparent,
                   buttoncolor_hover = col(GuiModel.ACCENT_SOFT), buttoncolor_active = col(GuiModel.ACCENT_SOFT),
                   labelcolor = col(GuiModel.INK), strokewidth = 0, padding = (14, 4, 0, 0))
        ls = [Label(left[10 + r, c], "", fontsize = 10, halign = :right, padding = (4, c == 4 ? 14 : 4, 0, 0)) for c in 2:4]
        push!(rows, (b, ls))
        on(b.clicks) do _
            r <= length(ds) || return
            playing[] = false
            theta[] = 0.0
            sel[] = ds[r]
        end
    end
    function refresh_table!()
        for (r, (b, ls)) in enumerate(rows)
            if r > length(ds)          # unused rows are blank (transparent, no stroke)
                b.label = ""
                foreach(l -> l.text = "", ls)
                continue
            end
            d = ds[r]
            b.label = d.id
            b.labelcolor = d === sel[] ? col(GuiModel.ACCENT) : col(GuiModel.INK)
            ls[1].text = d.sigma; ls[2].text = string(d.F); ls[3].text = fmt(d.theta_max)
        end
    end
    on(_ -> refresh_table!(), sel)
    isempty(ds) && (progress.text = "no designs yet: pick a pattern and press Design")

    # ---- centre: readouts strip, canvas, slider with coloured track
    ro = centre[1, 1] = GridLayout(tellwidth = false)
    badge_box = Box(ro[1, 1], color = col(GuiModel.ACCENT_SOFT), strokewidth = 0, cornerradius = 6,
                    width = 118, height = 22, halign = :left, alignmode = Outside(0, 0, 4, 4))
    badge = Label(ro[1, 1], "closed", fontsize = 11, color = col(GuiModel.ACCENT), halign = :left,
                  padding = (10, 10, 0, 0), tellwidth = false)
    readout(k, r, c) = (Label(ro[1, c], k, fontsize = 11, color = col(GuiModel.INK2), padding = (8, 2, 6, 6));
                        Label(ro[1, c+1], r, fontsize = 11, font = :bold, padding = (0, 4, 6, 6)))
    sT = lift(d -> d === nothing ? "—" : fmt(d.theta_max; d = 5) * " rad", sel)
    sE = lift(d -> d === nothing ? "—" : fmt(d.eps_max; d = 5) * " rad", sel)
    sF = lift(d -> d === nothing ? "—" : "$(d.F) / $(d.nsplit) / $(d.F - d.nsplit)", sel)
    sA = lift((d, t) -> d === nothing ? "—" : fmt(GuiModel.areal_expansion(d, t); d = 2) * "×", sel, theta)
    readout("exact Θmax", sT, 2); readout("certified εmax", sE, 4)
    readout("F / split / hinge", sF, 6); readout("expansion", sA, 8)

    ax = Axis(centre[2, 1], aspect = DataAspect(), backgroundcolor = col(GuiModel.PANEL))
    hidedecorations!(ax); hidespines!(ax)
    polys = Observable(Makie.Polygon[])
    pcols = Observable(Makie.RGBAf[])
    edgec = Observable(col(GuiModel.FACE_EDGE))
    poly!(ax, polys, color = pcols, strokecolor = edgec, strokewidth = 0.8)
    function refresh_canvas(d, t)
        if d === nothing
            polys[] = Makie.Polygon[]; pcols[] = Makie.RGBAf[]
            return
        end
        P = GuiModel.frame(d, t)
        isover = GuiModel.deploy_state(d, t) == GuiModel.over
        polys[] = [Makie.Polygon([Point2f(P[v, 1], P[v, 2]) for v in f]) for f in d.faces]
        pcols[] = [rgba(GuiModel.face_color(s), isover ? 0.55 : 1.0) for s in d.sense]
        edgec[] = isover ? col(GuiModel.BAD) : col(GuiModel.FACE_EDGE)
    end
    onany(refresh_canvas, sel, theta)
    on(sel) do d
        d === nothing && return
        x0, y0, x1, y1 = GuiModel.view_box(d)
        pad = 0.12 * max(x1 - x0, y1 - y0, 1e-9)
        limits!(ax, x0 - pad, x1 + pad, y0 - pad, y1 + pad)
    end
    onany(sel, theta) do d, t
        d === nothing && return
        s = GuiModel.deploy_state(d, t)
        badge.text = GuiModel.badge_text(s)
        if s == GuiModel.over
            badge_box.color = col(GuiModel.BAD_SOFT); badge.color = col(GuiModel.BAD)
        elseif s == GuiModel.usable
            badge_box.color = rgba(GuiModel.GOOD, 0.2); badge.color = col(GuiModel.GOOD)
        else
            badge_box.color = col(GuiModel.ACCENT_SOFT); badge.color = col(GuiModel.ACCENT)
        end
    end
    Label(centre[2, 1], "drag the slider · click a design on the left", fontsize = 10, color = col(GuiModel.INK2),
          halign = :right, valign = :bottom, alignmode = Outside(10), tellwidth = false, tellheight = false)

    bar = centre[3, 1] = GridLayout(tellwidth = false)
    Label(bar[1:2, 1], "Deploy θ", fontsize = 12, padding = (12, 6, 0, 0))
    # coloured track (usable = face-a, certified = accent, over = bad) above the slider
    track = Axis(bar[1, 2], height = 7, limits = (0, 1, 0, 1), backgroundcolor = :transparent)
    hidedecorations!(track); hidespines!(track)
    fu = lift(d -> d === nothing ? 0.0 : GuiModel.track_fractions(d)[1], sel)
    fc = lift(d -> d === nothing ? 0.0 : GuiModel.track_fractions(d)[2], sel)
    poly!(track, lift(u -> Rect2f(u, 0, 1 - u, 1), fu), color = rgba(GuiModel.BAD, 0.5))
    poly!(track, lift(u -> Rect2f(0, 0, u, 1), fu), color = col(GuiModel.FACE_A))
    poly!(track, lift(c -> Rect2f(0, 0, c, 1), fc), color = col(GuiModel.ACCENT))
    sl = Slider(bar[2, 2], range = range(0, pi, length = 721), startvalue = 0.0, color_active = col(GuiModel.ACCENT),
                color_active_dimmed = rgba(GuiModel.ACCENT, 0.5), color_inactive = col(GuiModel.LINE))
    sync = Ref(false)
    on(sl.value) do v
        sync[] && return
        playing[] = false
        theta[] = v
    end
    on(theta) do t            # external θ changes (play) move the slider without echo
        sync[] = true
        set_close_to!(sl, t)
        sync[] = false
    end
    Label(bar[1:2, 3], lift(t -> fmt(deg(t); d = 1) * "°", theta), fontsize = 12, font = :bold, padding = (6, 12, 0, 0))

    # ---- right: deployment, characterize, exports
    Label(right[1, 1], "Deployment", font = :bold, fontsize = 13, halign = :left, padding = (14, 0, 6, 14))
    play_btn = Button(right[2, 1], label = "Play deployment", buttoncolor = col(GuiModel.ACCENT), labelcolor = :white,
                      buttoncolor_hover = col(GuiModel.ACCENT), tellwidth = false, halign = :left, padding = (14, 0, 0, 0))
    on(playing) do p
        play_btn.label = p ? "Pause" : "Play deployment"
    end
    on(play_btn.clicks) do _
        d = sel[]
        d === nothing && return
        if playing[]
            playing[] = false
            return
        end
        playing[] = true
        start, stop = GuiModel.play_range(d, theta[])
        t0 = time()
        @async begin
            while playing[]
                u = (time() - t0) / PLAY_SECONDS
                theta[] = GuiModel.ease_theta(start, stop, u)
                u >= 1 && break
                sleep(1 / 60)
            end
            playing[] = false
        end
    end
    Label(right[3, 1], lift(d -> d === nothing ? "" : "$(d.id)  ·  $(d.kind)\nσ $(d.sigma)  ·  $(count(>(0), d.sense)) of $(d.F) faces turn +θ/2", sel),
          fontsize = 11, color = col(GuiModel.INK2), halign = :left, justification = :left, padding = (14, 0, 4, 4), tellwidth = false)

    Label(right[4, 1], "Characterize", font = :bold, fontsize = 13, halign = :left, padding = (14, 0, 16, 6))
    char_btn = Button(right[5, 1], label = "Characterize selected design", tellwidth = false, halign = :left, padding = (14, 0, 0, 0))
    char_out = Label(right[6, 1], "", fontsize = 10, color = col(GuiModel.INK2), halign = :left, justification = :left,
                     padding = (14, 14, 4, 4), tellwidth = false, word_wrap = true)
    on(char_btn.clicks) do _
        d = sel[]
        d === nothing && return
        try
            ch = GuiBackend.recharacterize(d)
            char_out.text = "Θmax $(fmt(ch.theta_max; d = 6)) rad · εmax $(fmt(ch.eps_max; d = 6)) rad · " *
                            (ch.certified ? "certified" : "not certified") * " at eps $(ch.certificate.eps)\n" *
                            "contacts $(length(ch.contacts)) · binding $(ch.binding) · inverted $(ch.n_inverted)"
        catch e
            char_out.text = sprint(showerror, e)
        end
    end

    Label(right[7, 1], "Export", font = :bold, fontsize = 13, halign = :left, padding = (14, 0, 16, 6))
    xsvg = Button(right[8, 1], label = "laser-cut SVG at current angle", tellwidth = false, halign = :left, padding = (14, 0, 0, 0))
    xsvg0 = Button(right[9, 1], label = "laser-cut SVG (closed)", tellwidth = false, halign = :left, padding = (14, 0, 0, 0))
    x3mf = Button(right[10, 1], label = "3MF solid at current angle", tellwidth = false, halign = :left, padding = (14, 0, 0, 0))
    xjson = Button(right[11, 1], label = "state JSON", tellwidth = false, halign = :left, padding = (14, 0, 0, 0))
    log = Label(right[12, 1], "", fontsize = 10, color = col(GuiModel.INK2), halign = :left, justification = :left,
                padding = (14, 14, 8, 8), tellwidth = false, word_wrap = true)
    function export_file(kind)
        d = sel[]
        d === nothing && return
        mkpath(EXPORT_DIR)
        t = kind == :svg0 ? 0.0 : theta[]
        name = d.id * "_theta" * string(round(Int, deg(t)))
        try
            if kind == :json
                path = joinpath(EXPORT_DIR, name * ".json")
                write(path, GuiModel.state_json(d, t))
            elseif kind == :mf
                path = joinpath(EXPORT_DIR, name * ".3mf")
                w = GuiBackend.export_3mf(d, t, path; avg_edge_mm = AVG_EDGE_MM, hinge_mm = HINGE_MM)
                isempty(w) || (log.text = join(w, "\n"))
            else
                path = joinpath(EXPORT_DIR, name * ".svg")
                if d.source isa GuiBackend.DesignSource
                    write(path, GuiBackend.export_svg_string(d, t; avg_edge_mm = AVG_EDGE_MM, hinge_mm = HINGE_MM))
                else
                    write(path, GuiModel.svg_string(d, t))   # frame-file design: preview only
                end
            end
            log.text = "wrote exports/" * basename(path)
        catch e
            log.text = "export failed: " * sprint(showerror, e)
        end
    end
    on(_ -> export_file(:svg), xsvg.clicks)
    on(_ -> export_file(:svg0), xsvg0.clicks)
    on(_ -> export_file(:mf), x3mf.clicks)
    on(_ -> export_file(:json), xjson.clicks)

    lg = right[13, 1] = GridLayout(halign = :left)
    for (i, (c, txt)) in enumerate(((GuiModel.FACE_A, "+θ/2"), (GuiModel.FACE_B, "−θ/2"),
                                    (GuiModel.ACCENT, "certified"), (GuiModel.BAD, "past Θmax")))
        Box(lg[i, 1], color = col(c), width = 10, height = 10, strokewidth = 0, cornerradius = 2)
        Label(lg[i, 2], txt, fontsize = 11, halign = :left, padding = (4, 0, 1, 1))
    end

    # ---- design jobs: start on the worker, poll from the render thread
    function start_design!(pattern::AbstractString, rep::Int, sigma_rule::AbstractString, method::AbstractString)
        job[] === nothing || return
        m = match(r"^corpus k1a #(\d+)$", pattern)
        f = if m !== nothing
            id = parse(Int, m.captures[1])
            prog -> GuiBackend.design_corpus(compute, id; sigma_rule = sigma_rule, method = method, progress = prog)
        else
            p = compute.patterns[findfirst(x -> x.name == pattern, compute.patterns)]
            prog -> GuiBackend.design_pattern(compute, p; rep_x = rep, rep_y = rep, sigma_rule = sigma_rule,
                                              method = method, progress = prog)
        end
        job[] = GuiBackend.start_job(f)
        design_btn.label = "Designing…"
        progress.text = "starting"
    end
    function finish_job!()
        j = job[]
        j === nothing && return
        prog, done, result = GuiBackend.snapshot(j)
        progress.text = prog
        done || return
        job[] = nothing
        design_btn.label = "Design"
        if result isa GuiModel.Design
            pushfirst!(ds, result)
            push!(compute.designs, result)
            length(ds) > TABLE_ROWS && resize!(ds, TABLE_ROWS)
            progress.text = "done: Θmax $(fmt(result.theta_max; d = 5)) rad, εmax $(fmt(result.eps_max; d = 5)) rad"
            playing[] = false
            theta[] = 0.0
            sel[] = result
        else
            progress.text = "failed: " * sprint(showerror, result)
        end
    end
    on(design_btn.clicks) do _
        pattern_menu.selection[] === nothing && return
        start_design!(pattern_menu.selection[], rep_slider.value[], sigma_menu.selection[], method_menu.selection[])
    end
    on(cancel_btn.clicks) do _
        job[] === nothing || GuiBackend.cancel!(job[])
    end
    poll = Timer(_ -> finish_job!(), 0.0; interval = 0.25)

    notify(sel)               # initial draw
    return fig, (; sel, theta, playing, designs = ds, job, start_design!, finish_job!, poll)
end

if abspath(PROGRAM_FILE) == @__FILE__
    o = parse_args(ARGS)
    compute = GuiBackend.ComputeBackend(PATTERNS_FILE, CORPUS_FILE)
    fig, st = build_app(GuiModel.FileBackend(o["data"]), compute)
    screen = display(fig)
    if o["design"] !== nothing      # scripted design (used for the screenshot)
        st.start_design!(o["design"], 2, o["sigma"], o["method"])
        while st.job[] !== nothing
            sleep(0.25)
        end
    end
    o["theta"] > 0 && (st.theta[] = o["theta"])
    if o["screenshot"] !== nothing
        sleep(0.5)
        save(o["screenshot"], fig)
        println("screenshot written to ", o["screenshot"])
    else
        wait(screen)
    end
end
