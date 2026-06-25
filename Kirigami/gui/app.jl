# gui/app.jl -- desktop explorer, phase 1 (viewer).
#
#   julia --project=Kirigami/gui Kirigami/gui/app.jl [patterns.json] [--theta rad] [--screenshot out.png]
#
# Three-column layout after web/app.html: pattern table (left), canvas (centre) with the
# θ slider and coloured track below it, actions/readouts (right).  All numbers come from
# gui/model.jl; this file only wires Observables to widgets.
using GLMakie
using Observables
import JSON
include(joinpath(@__DIR__, "model.jl"))
using .GuiModel

const GUI_DIR = @__DIR__
const DEFAULT_DATA = normpath(joinpath(GUI_DIR, "..", "..", "data", "deploy_frames.json"))
const EXPORT_DIR = joinpath(GUI_DIR, "exports")
const PLAY_SECONDS = 4.0

col(hex) = Makie.RGBf(parse(Makie.Colorant, hex))
rgba(hex, a) = Makie.RGBAf(col(hex), a)

fmt(x; d = 3) = string(round(x; digits = d))
deg(rad) = rad * 180 / pi

function parse_args(args)
    data = DEFAULT_DATA
    shot = nothing
    theta0 = 0.0
    i = 1
    while i <= length(args)
        if args[i] == "--screenshot"
            shot = args[i+1]; i += 2
        elseif args[i] == "--theta"
            theta0 = parse(Float64, args[i+1]); i += 2
        else
            data = args[i]; i += 1
        end
    end
    return data, shot, theta0
end

"""Build the figure.  Returns (fig, state) so a driver can screenshot or script it."""
function build_app(backend::GuiModel.AbstractBackend)
    ds, skipped = GuiModel.designs(backend)

    set_theme!(fontsize = 13, backgroundcolor = col(GuiModel.BG),
               textcolor = col(GuiModel.INK))
    fig = Figure(size = (1380, 900), figure_padding = 0)
    left = fig[1, 1] = GridLayout(tellheight = false, valign = :top)
    centre = fig[1, 2] = GridLayout()
    right = fig[1, 3] = GridLayout(tellheight = false, valign = :top)
    colsize!(fig.layout, 1, Fixed(300))
    colsize!(fig.layout, 3, Fixed(300))
    Box(fig[1, 1], color = col(GuiModel.PANEL), strokecolor = col(GuiModel.LINE), strokewidth = 1)
    Box(fig[1, 3], color = col(GuiModel.PANEL), strokecolor = col(GuiModel.LINE), strokewidth = 1)

    # ---- state
    sel = Observable{Union{Nothing,GuiModel.Design}}(isempty(ds) ? nothing : ds[1])
    theta = Observable(0.0)
    playing = Observable(false)

    # ---- left: pattern table
    Label(left[1, 1:4], "Kirigami Studio", font = :bold, fontsize = 17, halign = :left, padding = (14, 0, 6, 14))
    Label(left[2, 1:4], "Patterns", fontsize = 12, color = col(GuiModel.INK2), halign = :left, padding = (14, 0, 2, 2))
    hdr = ("design", "σ", "faces", "Θmax")
    for (c, h) in enumerate(hdr)
        Label(left[3, c], h, fontsize = 11, color = col(GuiModel.INK2), halign = c == 1 ? :left : :right,
              padding = (c == 1 ? 14 : 4, c == 4 ? 14 : 4, 2, 2))
    end
    rowbuttons = Button[]
    if isempty(ds)
        Label(left[4, 1:4], "no characterised designs in\n$(basename(DEFAULT_DATA))\n($(skipped) raw entries skipped)",
              fontsize = 11, color = col(GuiModel.BAD), halign = :left, justification = :left, padding = (14, 14, 8, 8))
    end
    for (i, d) in enumerate(ds)
        r = 3 + i
        b = Button(left[r, 1], label = d.id, halign = :left, fontsize = 10.5, buttoncolor = :transparent,
                   buttoncolor_hover = col(GuiModel.ACCENT_SOFT), buttoncolor_active = col(GuiModel.ACCENT_SOFT),
                   labelcolor = col(GuiModel.INK), strokewidth = 0, padding = (14, 4, 0, 0))
        push!(rowbuttons, b)
        Label(left[r, 2], d.sigma, fontsize = 10.5, halign = :right, padding = (4, 4, 0, 0))
        Label(left[r, 3], string(d.F), fontsize = 10.5, halign = :right, padding = (4, 4, 0, 0))
        Label(left[r, 4], fmt(d.theta_max), fontsize = 10.5, halign = :right, padding = (4, 14, 0, 0))
        on(b.clicks) do _
            playing[] = false
            theta[] = 0.0
            sel[] = d
        end
    end
    on(sel) do d              # highlight the selected row
        for (b, dd) in zip(rowbuttons, ds)
            b.labelcolor = dd === d ? col(GuiModel.ACCENT) : col(GuiModel.INK)
        end
    end

    # ---- centre: readouts strip, canvas, slider with coloured track, play
    ro = centre[1, 1] = GridLayout(tellwidth = false)
    # state badge: a coloured box under a label, in its own cell (a Box sharing the axis cell
    # would be drawn beneath the axis scene)
    badge_box = Box(ro[1, 1], color = col(GuiModel.ACCENT_SOFT), strokewidth = 0, cornerradius = 6,
                    width = 118, height = 22, halign = :left, alignmode = Outside(0, 0, 4, 4))
    badge = Label(ro[1, 1], "closed", fontsize = 11, color = col(GuiModel.ACCENT), halign = :left,
                  padding = (10, 10, 0, 0), tellwidth = false)
    readout(k, r, c) = (Label(ro[1, c], k, fontsize = 11, color = col(GuiModel.INK2), padding = (8, 2, 6, 6));
                        Label(ro[1, c+1], r, fontsize = 11, font = :bold, padding = (0, 4, 6, 6)))
    sT = lift(d -> d === nothing ? "—" : fmt(d.theta_max; d = 5) * " rad", sel)
    sE = lift(d -> d === nothing ? "—" : fmt(d.eps_max; d = 5) * " rad", sel)
    sF = lift(d -> d === nothing ? "—" : "$(d.F) / $(d.nsplit) / $(length(d.faces) == 0 ? 0 : d.F - d.nsplit)", sel)
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
    Label(centre[2, 1], "drag the slider · click a pattern on the left", fontsize = 10, color = col(GuiModel.INK2),
          halign = :right, valign = :bottom, alignmode = Outside(10), tellwidth = false, tellheight = false)

    bar = centre[3, 1] = GridLayout(tellwidth = false)
    Label(bar[1:2, 1], "Deploy θ", fontsize = 12, padding = (12, 6, 0, 0))
    # coloured track (usable = face-a, certified = accent, over = bad) drawn on a thin axis
    # directly above the slider, in slider-fraction coordinates 0..1
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
    thv = Label(bar[1:2, 3], lift(t -> fmt(deg(t); d = 1) * "°", theta), fontsize = 12, font = :bold, padding = (6, 12, 0, 0))

    # ---- right: actions, exports
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

    Label(right[4, 1], "Export", font = :bold, fontsize = 13, halign = :left, padding = (14, 0, 18, 6))
    xsvg = Button(right[5, 1], label = "SVG at current angle", tellwidth = false, halign = :left, padding = (14, 0, 0, 0))
    xsvg0 = Button(right[6, 1], label = "SVG (closed)", tellwidth = false, halign = :left, padding = (14, 0, 0, 0))
    xjson = Button(right[7, 1], label = "state JSON", tellwidth = false, halign = :left, padding = (14, 0, 0, 0))
    log = Label(right[8, 1], "", fontsize = 10, color = col(GuiModel.INK2), halign = :left, justification = :left,
                padding = (14, 14, 8, 8), tellwidth = false, word_wrap = true)
    function export_file(kind)
        d = sel[]
        d === nothing && return
        mkpath(EXPORT_DIR)
        t = kind == :svg0 ? 0.0 : theta[]
        name = d.id * "_theta" * string(round(Int, deg(t)))
        if kind == :json
            path = joinpath(EXPORT_DIR, name * ".json")
            write(path, GuiModel.state_json(d, t))
        else
            path = joinpath(EXPORT_DIR, name * ".svg")
            write(path, GuiModel.svg_string(d, t))     # TODO(export): Kirigami.export_svg
        end
        log.text = "wrote exports/" * basename(path)
    end
    on(_ -> export_file(:svg), xsvg.clicks)
    on(_ -> export_file(:svg0), xsvg0.clicks)
    on(_ -> export_file(:json), xjson.clicks)

    # legend
    lg = right[9, 1] = GridLayout(halign = :left)
    for (i, (c, txt)) in enumerate(((GuiModel.FACE_A, "+θ/2"), (GuiModel.FACE_B, "−θ/2"),
                                    (GuiModel.ACCENT, "certified"), (GuiModel.BAD, "past Θmax")))
        Box(lg[i, 1], color = col(c), width = 10, height = 10, strokewidth = 0, cornerradius = 2)
        Label(lg[i, 2], txt, fontsize = 11, halign = :left, padding = (4, 0, 1, 1))
    end

    notify(sel)               # initial draw
    return fig, (; sel, theta, playing, designs = ds)
end

if abspath(PROGRAM_FILE) == @__FILE__
    data, shot, theta0 = parse_args(ARGS)
    fig, st = build_app(GuiModel.FileBackend(data))
    theta0 > 0 && (st.theta[] = theta0)
    screen = display(fig)
    if shot !== nothing
        sleep(0.5)
        save(shot, fig)
        println("screenshot written to ", shot)
    else
        wait(screen)
    end
end
