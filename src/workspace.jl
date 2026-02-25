mutable struct BSWorkspace
    W::Matrix{Float64}
    wnorm2::Vector{Float64}
    s::Vector{Float64}
    s_ref::Vector{Float64}
    beta::Vector{Float64}
    dots::Vector{Float64}
    work::Vector{Float64}
end

function BSWorkspace(_m::Int, n::Int, kmax::Int)
    nrem = max(n - 1, 0)
    return BSWorkspace(
        zeros(Float64, kmax, n),
        zeros(Float64, n),
        zeros(Float64, n),
        zeros(Float64, n),
        zeros(Float64, nrem),
        zeros(Float64, nrem),
        zeros(Float64, max(n, 1)),
    )
end

function _require_workspace(
    workspace::Union{Nothing,BSWorkspace},
    m::Int,
    n::Int,
    kmax::Int,
)
    if workspace === nothing
        return BSWorkspace(m, n, kmax)
    end

    ws = workspace
    size(ws.W, 1) >= kmax || throw(ArgumentError("workspace.W has too few rows"))
    size(ws.W, 2) >= n || throw(ArgumentError("workspace.W has too few columns"))
    length(ws.wnorm2) >= n || throw(ArgumentError("workspace.wnorm2 too small"))
    length(ws.s) >= n || throw(ArgumentError("workspace.s too small"))
    length(ws.s_ref) >= n || throw(ArgumentError("workspace.s_ref too small"))
    length(ws.beta) >= max(n - 1, 0) || throw(ArgumentError("workspace.beta too small"))
    length(ws.dots) >= max(n - 1, 0) || throw(ArgumentError("workspace.dots too small"))
    length(ws.work) >= max(n, 1) || throw(ArgumentError("workspace.work too small"))

    # Reuse workspace buffers without eager full zeroing; kernel overwrites active regions.
    return ws
end
