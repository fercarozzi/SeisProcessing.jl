"""
    Ricker(; <keyword arguments>)

Create a Ricker wavelet.

# Arguments
* `dt=0.002`: sampling interval in secs.
* `f0=20.0`: central frequency in Hz.

# Examples
```julia
julia> w = Ricker(); plot(w);
julia> w = Ricker(dt=0.004, f0=20); plot(w);
```
# Reference
Sheriff, Robert, 2002, Encyclopedic Dictionary of Applied Geophysics, fourth
ed.: Society of Exploration Geophysicists. Geophysical Reference Series No. 13.
"""
function Ricker(; dt::Tf=0.002, f0::Tr=20.0)where {Tf<:AbstractFloat, Tr<:Real}

    nw = 2.0/(f0*dt)
    nc = floor(Int, nw/2)
    t = dt*collect(-nc:1:nc)
    b = [ (pi*f0*t[i])^2 for i in eachindex(t) ]
    w = [(1.0 - 2.0 *b[i]) * exp(-b[i]) for i in eachindex(t) ]

end
