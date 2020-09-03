"""
    Ormsby(; <keyword arguments>)

Create a Ormsby wavelet sampled every dt seconds with corner frequencies
defined by the vector f = [f1, f2, f3, f4]. The final wavelet is multiplied by
a Hamming window.

# Arguments
- `dt::AbstractFloat=0.002`: sampling interval in secs.
- `f::Vector{AbstractFloat}=[2.0, 10.0, 40.0, 60.0]`: corner frequencies in Hz.
      ^
    1 |     ***************
      |    *               *
      |   *                 *
      |  *                   *
      | *                     *
      -----------------------------> f
        f1  f2           f3  f4
# Example
```julia
julia> w = Ormsby(); plot(w);
```
"""
function Ormsby(; dt::Tf=0.002, f::Vector{Tf}=[2.0, 10.0, 40.0, 60.0]) where {Tf<:AbstractFloat}

    f1,f2,f3,f4 = f

    fc = (f2+f3)/2.0
    nw = 2.2/(fc*dt)
    nc = floor(Int, nw/2)
    t = dt*collect(-nc:1:nc)
    nw = 2*nc + 1
    a4 = (pi*f4)^2/(pi*(f4-f3))
    a3 = (pi*f3)^2/(pi*(f4-f3))
    a2 = (pi*f2)^2/(pi*(f2-f1))
    a1 = (pi*f1)^2/(pi*(f2-f1))

    u = [ a4 * (sinc(f4 * t[i] ))^2 - a3 * (sinc(f3 * t[i] ))^2 for i in eachindex(t) ]
    v = [ a2 * (sinc(f2 * t[i] ))^2 - a1 * (sinc(f1 * t[i] ))^2 for i in eachindex(t) ]

    w  = u - v
    Hm = Hamming(nw) 
    w  = [w[i] * Hm[i] for i in eachindex(w)]
    w  = w/maximum(w)

end
