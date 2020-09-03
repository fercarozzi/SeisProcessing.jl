"""
    Berlage(; <keyword arguments>)

Create a Berlage wavelet.

# Arguments
- `dt::AbstractFloat=0.002`: sampling interval in secs.
- `f0::Real=20.0`: central frequency in Hz.
- `m::Real=2`: exponential parameter of Berlage wavelet.
- `alpha::Real=180.0`: alpha parameter of Berlage wavelet in rad/secs.
- `phi0::Real=0.0`: phase rotation in radians.

# Example
```julia
julia> w = Berlage(); plot(w);
```
**Reference**
* Aldridge, David F., 1990, The berlage wavelet: GEOPHYSICS, 55, 1508--1511.
"""
function Berlage(; dt::AbstractFloat=0.002, f0::Real=20.0, m::Real=2, alpha::Real=180.0,
                 phi0::Real=0.0)

    nw = floor(Int, 2.0/(f0*dt))
    t = dt*collect(0:1:nw-1)
    w = [ (t[i] ^ m) * exp(-alpha*t[i]) * cos(2*pi*f0*t[i] + phi0) for i in eachindex(t) ]
    
    w = w/maximum(w)

end
