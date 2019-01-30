"""
    SeisBandPass(in; <keyword arguments>)

Apply a BandPass filter to an input seismic array,
Bandpass filter: flat passband without gain or attenuation, and attenuate all frequencies outside.

# Arguments
* `in`: input seismic array

# Keyword arguments
* `dt=0.004`: time sampling interval in secs.
* `fa=0`: lower frequency cut [Hz]
* `fb=0`: lower frequency banpass [Hz]
* `fc=60`: upper frequency bandpass [Hz]
* `fd=80`: upper frequency cut [Hz]

# Output
Bandpassed seismic data
# Example
```julia
julia> using SeisPlot
julia> d = SeisLinearEvents(dy=0.004); db = SeisBandPass(d,dy=0.004,fc=30,fd=45);
julia> SeisPlotFK(d,dy=0.004)
julia> SeisPlotFK(db,dy=0.004)
```
*Credits: Aaron Stanton,2017*

"""
function SeisBandPass(d;dt=0.004,fa=0,fb=0,fc=60,fd=80)


	nt = size(d,1)
	dn = reshape(d,nt,:)
	nx = size(dn,2)
	nf = iseven(nt) ? nt : nt + 1
	df = 1/nf/dt
	nw = round(Int,nf/2) + 1

	if(fd*dt*nf < nw)
		iw_max = round(Int,floor(fd*dt*nf))
	else
		iw_max = round(Int,floor(0.5/dt))
	end

	d = pad_first_axis(dn,nf)
	m = fft(d,1)/sqrt(size(d,1))
	if fa > 0.
		m[1,:] *= 0.
	end
	for iw=2:iw_max
		f = df*(iw-1)
		if (f<fa)
			m[iw,:] *= 0.
		elseif (f >= fa && f < fb)
			m[iw,:] *= (f-fa)/(fb-fa)
		elseif (f >= fb && f <= fc)
			m[iw,:] *= 1.
		elseif (f > fc && f <= fd)
			m[iw,:] *= 1. - (f-fc)/(fd-fc)
		else
			m[iw,:] *= 0.
		end
	end
	m[iw_max:end,:] .= 0.

	# symmetries
	for iw=nw+1:nf
		m[iw,:] = conj(m[nf-iw+2,:])
	end
	dn = real(bfft(m,1)/sqrt(size(m,1)))
	dout = dn[1:nt,1:nx];
	return reshape(dout,size(d));
end

function pad_first_axis(a,N1)
	n1 = size(a,1)
	nx = size(a[:,:],2)
	b = zeros(N1,nx)
	for ix = 1 : nx
		b[1:n1,ix] = a[:,ix]
	end
	return b
end