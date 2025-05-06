function klobuchar(lla::LLA,az::Number,el::Number,gpstime::Number, α::Vector, β::Vector)
    ## Klobuchar model for ionospheric delay

    # check that α β have correct lenghts
    length(α) == 4 || error("α must have length 4")
    length(β) == 4 || error("β must have length 4")

    ## 1. Compute the earth centered elevation angle in the semi-circle units
    ψ = 0.0137 / (el/pi + 0.11) - 0.022
    
    ## 2. Compute the subionospheric latitude in semi-circle units
    Φᵢ = lla.lat / 180 + ψ * cos(az)
    Φᵢ = clamp(Φᵢ, -0.416, 0.416) # clamp to [-0.416, 0.416]

    ## 3. Compute the subionospheric longitude in semi-circle units
    λᵢ = lla.lon /180 + ψ * sin(az) / cos(Φᵢ * pi)

    ## 4. Find geomagnetic latitude in semi-circle units
    Φₘ = Φᵢ + 0.064 * cos(λᵢ * pi - 1.617 * pi)

    ## 5. Find the local time in seconds
    t = 4.32e4 * λᵢ + gpstime
    t = mod(t, 86400)

    ## 6. Find the slant factor
    F = 1 + 16 * (0.53 - el / pi)^3

    ## 7. Find the ionospheric delay
    x = 2 * pi * (t - 50400) / sum([βₙ * Φₘ ^ (n-1) for (n,βₙ) in enumerate(β)])
    T_iono = F * (5e-9 + sum([αₙ * Φₘ ^ (n-1) for (n,αₙ) in enumerate(α)]) * (1 - x^2 /2 + x^4 /24))
end

function klobuchar(lat::Real,lon::Real,alt::Real,az::Real,el::Real,gpstime::Real, α::Vector{Real}, β::Vector{Real}, freq::Real)
    lla = LLA(lat, lon, alt)
    return klobuchar(lla,az,el,gpstime, α, β) * (1575.42e6/freq)^2
end

function klobuchar(lat::Real,lon::Real,alt::Real,az::Real,el::Real,gpstime::Real, α::Vector{Real}, β::Vector{Real})
    lla = LLA(lat, lon, alt)
    return klobuchar(lla,az,el,gpstime, α, β)
end

function klobuchar(lla::Vector{Real},az::Real,el::Real,gpstime::Real, α::Vector{Real}, β::Vector{Real})
    return klobuchar(LLA(lla...),az,el,gpstime, α, β)
end

function klobuchar(lla::Vector{Real},az::Real,el::Real,gpstime::Real, α::Vector{Real}, β::Vector{Real}, freq::Real)
    return klobuchar(lla,az,el,gpstime, α, β) * (1575.42e6/freq)^2
end

function klobuchar(lla::LLA,az::Real,el::Real,gpstime::Real, α::Vector{Real}, β::Vector{Real}, freq::Real)
    return klobuchar(lla,az,el,gpstime, α, β) * (1575.42e6/freq)^2
end

function secondofweek(t::DateTime)
    return 3600*hour(t) + 60*minute(t) + second(t) + 24*3600*(Dates.dayofweek(t) % 7)
end

function secondofweek(t::TimeDate)
    return 3600*hour(t) + 60*minute(t) + second(t) + 24*3600*(Dates.dayofweek(t) % 7)
end