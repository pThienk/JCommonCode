# Functions for getting the average shape of avalanches

function resize(vels::Vector, times::Vector, length::Int)

    outtimes = (0:(length-1) |> collect) ./ (length-1) 
    mytimes = (times .- min(times)) ./ (max(times) - min(times))
    
    outvels = [linear_interpolate(mytimes, vels, t) for t in outtimes]

    return outtimes, outvels
end

function find_centers(smin::Real, smax::Real, w::Real, nbins::Int)

    return [(smin/(1-w))*(((1-w)*smax/((1+w)*smin))^(i/(nbins-1))) for i in 0:(nbins-1)]
end

function find_max_width(smin::Real, smax::Real, nbins::Int)
    r = smax/smin
    max_width = (r^(1/nbins) - 1)/(r^(1/nbins) + 1)
    return max_width
end

"""
    Parameters
    ----------
    v : Event velocities, output from get_slips. Required.
    t : Event times, output from get_slips. Required.
    s : Event sizes, output from get_slips. Required.
    d : Event durations, output from get_slips. Required.
    centers: list, optional
        Set the bin centers manually, if you so desire. Otherwise, the program
        will determine the bin centers automatically.
    style : str, optional
        Set the type of averaging to do. Can be 'duration', 'size', or 'duration_shapes'.
        'size' = bin according to size
        'duration' = bin according to duration, then average.
        'duration_shapes' = bin according to duration, but normalize the velocity prior to averaging. Allows the shapes to be observed.
        The default is 'size'.
    width : float, optional
        The fractional half-width of each bin. Events are collected into a bin
        with a center value X_c if an event has X in the range 
        (1-width)*X_c < X < (1+width)*X_c. Default value is ~0.15.
        
        If value is set to None, tells the system to automatically choose
        the width to be maximum, which is not always the best.
        
    ***MAX WIDTH ALLOWED FOR BOUNDARIES S_min, S_max, and N bins:
        |    Let r = S_max/S_min
        |    width_max = (r^(1/N) - 1)/(r^(1/N) + 1)
        |    
        |    **BIN POSITIONS**
        |    Number each bin s_1, s_2, ... s_N. Furthermore, (1-w)*s_1 = s_min and (1+w)*s_N = smax.
        |    Then, each bin position is given by: s_1 = s_min/(1-w), s_(n-1) *(1+w)/(1-w) = s_n, 1 < n <= N.
        
        ***FINDING BIN CENTERS FOR AUTOMATIC BINNING with S_min, S_max, width w, and N bins:
            Number each bin s_1, s_2, ... s_N
            Let s_1 = s_min/(1-w) and s_N = s_max/(1+w)
            s_n = s_1*((s_N/s_1)**(n/N))
            s_n = s_min/(1-w)*((1-w)s_max/((1+w)s_min))**(n/(N-1))
        
    nbins : int, optional
        The number of bins to create. The default is 5.
    mylimits : list, optional
        A list which holds the minimum and maximum X to search for, [min(X), max(X)].
        X can either be duration or size, depending on the style chosen.
        The default is None.
        
    errbar_type: string, optional
        Either 'bootstrap' or 'std'. 'std' for standard deviation (faster but 
        less realistic for points near the edges) and 'bootstrap' for using
        a bootstrapping approach (slow, but far more accurate)
        
    ci: float, optional
        The confidence interval, between 0 and 1. Defaults to 0.95 for 95% CI

    Returns
    -------
    
    allt: 3-dimenisonal list.
        A (nbins)-long list of lists of lists which holds the individual time
        traces for the avalanches which make up each bin.
    allv: 3-dimensional list.
        A (nbins)-long list of lists of lists which holds the individual
        velocity profiles for the avalanches which make up each bin.
    avgt: list of lists of floats.
        Holds the time vectors for each of the averaged profiles.
    avgv: list of list of floats.
        Holds the velocity vector for each of the averaged profiles.
    avgstd: list of list of floats
        Holds the [lo,hi]  95% confidence intervals for each bin.
        As of v4, error bars are identified via bootstrapping.
    centers: list of floats.
        The bin centers. If binning by size, this is the center of each bin in size.
    width: float
        Fractional width width for each bin. Defined as (1-w)*x_c < x_c < (1+w)*x_c
        for center x_c and width w.

"""
function shapes(v, t, s, d; centers=nothing, style::Symbol=:size, width::Real=0.15, nbins::Int=5, limits::Tuple=(), errbar_type::Symbol=:bootstrap, ci=0.95)

    data = nothing
    if style == :size
        data = s
    elseif style == :duration || style == :duration_shapes
        data = d
    else
        @error "Please give a valid option for style! Should be size, duration, or duration_shapes."
        return nothing
    end

    if isempty(limits)
        limits = (min(data), max(data))
    end

    max_width = find_max_width(limits..., nbins)

    if width == -1
        width = max_width
    end

    if width > max_width
        @error "Width is too large! max_width = $(max_width), input width = $(width)"
    end

    if centers === nothing
        centers = find_centers(limits..., width, nbins)
    end

    allt = []
    allv = []
    avgt = []
    avgv = []
    avgstd = []

    @inbounds for center in centers 
        
        bt = [t[idx] - min(t[idx]) for (idx, val) in enumerate(data) if val >= center*(1-width) && val < center*(1+width)]
        bv = [v[idx] - min(v[idx]) for (idx, val) in enumerate(data) if val >= center*(1-width) && val < center*(1+width)]

        lens = length.(bt)
        idxmax = argmax(lens)

        at = zeros(lens[idxmax])
        av = zeros(lens[idxmax])

        toav_v = zeros(length(bt), lens[idxmax])
        for i in eachindex(bv) 
            toav_v[i, 1:length(bv[i])] .+= bv[i]
        end

        if style == :size
            at = bt[idxmax] .- min(bt[idxmax])
        else

            maxlen = length(bt[idxmax])
            tmpt = nothing
            tmpv = nothing
            for i in eachindex(bv)
                avv = bv[i]
                avt = bt[i]
                if style == :duration_shapes
                    avv ./= max(avv)
                end

                tmpt, tmpv = resize(avv, avt, maxlen)
                toav_v[i, :] .= tmpv

            end

            at = tmpt
            at .*= center
        end

        av = mean(toav_v, dims=1)

        astd = []
        if errbar_type == :bootstrap
            for i in 1:max(lens)
                bs_data = bootstrap(mean, toav_v[:, i], BasicSampling(9999))
                (bs_c, ) = confint(bs_data, BasicConfInt(ci))
                lo = av[i] - bs_c[3]
                hi = bs_c[2] - av[i]
                push!(astd, (lo, hi))
            end

        else
            for i in 1:max(lens)
                z = invlogcdf(Normal(), log((1 - ci) / 2))
                sigma = z*std(toav_v[:, i])
                push!(astd, (sigma, sigma))
            end
        end

        push!(allt, bt)
        push!(allv, bv)
        push!(avgt, at)
        push!(avgv, av)
        push!(avgstd, astd)
    end

    return allt, allv, avgt, avgv, avgstd, centers, width
end