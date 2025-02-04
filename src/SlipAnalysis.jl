# Rewritten code for Dahmen group's slip avalanches data statistical
# analysis. As with the Python version, the get_slips function accounts
# for both displacement-based ("SLIP"), and velocity-based ("SLIP rate")

"""
    Wrapper function for the vectorized and fixed version of get_slips.
    Extracts basic avalanche statistics from provided data.
    Please ready parameter descriptions carefully!

    Parameters
    ----------
    displacement: (Vector; REQUIRED)
        Time series data to be analyzed for avalanches IF data is some accumulated quantity.
        I.e., net displacement or accumulated stress over time.
        An avalanche in this perspective is a "SLIP", in the parlance of the Dahmen Group.
    velocity: (Vector; REQUIRED)
        Time series data to be analyzed for avalanches IF data is some quantity where at each time a new value is acquired.
        I.e., number of spins flipped in one timestep of the random-field Ising model (RFIM) or the number of cell failures in one timestep in a slip model.
        An avalanche in this perspective is a "slip RATE", in the parlance of the Dahmen Group.
    time: (Vector; REQUIRED)
        Time vector in data units.
        Defaults to an index array, i.e., an array ranging from 0 to N - 1, where N is the length of the input data.
    drops: (Bool; OPTIONAL)
        Default value is TRUE.
        Whether to scan the time series for drops in the data.
    threshold: (<:Real; OPTIONAL)
        Default value is 0.
        Number of standard deviations above the average velocity a fluctuation must be before it is considered an event.
        Recommend 0 for most applications.
        Setting this equal to -1 forces a zero-velocity threshold on velocity curve.
            This is useful for simulations since there's little to no noise to be mistaken for an avalanche.
    mindrop: (<:Real; OPTIONAL)
        Default value is 0.
        Minimum size required to be counted as an event.
        Recommend 0 when starting analysis, but should be set higher when i.e. data culling to find the true value of tau.
    threshtype: (String; OPTIONAL)
        Default value is 'median'.
        What type of threshold to use. Options:
        'median'
            -- Uses the median velocity instead of the mean velocity to quantify the average velocity.
            -- Works best in signals with many excursions (i.e., many avalanches) and is not sensitive to outliers.
            -- Threshold is calculated using median absolute deviation (MAD) with this option instead of standard deviation because it more accurately describes the dispersion of the noise fluctuations while ignoring the avalanches.
        'mean'
            -- This is the traditional method.
            -- The threshold is compared to the mean velocity, (displacement[end]-displacement[start])/(time[end]-time[start]).
            -- Threshold is calculated using the standard deviation with this setting.
            -- This method has some major issues in non-simulation environments, as the standard deviation is very sensitive to large excursions from the noise floor
        'sliding_median'
            -- Uses a sliding median of width window_length to obtain a sliding median estimate of the background velocity.
            -- Useful when the background rate of the process is not constant over the course of the experiment.
            -- Threshold is calculated using median absolute deviation with this option and follows the change in average velocity while accurately describing the dispersion of the noise.
    window_size: (Int; OPTIONAL)
        Default value is 1.
        The window size, in datapoints, to be used when calculating the sliding median.
        Should be set to be much longer than the length of an avalanche in your data.
        Jordan found setting window_size = 3% the length of the data (in one case) was good, but this value can be anywhere from just a few hundred datapoints (very short avalanches, many datapoints) to up to 10-20% of the total length of the signal (when the data are shorter but contain several very long avalanches).
        This is worth playing with!

    Returns
    -------
    [0] velocity: list of lists
        -- A list of avalanche velocity curves.
        -- E.g. v[9] will be the velocity curve for the 9th avalanche.
        -- Velocity curve will always begin and end with 0s because the velocity had to cross 0 both times in its run.
    [1] times: list of lists
        -- A list of avalanche time curves.
        -- E.g. t[9] will be the corresponding times at which v[9] velocity curve occured.
        -- The time curve will have the 0s in the velocity curve occur at t_start - ts/2 and t_end + ts/2 as a 0-order estimate of when the curve intersected with zero.
    [2] sizes: list of floats
        -- List of avalanche sizes.
        -- The size is defined as the amount the displacement changes while the velocity is strictly above the threshold.
        -- Each size is corrected by the background rate.
        -- That is, (background rate)*(duration) is removed from the event size.
    [3] durations: list of floats
        -- List of avalanche durations.
        -- The duration is the amount of time the velocity is strictly above the threshold.
        -- For example, an avalanche with velocity profile [0,1,2,3,2,1,0] has a duration of 5 timesteps.
    [4] st: list of ints
        -- List of avalanche start indices.
        -- E.g. st[9] is the index on the displacement where the 9th avalanche occurred.
        -- The start index is the first index that the velocity is above the threshold.
    [5] en: list of ints
        -- List of avalanche end indices.
        -- E.g. en[9] will be the index on the displacement where the 9th avalanche ends.
        -- The end index is the first index after the velocity is below the threshold (accounting for Python index counting).
"""
@inline function get_slips(; disp::Vector=[], vel::Vector=[], time::Vector=[], drops::Bool=true,
     threshold::Real=0, mindrop::Real=0, threshtype::String="median", window_size::Int=101)::Tuple
    
    # Alert the user if no data is given
    @assert !(isempty(disp) && isempty(vel)) "Give at least one of velocity or displacement!"

    # Set the threshtype value based on input or alert the user if no valid choice is given
    thrt::Int = 0
    if threshtype == "mean"
        thrt = 1
    elseif threshtype == "median"
        thrt = 2
    elseif threshtype == "sliding_median"
        thrt = 3
    else
        @assert (false) "Parameter threshtype must be one of the following: \n mean \n median \n sliding_median"
    end

    # Create the time vector if none is given
    if isempty(time)
        time = isempty(disp) ? eachindex(vel) |> collect : eachindex(disp) |> collect
    end

    # Make window_size the nearest odd number for compatability with sliding_median.
    window_size = window_size + (1 - window_size % 2)

    # Velocity input.
    # st and en indices are the start and end in velocity-land.
    # The time vector is shifted forward by 1 index so the displacement and time are same length.
    is_integrated::Bool = false
    if isempty(disp)
        # Displacement start index is index_velocity_begins + 1.
        # Displacement end index is index_velocity_ends + 1.
        disp = cumulative_trapz_int(time, vel)
        prepend!(disp, [disp[1],disp[1]])
        dt = median(diff(time))
        # Make time and displacement the same length.
        time = [time[1] ; time .+= dt]
        is_integrated = true
    end

    # Displacement input.
    # st and en indices are start and end in displacement-land.
    if isempty(vel)
        # Displacement start index is index_velocity_begins + 1.
        # Displacement end index is index_displacement_ends.
        vel = diff(disp) ./ diff(time)
    end

    # If looking at drops, invert the signal.
    disp = (-1)^(drops) .* disp
    vel = (-1)^(drops) .* vel

    return get_slips_core(displacement, velocity, time, threshold, mindrop, is_integrated, threshtype, window_size)
end

function get_slips_core(smoothed::Vector, deriv::Vector, time::Vector, threshhold::Real, mindrop::Real,
     is_integrated::Bool, threshtype::Real, window_size::Real)
    
    # We now take a numeric derivative and get an average
    diff_avg::Vector = [] # empty
    if threshtype == 1
        diff_avg = mean(deriv)
        diff_avg = [diff_avg for _ in eachindex(deriv)]
    # Do pure median
    # Works better at getting the true size
    elseif threshtype == 2
        diff_avg = median(deriv)
        diff_avg = [diff_avg for _ in eachindex(deriv)]
    else
        diff_avg = sliding_median(deriv, window_size)
        diff_avg[1:(window_size ÷ 2)] .= diff_avg[window_size ÷ 2]
        diff_avg[end-(window_size ÷ 2 - 1):end] .= diff_avg[end-(window_size ÷ 2 - 1)]
    end

    # We now set the minimum slip rate for it to be considered an avalanche
    min_diff::Real = 0 # []
    if threshhold == -1
        min_diff = 0 # zeros(length(deriv))
    else
        if threshtype == 1
            min_diff = diff_avg + (std(deriv)*threshhold)
        else
            # Works for both cases of sliding median and overall median.
            # New code using median average difference (see https://en.wikipedia.org/wiki/Median_absolute_deviation)
            # Estimates the standard deviation of the noisy part of the data by assuming noise is normally distributed
            # and using relationship between MAD and standard deviation for normal data.

            # Rewritter's Comment: Not sure where the magic number 1.4826 came from, tbh
            min_diff = diff_avg + (threshhold * mad(deriv .- diff_avg))
        end
    end

    # Now we see if a slip is occurring, i.e. if the derivative is above the minimum value.

    # Shift by one index forward
    slips::Vector = [0 ; diff(deriv .> min_diff)]
    # Velocity start index (first index above 0)
    velocity_index_begins::Vector{Int} = [i for i in eachindex(slips) if slips[i] == 1]
    # Velocity end index (last index above 0)
    velocity_index_ends::Vector{Int} = [i for i in eachindex(slips) if slips[i] == -1]

    # We must consider the case where we start or end on an avalanche.
    # This checks if this is case and if so makes the start or end of the data a start or end of an avalanche.
    if isempty(velocity_index_begins)
        velocity_index_begins = [1]
    end
    if isempty(velocity_index_ends)
        velocity_index_ends = [lastindex(time)]
    end
    if velocity_index_begins[end] >= velocity_index_ends[end]
        push!(velocity_index_ends, lastindex(time))
    end
    if velocity_index_begins[1] >= velocity_index_ends[1]
        prepend!(velocity_index_begins, [1])
    end

    # Correcting for if displacement is integrated or not.

    # In all cases, the displacement index is one more than velocity index begins.
    displacement_index_begins::Vector{Int} = velocity_index_begins
    # + is_integrated
    displacement_index_ends::Vector{Int} = velocity_index_ends

    # Reported start and end index different depending on if displacement or velocity is given.
    index_begins::Vector{Int} = []
    index_ends::Vector{Int} = []
    if is_integrated
        index_begins = velocity_index_begins
        index_ends = velocity_index_ends
    else
        index_begins = displacement_index_begins
        index_ends = displacement_index_ends
    end

    # Now we see if the drops were large enough to be considered an avalanche.
    # Avalanche duration calculated from signal input.
    # First-order approximation assuming the background velocity is constant throughout the avalanche.
    mindrop_correction::Vector = []
    if threshhold != -1
        mindrop_correction = diff_avg[index_begins] .* (time[index_ends .- is_integrated] .- time[index_begins])
    else
        mindrop_correction = zeros(length(index_begins))
    end
    
    index_av_begins::Vector{Int} = index_begins[mindrop .< (smoothed[displacement_index_ends] .- smoothed[displacement_index_begins] .- mindrop_correction)]
    index_av_ends::Vector{Int} = index_ends[mindrop .< (smoothed[displacement_index_ends] .- smoothed[displacement_index_begins] .- mindrop_correction)]

    # Finally we use these indices to get the durations and sizes of the events, accounting for diff().
    slip_durations::Vector = time[index_av_ends] - time[index_av_begins]
    dt = median(diff(time))

    # Mindrop correction term.
    # Term goes to 0 if avalanche is a single time step, leading to LeBlanc et al 2016 definition of size.
    duration_calculation = time[index_av_ends .- is_integrated] .- time[index_av_begins]

    # Increment size calculation by 1 if the displacement was integrated to account for cumtrapz()
    is_step = 1 .* ((slip_durations .<= dt) .* is_integrated .* (index_av_ends .< length(smoothed)))

    # Sizes are more accurately reported by only correcting for background rate, not the rate + threshold
    slip_sizes::Vector = smoothed[index_av_ends .+ is_step] .- smoothed[index_av_begins] .- diff_avg[index_av_begins] .* duration_calculation .* Int(threshhold != -1)
    # time_begins = time[index_av_begins]
    # time_ends = time[index_av_ends]
    time2 = 0.5 .* (time[1:end-1] .+ time[2:end])
    # Sampling time
    tsamp = median(diff(time2))

    velocity::Vector = []
    times::Vector = []
    for k ∈ eachindex(index_av_begins)
        st = index_av_begins[k]
        en = index_av_ends[k]
        mask = st:en-1 |> collect
        if st == en
            mask = st
        end

        # First-order approximation: assume the shape begins and ends at min_diff halfway between the start index and the preceeding index.
        curv = zeros(en - st + 2)
        curt = zeros(en - st + 2)
        curv[2:end-1] .= (deriv[mask] .- min_diff[st])
        curt[2:end-1] .= time2[mask]
        curt[1] = curt[2] - tsamp / 2
        curt[end] = curt[end-1] + tsamp / 2

        velocity = curv
        times = curt
    end

    return (velocity, times, slip_sizes, slip_durations, index_av_begins, index_av_ends)
end

"""
    Identical to get_slips_wrap & get_slips_core in terms of parameters and outputs, but specifically designed for velocity signals & to ensure no negative sizes.
    Use the trapezoidal rule + chopping off all parts of the signal less than the threshold (i.e. velocity < threshold*std(data) or whatever) to get a more accurate view of the size of an event.

    Parameters
    ----------
    time: (Vector; REQUIRED)
        Time vector in data units.
        Defaults to an index array, i.e., an array ranging from 0 to N - 1, where N is the length of the input data.
    velocity: (Vector; REQUIRED)
        Time series data to be analyzed for avalanches IF data is some quantity where at each time a new value is acquired.
        I.e., number of spins flipped in one timestep of the random-field Ising model (RFIM) or the number of cell failures in one timestep in a slip model.
        An avalanche in this perspective is a "slip RATE", in the parlance of the Dahmen Group.
    drops: (Bool; OPTIONAL)
        Default value is TRUE.
        Whether to scan the time series for drops in the data.
    threshold: (<:Real; OPTIONAL)
        Default value is 0.
        Number of standard deviations above the average velocity a fluctuation must be before it is considered an event.
        Recommend 0 for most applications.
        Setting this equal to -1 forces a zero-velocity threshold on velocity curve.
            This is useful for simulations since there's little to no noise to be mistaken for an avalanche.
    mindrop: (<:Real; OPTIONAL)
        Default value is 0.
        Minimum size required to be counted as an event.
        Recommend 0 when starting analysis, but should be set higher when i.e. data culling to find the true value of tau.
    threshtype: (String; OPTIONAL)
        Default value is 'median'.
        What type of threshold to use. Options:
        'median'
            -- Uses the median velocity instead of the mean velocity to quantify the average velocity.
            -- Works best in signals with many excursions (i.e., many avalanches) and is not sensitive to outliers.
            -- Threshold is calculated using median absolute deviation (MAD) with this option instead of standard deviation because it more accurately describes the dispersion of the noise fluctuations while ignoring the avalanches.
        'mean'
            -- This is the traditional method.
            -- The threshold is compared to the mean velocity, (displacement[end]-displacement[start])/(time[end]-time[start]).
            -- Threshold is calculated using the standard deviation with this setting.
            -- This method has some major issues in non-simulation environments, as the standard deviation is very sensitive to large excursions from the noise floor
        'sliding_median'
            -- Uses a sliding median of width window_length to obtain a sliding median estimate of the background velocity.
            -- Useful when the background rate of the process is not constant over the course of the experiment.
            -- Threshold is calculated using median absolute deviation with this option and follows the change in average velocity while accurately describing the dispersion of the noise.
    window_size: (Int; OPTIONAL)
        Default value is 1.
        The window size, in datapoints, to be used when calculating the sliding median.
        Should be set to be much longer than the length of an avalanche in your data.
        Jordan found setting window_size = 3% the length of the data (in one case) was good, but this value can be anywhere from just a few hundred datapoints (very short avalanches, many datapoints) to up to 10-20% of the total length of the signal (when the data are shorter but contain several very long avalanches).
        This is worth playing with!

    Returns
    -------
    [0] velocity: list of lists
        -- A list of avalanche velocity curves.
        -- E.g. v[9] will be the velocity curve for the 9th avalanche.
        -- Velocity curve will always begin and end with 0s because the velocity had to cross 0 both times in its run.
    [1] times: list of lists
        -- A list of avalanche time curves.
        -- E.g. t[9] will be the corresponding times at which v[9] velocity curve occured.
        -- The time curve will have the 0s in the velocity curve occur at t_start - ts/2 and t_end + ts/2 as a 0-order estimate of when the curve intersected with zero.
    [2] sizes: list of floats
        -- List of avalanche sizes.
        -- The size is defined as the amount the displacement changes while the velocity is strictly above the threshold.
        -- Each size is corrected by the background rate.
        -- That is, (background rate)*(duration) is removed from the event size.
    [3] durations: list of floats
        -- List of avalanche durations.
        -- The duration is the amount of time the velocity is strictly above the threshold.
        -- For example, an avalanche with velocity profile [0,1,2,3,2,1,0] has a duration of 5 timesteps.
    [4] st: list of ints
        -- List of avalanche start indices.
        -- E.g. st[9] is the index on the displacement where the 9th avalanche occurred.
        -- The start index is the first index that the velocity is above the threshold.
    [5] en: list of ints
        -- List of avalanche end indices.
        -- E.g. en[9] will be the index on the displacement where the 9th avalanche ends.
        -- The end index is the first index after the velocity is below the threshold (accounting for Python index counting).
"""
function get_slips_vel(; time::Vector{<:Real}, velocity::Vector{<:Real}, drops::Bool=true,
     threshold::Real=0, mindrop::Real=0, threshtype::String="median", window_size::Int=101)::Tuple

    stddev = mad
    avg = nothing # empty Any
    if threshtype == "median"
        avg = median
    end
    if threshtype == "mean"
        stddev = std
        avg = mean
    end

    cutoff_velocity = (avg(velocity) + stddev(velocity) * threshold * Int(threshold != -1)) .* ones(length(velocity))
    if threshtype == "sliding_median"
        cutoff_velocity = sliding_median(deriv, window_size)
        cutoff_velocity[1:(window_size ÷ 2)] .= cutoff_velocity[window_size ÷ 2]
        cutoff_velocity[end-(window_size ÷ 2 - 1):end] .= cutoff_velocity[end-(window_size ÷ 2 - 1)]
        cutoff_velocity = cutoff_velocity .+ stddev(velocity) * threshold * Int(threshold != -1)
    end

    # Treat the velocity by removing the trend such that its centered around zero.
    deriv::vector = []
    if drops
        deriv = -1 .* velocity
    else
        deriv = velocity .- cutoff_velocity
    end

    # Search for rises in the deriv curve.
    # Set all parts of the curve with velocity less than zero to be equal to zero.
    deriv[deriv .< 0] .= 0
    # Get the slips
    slips = [0 ; diff(deriv .> 0)]
    # Velocity start index (first index above 0)
    index_begins::Vector{Int} = [i for i in eachindex(slips) if slips[i] == 1]
    # Velocity end index (last index above 0)
    index_ends::Vector{Int} = [i for i in eachindex(slips) if slips[i] == -1]

    if isempty(index_begins)
        index_begins = [1]
    end
    if isempty(index_ends)
        index_ends = [lastindex(time)]
    end
    if index_begins[end] >= index_ends[end]
        push!(index_ends, lastindex(time))
    end
    if index_begins[1] >= index_ends[1]
        prepend!(index_begins, [1])
    end

    # Get the possible sizes
    possible_sizes::Vector = []
    possible_durations::Vector = []
    for i ∈ eachindex(index_begins) 
        st = index_begins[i]
        en = index_ends[i]
        trapz_st = max(st - 1, 0)
        trapz_en = min(en + 1, length(deriv))
        push!(possible_sizes, trapz(time[trapz_st:trapz_en], deriv[trapz_st:trapz_en]))
        push!(possible_durations, time[en] - time[st])
    end

    idxs::Vector{Int} = [i for i in eachindex(possible_sizes) if possible_sizes[i] < mindrop]
    sizes::Vector = possible_sizes[idxs]
    durations::Vector = possible_durations[idxs]
    index_av_begins::Vector{Int} = index_begins[idxs]
    index_av_ends::Vector{Int} = index_ends[idxs]

    time2 = 0.5 .* (time[1:end-1] .+ time[2:end])
    # Sampling time
    tsamp = median(diff(time2))
    push!(time2, time2[end] + tsamp)

    velocity_out::Vector = []
    times_out::Vector = []
    for k ∈ eachindex(index_av_begins)
        st = index_av_begins[k]
        en = index_av_ends[k]
        mask = st:en-1 |> collect
        if st == en
            mask = st
        end

        # First-order approximation: assume the shape begins and ends at min_diff halfway between the start index and the preceeding index.
        curv = zeros(en - st + 2)
        curt = zeros(en - st + 2)
        curv[2:end-1] .= deriv[mask]
        curt[2:end-1] .= time2[mask]
        curt[1] = curt[2] - tsamp / 2
        curt[end] = curt[end-1] + tsamp / 2

        velocity_out = curv
        times_out = curt
    end

    return (velocity_out, times_out, sizes, durations, index_av_begins, index_av_ends)
end