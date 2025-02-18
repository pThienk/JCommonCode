# File dedicated to the construction of CDF, AND CCDF (despite the name)
# of general (simulated or experimental) data's event "sizes" (to be defined)
# and their frequency.

"""
HISTORY
-----------------------------------------
Original CCDF code: Alan Long Jun 10 2016 \n
Vectorized & edited by Jordan on Dec 19 2023 \n
Edited by Ethan on Jun 12 2024 \n
Edited slightly by Jordan on June 20 2024 (make the inputs to the CCDF function required) \n
Julia rewrite by Porpun 2/3/2025
____________________________________________________________________________________________
This code takes data as a Vector and returns its complementary cumulative distribution function (CCDF).
- Returns two lists: histx and histy, the x and y values for the ccdf as a Tuple.
- If input lists have negative, NaN, or inf values, it will throw them out and proceed normally.
- Has two separate functions:
[1] ccdf(): This is the standard method. Takes in the full list of data.
[2] ccdf2(): This method allows you to send in the unique values of your data and their counts. Useful if you have a
             very large data array that is best kept as unique values & their counts. These do not need to be sorted,
             but they must have the same length and each data value must have its corresponding count at the same index.
- There are two inputs methods for both functions:
[1] 'scipy': CCDF is defined as P(X > x). The data is appended with a zero at the beginning such that
             P(X > 0) = 1 is the first pair in histx & histy. [2] 'dahmen': CCDF is defined as P(X >= x).
             The data is unchanged, so the first pair in histx & histy is
            P(X >= [smallest array value]) = 1.
"""
function ccdf(data::Vector{<:Real}; method="scipy")::Tuple
    
    @assert (!isempty(data)) "Data cannot be empty!"
    @assert (method == "scipy" || method == "dahmen") "Please choose between two methods: \'scipy\' or \'dahmen\'."

    histx = [value for value in data if (value > 0) && (isfinite(value))]

    counts = unique_count!(histx)

    cumulative_counts::Vector{Int} = []
    histy::Vector = []
    if method == "scipy"
        
        prepend!(histx, [0])
        cumulative_counts = accumulate(+, counts)

        total_count = cumulative_counts[end]
        histy = [1 ; 1 .- (cumulative_counts ./ total_count)]
    else

        prepend!(histx, [0])
        cumulative_counts = accumulate(+, counts)

        total_count = cumulative_counts[end]
        histy = 1 .- (cumulative_counts ./ total_count)
    end

    return (histx, histy)
end

"""
    Indentical to ccdf() in every aspect, except it returns the cumulative distribution
    instead of the complementary version.
"""
function cdf(data::Vector{<:Real}; method="scipy")::Tuple
    @assert (!isempty(data)) "Data cannot be empty!"
    @assert (method == "scipy" || method == "dahmen") "Please choose between two methods: \'scipy\' or \'dahmen\'."

    histx = [value for value in data if (value > 0) && (isfinite(value))]

    counts = unique_count!(histx)

    cumulative_counts::Vector{Int} = []
    histy::Vector = []
    if method == "scipy"
        
        prepend!(histx, [0])
        cumulative_counts = accumulate(+, counts)

        total_count = cumulative_counts[end]
        histy = [0 ; (cumulative_counts ./ total_count)]
    else

        prepend!(histx, [0])
        cumulative_counts = accumulate(+, counts)

        total_count = cumulative_counts[end]
        histy = (cumulative_counts ./ total_count)
    end

    return (histx, histy)
end