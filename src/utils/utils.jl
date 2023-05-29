#=
List of the utils functions
- pseudoticks
    Returns coordinates of the new ticks.

- sortnatural
    Natural sort a string vector accounting for numeric sorting.

=#


"""
    pseudotick(mytick::Vector{Float64}) => Vector{Float64}
    
Returns coordinates of the new ticks. It returns the midpoints 
between each consecutive pairs of value.

## Example

```julia
julia> vticks = [1,4,5,7] |> float;

julia> pseudoticks(vticks)
4-element Vector{Float64}:
 0.5
 2.5
 4.5
 6.0

```
"""
function pseudoticks(myticks::Vector{Float64})
    new_ticks = zeros(size(myticks,1))
    for i in 1:size(myticks,1)
        if i == 1
            new_ticks[i] = myticks[i]/2
        else
            new_ticks[i] = myticks[i-1] + (myticks[i] - myticks[i-1])/2
        end 
    end   
    return new_ticks
end

"""

sortnatural(x::Vector{String}) => Vector(::String)

Natural sort a string vector accounting for numeric sorting.

## Example
```julia
julia> myvec = ["12", "2", "1", "Y", "10", "X"];
julia> sort(myvec)
6-element Vector{String}:
 "1"
 "10"
 "12"
 "2"
 "X"
 "Y"
 julia> sortnatural(myvec)
 6-element Vector{String}:
  "1"
  "2"
  "10"
  "12"
  "X"
  "Y"
```
"""
function sortnatural(x::Vector{String})
    f = text -> all(isnumeric, text) ? Char(parse(Int, text)) : text
    sorter = key -> join(f(m.match) for m in eachmatch(r"[0-9]+|[^0-9]+", key))
    sort(x, by=sorter)
end
