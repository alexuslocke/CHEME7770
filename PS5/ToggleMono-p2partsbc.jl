include("PhasePortraitV2.jl")

# Function for a dual repression system without cooperativity
# x1: range of x1 values (i.e. A values)
# x2: range of x2 values (i.e. R values)
# We use `@.` to apply the calculations across all rows.
# Note that model parameters are specified within the function
# Returns computed (dx1/dt, dx2/dt) over the range of (x1, x2)
function toggleMono(ca, cr)
    ror = 1
    rr = 100
    roa = 100
    ra = 5000
    da = 30
    
    u = @. -da*ca + (roa+ra*ca^2)/(1+ca^2+cr^2)
    v = @. -cr + (ror+rr*ca^2)/(1+ca^2)

    return (u,v)
end

#Range of x1, x2 values
x1range = (0,200,25)          #Has the form (min, max, num points)
x2range = (0,100,25)          #Has the form (min, max, num points)
x₀ = ([1.0,10.0],)  #initial state vectors; a common must be included after the first array
tspan=(0.0,30.0)             #time span

#Call the phaseplot functon to construct the phase portrait
phaseplot(toggleMono, x1range, x2range, xinit=x₀, t=tspan, clines=true,
        norm=true, scale=0.25)
