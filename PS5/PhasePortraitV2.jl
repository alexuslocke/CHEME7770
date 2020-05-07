using DifferentialEquations
using Plots
gr(size=(500,500), dpi=300)  #use the gr backend; set resolution to 300 dpi

############################  PHASEPLOT FUNCTION ##############################


# Function for a construction of a phase portrait
# Output of phaseplot is "phase_portrait.png" saved in the current directory
# model: functon that returns dx1/dt and dx2/dt for the ODE system
# x1lim: boundaries for the x1-phase space; has the form (min, max, num points)
# x2lim: boundaries for the x1-phase space; has the form (min, max, num points)
# clines: set to true to display system nullclines
# xinit: tuple of arrays of initial conditions, [x1_init, x2_init]; a
# trajectory will be plotted for each initial condition array in xinit
# t: time span (t_start, t_end) for the trajectory; trajectory will be drawn
# for times from t_start to t_end
# norm: when set to true, (dx1/dt, dx2/dt) vectors will be normalized by their
# magnitudes
# scale: scaling factor for (dx1/dt, dx2/dt) vectors
function phaseplot(model, x1lim, x2lim; clines=true, xinit=(), t=(0.0,50.0),
    norm=true, scale=0.5)

    #Create a grid of points spanning the specified x1 and x2 ranges
    x1_spacing = (x1lim[2] - x1lim[1])/x1lim[3]
    x2_spacing = (x2lim[2] - x2lim[1])/x2lim[3]

    x1_range = x1lim[1]:x1_spacing:x1lim[2]
    x2_range = x2lim[1]:x2_spacing:x2lim[2]

    x, y = meshgrid(x1_range, x2_range)

    #Compute dx1/dt and dx2/dt for each (x1, x2) grid point. Normalize the length
    #of dx1/dt and dx2/dt
    u,v = model(x, y)

    #Create the phase portrait for the model using a quiver plot.
    if norm
        #Compute a scaling factor and normalize the vector length
        factor = min(x1_spacing, x2_spacing)*scale
        u_norm,v_norm = normalize(u,v).*factor
        #Plot
        plt1 = quiver(x, y, quiver=(u_norm, v_norm), c=:grey, xaxis=("x1",
            (x1lim[1], x1lim[2])), yaxis=("x2",(x2lim[1], x2lim[2])))
    else
        plt1 = quiver(x, y, quiver=(u.*scale, v.*scale), c=:grey, xaxis=("x1",
            (x1lim[1], x1lim[2])), yaxis=("x2",(x2lim[1], x2lim[2])))
    end

    #Plot the x1- and x2- nullclines unless the nullclines argument is
    #set to false.  Rather than solving for the set of (x1, x2)
    #values for which dx1/dt = 0 or dx2/dt = 0, we create a contour plot of
    #computed dx1/dt and dx2/dt values over our grid and display the level set
    #of zero values (i.e. using the levels=[0] argument in countour).  Note that
    #with inclusion of the exclamation point symbol after the countour fxn,
    #the plots will be constructed on the quiver plot rather than on new plots
    if clines
            contour!(x1_range, x2_range, u, lw=2, c=:red, levels=[0],
                    colorbar = :none)
            contour!(x1_range, x2_range, v, lw=3, linestyle=:dot, c=:red,
                    levels=[0], colorbar = :none, xaxis=("x1", (x1lim[1], x1lim[2])),
                    yaxis=("x2",(x2lim[1], x2lim[2])))
    end


    #Solve the model and plot the trajectories for any specified initial
    #conditions
    for i in xinit
        prob = ODEProblem(ODE_model!,i,t,model)
        sol = solve(prob)

        plot!(sol, vars=(1,2), lw=2, c=:black, legend = :none, xaxis=("x1",
            (x1lim[1], x1lim[2])), yaxis=("x2",(x2lim[1], x2lim[2])))
    end

    #Display and save the phase portrait
    savefig(plt1, "./phase_portrait.png")
    display(plt1)

    #Our work is done
    println("Construction of phase portrait completed!")

end

#################  ACCESSORY FUNCTIONS FOR PHASEPLOT  #########################

#Function analogous to Matlab meshgrid
#x: range of x values (i.e. min:spacing:max)
#y: range of y values
#Returns a grid of x and y values over the specified ranges
meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))


#Function to normalize vectors with components (xi,yi).
#x: array of x-vector components
#y: array of y-vector components
#@. is used to apply the calculations across all rows.
#Returns the normalized x- and y- vector components
function normalize(x,y)
    mag = @. hypot(x,y)
    return (x./mag, y./mag)
end

#Function to wrap the model function into a form that can be solved numerically
#with the DifferentialEquations.jl package.
function ODE_model!(du,u,model,t)
    du[1], du[2] = model(u[1], u[2])
end
