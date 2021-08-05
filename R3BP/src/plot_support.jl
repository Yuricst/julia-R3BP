"""
Plot heplers
"""



"""
    list_to_plotarray_3d(sollist, state_idx_x, state_idx_y, state_idx_z)

Function for plotting from list of states in 3D
"""
function sol_to_arrays(sollist)
    nsv = length(sollist[1])
    sols = zeros(nsv, length(sollist))
    for k = 1:length(sollist)
        for i = 1:nsv
            sols[i, k] = sollist[k][i]
        end
    end
    return sols
end


"""
Get circle cartesian coordinates
"""
function get_circle(center, radius, n=100)
    sols = zeros(2,n)
    thetas = LinRange(0.0, 2π, n)
    for i = 1:n
        sols[1,i] = center[1] + radius*cos(thetas[i])
        sols[2,i] = center[2] + radius*sin(thetas[i])
    end
    return sols
end


#
# function plot_manifold!(myplot, manif, c="navy", lw=0.25)
#     for branch in manif
#         branch_arr = sol_to_arrays(branch)
#         plot!(myplot, branch_arr[1,:], branch_arr[2,:], c=c, lw=lw)
#     end
# end
