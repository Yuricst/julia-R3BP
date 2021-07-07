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
