using LinearAlgebra
using Arpack
using SparseArrays
using Expokit

""" 
L = 4
J = 0.7
hx = 0.3
hz = 0
"""
# Define our Pauli Matrices

σx = [0 1; 1 0]
σz = [1 0; 0 -1]
I2 = Matrix{ComplexF64}(I, 2, 2)


# Define a tensor product function

function tensor(matrix_list)
    
    result = matrix_list[1]

    for i in 2:length(matrix_list)
        result = kron(result, matrix_list[i]) # Calculates tensor product for 2 matrices
    end

    return result
end

# Define our Hamiltonian 

function hamiltonian(L, J, hx, hz)

    dim = 2^L
    H = spzeros(ComplexF64, dim, dim) # create a sparce matrix for quicker calculations

    for i in 1:L
        # Nearest neighbour coupling
        σz_i   = [I2 for _ in 1:L] # Array of Identities
        σz_ip1 = [I2 for _ in 1:L]
        σz_i[i] = σz # change Identity at site i to σz, acts site i
        σz_ip1[mod1(i+1, L)] = σz # mod1 for periodicity

        H -= J * sparse(tensor(σz_i) * tensor(σz_ip1))

        # Transverse Field coupling

        σx_i = [I2 for _ in 1:L]
        σx_i[i] = σx

        H -= hx * sparse(tensor(σx_i))

        # Longitudal Field coupling
        # σz_i are already defined above
        H -= hz * sparse(tensor(σz_i))

    end
    return H
end

# for L = 2, J = 1, hx = hz = 0 should give diag matrix [-2, 2, 2, -2]
#println(hamiltonian(2, 1, 0, 0))
# which it does.

# Define our M Operator.
# e^(iBM_α), taking α = x and M_x is sum over all spins, but these just pauli matricies with constants

function Mx(L)
    dim = 2^L
    M = spzeros(ComplexF64, dim, dim)

    for i in 1:L
        σx_i   = [I2 for _ in 1:L] # Array of Identities
        σx_i[i] = σx # change Identity at site i to σz, acts site i

        M += sparse(tensor(σx_i))
    end
    return M
end



# Apply and Evolve; 
function pulse(ψ, B, M)
    #return expv(1im * B, M, ψ)
    return expmv(1im * B, M, ψ)
end

function evolve(H, t, ψ)
    return expmv(-1im * t, H, ψ)
end

# Calculate Ground state:
Ham = hamiltonian(2, 1, 0, 0)
ψ0 = eigs(Ham; nev=1, which=:SR)[2][:,1] # take smallest real eigenvalues eigenvector
#println(ψ0) # Should give a random state that is superpostion between |00> and |11> states. 


# Hit our gs with the pulses and evolve

function two_pulse_evolution(H, M, ψ0, t, τ, B0, Bτ)

    ψ1 = pulse(ψ0, B0, M) # Hit with first pulse B0 at t' = 0
    ψ2 = evolve(H, τ, ψ1) # evolve for a time τ
    ψ3 = pulse(ψ2, Bτ, M) # Hit with second pulse Bτ at t' = τ
    ψ4 = evolve(H, t, ψ3) # evolve until measurment time t' = t + τ

    return ψ4

end

# call 3 pulse function M to mimic paper, with time between pulses t1, t2, t3. First pulse hits at t' = 0
function final_state(H, M, ψ0, t1, t2, t3, B, B0, Bτ)

    ψ1 = pulse(ψ0, B, M) # Hit with first pulse at t' = 0
    ψ2 = evolve(H, t1, ψ1) # evolve for a time t1
    ψ3 = pulse(ψ2, B0, M) # Hit with second pulse at t' = t1
    ψ4 = evolve(H, t2, ψ3) # evolve for a time t2
    ψ5 = pulse(ψ4, Bτ, M) # Hit with third pulse at t' = t1 + t2
    ψ6 = evolve(H, t3, ψ5) # evolve until measurment at t' = t1 + t2 + t3

    return ψ6
end


# Take inner product 

function mx(ψ, Malpha, L)
    expectation =  real(dot(ψ, Malpha * ψ)) / dot(ψ, ψ) # Take expectation and normalise to be safe
    return expectation/L # divide by L to get per site magnetization
end

# Calculate χ^3 using central difference method
# Looking at t1 = 0
function chi3(L, t1, t2, t3, B, J, hx, hz)

    M = Mx(L)
    hamil = hamiltonian(L, J, hx, hz)
    ψ0 = eigs(hamil; nev=1, which=:SR)[2][:,1] # take smallest real eigenvalues eigenvector
    
    if t1 == 0
        # For different Bs, only play affect here so determine different values
        ψpp = final_state(hamil, M, ψ0, t1, t2, t3, B, B, B)
        ψpz = final_state(hamil, M, ψ0, t1, t2, t3, B, B, 0)
        ψpm = final_state(hamil, M, ψ0, t1, t2, t3, B, B, -B)
        ψmp = final_state(hamil, M, ψ0, t1, t2, t3, B, -B, B)
        ψmz = final_state(hamil, M, ψ0, t1, t2, t3, B, B, 0)
        ψmm = final_state(hamil, M, ψ0, t1, t2, t3, B, -B, -B)


        deriv = (
            mx(ψpp, M, L) - 2*mx(ψpz, M, L) + mx(ψpm, M, L) -
            mx(ψmp, M, L) + 2*mx(ψmz, M, L) - mx(ψmm, M, L)
        ) / (2 * B^3)

        return deriv

    else

        ψpp = final_state(hamil, M, ψ0, t1, t2, t3, B, B, B)
        ψzp = final_state(hamil, M, ψ0, t1, t2, t3, B, 0, B)
        ψpm = final_state(hamil, M, ψ0, t1, t2, t3, B, B, -B)
        ψmp = final_state(hamil, M, ψ0, t1, t2, t3, B, -B, B)
        ψzm = final_state(hamil, M, ψ0, t1, t2, t3, B, B, 0)
        ψmm = final_state(hamil, M, ψ0, t1, t2, t3, B, -B, -B)


        deriv = (
            mx(ψpp, M, L) - 2*mx(ψzp, M, L) - mx(ψpm, M, L) +
            mx(ψmp, M, L) + 2*mx(ψzm, M, L) - mx(ψmm, M, L)
        ) / (2 * B^3)

        return deriv
    end
end

# test to see if runs
println(chi3(2, 0, 25, 25, 0.001, 0.7, 0.3, 0)) # this has a limit on L due to memory size.
println(chi3(2, 25, 0, 25, 0.001, 0.7, 0.3, 0))
