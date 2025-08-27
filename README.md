A julia file to "recreate" the main ideas behind this paper https://doi.org/10.1103/PhysRevB.110.134443.
This paper studied the non linear magnetic response using 2D coherent spectroscopy. I followed the outline of this paper and looked at the 1D Transverse Field Ising Model (TFIM)

# Hamiltonian
The model we consider is the 1D-TFIM in the presence of
a longitudinal field with the Hamiltonian

$$H = -J\sum_{i=1}^L\sigma_i^z\sigma_{i+1}^z - h_x\sum_{i=1}^L\sigma_i^x - h_z\sum_{i=1}^L\sigma_i^z$$

where $\sigma_i^{\alpha}$ $(\alpha = x,y,z)$ are Pauli matrices at site $i$, $J>0$ is the nearest neighbour 
ferromagnetic coupling $h_x$ is the transverse field and $h_z$ is the longitudinal field.

# Pulses
The goal is to measure the nonlinear susceptibilites that arise due to different pulses at specific times
$t{\prime}$. The code is aimed to be general-ish with the pulse function, called `final_state`. This function takes in an
initial state, usually a ground state, and hits the state with pulses of strenght $B_0, B_{t_1}$ and $B_{t_2}$ at
$t{\prime} = 0, t{\prime} = t_1, t{\prime} = t_1 + t_2$ respectively, before the magnetisation is measure at a time 
$t{\prime} = t_1 + t_2 + t_3$.The function takes in the time difference between pulses and measurement, $t_i$ and allows for choice in pulse
strength. So what does this function actaully do in practice. 

Using `Expokit.jl` we hit our intial state with a pulse $B_0$, evolve for $t_1$, hit with our second pulse 
$B_{t_1}$, evolve for $t_2$, hit with pulse 3 $B_{t_2}$ and evolve for $t_3$, when measurement occurs. Pulse
action on an arbitrary state $\ket{\psi}$ is 

$$\ket{\psi \prime} = \exp{(iB_tM^{\alpha})}\ket{\psi}$$

where $M^{\alpha}(t) = \sum_iS^{\alpha}_i(t)$ is the $\alpha$ component of the spin operator on site i at time
t.

And evolution is done using the standard time evolution operator

$$U(t) = \exp{(-iHt)}$$

Putting all this together, for a 3 pulse system, our final state at measurment is:

$$ \ket{\psi(t_1, t_2, t_3)} = U(t_3)e^{iB_{t_2}M^{\alpha}}U(t_2)e^{iB_{t_1}M^{\alpha}}U(t_1)e^{iB_{t_0}M^{\alpha}}\ket{\psi} $$

# Magnetization and Susceptibilites

We define magnetization per site as:

$$ m^{\alpha}(t_1, t_2, t_3, B_0, B_{t_1}, B_{t_2}) = \frac{1}{L}\bra{\psi(t_1, t_2, t_3)}M^{\alpha}
\ket{\psi(t_1, t_2, t_3)} $$

To mimic real life experiments and for ease of calculations we often take the semi-impulsive limit by having 2 pulses
occur at the same time with the same strength. Therefore we have $t$ for the time between the final pulse and measurment and $\tau$ for the time
between intermediate pulses. This now leads to 2 choices for when our pulses can occur together. We can have the first two pulses
occur together, $t{\prime} = t_1 = 0, B_0 = B_{t_1}$ or the second and third pulses occur together $t_2 = 0, B_{\tau}$. 

We can now define the non linear susceptibilites as follows:

$$ \chi_{\alpha \alpha \alpha \alpha} ^{(3)}(t, \tau, 0) = \frac{\partial^3m^{\alpha}(t, \tau, B_{\tau}, B_{0})}{\partial B_{\tau} \partial B_0^2} $$

$$ \chi_{\alpha \alpha \alpha \alpha} ^{(3)}(t, 0, \tau) = \frac{\partial^3m^{\alpha}(t, \tau, B_{\tau}, B_{0})}{\partial B_{\tau}^2 \partial B_0} $$

These derivatives are evaluated at $B_0 = B_{\tau} = 0$ and are computed numerically using the central difference method.

# Code Woes 

Despite being pretty cool and julia being a fast language, this code needs to be done on some sort of monster computer for any
reasonable results. This is because Exact Diagonalisation is computationally intense. If you want to see cool figures, plots
or what this stuff actually means, check out the original paper.
