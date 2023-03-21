### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	
	using VoronoiFVM
	using ExtendableGrids
	using GridVisualize
	using LinearAlgebra
	using StaticArrays
	using LessUnitful
	using PlutoVista	
	using PlutoUI
	using Colors

	using GasDiffusionElectrodes
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 3ebfa38c-6b69-4108-883b-2c54a453f1c3
html"<button onclick='present()'>present</button>"

# ╔═╡ 863c9da7-ef45-49ad-80d0-3594eca4a189
PlutoUI.TableOfContents(title="Gas Diffusion Electrode")

# ╔═╡ 0c8a42a7-3798-4c64-8455-cab89dd30128
md"""
# Introduction
The following establishes (WIP) a 1D, isothermal, steady state model for gas diffusion electrodes (GDE) used in $\text{CO}$ / $\text{CO}_2$ electrolysis.

This file is part of a public [repository](https://github.com/DavidBrust/GasDiffusionElectrodes) in collaboration with [J. Fuhrmann](https://github.com/j-fu), the main author of the numerical software ([VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl)) used throughout.

The approach follows:

1. [Weng, L. C., et al. (2018).](https://pubs.rsc.org/en/content/articlehtml/2018/cp/c8cp01319e) (A. Z. Weber group, JCAP, Open Access) 
1. [Divisek, J., et al. (2003).](https://iopscience.iop.org/article/10.1149/1.1572150) (FZJ, J. Fuhrmann et al, Open Access)
"""

# ╔═╡ e4b5de7f-be77-426f-9bea-20fc9db02001
md"""
# Modelling Domain
	
Focus on modelling the regions of diffusion medium (DM) and catalyst layer (CL).

$(Resource("https://pubs.rsc.org/image/article/2018/CP/c8cp01319e/c8cp01319e-f1_hi-res.gif"))

"""

# ╔═╡ 2bb77cb7-7323-4539-8e80-d84f9ae53f38
md"""
# Discretization grid
- region 1: catalyst layer (__CL__), phases: __gas__ and __liquid__ and __solid__ phase, species: gas and ionic species

- region 2: diffusion medium (__DM__), phases: __gas__ and __solid__ phase, species: gas phase speceis
$(Resource("https://pubs.rsc.org/image/article/2018/CP/c8cp01319e/c8cp01319e-f3_hi-res.gif"))
"""

# ╔═╡ 2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
md"""
# Governing system
Describe the transport of mass (chemical species) through the porous diffusion medium to facilitate the electro-chemical reaction in the catalyst layer.

Within both porous regions (DM, CL), the existance of a liquid water phase adjacent to a gaseous phase is assumed (gas phase g/ liquid water phase w).
"""

# ╔═╡ aa5facd4-ffac-4b24-8868-d94b1e7b615d
md"""
## Charge transfer reaction
For Ag catalyst at the cathode, the predominant CO2 reduction product is CO. Therefore only CO2r to CO is considered. Literature on CO2r suggests the alkaline reaction to take place, where water (H₂O) is the proton source.
The cathodic reaction then reads:
```math
\text{CO}_2 + \text{H}_2\text{O} + 2 \text{e}^- \rightarrow \text{CO} + 2 \text{OH}^-
```
"""

# ╔═╡ 6c37cb58-d675-4025-a449-13b67bf56899
md"""
## Porosity
The void fraction (fraction of empty space) in porous materials is quantified by the porosity $\phi$.
```math
	\phi = \frac{\text{void volume of porous material}}{\text{total volume of porous material}}
```
"""

# ╔═╡ 6ac23f0b-6b5a-4ac1-9ff4-3d5648fd6901
md"""
The volume void fraction remains when subtraction volume fractions occupied by the solid carbon $\epsilon_{\text s}$ (diffusion medium / catalyst layer) and polymer electrolyte $\epsilon_{\text e}$ (catalyst layer / membrane) phases. 
"""

# ╔═╡ 53a8751a-004b-489d-a4bf-2a2e0ac8e7c8
md"""
```math
	\phi = 1 - \epsilon_{\text s} - \epsilon_{\text e}
```
"""

# ╔═╡ 875448c5-e6cc-4c8e-b53a-6d7411828ebb
function poros(i,data)
	1.0 - data.ϵ_s[i] - data.ϵ_e[i]
end

# ╔═╡ 3ecac2d7-19b5-4528-86c2-5b8340a57841
md"""
When different phases are involved, the available pore space is shared among them (gas phase g/ liquid water phase w). 
"""

# ╔═╡ 17189a30-0a59-470b-9f59-2fa354e79eca
md"""
```math
	\phi_{\text{g/w}} = \frac{\text{volume of phase g/w}}{\text{total volume of porous material}}
```
"""

# ╔═╡ 0f4e8933-d570-484b-8562-a83f2425d9c1
md"""
## Capillary pressure
In case of "wetting" liquid in contact with solid interface, through the intermolecular and interfacial forces, the interface assumes a form that minimizes the total potential energy. The liquid water phase shows the tendency to creep up through pores. This phenomenon is described by capillary pressure that establishes across the phase boundary gas/liquid.

$(Resource("https://upload.wikimedia.org/wikipedia/commons/4/43/CapillaryRise.png"))
"""

# ╔═╡ d1a185dd-7eae-4f10-be33-b4cc845d7b61
md"""
```math
	p_{\text c} = p_{\text{gas, tot} } - p_{\text w}
```
"""

# ╔═╡ a0f452b6-5ca4-40aa-b29a-d8d82fcdbc58
md"""
## Water Saturation and Effective Permeability
The fraction of pore space occupied by each phase (gas phase g/ liquid water phase w) can be quantified by saturation $S_{\text{g/w}}$
"""

# ╔═╡ 34000852-34bd-4bb9-a21c-b7960df63cda
md"""
```math
	S_{\text{g/w}} = \frac{\phi_{\text{g/w}}}{\phi}
```

```math
	\sum_{i=1}^{n_\text{phases}} S_i = 1
```
"""

# ╔═╡ 75d8a575-9899-4110-b34c-3889e1266809
md"""
The water saturation (fraction of pore space occupied by liquid water) is a function of capillary pressure.

Below is an exemplary saturation curve used herein. __This curve describes the wettability properties of the porous medium.__ The relationship must be measured experimentally for the materials used in the electrolyser system. (hydrophobic, carbon based gas diffusion medium)
"""

# ╔═╡ 06b2b42c-f45f-4c01-a220-cf8287d93184
# effective saturation (liquid water phase)
function sat_eff(pc, data)
	B_prime(x) = -(exp(x)*(x-1)+1)/(exp(x)-1)^2	
	#δpc = pc-data.p12	
	se = -B_prime(data.α*(pc-data.p12))
	if isnan(se)
		return 0.5
	else
		return se
	end	
end;

# ╔═╡ 182e4a27-03c7-43e3-9b22-4a8a0482da72
let
	x = collect(range(-2,2,length=101))
	#y = @. sat_eff(x*ufac"bar", ModelData())
	p = plot(xlabel="Capillary pressure / bar", ylabel="Water saturation / -", axisfontsize=16, tickfontsize=16, limits=(0,1), legend=:rt)
	αs = collect((1:2:9)*1.0e-5)
	for α in αs
		y = [sat_eff(x_*ufac"bar", (α=α, p12=0.0)) for x_ in x]
		plot!(p, x,y, label="α=$(α)")
	end
	p
end

# ╔═╡ f922d570-9564-4122-bfd7-3224c17b4188
md"""
The effective relative permeability $k_{\text{rg/w}}$, governs the proportionality between volume flux and driving pressure gradient of fluid flow through porous media (Darcy law).
"""

# ╔═╡ 42e49d99-e68b-44a6-9ed5-ed96b9c84dc4
md"""
```math
	u_{\text{g/w}} = - \frac{k_{\text f} k_{\text{rg/w}}}{\mu_{{\text{g/w}}}} \nabla p_{\text{g/w}}
```
"""

# ╔═╡ 2a860dc8-c798-4d27-8bc7-4d0e40cbe62f
md"""
The relative permeability $k_{\text{rg/w}}$ is a function of phase saturation. For numerical stability, an additional __residual permeability__ parameter $k_{\text{res}}$, one for each phase, is introduced which regularizes the equations (prevents the relative permeabilities from vanishing)

Without the residual permeability parameter, situations could arise where relative permeabilities can become zero. Depending on the choice of the saturation parameter $\alpha$, constant saturation (e.g. complete water saturation) could occur in part of the domain, leading to relative permeabilities of zero which in turn leads to vanishing equations.
"""

# ╔═╡ e8755eed-0d23-45f2-b7ac-1e33e772fd95
md"""
Liquid water phase relative permeability:
```math
	k_{\text{r,w}} = k_{\text{res,w}} + (1-k_{\text{res,w}}) S_{\text{w}}^{\gamma}
```
"""

# ╔═╡ 194c9e3e-4e84-4f9e-bc28-624c1959dcf9
md"""
Gas phase relative permeability:
```math
	k_{\text{r,g}} = k_{\text{res,g}} + (1-k_{\text{res,g}}) (1-S_{\text{w}}^{\gamma})
```
"""

# ╔═╡ e8572b3b-3527-43e6-8f37-603674570667
function kr(se,data)
	krw = data.kres_w + (1.0-data.kres_w) * se^data.γw # liquid water phase w
	krg = data.kres_g + (1.0-data.kres_g) * (1-se^data.γw) # gas phase g
	krg,krw
end

# ╔═╡ cdbe0982-0e0a-4f3c-b3ad-618dc27c0ed4
md"""
A consequence of the residual permeability parameter discussed above is that the relative permeability for either phase does not vanish. To see this zoom into the graph.
"""

# ╔═╡ eab9ff54-9322-4794-b0a9-de9ac11d78f2
md"""
## Liquid Water Phase
"""

# ╔═╡ 17593695-c27f-42b3-b326-cb3fb0242f20
md"""
The liquid water phase is described in terms of the liquid water pressure $p_{\text w}$. Its transport equation takes the form
"""

# ╔═╡ 5cf89fb0-1370-4600-9bbf-a10aaf49adfa
md"""
```math
\begin{align}

\nabla \cdot \mathbf{J}_{\text w} &= R_{\text w} \\
\mathbf{J}_{\text w} &= - \frac{\rho_{\text w} k_{\text f} k_{\text{rw}}}{M_{\text w} \mu_{\text w}} \nabla p_{\text w}

\end{align}

```
"""

# ╔═╡ 1843685f-233e-457b-b57e-30a573374325
md"""
With species molar water fluxes $J_{\text w}$ in units of $\frac{\text{mol}}{\text{m}^2 \text{s}}$ and liquid water source terms $R_{\text w}$ from chemical reactions or phase transitions in units of $\frac{\text{mol}}{\text{m}^3 \text{s}}$.

In future interations, the electro-osmotic water drag resulting from the transport of ions in the catalyst layer can be considered.
"""

# ╔═╡ 574329a5-51eb-4c07-8e31-be9b9e2f5209
md"""
## Dissolved (Solution Phase) Species

Dissolved species ($\text{CO}_2$) occuring in the system are described via molar concentrations $c_i$ ($\frac{\text{mol}}{\text{m}^3}$).

The dissovled species fluxes are induced by a volumetric water flux and dispersion effects ([Divisek, J., et al. (2003)] (https://iopscience.iop.org/article/10.1149/1.1572150)):
"""

# ╔═╡ 3e8661e2-fccc-4c0a-be8c-772f9f8e3bcc
md"""
```math
\begin{align}

\nabla \cdot \mathbf{J}_{\text{sol},i } &= R_{\text{sol},i} \\
\mathbf{J}_{\text{sol},i} &= - D^{\text d}_i \phi s_{\text w} \nabla c_i + c_i \frac{M_{\text w}}{\rho_{\text w}} \mathbf{J}_{\text w}

\end{align}

```
"""

# ╔═╡ f850ff71-da70-4ec9-9f2a-4080585ad0ce
md"""
where $D^{\text d}_i, \phi$ and $s_{\text w}$ denote the dispersion coefficient, porosity (volume void fraction) and water saturation respectively.
"""

# ╔═╡ 06e742d6-280b-4ab7-b8e3-09453bb3b688
md"""
## Gas Phase
At this point, tranport in the gas phase is limited to multi-component (Stefan-Maxwell) diffusion. Currently the interaction of gas phase transport with the surrounding porous medium is not taken into account. This can be handled by incorporating e.g. the dusty-gas-model or the mean transport pore model.
"""

# ╔═╡ f4dcde90-6d8f-4b17-b4ec-367d2372637f
md"""
### Species Mass balances
For a system consisting of $N$ gas phase species, the gas phase composition is expressed via partial pressures $p_{\alpha}, \alpha = 1 ... N$.
"""

# ╔═╡ 3703afb0-93c4-4664-affe-b723758fb56b
md"""
The transport equation for each species $\alpha$ ($\alpha$ = CO, CO2, H2, H2O ... )in the gas phase takes the form
"""

# ╔═╡ 362b9cfc-2de0-4fb3-9110-408f6ab0fb45
md"""
```math
\nabla \cdot \mathbf{J}_{\alpha} = R_{\alpha}
~,
\alpha = 1 ... N
```
"""

# ╔═╡ dd21d46b-8d8a-443f-9e22-d2cf95fe5107
md"""
With species molar fluxes $J_{\alpha}$ in units of $\frac{\text{mol}}{\text{m}^2 \text{s}}$ and species source terms $R_{\alpha}$ from chemical reactions or phase transitions in units of $\frac{\text{mol}}{\text{m}^3 \text{s}}$.
"""

# ╔═╡ eee8e226-b928-4468-af11-37a55108b35a
md"""
The influence of the liquid water phase on the transport of gas phase species can be incorporated via the relative permeability:
"""

# ╔═╡ 932d34bc-bf11-425c-8f9c-e9c516bcc730
md"""
```math
\nabla \cdot \left( k_{\text f} k_{\text{rg}} \mathbf{J}_{\alpha} \right) = R_{\alpha}

```
"""

# ╔═╡ c4e08410-3d1a-4921-aa5c-2f232daaa07a
md"""
### Multicomponent gas species diffusion

Express the Stefan-Maxwell diffusion equations in Matrix-Vector notation, notation taken from __R. Taylor and R. Krishna__, Multicomponent mass transfer. 

"""

# ╔═╡ 649940bd-faa8-4b95-bf62-d82f06935c64
md"""
For ideal gases, the relationship between partial pressures $p_{\alpha}$ and mole fractions $y_{\alpha}$ is

```math
\left. y_{\alpha} = p_{\alpha} \middle/ \sum_{\alpha} p_{\alpha} \right.
```
"""

# ╔═╡ 313cd88a-1497-4d3a-b09a-9b98a6dad9c2
md"""
In a multi-component system of $N$ species, there are only $N-1$ independent equations. The system is completed by the constraint
```math
\sum_{\alpha} y_{\alpha} = 1.
```
The species flux of the last component $J_N$ can be eliminated via
```math
J_N = -\sum_{\alpha=1}^{N-1} J_{\alpha}.
```
"""

# ╔═╡ 8528e15f-cce7-44d7-ac17-432f92cc5f53
md"""
Here we look at ideal gas species at constant temperature $T$. The driving force for transport of species $\alpha$ relative to the other species then is $\nabla p_{\alpha}$. The force-flux relations (without influence of liquid water phase on gas phase transport) read:
"""

# ╔═╡ 19a46c60-d7b5-45da-a3c5-3687316b2a2c
md"""
```math
\begin{align}
\frac{1}{RT} \nabla \mathbf{p} &= - \mathbf{M} ~\mathbf{J} \\
\frac{1}{RT} \nabla p_{\alpha} &= - M_{\alpha \alpha} J_{\alpha} - \sum_{\beta=1 \atop \beta \neq \alpha}^{N-1} M_{\alpha \beta} J_{\beta}
\end{align}
```
"""

# ╔═╡ d5987ad8-e71b-4d40-b1c1-88142558265d
md"""
Equation numbers refer to R. Taylor and R. Krishna, Multicomponent mass transfer.
The coefficients $M_{\alpha \alpha}$ and $M_{\alpha \beta}$ are defined by (2.1.21, 2.1.22):

"""

# ╔═╡ d2adc467-d366-4461-9a00-9832e7d5a737
md"""
```math
	M_{\alpha \alpha} = \frac{x_{\alpha}}{D_{\alpha N}} + \sum_{\gamma=1 \atop \gamma \neq \alpha}^{N} \frac{x_{\gamma}}{D_{\alpha \gamma}}
```
```math
M_{\alpha \beta} = -x_{\alpha} \left ( \frac{1}{D_{\alpha \beta}} - \frac{1}{D_{\alpha N}} \right)
```
"""

# ╔═╡ b6381008-0280-404c-a86c-9c9c3c9f82eb
function M_matrix(y,data)
	n=data.ng
	#D=Dmatrix(n, data.Dg)
	D=data.D
	M=zeros(eltype(y), n-1, n-1)
	for i=1:(n-1)
		M[i,i] = y[i]/D[i,n]
		for j=1:(n-1)
			if j != i
				M[i,j] = - y[i]*(1/D[i,j]-1/D[i,n])
			end
		end
		for k=1:n
			if k != i
				M[i,i] += y[k]/D[i,k]
			end
		end
	end
	#-1.0*M./data.ct
	-1.0*M
end

# ╔═╡ e23d85e3-c0fc-4078-b221-5d3818b193e6
md"""
## Electron transport
It is assumed that electrons are transported through the solid graphite phase. Further assuming the raltionship $\epsilon \approx 1/ \tau$ holds between the solid volume fraction of porous conductor material $\epsilon$ and the turtuosity $\tau$. Then, within the mean-transport pore model ([Divisek, J., et al. (2003)](https://iopscience.iop.org/article/10.1149/1.1572150)), the transport parameter $\psi$ for porous media transport is given by $\psi=\epsilon / \tau=\epsilon^2$.

For electron transport we can use the following potential equation:
"""

# ╔═╡ 1cdfd4f3-38be-4547-bcaa-19d01c181c78
md"""
```math
\begin{align}

\nabla \cdot \mathbf{J}_{\text{e}^-} &= R_{\text{e}^-} \\
\mathbf{J}_{\text{e}^-} &= - \frac{1}{F} \epsilon_{\text s}^2 \sigma_{\text s} \nabla \phi_{\text{e}^-}

\end{align}

```
"""

# ╔═╡ 5e3eb099-f274-4faa-b5db-fe2585afe5a2
md"""
## Ion transport
(WIP) simplified treatment of hydroxide ion (OH⁻) transport. It is assumed that OH⁻ are transported through the solid (polymere) electrolyte (PiperION, carbonate AEM). The same assumptions regarding porous media effictive conductivities as for electron transport applies.
We then can write for ion (OH⁻) transport:
"""

# ╔═╡ 96b04241-9d2d-4eea-9486-384ecc14c7a9
md"""
```math
\begin{align}

\nabla \cdot \mathbf{J}_{\text{OH}^-} &= R_{\text{OH}^-} \\
\mathbf{J}_{\text{OH}^-} &= - \frac{1}{F} \epsilon_{\text s}^2 \sigma_{\text s} \nabla \phi_{\text{OH}^-}

\end{align}

```
"""

# ╔═╡ b46cc553-9cec-46bd-8568-6e1c53692492
md"""
# System Setup and Solution
"""

# ╔═╡ c7ec471a-f1eb-4ed2-b29a-ba805ceda283
md"""
## Boundary Conditions
$(Resource("https://pubs.rsc.org/image/article/2018/CP/c8cp01319e/c8cp01319e-f3_hi-res.gif"))
"""

# ╔═╡ 2f3352b5-4ad4-4718-bc86-6e545083251f
md"""
## Results
"""

# ╔═╡ fcaf66aa-7279-4532-8685-2e0abeae7901
md"""
Catalyst layer thickness: 10 μm

Diffusion Medium thickness: 325 μm

X-Axis: through plane coordinate in meters

Y-Axis: Species partial pressure (gas phase) / Liquid water pressure (liquid phase)

"""

# ╔═╡ 62dc6e42-a45b-446b-b075-98117c32c59e
md"""
Liquid water saturation / -
"""

# ╔═╡ cbb76fa1-25c7-4d63-833d-e6bf88428f8e
md"""
Dissolved species concentrations / mol m⁻³
"""

# ╔═╡ 529cf3b7-0d42-450f-b154-4ac5e7829c7c
md"""
Charged species potentials / V vs. SHE

OH⁻ (hydroxide anion) exists only in the catalyst layer region, where there is polymere electrolyte (ionomer)
"""

# ╔═╡ 455d2b27-154d-41ff-b96d-f41d77a520e1
md"""
## Solver Setup
"""

# ╔═╡ ec327c35-528a-43b1-b618-a51a85c3b3f5
md"""
## Mole Balance Check
Integrate fluxes crossing the boundaries:

The `integrate` method with a test function parameter returns a value for each species, the sign convention assumes that species leaving the domain lead to negative values.


Reaction = In/Outflow
"""

# ╔═╡ c08f2ce6-49e2-41b3-be12-da371977d2ad
md"""
Due to convention in VoronoiVFM.jl the reaction term is on the left hand side of the governing equations. Reaction terms that act as sources (produce species) therefore bear a negative sign (e.g. CO), whereas sink terms (consume species) have a positve sign (e.g. CO2).
"""

# ╔═╡ 1fcb320c-4ca6-47f5-876b-78b44c7484d9
md"""
# Model data
"""

# ╔═╡ ef134982-17d7-4dbc-ac6e-fe957572fcb2
md"""
# Auxiliary
"""

# ╔═╡ 50ac807e-7656-45e1-8414-0d9317a2d3c5
function means(u,data)
	nspec=data.ng+1	# gas phase species + liquid water
	means=zeros(eltype(u),nspec)
	for i=1:nspec
		means[i] = 0.5*(u[i,1]+u[i,2])
	end
	means
end

# ╔═╡ f3dc68cb-f9c3-4e85-892f-34358831e3eb
function Dmatrix(n, D)
	v = zeros(eltype(D),n,n)
	k=1
	for i=1:(n-1)
		for j=(i+1):n
			v[j,i] = D[k]
			k +=1
		end
	end
	Symmetric(v, :L)
end

# ╔═╡ 59bb0503-52d8-4a76-a3f3-ce0fdcfcb170
@Base.kwdef mutable struct ModelData
# diffusion medium thickness
L_DM::Float64 			= 325.0*ufac"μm"	
# catalyst layer thickness
L_CL::Float64 			= 5.0*ufac"μm"
	
# number of gas phase species
ng::Int64		 		= 4
# 1:CO, 2:H2, 3:CO, 4:N2
gn::Dict{Int, String} 	= Dict(1=>"CO", 2=>"H2", 3=>"CO2", 4=>"N2")
# stochiometric coeff of gas species for CT reaction
nug::Vector{Int} 		= [1,0,-1,0] # this will be a matrix for multiple reactions
# stochiometric coeff of charged species for CT reactions
nuc::Vector{Int} 		= [-2,2] # see dict cn for charged species indices
# number of transferred e- per CT reaction turnover
nCT::Int64 		= 2
	
# gas density
ρg::Float64 			= 1.0*ufac"kg/m^3"
# gas dynamic viscosity (use air at 25 °C)
μ_g::Float64 			= 18.37e-6*ufac"Pa*s"
# ambient pressure (absolute)
p_amb::Float64 			= 1.0*ufac"bar"
# pressure scaling parameter
pscale::Float64 		= 1.0*ufac"GPa"

# system temperature (isotherm)
T::Float64 				= 298.15*ufac"K"

# liquid water density	
ρw::Float64 			= 997.0*ufac"kg/m^3"
# liquid water molar mass
Mw::Float64 			= 18.015*ufac"g/mol"
# liquid water dynamic viscosity (at 25 °C)
μw::Float64 			= 0.00089*ufac"Pa*s"
	
# graphite volume fractions of DM and CL (and membrane M) regions
# region 1: CL
# region 2: DM
ϵ_s::Vector{Float64} 	= [0.2, 0.15] # from Divisek
# polymer electrolyte (Nafion/PiperION) volume fractions of DM and CL (and membrane M) regions
ϵ_e::Vector{Float64} 	= [0.21, 1.0e-3] # from Divisek
	
# gas permeability in porous medium
κ::Float64 				= 1.34e-12*ufac"m^2"


# gas phase diffusivity coefficients
# from Weng et. al. (2018)
Dg::Vector{Float64} 	= [0.743,0.152,0.202,0.646,0.779,0.165]*ufac"cm^2/s"
D::Symmetric{Float64, Matrix{Float64}} = Dmatrix(ng, Dg)

# gas phase species molar mass
Mg::Vector{Float64} 	= fill(10.0*ufac"g/mol",ng)

# number of dissolved phase species
nsol::Int64				= 1
#nsol::Int64				= 2
# names of dissolved phase species
soln::Dict{Int, String} = Dict(1=>"CO2")
#soln::Dict{Int, String} = Dict(1=>"CO2", 2=>"OH-")

# number of charged species
nc::Int64 				= 2
# names of charged species
cn::Dict{Int, String} = Dict(1=>"e-", 2=>"OH-")
# charge numbers
z::Vector{Int}   		= [-1,-1]

# reference voltage
E_ref::Float64   		= 0.0*ufac"V"

# solid phase (graphite) electrical conductivity
σ_s::Float64 	= 5000.0*ufac"S/m"

# solid electrolyte phase (ionomer) electrical conductivity
σ_e::Float64 	= 10.0*ufac"S/m"

# volume specific active surface area
#av::Float64 	= 1.0*ufac"1/m"
av::Float64 	= 1.5e7*ufac"1/m"

# quatities for saturation curve
# saturation parameter
α::Float64 		= 5.0e-5*ufac"1/m"
# sturation curve shift
p12::Float64 	= 0.0*ufac"Pa"
# wettable phase permeability exponent
γw::Float64 	= 1.0

# hydraulic permiability
kf::Float64 	= 5.0e-18*ufac"m^2"
# residual permeabilities for liquid water and gas phases
kres_w::Float64 	= 1.0e-3
kres_g::Float64 	= 1.0e-3

# physical constants
e::Float64 = ph"e"
F::Float64 = ph"e"*ph"N_A"
R::Float64 = ph"k_B"*ph"N_A"
end

# ╔═╡ d671f4ca-aa56-481e-b8e5-41f2561feb53
let
	l=101
	x = collect(range(-2,2,length=l))
	krg = zeros(Float64, l)
	krw = zeros(Float64, l)
	
	pg = plot(xlabel="Capillary pressure / bar", ylabel="Gas phase Relative permeability / -", axisfontsize=16, tickfontsize=16, limits=(0,1), legend=:none)

	pw = plot(xlabel="Capillary pressure / bar", ylabel="Liq. water phase Relative permeability / -", axisfontsize=16, tickfontsize=16, limits=(0,1), legend=:none)
	
	αs = collect((1:2:9)*1.0e-5)
	for α in αs
		se = [sat_eff(x_*ufac"bar", (α=α, p12=0.0)) for x_ in x]
		for (i,se) in enumerate(se)
			krg[i], krw[i] = kr(se,ModelData())
		end
		plot!(pg, x,krg, label="α=$(α)")
		plot!(pw, x,krw, label="α=$(α)")
	end
	pg, pw
end

# ╔═╡ 3f643eb7-ff5f-4e88-8a57-5c559687ad1b
begin
	const ngas=ModelData().ng
	const nsol=ModelData().nsol
	const nc=ModelData().nc
	
	# variable indices	
	const iw=ngas+1 	# index liq. water
	const isol=iw+1 	# index first solution phase species
	const ic=isol+nsol 	# index first charged species
	

	# boundary and region numbers
	const Γ_EL_CL=1
	const Γ_DM_GC=2
	const Γ_CL_DM=3
	# bulk regions
	# CL = 1
	# DM = 2
end;

# ╔═╡ fa2002f1-1b35-4fe8-bbd1-f9fff62acbdf
function grid2(;L_CL=5.0*ufac"μm",L_DM=325.0*ufac"μm",nref=0)
	X_CL=range(0,L_CL,length=10*2^nref+1)
	X_DM=range(L_CL,L_CL+L_DM,length=10*2^nref+1)
	grid=simplexgrid(glue(X_CL,X_DM))

	cellmask!(grid,[L_CL],[L_CL+L_DM],2)
    bfacemask!(grid,[L_CL], [L_CL],Γ_CL_DM)
	
    grid
end;

# ╔═╡ 7068922e-c797-416d-ad8e-197d80e3af1a
gridplot(grid2(L_CL=50.0*ufac"μm"),resolution=(600,200),legend=:rt,show=true)

# ╔═╡ f9af078f-6262-4bee-9b1d-04cdf2e9adec
grid=grid2(L_CL=10*ufac"μm",nref=1)

# ╔═╡ 29d66705-3d9f-40b1-866d-dd3392a1a268
function bcond(f,u,bnode,data)
	# gas phase
	boundary_dirichlet!(f,u,bnode,1,Γ_DM_GC,0.0*data.p_amb) # p_CO
	boundary_dirichlet!(f,u,bnode,2,Γ_DM_GC,0.0*data.p_amb) # p_H2
	boundary_dirichlet!(f,u,bnode,3,Γ_DM_GC,1.0*data.p_amb) # p_CO2
	boundary_dirichlet!(f,u,bnode,4,Γ_DM_GC,0.0*data.p_amb) # p_N2

	# liquid water
	boundary_dirichlet!(f,u,bnode,iw,Γ_DM_GC,1.0*data.p_amb) # p_w
	boundary_dirichlet!(f,u,bnode,iw,Γ_EL_CL,0.8*data.p_amb) # p_w	

	# dissolved species
	boundary_dirichlet!(f,u,bnode,isol,Γ_DM_GC,1.0) # c_CO2
	#boundary_dirichlet!(f,u,bnode,isol+1,Γ_DM_GC,1.0) # c_OH-

	boundary_dirichlet!(f,u,bnode,isol,Γ_EL_CL,0.8) # c_CO2
	#boundary_dirichlet!(f,u,bnode,isol+1,Γ_EL_CL,0.8) # c_OH-

	# electron potential (V vs SHE)
	boundary_dirichlet!(f,u,bnode,ic,Γ_DM_GC,-1.0) # ϕe- V vs SHE

	# hydroxide potential (V vs SHE)
	boundary_dirichlet!(f,u,bnode,ic+1,Γ_EL_CL,0.0) # ϕOH- V vs SHE
end

# ╔═╡ 039a6a6b-99ce-49d4-a967-928b55c5e2e2
function yα(p,α)
	p[α]/sum(p)
end

# ╔═╡ ed7941c4-0485-4d84-ad5b-383eb5cae70a
function flux(f,u,edge,data)
	ngas=data.ng
	#p_mean,y_mean=py_means(u,data)
	um = means(u,data)
	F=zeros(eltype(u), ngas-1)
	for i=1:(ngas-1)
		#F[i] = u[i,1]-u[i,2]
		F[i] = (u[i,1]-u[i,2])/(data.R*data.T)
	end
    
	# computation of fluxes J
	
	# total gas pressure
	pt = sum(um[1:ngas])
	yαm = [yα(um[1:ngas],α) for α in 1:ngas] # gas phase species mole fractions
	J = -M_matrix(yαm,data) \ F

	
	pc = pt - um[iw] # capillary pressure
	se = sat_eff(pc, data)
	
	# gas phase species fluxes, closing via ∑pi = p_T
	# f[ngas] via reaction: sum xi = 1
	
	# relative permeability considers pore blocking by liquid water phase
	krg, krw = kr(se,data)	
	f[1:(ngas-1)] = J * krg
	
	# liquid water flux
	f[iw] = data.ρw/data.Mw/data.μw*data.kf*krw*(u[iw,1]-u[iw,2])

	# dissolved species flux
	ϕ = poros(edge.region,data)

	for i=1:nsol
		i_=isol+i-1
	#for i=isol:(isol+nsol-1)
		D_dis = 1.0e-6
		#bp,bm=fbernoulli_pm(f[iw]*data.Mw/data.ρw/(D_dis*ϕ*se))		
        #f[i]=D_dis*ϕ*se*(bm*u[i,1]-bp*u[i,2])
		# use relative permeability for liquid water instead of saturation: it includes regularization
		bp,bm=fbernoulli_pm(f[iw]*data.Mw/data.ρw/(D_dis*ϕ*krw))		
        f[i_]=D_dis*ϕ*krw*(bm*u[i_,1]-bp*u[i_,2])
	end

	# electron flux (mol e- m-2 s-1)
	ie=ic	
	f[ie] = 1/data.F*data.ϵ_s[edge.region]^2*data.σ_s*(u[ie,1]-u[ie,2])


	# hydroxide flux (mol OH- m-2 s-1)
	iH=ic+1
	f[iH] = 1/data.F*data.ϵ_e[edge.region]^2*data.σ_e*(u[iH,1]-u[iH,2])	

end

# ╔═╡ 71e05017-adb5-430f-8477-10b62eff8ba7
function R_CT(u, data::ModelData)
 
	# current density of CT reaction, mA cm-2, reaction rate limited by availability of reactant CO2
	2.5*ufac"mA/cm^2"*yα(u,3)
	
	# species volumetric molar source term from CT reaction, eq. (26)
	#@. data.nug*data.av*rr/data.nCT/data.F	
end

# ╔═╡ 78cf4646-c373-4688-b1ac-92ed5f922e3c
function reaction(f,u,node,data)
	# volumetric charge-transfer reaction in CL region
	if node.region == 1
		rr=R_CT(u,data)
		for i=1:(ngas-1)
			# reaction term is on the left hand side, source term < 0, sink > 0
			f[i]= -rr*data.nug[i]*data.av/data.nCT/data.F	
		end

		for i=1:nc
		 	i_=ic+i-1
			f[i_]= -rr*data.nuc[i]*data.av/data.nCT/data.F	
		end
		
	end
	
	# isobaric conditions
	f[ngas]=log(sum(u[1:ngas])/data.p_amb)
	
	# one would implement a pressure drop in total pressure here
	
end

# ╔═╡ 333b5c80-259d-47aa-a441-ee7894d6c407
begin
	sys=VoronoiFVM.System( 	grid;
							data=ModelData(),
							flux=flux,
							reaction=reaction,
							bcondition=bcond
							)
	
	
	#enable_species!(sys; species=collect(1:(ngas+1)), regions=[1,2])
	
	# include liquid water, dissolved species CO2 and charged species e-, OH-
	#enable_species!(sys; species=collect(1:(ic+1)), regions=[1,2])
	enable_species!(sys; species=collect(1:ic), regions=[1,2])
	# charged species OH- only exists in CL
	enable_species!(sys; species=ic+1, regions=[1])

	inival=unknowns(sys)
	inival[:,:].=1
	inival[1:ngas,:].=ModelData().p_amb/ngas


	
	sol=solve(inival,sys;)
end;

# ╔═╡ 4b41d985-8ebc-4cab-a089-756fce0d3060
let
	# visualization
	vis=GridVisualizer(resolution=(400,200))
	
	c1 = colorant"red"
	c2 = colorant"blue"
	cols=range(c1, stop=c2, length=ngas+1)

	for i in 1:ngas
		scalarplot!(vis, grid, sol[i,:], color=cols[i], label="p$(ModelData().gn[i])", clear=false, legend=:ct)
	end

	scalarplot!(vis, grid, sol[iw,:], color=cols[iw], label="pw", clear=false, legend=:rc, title="Gas & liq. water pressures")

	reveal(vis)
end

# ╔═╡ a631d89a-c22d-43fa-9e64-e1dd6018615e
let
# plot water saturation
	vis=GridVisualizer(resolution=(400,200))
	
	pt = vec(sum(sol[1:ngas,:],dims=1))
	pw = sol[iw,:]
	pc = pt -pw
	se = @. sat_eff(pc, ModelData())
	
	scalarplot!(vis, grid, se, label="water saturation", limits=(0,1), legend=:ct)
	reveal(vis)
end

# ╔═╡ a597bc50-eed7-486a-a434-352752172508
let
	# visualization
	vis=GridVisualizer(resolution=(400,200))
	
	c1 = colorant"red"
	c2 = colorant"blue"
	if nsol > 1
		cols=range(c1, stop=c2, length=nsol)
	else
		cols=[c1]
	end

	for i in 1:nsol
		scalarplot!(vis, grid, sol[isol+i-1,:], color=cols[i], label="c_$(ModelData().soln[i])", clear=false, legend=:ct)
	end

	reveal(vis)
end

# ╔═╡ cb79b54e-9114-49bc-a669-2f1facee6e7f
let
	# plot charged species potentials
	vis=GridVisualizer(resolution=(400,200))
	d=ModelData()

	c1 = colorant"red"
	c2 = colorant"blue"
	if nc > 1
		cols=range(c1, stop=c2, length=nc)
	else
		cols=[c1]
	end

	for i in 1:nc
		if d.cn[i] == "OH-"
			sg=subgrid(grid,[1],)
			scalarplot!(vis, sg, sol[ic+i-1,:], color=cols[i], label="ϕ_$(d.cn[i])", clear=false, legend=:ct)
		else
			scalarplot!(vis, grid, sol[ic+i-1,:], color=cols[i], label="ϕ_$(d.cn[i])", clear=false, legend=:ct)
		end
	end
	reveal(vis)
end

# ╔═╡ 11979bc8-5c15-458e-a408-94f73f6d4a89
begin
	# CO, H2 outflow over Γ_DM_GC boundary
	tf=VoronoiFVM.TestFunctionFactory(sys)

	Γ_where_T_equal_1=[Γ_DM_GC]
	Γ_where_T_equal_0=[Γ_EL_CL]
	T=testfunction(tf,Γ_where_T_equal_0,Γ_where_T_equal_1)
	
	I=integrate(sys,T,sol)
end

# ╔═╡ c76159c7-aed5-44a4-9591-0d0a945ea2fe
R=integrate(sys,reaction,sol)

# ╔═╡ 6795a6c3-44eb-43f3-8969-1d012be9561e
let
	d=ModelData()
	R[7]*d.F
end

# ╔═╡ 6c667865-a5bc-4348-b3f8-8dfe6b3668a1
function ysum(p,data)
	ysum=zero(eltype(p))
	for α=1:data.ng
		ysum+=yα(p,α)
	end
	ysum
end

# ╔═╡ Cell order:
# ╟─3ebfa38c-6b69-4108-883b-2c54a453f1c3
# ╟─863c9da7-ef45-49ad-80d0-3594eca4a189
# ╟─11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
# ╟─0c8a42a7-3798-4c64-8455-cab89dd30128
# ╟─e4b5de7f-be77-426f-9bea-20fc9db02001
# ╟─2bb77cb7-7323-4539-8e80-d84f9ae53f38
# ╠═7068922e-c797-416d-ad8e-197d80e3af1a
# ╠═fa2002f1-1b35-4fe8-bbd1-f9fff62acbdf
# ╟─2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
# ╟─aa5facd4-ffac-4b24-8868-d94b1e7b615d
# ╟─6c37cb58-d675-4025-a449-13b67bf56899
# ╟─6ac23f0b-6b5a-4ac1-9ff4-3d5648fd6901
# ╟─53a8751a-004b-489d-a4bf-2a2e0ac8e7c8
# ╠═875448c5-e6cc-4c8e-b53a-6d7411828ebb
# ╟─3ecac2d7-19b5-4528-86c2-5b8340a57841
# ╟─17189a30-0a59-470b-9f59-2fa354e79eca
# ╟─0f4e8933-d570-484b-8562-a83f2425d9c1
# ╟─d1a185dd-7eae-4f10-be33-b4cc845d7b61
# ╟─a0f452b6-5ca4-40aa-b29a-d8d82fcdbc58
# ╟─34000852-34bd-4bb9-a21c-b7960df63cda
# ╟─75d8a575-9899-4110-b34c-3889e1266809
# ╠═182e4a27-03c7-43e3-9b22-4a8a0482da72
# ╠═06b2b42c-f45f-4c01-a220-cf8287d93184
# ╟─f922d570-9564-4122-bfd7-3224c17b4188
# ╟─42e49d99-e68b-44a6-9ed5-ed96b9c84dc4
# ╟─2a860dc8-c798-4d27-8bc7-4d0e40cbe62f
# ╟─e8755eed-0d23-45f2-b7ac-1e33e772fd95
# ╟─194c9e3e-4e84-4f9e-bc28-624c1959dcf9
# ╠═e8572b3b-3527-43e6-8f37-603674570667
# ╟─cdbe0982-0e0a-4f3c-b3ad-618dc27c0ed4
# ╟─d671f4ca-aa56-481e-b8e5-41f2561feb53
# ╟─eab9ff54-9322-4794-b0a9-de9ac11d78f2
# ╟─17593695-c27f-42b3-b326-cb3fb0242f20
# ╟─5cf89fb0-1370-4600-9bbf-a10aaf49adfa
# ╟─1843685f-233e-457b-b57e-30a573374325
# ╟─574329a5-51eb-4c07-8e31-be9b9e2f5209
# ╟─3e8661e2-fccc-4c0a-be8c-772f9f8e3bcc
# ╟─f850ff71-da70-4ec9-9f2a-4080585ad0ce
# ╟─06e742d6-280b-4ab7-b8e3-09453bb3b688
# ╟─f4dcde90-6d8f-4b17-b4ec-367d2372637f
# ╟─3703afb0-93c4-4664-affe-b723758fb56b
# ╟─362b9cfc-2de0-4fb3-9110-408f6ab0fb45
# ╟─dd21d46b-8d8a-443f-9e22-d2cf95fe5107
# ╟─eee8e226-b928-4468-af11-37a55108b35a
# ╟─932d34bc-bf11-425c-8f9c-e9c516bcc730
# ╟─c4e08410-3d1a-4921-aa5c-2f232daaa07a
# ╟─649940bd-faa8-4b95-bf62-d82f06935c64
# ╟─313cd88a-1497-4d3a-b09a-9b98a6dad9c2
# ╟─8528e15f-cce7-44d7-ac17-432f92cc5f53
# ╟─19a46c60-d7b5-45da-a3c5-3687316b2a2c
# ╟─d5987ad8-e71b-4d40-b1c1-88142558265d
# ╟─d2adc467-d366-4461-9a00-9832e7d5a737
# ╠═b6381008-0280-404c-a86c-9c9c3c9f82eb
# ╟─e23d85e3-c0fc-4078-b221-5d3818b193e6
# ╟─1cdfd4f3-38be-4547-bcaa-19d01c181c78
# ╟─5e3eb099-f274-4faa-b5db-fe2585afe5a2
# ╟─96b04241-9d2d-4eea-9486-384ecc14c7a9
# ╟─b46cc553-9cec-46bd-8568-6e1c53692492
# ╟─c7ec471a-f1eb-4ed2-b29a-ba805ceda283
# ╠═29d66705-3d9f-40b1-866d-dd3392a1a268
# ╠═f9af078f-6262-4bee-9b1d-04cdf2e9adec
# ╟─2f3352b5-4ad4-4718-bc86-6e545083251f
# ╟─fcaf66aa-7279-4532-8685-2e0abeae7901
# ╟─4b41d985-8ebc-4cab-a089-756fce0d3060
# ╟─62dc6e42-a45b-446b-b075-98117c32c59e
# ╟─a631d89a-c22d-43fa-9e64-e1dd6018615e
# ╟─cbb76fa1-25c7-4d63-833d-e6bf88428f8e
# ╟─a597bc50-eed7-486a-a434-352752172508
# ╟─529cf3b7-0d42-450f-b154-4ac5e7829c7c
# ╟─cb79b54e-9114-49bc-a669-2f1facee6e7f
# ╟─455d2b27-154d-41ff-b96d-f41d77a520e1
# ╠═333b5c80-259d-47aa-a441-ee7894d6c407
# ╠═3f643eb7-ff5f-4e88-8a57-5c559687ad1b
# ╠═ed7941c4-0485-4d84-ad5b-383eb5cae70a
# ╠═78cf4646-c373-4688-b1ac-92ed5f922e3c
# ╠═71e05017-adb5-430f-8477-10b62eff8ba7
# ╟─ec327c35-528a-43b1-b618-a51a85c3b3f5
# ╠═11979bc8-5c15-458e-a408-94f73f6d4a89
# ╠═6795a6c3-44eb-43f3-8969-1d012be9561e
# ╟─c08f2ce6-49e2-41b3-be12-da371977d2ad
# ╠═c76159c7-aed5-44a4-9591-0d0a945ea2fe
# ╟─1fcb320c-4ca6-47f5-876b-78b44c7484d9
# ╠═59bb0503-52d8-4a76-a3f3-ce0fdcfcb170
# ╟─ef134982-17d7-4dbc-ac6e-fe957572fcb2
# ╠═50ac807e-7656-45e1-8414-0d9317a2d3c5
# ╠═6c667865-a5bc-4348-b3f8-8dfe6b3668a1
# ╠═f3dc68cb-f9c3-4e85-892f-34358831e3eb
# ╠═039a6a6b-99ce-49d4-a967-928b55c5e2e2
