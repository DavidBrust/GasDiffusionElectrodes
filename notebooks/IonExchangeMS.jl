### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 23dc6930-c354-11ed-110a-235339eca423
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	
	using VoronoiFVM
	using ExtendableGrids
	using GridVisualize
	using LinearAlgebra

	using LessUnitful
	using PlutoVista	
	using PlutoUI
	using Colors

	using GasDiffusionElectrodes
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 5ede4875-1dbf-4930-9538-39b10d4bd377
PlutoUI.TableOfContents(title="Ion Exchange Maxwell-Stefan")

# ╔═╡ 1593e889-7b4f-40f9-b6e8-ef0564381c81
md"""
__Wesselingh, J. A., P. Vonk and G. Kraaijeveld (1995)__. "Exploring the Maxwell-Stefan description of ion exchange." The Chemical Engineering Journal and the Biochemical Engineering Journal 57(2): 75-89.
"""

# ╔═╡ 711a3678-e3c0-4193-a6cf-0d525ddeff91
md"""
# Continuity equation
For the species concentrations we assume a continuity equation involving the particle fluxes $N_i$ and species source/sink terms via electro-chemical reactions $R_i$ (for steady state):

```math
\nabla\cdot N_i + R_i= 0
```
"""

# ╔═╡ b55b90f2-bace-485c-a745-3fca9871152f
md"""
# Maxwell-Stefan Formulation
"""

# ╔═╡ 47fea77b-d605-40ac-b436-078fd91294fe
md"""
For a system consisting of $\nu$ species (ions and solvent, usually water), $\nu -1$ Maxwell-Stefan equations for multi-component mass transport can formulated:
"""

# ╔═╡ 4a33c6a3-ad59-4efc-997e-09ca780e5fd2
md"""
```math
	-RT \nabla \ln (\gamma x)_i - Fz_i \nabla \phi - v_i \nabla p = \sum_{i=1 \atop j \neq i}^{\nu} x_j\frac{RT}{D_{ij}}(v_i -v_j)
```
"""

# ╔═╡ 2aefd6b8-f35c-419e-a2c4-83ad44f1be77
md"""
From the article below, the following version is applied (neglecting the pressure gradient):
"""

# ╔═╡ e7177289-5c00-44f0-9ad3-92681d02c05d
md"""
## Aqueous Electrolyte
"""

# ╔═╡ 67a38c46-a952-473d-864e-53ca99b16174
md"""
__Kraaijeveld, G., V. Sumberova, S. Kuindersma and H. Wesselingh (1995)__. "Modelling electrodialysis using the Maxwell-Stefan description." The Chemical Engineering Journal and the Biochemical Engineering Journal 57(2): 163-176.
"""

# ╔═╡ a4bb4cff-e933-45ac-8ad2-c2f69328f62d
md"""
```math
-\frac{c_i}{RT} \nabla \mu_i - \frac{z_i c_iF}{RT}\nabla \phi = \sum_{i=1 \atop j \neq i}^{\nu} \frac{c_i c_j}{c_{\text{tot}} D_{ij}}(v_i-v_j)
```
"""

# ╔═╡ 99faa887-0c70-4834-80e0-9f4e39b1b00a
md"""
When we want to consider different species molar volumes $v_i$ or ion solvation effects (via $\kappa_i$), a formulation based on molar densities ("concentrations") in $\frac{\text{mol}}{\text{dm}^3}$ is preferable (see above )

"""

# ╔═╡ 87c2c676-9bf0-4716-8fe4-e44a0f884412
md"""
Relation between average species velocity and species flux:
"""

# ╔═╡ 72514d88-d17b-4dfe-8921-623aa6e3b42f
md"""
```math
N_i = c_i v_i
```
"""

# ╔═╡ 3ede7aa4-d924-4825-b0cf-f5f9b52c5e3b
md"""
Electro-neutrality condition
"""

# ╔═╡ 21d032ab-3de3-400a-b6fb-2c9805e4916b
md"""
```math
\sum_{i=1}^{\nu} z_i c_i = 0
```
"""

# ╔═╡ ac882611-643c-4c94-9acd-38de7566498e
md"""
Volume-fixed reference frame (incompressible liquid electrolyte):
"""

# ╔═╡ 0dbfa843-7b38-41bf-aa4c-14d528bbf75f
md"""
```math
\sum_{i=1}^{\nu} c_i v_i = 1
```
"""

# ╔═╡ 25ba905a-3342-4139-935d-ebcf0ad4be56
md"""
## Free-Solution Diffusivities
"""

# ╔═╡ 10f5b322-e52b-45b7-ae60-a9f5d82b5f99
md"""
$(LocalResource("../img/FreeSolDiff.png"))
Table 4 from __Kraaijeveld, G., V. Sumberova, S. Kuindersma and H. Wesselingh (1995)__. "Modelling electrodialysis using the Maxwell-Stefan description." The Chemical Engineering Journal and the Biochemical Engineering Journal 57(2): 163-176.
"""

# ╔═╡ e93d0059-0b6b-4911-a148-0c244ea010c4
md"""
## Ion Exchange Membranes
The Maxwell-Stefan (MS) approach lends itself for the description of ion exchange phenomena in polymere ion exchange membranes. The polymere __m__embrane possess fixed background charges. The membrane and its fixed background charges are therefore treated as an additional species denoted "m", which is stationary, $N_{\text m} = 0$.
$(LocalResource("../img/IonExchangeMemSpecies.png"))

Fig. 13 (coloring added) from __Wesselingh, J. A., P. Vonk and G. Kraaijeveld (1995)__. "Exploring the Maxwell-Stefan description of ion exchange." The Chemical Engineering Journal and the Biochemical Engineering Journal 57(2): 75-89.
"""

# ╔═╡ 27b0be8c-e57f-4396-b132-b4ddd64a1892
md"""
Not all component pairs and therefore diffusivities are equally important, depending on the situation.
"""

# ╔═╡ 878dbe91-b305-4155-a98d-c4aced05e8f3
md"""
__Key challenge__: determination of the different friction coefficients (inverse of diffusivites) among the ions and the membrane involved.
"""

# ╔═╡ 973a9c4e-b7c3-4640-96aa-b8f09ebe0b28
md"""
Different experiments and combinations thereof for the calculation of friction coefficients of ions in polymere membranes are presented in:
 __Wesselingh, J. A., P. Vonk and G. Kraaijeveld (1995)__. "Exploring the Maxwell-Stefan description of ion exchange." The Chemical Engineering Journal and the Biochemical Engineering Journal 57(2): 75-89.
$(LocalResource("../img/ExperimentsDetermineFrictionCoeffs.png"))
"""

# ╔═╡ 49d41a14-9f28-4130-bb6b-a03159336ef1
md"""
# Simplified case
"""

# ╔═╡ 22f39a64-bfb9-464e-bd2e-24fc72595749
md"""
When we neglect ion solvation effects, and assume equal partial molar volumes $v_i$ for each species, the total concentration $c_{\text{tot}}$ is constant throughout the solution. Therefore $x_i=\frac{c_i}{c_{\text{tot}}}$. Then the composition can be expressed via mole fractions $x_i$:
"""

# ╔═╡ 6d6d0af1-4185-43db-a5a6-2faad0592668
md"""
```math
- \left( \nabla x_i + x_i \frac{F z_i}{RT} \nabla \phi \right )= \sum_{i=1 \atop j \neq i}^{\nu} \frac{1}{D_{ij} c}(x_j N_i - x_i N_j)
```
"""

# ╔═╡ d7e6e457-9d3f-4081-b2f7-060f85f18b04
md"""
## Grid
"""

# ╔═╡ a6073864-e1ac-4387-a1b4-9d6cb49dcf3d
md"""
Simple 1-D domain.
"""

# ╔═╡ c7a7fc8b-7364-48ea-ba1a-892c7ace96f5
function domain(data; nref=0)
	w=data.width
	hw=w/10.0*2.0^(-nref)
	W=collect(0:hw:w)	    
	simplexgrid(W)	
end

# ╔═╡ bd558b01-7085-4f61-92c2-c05c4e167b8f
md"""
## System
"""

# ╔═╡ a0d45d79-19fd-405b-bf83-e5f4f0a51cd7
md"""
The numerical fluxes $\textbf{J}$ are computed using the backslash operator (Julia solver for linear systems) 

$\textbf{J} 
= 
M
\ \backslash\ 
\textbf{F}$

"""

# ╔═╡ 1afba52b-5435-4514-a668-dbe9ca8b1eb4
function reaction(f,u,node,data)
	(;ni,z,iϕ,c_tot,Fa) = data
	# ∑xi = 1
	f[ni] = sum(u[1:ni]) - 1.0
	
	# electro neutrality
	sumq = zero(eltype(u))
	for i=1:ni
		#sumq += u[i]*cw*z[i]*Fa
		sumq += u[i]*z[i]
	end
	f[iϕ] = -sumq*Fa*c_tot
	#f[iϕ] = sumq
	#f[iϕ] = zero(eltype(u))
end

# ╔═╡ 287d66ec-bfea-4d87-9b73-d54727d85236
function bcond(f,u,bnode,data)
	(;ni,iϕ,cHp_bulk,cClm_bulk,c_tot) = data

	# left = bulk boundary condition
	boundary_dirichlet!(f,u,bnode,1,1, cHp_bulk/c_tot) # H+ bulk / left boundary
	boundary_dirichlet!(f,u,bnode,2,1, cClm_bulk/c_tot) # Cl- bulk / left boundary

		
	# right = interface of interest: eg. working electrode / ion exchange membrane
	boundary_dirichlet!(f,u,bnode,1,2,0.0) # H+ right
	#boundary_dirichlet!(f,u,bnode,2,2,0.0) # Cl- right

	boundary_dirichlet!(f,u,bnode,iϕ,1,0.0) # bulk potential ϕ
	#boundary_dirichlet!(f,u,bnode,iϕ,2,0.0) # potential ϕ	

	
	#boundary_dirichlet!(f,u,bnode,ispec,ireg,val)
end

# ╔═╡ 3d2b9534-ac73-45b9-a774-8c529af64862
md"""
## Fluxes Through Boundaries
"""

# ╔═╡ c1adc446-4304-4524-9010-9387e1d27d9e
md"""
## Model Data
"""

# ╔═╡ 9d2ec456-c931-48de-b268-84161be712c1
Base.@kwdef mutable struct ModelData

	ions::Vector{String} = ["H+", "Cl-", "w"]
	z::Vector{Float64} = [1.0, -1.0, 0.0]
	ni::Int64 = length(ions)
	iϕ::Int64 = ni+1

	embed_p::Float64 = 1.0e-2
	
	rhow::Float64 = 1000.0*ufac"kg/m^3"
	Mw::Float64 = 18.015*ufac"g/mol"
	cw::Float64 = rhow/Mw*ufac"mol/m^3"

	# assume: molal partial volumes are the same for all species
	# v_H+ = v_Cl- = v_w -> c_tot = const

	cHp_bulk::Float64=1.0*ufac"mol/dm^3"
	cClm_bulk::Float64=cHp_bulk

	c_tot::Float64 = cw+cHp_bulk+cClm_bulk

	Dij::Vector{Float64} = [Inf, 9.3, 2.0] .* 1e-9 *ufac"m^2/s"
	#Dij::Vector{Float64} = [9.3, 9.3, 2.0] .* 1e-9 *ufac"m^2/s"
	width::Float64 = 1.0*ufac"cm"
	# width::Float64 = 20.0*ufac"nm"
	
	T::Float64 = 298.15*ufac"K"
	p::Float64 = 1.0*ufac"bar"

	# Faradays constant
	const Fa::Float64 = ph"N_A*e"*ufac"C/mol"
end

# ╔═╡ 77798783-8eac-4e74-afb8-06f4d4e06a58
let
	vis=GridVisualizer()
	gridplot!(vis, domain(ModelData()), legend=:best, show=true)	
end

# ╔═╡ cf5b6bf6-c0af-47be-beba-d63586df4485
function D_matrix_ion(data)
	ni=data.ni
	D = zeros(eltype(data.Dij),ni,ni)
	cnt = 1
	for i=1:(ni-1)
		for j=(i+1):ni
			D[j,i] = data.Dij[cnt]
			cnt +=1
		end
	end
	Symmetric(D, :L)
end

# ╔═╡ a9870a11-d529-414d-a745-cc8c24f04e06
function M_matrix_ion(data,x)
	
	ni=data.ni
	D=D_matrix_ion(data)
	M=zeros(eltype(x), ni, ni)
	for i=1:ni
		for j=1:ni
			if j != i
				M[i,i] = x[j]/D[i,j]
				M[i,j] = -x[i]/D[i,j]
			end
		end	
	end
	#c = p/(ph"R"*T)
	M / data.c_tot	
end

# ╔═╡ c8881b90-b33b-4c2c-aff4-ca0377fda5b2
function flux(f,u,edge,data)
	(;ni,iϕ,z,cw,Fa,T) = data
	
	xbar = zeros(eltype(u), ni)
	F = zeros(eltype(u), ni)

	δϕ = u[iϕ,1]-u[iϕ,2]
	
	for i=1:ni
		xbar[i] = 0.5*(u[i,1]+u[i,2])

		Felec = Fa*z[i]/(ph"R"*T)*δϕ
		bp,bm = fbernoulli_pm(Felec)
		
		F[i] = (bm*u[i,1]-bp*u[i,2])
	end	
	
	# computation of fluxes J
	
	J = M_matrix_ion(data,xbar) \ F
			
	f[1:(ni-1)] = J

	#f[iϕ] = δϕ

end

# ╔═╡ ad43d343-7e02-4cd7-adcf-5e1f483418dc
function main(;data=ModelData())
	(;ni,iϕ,cHp_bulk,cClm_bulk,cw,c_tot) = data
	grid=domain(data, nref=5)

	 function pre(sol,par)
	 	ng=data.ng
	 	iT=data.iT

	 end
	
	sys=VoronoiFVM.System( 	grid;
							data=data,
							flux=flux,
							reaction=reaction,
							bcondition=bcond
							)

	
	enable_species!(sys; species=collect(1:(ni+1))) # aqueous e'lyte species + elec. pot
	inival=unknowns(sys)
	inival[1,:] .= cHp_bulk/c_tot
	inival[2,:] .= cClm_bulk/c_tot
	inival[3,:] .= cw/c_tot
	inival[iϕ,:] .= 0.0
	
	
	sol=solve(sys;inival)

	sol,sys,grid,data
end

# ╔═╡ be989d1d-fd95-4a1a-a6a4-429cbaac4939
sol,sys,grid,data = main();

# ╔═╡ 9ed7cecd-19d4-4621-9593-29689d2bfd23
let
	(;ni,ions,iϕ) = data
	
	c1 = colorant"red"
	c2 = colorant"blue"
	cols=range(c1, stop=c2, length=ni+1)
	
	vis=GridVisualizer()
	
	for i=1:ni
		scalarplot!(vis, grid, sol[i,:], color=cols[i], label="$(ions[i])", clear=false, legend=:best)
	end

	scalarplot!(vis, grid, sol[iϕ,:], color=cols[iϕ], label="ϕ", clear=false, legend=:best)

	reveal(vis)
end

# ╔═╡ c29d1156-7941-4091-b94e-94569786dbb0
let
	tf=VoronoiFVM.TestFunctionFactory(sys)
	Γ_where_T_equal_1=[2]
	Γ_where_T_equal_0=[1]
	T=testfunction(tf,Γ_where_T_equal_0,Γ_where_T_equal_1)
	I=integrate(sys,T,sol)
end

# ╔═╡ Cell order:
# ╠═23dc6930-c354-11ed-110a-235339eca423
# ╟─5ede4875-1dbf-4930-9538-39b10d4bd377
# ╠═1593e889-7b4f-40f9-b6e8-ef0564381c81
# ╟─711a3678-e3c0-4193-a6cf-0d525ddeff91
# ╟─b55b90f2-bace-485c-a745-3fca9871152f
# ╟─47fea77b-d605-40ac-b436-078fd91294fe
# ╟─4a33c6a3-ad59-4efc-997e-09ca780e5fd2
# ╟─2aefd6b8-f35c-419e-a2c4-83ad44f1be77
# ╟─e7177289-5c00-44f0-9ad3-92681d02c05d
# ╟─67a38c46-a952-473d-864e-53ca99b16174
# ╟─a4bb4cff-e933-45ac-8ad2-c2f69328f62d
# ╟─99faa887-0c70-4834-80e0-9f4e39b1b00a
# ╟─87c2c676-9bf0-4716-8fe4-e44a0f884412
# ╟─72514d88-d17b-4dfe-8921-623aa6e3b42f
# ╟─3ede7aa4-d924-4825-b0cf-f5f9b52c5e3b
# ╟─21d032ab-3de3-400a-b6fb-2c9805e4916b
# ╟─ac882611-643c-4c94-9acd-38de7566498e
# ╟─0dbfa843-7b38-41bf-aa4c-14d528bbf75f
# ╟─25ba905a-3342-4139-935d-ebcf0ad4be56
# ╟─10f5b322-e52b-45b7-ae60-a9f5d82b5f99
# ╟─e93d0059-0b6b-4911-a148-0c244ea010c4
# ╟─27b0be8c-e57f-4396-b132-b4ddd64a1892
# ╟─878dbe91-b305-4155-a98d-c4aced05e8f3
# ╟─973a9c4e-b7c3-4640-96aa-b8f09ebe0b28
# ╟─49d41a14-9f28-4130-bb6b-a03159336ef1
# ╟─22f39a64-bfb9-464e-bd2e-24fc72595749
# ╟─6d6d0af1-4185-43db-a5a6-2faad0592668
# ╟─d7e6e457-9d3f-4081-b2f7-060f85f18b04
# ╟─a6073864-e1ac-4387-a1b4-9d6cb49dcf3d
# ╠═c7a7fc8b-7364-48ea-ba1a-892c7ace96f5
# ╠═77798783-8eac-4e74-afb8-06f4d4e06a58
# ╠═bd558b01-7085-4f61-92c2-c05c4e167b8f
# ╠═c8881b90-b33b-4c2c-aff4-ca0377fda5b2
# ╟─a0d45d79-19fd-405b-bf83-e5f4f0a51cd7
# ╠═1afba52b-5435-4514-a668-dbe9ca8b1eb4
# ╠═287d66ec-bfea-4d87-9b73-d54727d85236
# ╠═ad43d343-7e02-4cd7-adcf-5e1f483418dc
# ╠═be989d1d-fd95-4a1a-a6a4-429cbaac4939
# ╠═9ed7cecd-19d4-4621-9593-29689d2bfd23
# ╟─3d2b9534-ac73-45b9-a774-8c529af64862
# ╠═c29d1156-7941-4091-b94e-94569786dbb0
# ╠═c1adc446-4304-4524-9010-9387e1d27d9e
# ╠═9d2ec456-c931-48de-b268-84161be712c1
# ╠═cf5b6bf6-c0af-47be-beba-d63586df4485
# ╠═a9870a11-d529-414d-a745-cc8c24f04e06
