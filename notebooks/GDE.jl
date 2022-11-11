### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# ╔═╡ beeb281b-c15b-4e9f-beae-00a2fd2d7104
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	
	using VoronoiFVM
	using ExtendableGrids
	using GridVisualize
	using LessUnitful
	using PlutoVista	
	using PlutoUI
	using Colors
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 99d7a0bf-6ce2-403f-b382-00e726cbd762
TableOfContents()

# ╔═╡ 7dde1a15-7bd1-4e3f-a181-68040a9f51f2
md"""
# Gas diffusion electrode
Pluto notebook aiming at recreating published models for gas diffusion electrodes (GDEs) using the VoronoiFVM.jl package.

In particular, the 1D, isothermal, steady state model of a GDE used in CO2 electrolysis is to be reimplemented, presented in:

[Weng, L. C., et al. (2018).](https://pubs.rsc.org/en/content/articlehtml/2018/cp/c8cp01319e)
(Open Access)

$(Resource("https://pubs.rsc.org/image/article/2018/CP/c8cp01319e/c8cp01319e-f1_hi-res.gif"))

"""

# ╔═╡ e1a50f4f-158d-4ea7-8436-c7b0e5ee0865
md"""
# Model data
"""

# ╔═╡ 3066fadb-ae87-4672-9fc4-73a90e3a1e86
@Base.kwdef mutable struct ModelData
	# number of gas phase species
	ng::Int64		 		= 3
	
	# 1:CO, 2:H2, 3:CO
	gn::Dict{Int, String} 	= Dict(1=>"CO", 2=>"H2", 3=>"CO2")
	
	# stochiometric coeff for CT reaction
	nug::Vector{Int} 	= [1,1,-2]

	# number of transferred e- per CT reaction turnover
	nCT::Vector{Int} 		= [2]
	
	# gas density
	ρg::Float64 			= 1*ufac"kg/(m^3)"

	# gas velocity, convective transport, different for the 2 regions
	vg::Vector{Float64}		= 0*[1.7e-5,0.35]*ufac"m/s"
	#vg::Vector{Float64}		= -1e-5*[1.7e-5,0.35]*ufac"m/s"

	# gas phase diffusivity coefficients
	Dg::Vector{Float64} 	= fill(1.0*ufac"cm^2/s",ng)
	
	# gas phase species molar mass
	Mg::Vector{Float64} 	= fill(10.0*ufac"g/mol",ng)

	# number of aqueous phase species
	ni::Int64 				= 2

	# ion charge numbers
	z::Vector{Int}   		= [-1,1]
	
	# reference voltage
	E_ref::Float64   		= 0.0*ufac"V"

	# solid phase electrical conductivity
	σs::Float64 	= 200.0*ufac"S/m"

	# volume specific active surface area
	#av::Float64 	= 1.0*ufac"1/m"
	av::Float64 	= 1.5e7*ufac"1/m"
	
	e::Float64 = ph"e"
	F::Float64 = ph"e"*ph"N_A"

end

# ╔═╡ a9fa57c9-2ff4-4086-9a11-ab506c4531d6
md"""
### Discretization grid
- region 1: catalyst layer (__CL__), phases: __gas__ and __liquid__ and __solid__ phase, species: gas and ionic species

- region 2: diffusion medium (__DM__), phases: __gas__ and __solid__ phase, species: gas phase speceis
$(Resource("https://pubs.rsc.org/image/article/2018/CP/c8cp01319e/c8cp01319e-f3_hi-res.gif"))
"""

# ╔═╡ 49d6cfe9-2aaa-4f53-a34b-7e5eb72a001f
function grid2(;L_CL=5.0*ufac"μm",L_DM=325.0*ufac"μm",nref=0)
	X_CL=range(0,L_CL,length=10*2^nref+1)
	X_DM=range(L_CL,L_DM,length=10*2^nref+1)
	grid=simplexgrid(glue(X_CL,X_DM))

	cellmask!(grid,[L_CL],[L_DM],2)
    #bfacemask!(grid,[L_CL], [L_CL],2)
	#bfacemask!(grid,[L_DM], [L_DM],3)
    grid
end

# ╔═╡ 9b18f479-3141-4e84-a465-8f5542c9dbe2
md"""
### Boundary regions
Boundary region numbers are assigned to the interfaces between neighboring regions in the electrolyser cell.

The regions electrolyte channel and the gas channel are not part of this model but will be incorporated via boundary conditions:
1) electrolyte channel / catalyst layer: $\Gamma_{\text{EL,CL}} = 1$
1) diffusion medium / gas channel: $\Gamma_{\text{DM,GC}} = 2$
"""
#1) catalyst layer / diffusion medium: $\Gamma_{\text{CL,DM}} = 3$

# ╔═╡ 74444ac9-748a-48d4-b693-c16f53ae188f
begin
	const Γ_EL_CL=1
	const Γ_DM_GC=2
	#const Γ_CL_DM=3

	# bulk regions
	# CL = 1
	# DM = 2
end;

# ╔═╡ 41dffd5d-e4e5-4601-9e95-dc7598509f3f
gridplot(grid2(L_CL=50.0*ufac"μm"),resolution=(600,200),legend=:rt,show=true)

# ╔═╡ 69992e99-7485-4af1-8594-a254e541efe4
md"""
# System of Equations
"""

# ╔═╡ b83fc5d5-083d-401a-a78b-0c136297f752
md"""

## Gas phase species transport

Convective-diffusive gas species transport takes place in the diffusion medium (DM) and catalyst layer (CL). The variable describing the gas species distribution is the mass fraction $\omega_i$. (as in Weng, L. C., et al. 2018)
```math
\begin{aligned}
\nabla\cdot \vec N_i &= 0 \\
\vec N_i &= - \rho_{\text{g}} \left( D^{\text{eff}}_i \vec\nabla\omega_i- \omega_i \left[\vec u+\frac{D^{\text{eff}}_i}{M_{\text{n}}} \vec \nabla M_{\text{n}} \right] \right)
\end{aligned}
```

In the simplified case consider charge transfer reactions as a boundary condition.
```math
R_{\text{CT},i} = -M_i \sum_k \nu_{i,k} \frac{I_k}{n_k F}
```

"""

# ╔═╡ b03473a8-bda9-421d-a318-5073ab8a8a34
md"""

## Electron transport (solid phase $\phi_{\text{s}}$)

Equations from __Weng, L. C., et al. (2018).__

### Diffusion Medium
For the region DM (diffusion medium) where solid and gas phases are present, the solid phase potential $\phi_{\text{s}}$ is described via charge conservation and Ohm's law. No charge-transfer reactions occur in the DM.

```math
\begin{aligned}
\nabla\cdot \vec i_{\text{s}} &= 0 \\
\vec i_{\text{s}} &= - \sigma_{\text{s}} \vec \nabla\phi \\
\end{aligned}
```

### Catalyst Layer
For the region CL (catalyst layer) where (for now, simplified) solid and liquid phases are present, the solid phase potential $\phi_{\text{s}}$ is described via Ohm's law.

In the CL, $k$ charge-transfer reactions occur (measured by currentat density $i_k$), there is an reaction term, that causes the transport of charges (electrons) and therefore gradients in  $\phi_{\text{s}}$.

```math

- \nabla \cdot \left ( \sigma_{\text{s,CL}} \vec \nabla \phi_{\text{s}} \right) = - a_{\text{v}} \sum_k i_k

```


"""

# ╔═╡ 01da1e64-e500-4d6b-9102-9de56aa874a9
md"""
## Electrolyte potential (aqueous/liquid phase $\phi_{\text{l}}$)

### Catalyst layer
Assuming the DM is dry, the aqueous/liquid phase is only present in the CL.

Enforcing electro-neutrality, the Poisson equation for the potential in the aqueous electrolyte takes the following form:

```math
\nabla\cdot \left( \epsilon \vec \nabla \phi_{\text{l}} \right) = F \rho_{\text{l}} \sum_j \frac{\omega_j z_j}{M_j} = 0
```
where $F$ is Faraday's constant and $\rho_{\text{l}}$ is the liquid phase density. 

This simplifies to the implicit equation for $\phi_{\text{l}}$:
```math
\sum_j \frac{\omega_j z_j}{M_j} = 0
```
"""

# ╔═╡ 021fa393-10aa-4339-bf90-6e818e46b299
md"""
## Ionic species transport

### Catalyst layer
Ionic species only exist in the aqueous/liquid phase which is only present in the CL.
The transport is described by the Nernst-Planck flux, also commonly called drift-diffusion flux.

The variable describing the ionic species distribution is the mass fraction $\omega_j$. (as in Weng, L. C., et al. 2018)

__Volumetric__ charge transfer reactions $R_{\text{CT},j}$ occur in the CL.

```math
\begin{aligned}
\nabla\cdot \vec N_j &= - R_{\text{CT},j} \\

\vec N_j &= - \rho_{\text{l}} \left( D^{\text{eff}}_j \vec\nabla\omega_j - \omega_j  z_j u_j F \vec \nabla \phi_{\text{l}}  \right) \\

R_{\text{CT},j} &= M_j a_{\text{v}} \sum_k \nu_{j,k} \frac{I_k}{n_k F}

\end{aligned}
```
"""

# ╔═╡ b522b245-be17-444c-92eb-bc1c64e59706
md"""

## Charge transfer reaction
In the simplified case only consider the CO evolution reaction:
```math
CO_2 + H_2O + 2 e^-\rightarrow CO + 2 OH^-
```
"""

# ╔═╡ 393f5b04-9368-4fb0-a64a-de3a41410fcf
begin
	const ngas=ModelData().ng
	const nion=ModelData().ni
	
	# variable indices	
	const igas=1
	const iion=igas+ngas
	const iϕs=iion+nion
	const iϕl=iϕs+1
end;

# ╔═╡ d857b437-e639-45bb-a64c-b3153eb12f53
# charge transfer reaction rate
function R_CT(u, data::ModelData)
 
	# current density of CT reaction, mA cm-2, reaction rate limited by availability of reactant CO2
	rr=500*ufac"mA/cm^2"*u[3]
	
	
	# species mass flux from CT reaction, eq. (26)
	@. data.Mg*data.nug*data.av*rr/data.nCT/data.F
	
end

# ╔═╡ fe3e0095-31b8-4921-ab1a-da4cb42b9302
function flux(f,u,edge,data)
		
	# gas phase species flux eqs. (7,8)
	ngas=data.ng
	nion=data.ni
	
	
	sumw_M_k=0.0
	sumw_M_l=0.0
	for i=1:ngas
		sumw_M_k+=u[i,1]/data.Mg[i]
		sumw_M_l+=u[i,2]/data.Mg[i]
	end
	Mn_k = 1/sumw_M_k
	Mn_l = 1/sumw_M_l
	#δM = 2/(Mn_k+Mn_l)*(Mn_k-Mn_l)
	δM = log(Mn_k)-log(Mn_l)

	for i=1:(ngas-1)
	
		#bp, bm = fbernoulli_pm(data.vg/data.Dg[i])
		#bp, bm = fbernoulli_pm(data.vg/data.Dg[i] + δM)
		bp, bm = fbernoulli_pm(data.vg[edge.region]/data.Dg[i] + δM)
		
		f[i]=data.Dg[i]*(bm*u[i,1]-bp*u[i,2])
		#f[i]=data.Dg[i]*(u[i,1]-u[i,2])
	end
	
	# aqueous phase species flux
	for i=iion:(iion+nion-1)
		f[i]=u[i,1]-u[i,2]
	end
	
	# solid phase potential / DM + CL

	f[iϕs]=data.σs*(u[iϕs,1]-u[iϕs,2])

	# aqueous phase potential / CL
	f[iϕl]=u[iϕl,1]-u[iϕl,2]


end

# ╔═╡ 97879adc-8493-4b51-9071-82c7d5151d01
function bcond(f,u,bnode,data)
	
	boundary_dirichlet!(f,u,bnode,1,Γ_DM_GC,0.01) # CO
	boundary_dirichlet!(f,u,bnode,2,Γ_DM_GC,0.01) # H2

	for i=iion:(iion+nion-1) # define bc for all ion species
		boundary_dirichlet!(f,u,bnode,i,Γ_EL_CL,1)
	end
	boundary_dirichlet!(f,u,bnode,iϕs,Γ_DM_GC,1) # spec, boundary(left), value
	boundary_dirichlet!(f,u,bnode,iϕl,Γ_EL_CL,1)
	

end

# ╔═╡ 5d7eb635-14cb-46df-9b6c-4c402596889e
function reaction(f,u,node,data)

	# volumetric charge-transfer reaction in CL region
	if node.region == 1
		rr=R_CT(u,data)
		for i=1:(ngas-1)
			f[i]=-rr[i]
		end
	end
	
	# mole/mass fractions sum to 1 for last species nspec
	#f[nspec]=sum(u[1:nspec])-1
	
	f[ngas]=log(sum(u[1:ngas]))
end

# ╔═╡ 428a8c4a-9786-44a7-89be-9c2954f8a5d9
grid=grid2(L_CL=50*ufac"μm",nref=0)

# ╔═╡ 1e957121-f698-481b-86dc-b61bd1b7f874
begin
	# create system + solution
	sys=VoronoiFVM.System( 	grid;
							data=ModelData(),
							flux=flux,
							reaction=reaction,
							bcondition=bcond
							)
	enable_species!(sys; species=collect(igas:(igas+ngas-1)), regions=[1,2])
	enable_species!(sys; species=collect(iion:(iion+nion-1)), regions=[1])
	enable_species!(sys; species=[iϕs], regions=[1,2])
	enable_species!(sys; species=[iϕl], regions=[1])

	
	inival=unknowns(sys)
	inival[:,:].=1
	inival[igas:(igas+ngas-1),:].=1/ngas
	
	sol=solve(inival,sys)
end;

# ╔═╡ 3f938444-c12d-4853-9eab-6658924fed1f
let
	# CO, H2 outflow over Γ_DM_GC boundary
	tf=VoronoiFVM.TestFunctionFactory(sys)

	Γ_where_T_equal_1=[Γ_DM_GC]
	Γ_where_T_equal_0=[Γ_EL_CL]
	T=testfunction(tf,Γ_where_T_equal_0,Γ_where_T_equal_1)
	
	I=integrate(sys,T,sol)
end

# ╔═╡ 25bdfa94-e0b4-408e-b690-decf383eba2e
R=integrate(sys,reaction,sol)

# ╔═╡ b5767de1-bb9c-4ef3-a0f3-78bc2b1fbd07
let
	# visualization
	vis=GridVisualizer(layout=(2,1), resolution=(600,600))
	
	c1 = colorant"blue"
	c2 = colorant"red"
	cols=range(c1, stop=c2, length=ngas)

	for i in 1:ngas
		scalarplot!(vis[1,1], grid, sol[i,:], color=cols[i], label="$(ModelData().gn[i])", clear=false, legend=:ct)
	end

	for i in iion:(iion+nion-1)
		scalarplot!(vis[2,1], grid, sol[i,:], color=cols[i-iion+1], label="ion$(i)", clear=false, legend=:ct)
	end
	
	#scalarplot!(vis[2,1], grid, sol[iϕs,:], color=cols[1], label="ϕs", clear=false, legend=:ct)
	reveal(vis)

end

# ╔═╡ Cell order:
# ╟─99d7a0bf-6ce2-403f-b382-00e726cbd762
# ╠═beeb281b-c15b-4e9f-beae-00a2fd2d7104
# ╟─7dde1a15-7bd1-4e3f-a181-68040a9f51f2
# ╟─e1a50f4f-158d-4ea7-8436-c7b0e5ee0865
# ╠═3066fadb-ae87-4672-9fc4-73a90e3a1e86
# ╟─a9fa57c9-2ff4-4086-9a11-ab506c4531d6
# ╠═49d6cfe9-2aaa-4f53-a34b-7e5eb72a001f
# ╟─9b18f479-3141-4e84-a465-8f5542c9dbe2
# ╠═74444ac9-748a-48d4-b693-c16f53ae188f
# ╠═41dffd5d-e4e5-4601-9e95-dc7598509f3f
# ╟─69992e99-7485-4af1-8594-a254e541efe4
# ╟─b83fc5d5-083d-401a-a78b-0c136297f752
# ╟─b03473a8-bda9-421d-a318-5073ab8a8a34
# ╟─01da1e64-e500-4d6b-9102-9de56aa874a9
# ╟─021fa393-10aa-4339-bf90-6e818e46b299
# ╟─b522b245-be17-444c-92eb-bc1c64e59706
# ╠═393f5b04-9368-4fb0-a64a-de3a41410fcf
# ╠═d857b437-e639-45bb-a64c-b3153eb12f53
# ╠═fe3e0095-31b8-4921-ab1a-da4cb42b9302
# ╠═97879adc-8493-4b51-9071-82c7d5151d01
# ╠═5d7eb635-14cb-46df-9b6c-4c402596889e
# ╠═1e957121-f698-481b-86dc-b61bd1b7f874
# ╠═3f938444-c12d-4853-9eab-6658924fed1f
# ╠═25bdfa94-e0b4-408e-b690-decf383eba2e
# ╠═428a8c4a-9786-44a7-89be-9c2954f8a5d9
# ╠═b5767de1-bb9c-4ef3-a0f3-78bc2b1fbd07
