### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 4eebb4d0-4ee9-11ed-0892-05b637938ba7
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

# ╔═╡ 76b3f5de-efb6-4837-bfa6-c26f7f3fb215
TableOfContents()

# ╔═╡ bd3af7ca-e2fe-4fce-bd0c-88bc51ae046b
abstract type AbstractElectrolyteData end

# ╔═╡ b48f424e-a9b7-434c-b5b0-5f953de68664
md"""
# The model
"""

# ╔═╡ dc1d238c-a04a-477f-a8fb-5b13c6684ff2
md"""
The envisioned model covers the regions diffusion medium (DM), catalyst layer (CL) and membrane (Mem) with phases and species as described in the table.
"""

# ╔═╡ 4bf59f98-0f54-499d-a081-41f456fa5719
md"""
$\begin{array}{c|c|p{}|c}
\hline
\text{Region} & \text{Phases} & \text{Species} & \text{Variables}\\
\hline
\text{DM} & \text{s, g, l} & \text{CO}^{\text{(g)}}, \text{CO}^\text{(g)}_{2}, \text{H}_{2}\text{O}^\text{(g)}, \text{H}_{2}\text{O}^\text{(l)}, \text{e}^- & p_{\text{CO}}, p_{\text{CO}_2}, p_{\text{H}_2\text{O}}, p_{\text{w,L}}, \phi_{\text{s}}\\

\text{CL} & \text{s, g, l, ion} & \text{CO}^{\text{(g)}}, \text{CO}^\text{(g)}_{2}, \text{H}_{2}\text{O}^\text{(g)}, \text{H}_{2}\text{O}^\text{(l)}, \text{e}^-, \text{H}_{2}\text{O}^\text{(aq)}, \text{CO}_{2}^\text{(aq)}, \text{OH}^{-\text{(aq)}} & p_{\text{CO}}, p_{\text{CO}_2}, p_{\text{H}_2\text{O}}, p_{\text{w,L}}, \phi_{\text{s}}, c_{\text{H}_2\text{O}}, c_{\text{CO}_2}, c_{\text{OH}^-}, \phi_{\text{ion}}\\

\text { Mem } & \text{l, ion} & \text{H}_{2}\text{O}^\text{(l)}, \text{H}_{2}\text{O}^\text{(aq)}, \text{CO}_{2}^\text{(aq)}, \text{OH}^{-\text{(aq)}} & p_{\text{w,L}}, c_{\text{H}_2\text{O}}, c_{\text{CO}_2}, c_{\text{OH}^-}, \phi_{\text{ion}}\\
\hline
\end{array}$
"""

# ╔═╡ 01a76a61-acc2-463d-8a5a-b8e22392a50f
md"""
Charged particles in a bounded domain ``\Omega``:
- ``N`` Particle densities (concentration ``c_i(x)``)
- Charge numbers ``z_i``
- Electrolyte potential resulting from densities of charged particles and enforcement of charge neutrality ``\phi``
"""

# ╔═╡ 8abace6c-f126-4de3-b92d-484fa31811b9
md"""
## Electrostatics
"""

# ╔═╡ bbb39bea-26a9-45b0-8528-fecd6636a28e
md"""
The charge neutrality constraint is used to calculate the electrolyte potential:

`` 0 = \sum_{i=1}^N z_i c_i``


"""

# ╔═╡ 1c25bd96-09dc-4593-8616-f6ad2b8c6746
md"""
### Continuity equation

For the particle densities we assume a continuity equation  involving the particle fluxes $N_i$:

``\partial_t c_i + \nabla\cdot N_i = R_i``
"""

# ╔═╡ b674382d-1c6c-4b8e-891c-2db213ed2bcd
md"""
### Which forces make the particles move ?
- concentration gradients  $\nabla c_i$
- electric field  $\nabla \phi$

So we can formulate in the simplest case:

``N_i = - D_i (\nabla  c_i + \frac{e}{k_BT} z_i c_i \nabla \phi``)
"""

# ╔═╡ 0ba193f4-fb72-4a4c-9f35-f46e948df2d7
md"""
We derived the simplest case of the Poisson-Nernst-Plank system (as electrochemiststs say) or the van Roosbroeck system (as semiconductor people say).
"""

# ╔═╡ 1ecc9a61-5c98-46a1-ba03-14095b33b20c
md"""
In the sequel, for simplicity, we set all constants $\varepsilon, e, D_i, k_B, T=1$


We arrive at the system

``
\begin{aligned}
0 &= \sum z_i c_i\newline
\partial_t c_i  - \nabla(\nabla c_i - z_i c_i \nabla \phi)  &= R_i
\end{aligned}
``

This system is nonlinear, convection dominant in places where ``\nabla\phi`` is large. We want the concentrations to be nonnegative and to be consistent to the thermodynamic properties behind this problem (not explained here).
"""

# ╔═╡ 932958ba-a7a2-4100-addb-851816b45b24
md"""
### Creating a discretization grid

For this purpose we use the package ExtendableGrids which allows to manage
simplex grid in 1/2/3D.
"""

# ╔═╡ 2bd278e7-898d-46c3-b9eb-5d4edcb060bb
function grid2(;nref=0,L_CL=5.0*ufac"nm")
	hmin=1.0e-1*ufac"nm"*2.0^(-nref)
    hmax=1.0*ufac"nm"*2.0^(-nref)
    L=20.0*ufac"nm"
    X=geomspace(0,L,hmin,hmax)
    grid=simplexgrid(X)
	
	cellmask!(grid,[0],[L_CL],2)
	grid
end;

# ╔═╡ 0dd23774-4f69-44af-923d-0e7eaefe6bda
begin
	grid=grid2()
 	gridplot(grid,resolution=(600,200),legend=:rt,show=true)
end

# ╔═╡ 032f7e2b-9bd2-4ca0-ba41-edadb2ada2ee
function charge(u,data::AbstractElectrolyteData)
    q=zero(eltype(u))
    for ic=1:data.nc
        q+=u[ic] * data.z[ic]
    end
    q*data.F
end

# ╔═╡ 0c414c82-6801-4b8d-8684-e1ab756c30a5
function sflux(ic,dϕ,ck,cl,data)
    (;D,z,F,RT) = data
    #bp, bm = fbernoulli_pm(z[ic] * dϕ*F/RT  + dμex(βk,βl,electrolyte)/RT)
	# simplified case for β=1 -> dμex = 0
	bp, bm = fbernoulli_pm(z[ic] * dϕ*F/RT)
    D[ic] * (bm*ck - bp*cl)
end

# ╔═╡ eb5df477-a000-46fc-8518-47d8147df2c2
md"""
# Physics functions
"""

# ╔═╡ a78ebcf6-eab5-4b54-bb7e-3f6bee5b0795
md"""
PNP: Poisson-Nernst-Planck
"""

# ╔═╡ 7e0fe75c-b3f9-4cbc-9797-1b2253c31a2b
function pnpflux(f,u,edge,data)
	(;iϕ, v0, v, nc) = data
	dϕ = u[iϕ,1]-u[iϕ,2]

	for ic = 1:nc        
        ## Regularize ck,cl so they don't become zero
        ck,cl=u[ic,1]+data.epsreg,u[ic,2]+data.epsreg

    	#excess chemical potential (SEDAN) scheme
		#f[ic]=sflux(ic,dϕ,ck,cl,βk,βl,bar_ck, bar_cl,electrolyte)
		f[ic]=sflux(ic,dϕ,ck,cl,data)
    end
end

# ╔═╡ f097c058-977f-48cb-bad9-0a04319379f3
md"""
Bulk electrolyte boundary condition
"""

# ╔═╡ fd07e120-3312-462c-8efa-500f602584fa
function bulkbcondition(f,u,bnode,data;region=data.Γ_bulk)
    (;iϕ,nc,ϕ_bulk,c_bulk) = data
    if bnode.region==region
        boundary_dirichlet!(f,u,bnode;species=iϕ,region,value=ϕ_bulk)
        #boundary_dirichlet!(f,u,bnode;species=ip,region,value=p_bulk)
        for ic=1:nc
            boundary_dirichlet!(f,u,bnode;species=ic,region,value=data.c_bulk[ic])
        end
    end
end

# ╔═╡ d784e1d2-17fa-428f-8963-dfd0ddb60461
function grid1(;nref=0)
	hmin=1.0e-1*ufac"nm"*2.0^(-nref)
    hmax=1.0*ufac"nm"*2.0^(-nref)
    L=20.0*ufac"nm"
    X=geomspace(0,L,hmin,hmax)
    grid=simplexgrid(X)
end

# ╔═╡ f0d7f04a-0b9f-4be7-bb42-8b97ce0a8676
begin
	const ihplus=1
    const ife2 = 2
    const ife3 = 3
    const iso4 = 4
	const F=ph"N_A*e"
end;

# ╔═╡ 09847b05-b5ef-41fd-ae23-1762c281d69a
function pnpflux_mem(f,u,edge,data)
	(;iϕ, v0, v, nc) = data
	dϕ = u[iϕ,1]-u[iϕ,2]

	for ic = 1:nc
		if ic == iso4
			f[ic] = 0 #fixed species in polymere matrix
		else
	        ## Regularize ck,cl so they don't become zero
	        ck,cl=u[ic,1]+data.epsreg,u[ic,2]+data.epsreg
	
	    	#excess chemical potential (SEDAN) scheme
			#f[ic]=sflux(ic,dϕ,ck,cl,βk,βl,bar_ck, bar_cl,electrolyte)
			f[ic]=sflux(ic,dϕ,ck,cl,data)
		end
    end
end

# ╔═╡ d9b34916-597f-463a-86b9-f9b2b1856c05
function pnpunknowns(sys)
    (; iϕ,nc,c_bulk) = sys.physics.data
    u=unknowns(sys)
    @views u[iϕ,:] .= 0
    #@views u[ip,:] .= 0
    for ic=1:nc
        @views u[ic,:] .= c_bulk[ic]
    end
    u
end

# ╔═╡ 5a568ec5-9770-45d3-b95c-db7ba6605110
function ivsweep(sys;potentials=[0.0:0.005:1.0...]*ufac"V",ispec=1,solver_kwargs...)
		
	sols = []
	vs = zeros(0)
	is = zeros(0)
	
	data=sys.physics.data
	data.ϕ_we=0
    control=SolverControl(;solver_kwargs...)
    iϕ=data.iϕ

	inival0 = solve(sys;inival=pnpunknowns(sys))
    inival=copy(inival0)
    sol=copy(inival0)

	for ϕ in potentials
		data.ϕ_we=ϕ
        solve!(sol, inival, sys; control)

		inival .= sol
		I=-integrate(sys,sys.physics.breaction,sol; boundary=true)[:,data.Γ_we]

		push!(vs, ϕ)
		push!(is, I[ispec]*F)
		push!(sols, copy(sol))
	end
	vs,is,sols
end

# ╔═╡ f9178a7d-2b2b-4ad1-9ee3-941f4e2f27a6
md"""
# Auxiliary
"""

# ╔═╡ 4aaaf629-b265-45b2-b449-99379c3e3cd5
rlog(x,data::AbstractElectrolyteData)=rlog(x,eps=data.epsreg);

# ╔═╡ 75392e35-c597-47b2-8698-cc7951529fa8
function rlog(x;eps=1.0e-20)
    if x<eps
        return log(eps)+(x-eps)/eps
    else
        return log(x)
    end
end;

# ╔═╡ 5dae223d-951a-4746-b5af-27c7a28a42a0
function rexp(x;trunc=500.0)
    if x<-trunc
        1.0/rexp(-x;trunc)
    elseif x<=trunc
        exp(x)
    else
        exp(trunc)*(x-trunc+1)
    end
end;

# ╔═╡ 8e2b1dd6-af4d-4b83-987b-dda484e01b0c
rrate(R0,β,A)=R0*(rexp(-β*A) - rexp((1-β)*A));

# ╔═╡ bc02abf0-f64d-4248-8c7e-5c8ac7de9358
md"""
Calculate chemical potential of species with concentration c
```math
        μ = v(p-p_{ref}) + RT\log \frac{c}{\bar c}
```
"""

# ╔═╡ 44e785f4-74fe-4fb0-a0fe-5ee205c5bcbe
chemical_potential(c, barc, p, v, data)=rlog(c/barc,data)*data.RT+v*data.pscale*(p-data.p_bulk);

# ╔═╡ 903d9fe7-ddbd-4b42-8816-d216e46b6bbf
md"""
Calculate solvent concentration ``c_0`` and summary concentration ``\bar c`` from vector of concentrations `c`
using the incompressibility constraint (assuming ``κ_0=0``):
```math
 \sum_{i=0}^N c_i (v_i + κ_iv_0) =1
```

This gives

```math
 c_0v_0=1-\sum_{i=1}^N c_i (v_i+ κ_iv_0)
```

```math
c_0= 1/v_0 - \sum_{i=1}^N c_i(\frac{v_i}{v_0}+κ)
```

Then we can calculate 
```math
 \bar c= \sum_{i=0}^N c_i
```

"""

# ╔═╡ aa30627a-8359-40d4-8c20-34ea0f7ab775
md"""
Calculate relative (wrt. solvent) molar volume of i-th species ``v_{i,rel}=κ_i+\frac{v_i}{v_0}``.
"""

# ╔═╡ 9dd8e4a5-391e-4414-ad47-5247c2803101
vrel(ic,data)=data.v[ic]/data.v0+data.κ[ic]

# ╔═╡ 2fa560d3-c9b6-4ae6-9c99-cbb6ff8a98ec
function c0_barc(c, data)
    c0 = one(eltype(c)) / data.v0
    barc = zero(eltype(c))
    for ic = 1:data.nc
        barc += c[ic]
        c0 -= c[ic] * vrel(ic,data)
    end
    barc += c0
    c0+data.epsreg, barc+data.epsreg
end;

# ╔═╡ 58498779-d01f-451b-8a65-634baa2975ae
function pnpreaction(f,u,node,data)
	(;iϕ, nc,ϕ_we,ip,iϕ,R0,Δg,β,p_bulk,v,κ,v0,RT)=data
    ## Charge density
    f[iϕ] = -charge(u,data)

	# catalyst layer, CL=2
	if node.region==2
		c0,barc=c0_barc(u,data)
		μfe2=chemical_potential(u[ife2], barc, p_bulk, v[ife2]+κ[ife2]*v0, data)
		μfe3=chemical_potential(u[ife3], barc, p_bulk, v[ife3]+κ[ife3]*v0, data)
		r=rrate(R0,β,(μfe2 - μfe3 + Δg - F*(u[iϕ]-ϕ_we))/RT)
		f[ife2]-=r
		f[ife3]+=r
	end
end

# ╔═╡ 744cae5e-ec70-4ae3-9283-1c1d96df729f
function pnpreaction_mem(f,u,node,data)
	(;iϕ, nc,ϕ_we,ip,iϕ,R0,Δg,β,p_bulk,v,κ,v0,RT)=data
    ## Charge density
    f[iϕ] = -charge(u,data)
	f[iso4] = u[iso4] - data.c_bulk[iso4]

	# catalyst layer, CL=2
	if node.region==2
		c0,barc=c0_barc(u,data)
		μfe2=chemical_potential(u[ife2], barc, p_bulk, v[ife2]+κ[ife2]*v0, data)
		μfe3=chemical_potential(u[ife3], barc, p_bulk, v[ife3]+κ[ife3]*v0, data)
		r=rrate(R0,β,(μfe2 - μfe3 + Δg - F*(u[iϕ]-ϕ_we))/RT)
		f[ife2]-=r
		f[ife3]+=r
	end
end

# ╔═╡ 20f230b0-bc12-4c38-a446-ce7d8607896f
Base.@kwdef mutable struct ElectrolyteData <: AbstractElectrolyteData
    "Number of ionic species."
    nc::Int=2
    
    "Potential index in species list."
    iϕ::Int=nc+1
    
    "Pressure index in species list"
    ip::Int=nc+2
    
    "Mobility coefficient"
    D::Vector{Float64}=fill(2.0e-9*ufac"m^2/s",nc) 

    "Charge numbers of ions"
    z::Vector{Int}=[ (-1)^(i-1) for i=1:nc]
    
    "Molar weight of solvent"
    M0::Float64=18.0153*ufac"g/mol"
    
    "Molar weight of ions"
    M::Vector{Float64}=fill(M0,nc)

    "Molar volume of solvent"
    v0::Float64=1/(55.4*ufac"M")

    "Molar volumes of ions"
    #v::Vector{Float64}=fill(v0,nc)
	v::Vector{Float64}=fill(0.0,nc)
	
	"Charge transfer reaction rate constant"
	R0::Float64=1.0e-6*ufac"mol/(cm^2*s)"

	"Change in gibbs free energy of charge transfer reaction"
	Δg::Float64=0.0*ufac"J/mol"

	"Charge transfer coefficient of Butler-Volmer equation"
    β::Float64=0.5
	
    "Solvation numbers"
    κ::Vector{Float64}=fill(10.0,nc)
    
    "Bulk ion concentrations"
    c_bulk::Vector{Float64}=fill(0.1*ufac"M",nc)
    
    "Bulk voltage"
    ϕ_bulk::Float64=0.0*ufac"V"
    
    "Bulk pressure"
    p_bulk::Float64=0.0*ufac"Pa"

    "Bulk boundary number"
    Γ_bulk::Int=2

    "Working electrode voltage"
    ϕ_we::Float64=0.0*ufac"V"
    
    "Working electrode  boundary number"
    Γ_we::Int=1
    
    "Temperature"
    T::Float64=(273.15+25)*ufac"K"
    
    "Molar gas constant scaled with temperature"
    RT::Float64=ph"R"*T
    
    "Faraday constant"
    F::Float64=ph"N_A*e"
    
    "Dielectric permittivity of solvent"
    ε::Float64=78.49
    
    "Dielectric permittivity of vacuum"
    ε_0::Float64=ph"ε_0"
    
    "Pressure scaling factor"
    pscale::Float64=1.0e9
    
    "Local electroneutrality switch"
    #eneutral::Bool=true
    
    """
    [Flux caculation scheme](@id fluxes)
    This allows to choose between
    - `:μex` (default): excess chemical potential (SEDAN) scheme, see [`sflux`](@ref)
    - `:act` : scheme based on reciprocal activity coefficients, see [`aflux`](@ref)
    - `:cent` : central scheme, see [`cflux`](@ref).
    """
    #scheme::Symbol=:μex
    
    """
    Regularization parameter used in [`rlog`](@ref)
    """
    epsreg::Float64=1.0e-20
end

# ╔═╡ a02a4def-e2dc-4c9e-98b8-75524963e242
function PNPSystem(grid,bcondition;data=ElectrolyteData(),kwargs...)
    sys=VoronoiFVM.System(grid;
                          data=data,
                          flux=pnpflux,
                          reaction=pnpreaction,
                          #storage=pnpstorage,
                          bcondition,
                          species=[1:data.nc..., data.iϕ],
						  regions=[1,2],
                          kwargs...
                          )
end

# ╔═╡ 2299256d-ba0c-45b3-a331-1f75fcc02c45
function Setup(;nref=0,R0=1.0e-2,κ=10.0,kwargs...)
	@local_phconstants N_A e R ε_0
    F=N_A*e
    @local_unitfactors cm μF mol dm s mA A nm

	defaults=(; max_round=3,
              tol_round=1.0e-9,
              verbose=false,
              tol_relative=1.0e-8,
              tol_mono=1.0e-10)
	
	kwargs=merge(defaults, kwargs)


	#grid=grid1()
	grid=grid2()


    R0=R0*ufac"mol/(cm^2*s)"

    ihplus=1
    ife2 = 2
    ife3 = 3
    iso4 = 4
	
	#Working electrode boundary condition
	function halfcellbc(f,u,bnode,data)
		(;nc,Γ_we,Γ_bulk,ϕ_we,ip,iϕ,v,v0,RT)=data
		bulkbcondition(f,u,bnode,data;region=Γ_bulk)
		#if bnode.region==Γ_we
		#	Δg = data.Δg
    	#	β = data.β
		#	c0,barc=c0_barc(u,data)
		#	μfe2=chemical_potential(u[ife2], barc, u[ip], v[ife2]+κ*v0, data)
		#	μfe3=chemical_potential(u[ife3], barc, u[ip], v[ife2]+κ*v0, data)
		#	r=rrate(R0,β,(μfe2 - μfe3 + Δg - F*(u[iϕ]-ϕ_we))/RT)
		#	f[ife2]-=r
		#	f[ife3]+=r
		#end
		nothing
	end

    data=ElectrolyteData(;nc=4,
						 z=[1,2,3,-2],
						 κ=fill(κ,4),
                         Γ_we=1,
                         Γ_bulk=2,
						 R0=R0,
                         )


	data.c_bulk[ihplus]=1.0*mol/dm^3
    data.c_bulk[ife2]=0.1*mol/dm^3
    data.c_bulk[ife3]=0.1*mol/dm^3
    data.c_bulk[iso4]=0.75*mol/dm^3

	@assert isapprox(data.c_bulk'*data.z,0, atol=1.0e-12)

	PNPSystem(grid,halfcellbc;data,kwargs...)


end

# ╔═╡ bf4838cd-390b-4706-8903-a35be1915600
function PNPSystem_mem(grid,bcondition;data=ElectrolyteData(),kwargs...)
    sys=VoronoiFVM.System(grid;
                          data=data,
                          flux=pnpflux_mem,
                          reaction=pnpreaction_mem,
                          bcondition,
                          species=[1:data.nc..., data.iϕ],
						  regions=[1,2],
                          kwargs...
                          )
end

# ╔═╡ 44c8d993-1b4f-4e72-ba29-2261d5e6265b
function Setup_mem(;nref=0,R0=1.0e-2,κ=10.0,kwargs...)
	@local_phconstants N_A e R ε_0
    F=N_A*e
    @local_unitfactors cm μF mol dm s mA A nm

	defaults=(; max_round=3,
              tol_round=1.0e-9,
              verbose=false,
              tol_relative=1.0e-8,
              tol_mono=1.0e-10)
	
	kwargs=merge(defaults, kwargs)
    R0=R0*ufac"mol/(cm^2*s)"

    ihplus=1
    ife2 = 2
    ife3 = 3
    iso4 = 4
	
	#Working electrode boundary condition
	function halfcellbc(f,u,bnode,data)
		(;nc,Γ_we,Γ_bulk,ϕ_we,ip,iϕ,v,v0,RT)=data
		bulkbcondition(f,u,bnode,data;region=Γ_bulk)
		nothing
	end

    data=ElectrolyteData(;nc=4,
						 z=[1,2,3,-2],
						 κ=fill(κ,4),
                         Γ_we=1,
                         Γ_bulk=2,
						 R0=R0,
                         )

	data.c_bulk[ihplus]=1.0*mol/dm^3
    data.c_bulk[ife2]=0.1*mol/dm^3
    data.c_bulk[ife3]=0.1*mol/dm^3
    data.c_bulk[iso4]=0.75*mol/dm^3

	@assert isapprox(data.c_bulk'*data.z,0, atol=1.0e-12)

	grid=grid2()
	PNPSystem_mem(grid,halfcellbc;data,kwargs...)
end

# ╔═╡ d72f83c5-ea94-476f-8d9e-ea942e41f851
let
	#sys=Setup()
	sys=Setup_mem()
	data=sys.physics.data
	volts,currs,sols = ivsweep(sys;ispec=ife2)
	
	tsol=VoronoiFVM.TransientSolution(sols,volts)

    for it=1:length(tsol.t)
        tsol.u[it][ihplus,:]/=ufac"mol/dm^3"
		tsol.u[it][ife2,:]/=ufac"mol/dm^3"
        tsol.u[it][ife3,:]/=ufac"mol/dm^3"
		tsol.u[it][iso4,:]/=ufac"mol/dm^3"
    end

	
    vis=GridVisualizer(;resolution=(900,400),layout=(2,3),clear=true)
    #aspect=3.5*20.0/(tsol.t[end]-tsol.t[begin])
	
    #scalarplot!(vis[1,1],F*currs/(ufac"mA/cm^2"),volts,markershape=:none, title="IV",xlabel="I",ylabel="ϕ")
	
    scalarplot!(vis[1,1],sys,sols[end];species=ife2,title="Fe2+",colormap=:summer,ylabel="c / mol dm⁻³")

    scalarplot!(vis[1,2],sys,sols[end];species=ife3,title="Fe3+",colormap=:summer,ylabel="c / mol dm⁻³")
	
    scalarplot!(vis[1,3],sys,sols[end];species=data.iϕ,title="ϕ",colormap=:bwr,ylabel="ϕ")
	
	scalarplot!(vis[2,1],sys,sols[end];species=ihplus,title="H+",colormap=:summer,ylabel="c / mol dm⁻³")

	scalarplot!(vis[2,2],sys,sols[end];species=iso4,title="SO₄²⁻",colormap=:summer,ylabel="c / mol dm⁻³")
    
	reveal(vis)
end

# ╔═╡ Cell order:
# ╟─76b3f5de-efb6-4837-bfa6-c26f7f3fb215
# ╠═4eebb4d0-4ee9-11ed-0892-05b637938ba7
# ╟─bd3af7ca-e2fe-4fce-bd0c-88bc51ae046b
# ╟─b48f424e-a9b7-434c-b5b0-5f953de68664
# ╟─dc1d238c-a04a-477f-a8fb-5b13c6684ff2
# ╟─4bf59f98-0f54-499d-a081-41f456fa5719
# ╟─01a76a61-acc2-463d-8a5a-b8e22392a50f
# ╟─8abace6c-f126-4de3-b92d-484fa31811b9
# ╟─bbb39bea-26a9-45b0-8528-fecd6636a28e
# ╟─1c25bd96-09dc-4593-8616-f6ad2b8c6746
# ╟─b674382d-1c6c-4b8e-891c-2db213ed2bcd
# ╟─0ba193f4-fb72-4a4c-9f35-f46e948df2d7
# ╟─1ecc9a61-5c98-46a1-ba03-14095b33b20c
# ╟─932958ba-a7a2-4100-addb-851816b45b24
# ╠═2bd278e7-898d-46c3-b9eb-5d4edcb060bb
# ╠═0dd23774-4f69-44af-923d-0e7eaefe6bda
# ╠═032f7e2b-9bd2-4ca0-ba41-edadb2ada2ee
# ╠═0c414c82-6801-4b8d-8684-e1ab756c30a5
# ╟─eb5df477-a000-46fc-8518-47d8147df2c2
# ╟─a78ebcf6-eab5-4b54-bb7e-3f6bee5b0795
# ╠═7e0fe75c-b3f9-4cbc-9797-1b2253c31a2b
# ╠═09847b05-b5ef-41fd-ae23-1762c281d69a
# ╠═58498779-d01f-451b-8a65-634baa2975ae
# ╠═744cae5e-ec70-4ae3-9283-1c1d96df729f
# ╟─f097c058-977f-48cb-bad9-0a04319379f3
# ╠═fd07e120-3312-462c-8efa-500f602584fa
# ╠═d784e1d2-17fa-428f-8963-dfd0ddb60461
# ╠═2299256d-ba0c-45b3-a331-1f75fcc02c45
# ╠═44c8d993-1b4f-4e72-ba29-2261d5e6265b
# ╠═f0d7f04a-0b9f-4be7-bb42-8b97ce0a8676
# ╠═5a568ec5-9770-45d3-b95c-db7ba6605110
# ╠═d72f83c5-ea94-476f-8d9e-ea942e41f851
# ╠═a02a4def-e2dc-4c9e-98b8-75524963e242
# ╠═bf4838cd-390b-4706-8903-a35be1915600
# ╠═d9b34916-597f-463a-86b9-f9b2b1856c05
# ╟─f9178a7d-2b2b-4ad1-9ee3-941f4e2f27a6
# ╠═4aaaf629-b265-45b2-b449-99379c3e3cd5
# ╠═75392e35-c597-47b2-8698-cc7951529fa8
# ╠═5dae223d-951a-4746-b5af-27c7a28a42a0
# ╠═8e2b1dd6-af4d-4b83-987b-dda484e01b0c
# ╟─bc02abf0-f64d-4248-8c7e-5c8ac7de9358
# ╠═44e785f4-74fe-4fb0-a0fe-5ee205c5bcbe
# ╟─903d9fe7-ddbd-4b42-8816-d216e46b6bbf
# ╠═2fa560d3-c9b6-4ae6-9c99-cbb6ff8a98ec
# ╟─aa30627a-8359-40d4-8c20-34ea0f7ab775
# ╠═9dd8e4a5-391e-4414-ad47-5247c2803101
# ╠═20f230b0-bc12-4c38-a446-ce7d8607896f
