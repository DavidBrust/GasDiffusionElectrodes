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
"""

# ╔═╡ Cell order:
# ╠═23dc6930-c354-11ed-110a-235339eca423
# ╟─5ede4875-1dbf-4930-9538-39b10d4bd377
# ╟─1593e889-7b4f-40f9-b6e8-ef0564381c81
# ╟─711a3678-e3c0-4193-a6cf-0d525ddeff91
# ╟─b55b90f2-bace-485c-a745-3fca9871152f
# ╟─47fea77b-d605-40ac-b436-078fd91294fe
# ╟─4a33c6a3-ad59-4efc-997e-09ca780e5fd2
# ╟─2aefd6b8-f35c-419e-a2c4-83ad44f1be77
# ╟─67a38c46-a952-473d-864e-53ca99b16174
# ╟─a4bb4cff-e933-45ac-8ad2-c2f69328f62d
# ╟─87c2c676-9bf0-4716-8fe4-e44a0f884412
# ╟─72514d88-d17b-4dfe-8921-623aa6e3b42f
# ╟─3ede7aa4-d924-4825-b0cf-f5f9b52c5e3b
# ╟─21d032ab-3de3-400a-b6fb-2c9805e4916b
# ╟─ac882611-643c-4c94-9acd-38de7566498e
# ╟─0dbfa843-7b38-41bf-aa4c-14d528bbf75f
# ╟─25ba905a-3342-4139-935d-ebcf0ad4be56
# ╟─10f5b322-e52b-45b7-ae60-a9f5d82b5f99
