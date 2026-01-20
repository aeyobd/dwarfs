### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./style.jl")

# ╔═╡ 2c7fe797-4acd-4390-a5df-3a8584c9ee20
import FileIO

# ╔═╡ e8e2eb62-46ed-4367-b6d5-79dababe209e
img_dark_matter = FileIO.load("resources/Aq-fullbox_1024x1024.jpg")[:, 512:end]

# ╔═╡ 59610717-6dae-402e-8386-c234cd2d5799
function fit_image_to_rect(img, target_width, target_height)
    img_height, img_width = size(img)
    img_aspect = img_width / img_height
    target_aspect = target_width / target_height
    
    # Determine crop dimensions to match target aspect ratio
    if img_aspect > target_aspect
        # Image is wider - crop width
        new_width = round(Int, img_height * target_aspect)
        new_height = img_height
        x_offset = (img_width - new_width) ÷ 2
        y_offset = 0
    else
        # Image is taller - crop height
        new_width = img_width
        new_height = round(Int, img_width / target_aspect)
        x_offset = 0
        y_offset = (img_height - new_height) ÷ 2
    end
    
    # Crop image (centered crop)
    cropped = img[y_offset+1:y_offset+new_height, x_offset+1:x_offset+new_width]

	return cropped
end



# ╔═╡ 3612882d-3824-4dea-ad9e-ce3da30d14f8
import ImageTransformations

# ╔═╡ 2def5dfe-41be-4832-8467-e69eea3b3e02
img_ordinary_matter = FileIO.load("resources/2560px-Webb's_First_Deep_Field.jpg")

# ╔═╡ d439795d-912b-450f-a5bd-3c38de665a4c
image(0..0.1, 0..0.2, rotr90(fit_image_to_rect(img_ordinary_matter, 0.1, 0.2)), axis=(; aspect=DataAspect()))

# ╔═╡ dd519973-9e0a-47bd-8bbc-fa368cff27aa
@savefig "LCDM" let
	# Create figure
	fig = Figure(size=(668, 668), figure_padding=0)
	ax = Axis(fig[1, 1], aspect=DataAspect())
	
	# Golden ratio
	φ = (1 + 1) / 2
	
	# Overall rectangle dimensions (width = φ, height = 1 for golden ratio)
	width = φ
	height = 1.0
	total_area = width * height
	
	# Area proportions
	areas = [0.68, 0.27, 0.05]
	
	# Calculate rectangle dimensions
	# Using a simple layout: stack them horizontally
	area1 = areas[1] * total_area
	area2 = areas[2] * total_area
	area3 = areas[3] * total_area
	
	# All rectangles have the same height
	w1 = area1 / height
	w2 = width - w1

	h1 = area2 / (area2 + area3)
	
	# Define rectangles as (x_start, y_start, width, height)
	rectangles = [
	    (0, 0, w1, height),
	    (w1, 0, w2, h1),
	    (w1, h1, w2, height - h1)
	]

	labels = ["dark energy", "dark\nmatter", "ordinary\nmatter"]
	colors = [COLORS[5], :transparent, :transparent]
	# Draw rectangles
	for (i, (x, y, w, h)) in enumerate(rectangles)
		@info w, h

		if i == 2
			image!(ax, x..(x+w), y..(y+h), rotr90(fit_image_to_rect(img_dark_matter, w, h)))
		end

		if i == 3
			image!(ax, x..(x+w), y..(y+h), rotr90(fit_image_to_rect(img_ordinary_matter, w, h)))
		end
		
	    poly!(ax, 
	          Point2f[(x, y), (x+w, y), (x+w, y+h), (x, y+h)],
	          color=colors[i],
	          strokecolor=:white,
	          strokewidth=theme(:linewidth))



	    
	    # Add percentage labels
	    text!(ax, x + w/2, y + h/2, 
	          text="$(Int(round(areas[i]*100)))% $(labels[i])",
	          align=(:center, :center),
	          fontsize=theme(:fontsize)[] * 0.8,
	          color=:white,
	          font=:bold)
	end
	
	# Set axis limits and hide decorations
	xlims!(ax, -0.0, width + 0.0)
	ylims!(ax, -0.0, height + 0.0)
	hidedecorations!(ax)
	hidespines!(ax)
	resize_to_layout!()
	
	fig
end

# ╔═╡ 1e6bea19-8010-4549-9f9c-0de6b5d4051a


# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═2c7fe797-4acd-4390-a5df-3a8584c9ee20
# ╠═e8e2eb62-46ed-4367-b6d5-79dababe209e
# ╠═59610717-6dae-402e-8386-c234cd2d5799
# ╠═3612882d-3824-4dea-ad9e-ce3da30d14f8
# ╠═2def5dfe-41be-4832-8467-e69eea3b3e02
# ╠═d439795d-912b-450f-a5bd-3c38de665a4c
# ╠═dd519973-9e0a-47bd-8bbc-fa368cff27aa
# ╠═1e6bea19-8010-4549-9f9c-0de6b5d4051a
