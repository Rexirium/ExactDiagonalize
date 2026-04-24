using CairoMakie

include("haldane.jl")

function plot_Haldane_spectrum(Lx::Int, nky::Int, t1::Real, t2::Number; m2::Real=0.0) 
    kys = range(0, 2π, nky)
    spectra = Matrix{Float64}(undef, Lx, nky)

    H = zeros(Lx, Lx)
    for (i, ky) in enumerate(kys)
        updateHaldaneHamiltonian!(H, Lx, ky, t1, t2; m2 = m2)
        spectra[:, i] = eigvals(H)
    end

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], title="Haldane Spectrum", 
        xlabel=L"k_y a_y", ylabel=L"E", 
        xticks=(0 : π/2 : 2π, [L"0", L"π/2", L"π", L"3π/2", L"2π"])
    )

    for j in 1:Lx
        lines!(ax, kys, spectra[j, :], color=:black, linewidth=1)
    end
    vlines!(ax, [π / 2], color=:red, linestyle=:dash)
    vlines!(ax, [3π/2], color=:red, linestyle=:dash)
    fig
end

function plot_Haldane_edgestates(Lx::Int, kys::Vector{Float64}, t1::Real, t2::Number; m2=0.0) 
    fig = Figure(size=(800, 900))

    for (row, kyc) in enumerate(kys)
        H = makeHaldaneHamiltonian(Lx, kyc * π, t1, t2; m2=m2)
        edgestates = get_edgestates(H)

        ax = Axis(fig[row, 1], 
            title="Haldane Edge States at k_y a = $(kyc)π", 
            xlabel=L"j", ylabel=L"|ψ_j|^2", 
            xticks = 0 : 10 : Lx
        )
        deg = size(edgestates, 2)
        for i in 1:deg
            barplot!(ax, 1 : Lx, abs2.(edgestates[:, i]), label="Edge State $i")
        end
        axislegend(ax; position=:ct)
    end
    fig   
end

let 
    set_theme!(Axis=(
        xtickalign = 1,
        ytickalign = 1,
        xticklabelsize = 16,
        yticklabelsize = 16,
        xlabelsize = 18,
        ylabelsize = 18,
    ))

    fig = plot_Haldane_spectrum(50, 201, 1.0, 0.2im)
    
end