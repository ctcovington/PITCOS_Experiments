using Distributions
using Random
using PyPlot

using Pkg 
Pkg.add(path="/Users/christiancovington/.julia/dev/PITCOS")
using PITCOS
mkpath("../../results/power_calcs")

# Set logging level to only show warnings and above
using Logging
Logging.disable_logging(LogLevel(Logging.Info))  # Only show warnings and above


# PITCOS-based methods
# automatically set correction to false because we are showing false and true in the plots
PITOS(x) = pitcos(x, index_sets = [[i] for i in 1:length(x)], correction=false)
PITCOS_1(x) = pitcos(x, index_sets = [1:length(x)], correction=false)
HP(x) = pitcos(x, method="Harmonic", correction=false)

include("gof_statistics.jl")
mkpath("../../results/nulls")
colors = ["#3b7c70", "#ce9642", "#898e9f", "#3b3a3e"]

function plot_ecdfs(no_correction_ps, ps, label, fn)
    figure(figsize=(8, 5))
    clf();

    # First CDF: No correction
    subplot(1, 2, 1)
    sorted_no_corr = sort(no_correction_ps)
    n_sims = length(sorted_no_corr)
    plot(sorted_no_corr, range(1/n_sims, 1, length=n_sims), color=colors[1])
    plot([0, 1], [0, 1], "k--", linewidth=1)  # reference line
    xlim(0,1.005)
    ylim(0,1)
    xlabel(raw"$\alpha$", fontsize = 20)
    ylabel("cdf", fontsize = 20)
    title("$(label)\nno correction", fontsize = 20)

    # Second CDF: With correction
    subplot(1, 2, 2)
    sorted_corr = sort(ps)
    plot(sorted_corr, range(1/n_sims, 1, length=n_sims), color=colors[2])
    plot([0, 1], [0, 1], "k--", linewidth=1)  # reference line
    xlim(0,1.005)
    ylim(0,1)
    xlabel(raw"$\alpha$", fontsize = 20)
    ylabel("cdf", fontsize = 20)
    title("$(label)\nwith correction", fontsize = 20)

    tight_layout()
    savefig("../../results/nulls/$(fn)_ecdfs.png", dpi=200)
end

function plot_truncated_ecdfs(no_correction_ps, ps, label, fn, t)
    figure(figsize=(8, 5))
    clf();

    # First CDF: No correction
    subplot(1, 2, 1)
    sorted_no_corr = sort(no_correction_ps)
    n_sims = length(sorted_no_corr)
    plot(sorted_no_corr, range(1/n_sims, 1, length=n_sims), color=colors[1])
    plot([0, t], [0, t], "k--", linewidth=1)  # reference line
    xlim(0, t)
    ylim(0, max(t, findfirst((sorted_no_corr.>=0.1))/n_sims))
    xlabel(raw"$\alpha$", fontsize = 20)
    ylabel("cdf", fontsize = 20)
    title("$(label)\nno correction", fontsize = 20)

    # Second CDF: With correction
    subplot(1, 2, 2)
    sorted_corr = sort(ps)
    plot(sorted_corr, range(1/n_sims, 1, length=n_sims), color=colors[2])
    plot([0, t], [0, t], "k--", linewidth=1)  # reference line
    xlim(0, t)
    ylim(0, max(t, findfirst((sorted_corr.>=0.1))/n_sims))
    xlabel(raw"$\alpha$", fontsize = 20)
    ylabel("cdf", fontsize = 20)
    title("$(label)\nwith correction", fontsize = 20)

    tight_layout()
    savefig("../../results/nulls/$(fn)_trunc_ecdfs.png", dpi=200)
end

function plot_hists(no_correction_ps, ps, label, fn)
    figure(figsize=(8, 5))
    clf();

    # First CDF: No correction
    subplot(2, 1, 1)
    hist(no_correction_ps, bins=100, alpha=1, density = true, color = colors[1], edgecolor="black")
    xlabel("p-value")
    ylabel("density")
    title("$(label)\nno correction")

    # Second CDF: With correction
    subplot(2, 1, 2)
    hist(ps, bins=100, alpha=1, density = true, color = colors[2], edgecolor="black")
    xlabel("p-value")
    ylabel("density")
    title("$(label)\nwith correction")

    tight_layout()
    savefig("../../results/nulls/$(fn)_hists.png", dpi=200)
end

n = 30;
n_sims = 10^5;

methods = [AD, KS, SW, NB, PITOS, PITCOS_1, HP];
method_labels = ["Anderson-Darling",
            "Kolmogorov-Smirnov",
            "Shapiro-Wilk",
            "Neyman-Barton",
            raw"$\mathtt{PITOS}$",
            raw"$\mathtt{PITCOS}_1$",
            raw"$\mathtt{PITCOS}_{\mathrm{Harmonic}}$" #,
];
fig_names = ["AD", "KS", "SW", "NB", "PITOS", "PITCOS_1", "PITCOS_Harmonic"]
# fig_names = ["PITOS", "PITCOS_1", "PITCOS_Pairs", "PITCOS_Scan", "PITCOS_Harmonic", "Combined Test"]

for (method, label, fn) in zip(methods, method_labels, fig_names)
    no_correction_ps = zeros(n_sims);
    ps = zeros(n_sims);

    for i in 1:n_sims
        # sample from Unif(0,1)
        x = rand(Beta(1,1), n)

        # get p-values with/without correction
        p = method(x)
        no_correction_ps[i] = p 
        ps[i] = min(1, 1.15*p)
    end
    
    plot_ecdfs(no_correction_ps, ps, label, fn)
    plot_hists(no_correction_ps, ps, label, fn)
    plot_truncated_ecdfs(no_correction_ps, ps, label, fn, 0.1)
end 