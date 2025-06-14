using Distributions
using Random
using PyPlot

include("gof_statistics.jl")
mkpath("results/nulls")
colors = ["#3b7c70", "#ce9642", "#898e9f", "#3b3a3e"]

function plot_ecdfs(no_correction_ps, ps, fn)
    figure(figsize=(8, 5))
    clf();

    # no correction
    subplot(1, 2, 1)
    sorted_no_corr = sort(no_correction_ps)
    n_sims = length(sorted_no_corr)
    plot(sorted_no_corr, range(1/n_sims, 1, length=n_sims), color=colors[1])
    plot([0, 1], [0, 1], "k--", linewidth=1)  # reference line
    xlim(0,1.005)
    ylim(0,1)
    xlabel(raw"$\alpha$", fontsize = 20)
    ylabel("cdf", fontsize = 20)
    title("no correction", fontsize = 20)

    # with correction
    subplot(1, 2, 2)
    sorted_corr = sort(ps)
    plot(sorted_corr, range(1/n_sims, 1, length=n_sims), color=colors[2])
    plot([0, 1], [0, 1], "k--", linewidth=1)  # reference line
    xlim(0,1.005)
    ylim(0,1)
    xlabel(raw"$\alpha$", fontsize = 20)
    ylabel("cdf", fontsize = 20)
    title("with correction", fontsize = 20)

    tight_layout()
    savefig("results/nulls/$(fn)_ecdfs.png", dpi=200)
end

function plot_truncated_ecdfs(no_correction_ps, ps, fn, t)
    figure(figsize=(8, 5))
    clf();

    # no correction
    subplot(1, 2, 1)
    sorted_no_corr = sort(no_correction_ps)
    n_sims = length(sorted_no_corr)
    plot(sorted_no_corr, range(1/n_sims, 1, length=n_sims), color=colors[1])
    plot([0, t], [0, t], "k--", linewidth=1)  # reference line
    xlim(0, t)
    ylim(0, max(t, findfirst((sorted_no_corr.>=0.1))/n_sims))
    xlabel(raw"$\alpha$", fontsize = 20)
    ylabel("cdf", fontsize = 20)
    title("no correction", fontsize = 20)

    # with correction
    subplot(1, 2, 2)
    sorted_corr = sort(ps)
    plot(sorted_corr, range(1/n_sims, 1, length=n_sims), color=colors[2])
    plot([0, t], [0, t], "k--", linewidth=1)  # reference line
    xlim(0, t)
    ylim(0, max(t, findfirst((sorted_corr.>=0.1))/n_sims))
    xlabel(raw"$\alpha$", fontsize = 20)
    ylabel("cdf", fontsize = 20)
    title("with correction", fontsize = 20)

    tight_layout()
    savefig("results/nulls/$(fn)_trunc_ecdfs.png", dpi=200)
end

function plot_hists(no_correction_ps, ps, fn)
    figure(figsize=(8, 5))
    clf();

    # no correction
    subplot(2, 1, 1)
    hist(no_correction_ps, bins=100, alpha=1, density = true, color = colors[1], edgecolor="black")
    xlabel("p-value")
    ylabel("density")
    title("no correction")

    # with correction
    subplot(2, 1, 2)
    hist(ps, bins=100, alpha=1, density = true, color = colors[2], edgecolor="black")
    xlabel("p-value")
    ylabel("density")
    title("with correction")

    tight_layout()
    savefig("results/nulls/$(fn)_hists.png", dpi=200)
end

n = 30;
n_sims = 10^5;

HP(x) = harmonic_PITCOS(x, false)

methods = [HP];
method_labels = [raw"$\mathtt{PITCOS}_{\mathrm{Harmonic}}$" #,
];
fig_names = ["PITCOS_Harmonic"]

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
    
    plot_ecdfs(no_correction_ps, ps, fn)
    plot_hists(no_correction_ps, ps, fn)
    plot_truncated_ecdfs(no_correction_ps, ps, fn, 0.1)
end 