using HypothesisTests
using Distributions
using PyPlot
using SpecialFunctions
using StatsBase
using Random
using PyCall
using JLD2
@pyimport matplotlib.lines as mlines

include("gof_statistics.jl")
using Pkg 
Pkg.add(path="/Users/christiancovington/.julia/dev/PITCOS")
using PITCOS
mkpath("../../results/power_calcs")

# Set logging level to only show warnings and above
using Logging
Logging.disable_logging(LogLevel(Logging.Info))  # Only show warnings and above


# PITCOS-based methods
PITOS(x) = pitcos(x, index_sets = [[i] for i in 1:length(x)], correction=true)
PITCOS_1(x) = pitcos(x, index_sets = [1:length(x)], correction=false)
HP(x) = pitcos(x, method="Harmonic", correction=true)

rand_B_1_1(n) = rand(Beta(1,1),n)
rand_B_1_2(n) = rand(Beta(1,2),n)
rand_B_2_2(n) = rand(Beta(2,2),n)
rand_B_p5_p5(n) = rand(Beta(0.5,0.5),n)
rand_B_p8_p8(n) = rand(Beta(0.8,0.8),n)
rand_B_p9_p9(n) = rand(Beta(0.9,0.9),n)
rand_B_1p1_1p1(n) = rand(Beta(1.1,1.1),n)
rand_B_1p2_1p2(n) = rand(Beta(1.2,1.2),n)
rand_B_1p3_p8(n) = rand(Beta(1.3,0.8),n)
rand_B_p7_1(n) = rand(Beta(0.7,1),n)
rand_B_p9_1(n) = rand(Beta(0.9,1),n)
rand_B_1p1_1(n) = rand(Beta(1.1,1),n)

rand_mix_spike_right(n) = (z=rand(Bernoulli(0.9),n); z.*rand(Beta(1,1),n) .+ (1 .- z).*rand(Beta(20,4),n))
rand_no_edges(n) = rand(Beta(1,1),n)*0.9 .+ 0.05
rand_no_middle(n) = (z=rand(Bernoulli(0.5),n); z.*rand(Uniform(0,0.45),n) .+ (1 .- z).*rand(Uniform(0.55,1),n))
rand_outliers(n) = (x = rand(Beta(1,1),n); m = ceil(Int,n/50); x[1:m] = x[1:m]*1e-5; x)
rand_t10(n) = cdf.(Normal(), rand(TDist(10), n))
rand_laplace(n) = cdf.(Normal(), rand(Laplace(),n))
discrete_50(n) = rand(Categorical(fill(1/50, 50)), n) ./ 50

# list alternative hypotheses
Rs = [rand_B_1_1,
    #   rand_B_1_2,
      rand_B_2_2,
    #   rand_B_p8_p8,
    #   rand_B_1p2_1p2,
      rand_B_1p3_p8,
      rand_B_p9_1,
      rand_B_1p1_1,
      rand_mix_spike_right,
      rand_no_edges,
      rand_no_middle,
      discrete_50,
      rand_outliers,
      rand_t10,
      rand_laplace
];
R_labels = ["Uniform(0,1)",
            # "Beta(1,2)",
            "Beta(2,2)",
            # "Beta(0.8,0.8)",
            # "Beta(1.2,1.2)",
            "Beta(1.3,0.8)",
            "Beta(0.9,1)",
            "Beta(1.1,1)",
            "Spike Right",
            "No Edges",
            "No Middle",
            "Discrete",
            "Outliers",
            raw"t$_{10}$",
            "Laplace"
];

nRs = length(Rs);

# list tests
Ts = [AD, KS, SW, NB, PITOS, PITCOS_1, HP];
T_labels = ["Anderson-Darling",
            "Kolmogorov-Smirnov",
            "Shapiro-Wilk",
            "Neyman-Barton",
            raw"$\mathtt{PITOS}$",
            raw"$\mathtt{PITCOS}_1$",
            raw"$\mathtt{PITCOS}_{\mathrm{Harmonic}}$" #,
];
nTs = length(Ts);

# list sample sizes
ns = [20, 40, 60, 80, 100, 200, 400];
nNs = length(ns);

# set alpha level and number of simulations
alpha = 0.05;
n_sims = 10^5;

power = zeros( nTs, nRs, nNs );
cb_palette = [
    "#E41A1C",  # Red
    "#377EB8",  # Blue
    "#4DAF4A",  # Green
    "#984EA3",  # Purple
    "#FF7F00",  # Orange
    "#00CED1",  # Teal
    "#A65628"   # Brown
]

for (i_R, R) in enumerate(Rs)
    for (i_n, n) in enumerate(ns)
        println("Computing power for $(R_labels[i_R]) with n = $n ...")
        
        x = [R(n) for rep = 1:n_sims]
        for (i_T, T) in enumerate(Ts)
            power[i_T, i_R, i_n] = mean([(T(x[rep]) .< alpha) for rep = 1:n_sims])
        end
    end
    # Plot to individual file
    figure(1); clf()
    for (i_T, T) in enumerate(Ts)
        plot(ns, vec(power[i_T, i_R, :]), label = T_labels[i_T], color=cb_palette[i_T])
    end
    xlabel("n")
    ylabel("power")
    legend()
    savefig("../../results/power_calcs/power-$(R_labels[i_R]).png", dpi=200)
end

# generate vectors of samples from each R
mkpath("../../results/power_calcs/samples")
if !isfile("../../results/power_calcs/samples.jld2")
    Rs_samples = [R(10^7) for R in Rs]
    @save "../../results/power_calcs/samples.jld2" Rs_samples
else 
    Rs_samples = load("../../results/power_calcs/samples.jld2")
    Rs_samples = Rs_samples["Rs_samples"]
end

# Create a new figure
figure(figsize=(8, 36))
clf()

n_row = 8
n_col = 3

# Define height ratios: top 4 rows = 1, bottom 4 rows = 1.5
row_heights = vcat(fill(1.0, 4), fill(1.5, 4))  # [1, 1, 1, 1, 1.5, 1.5, 1.5, 1.5]

gs = matplotlib.gridspec.GridSpec(n_row, n_col, height_ratios=row_heights)

for (i_R, R) in enumerate(Rs)
    row_hist = div(i_R-1, n_col) + 1
    col = mod1(i_R, n_col)
    ax1 = subplot(gs[row_hist, col])
    hist(Rs_samples[i_R], bins=100, color="blue")
    xlim(-0.01, 1.01)
    xlabel("")
    xticks([])
    ylabel("")
    yticks([])
    title(R_labels[i_R], fontsize=8, y=1.05, pad=2)  # Increased y and added padding

    # Bottom 4 rows: rows 4–7 used for power curves (i_R + 12 in 13–24)
    row_power = row_hist + 4
    ax2 = subplot(gs[row_power, col])
    for (i_T, T) in enumerate(Ts)
        plot(ns, vec(power[i_T, i_R, :]),
             label=T_labels[i_T],
             color=cb_palette[i_T],
             lw=0.8)
    end

    title(R_labels[i_R], fontsize=8, y=1.05, pad=2)  # Increased y and added padding

    if row_power == 8  
        xlabel("n", fontsize=6)
    else
        xticks([])
    end
    
    if col == 1
        ylabel("power", fontsize=6)
    end
    tick_params(labelsize=4)
end

tight_layout(pad=0.3, w_pad=0.2, h_pad=0.2)
legend_lines = [mlines.Line2D([], [], color=cb_palette[i], label=T_labels[i]) for i in 1:length(Ts)]
figlegend(legend_lines, T_labels, loc="lower center", ncol=3, fontsize=10, bbox_to_anchor=(0.5, -0.02))
subplots_adjust(bottom=0.1)

savefig("../../results/power_calcs/data-power-grid.png", dpi=300, bbox_inches="tight")


# # Create a new figure 
# figure(figsize=(8, 36)) 
# clf();

# n_row = 8
# n_col = 3

# for (i_R, R) in enumerate(Rs)
#     # Plot histogram of R(n)
#     subplot(n_row, n_col, i_R)
#     hist(Rs_samples[i_R], bins=100, color="blue")
#     xlim(-0.01, 1.01)
#     xlabel("")
#     xticks([])
#     ylabel("")
#     yticks([])
#     title(R_labels[i_R], fontsize=8, y=1.05, pad=2)  # Increased y and added padding  # lower title to avoid overlap

#     # Plot power curves
#     subplot(n_row, n_col, i_R + 12)
#     for (i_T, T) in enumerate(Ts)
#         plot(ns, vec(power[i_T, i_R, :]), 
#              label=T_labels[i_T], 
#              color=cb_palette[i_T],
#              lw = 0.8
#              )
#     end
#     title(R_labels[i_R], fontsize=8, y=1.05, pad=2)  # Increased y and added padding
#     if i_R + 12 >= (24-n_col)
#         xlabel("n", fontsize=6)
#     else
#         xticks([])
#     end
#     if i_R % n_col == 1  # Only add y-label for the first column
#         ylabel("power", fontsize=6)
#     end
#     tick_params(labelsize=4)
# end
# tight_layout(pad=0.3, w_pad=0.2, h_pad=0.2)
# legend_lines = [mlines.Line2D([], [], color=cb_palette[i], label=T_labels[i]) for i in 1:length(Ts)]
# figlegend(legend_lines, T_labels, loc="lower center", ncol=3, fontsize=10, bbox_to_anchor=(0.5, -0.02))
# subplots_adjust(bottom=0.1)  # reserve space

# savefig("../../results/power_calcs/data-power-grid.png", dpi=300, bbox_inches="tight")

# #=
# experimenting w new figure
# =#
# # Create a figure with one row per distribution
# n_distributions = length(Rs)
# fig, axes = plt.subplots(n_distributions, 2, figsize=(12, 5*n_distributions), 
#                        gridspec_kw={"width_ratios": [1, 2]},  # Make right column wider
#                        sharey=False)  # Don't share y-axes

# for (i_R, R) in enumerate(Rs)

# end

# for (i_R, R) in enumerate(Rs)
#     # Left column: Histogram
#     ax1 = axes[i_R, 1]
#     ax1.hist(Rs_samples[i_R], bins=50, color="blue", alpha=0.7, density=true)
#     ax1.set_title(R_labels[i_R], fontsize=12, pad=10)
#     ax1.set_xticks([])
#     ax1.grid(True, alpha=0.3)
#     ax1.set_ylabel("density", fontsize=10)
    
#     # Right column: Power curves
#     ax2 = axes[i_R, 2]  # Right column
#     for (i_T, T) in enumerate(Ts)
#         ax2.plot(ns, vec(power[i_T, i_R, :]), 
#                 label=T_labels[i_T], 
#                 color=cb_palette[i_T],
#                 lw=2.0,  # Thicker lines
#                 alpha=0.9)  # Slightly more opaque
#     end
    
#     # Customize axes
#     ax2.set_ylim(0, 1.05)  # Slight padding at top
#     ax2.grid(true, alpha=0.3)
#     ax2.set_ylabel("power", fontsize=10)
    
#     # Add x-label to bottom plot only
#     if i_R == n_distributions - 1  # 0-based indexing
#         ax2.set_xlabel("n", fontsize=10)
#     else
#         ax2.set_xticklabels([])  # Hide x-tick labels for non-bottom plots
#     end
    
#     # Add legend to first subplot only to save space
#     if i_R == 0
#         ax2.legend(fontsize=9, loc="upper left", bbox_to_anchor=(1.02, 1),
#                   borderaxespad=0., frameon=False)
#     end
# end

# # Adjust layout with more padding
# plt.subplots_adjust(hspace=0.4, wspace=0.4)  # More space between subplots

# # Save with high resolution and tight bounding box
# savefig("../../results/power_calcs/improved_power_grid.png", 
#         dpi=300, 
#         bbox_inches="tight",
#         facecolor="white")
