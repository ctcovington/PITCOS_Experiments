using HypothesisTests
using Distributions
using PyPlot
using SpecialFunctions
using StatsBase
using Random
using PyCall
@pyimport matplotlib.lines as mlines

include("gof_statistics.jl")
mkpath("results/power_calcs")


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
      rand_B_1_2,
      rand_B_2_2,
      rand_B_p8_p8,
      rand_B_1p2_1p2,
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
            "Beta(1,2)",
            "Beta(2,2)",
            "Beta(0.8,0.8)",
            "Beta(1.2,1.2)",
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
Ts = [AD, KS, SW, NB, harmonic_PITCOS];
T_labels = ["Anderson-Darling",
            "Kolmogorov-Smirnov",
            "Shapiro-Wilk",
            "Neyman-Barton",
            raw"$\mathtt{PITCOS}_{\mathrm{Harmonic}}$" #,
];
nTs = length(Ts);

# list sample sizes
ns = [20, 40, 60, 80, 100];
nNs = length(ns);

# set alpha level and number of simulations
alpha = 0.05;
n_sims = 10^5;

power = zeros( nTs, nRs, nNs );
cb_palette = ["#332288",  # dark blue
          "#88CCEE",  # light blue
          "#117733",  # green
          "#DDCC77",  # yellow
          "#CC6677"#,  # red
        ];

# calculate power and save to individual files
for (i_R, R) in enumerate(Rs)
    for (i_n, n) in enumerate(ns)
        println("Computing power for $(R_labels[i_R]) with n = $n ...")
        
        x = [R(n) for rep = 1:n_sims]
        for (i_T, T) in enumerate(Ts)
            power[i_T, i_R, i_n] = mean([(T(x[rep]) .< alpha) for rep = 1:n_sims])
        end
    end

    figure(1); clf()
    for (i_T, T) in enumerate(Ts)
        plot(ns, vec(power[i_T, i_R, :]), label = T_labels[i_T], color=cb_palette[i_T])
    end
    xlabel("n")
    ylabel("power")
    legend()
    savefig("results/power_calcs/power-$(R_labels[i_R]).png", dpi=200)
end

# create main figure
clf();
for (i_R, R) in enumerate(Rs)
    # Plot histogram of R(n)
    subplot(6, 5, i_R)
    hist(R(10^7), bins=100, color="blue")
    xlim(-0.01, 1.01)
    xlabel("")
    xticks([])
    ylabel("")
    yticks([])
    title(R_labels[i_R], fontsize=8, y=0.9)  # lower title to avoid overlap

    # plot power curves
    subplot(6, 5, i_R + 15)
    for (i_T, T) in enumerate(Ts)
        plot(ns, vec(power[i_T, i_R, :]), 
             label=T_labels[i_T], 
             color=cb_palette[i_T],
             lw = 0.8
             )
    end
    title(R_labels[i_R], fontsize=8, y=0.9)
    if i_R + 15 > 25  # i_R + 15 in 26–30 => last row
        xlabel("n", fontsize=6)
    else
        xticks([])
    end
    ylabel("power", fontsize=6)
    tick_params(labelsize=4)
end
tight_layout(pad=0.3, w_pad=0.2, h_pad=0.2)
legend_lines = [mlines.Line2D([], [], color=cb_palette[i], label=T_labels[i]) for i in 1:length(Ts)]
figlegend(legend_lines, T_labels, loc="lower center", ncol=3, fontsize=10, bbox_to_anchor=(0.5, -0.02))
subplots_adjust(bottom=0.2)

savefig("results/power_calcs/data-power-grid.png", dpi=300, bbox_inches="tight")