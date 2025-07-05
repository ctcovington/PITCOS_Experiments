using HypothesisTests
using Distributions
using Random

# Kolmogorov-Smirnov
KS(x) = pvalue(HypothesisTests.ExactOneSampleKSTest(x, Uniform(0,1)))

# Anderson-Darling
AD(x) = pvalue(HypothesisTests.OneSampleADTest(x, Uniform(0,1)))

# Shapiro-Wilk
SW(x) = pvalue(HypothesisTests.ShapiroWilkTest(quantile.(Normal(0,1),x)))

# Neyman-Barton
# (This is N_2 in Blinov and Lemeshko (2014).)
function NB(x)
    n = length(x)
    y = x .- 0.5
    V1 = sqrt(n)*mean(2*sqrt(3)*y)
    V2 = sqrt(n)*mean(sqrt(5)*(6*y.^2 .- 0.5))
    t = V1^2 + V2^2
    p = ccdf(Chisq(2),t)
    return p
end

# # PITOS
# cauchy_combination(pvalues) = ccdf(Cauchy(), mean(tan.(pi*(0.5 .- pvalues))))

# function PITOS(x, correction)
#     if correction
#         cf = 1.15
#     else
#         cf = 1
#     end
#     xo = sort(x)
#     n = length(xo)
#     j = (1:n)
#     u = [ccdf(Binomial(n,xo[j]), j-1) for j=1:n]
#     p = 2*min.(u, 1 .- u)
#     t = min(1, cf*cauchy_combination(p))
#     return t
# end

# function indexed_PITCOS(xo, n, indices)
#     m = length(indices)
#     u = zeros(m)

#     u[1] = cdf( Beta(indices[1], n - indices[1] + 1), xo[indices[1]] ) # map through marginal distribution
#     for i in 1:(m-1)
#         # map through cdf for C_{(i+1)} | C_{(i)}, where C = X[indices]
#         if indices[i] < indices[i+1]
#             # if consecutive xo values are identical, default to extreme u value
#             if xo[indices[i+1]] == xo[indices[i]]
#                 u[i+1] = 0
#             else
#                 u[i+1] = cdf( Beta(indices[i+1] - indices[i], n - indices[i+1] + 1), 
#                             (xo[indices[i+1]] - xo[indices[i]])/(1-xo[indices[i]]) 
#                             )
#             end
#         else
#             if xo[indices[i+1]] == xo[indices[i]]
#                 u[i+1] = 0
#             else
#                 u[i+1] = cdf( Beta(indices[i+1], indices[i] - indices[i+1]), 
#                             xo[indices[i+1]] / xo[indices[i]] 
#                             )
#             end
#         end
#     end
#     ps = 2*min.(u, 1 .- u);
#     p = cauchy_combination(ps);
#     return p;
# end

# function PITCOS_1(x, correction=false)
#     xo = sort(x)
#     n = length(x)
#     p = indexed_PITCOS(xo, n, 1:n)
#     return p
# end

# function pairs_PITCOS(x, correction=true)
#     if correction
#         cf = 1.15
#     else
#         cf = 1
#     end
#     xo = sort(x)
#     n = length(x);
#     m = Int64(n*(n+1)/2);
#     ps = zeros(m);
#     pairs = [(i, j) for i in 1:n for j in i+1:n];
#     for i in 1:n
#         ps[i] = indexed_PITCOS(xo, n, i);
#     end
#     for (ip,pair) in enumerate(pairs)
#         ps[n+ip] = indexed_PITCOS(xo, n, pair);
#     end

#     p = min(1, cf*cauchy_combination(ps)); 
#     return p
# end

# function scan_PITCOS(x, correction=true)
#     if correction
#         cf = 1.15
#     else
#         cf = 1
#     end
#     xo = sort(x)
#     n = length(x);
#     m = n-1;
#     ps = []
#     for k in 1:m
#         for j in 1:(min(k,n-k))
#             v = collect(j:k:n)
#             push!(ps, indexed_PITCOS(xo, n, v))
#         end
#     end
#     p = min(1, cf*cauchy_combination(ps))
#     return p
# end

# function harmonic_PITCOS(x, correction=true)
#     if correction
#         cf = 1.15
#     else
#         cf = 1
#     end
#     xo = sort(x)
#     n = length(x);
#     m = n-1;
#     ps = zeros(m)
#     for k in 1:m
#         s = rand(1:min(k, n-k))
#         v = collect(s:k:n)
#         ps[k] = indexed_PITCOS(xo, n, v)
#     end
#     p = min(1, cf*cauchy_combination(ps))
#     return p
# end

# function comb_test(x, correction=true)
#     if correction
#         cf = 1.15
#     else
#         cf = 1
#     end
#     p1 = harmonic_PITCOS(x, false);
#     p2 = AD(x);
#     p3 = NB(x);
#     ps = [p1,p2,p3]
#     p = min(1, cf*cauchy_combination(ps))
#     return p
# end