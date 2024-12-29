using LinearAlgebra, Plots, Random

function randomVector()
    θ = 2*π*rand()
    return [cos(θ), sin(θ)]
end

generateGridPoints(xlim, ylim) = [[x, y] for x in 0:xlim, y in 0:ylim]

function gridLines(gridPoints)
    n, m = size(gridPoints)
    verLines = Vector{Vector{}}(undef, n*(m-1))
    for i in 1:n, j in 1:m-1
        verLines[(i-1)*(n-1) + j] = [Tuple(gridPoints[i, j]), Tuple(gridPoints[i, j+1])]
    end

    horLines = Vector{Vector{}}(undef, m*(n-1))
    for i in 1:n-1, j in 1:m
        horLines[(j-1)*(m-1) + i] = [Tuple(gridPoints[i, j]), Tuple(gridPoints[i+1, j])]
    end

    skewLines = Vector{Vector{}}(undef, 2*sum(1:m-2)+m-1)
    c = 1
    for i in 2:m-1, j in 1:m-i
        skewLines[c] = [Tuple(gridPoints[i+j-1, j]), Tuple(gridPoints[i+j, j+1])]
        skewLines[length(skewLines)-c+1] = [Tuple(gridPoints[n-(i+j-1)+1, m-j+1]), Tuple(gridPoints[n-(i+j)+1, m-(j+1)+1])]
        c += 1
    end
    for j in 1:m-1
        skewLines[c-1+j] = [Tuple(gridPoints[j, j]), Tuple(gridPoints[j+1, j+1])]
    end
    return [verLines; horLines; skewLines] 
end

function plotGrid!(gridPoints, color=1)
    lines = gridLines(gridPoints)
    for line in lines
        plot!(line, legend=false, c=color)
    end
    plot!()
end

linearInterpolation(a0, a1, w) = (a1 - a0) * w + a0

function smoothstep(x; n=1)
    if x <= 0
        return 0
    elseif x <= 1
        #sum([binomial(n + k, k) * binomial(2n + 1, n - k) .* (-x)^k for k in 0:n])
        #sum(binomial.(n .+ 0:n, 0:n) .* binomial.(2n + 1, n .- 0:n) .* (-x).^(0:n))
        return x^(n+1) * sum([binomial(n + k, k) * binomial(2n + 1, n - k) .* (-x)^k for k in 0:n])
    else
        return 1
    end
end

smoothstepInterploation(a0, a1, w; n=1) = (a1 - a0) * smoothstep(w, n=n) + a0

smoothstep1(x) = 3x^2 - 2x^3
smoothstepInterploation1(a0, a1, w) = (a1 - a0) * smoothstep1(w) + a0
smoothstep2(x) = 6x^5 - 15x^4 + 10x^3
smoothstepInterploation2(a0, a1, w) = (a1 - a0) * smoothstep2(w) + a0

gridCorners(grid, cell) = grid[cell[1]:cell[1]+1, cell[2]:cell[2]+1]

function plotVector!(v0, v; color=1, lw=1, scale=1)
    x0, y0 = v0
    x, y = scale*v
    plot!([x0, x0 + x], [y0, y0 + y], color=color, lw=lw, arrow=true, label=false)
end

function plotGridVectors!(gridVectors; color=1, lw=1, scale=1)
    n, m = size(gridVectors)
    for i in 1:n, j in 1:m
        plotVector!([i-1, j-1], gridVectors[i, j], color=color, lw=lw, scale=scale)
    end
    plot!()
end

function perlinNoise(x, y, gridVectors; interpolationFunction = smoothstepInterploation2)
    i, j = floor(Int64, x), floor(Int64, y)
    dx, dy = x - i, y - j
    g00 = gridVectors[i+1, j+1]
    g01 = gridVectors[i+1, j+2]
    g10 = gridVectors[i+2, j+1]
    g11 = gridVectors[i+2, j+2]

    n00 = g00 ⋅ [dx, dy]
    n01 = g01 ⋅ [dx, dy-1]
    n10 = g10 ⋅ [dx-1, dy]
    n11 = g11 ⋅ [dx-1, dy-1]

    nx1 = interpolationFunction(n00, n10, dx)
    nx2 = interpolationFunction(n01, n11, dx)

    return interpolationFunction(nx1, nx2, dy)
end

function generatePerlinNoise(gridVectors, N; interpolationFunction=smoothstepInterploation2)
    n, m = size(gridVectors)
    xRange = LinRange(0, n-1.01, N)
    yRange = LinRange(0, m-1.01, N)
    noise = Matrix{Float64}(undef, N, N)
    for (i, x) in enumerate(xRange), (j, y) in enumerate(yRange)
        noise[i, j] = perlinNoise(x, y, gridVectors, interpolationFunction=interpolationFunction)
    end
    return noise
end

##
function perlinProds(x, y, plotGridVectors)
    i, j = floor(Int64, x), floor(Int64, y)
    dx, dy = x - i, y - j
    g00 = gridVectors[i+1, j+1]
    g01 = gridVectors[i+1, j+2]
    g10 = gridVectors[i+2, j+1]
    g11 = gridVectors[i+2, j+2]

    n00 = g00 ⋅ [dx, dy]
    n01 = g01 ⋅ [dx, dy-1]
    n10 = g10 ⋅ [dx-1, dy]
    n11 = g11 ⋅ [dx-1, dy-1]
    return n00, n01, n10, n11
end

function generatePerlinProds(gridVectors, N)
    n, m = size(gridVectors)
    xRange = LinRange(0, n-1.01, N)
    yRange = LinRange(0, m-1.01, N)
    noise00 = Matrix{Float64}(undef, N, N)
    noise01 = Matrix{Float64}(undef, N, N)
    noise10 = Matrix{Float64}(undef, N, N)
    noise11 = Matrix{Float64}(undef, N, N)
    for (i, x) in enumerate(xRange), (j, y) in enumerate(yRange)
        noise00[i, j], noise01[i, j], noise10[i, j], noise11[i, j] = perlinProds(x, y, gridVectors)
    end
    return noise00, noise01, noise10, noise11
end

