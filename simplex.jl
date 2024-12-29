using Plots, LinearAlgebra, BenchmarkTools
include("perlin.jl")

generateGridPoints(xlim, ylim) = [[x, y] for x in 0:xlim, y in 0:ylim]

function gridLines(gridPoints; across = true)
    n, m = size(gridPoints)
    verLines = Vector{Vector{}}(undef, n*(m-1))
    for i in 1:n, j in 1:m-1
        verLines[(i-1)*(n-1) + j] = [Tuple(gridPoints[i, j]), Tuple(gridPoints[i, j+1])]
    end

    horLines = Vector{Vector{}}(undef, m*(n-1))
    for i in 1:n-1, j in 1:m
        horLines[(j-1)*(m-1) + i] = [Tuple(gridPoints[i, j]), Tuple(gridPoints[i+1, j])]
    end
    if !across
        return [verLines; horLines]    
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

function plotGrid!(gridPoints, color=1; across=true)
    lines = gridLines(gridPoints, across=across)
    for line in lines
        plot!(line, legend=false, c=color)
    end
    plot!()
end

const US = (sqrt(3)-1)/2
const S = -(3-sqrt(3))/6

skewPoint(P, A) = P .+ sum(P)*A
skewPointM(P, M) = M*P

function skewMatrix(angle)
    M = [cos(-angle) cos(pi/2+angle);
        sin(-angle) sin(pi/2+angle)]
    return M
end

function borderCorners(gridPoints, A)
    n, m = size(gridPoints)
    corners = gridPoints[[1, n, 1, n], [1, 1, m, m]]
    return skewPoint.(corners, A)
end

function safeRange(gridPoints, A)
    corners = borderCorners(gridPoints, A)
    return [corners[1][1], corners[end][1], corners[1][2], corners[end][2]]
end


function simplexNoise(x, y, gridVectors)
    i, j = floor.(Int64, skewPoint([x, y], US))
    xCell, yCell = skewPoint([i, j], S)
    dx0, dy0 = x - xCell, y - yCell
    
    di, dj = 1, 0
    if dy0 > dx0
        di, dj = 0, 1
    end

    dx1 = dx0 - di - S
    dy1 = dy0 - dj - S
    dx2 = dx0 - 1 - 2S
    dy2 = dy0 - 1 - 2S

    g0 = gridVectors[i+1, j+1]
    g1 = gridVectors[i+1+di, j+1+dj]
    g2 = gridVectors[i+2, j+2]

    t0 = 0.5 - dx0^2 - dy0^2
    t1 = 0.5 - dx1^2 - dy1^2
    t2 = 0.5 - dx2^2 - dy2^2
    if t0 < 0
        n0 = 0
    else      
        n0 = t0^4*dot(g0, [dx0, dy0])
    end
    if t1 < 0
        n1 = 0
    else      
        n1 = t1^4*dot(g1, [dx1, dy1])
    end
    if t2 < 0
        n2 = 0
    else      
        n2 = t2^4*dot(g2, [dx2, dy2])
    end
    return 100*(n0 + n1 + n2)
end

function generateSimplexNoise(size, N)
    m = ceil(Int64, skewPoint([size, size], US)[1])
    xs = 0:m
    ys = 0:m

    gridVectors = [randomVector() for x in xs, y in ys]
    
    xRange = LinRange(0, size, N)
    yRange = LinRange(0, size, N)   

    noise = Matrix{Float64}(undef, N, N)
    for (i, x) in enumerate(xRange), (j, y) in enumerate(yRange)
        noise[i, j] = simplexNoise(x, y, gridVectors)
    end
    return noise, m
end
