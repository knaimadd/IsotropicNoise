include("perlin.jl")
include("simplex.jl")
include("trueMonotile.jl")

# Function that saves all needed information of tessellation into discrete data structures
function kiteMatrix(n, type)
    kites = generateProperTiling(n, type)
    skewed = getSkewedPoints.(kites)
    x, y = skewed[1][:,1]
    i, j = floor(Int64, x), floor(Int64, y)
    trans = [x-i; y-j]

    xInds = Vector{Int64}(undef, length(kites))
    yInds = Vector{Int64}(undef, length(kites))
    rots = Vector{Int64}(undef, length(kites))
    ids = Vector{Int64}(undef, length(kites))
    parents = Vector{Int64}(undef, length(kites))
    mirrors = Vector{Bool}(undef, length(kites))
    for (k, s) in enumerate(skewed)
        m = mean(s .- trans, dims=2)
        xInds[k], yInds[k] = floor.(Int64, m)
        rots[k] = mod1(round(Int64, 3/pi*mod2pi(kites[k].rotation-pi/12)), 6)
        ids[k] = kites[k].id
        parents[k] = kites[k].parent
        mirrors[k] = kites[k].mirrored
    end
    xMin, xMax = extrema(xInds)
    yMin, yMax = extrema(yInds)

    kM = Matrix{Vector{Int64}}(undef, xMax-xMin+1, yMax-yMin+1)
    pM = Matrix{Vector{Int64}}(undef, xMax-xMin+1, yMax-yMin+1)
    mM = Matrix{Vector{Bool}}(undef, xMax-xMin+1, yMax-yMin+1)
    for I in eachindex(kM)
        kM[I] = zeros(6)
        pM[I] = zeros(6)
        mM[I] = zeros(6)
    end
    for k in 1:length(xInds)
        kM[xInds[k].-xMin+1, yInds[k].-yMin+1][[5, 3, 6, 2, 4, 1][rots[k]]] = ids[k]
        pM[xInds[k].-xMin+1, yInds[k].-yMin+1][[5, 3, 6, 2, 4, 1][rots[k]]] = parents[k]
        mM[xInds[k].-xMin+1, yInds[k].-yMin+1][[5, 3, 6, 2, 4, 1][rots[k]]] = mirrors[k]
    end
    return kM, pM, mM, trans, xMin, yMin
end

# Generating random vectors table
function randomKiteGridVectors(kM)
    n, m = size(kM)
    vecs = Matrix{Vector{Vector{Float64}}}(undef, n, m)
    for i in 1:n, j in 1:m
        vecs[i, j] = [randomVector() for _ in 1:6]
    end
    return vecs
end

# Generating random vecotrs to each hat tile
function randomParentVectors(pM)
    m = findmax(x->findmax(x), pM)[1][1]
    return [randomVector() for _ in 1:m]
end

angels = [5pi/12, -pi/4, -11pi/12, pi/12, 3pi/4, -7pi/12] # table of all possible rotations of kite tile on proper tiling

# Function that evaluates kite noise at (x, y) given random vectors and inormation from tessellation  
function kiteNoise(x, y, gridVectors, trans, xMin, yMin)
    i, j = floor.(Int64, skewPoint([x, y], US) - trans)
    xCell, yCell = skewPoint([i, j] + trans, S)
    dx0, dy0 = x - xCell, y - yCell
        
    di, dj = 1, 0
    if dy0 > dx0
        di, dj = 0, 1
    end

    dx1 = dx0 - di - S
    dy1 = dy0 - dj - S
    dx2 = dx0 - 1 - 2S
    dy2 = dy0 - 1 - 2S

    k = findmin([sqrt(dx0^2 + dy0^2), sqrt(dx1^2 + dy1^2), sqrt(dx2^2 + dy2^2)])[2]
    cellInd = 3di + k
    dx0 = [dx0, dx1, dx2][k]
    dy0 = [dy0, dy1, dy2][k]
    a = angels[cellInd]
    dx1, dy1 = [dx0, dy0] - [cos(a+pi/6), sin(a+pi/6)]*sqrt(6)/6
    dx2, dy2 = [dx0, dy0] - [cos(a), sin(a)]*sqrt(2)/3
    dx3, dy3 = [dx0, dy0] - [cos(a-pi/6), sin(a-pi/6)]*sqrt(6)/6

    if cellInd == 1
        g0 = gridVectors[i - xMin + 1, j - yMin + 1][1]
        g1 = gridVectors[i - xMin + 1, j - yMin + 1][2]
        g2 = gridVectors[i - xMin + 1, j - yMin + 1][3]
        g3 = gridVectors[i - xMin + 1, j - yMin + 1][4]
    elseif cellInd == 2
        g0 = gridVectors[i - xMin + 1, j - yMin + 2][1]
        g1 = gridVectors[i - xMin + 1, j - yMin + 2][6]
        g2 = gridVectors[i - xMin + 1, j - yMin + 1][3]
        g3 = gridVectors[i - xMin + 1, j - yMin + 1][2]
    elseif cellInd == 3
        g0 = gridVectors[i - xMin + 2, j - yMin + 2][1]
        g1 = gridVectors[i - xMin + 1, j - yMin + 1][4]
        g2 = gridVectors[i - xMin + 1, j - yMin + 1][3]
        g3 = gridVectors[i - xMin + 1, j - yMin + 2][6]
    elseif cellInd == 4
        g0 = gridVectors[i - xMin + 1, j - yMin + 1][1]
        g1 = gridVectors[i - xMin + 1, j - yMin + 1][4]
        g2 = gridVectors[i - xMin + 1, j - yMin + 1][5]
        g3 = gridVectors[i - xMin + 1, j - yMin + 1][6]
    elseif cellInd == 5
        g0 = gridVectors[i - xMin + 2, j - yMin + 1][1]
        g1 = gridVectors[i - xMin + 1, j - yMin + 1][6]
        g2 = gridVectors[i - xMin + 1, j - yMin + 1][5]
        g3 = gridVectors[i - xMin + 2, j - yMin + 1][2]
    elseif cellInd == 6
        g0 = gridVectors[i - xMin + 2, j - yMin + 2][1]
        g1 = gridVectors[i - xMin + 2, j - yMin + 1][2]
        g2 = gridVectors[i - xMin + 1, j - yMin + 1][5]
        g3 = gridVectors[i - xMin + 1, j - yMin + 1][4]
    end

    t0 = (1/6 - dx0^2 - dy0^2)*6
    t1 = (2/243 - 1/18*(dx1*cos(a+pi/6) + dy1*sin(a+pi/6))^2 - 4/27*(dx1*sin(a+pi/6) - dy1*cos(a+pi/6))^2)*243/2
    t2 = (1/18 - dx2^2 - dy2^2)*18
    t3 = (2/243 - 1/18*(dx3*cos(a-pi/6) + dy3*sin(a-pi/6))^2 - 4/27*(dx3*sin(a-pi/6) - dy3*cos(a-pi/6))^2)*243/2
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
    if t3 < 0
        n3 = 0
    else
        n3 = t3^4*dot(g3, [dx3, dy3])
    end
    return 8*(n0 + n1 + n2 + n3)
end

# Auxilary function that generates kite noise on field
function generateKiteNoise(M, N, gridVectors, trans, xMin, yMin)
    xs = LinRange(0, M, N)
    ys = LinRange(0, M, N)

    noise = Matrix{Float64}(undef, N, N)
    for (i, x) in enumerate(xs), (j, y) in enumerate(ys)
        noise[i, j] = kiteNoise(x, y, gridVectors, trans, xMin, yMin)
    end
    return noise 
end

# Function that evaluates index of kite tile at point (x, y) 
function hatHeights(x, y, kM, trans, xMin, yMin)
    i, j = floor.(Int64, skewPoint([x, y], US) - trans)
    xCell, yCell = skewPoint([i, j] + trans, S)
    dx0, dy0 = x - xCell, y - yCell
        
    di, dj = 1, 0
    if dy0 > dx0
        di, dj = 0, 1
    end

    dx1 = dx0 - di - S
    dy1 = dy0 - dj - S
    dx2 = dx0 - 1 - 2S
    dy2 = dy0 - 1 - 2S

    k = findmin([sqrt(dx0^2 + dy0^2), sqrt(dx1^2 + dy1^2), sqrt(dx2^2 + dy2^2)])[2]
    cellInd = 3di + k
    return kM[i - xMin + 1, j - yMin + 1][cellInd]
end

# Function that evaluates hat noise at point (x, y) given random vecotrs table and information about tessellation
function hatNoise(x, y, gridVectors, kM, pM, mM, pVectors, trans, xMin, yMin)
    i, j = floor.(Int64, skewPoint([x, y], US) - trans)
    xCell, yCell = skewPoint([i, j] + trans, S)
    dx0, dy0 = x - xCell, y - yCell
        
    di, dj = 1, 0
    if dy0 > dx0
        di, dj = 0, 1
    end

    dx1 = dx0 - di - S
    dy1 = dy0 - dj - S
    dx2 = dx0 - 1 - 2S
    dy2 = dy0 - 1 - 2S

    k = findmin([sqrt(dx0^2 + dy0^2), sqrt(dx1^2 + dy1^2), sqrt(dx2^2 + dy2^2)])[2]
    cellInd = 3di + k
    dx0 = [dx0, dx1, dx2][k]
    dy0 = [dy0, dy1, dy2][k]
    a = angels[cellInd]
    dx1, dy1 = [dx0, dy0] - [cos(a+pi/6), sin(a+pi/6)]*sqrt(6)/6
    dx2, dy2 = [dx0, dy0] - [cos(a), sin(a)]*sqrt(2)/3
    dx3, dy3 = [dx0, dy0] - [cos(a-pi/6), sin(a-pi/6)]*sqrt(6)/6

    kInd = kM[i - xMin + 1, j - yMin + 1][cellInd]
    mInd = mM[i - xMin + 1, j - yMin + 1][cellInd]

    if (kInd in [1, 4, 6, 7] && !mInd) || (kInd in [2, 3, 5, 8] && mInd)
        dxp, dyp = [dx0, dy0] - [cos(a-pi/6), sin(a-pi/6)]*sqrt(2)/6
        tp = 1 - (dxp*cos(a-2pi/3) + dyp*sin(a-2pi/3))^2/1.3225*18 - (dxp*sin(a-2pi/3) - dyp*cos(a-2pi/3))^2/(4-2sqrt(3))*18
    else
        dxp, dyp =  [dx0, dy0] - [cos(a+pi/6), sin(a+pi/6)]*sqrt(2)/6
        tp = 1 - (dxp*cos(a+2pi/3) + dyp*sin(a+2pi/3))^2/1.3225*18 - (dxp*sin(a+2pi/3) - dyp*cos(a+2pi/3))^2/(4-2sqrt(3))*18
    end
    if cellInd == 1
        g0 = gridVectors[i - xMin + 1, j - yMin + 1][1]
        g1 = gridVectors[i - xMin + 1, j - yMin + 1][2]
        g2 = gridVectors[i - xMin + 1, j - yMin + 1][3]
        g3 = gridVectors[i - xMin + 1, j - yMin + 1][4]
    elseif cellInd == 2
        g0 = gridVectors[i - xMin + 1, j - yMin + 2][1]
        g1 = gridVectors[i - xMin + 1, j - yMin + 2][6]
        g2 = gridVectors[i - xMin + 1, j - yMin + 1][3]
        g3 = gridVectors[i - xMin + 1, j - yMin + 1][2]
    elseif cellInd == 3
        g0 = gridVectors[i - xMin + 2, j - yMin + 2][1]
        g1 = gridVectors[i - xMin + 1, j - yMin + 1][4]
        g2 = gridVectors[i - xMin + 1, j - yMin + 1][3]
        g3 = gridVectors[i - xMin + 1, j - yMin + 2][6]
    elseif cellInd == 4
        g0 = gridVectors[i - xMin + 1, j - yMin + 1][1]
        g1 = gridVectors[i - xMin + 1, j - yMin + 1][4]
        g2 = gridVectors[i - xMin + 1, j - yMin + 1][5]
        g3 = gridVectors[i - xMin + 1, j - yMin + 1][6]
    elseif cellInd == 5
        g0 = gridVectors[i - xMin + 2, j - yMin + 1][1]
        g1 = gridVectors[i - xMin + 1, j - yMin + 1][6]
        g2 = gridVectors[i - xMin + 1, j - yMin + 1][5]
        g3 = gridVectors[i - xMin + 2, j - yMin + 1][2]
    elseif cellInd == 6
        g0 = gridVectors[i - xMin + 2, j - yMin + 2][1]
        g1 = gridVectors[i - xMin + 2, j - yMin + 1][2]
        g2 = gridVectors[i - xMin + 1, j - yMin + 1][5]
        g3 = gridVectors[i - xMin + 1, j - yMin + 1][4]
    end
    pInd = pM[i - xMin + 1, j - yMin + 1][cellInd]
    gp = pVectors[pInd]

    t0 = (1/6 - dx0^2 - dy0^2)*6
    t1 = (2/243 - 1/18*(dx1*cos(a+pi/6) + dy1*sin(a+pi/6))^2 - 4/27*(dx1*sin(a+pi/6) - dy1*cos(a+pi/6))^2)*243/2
    t2 = (1/18 - dx2^2 - dy2^2)*18
    t3 = (2/243 - 1/18*(dx3*cos(a-pi/6) + dy3*sin(a-pi/6))^2 - 4/27*(dx3*sin(a-pi/6) - dy3*cos(a-pi/6))^2)*243/2
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
    if t3 < 0
        n3 = 0
    else
        n3 = t3^4*dot(g3, [dx3, dy3])
    end
    if tp < 0
        np = 0
    else
        np = tp^4*dot(gp, [dxp, dyp])
    end
    return 7*(n0 + n1 + n2 + n3 + np)
    return xCell, yCell
end

# Auxilary function for generating hat noise on field
function generateHatNoise(M, N, gridVectors, kM, pM, mM, pVectors, trans, xMin, yMin)
    xs = LinRange(0, M, N)
    ys = LinRange(0, M, N)

    noise = Matrix{Float64}(undef, N, N)
    for (i, x) in enumerate(xs), (j, y) in enumerate(ys)
        noise[i, j] = hatNoise(x, y, gridVectors, kM, pM, mM, pVectors, trans, xMin, yMin)
    end
    return noise 
end


