using Plots, Statistics

# Structs that are geometric figures

mutable struct Metatile
    type::Char
    a::Float64 #
    b::Float64 # parameters of monotiles
    c::Float64 #
    translation::Vector{Float64}
    rotation::Float64
    scale::Float64
end

mutable struct HatTile
    translation::Vector{Float64}
    rotation::Float64
    scale::Float64
    mirrored::Bool
    parentType::Char
end

mutable struct KiteTile
    translation::Vector{Float64}
    rotation::Float64
    scale::Float64
    parent::Int64
    id::Int64
    mirrored::Bool
end

# Function for rotation
function rotatePoints(points, angle)
    R = [cos(angle) -sin(angle);
         sin(angle)  cos(angle)]
    return R*points
end

# Function for mirroring
function mirrorPoints(points)
    T = [-1 0;
          0 1]
    return T*points
end

# Methods of getting points of vertex of diffrent tiles

function getPoints(tile::Metatile)
    if tile.type == 'H'
        short = tile.c - tile.b
        long = tile.a

        p0 = [0, 0]
        p1 = [-1/2, sqrt(3)/2]*short
        p2 = p1 + [1/2, sqrt(3)/2]*long
        p3 = p2 + [1, 0]*short
        p4 = p3 + [1/2, -sqrt(3)/2]*long
        p5 = p4 + [-1/2, -sqrt(3)/2]*short
        return (tile.translation .+ rotatePoints([p0 p1 p2 p3 p4 p5], tile.rotation))*tile.scale
    elseif tile.type == 'T'
        side = tile.a + tile.b - tile.c

        p0 = [0, 0]
        p1 = [1/2, sqrt(3)/2]*side
        p2 = [1, 0]*side
        return (tile.translation .+ rotatePoints([p0 p1 p2], tile.rotation))*tile.scale
    elseif tile.type == 'P'
        a = tile.a
        b = tile.b

        p0 = [0, 0]
        p1 = [-1/2, sqrt(3)/2]*b
        p2 = p1 + [1, 0]*a
        p3 = [1, 0]*a
        return (tile.translation .+ rotatePoints([p0 p1 p2 p3], tile.rotation))*tile.scale
    elseif tile.type == 'F'
        a = tile.a
        b = tile.b
        c = tile.c
        dx = -b/2 + a - c
        dy = sqrt(3)/2*b
        e = sqrt(dx^2+dy^2)*sqrt(3)/3
        γ = atan(dy, -dx)

        p0 = [0, 0]
        p1 = [-1/2, sqrt(3)/2]*b
        p2 = p1 + [1, 0]*a
        p4 = [1, 0]*c
        p3 = p4 + [cos(5pi/6-γ), sin(5pi/6-γ)]*e
        return (tile.translation .+ rotatePoints([p0 p1 p2 p3 p4], tile.rotation))*tile.scale
    else
        throw("No Metatile of type $(tile.type)")
    end
end

function getPoints(tile::HatTile)
    basePoints = [0      -1/2       1        1        2        5/2        4       4       5      11/2 4       5/2 2 1;
                  0 sqrt(3)/2 sqrt(3) 2sqrt(3) 2sqrt(3) 5sqrt(3)/2 2sqrt(3) sqrt(3) sqrt(3) sqrt(3)/2 0 sqrt(3)/2 0 0]
    if tile.mirrored
        basePoints = mirrorPoints(basePoints)
    end
    return (tile.translation .+ rotatePoints(basePoints, tile.rotation))*tile.scale
end

function getPoints(tile::KiteTile)
    p0 = [0, 0]
    p1 = [0, sqrt(3)]
    p2 = [1, sqrt(3)]
    p3 = [3/2, sqrt(3)/2]
    return (tile.translation .+ rotatePoints([p0 p1 p2 p3], tile.rotation))*tile.scale
end

# Methods of plotting diffrent tiles 
                    
function plotTile!(tile::Metatile; alpha=0.4, lw=1, lc=:black)
    colorDict = Dict{Char, Symbol}('H' => :red, 'T' => :yellow, 'P' => :green, 'F' => :blue)
    points = getPoints(tile)
    plot!(points[1, :], points[2, :], seriestype=[:shape], alpha=alpha, color=colorDict[tile.type], legend=false, lw=lw, lc=lc)
end

function plotTile!(tile::Metatile, color; alpha=0.4, lw=1, lc=:black)
    points = getPoints(tile)
    plot!(points[1, :], points[2, :], seriestype=[:shape], alpha=alpha, color=color, legend=false, lw=lw, lc=lc)
end

function plotLine!(tile::Metatile, color; lw=1, lc=:black)
    points = getPoints(tile)
    plot!(points[1, :], points[2, :], seriestype=[:shape], color=color, legend=false, lw=lw, lc=lc)
end

function plotTile!(tile::HatTile, color; alpha=0.5, lw=1, lc=:black)
    points = getPoints(tile)
    plot!(points[1, :], points[2, :], seriestype=[:shape], alpha=alpha, color=color, legend=false, lw=lw, lc=lc)
end

function plotLine!(tile::HatTile; lw=1, lc=:black)
    points = getPoints(tile)
    plot!(points[1, :], points[2, :], seriestype=[:shape], color=false, legend=false, lw=lw, lc=lc)
end

function plotTile!(tile::HatTile; alpha=0.5, lw=1, lc=:black)
    colorDict = Dict{Char, Symbol}('H' => :red, 'T' => :yellow, 'P' => :green, 'F' => :blue)
    points = getPoints(tile)
    plot!(points[1, :], points[2, :], seriestype=[:shape], alpha=alpha, color=colorDict[tile.parentType], legend=false, lw=lw, lc=lc)
end

function plotTile!(tile::KiteTile, color; alpha=0.5, lw=1, lc=:black)
    points = getPoints(tile)
    plot!(points[1, :], points[2, :], seriestype=[:shape], alpha=alpha, color=color, legend=false, lw=lw, lc=lc)
end

function plotTile!(tile::KiteTile; alpha=0.5)
    points = getPoints(tile)
    plot!(points[1, :], points[2, :], seriestype=[:shape], alpha=alpha, color=tile.id, legend=false)
end

function plotLine!(tile::KiteTile; lw=1, lc=:black)
    points = getPoints(tile)
    plot!(points[1, :], points[2, :], seriestype=[:shape], color=false, legend=false, lw=lw, lc=lc)
end


# Function that calculates values of a, b, c parameters on next iteration                    
function calculateNextValues(a, b, c)
    ## part of F Metatile substitution
    tile0 = Metatile('F', a, b, c, [0, 0], pi, 1)
    
    trans1 = [1/2, -sqrt(3)/2]*c
    tile1 = Metatile('F', a, b, c, trans1, 0, 1)
    
    keyPoint0 = getPoints(tile0)[:, 4]
    keyPoint1 = getPoints(tile1)[:, 4]

    cix, ciy = keyPoint0 - keyPoint1
    ci = sqrt(cix^2 + ciy^2)
    αi = atan(ciy, -cix)
    
    ## part of H Metatile substitution
    tile2 = Metatile('F', a, b, c, [0, 0], 0, 1)

    trans3 = [-1/2, sqrt(3)/2]*c + [-1, 0]*a + [1/2, -sqrt(3)/2]*b
    tile3 = Metatile('F', a, b, c, trans3, 2pi/3, 1)

    keyPoint2 = getPoints(tile2)[:, 4]
    keyPoint3 = getPoints(tile3)[:, 4]
    
    zix, ziy = keyPoint3 - keyPoint2
    zi = sqrt(zix^2 + ziy^2)
    βi = atan(ziy, -zix)
    bi = ci - 2sin(βi-αi)*zi/sqrt(3)
    ai = 2sin(pi/3+αi-βi)*zi/sqrt(3)
    #ai = sqrt(bi^2+zi^2-2bi*zi*cos(βi-αi))
    
    return ai, bi, ci, αi
end

# Substitution of metatiles with new metatiles
# a, b, c are previously calculated value tile is supertile of tiles based on a, b, c
function metatileSubstitute(supertile, a, b, c, αi)
    dx = -b/2 + a - c
    dy = sqrt(3)/2*b
    e = sqrt(dx^2+dy^2)*sqrt(3)/3
    γ = atan(dy, -dx)
    if supertile.type == 'H'
        defRot = supertile.rotation+αi
        defScale = supertile.scale
        aT(angle) = [cos(defRot+angle), sin(defRot+angle)]
        defTrans = supertile.translation - aT(2pi/3)*(c-b) - aT(pi)*(c-b) - aT(-2pi/3)*a - aT(3pi/2-γ)*e + aT(2pi/3-αi)*(supertile.c-supertile.b)

        tile0 = Metatile('H', a, b, c, defTrans, defRot, defScale)

        trans1 = aT(pi)*a + aT(-2pi/3)*(c-b) + aT(-pi/3)*(a+b-c)
        tile1 = Metatile('H', a, b, c, defTrans + trans1, defRot, defScale)

        trans2 = aT(0)*(a+b-c)
        tile2 = Metatile('H', a, b, c, defTrans + trans2, defRot-2pi/3, defScale)

        tile3 = Metatile('T', a, b, c, defTrans, defRot-pi/3, defScale)

        trans4 = aT(2pi/3)*(c-b)
        tile4 = Metatile('P', a, b, c, defTrans + trans4, defRot+pi/3, defScale)

        trans5 = aT(0)*a
        tile5 = Metatile('P', a, b, c, defTrans + trans5, defRot-pi/3, defScale)

        trans6 = trans1 + aT(-pi/3)*b
        tile6 = Metatile('P', a, b, c, defTrans + trans6, defRot, defScale)

        trans7 = trans4 + aT(pi/3)*a + aT(pi)*b
        tile7 = Metatile('F', a, b, c, defTrans + trans7, defRot, defScale)

        trans8 = trans5 + aT(pi/3)*c
        tile8 = Metatile('F', a, b, c, defTrans + trans8, defRot+2pi/3, defScale)

        trans9 = trans4 + aT(pi)*c
        tile9 = Metatile('F', a, b, c, defTrans + trans9, defRot-2pi/3, defScale)

        trans10 = trans6
        tile10 = Metatile('F', a, b, c, defTrans + trans10, defRot+2pi/3, defScale)

        trans11 = trans5 + aT(-pi/3)*a + aT(pi/3)*b
        tile11 = Metatile('F', a, b, c, defTrans + trans11, defRot-2pi/3, defScale)

        trans12 = trans2 + aT(-2pi/3)*a + aT(-pi/3)*c
        tile12 = Metatile('F', a, b, c, defTrans + trans12, defRot, defScale)

        return [tile0, tile1, tile2, tile3, tile4, tile5, tile6, tile7, tile8, tile9, tile10, tile11, tile12]
    elseif supertile.type == 'T'
        defRot = supertile.rotation+αi
        defScale = supertile.scale
        aT1(angle) = [cos(defRot+angle), sin(defRot+angle)]
        defTrans = supertile.translation - aT1(2pi/3)*(c-b) - aT1(pi)*c + aT1(5pi/6-γ)*e - aT1(-αi)*(supertile.c-supertile.b)

        tile0 = Metatile('H', a, b, c, defTrans, defRot, defScale)

        trans1 = aT1(0)*a + aT1(pi/3)*(c-b) + aT1(2pi/3)*a
        tile1 = Metatile('P', a, b, c, defTrans + trans1, defRot-pi/3, defScale)

        trans2 = aT1(-pi/3)*b
        tile2 = Metatile('P', a, b, c, defTrans + trans2, defRot, defScale)
        
        trans3 = aT1(2pi/3)*(c-b)
        tile3 = Metatile('P', a, b, c, defTrans + trans3, defRot+pi/3, defScale)

        trans4 = aT1(0)*a
        tile4 = Metatile('F', a, b, c, defTrans + trans4, defRot-pi/3, defScale)

        trans5 = trans1
        tile5 = Metatile('F', a, b, c, defTrans + trans5, defRot+pi/3, defScale)

        trans6 = trans3
        tile6 = Metatile('F', a, b, c, defTrans + trans6, defRot+pi, defScale)

        return [tile0, tile1, tile2, tile3, tile4, tile5, tile6]
    elseif supertile.type == 'P'
        defRot = supertile.rotation+αi
        defScale = supertile.scale
        aT2(angle) = [cos(defRot+angle), sin(defRot+angle)]
        defTrans = supertile.translation - aT2(2pi/3)*(c-b) - aT2(pi)*c + aT2(5pi/6-γ)*e

        tile0 = Metatile('H', a, b, c, defTrans, defRot, defScale)

        trans1 = aT2(2pi/3)*(c-b) + aT2(pi)*b
        tile1 = Metatile('H', a, b, c, defTrans + trans1, defRot+pi/3, defScale)

        trans2 = trans1 + aT2(pi/3)*a + aT2(2pi/3)*(c-b) + aT2(pi)*a
        tile2 = Metatile('P', a, b, c, defTrans + trans2, defRot, defScale)

        trans3 = aT2(-pi/3)*b
        tile3 = Metatile('P', a, b, c, defTrans + trans3, defRot, defScale)
        
        trans4 = aT2(2pi/3)*(c-b)
        tile4 = Metatile('P', a, b, c, defTrans + trans4, defRot+pi/3, defScale)
        
        trans5 = aT2(0)*a + aT2(pi/3)*c
        tile5 = Metatile('F', a, b, c, defTrans + trans5, defRot+2pi/3, defScale)
        
        trans6 = trans1 + aT2(pi/3)*a
        tile6 = Metatile('F', a, b, c, defTrans + trans6, defRot, defScale)

        trans7 = trans2
        tile7 = Metatile('F', a, b, c, defTrans + trans7, defRot+2pi/3, defScale)

        trans8 = aT2(0)*a
        tile8 = Metatile('F', a, b, c, defTrans + trans8, defRot-pi/3, defScale)

        trans9 = trans4
        tile9 = Metatile('F', a, b, c, defTrans + trans9, defRot+pi, defScale)

        trans10 = trans2 + aT2(-2pi/3)*c
        tile10 = Metatile('F', a, b, c, defTrans + trans10, defRot-pi/3, defScale)

        return [tile0, tile1, tile2, tile3, tile4, tile5, tile6, tile7, tile8, tile9, tile10]
    elseif supertile.type == 'F'
        defRot = supertile.rotation+αi+pi
        defScale = supertile.scale
        aT3(angle) = [cos(defRot+angle), sin(defRot+angle)]
        defTrans = supertile.translation + aT3(pi)*a + aT3(-2pi/3)*(c-b) + aT3(-pi/3)*a - aT3(5pi/6-γ)*e
        
        tile0 = Metatile('H', a, b, c, defTrans, defRot, defScale)

        trans1 = aT3(2pi/3)*(c-b) + aT3(pi)*b
        tile1 = Metatile('H', a, b, c, defTrans + trans1, defRot+pi/3, defScale)

        trans2 = aT3(-pi/3)*b
        tile2 = Metatile('P', a, b, c, defTrans + trans2, defRot, defScale)

        trans3 = aT3(2pi/3)*(c-b)
        tile3 = Metatile('P', a, b, c, defTrans + trans3, defRot+pi/3, defScale)

        trans4 = aT3(0)*a + aT3(pi/3)*c
        tile4 = Metatile('F', a, b, c, defTrans + trans4, defRot+2pi/3, defScale)
        
        trans5 = trans1 + aT3(pi/3)*a
        tile5 = Metatile('F', a, b, c, defTrans + trans5, defRot, defScale)
        
        trans6 = trans5 + aT3(2pi/3)*c
        tile6 = Metatile('F', a, b, c, defTrans + trans6, defRot+pi, defScale)

        trans7 = trans1 + aT3(pi)*(c-b) + aT3(2pi/3)*a + aT3(-2pi/3)*b
        tile7 = Metatile('F', a, b, c, defTrans + trans7, defRot+pi/3, defScale)
        
        trans8 = aT3(0)*a
        tile8 = Metatile('F', a, b, c, defTrans + trans8, defRot-pi/3, defScale)

        trans9 = trans3
        tile9 = Metatile('F', a, b, c, defTrans + trans9, defRot+pi, defScale)

        trans10 = trans7
        tile10 = Metatile('F', a, b, c, defTrans + trans10, defRot-pi/3, defScale)
        return [tile0, tile1, tile2, tile3, tile4, tile5, tile6, tile7, tile8, tile9, tile10]
    else
        throw("No Metatile of type $(tile.type)")
    end
end

# Function that generates metatile tiling
function generateHTPFTiling(n, type, translation=[0, 0], rotation=0, scale=1)
    as = Vector{Float64}(undef, n+1)
    bs = Vector{Float64}(undef, n+1)
    cs = Vector{Float64}(undef, n+1)
    αs = Vector{Float64}(undef, n+1)
    as[1], bs[1], cs[1] = 8, 4, 6
    αs[1] = 0
    for i in 2:n+1
        as[i], bs[i], cs[i], αs[i] = calculateNextValues(as[i-1], bs[i-1], cs[i-1])
    end

    supertile = Metatile(type, as[end], bs[end], cs[end], translation, -sum(αs) + rotation, scale)
    tiles = Metatile[supertile]
    for i in n+1:-1:2
        tempTiles = Metatile[]
        for tile in tiles
            push!(tempTiles, metatileSubstitute(tile, as[i-1], bs[i-1], cs[i-1], αs[i])...)
        end
    tiles = tempTiles
    end
    
    return tiles 
end


# Substitution of metatiles into hat (Einstein tiles)
function hatSubstitute(tile::Metatile)
    if tile.type == 'H'
        calcOffset0 = 4*[cos(tile.rotation), sin(tile.rotation)]
        hat0 = HatTile(tile.translation + calcOffset0, tile.rotation+2pi/3, tile.scale, false, 'H')

        calcOffset1 = 4*[cos(tile.rotation), sin(tile.rotation)]
        hat1 = HatTile(tile.translation + calcOffset1, tile.rotation, tile.scale, false, 'H')

        calcOffset2 = (1*[cos(tile.rotation), sin(tile.rotation)] + 3sqrt(3)*[cos(tile.rotation+pi/2), sin(tile.rotation+pi/2)])
        hat2 = HatTile(tile.translation + calcOffset2, tile.rotation+pi, tile.scale, true, 'H')

        calcOffset3 = (7*[cos(tile.rotation), sin(tile.rotation)] + 3sqrt(3)*[cos(tile.rotation+pi/2), sin(tile.rotation+pi/2)])
        hat3 = HatTile(tile.translation + calcOffset3, tile.rotation+2pi/3, tile.scale, false, 'H')

        hats = [hat0, hat1, hat2, hat3]
    elseif tile.type == 'T'
        calcOffset0 = 6*[cos(tile.rotation), sin(tile.rotation)]
        hat0 = HatTile(tile.translation + calcOffset0, tile.rotation+2pi/3, tile.scale, false, 'T')
        
        hats = [hat0]
    elseif tile.type == 'P'
        calcOffset0 = (6*[cos(tile.rotation), sin(tile.rotation)] + 4*[cos(tile.rotation+2pi/3), sin(tile.rotation+2pi/3)])
        hat0 = HatTile(tile.translation + calcOffset0, tile.rotation+pi, tile.scale, false, 'P')

        calcOffset1 = 8*[cos(tile.rotation), sin(tile.rotation)]
        
        hat1 = HatTile(tile.translation + calcOffset1, tile.rotation+2pi/3, tile.scale, false, 'P')
        hats = [hat0, hat1]
    elseif tile.type == 'F'
        calcOffset0 = 2*[cos(tile.rotation), sin(tile.rotation)]
        hat0 = HatTile(tile.translation + calcOffset0, tile.rotation, tile.scale, false, 'F')

        calcOffset1 = 4*[cos(tile.rotation+2pi/3), sin(tile.rotation+2pi/3)]
        hat1 = HatTile(tile.translation + calcOffset1, tile.rotation-pi/3, tile.scale, false, 'F')

        hats = [hat0, hat1]
    else
        throw("No Monotile fo type $(tile.type)")
    end
    return hats
end

# Function to generate hat tiling
function generateHatTiling(n, type, translation=[0, 0], rotation=0,  scale=1)
    tiles = generateHTPFTiling(n, type, translation, rotation, scale)
    hats = HatTile[]
    for tile in tiles
        push!(hats, hatSubstitute(tile)...)
    end
    return hats
end

# Substitution of hat tiles into kites
function kiteSubstitute(tile::HatTile, i)
    if tile.mirrored
        defRotation = tile.rotation
        defScale = tile.scale
        aT(angle) = [cos(defRotation+angle), sin(defRotation+angle)]
        defTrans = tile.translation + aT(pi/2)*sqrt(3) + aT(pi)

        kite0 = KiteTile(defTrans, defRotation+pi/3, defScale, i, 1, true)

        trans1 = aT(5pi/6)*2sqrt(3)
        kite1 = KiteTile(defTrans + trans1, defRotation-pi/3, defScale, i, 2, true)

        kite2 = KiteTile(defTrans, defRotation+2pi/3, defScale, i, 3, true)

        trans3 = trans1
        kite3 = KiteTile(defTrans + trans3, defRotation-2pi/3, defScale, i, 4, true)

        kite4 = KiteTile(defTrans, defRotation-2pi/3, defScale, i, 5, true)

        kite5 = KiteTile(defTrans, defRotation-pi, defScale, i, 6, true)

        trans6 = aT(-5pi/6)*2sqrt(3)
        kite6 = KiteTile(defTrans + trans6, defRotation, defScale, i, 7, true)

        trans7 = trans6
        kite7 = KiteTile(defTrans + trans7, defRotation+pi/3, defScale, i, 8, true)

        return [kite0, kite1, kite2, kite3, kite4, kite5, kite6, kite7]
    else
        defRotation = tile.rotation
        defScale = tile.scale
        aT1(angle) = [cos(defRotation+angle), sin(defRotation+angle)]
        defTrans = tile.translation + aT1(pi/2)*sqrt(3) + aT1(0)

        kite0 = KiteTile(defTrans, defRotation, defScale, i, 1, false)

        trans1 = aT1(pi/6)*2sqrt(3)
        kite1 = KiteTile(defTrans + trans1, defRotation+2pi/3, defScale, i, 2, false)

        kite2 = KiteTile(defTrans, defRotation-pi/3, defScale, i, 3, false)

        trans3 = trans1
        kite3 = KiteTile(defTrans + trans3, defRotation+pi, defScale, i, 4, false)

        kite4 = KiteTile(defTrans, defRotation+pi, defScale, i, 5, false)

        kite5 = KiteTile(defTrans, defRotation-2pi/3, defScale, i, 6, false)

        trans6 = aT1(-pi/6)*2sqrt(3)
        kite6 = KiteTile(defTrans + trans6, defRotation+pi/3, defScale, i, 7, false)

        trans7 = trans6
        kite7 = KiteTile(defTrans + trans7, defRotation, defScale, i, 8, false)

        return [kite0, kite1, kite2, kite3, kite4, kite5, kite6, kite7]
    end
end

# Function that generates kite tiling
function generateKiteTiling(n, type, translation=[0, 0], rotation=0,  scale=1)
    hats = generateHatTiling(n, type, translation, rotation, scale)
    kites = KiteTile[]
    for (i, hat) in enumerate(hats) 
        push!(kites, kiteSubstitute(hat, i)...)
    end
    return kites
end

# Function that generates kite tiling that is inscribed in regular triangle tiling used in Simplex Noise
generateProperTiling(n, type) = generateKiteTiling(n, type, [0, 0], pi/12, sqrt(2)/6)

# Auxilary function to skew points of kite tile
function getSkewedPoints(tile::KiteTile)
    points = getPoints(tile)
    M = [1+sqrt(3) 1-sqrt(3);
    1-sqrt(3) 1+sqrt(3)]/(2sqrt(3))
    return M^-1 * points
end

function plotSkewedTile!(tile::KiteTile)
    points = getSkewedPoints(tile)
    plot!(points[1,:], points[2,:], seriestype=[:shape], color=false, legend=false)
end

# Function that given vector of kite tiles deletes repeated tiles
function filterRepeats(kites)
    skewed = getSkewedPoints.(kites)
    x, y = skewed[1][:,1]
    i, j = floor(Int64, x), floor(Int64, y)
    trans = [x-i; y-j]

    notRepeat = Tuple{Int64, Int64, Int64}[]
    notRepeatInds = Int64[]
    for (k, s) in enumerate(skewed)
        m = mean(s .- trans, dims=2)
        i, j = floor.(Int64, m)
        rot = mod1(round(Int64, 3/pi*mod2pi(kites[k].rotation-pi/12)), 6)
        if !((i, j, rot) in notRepeat)
            push!(notRepeat, (i, j, rot))
            push!(notRepeatInds, k)
        end
    end
    return kites[notRepeatInds]
end


