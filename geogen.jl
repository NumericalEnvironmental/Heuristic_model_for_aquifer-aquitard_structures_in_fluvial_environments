# **********************************************************************************************
#
# GeoGen.jl - a geologically-oriented random field generator
#
# by Walt McNab
#
# Key concepts:
# (1) Interwoven, meandering stream channels create self-correlation structures
# (2) Streams deposit high-conductivity materials and cut older structures
# (3) Low-conductivity materials are deposited in overbank areas and do not cut underlying structures
# (4) Grid cells are assigned a network of other cells with which they are correlated
# (5) A virtual nearby point technique is used to generate values, which gradually become less important as grid fills
# (6) Correlated structures are by x-y layer; no z-direction vertical correlation is assumed
# (7) Streams follow a preferred direction, generally aligned with x-axis
# (8) Domain is divided into formations by faults with arbitrary strike and dip orientation
#
# **********************************************************************************************


using DelimitedFiles


#####################
#
# data structs
#
#####################


mutable struct Cell
    noticed::Bool                   # boolean flag as to whether cell has yet been processed
    x::Float64                      # coordinates of cell's center point within layer
    y::Float64
    z::Float64
    logK::Float64                   # hydraulic conductivity
    neighborList::Array{Int64, 1}   # an array of cell indices corresponding to search neighborhood
end


mutable struct Model
    xLength::Float64                # spatial extent of model
    yLength::Float64
    zLength::Float64
    nx::Int64                       # spatial extent of model
    ny::Int64
    nz::Int64
    dx::Float64                     # discretization
    dy::Float64
    dz::Float64
    aMin::Float64                   # max & min search distances for cells along a stream (i.e., reach and width factors)
    aMax::Float64
    bMin::Float64
    bMax::Float64
    downstream::Float64             # downstream (x-axis) random-walk weighting; 1 = neutral
    minStreamFrac::Float64          # min & max proportions of stream cells (seeds) at upgradient boundary
    maxStreamFrac::Float64
    minLogK::Float64                # min & max log hydraulic conductivity corresponding to stream flow intensity scale (0 to 1)
    maxLogK::Float64
    logKStdev::Float64              # universal standard deviation for (log) K
    searchLK::Float64               # search range for overbank cells
    logKMeanLK::Float64             # log mean hydraulic conductivity for overbank (low-K) sediments
    logKstdevLK::Float64            # log hydraulic conductivity standard deviation for overbank (low-K) sediments
    expF::Float64                   # inverse distance weighting exponent
    smoothF::Float64                # inverse distance weighting smoothing factor (employ a minimum to prevent divided by zero)
    distCorPt::Float64              # distance to virtual point when positing grid cell hydraulic conductivity values
    vertDrift::Float64              # random-walk drift factor to address variability in vertical hydraulic conductivity
end


mutable struct Fault
    d::Float64                      # fault thickness
    logK::Float64                   # fault hydraulic conductivity
    A::Float64                      # fault plane equation coefficients: Ax + By + Cz + D = 0
    B::Float64
    C::Float64
    D::Float64
end


#####################
#
# utility functions
#
#####################


function Distance(x1::Float64, y1::Float64, x2::Float64, y2::Float64)::Float64
    # distance between locations (x1, y1) and (x2, y2)
    return sqrt((x1-x2)^2 + (y1-y2)^2)
end


function InverseDist(x::Float64, y::Float64, xSet::Array{Float64, 1}, ySet::Array{Float64, 1}, vSet::Array{Float64, 1}, expF::Float64, smoothF::Float64)::Float64
    # interpolation of value at a point using inverse distance weighting (xSet/ySet/v are assumed to have been screened for search distance)
    denom = 0.
    numer = 0.
    for i = 1:length(xSet)
        dist = Distance(x, y, xSet[i], ySet[i])
        h = sqrt(dist + smoothF)
        denom += vSet[i]/h^expF
        numer += 1.0/h^expF
    end
    return denom/numer
end


function FindIndex(x::Float64, y::Float64, model::Model)::Int64
    # find the index number of nearest cell object associated with location (x, y)
    row = round(Int, y/model.dy + 0.5)
    col = round(Int, x/model.dx + 0.5)
    return (row-1)*model.nx + col
end


function FindCoord(index::Int64, model::Model)
    # return center point coordinates for cell[index]
    row = round(Int, index/model.nx - 0.5)
    col = remainder
    return row*(model.dy-0.5), y
end


function FindNeighbors(group1::Array{Int64, 1}, group2::Array{Int64, 1}, cell::Array{Cell, 1}, maxDist::Float64)
    # return an array of cell pairs (by index number) separated by <= maxDist
    neighbors = []
    for i = 1:length(group1)
        for j = i:length(group2)
            if (group1[i] != group2[j]) && Distance(cell[group1[i]].x, cell[group1[i]].y, cell[group2[i]].x, cell[group2[i]].y) <= maxDist
                push!(neighbors, [group1[i], group2[j]])
            end
        end
    end
    return neighbors
end


function cellLinks(cell::Array{Cell, 1}, streamCell::Array{Int64, 1}, neighbors)::Array{Cell, 1}
    # clear neighbors lists for all stream cells (newer streams overwrite older ones)
    for index in streamCell
        empty!(cell[index].neighborList)
        cell[index].noticed = false
    end
    # cross-reference all cell (pre-selected) neighbors within search distance
    for pair in neighbors
        push!(cell[pair[1]].neighborList, pair[2])
        push!(cell[pair[2]].neighborList, pair[1])
    end
    return cell
end


function BoundedWalk(lowerBound::Float64, upperBound::Float64, numSteps::Int64, stepSize::Float64)::Array{Float64, 1}
    # random walk between two boundaries
    wander = Float64[]
    point = RandBetween(lowerBound, upperBound)
    push!(wander, point)
    for i = 1:numSteps
        shift = RandBetween(-stepSize, stepSize)
        if point + shift < lowerBound
            point = 2*lowerBound - (point+shift)
        elseif point + shift > upperBound
            point = 2*upperBound - (point+shift)
        else
            point += shift
        end
        push!(wander, point)
    end
    return wander
end


function BoundedRandNormal(lowerBound::Float64, upperBound::Float64, meanVal::Float64, stdev::Float64)::Float64
    # selection of a location from a bounded normal distribution
    # avoids requiring installation of Distributions package (with truncated normal function)
    rVal = lowerBound - 1.
    while ((rVal < lowerBound) || (rVal > upperBound))
        rVal = randn()*stdev + meanVal
    end
    return rVal
end


function RandBetween(a::Float64, b::Float64)::Float64
    rVal = a + rand() * (b - a)
    return rVal
end


##############################################
#
# fault/plane geometry functions
#
##############################################

        
function NormPlane(phi::Float64, theta::Float64)
    # calculate the vector normal to a plane of strike theta and dip phi in 3-space
    # see http://eqseis.geosc.psu.edu/~cammon/HTML/UsingMATLAB/PDF/ML1%20FaultNormals.pdf
    xn = sin(phi) * cos(theta)    
    yn = -sin(phi) * sin(theta)
    zn = -cos(phi)
    return xn, yn, zn
end


function WithinFault(x::Float64, y::Float64, z::Float64, fault::Fault)::Bool
    # determine if point (x, y, z) is within the fault zone
    interior = (DistPtPlane(x, y, z, fault) <= 0.5*fault.d)
    return interior         # return result as boolean flag
end

        
function DistPtPlane(x::Float64, y::Float64, z::Float64, fault::Fault)::Float64
    # normal distance between point (x, y, z) and fracture plane
    # see College Calculus with Analytic Geometry, 3rd Edition, Protter and Morrey (1982); p. 410
    num = abs(fault.A*x + fault.B*y + fault.C*z + fault.D)
    den = sqrt(fault.A^2 + fault.B^2 + fault.C^2)
    return num/den
end


function ProcessFmts(lower::Array{Cell, 1}, upper::Array{Cell, 1}, fault::Fault)::Array{Cell, 1}
    # create composite model, dividing log K distribution above/below fault plane 
    println("Processing fault ...")
    composite = Cell[]
    for (i, cUp) in enumerate(upper)
        if WithinFault(cUp.x, cUp.y, cUp.z, fault)==true                # cell is within fault plane zone itself
            cUp.logK = fault.logK
            push!(composite, cUp)
        elseif cUp.z > (-fault.A*cUp.x - fault.B*cUp.y - fault.D) / fault.C     # cell is above fault plane
            push!(composite, cUp)       
        else                                                            # cell is below fault plane
            push!(composite, lower[i])      
        end
    end
    return composite
end

        
##############################################
#
# streambed conductivity generation functions
#
##############################################


function Overbank(cell::Array{Cell, 1}, model::Model)
    # (1) identify overbank cells
    overbank = Int64[]
    for (i, c) in enumerate(cell)
        if c.noticed == false
            push!(overbank, i)
        end
    end
    # (2) find all neighbor pairs within search distance
    neighbors = FindNeighbors(overbank, overbank, cell, model.searchLK)
    # (3) update neighbor lists for cells in overbank material
    cell = cellLinks(cell, Int64[], neighbors)
    return cell, overbank
end


function Meander(model::Model, origin::Float64, density::Float64)::Array{Int64, 1}
    # create stream pathway
    xStream = Float64[]                             # arrays to hold meandering stream sample points
    yStream = Float64[]
    d = model.xLength/20.                           # displacement distance for sample point generation
    x = 0.5 * model.dx                              # stream starts along x = 0 boundary
    y = BoundedRandNormal(0., model.yLength, origin, density*model.yLength)

    # first loop creates random (x, y) sample points for stream centerline
    while (x < model.xLength)
        push!(xStream, x)
        push!(yStream, y)
        theta = -pi/2 + rand()*pi                  # pick random downstream direction
        x += d * cos(theta) * model.downstream
        y += d * sin(theta)
    end
    # second loop interpolates stream y-values at each grid cell center point
    pathway = Int64[]
    dummyY = zeros(length(xStream))     # inverse distance modeling is 1-D here ...
    for i = 1:model.nx
        xc = (i-0.5)*model.dx
        yc = InverseDist(xc, 0., xStream, dummyY, yStream, model.expF, model.smoothF)
        if ((yc > 0.) && (yc < model.yLength))
            cellIndex = FindIndex(xc, yc, model)    # stream can meander in and out of model domain
            push!(pathway, cellIndex)
        end
    end
    return pathway      # 'pathway' is an array of cell indices (not coordinates) representing stream center line
end


function StreamBed(pathway::Array{Int64, 1}, model::Model)
    # find all cells (by index) belonging to streambed delineated by 'pathway'
    streamCell = Int64[]
    searchGroup = Int64[]
    width = RandBetween(model.bMin, model.bMax)
    cellStep = round(Int64, 0.5*width/model.dy)
    meanSearch = 0.5*(model.aMin + model.aMax)
    for (i, index) in enumerate(pathway)
        groupIndex = round(Int64, floor(i*model.dx/meanSearch)) + 1
        push!(streamCell, index)
        push!(searchGroup, groupIndex)
        for j = 1:length(cellStep)              # add in cells spanning y-direction, both up and down
            cellUp = index + j*model.nx
            if cellUp <= model.nx * model.ny
                push!(streamCell, cellUp)
                push!(searchGroup, groupIndex)
            end
            cellDown = index - j*model.nx
            if cellDown > 0
                push!(streamCell, cellDown)
                push!(searchGroup, groupIndex)
            end
        end
    end
    return streamCell, searchGroup
end


function StreamNeighbors(streamCell, searchGroup, cell::Array{Cell, 1}, model::Model)
    neighbors = []
    cellGroup = []
    for i = 1:maximum(searchGroup)
        push!(cellGroup, Int64[])       # first index is the cell group number, the second will be the cell number
    end
    # delineate search groups from among streamCell array
    for (i, groupIndex) in enumerate(searchGroup)
        push!(cellGroup[groupIndex], streamCell[i])
    end
    searchLength = RandBetween(model.aMin, model.aMax)
    # search within stream cell segment groups for pairs < searchLength
    for group in cellGroup
        cellNeighbors = FindNeighbors(group, group, cell, searchLength)
        neighbors = vcat(neighbors, cellNeighbors)
    end
    # search among adjacent cell segment groups for pairs < searchLength
    if length(cellGroup) > 1
        for i in length(cellGroup)-1
            cellNeighbors = FindNeighbors(cellGroup[i], cellGroup[i+1], cell, searchLength)
            neighbors = vcat(neighbors, cellNeighbors)
        end
    end
    return neighbors        # return list of neighbors within search range
end


function AssignFieldVal(remainingCells::Array{Int64, 1}, cell::Array{Cell, 1}, model::Model, meanVal::Float64, stdevVal::Float64)
    # select a random member of remainingCells
    rIndex = round(Int64, rand()*length(remainingCells) + 0.5)
    rMember = remainingCells[rIndex]
    # find the member's 'noticed' FindNeighbors; note x, y, and value
    xCheckPoint = Float64[]
    yCheckPoint = Float64[]
    vCheckPoint = Float64[]
    for nearCell in cell[rMember].neighborList
        if cell[nearCell].noticed == true
            push!(xCheckPoint, cell[nearCell].x)
            push!(yCheckPoint, cell[nearCell].y)
            push!(vCheckPoint, cell[nearCell].logK)
        end
    end
    # posit virtual point
    rVal = randn()*stdevVal + meanVal
    push!(xCheckPoint, cell[rMember].x + model.distCorPt)
    push!(yCheckPoint, cell[rMember].y)
    push!(vCheckPoint, rVal)
    # estimate field (e.g., hydraulic conductivity) by inverse distance weighting
    cell[rMember].logK = InverseDist(cell[rMember].x, cell[rMember].y, xCheckPoint, yCheckPoint, vCheckPoint,
        model.expF, model.smoothF)
    cell[rMember].noticed = true                # mark cell as noticed
    deleteat!(remainingCells, rIndex)           # remove the random member from remainingCells list
    return cell, remainingCells
end


##################################
#
# input, output & setup functions
#
##################################


function ReadFault()::Fault
    # fault parameters
    data = readdlm("fault.txt", '\t', header=false)
    phi = Float64(data[1, 2]) * pi/180
    theta = Float64(data[2, 2]) * pi/180
    yIntercept = Float64(data[3, 2])
    d = Float64(data[4, 2])
    logK = Float64(data[5, 2])
    A, B, C = NormPlane(phi, theta)         # normal vector components/equation coefficient to fracture plane
    D = -B * yIntercept                     # complete equation for fault plane (running through y-axis intercept)
    fault = Fault(d, logK, A, B, C, D)      # see http://mathonline.wikidot.com/point-normal-form-of-a-plane
    println("Read and computed fault properties.")
    return fault
end


function ModelParams(fmtName::String)::Model
    # read various model parameters from file
    data = readdlm(fmtName * ".txt", '\t', header=false)
    xLength = Float64(data[1, 2])
    yLength = Float64(data[2, 2])
    zLength = Float64(data[3, 2])
    nx = Int64(data[4, 2])
    ny = Int64(data[5, 2])
    nz = Int64(data[6, 2])
    aMin = Float64(data[7, 2])
    aMax = Float64(data[8, 2])
    bMin = Float64(data[9, 2])
    bMax = Float64(data[10, 2])
    downstream = Float64(data[11, 2])
    minStreamFrac = Float64(data[12, 2])
    maxStreamFrac = Float64(data[13, 2])
    minLogK = Float64(data[14, 2])
    maxLogK = Float64(data[15, 2])
    logKStdev = Float64(data[16, 2])
    searchLK = Float64(data[17, 2])
    logKMeanLK = Float64(data[18, 2])
    logKstdevLK = Float64(data[19, 2])
    expF = Float64(data[20, 2])
    smoothF = Float64(data[21, 2])
    distCorPt = Float64(data[22, 2])
    vertDrift = Float64(data[23, 2])
    dx = xLength/nx
    dy = yLength/ny
    dz = zLength/nz
    model = Model(xLength, yLength, zLength, nx, ny, nz, dx, dy, dz, aMin, aMax, bMin, bMax,
        downstream, minStreamFrac, maxStreamFrac, minLogK, maxLogK, logKStdev, searchLK,
        logKMeanLK, logKstdevLK, expF, smoothF, distCorPt, vertDrift)
    println("\tRead in model parameters.")
    return model
end


function WriteCells(cell::Array{Cell, 1}, fmtName::String)
    # write property distribution to file
    fname = fmtName * ".csv"
    csvfile = open(fname,"w")
    line_out = "x" * "," * "y" * "," * "z" * "," * "log_K"
    println(csvfile, line_out)
    for ce in cell
        line_out = string(ce.x) * "," * string(ce.y) * "," * string(ce.z) * "," * string(ce.logK)
        println(csvfile,line_out)
    end
    close(csvfile)
    println("Wrote output file for " * fmtName * ".")
end


function SpawnCells(model::Model, iLayer::Int64)::Array{Cell, 1}
    # create cell objects
    cell = Cell[]
    for j = 1:model.ny, i = 1:model.nx
        x = (i-0.5) * model.dx
        y = (j-0.5) * model.dy
        z = (iLayer-0.5) * model.dz
        push!(cell, Cell(false, x, y, z, 0.0, Int64[]))    # cells will be stacked later to create layers/z-values
    end
    return cell
end


######################
#
# control functions
#
######################


function Formation(fmtName::String)::Array{Cell, 1}

    # fill out the model grid based on model parameters given in 'fmtName'
    println("Working on formation " * fmtName * " ...")
    model = ModelParams(fmtName)            # read parameters for this formation

    cellsAll = Cell[]                           #

    # assumptions regarding heterogeneity:
    # (1) stream origins in each layer associated with a focus point & associated zone
    # (2) focus point drifts with layer shift (random walk)
    # (3) intensity multiplier, and number of streams, associated with layer intensity
    # (4) layer intensity drifts with layer shift (random walk)
    # (5) layer intensity quantified by stream fraction min & max

    # stream generation across all layers
    streamDensityStep = (model.maxStreamFrac - model.minStreamFrac) / (model.vertDrift * model.nz)
    streamOriginStep = model.yLength / (0.5*model.nz)
    streamDensity = BoundedWalk(model.minStreamFrac, model.maxStreamFrac, model.nz-1, streamDensityStep)
    streamOrigin = BoundedWalk(0., model.yLength, model.nz-1, streamOriginStep)

    for k = 1:model.nz                                  # for each layer

        cell = SpawnCells(model, k)                     # create cell objects
        numStreams = round(Int64, streamDensity[k] * model.ny)
        streamIntensity = sort(rand(numStreams) * streamDensity[k])                 # used to assign K values and overwrite weaker streams

        println("\t", "Processing ", numStreams, " streams in Layer ", k)
        for level in streamIntensity
            pathway = Meander(model, streamOrigin[k], streamDensity[k])             # posit stream centerline pathway
            streamCell, searchGroup = StreamBed(pathway, model)                     # designate cells that are part of stream
            if length(searchGroup) > 0
                neighbors = StreamNeighbors(streamCell, searchGroup, cell, model)   # find stream cells within search distance
                cell = cellLinks(cell, streamCell, neighbors)                       # cross-reference cell neighbors lists
                logKMean = model.minLogK + level*(model.maxLogK - model.minLogK)
                # posit values for cells until remaining streamCell array goes to zero
                while length(streamCell) > 0
                    cell, streamCell = AssignFieldVal(streamCell, cell, model, logKMean, model.logKStdev)
                end
            end
        end

        # process overbank materials
        println("\t", "Processing overbank materials in Layer ", k)
        cell, overbankCell = Overbank(cell, model)
        while length(overbankCell) > 0
            cell, overbankCell = AssignFieldVal(overbankCell, cell, model, model.logKMeanLK, model.logKstdevLK)
        end

        # stack cells onto full 3-D array
        cellsAll = vcat(cellsAll, cell)

    end

    return cellsAll

end


function GeoGen()

    ### main routine ###

    # 1-formation format, or 2-formation, the latter separated by fault or unconformity
    fmtList = ["fluvial_1", "low_energy"]
    if length(fmtList) == 2
        # two formations; first one listed is placed above fault
        upper = Formation(fmtList[1])
        lower = Formation(fmtList[2])
        WriteCells(upper, fmtList[1])
        WriteCells(lower, fmtList[2])           
        fault = ReadFault()
        composite = ProcessFmts(lower, upper, fault)
    else
        # only a single formation; first-listed formation is the default value
        composite = Formation(fmtList[1])
    end

    WriteCells(composite, "composite")              # write output
    println("Finished.")

end


GeoGen()        ### run main routine ###

