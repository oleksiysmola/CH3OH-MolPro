using Printf

stretchesSpacing::Vector{Float64} = [0.0000, 0.035, -0.035, 0.04, -0.04, 0.045, -0.045, 0.05, -0.05, 0.055, -0.055, 0.06, -0.06]
angleSpacing::Vector{Float64} = [0.0000, 1.7500, -1.7500, 2.000, 2.000, 2.25, -2.25, 2.50, -2.50, 3.75, -3.75, 4.00, -4.00]
dihedralSpacing::Vector{Float64} = [0.0000, 1.7500, -1.7500, 2.000, 2.000, 2.25, -2.25, 2.50, -2.50, 3.75, -3.75, 4.00, -4.00]
torsionSpacing::Vector{Float64} = [0.0000, 2.100, 2.200, 2.30, 2.40, 2.5, 2.6, 2.7, 2.8, 2.9, 3.00, 3.2, 3.4, 3.6, 3.8, 3.00]

stretchesGrid::Int64 = size(stretchesSpacing)[1]
angleGrid::Int64 = size(angleSpacing)[1]
dihedralGrid::Int64 = size(angleSpacing)[1]
torsionGrid::Int64 = size(torsionSpacing)[1]
convertToRadians::Float64 = 2*pi/360

function EqCH(tau::Float64)::Float64
    tau = tau*convertToRadians
    a0::Float64 = 1.08902324e+00
    a1::Float64 = 2.71787366e-03
    a2::Float64 = -2.16000169e-03
    a3::Float64 = -2.55337930e-04
    rCH::Float64 = a0 + a1*cos(tau) + a2*cos(2*tau) + a3*cos(3*tau)
    return rCH
end

function EqaHCO(tau::Float64)::Float64
    tau = tau*convertToRadians
    a0::Float64 = 1.10243240e+02 
    a1::Float64 = 2.29302542e+00
    a2::Float64 = -7.16664340e-01
    a3::Float64 = 7.93169192e-02
    a4::Float64 = 3.31054618e-02
    rHCO::Float64 = a0 + a1*cos(tau) + a2*cos(2*tau) + a3*cos(3*tau) + a4*cos(4*tau)
    return rHCO
end

d1::Float64 =  61.43364279
d2::Float64 = 180.00000000
d3::Float64 = 298.56635721
SaEq::Float64 = -1.7558466544600673
SbEq::Float64 =  3.04121561582463
tau::Float64 = 60.00

ahh1::Float64 = tau+1.0/3.0*sqrt(2.0)*SbEq
ahh2::Float64 = 120.0+tau-1.0/6.0*sqrt(2.0)*SbEq-1.0/6.0*sqrt(6.0)*SaEq
ahh3::Float64 = 240.0+tau-1.0/6.0*sqrt(2.0)*SbEq+1.0/6.0*sqrt(6.0)*SaEq

rCOeq::Float64 =                1.42077677
rOHeq::Float64=                 0.96013932
aCOHeq::Float64=              108.12930637
aHH1eq::Float64=               61.43364279
aHH2eq::Float64=              180.00000000
aHH3eq::Float64=              298.56635721
tauEq::Float64 = 60.00000
rCH1eq::Float64=             EqCH(tauEq)         #1.09108970
rCH2eq::Float64=             EqCH(tauEq + 120)   #1.08555104
rCH3eq::Float64=             EqCH(tauEq + 240)    #1.09108970
aOCH1eq::Float64=            EqaHCO(tauEq)        #111.95221297
aOCH2eq::Float64=            EqaHCO(tauEq + 120)  #106.58134561
aOCH3eq::Float64=            EqaHCO(tauEq + 240)  #111.95221297
equilibriumGrid::Vector{Float64} = [rCOeq, rOHeq, rCH1eq, rCH2eq, rCH3eq, aCOHeq, aOCH1eq, aOCH2eq, aOCH3eq, aHH1eq, aHH2eq, aHH3eq]

function PrintGeometry(point::Int64, grid::Vector{Float64})
    @printf("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", point, grid[1], grid[2], grid[3], grid[4], grid[5], grid[6], grid[7], grid[8], grid[9], grid[10], grid[11], grid[12])
end

function SubmitJob(point::Int64, grid::Vector{Float64})
    submission::Cmd = `qsub -e CH3OH_2D_MEP_AdaptiveGrid_$(point).e -o CH3OH_2D_MEP_AdaptiveGrid_$(point).o -l h_rt="11:59:00" GenerateMolproScript2D.csh $(point) $(grid[1]) $(grid[2]) $(grid[3]) $(grid[4]) $(grid[5]) $(grid[6]) $(grid[7]) $(grid[8]) $(grid[9]) $(grid[10]) $(grid[11]) $(grid[12])`
    run(submission)
end

point::Int64 = 125
grids::Vector{Vector{Float64}} = []

for i in 1:5
    for k in i+1:5
        for j in 2:stretchesGrid
            displacementVector::Vector{Float64} = zeros(12)
            displacementVector[i] = stretchesSpacing[j]
            for l in 2:stretchesGrid
                displacementVector[k] = stretchesSpacing[l]
                grid::Vector{Float64} = equilibriumGrid + displacementVector
                push!(grids, grid)
                PrintGeometry(point, grid)
                global point = point + 1
            end
        end
    end
    for k in 6:9
        for j in 2:stretchesGrid
            displacementVector::Vector{Float64} = zeros(12)
            displacementVector[i] = stretchesSpacing[j]
            for l in 2:angleGrid
                displacementVector[k] = angleSpacing[l]
                grid::Vector{Float64} = equilibriumGrid + displacementVector
                push!(grids, grid)
                PrintGeometry(point, grid)
                global point = point + 1
            end
        end
    end
    for k in 10:11
        for j in 2:stretchesGrid
            displacementVector::Vector{Float64} = zeros(12)
            displacementVector[i] = stretchesSpacing[j]
            for l in 2:dihedralGrid
                symmeterisedDihedrals::Vector{Float64} = [-1.7558466544600673, 3.04121561582463]
                symmeterisedDihedrals[k - 9] = symmeterisedDihedrals[k - 9] + dihedralSpacing[l]
                grid::Vector{Float64} = equilibriumGrid + displacementVector
                grid[10] = tauEq+1.0/3.0*sqrt(2.0)*symmeterisedDihedrals[2]
                grid[11] = 120.0+tauEq-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]-1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
                grid[12] = 240.0+tauEq-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]+1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
                PrintGeometry(point, grid)
                push!(grids, grid)
                global point = point + 1
            end
        end
    end
    for j in 2:stretchesGrid
        displacementVector::Vector{Float64} = zeros(12)
        displacementVector[i] = stretchesSpacing[j]
        for l in 2:torsionGrid
            displacementVector[10] = torsionSpacing[l]
            displacementVector[11] = torsionSpacing[l]
            displacementVector[12] = torsionSpacing[l]
            grid::Vector{Float64} = equilibriumGrid + displacementVector
            grid[3] = EqCH(tauEq + torsionSpacing[l])
            grid[4] = EqCH(tauEq + torsionSpacing[l] + 120)
            grid[5] = EqCH(tauEq + torsionSpacing[l] + 240)
            grid[7] = EqaHCO(tauEq + torsionSpacing[l])
            grid[8] = EqaHCO(tauEq + torsionSpacing[l] + 120)
            grid[9] = EqaHCO(tauEq + torsionSpacing[l] + 240)
            PrintGeometry(point, grid)
            push!(grids, grid)
            global point = point + 1
        end
    end
end

for i in 6:9
    for k in i+1:9
        for j in 2:angleGrid
            displacementVector::Vector{Float64} = zeros(12)
            displacementVector[i] = angleSpacing[j]
            for l in 2:angleGrid
                displacementVector[k] = angleSpacing[l]
                grid::Vector{Float64} = equilibriumGrid + displacementVector
                push!(grids, grid)
                PrintGeometry(point, grid)
                global point = point + 1
            end
        end
    end
    for k in 10:11
        for j in 2:angleGrid
            displacementVector::Vector{Float64} = zeros(12)
            displacementVector[i] = angleSpacing[j]
            for l in 2:dihedralGrid
                symmeterisedDihedrals::Vector{Float64} = [-1.7558466544600673, 3.04121561582463]
                symmeterisedDihedrals[k - 9] = symmeterisedDihedrals[k - 9] + dihedralSpacing[l]
                grid::Vector{Float64} = equilibriumGrid + displacementVector
                grid[10] = tauEq+1.0/3.0*sqrt(2.0)*symmeterisedDihedrals[2]
                grid[11] = 120.0+tauEq-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]-1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
                grid[12] = 240.0+tauEq-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]+1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
                PrintGeometry(point, grid)
                push!(grids, grid)
                global point = point + 1
            end
        end
    end
    for j in 2:angleGrid
        displacementVector::Vector{Float64} = zeros(12)
        displacementVector[i] = angleSpacing[j]
        for l in 2:torsionGrid
            displacementVector[10] = torsionSpacing[l]
            displacementVector[11] = torsionSpacing[l]
            displacementVector[12] = torsionSpacing[l]
            grid::Vector{Float64} = equilibriumGrid + displacementVector
            grid[3] = EqCH(tauEq + torsionSpacing[l])
            grid[4] = EqCH(tauEq + torsionSpacing[l] + 120)
            grid[5] = EqCH(tauEq + torsionSpacing[l] + 240)
            grid[7] = EqaHCO(tauEq + torsionSpacing[l])
            grid[8] = EqaHCO(tauEq + torsionSpacing[l] + 120)
            grid[9] = EqaHCO(tauEq + torsionSpacing[l] + 240)
            PrintGeometry(point, grid)
            push!(grids, grid)
            global point = point + 1
        end
    end
end

for i in 10:11
    for k in i+1:11
        for j in 2:dihedralGrid
            displacementVector::Vector{Float64} = zeros(12)
            symmeterisedDihedrals::Vector{Float64} = [-1.7558466544600673, 3.04121561582463]
            symmeterisedDihedrals[i - 9] = symmeterisedDihedrals[i - 9] + dihedralSpacing[j]
            for l in 2:dihedralGrid
                symmeterisedDihedrals[k - 9] = symmeterisedDihedrals[k - 9] + dihedralSpacing[l]
                grid::Vector{Float64} = equilibriumGrid + displacementVector
                grid[10] = tauEq+1.0/3.0*sqrt(2.0)*symmeterisedDihedrals[2]
                grid[11] = 120.0+tauEq-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]-1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
                grid[12] = 240.0+tauEq-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]+1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
                PrintGeometry(point, grid)
                push!(grids, grid)
                global point = point + 1
            end
        end
    end
    for j in 2:dihedralGrid
        displacementVector::Vector{Float64} = zeros(12)
        symmeterisedDihedrals::Vector{Float64} = [-1.7558466544600673, 3.04121561582463]
        symmeterisedDihedrals[i - 9] = symmeterisedDihedrals[i - 9] + dihedralSpacing[j]
        for l in 2:torsionGrid
            grid::Vector{Float64} = equilibriumGrid + displacementVector
            grid[3] = EqCH(tauEq + torsionSpacing[l])
            grid[4] = EqCH(tauEq + torsionSpacing[l] + 120)
            grid[5] = EqCH(tauEq + torsionSpacing[l] + 240)
            grid[7] = EqaHCO(tauEq + torsionSpacing[l])
            grid[8] = EqaHCO(tauEq + torsionSpacing[l] + 120)
            grid[9] = EqaHCO(tauEq + torsionSpacing[l] + 240)
            grid[10] = tauEq+torsionSpacing[l]+1.0/3.0*sqrt(2.0)*symmeterisedDihedrals[2]
            grid[11] = 120.0+tauEq+torsionSpacing[l]-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]-1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
            grid[12] = 240.0+tauEq+torsionSpacing[l]-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]+1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
            PrintGeometry(point, grid)
            push!(grids, grid)
            global point = point + 1
        end
    end
end