using Printf

stretchesSpacing::Vector{Float64} = [0.0000, 0.005, -0.005, 0.01, -0.01, 0.015, -0.015, 0.02, -0.02, 0.025, -0.025, 0.03, -0.03]
angleSpacing::Vector{Float64} = [0.0000, 0.2500, -0.2500, 0.500, -0.500, 0.75, -0.75, 1.00, -1.00, 1.25, -1.25, 1.50, -1.50]
dihedralSpacing::Vector{Float64} = [0.0000, 0.2500, -0.2500, 0.500, -0.500, 0.75, -0.75, 1.00, -1.00, 1.25, -1.25, 1.50, -1.50]
torsionSpacing::Vector{Float64} = [0.0000, 0.100, 0.200, 0.30, 0.40, 0.5, 0.6, 0.7, 0.8, 0.9, 1.00, 1.2, 1.4, 1.6, 1.8, 2.00, 3.00, 4.00, 5.00]

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

function PrintGeometry(grid::Vector{Float64})
    @printf("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", grid[1], grid[2], grid[3], grid[4], grid[5], grid[6], grid[7], grid[8], grid[9], grid[10], grid[11], grid[12])
end
PrintGeometry(equilibriumGrid)

function SubmitJob(point::Int64, grid::Vector{Float64})
    submission::Cmd = `qsub -e CH3OH_1D_MEP_AdaptiveGrid_$(point).e -o CH3OH_1D_MEP_AdaptiveGrid_$(point).o -l h_rt="11:59:00" GenerateMolproScript1D.csh $(point) $(grid[1]) $(grid[2]) $(grid[3]) $(grid[4]) $(grid[5]) $(grid[6]) $(grid[7]) $(grid[8]) $(grid[9]) $(grid[10]) $(grid[11]) $(grid[12])`
    run(submission)
end

point::Int64 = 1
SubmitJob(point, equilibriumGrid)

for i in 1:4
    for j in 2:stretchesGrid
        global point = point + 1
        displacementVector::Vector{Float64} = zeros(12)
        displacementVector[i] = stretchesSpacing[j]
        grid::Vector{Float64} = equilibriumGrid + displacementVector
        SubmitJob(point, grid)
        PrintGeometry(grid)
    end
end
for i in 6:8
    for j in 2:angleGrid
        global point = point + 1
        displacementVector::Vector{Float64} = zeros(12)
        displacementVector[i] = angleSpacing[j]
        grid::Vector{Float64} = equilibriumGrid + displacementVector
        SubmitJob(point, grid)
        PrintGeometry(grid)
    end
end
for i in 10:11
    for j in 2:dihedralGrid
        global point = point + 1
        symmeterisedDihedrals::Vector{Float64} = [-1.7558466544600673, 3.04121561582463]
        symmeterisedDihedrals[i - 9] = symmeterisedDihedrals[i - 9] + dihedralSpacing[j]
        displacementVector::Vector{Float64} = zeros(12)
        grid::Vector{Float64} = equilibriumGrid + displacementVector
        grid[10] = tau+1.0/3.0*sqrt(2.0)*symmeterisedDihedrals[2]
        grid[11] = 120.0+tau-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]-1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
        grid[12] = 240.0+tau-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]+1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
        
        SubmitJob(point, grid)
        PrintGeometry(grid)
    end
end
for j in 2:torsionGrid
    global point = point + 1
    displacementVector::Vector{Float64} = zeros(12)
    displacementVector[10] = torsionSpacing[j]
    displacementVector[11] = torsionSpacing[j]
    displacementVector[12] = torsionSpacing[j]
    grid::Vector{Float64} = equilibriumGrid + displacementVector
    grid[3] = EqCH(tauEq + torsionSpacing[j])
    grid[4] = EqCH(tauEq + torsionSpacing[j] + 120)
    grid[5] = EqCH(tauEq + torsionSpacing[j] + 240)
    grid[7] = EqaHCO(tauEq + torsionSpacing[j])
    grid[8] = EqaHCO(tauEq + torsionSpacing[j] + 120)
    grid[9] = EqaHCO(tauEq + torsionSpacing[j] + 240)
    SubmitJob(point, grid)
    PrintGeometry(grid)
end