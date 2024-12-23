using Printf
using Random

stretchesSpacing::Vector{Float64} = [0.0000, -0.025,0.025, 0.050, -0.0500, -0.1000, 0.1000, 0.12500, -0.125000, 0.1500, -0.1500, 0.200, 
   -0.200, 0.300, -0.300, 0.400, 0.500, 0.600, 0.700]
angleSpacing::Vector{Float64} = [0.0000, 1.000, -1.000, 2.500, -2.5000, 5.0000, -5.0000, 7.500, -7.500, 10.0000, -10.0000, 20.0000, 
    -20.0000, 30.0000, -30.0000, 50.0000, -50.0000, 60.0000]
dihedralSpacing::Vector{Float64} = [0.0000, 1.000, -1.000, 2.5000, -2.5000, 5.0000, -5.0000, 7.500, -7.500, 10.0000, -10.0000, 40.0000, 
    -40.0000, 60.0000, -60.0000, 80.0000, -80.0000, 100.0000]
torsionSpacing::Vector{Float64} = [0.0000, 5.0000, 10.0000, 15.0000, 20.00000, 25.000, 30.0000, 35.0000, 40.0000, 45.0000, 50.0000, 55.0000, 60.0000]

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
    return 1.08924348 #rCH
end

function EqaHCO(tau::Float64)::Float64
    tau = tau*convertToRadians
    a0::Float64 = 1.10243240e+02 
    a1::Float64 = 2.29302542e+00
    a2::Float64 = -7.16664340e-01
    a3::Float64 = 7.93169192e-02
    a4::Float64 = 3.31054618e-02
    rHCO::Float64 = a0 + a1*cos(tau) + a2*cos(2*tau) + a3*cos(3*tau) + a4*cos(4*tau)
    return 110.16192385 # rHCO
end

d1::Float64 = 60.00#  61.43364279
d2::Float64 = 180.00# 180.00000000
d3::Float64 = 300.00# 298.56635721
SaEq::Float64 = 0 #-1.7558466544600673
SbEq::Float64 = 0 # 3.04121561582463
tau::Float64 = 60.00

ahh1::Float64 = tau+1.0/3.0*sqrt(2.0)*SbEq
ahh2::Float64 = 120.0+tau-1.0/6.0*sqrt(2.0)*SbEq-1.0/6.0*sqrt(6.0)*SaEq
ahh3::Float64 = 240.0+tau-1.0/6.0*sqrt(2.0)*SbEq+1.0/6.0*sqrt(6.0)*SaEq

rCOeq::Float64 =                1.42077677
rOHeq::Float64=                 0.96013932
aCOHeq::Float64=              108.12930637
aHH1eq::Float64=              60.000# 61.43364279
aHH2eq::Float64=              180.000#180.00000000
aHH3eq::Float64=              300.000#298.56635721
tauEq::Float64 = 60.00000
rCH1eq::Float64=             EqCH(tauEq)         #1.09108970    
rCH2eq::Float64=             EqCH(tauEq + 120)   #1.08555104    
rCH3eq::Float64=             EqCH(tauEq + 240)   #1.09108970    
aOCH1eq::Float64=            EqaHCO(tauEq)       #111.95221297  
aOCH2eq::Float64=            EqaHCO(tauEq + 120) #106.58134561  
aOCH3eq::Float64=            EqaHCO(tauEq + 240) #111.95221297  
equilibriumGrid::Vector{Float64} = [rCOeq, rOHeq, rCH1eq, rCH2eq, rCH3eq, aCOHeq, aOCH1eq, aOCH2eq, aOCH3eq, aHH1eq, aHH2eq, aHH3eq]

function PrintGeometry(point::Int64, grid::Vector{Float64})
    @printf("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", point, grid[1], grid[2], grid[3], grid[4], grid[5], grid[6], grid[7], grid[8], grid[9], grid[10], grid[11], grid[12])
end

function SubmitJob(point::Int64, grid::Vector{Float64})
    submission::Cmd = `qsub -e CH3OH_MonteCarlo_MEP_AdaptiveGrid_$(point).e -o CH3OH_MonteCarlo_MEP_AdaptiveGrid_$(point).o -l h_rt="11:59:00" GenerateMolproScript2D.csh $(point) $(grid[1]) $(grid[2]) $(grid[3]) $(grid[4]) $(grid[5]) $(grid[6]) $(grid[7]) $(grid[8]) $(grid[9]) $(grid[10]) $(grid[11]) $(grid[12])`
    run(submission)
end

point::Int64 = 37710
grids::Vector{Vector{Float64}} = []

numberOfGridPoints4D::Int64 = 7500
numberOfGridPoints5D::Int64 = 7500
numberOfGridPoints6D::Int64 = 7500
degreesOfFreedom::Int64 = 12
equilibriumProbabilities::Vector{Float64} = [0.3, 0.3, 0.3, 0.3, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0]

function GenerateMonteCarloGrid(dimension::Int64, numberOfGridPoints::Int64, maxGridRange::Float64)
    generatedPoint::Int64 = 1
    while generatedPoint <= numberOfGridPoints
        chosenCoordinates::Vector{Int64} = zeros(dimension)
        equilibriumProbabilitiesOfChosenCoordinates::Vector{Float64} = zeros(dimension)
        currentCoordinate::Int64 = 1
        while currentCoordinate <= dimension
            newCoordinate::Int64 = rand(1:degreesOfFreedom)
            if newCoordinate in chosenCoordinates
                continue
            else
                chosenCoordinates[currentCoordinate] = newCoordinate
                equilibriumProbabilitiesOfChosenCoordinates[currentCoordinate] = equilibriumProbabilities[newCoordinate]
                currentCoordinate += 1
            end
        end
        probabilitiesOfChosenCoordinates::Vector{Float64} = rand(dimension)
        if prod(abs.(probabilitiesOfChosenCoordinates .- equilibriumProbabilitiesOfChosenCoordinates)) > maxGridRange^dimension
            continue
        else
            grid::Vector{Float64} = equilibriumGrid .+ 0.0
            dihedralDisplacement::Vector{Float64} = [0.0, 0.0, 0.0]
                for j in 1:size(chosenCoordinates)[1]
                    if 1 <= chosenCoordinates[j] < 6
                        stretchDisplacement::Float64 = 0
                        if probabilitiesOfChosenCoordinates[j] < 0.3
                            # PDF p(x) = m1 x + c  (values given not normalised)
                            # c = 0.3
                            # m1 = 1
                            stretchDisplacement = -0.3+sqrt(0.3^2 - 2*(9/200 - 3*probabilitiesOfChosenCoordinates[j]/20))
                            grid[chosenCoordinates[j]] += stretchDisplacement
                        else
                            # PDF p(x) = m2 x + c  (values given not normalised)
                            # c = 0.3
                            # m2 = -3/7
                            stretchDisplacement = 7*(0.3-sqrt(0.3^2 + 6*(9/200 - 3*probabilitiesOfChosenCoordinates[j]/20)/7))/3
                            grid[chosenCoordinates[j]] += stretchDisplacement
                        end
                    elseif 6 <= chosenCoordinates[j] < 12
                        angleDisplacement::Float64 = real((25*(1- sqrt(3)*1im+(1+sqrt(3)*1im )*(-1+2*probabilitiesOfChosenCoordinates[j]+2*sqrt(Complex((-1+probabilitiesOfChosenCoordinates[j])*probabilitiesOfChosenCoordinates[j])))^(2/3)))/(-1+2*probabilitiesOfChosenCoordinates[j]+2*sqrt((Complex(-1+probabilitiesOfChosenCoordinates[j]))*probabilitiesOfChosenCoordinates[j]))^(1/3))
                        if chosenCoordinates[j] < 10
                            grid[chosenCoordinates[j]] += angleDisplacement
                        else
                            dihedralDisplacement[chosenCoordinates[j] - 9] += angleDisplacement
                        end
                    else
                        torsionDisplacement::Float64 = -log(1 - probabilitiesOfChosenCoordinates[j]*(1 - exp(-60*log(2)/100)))*100/log(2)
                        dihedralDisplacement[3] += torsionDisplacement
                    end
                end
            grid[10] = tauEq+dihedralDisplacement[3]+1.0/3.0*sqrt(2.0)*dihedralDisplacement[2]
            grid[11] = 120.0+tauEq+dihedralDisplacement[3]-1.0/6.0*sqrt(2.0)*dihedralDisplacement[2]-1.0/6.0*sqrt(6.0)*dihedralDisplacement[1]
            grid[12] = 240.0+tauEq+dihedralDisplacement[3]-1.0/6.0*sqrt(2.0)*dihedralDisplacement[2]+1.0/6.0*sqrt(6.0)*dihedralDisplacement[1]
            PrintGeometry(point, grid)
            generatedPoint += 1
            global point += 1
        end
    end
end

GenerateMonteCarloGrid(4, numberOfGridPoints4D, 0.1)
GenerateMonteCarloGrid(5, numberOfGridPoints5D, 0.5)
GenerateMonteCarloGrid(6, numberOfGridPoints6D, 0.25)