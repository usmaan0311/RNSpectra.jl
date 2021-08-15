using RNSpectra
using Test
using DataFrames,DelimitedFiles,Glob

@testset "RNSpectra.jl" begin
    f=Glob.glob("*.txt")
    dx=readdlm(f[1],',',skipstart=2)
    df=DataFrame(dx,:auto)
    Abn=df[:,2]
    λ = df[:,1]
    En = 1240*λ.^-1
    t=798e-7 # nm
    α = log(10)*Abn/t
    RNSpectra.toyozawa(α,En)

end
