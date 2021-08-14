module RNSpectra
using LsqFit, QuadGK, Plots



function ElliottToyozawa(α,En,p0=[1,0.08,5.1,0.2,0.2])

    αf(E,p) = ( p[1].*sqrt( abs(p[2]) )./(E.*(2π) )).*(sum(m-> (2*p[2]/(m^3))*(1/p[4])*exp.(-((E .- (p[3] - p[2]/(m^2))).^2)/(2*(p[4]^2))),1:1000  ) .+ (1/p[5]).*quadgk(x-> exp.(-(x .- E).^2/(2*p[5].^2))./(1 .- exp.(-2π*sqrt( abs(p[2])./(abs(x .- p[3]) ))) ),p[3],Inf,rtol=1e-6)[1]  )
    αₙₒᵣₘ = α/maximum(α)
    p0=[1,0.08,5.1,0.2,0.2] # [c0 R Eg Γm Γc]
    fit=curve_fit(αf,vec(En),αₙₒᵣₘ,p0)
    cf=coef(fit)
    J=fit.jacobian
    R=fit.resid

    return vec(cf), J, R, vec(αf(En,cf))
end

function ElliottPlot(En,α,kn)
        αₙₒᵣₘ = α/maximum(α)
        plotly()
        plot()
        p1=plot(vec(En),vec(αₙₒᵣₘ)*maximum(α))
        p2=scatter!(vec(En),vec(kn*maximum(α)),dpi=1200,legend=:false)
        return p2
        
end
    

export ElliottToyozawa, ElliottPlot

end
