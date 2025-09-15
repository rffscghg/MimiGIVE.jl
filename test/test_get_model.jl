module TestGetModel

using MimiGIVE
using Test
using MimiMooreEtAlAgricultureImpacts

import MimiGIVE: get_model, compute_scc

# function get_model(; 
#     Agriculture_gtap::String = "midDF",
#     socioeconomics_source::Symbol = :RFF,
#     SSP_scenario::Union{Nothing, String} = nothing,       
#     RFFSPsample::Union{Nothing, Int} = nothing,
#     Agriculture_floor_on_damages::Bool = true,
#     Agriculture_ceiling_on_benefits::Bool = false,
#     vsl::Symbol= :epa
# )

##------------------------------------------------------------------------------
## API - test that run without error
##------------------------------------------------------------------------------

m = get_model()
run(m)

# RFF socioeconomics
for Agriculture_gtap in ["AgMIP_AllDF", "AgMIP_NoNDF", "highDF", "lowDF", "midDF"]
    for RFFSPsample in [1, 2]
        for Agriculture_floor_on_damages in [true, false]
            for Agriculture_ceiling_on_benefits in [true, false]
                for vsl in [:epa, :fund]
                    get_model(; Agriculture_gtap=Agriculture_gtap,
                        socioeconomics_source=:RFF,
                        RFFSPsample=RFFSPsample,
                        Agriculture_floor_on_damages=Agriculture_floor_on_damages,
                        Agriculture_ceiling_on_benefits=Agriculture_ceiling_on_benefits,
                        vsl=vsl)
                end
            end
        end
    end
end

# SSP socioeconomics
for Agriculture_gtap in ["AgMIP_AllDF", "AgMIP_NoNDF", "highDF", "lowDF", "midDF"]
    for SSP_scenario in ["SSP126", "SSP245", "SSP370", "SSP585"]
        for Agriculture_floor_on_damages in [true, false]
            for Agriculture_ceiling_on_benefits in [true, false]
                for vsl in [:epa, :fund]
                    get_model(; Agriculture_gtap=Agriculture_gtap,
                        socioeconomics_source=:SSP,
                        SSP_scenario=SSP_scenario,
                        Agriculture_floor_on_damages=Agriculture_floor_on_damages,
                        Agriculture_ceiling_on_benefits=Agriculture_ceiling_on_benefits,
                        vsl=vsl)
                end
            end
        end
    end
end

# some errors
@test_throws ErrorException get_model(; socioeconomics_source=:SSP) # missing SSP scenario
@test_throws ErrorException get_model(; socioeconomics_source=:foo) # not a legal SSP option
@test_throws ErrorException get_model(; Agriculture_gtap="foo") # not a legal gtap spec option
@test_throws ErrorException get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP8") # not a legal SSP scenario option
@test_throws ErrorException get_model(; vsl=:foo) # not a legal vsl option

##------------------------------------------------------------------------------
## keyword arguments and values
##------------------------------------------------------------------------------

# Agriculture GTAP Parameter (Agriculture_gtap)
sccs = []
agcosts = []

for Agriculture_gtap in ["AgMIP_AllDF", "AgMIP_NoNDF", "highDF", "lowDF", "midDF"]
    m = get_model(; Agriculture_gtap=Agriculture_gtap)
    run(m)

    append!(sccs, compute_scc(m, year=2020))
    append!(agcosts, sum(skipmissing(m[:Agriculture, :agcost])))
    gtap_idx = findfirst(isequal(Agriculture_gtap), MimiMooreEtAlAgricultureImpacts.gtaps)

    @test m[:Agriculture, :gtap_df] == MimiMooreEtAlAgricultureImpacts.gtap_df_all[:, :, gtap_idx]
end

@test allunique(sccs)
@test allunique(agcosts)
@test agcosts[4] > agcosts[5] > agcosts[3] # lowDF > midDF > highDF
@test sccs[4] > sccs[5] > sccs[3]  # lowDF > midDF > highDF

# socioeconomics_source and SSP_scenario and RFFSPsample
sccs = []
co2_emissions = []
gdp = []
pop = []

for id in [1, 2, 3]
    m_rff = get_model(; RFFSPsample=id)
    run(m_rff)

    append!(sccs, compute_scc(m_rff, year=2020))
    push!(co2_emissions, m_rff[:Socioeconomic, :co2_emissions])
    push!(gdp, m_rff[:Socioeconomic, :gdp_global])
    push!(pop, m_rff[:Socioeconomic, :population_global])

    @test(m_rff[:Socioeconomic, :id] == id)
end

for ssp in ["SSP126", "SSP245", "SSP370", "SSP585"]
    m_ssp = get_model(; socioeconomics_source=:SSP, SSP_scenario=ssp)
    run(m_ssp)

    append!(sccs, compute_scc(m_ssp, year=2020))
    push!(co2_emissions, m_ssp[:Socioeconomic, :co2_emissions])
    push!(gdp, m_ssp[:Socioeconomic, :gdp_global])
    push!(pop, m_ssp[:Socioeconomic, :population_global])

    @test(m_ssp[:Socioeconomic, :SSP] == ssp[1:4])
    @test(m_ssp[:Socioeconomic, :emissions_scenario] == ssp)
end

@test allunique(sccs)
for i in 1:length(gdp), j in 1:length(gdp) # equivalent to allunique for two arrays
    if i !== j
        @test gdp[i] !== gdp[j]
        @test pop[i] !== pop[j]
        @test co2_emissions[i] !== co2_emissions[j]
    end
end

@test compute_scc(get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP585"), year=2020) > compute_scc(get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP126"), year=2020)
@test compute_scc(get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP245"), year=2020) > compute_scc(get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP126"), year=2020)

# vsl
m_epa = get_model(vsl=:epa)
m_fund = get_model(vsl=:fund)
run(m_epa)
run(m_fund)
@test skipmissing(m_epa[:VSL, :vsl]) !== skipmissing(m_fund[:VSL, :vsl])

end # module
