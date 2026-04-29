module MimiGIVE

using Mimi

# External Components
using MimiFAIRv1_6_2 # Climate module
using MimiSSPs # SSP socioeconomic projections and emissions
using MimiRFFSPs # RFF socioeconomic projections and emissions
using MimiBRICK # Sea Level Rise
using Mimi_NAS_pH # Ocean PH
using MimiCIAM # Sea Level Rise Damages
using MimiMooreEtAlAgricultureImpacts # Agriculture Damages

# Constants
# (10/25/2021) BEA Table 1.1.9, line 1 GDP annual values as linked here: https://apps.bea.gov/iTable/iTable.cfm?reqid=19&step=3&isuri=1&select_all_years=0&nipa_table_list=13&series=a&first_year=2005&last_year=2020&scale=-99&categories=survey&thetable=
const pricelevel_2010_to_2005 = 87.504 / 96.166
const pricelevel_2005_to_2020 = 113.648 / 87.504
const pricelevel_1995_to_2005 = 87.504 / 71.823
const pricelevel_2006_to_2005 = 87.504 / 90.204
const pricelevel_2011_to_2005 = 87.504 / 98.164

# Utilites
include("utils/utils.jl")
include("utils/lsl_downscaling.jl")
include("utils/ciam_params.jl")

# Local Helper Components
include("components/Agriculture_helper_components.jl")
include("components/PerCapitaGDP.jl")
include("components/VSL.jl")
include("components/DamageAggregator.jl")
include("components/netconsumption.jl")
include("components/identity.jl")
include("components/GlobalTempNorm.jl")
include("components/OceanHeatAccumulator.jl")
include("components/GlobalSLRNorm.jl")
include("components/Damages_RegionAggregatorSum.jl")

# Local Damage Components
include("components/energy_damages.jl")
include("components/cromar_mortality_damages.jl")
include("components/dice2016R2_damages.jl")
include("components/howard_sterner_damages.jl")

# Primary API
include("main_model.jl")
include("main_mcs.jl")
include("main_ciam.jl")
include("scc.jl")

end
