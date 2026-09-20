using Pkg

Pkg.add("DataDeps")

using DataDeps

ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

register(DataDep(
    "rffsps_v5",
    "RFF SPs version v5",
    "https://zenodo.org/record/6016583/files/rffsps_v5.7z",
    "a39b51d7552d198123b1863ea5a131533715a7a7c9ff6fad4b493594ece49497",
    post_fetch_method=unpack
))

# How many files each directory of the unpacked archive has.
#
# DataDeps checksums the archive between `fetch_method` and `post_fetch_method`, and that is
# the whole of its integrity checking: nothing validates the tree `unpack` leaves behind, and
# nothing re-validates it on later loads (oxinabox/DataDeps.jl#6, open since 2017). `resolve`
# does check that an inner path is readable -- `datadep"name/some/file"` aborts under `CI`
# when it is not -- but only for the one file named, and it does not re-fetch.
#
# So a depot cache holding a half-written copy of the deposit (a restore or a save that did
# not finish, or a job killed mid-unpack) leaves `datadep"rffsps_v5"` perfectly satisfied by
# the directory existing, and the truncated tree is then carried forward into the next run's
# cache. What a test item sees instead is a bare `SystemError: opening file
# .../pop_income/rffsp_pop_income_run_6546.feather`, thousands of lines into an activation
# that reported success.
#
# Hence this check, on the file counts of the one pinned Zenodo record the sha256 above
# admits. Deleting and re-resolving is the same remedy DataDeps offers interactively when
# `resolve` finds a file it cannot read; it just cannot reach that prompt from here.
const RFFSPS_V5_CONTENTS = Dict(
    "pop_income" => 10000,
    "death_rates" => 1000,
    "emissions" => 3,
    "sample_numbers" => 2,
    "ypc1990" => 1,
)

# Nothing, or a one-line description of the first thing that is wrong.
function rffsps_v5_defect(dir)
    for (subdir, expected) in sort!(collect(RFFSPS_V5_CONTENTS))
        path = joinpath(dir, subdir)
        isdir(path) || return "$subdir/ is missing"
        # Dotfiles are not part of the archive; macOS runners like to leave them anyway.
        found = count(f -> !startswith(f, "."), readdir(path))
        found == expected || return "$subdir/ has $found files, expected $expected"
    end
    return nothing
end

dir = datadep"rffsps_v5"
defect = rffsps_v5_defect(dir)

if defect !== nothing
    @warn "The rffsps_v5 data deposit is incomplete ($defect); deleting it and fetching again"
    rm(dir, recursive=true)
    dir = datadep"rffsps_v5"
    defect = rffsps_v5_defect(dir)
    defect === nothing || error("rffsps_v5 is still incomplete after re-fetching it: $defect")
end
