using CSV
using DataFrames
using Dates

# ---------------------------
# Load source datasets
# ---------------------------
adsl = CSV.read("adsl.csv", DataFrame)
adrs = CSV.read("adrs_onco.csv", DataFrame)

# If any date columns are strings, uncomment and adjust these:
# for col in [:RANDDT, :LSTALVDT]  # ADSL dates
#     adsl[!, col] = Date.(adsl[!, col])
# end
# for col in [:ADT]  # ADRS dates
#     adrs[!, col] = Date.(adrs[!, col])
# end

# ---------------------------
# Event/Censor source constructors (store filter as a Function)
# ---------------------------
struct SourceDef
    type::Symbol              # :event or :censor
    dataset_name::Symbol      # :adsl or :adrs
    filter::Union{Nothing,Function}
    date::Symbol
    set_values_to::Dict{Symbol,Any}
end

event_source(; dataset_name, filter=nothing, date, set_values_to=Dict()) =
    SourceDef(:event, Symbol(dataset_name), filter, date, set_values_to)

censor_source(; dataset_name, filter=nothing, date, set_values_to=Dict()) =
    SourceDef(:censor, Symbol(dataset_name), filter, date, set_values_to)

# ---------------------------
# Define events and censoring sources (filters as row -> Bool)
# ---------------------------
death_event = event_source(
    dataset_name = "adrs",
    filter = r -> r.PARAMCD == "DEATH" && r.AVALC == "Y" && r.ANL01FL == "Y",
    date = :ADT,
    set_values_to = Dict(:EVNTDESC => "Death", :SRCDOM => "ADRS", :SRCVAR => "ADT")
)

pd_event = event_source(
    dataset_name = "adrs",
    filter = r -> r.PARAMCD == "PD" && r.ANL01FL == "Y",
    date = :ADT,
    set_values_to = Dict(:EVNTDESC => "Progressive Disease", :SRCDOM => "ADRS", :SRCVAR => "ADT")
)

lastalive_censor = censor_source(
    dataset_name = "adsl",
    # no filter
    date = :LSTALVDT,
    set_values_to = Dict(
        :EVNTDESC => "Last Known Alive",
        :CNSDTDSC => "Last Known Alive Date",
        :SRCDOM   => "ADSL",
        :SRCVAR   => "LSTALVDT"
    )
)

lasta_censor = censor_source(
    dataset_name = "adrs",
    filter = r -> r.PARAMCD == "LSTA" && r.ANL01FL == "Y",
    date = :ADT,
    set_values_to = Dict(
        :EVNTDESC => "Progression Free Alive",
        :CNSDTDSC => "Last Tumor Assessment",
        :SRCDOM   => "ADRS",
        :SRCVAR   => "ADT"
    )
)

rand_censor = censor_source(
    dataset_name = "adsl",
    # no filter
    date = :RANDDT,
    set_values_to = Dict(
        :EVNTDESC => "Randomization Date",
        :CNSDTDSC => "Randomization Date",
        :SRCDOM   => "ADSL",
        :SRCVAR   => "RANDDT"
    )
)

# ---------------------------
# Helper: index source datasets by USUBJID for faster lookup
# ---------------------------
function index_by_usubjid(df::DataFrame)
    g = groupby(df, :USUBJID)
    Dict(key => d for (key, d) in pairs(g))
end

idx_adsl = index_by_usubjid(adsl)
idx_adrs = index_by_usubjid(adrs)

# ---------------------------
# Core TTE derivation matching {admiral} logic:
#   - Collect ALL candidate events and censors per subject
#   - Pick the earliest date overall
#   - If earliest is event => CNSR=0 else CNSR=1
# ---------------------------
function derive_param_tte(; 
    dataset_adsl::DataFrame,
    start_date::Symbol,
    event_conditions::Vector{SourceDef},
    censor_conditions::Vector{SourceDef},
    source_datasets::Dict{Symbol,DataFrame},
    set_values_to::Dict{Symbol,Any}
)
    out_cols = [:STUDYID, :USUBJID, :ADT, :EVNTDESC, :SRCDOM, :SRCVAR, :CNSR, :CNSDTDSC, :STARTDT, :PARAMCD, :PARAM]
    results = DataFrame([c => Any[] for c in out_cols]...)

    # dataset indexes by USUBJID for quick access
    adsl_idx = index_by_usubjid(source_datasets[:adsl])
    has_adrs = haskey(source_datasets, :adrs)
    adrs_idx = has_adrs ? index_by_usubjid(source_datasets[:adrs]) : Dict{Any,DataFrame}()

    for row in eachrow(dataset_adsl)
        usubjid = row.USUBJID
        studyid = row.STUDYID
        startdt  = row[start_date]

        candidates = DataFrame(
            STUDYID = String[], USUBJID = String[], ADT = Date[],
            EVNTDESC = String[], SRCDOM = String[], SRCVAR = String[],
            CNSR = Int[], CNSDTDSC = Union{Missing,String}[], STARTDT = Date[],
            PARAMCD = String[], PARAM = String[],
            _TYPE = Symbol[] # :event or :censor, internal
        )

        # Helper to add a row
        function push_candidate!(dt::Date, desc::String, srcdom::String, srcvar::String, cnsr::Int, cnsdesc::Union{Missing,String})
            push!(candidates, (
                string(studyid), string(usubjid), dt,
                desc, srcdom, srcvar,
                cnsr, cnsdesc, startdt,
                String(set_values_to[:PARAMCD]), String(set_values_to[:PARAM]),
                cnsr == 0 ? :event : :censor
            ))
        end

        # Collect events
        for e in event_conditions
            ds_idx = e.dataset_name == :adsl ? adsl_idx : adrs_idx
            if haskey(ds_idx, usubjid)
                dfsub = ds_idx[usubjid]
                # apply filter if present
                if e.filter !== nothing
                    dfsub = filter(e.filter, dfsub)
                end
                for r in eachrow(dfsub)
                    dt = r[e.date]
                    if !ismissing(dt)
                        push_candidate!(dt, e.set_values_to[:EVNTDESC], e.set_values_to[:SRCDOM], e.set_values_to[:SRCVAR], 0, missing)
                    end
                end
            end
        end

        # Collect censors
        for c in censor_conditions
            ds_idx = c.dataset_name == :adsl ? adsl_idx : adrs_idx
            if haskey(ds_idx, usubjid)
                dfsub = ds_idx[usubjid]
                if c.filter !== nothing
                    dfsub = filter(c.filter, dfsub)
                end
                for r in eachrow(dfsub)
                    dt = r[c.date]
                    if !ismissing(dt)
                        cnsdesc = get(c.set_values_to, :CNSDTDSC, missing)
                        push_candidate!(dt, c.set_values_to[:EVNTDESC], c.set_values_to[:SRCDOM], c.set_values_to[:SRCVAR], 1, cnsdesc)
                    end
                end
            end
        end

        if nrow(candidates) == 0
            continue  # no TTE record for this subject
        end

        # pick earliest overall (event or censor)
        i = argmin(candidates.ADT)
        chosen = candidates[i, :]
        # (Optional) enforce that ADT must be ≥ STARTDT like {admiral} typical usage.
        # If needed: skip or clamp. Here we keep as-is to mirror user script behavior.

        push!(results, chosen[!, out_cols])
    end

    return results
end

# ---------------------------
# Derive OS
# ---------------------------
adtte_os = derive_param_tte(
    dataset_adsl = adsl,
    start_date = :RANDDT,
    event_conditions = [death_event],
    censor_conditions = [lastalive_censor, rand_censor],
    source_datasets = Dict(:adsl => adsl, :adrs => adrs),
    set_values_to = Dict(:PARAMCD => "OS", :PARAM => "Overall Survival")
)

# ---------------------------
# Derive PFS
# ---------------------------
adtte_pfs = derive_param_tte(
    dataset_adsl = adsl,
    start_date = :RANDDT,
    event_conditions = [pd_event, death_event],
    censor_conditions = [lasta_censor, rand_censor],
    source_datasets = Dict(:adsl => adsl, :adrs => adrs),
    set_values_to = Dict(:PARAMCD => "PFS", :PARAM => "Progression-Free Survival")
)

# Combine OS and PFS like the R pipeline
adtte = vcat(adtte_os, adtte_pfs)

# ---------------------------
# Derive AVAL (duration in days)
# ---------------------------
function derive_vars_duration(df::DataFrame; new_var::Symbol, start_date::Symbol, end_date::Symbol)
    df = copy(df)
    df[!, new_var] = Dates.value.(df[!, end_date] .- df[!, start_date])
    return df
end

adtte_aval = derive_vars_duration(adtte; new_var=:AVAL, start_date=:STARTDT, end_date=:ADT)

# ---------------------------
# Derive ASEQ (within STUDYID, USUBJID ordered by PARAMCD)
# ---------------------------
function derive_var_obs_number(df::DataFrame; by_vars::Vector{Symbol}, order::Vector{Symbol})
    df = sort(df, vcat(by_vars, order))
    df[!, :ASEQ] = 0
    for g in groupby(df, by_vars)
        g[!, :ASEQ] = 1:nrow(g)
    end
    return df
end

adtte_aseq = derive_var_obs_number(adtte_aval; by_vars=[:STUDYID, :USUBJID], order=[:PARAMCD])

# ---------------------------
# Merge ADSL variables
# ---------------------------
adtte_adsl = leftjoin(adtte_aseq, adsl, on=[:STUDYID, :USUBJID])

# ---------------------------
# (Stub) Metadata & export
# ---------------------------
adtte_final = adtte_adsl
# Optional XPT export (uncomment if you have XPTFiles.jl)
# using XPTFiles
# XPTFiles.write_xpt("adtte.xpt", adtte_final)
