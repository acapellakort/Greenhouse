# =============================================================================
# train_sac_multiseed.jl -- repeat the warm-started SAC over several seeds and
#   report SAC mean +/- std vs the (seed-independent) CEM constant baseline.
# =============================================================================
# Long run: each seed is a full 24k-step training + 31-year eval (~12-15 min).
#   julia -t auto --project=. test/train_sac_multiseed.jl
include(joinpath(@__DIR__, "train_sac.jl"))     # defines run_seed(...) (guarded: no auto-run)
using Printf, Statistics

const SEEDS = [11, 22, 33]                       # 3 seeds (~45 min); add more for tighter bars

@printf("multi-seed SAC: %d seeds x %d steps\n", length(SEEDS), SAC().total_steps); flush(stdout)
res = NamedTuple[]
for (i, s) in enumerate(SEEDS)
    @printf("--- seed %d  (%d/%d) ---\n", s, i, length(SEEDS)); flush(stdout)
    r = run_seed(s; verbose=false, save=(i==1))
    push!(res, r)
    @printf("  held-out %.2f | all %.2f | tuned %.2f | worst %.2f\n",
            r.sac_held, r.sac_all, r.sac_tuned, r.sac_worst); flush(stdout)
end

ms(f) = (v = Float64[f(r) for r in res]; (mean(v), length(v) > 1 ? std(v) : 0.0))
println("\n==== SAC across $(length(SEEDS)) seeds (mean +/- std) vs CEM constant ====")
@printf("%-16s %18s %10s\n", "split", "SAC mean+/-std", "CEM")
for (name, fs, cemv) in (("tuned 4yr",     r->r.sac_tuned, res[1].cem_tuned),
                         ("held-out 27yr", r->r.sac_held,  res[1].cem_held),
                         ("all 31yr",      r->r.sac_all,   res[1].cem_all),
                         ("worst yr",      r->r.sac_worst, res[1].cem_worst))
    m, sd = ms(fs)
    @printf("%-16s   %7.2f +/- %4.2f     %7.2f\n", name, m, sd, cemv)
end
@printf("\nSAC held-out beats CEM in %d / %d seeds\n",
        count(r->r.sac_held > res[1].cem_held, res), length(res))
