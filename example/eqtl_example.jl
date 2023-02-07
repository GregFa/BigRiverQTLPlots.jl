
using RecipesBase, Plots, Plots.PlotMeasures, ColorSchemes
include(joinpath(bulklmmdir, "..", "plot_utils", "visuals_utils.jl"));

multiple_results_allTraits = scan_lite_multivar(pheno_processed, geno_processed, kinship, Threads.nthreads());
gmap_file = joinpath(bulklmmdir,"..","data/bxdData/gmap.csv");
gInfo = CSV.read(gmap_file, DataFrame);
phenocovar_file = joinpath(bulklmmdir,"..","data/bxdData/phenocovar.csv");
pInfo = CSV.read(phenocovar_file, DataFrame);

plot_eQTL(multiple_results_allTraits, pheno, gInfo, pInfo; thr = 5.0)