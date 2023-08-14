
using BigRiverQTLPlots
using BulkLMM
using Random, Statistics
using Plots

##############
# BXD spleen #
##############

########
# Data #
########
bulklmmdir = dirname(pathof(BulkLMM));

gmap_file = joinpath(bulklmmdir, "..", "data", "bxdData", "gmap.csv");
gInfo = BulkLMM.CSV.read(gmap_file, BulkLMM.DataFrames.DataFrame);

phenocovar_file = joinpath(bulklmmdir, "..", "data", "bxdData", "phenocovar.csv");
pInfo = BulkLMM.CSV.read(phenocovar_file, BulkLMM.DataFrames.DataFrame);

pheno_file = joinpath(bulklmmdir, "..", "data", "bxdData", "spleen-pheno-nomissing.csv");
pheno = BulkLMM.DelimitedFiles.readdlm(pheno_file, ',', header = false);
pheno_processed = pheno[2:end, 2:(end-1)] .* 1.0; # exclude the header, the first (transcript ID)and the last columns (sex)

geno_file = joinpath(bulklmmdir, "..", "data", "bxdData", "spleen-bxd-genoprob.csv")
geno = BulkLMM.DelimitedFiles.readdlm(geno_file, ',', header = false);
geno_processed = geno[2:end, 1:2:end] .* 1.0;

#################
# Preprocessing #
#################
traitID = 1112;
pheno_y = reshape(pheno_processed[:, traitID], :, 1);

###########
# Kinship #
###########
kinship = calcKinship(geno_processed);


########
# Scan #
########

@time single_results_perms = scan(
	pheno_y,
	geno_processed,
	kinship;
	permutation_test = true,
	nperms = 1000,
);

@time single_results = scan(
	pheno_y,
	geno_processed,
	kinship,
);

thrs = BigRiverQTLPlots.perms_thresholds(single_results_perms.L_perms, [0.90, 0.95]);
println(thrs)

x, y, vecSteps, v_chr_names = get_plot_QTL_inputs(single_results.lod, gInfo)

########
# Plot #
########

mycolor = ["#1f78b4", "#a6cee3", "#756bb1", "#bcbddc"]

plot_manhattan(single_results, gInfo, 
markersize = 3,
markercolor= manahattancolor(single_results, mycolor[2] , mycolor[1]),
)


plot_QTL(single_results_perms, gInfo)
plot_QTL(single_results_perms, gInfo)
plot_manhattan(single_results_perms, gInfo, significance = [])

scatter(x, y,  
	grid = false,
	legend = false, 
	xticks = (pseudoticks(vecSteps[2:end]), v_chr_names),
	markercolor= ifelse.(is_chr_odd , mycolor[3] , mycolor[4]),
	markersize = 3,
	markerstrokewidth = 0.1,

	)

	v_chr_names


idx_Inf = vcat([1],findall(y.==Inf))
n = length(y)
is_chr_odd = repeat([true], n)

for i in 1:length(idx_Inf)-1
	is_chr_odd[idx_Inf[i]:idx_Inf[i+1]].= isodd(i)
end
# is_chr_odd[idx_Inf[2:end]] .= Inf

scatter(x, y,
markercolor= ifelse.(is_chr_odd , "blue" , "red")
)




length(is_chr_odd)










plot_QTL(single_results.lod, gInfo, seriestype = :scatter,  marker = 0.5)


plot_QTL(single_results_perms, gInfo, 
seriestype = :scatter,  markersize = 0.5,
thresholds = [0.90, 0.95])




plot(single_results.lod, (seriestype = :scatter,  markersize = 0.5))