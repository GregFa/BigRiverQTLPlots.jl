#=
List of the utils functions
- get_eQTL_accMb
    Returns pseudo "acculumated" loci for plotting, 
    according to the phenotype and genotype annotations.

- eQTLplot
    Generates a scatter plot for eQTL analysis.
        
=#


"""
get_eQTL_accMb(mlodmax::Matrix{Float64}, dfpInfo::DataFrame, dfgInfo::DataFrame; 
               chrColname::String = "Chr", mbColname::String = "Mb", threshold::Float64 = 5.0)=> 


Returns pseudo "acculumated" loci for plotting, according to the phenotype and genotype annotations.

## Arguments
- `multiLODs` is the matrix containing the LOD's values.
- `dfpInfo` is a dataframe containing the phenotype information such as probeset, chromosomes names and Mb distance
- `dfgInfo` is a dataframe containing the genotype information such as locus, cM distance, chromosomes names and Mb distance  
- `chrColname` column name containing the chromosomes information, default is `"Chr"`, in the dataframes. 
- `mbColname` column name containing the Mb distance information, default is `"Mb"`, in the dataframes.
- `threshold` is the LOD threshold value, default is `5.0``.

## Output
- `pheno_gmap_lod.acc_geno_mb` vector contains pseudo "accumulated" loci (mB) based on genotype.
- `pheno_gmap_lod.acc_phenocovar_mb` vector contains pseudo "accumulated" loci (mB) based on phenotype.
- `pheno_gmap_lod.maxlod` vector contains the maximum LOD scores above the threshold value `threshold`.  
- `steps` contains the accumulated version of all maximimum loci above the threshold value `threshold`.
- `vChrNames` contains the chromosome names.

"""
function get_eQTL_accMb(multiLODs::Matrix{Float64}, dfpInfo::DataFrame, dfgInfo::DataFrame; 
                        chrColname::String = "Chr", mbColname::String = "Mb", threshold::Float64 = 5.0)
    
    # find the maximum LOD among each columns (gInfo) for each rows (pInfo)
    maxLODs_allTraits = mapslices(x -> findmax(x), multiLODs; dims = 1);
    maxLODs_allTraits = reduce(vcat, vec(maxLODs_allTraits));
    
    # mlodmax is the matrix containing the maximum value of LOD score of each 
    # phenotype and its corresponding index
    mlodmax = Array{Float64, 2}(undef, size(multiLODs, 2), 2);

    for i in 1:size(multiLODs, 2)
        mlodmax[i, 1] = maxLODs_allTraits[i][2];
        mlodmax[i, 2] = maxLODs_allTraits[i][1];
    end
    
    # match chromosomes in pheno dataframe according to chromosomes list in geno dataframe
    dfpInfo_filtered = copy(dfpInfo) #match_chrs_pheno_to_geno(dfpInfo, dfgInfo)

    # prepare filtered pheno dataframe to compute accumulated mb distance for plotting
    rename!(dfpInfo_filtered, Dict(Symbol(mbColname) => :phenocovar_mb));
    rename!(dfpInfo_filtered, Dict(Symbol(chrColname) => :phenocovar_chr)); 
    dfpInfo_filtered[:, :acc_phenocovar_mb] .= dfpInfo_filtered.phenocovar_mb;

    # get a copy of genotype info dataframe
    gmap = copy(dfgInfo) 
    gmap.acc_geno_mb .= gmap[!, mbColname]
    rename!(gmap, Dict(Symbol(mbColname) => :geno_mb));
    rename!(gmap, Dict(Symbol(chrColname) => :geno_chr)); 
    
    # get unique chromosomes names from genotype info dataframe
    vChrNames = unique(gmap.geno_chr);
    vChrNames = sortnatural(String.(vChrNames))

    # initiate steps vector
    steps = Array{Float64}(undef, length(vChrNames))

    # if more than one chromosome
    if length(vChrNames) > 1 
        # compute accumulated mb distance
        for i in 1:length(vChrNames)-1 
            
            # get temp matrix based on a chromosome name
            phenotemp = filter(:phenocovar_chr => x-> x == vChrNames[i], dfpInfo_filtered)
            genotemp = filter(:geno_chr => x-> x == vChrNames[i], gmap)

            # get maximum distance for this chromosome
            steps[i] = max(
                        maximum(genotemp.acc_geno_mb),
                        maximum(phenotemp.acc_phenocovar_mb) 
                    )

            # calculate the accumulated distance for the next chromosome
            nextchr_phenotemp = view(dfpInfo_filtered, dfpInfo_filtered.phenocovar_chr .== vChrNames[i+1], :)
            nextchr_genotemp = view(gmap, gmap.geno_chr .== vChrNames[i+1], :) 

            nextchr_phenotemp.acc_phenocovar_mb .= nextchr_phenotemp.phenocovar_mb .+ steps[i]
            nextchr_genotemp.acc_geno_mb .= nextchr_genotemp.geno_mb .+ steps[i]
        end
    end

    steps[end] = max(
                  maximum(dfpInfo_filtered.acc_phenocovar_mb), 
                  maximum(gmap.acc_geno_mb)
                )

    # get the corresponding index of the maximum values of LOD score of each phenotype
    idxlodmax = trunc.(Int, mlodmax[:,1])
    # concatenate results and gmap info
    pheno_gmap_lod = hcat(dfpInfo_filtered, gmap[idxlodmax,:], DataFrame(idx = idxlodmax, maxlod = mlodmax[:,2]));

    # filter results according to the chromosome names in geno file
    pheno_gmap_lod = filter(:phenocovar_chr => in(vChrNames), pheno_gmap_lod);

    # filter according to LOD threshold
    pheno_gmap_lod = filter(row -> row.maxlod > threshold, pheno_gmap_lod);

    return  pheno_gmap_lod.acc_geno_mb, pheno_gmap_lod.acc_phenocovar_mb,  pheno_gmap_lod.maxlod, steps, vChrNames
end


"""
ploteQTL(multiLODs::Array{Float64, 2}, dfpInfo::DataFrame, dfgInfo::DataFrame;
         threshold::Float64 = 5.0, kwargs...)

Generates a scatter plot for eQTL analysis.

## Arguments
- `multiLODs` is the matrix containing the LOD's values.
- `dfpInfo` is a dataframe containing the phenotype information such as probeset, chromosomes names and Mb distance.
- `dfgInfo` is a dataframe containing the genotype information such as locus, cM distance, chromosomes names and Mb distance. 
- `threshold` is the LOD threshold value, default is `5.0``.

"""
function ploteQTL(multiLODs::Array{Float64, 2}, dfpInfo::DataFrame, dfgInfo::DataFrame;
                  threshold::Float64 = 5.0, kwargs...)

    # get coordinates ready for plotting
    x, y, z, mysteps, mychr = get_eQTL_accMb(
                        multiLODs, 
                        dfpInfo,
                        dfgInfo;
                        threshold = threshold,
                        kwargs...
                    )

    eqtlplot(x, y, z, mysteps, mychr, kwargs...)

end
