# BigRiverQTLPlots.jl

`BigRiverQTLPlots.jl` is a versatile plotting package built in the Julia programming language. The package consists of specific plotting recipes, designed to streamline data visualization and enhance the process of statistical analysis in genetic studies. It's particularly suited for QTL (Quantitative Trait Loci) and eQTL (expression Quantitative Trait Loci) analyses, providing key features for clear and intuitive data representation.

## Features
Here are the main functions provided by BigRiverQTLPlots.jl:

- `plot_QTL()`: This function generates plots for LOD (logarithm of odds) scores with respect to marker positions. It's useful in viewing one genome scan result, or even multiple genome scan results, on a single plot. This helps to provide a broad overview of the QTL landscape.

- `plot_eQTL()`: This function is specifically designed to plot eQTL analysis results, assisting in visualization of genetic associations with gene expression levels.

## Installation
To install `BigRiverQTLPlots.jl`, you can use Julia's package manager. Here is the command:

```julia
using Pkg
Pkg.add("BigRiverQTLPlots")

```

## Usage
After installing `BigRiverQTLPlots.jl`, you can include it in your Julia script using the following command:

```julia
using BigRiverQTLPlots
```

From there, you can start using `plot_QTL`, `plot_eQTL`, and `confplot` to plot your data. For example:

```julia
# Assuming `single_results_perms` are restulting lod scores
# gInfo contains genotype information  
plot_QTL(single_results_perms, gInfo)


# Assuming `lod_scores` are in multipletraits_results
# pInfo contains phenotype information
# gInfo contains genotype information  
# thresh is your LOD threshold value
plot_eQTL(multipletraits_results, pInfo, gInfo; threshold = 5.0)
plot_QTL(single_results_perms, gInfo)


# Assuming `confidence_data` is your data
confplot(confidence_data)
```

## Contribution
Contributions to BigRiverQTLPlots.jl are welcome and appreciated. If you'd like to contribute, please fork the repository and make changes as you'd like. If you have any questions or issues, feel free to open an issue on the repository.

## License
`BigRiverQTLPlots.jl` is licensed under the MIT license. For more information, please refer to the LICENSE file in the repository.

## Support
If you have any problems or questions using `BigRiverQTLPlots.jl`, please open an issue on the GitHub repository. We'll be happy to help!