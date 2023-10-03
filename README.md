# SSIEEG simulations
Run Simbio simulations to test leadfield differences on a skull that has burrholes.

### Requirements
- [fieldtrip](https://github.com/fieldtrip/fieldtrip)
- [RainCloudPlots](https://github.com/RainCloudPlots/RainCloudPlots)
- [Robust Statistical Toolbox](https://github.com/CPernet/Robust_Statistical_Toolbox)

### Analyses
1. [fem_analysis_mni.mlx](/fem_analysis_mni.mlx)
   - Build head models
   - Compute leadfield differences at standard burrhole sizes

2. [fem_analysis_runmany.m](/fem_analysis_runmany.m)
   - Run leadfield computations at varying burrhole sizes

3. [fem_analysis_plottingmany.m](/fem_analysis_plottingmany.m)
   - Plot results from `fem_analysis_runmany.m`
   
4. [fem_distancetohole.mlx](/fem_distancetohole.mlx)
   - Simulate leadfield differences as a function of distance to burrhole
   
