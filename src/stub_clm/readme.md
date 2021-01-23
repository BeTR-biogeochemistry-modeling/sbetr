# stub_clm

Data structure for coupling betr to elm/clm
|Source File | Description|
|------|------|
|CNCarbonFluxType.F90          |  Column carbon flux data type |
|PfCarbonFluxType.F90    | Patch carbon flux data type|
|CNCarbonStateType.F90        |Column carbon flux data type|
|CNDecompCascadeConType.F90     |Parameter structure defines basic parameters for decomposition|
|CNNitrogenFluxType.F90       |Column nitrogen flux data type|
|PfNitrogenFluxType.F90         |Patch nitrogen flux data type|
|CNNitrogenStateType.F90        |Column nitrogen state data type|
|PfNitrogenStateType.F90        |Patch nitrogen state data type|
|CNSharedParamsMod.F90          |Shared data strcuture for CN and cent-like bgc|
|CNStateType.F90               |Column data type for CN bgc|
|CanopyStateType.F90         |Canopy related variables|
|ChemStateType.F90            |Soil chemical state variables|
|ColumnType.F90        |Variables define the column type|
|EcophysConType.F90             |Ecosystem level parameters|
|FuncPedotransferMod.F90        |Different pedotransfer functions|
|GridcellType.F90               |Variables define the grid type|
|LandunitType.F90         |Variables define the land unit type|
|PatchType.F90              |Variables define the patch type|
|lnd2atmType.F90            |Variables for land-atmosphere exchange|
|PhosphorusFluxType.F90       |Column phosphorus flux data type|
|PfPhosphorusFluxType.F90     |Patch phosphorus flux data type|                  
|PhosphorusStateType.F90        |Column phosphorus state data type|
|QSatMod.F90                   |Subroutines for vapor pressure calculation|
|SoilHydrologyType.F90          |Column level hydroloy related variables|
|SoilStateType.F90           |Soil properties|
|SoilWaterRetentionCurveClappHornberg1978Mod.F90 |Soil water retention curve based on Clapp Hornberg 1978 paper|
|SoilWaterRetentionCurveFactoryMod.F90    |Module interface for flexible soil water retention curve parameterization|
|SoilWaterRetentionCurveMod.F90   |Base type for soil water retention curve parameterization|
|SurfaceResistanceMod.F90       |Soil resistance for air-soil gas exchange|
|TemperatureType.F90            |Column soil temperature data type|
|PfTemperatureType.F90          |Patch soil temperature data type|
|TridiagonalMod.F90            |Tridiagnoal solver|
|WaterFluxType.F90           |Column water flux data type|
|PfWaterFluxType.F90         |Patch water flux data type|
|WaterStateType.F90        |Water state variables|
|accumulMod.F90                |Code to support temporal average of user defined variables|
|atm2lndType.F90          |Atmospheric forcing type|
|clm_initializeMod.F90          |Module to initialize clm|
|clm_instMod.F90            |Module to initialize clm data types|
|clm_time_manager.F90       |Module to time stamp management|
|clm_varcon.F90            |Define mathematical constants|
|clm_varctl.F90          |Define logical switches that control the model|
|clm_varpar.F90             |Define parameters related to model data structure|
|clmgridMod.F90               |Define and initialize clm grid|
|column_varcon.F90           |Column constant variables|
|debuginfoMod.F90            |Subroutines to assist model debug|
|decompMod.F90                |Data structure to support parallelization|
|PlantMicKineticsMod.F90         |Kinetic data type to link the big leaf module with soil bgc|
|fileutils.F90                 |Subroutines for file management|
|getdatetime.F90             |Date time subroutine|
|histFileMod.F90                |History file data structure|
|landunit_varcon.F90            |Land unit constants|
|pftvarcon.F90                 |Patch constants|
|restUtilMod.F90               |Subroutines to support restart file I/O|
|soilorder_varcon.F90          |Soil order parameters|
|spmdMod.F90                 |Mpi related variables|
|subgridAveMod.F90   |Subroutines to do subgrid average                 |
