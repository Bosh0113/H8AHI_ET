import numpy as np
from numpy.ma import exp, log

# Priestley-Taylor coefficient alpha
PRIESTLEY_TAYLOR_ALPHA = 1.26
BETA = 1.0
PSYCHROMETRIC_GAMMA = 0.0662 # Pa/K # http://www.fao.org/docrep/x0490e/x0490e07.htm
KRN = 0.6
KPAR = 0.5
DEFAULT_FLOOR_SATURATION_VAPOR_PRESSURE=1.


def filter_bad_values(matrix, lower_bound, upper_bound):
    matrix[matrix < lower_bound] = np.nan
    matrix[matrix > upper_bound] = np.nan
    
    return matrix


# saturation vapor pressure in kPa from air temperature in celsius
def saturation_vapor_pressure_from_air_temperature(air_temperature):
    # SVP_BASE = 0.611
    SVP_BASE = 0.61121
    # SVP_MULT = 17.27
    SVP_MULT = 17.502
    # SVP_ADD = 237.7
    SVP_ADD = 240.97
    svp = SVP_BASE * exp((air_temperature * SVP_MULT) / (air_temperature + SVP_ADD))
    
    return svp


# calculate Soil-Adjusted Vegetation Index from Normalized Difference Vegetation Index
# using linear relationship
def savi_from_ndvi(ndvi):
    SAVI_MULT = 0.45
    SAVI_ADD = 0.132
    savi = ndvi * SAVI_MULT + SAVI_ADD
    
    return savi


# calculate fraction of absorbed photosynthetically active radiation
# from Soil-Adjusted Vegetation Index using linear relationship
def fAPAR_from_savi(savi):
    FAPAR_MULT = 1.3632
    FAPAR_ADD = -0.048
    
    fAPAR = savi * FAPAR_MULT + FAPAR_ADD
    
    return fAPAR


# calculate fraction of absorbed photosynthetically active radiation
# from Normalized Difference Vegetation Index
def fAPAR_from_ndvi(ndvi):
    savi = savi_from_ndvi(ndvi)
    fAPAR = fAPAR_from_savi(savi)
    
    return fAPAR


# calculate slope of saturation vapor pressure to air temperature
# in pascals over kelvin
def delta_from_air_temperature(air_temperature):
    # return 4098 * (0.6108 * exp(17.27 * air_temperature / (237.7 + air_temperature))) / (air_temperature + 237.3) ** 2
    return 240.97*17.502*(0.61121 * exp((air_temperature * 17.502) / (air_temperature + 240.97))) / (air_temperature + 240.97) ** 2


# cut max and min values to defined values
def enforce_boundaries(matrix, lower_bound, upper_bound):
    matrix[matrix < lower_bound] = np.nan
    matrix[matrix > upper_bound] = np.nan
    
    return matrix


# calculate fraction of intercepted photosynthetically active radiation
# from Normalized Difference Vegetation Index
def fIPAR_from_ndvi(ndvi):
    FIPAR_ADD = -0.05
    fIPAR = ndvi + FIPAR_ADD
    
    return fIPAR


# calculate the relative stress from Temperature
# this is for high temperatures
def fT_fun(Tmax,Topt):
    '''plant temperature constraint '''
    fT = np.exp(-((Tmax-Topt)/Topt)**2)
    fT[Tmax<-5]=0.05; # <---- add a cold temperature constraint
    return fT


# theory trying to capture temp at peak photosynthesis (wet, green, high vpd, high radiation)
def Topt_fun(RN,TMAX,SAVI,VPD):
    ''' returns the optimal temperature
        inputs 14 day averages of: RN, TACTUAL, SAVI, VPD'''
    topt_mask = ~np.isnan(np.array(RN)) & ~np.isnan(TMAX) & ~np.isnan(SAVI) & ~np.isnan(VPD)
    
    PHEN_RAW = np.ones(np.shape(RN)); PHEN_RAW[:]=np.nan

    VPD_g1 = VPD
    VPD_g1[VPD<0.5]=0.5
    PHEN_RAW[topt_mask] = RN[topt_mask]*TMAX[topt_mask]*SAVI[topt_mask]/VPD_g1[topt_mask]

    max_idx = PHEN_RAW.argmax()
    T_opt = TMAX[max_idx]
    
    if np.isnan(T_opt):
        T_opt = np.nanmax(np.nanmax(TMAX))-5

    return T_opt


def ptjpl_area(r_net_input_array, 
               rh_input_array, 
               ta_input_array, 
               ndvi_input_array, verbose=True, floor_saturation_vapor_pressure=DEFAULT_FLOOR_SATURATION_VAPOR_PRESSURE):
    """
        :AA: is a dataframe from where each variable listed below is extracted
        I have attached a csv file containing what each variable name is provided.  

        :param air_temperature_K:
            Numpy matrix of air temperature near the surface in kelvin
        :param air_temperature_mean_K:
            Numpy matrix of average air temperature near the surface in kelvin
        :param ndvi_mean:
            Numpy matrix of average Normalized Difference Vegetation Index
        :param net_radiation:
            Numpy matrix of instantaneous net radiation in watts per square meter
        :param daily_radiation:
            Numpy matrix of daily net radiation in watts per square meter

        :param water_vapor_pressure_mean_Pa:
            Numpy matrix of average vapor pressure in pascals
        :param optimum_temperature:
            Numpy matrix of phenologically optimum temperature
        :param fAPARmax:
            Numpy matrix of maximum fraction of photosynthetically active radiation

        :param verbose:
            Flag to output activity to console
        :param floor_saturation_vapor_pressure:
            Option to floor calculation of saturation vapor pressure at 1 to avoid anomalous output

        :return:
            Dataframe with:
                evapotranspiration, 
                PTJPL original model potential_evapotranspiration, 
                T and S and I daily_evapotranspiration, 
                PTJPL original scaled with EF and adiation
        """
    # air_temperature_K =       np.array(AA.TA)          # (K)
    air_temperature_K =       ta_input_array        # (K)
    # air_temperature_mean_K =  np.array(AA.TA_day_mean) # (K)
    # ndvi_mean=                np.array(AA.NDVI_day_mean)               # (0-1.0)
    ndvi_mean=                ndvi_input_array              # (0-1.0)
    # net_radiation=            np.array(AA.NETRAD);            # (W/m2)
    net_radiation=            r_net_input_array           # (W/m2)
    # daily_radiation=          np.array(AA.NETRAD_day)         # (W/m2)
    
    # RH =                      np.array(AA.RH/100.)         # (%), Can run with vapor pressure see comments out below
    RH =                      rh_input_array/100.         # (%), Can run with vapor pressure see comments out below

    #     water_vapor_pressure_mean_Pa =np.array(AA.VPD_day_max)*100# set VPD of input dataset to Pa
    #     optimum_temperature=18; # use the function Topt_fun to return this value [ constrain to +5 o C, if below a certain value, set fT to 1.0]
    #     fAPARmax=0.65 # grab this value from MODIS timeseries dataset

    # convert temperatures from kelvin to celsius
    air_temperature = air_temperature_K - 273.15
    # AA['TA_C']= air_temperature
    # air_temperature_mean = air_temperature_mean_K - 273.15
    
    #     scale water vapor pressure from Pa to kPa
    #     water_vapor_pressure_mean = water_vapor_pressure_mean_Pa * 0.001
    
    # calculate surface wetness values
    if verbose:
        print('calculating surface wetness values [%]')
        print('calculating vapor pressure deficit [kPa]')
    
    # calculate saturation vapor pressure in kPa from air temperature in celcius
    # saturation_vapor_pressure = saturation_vapor_pressure_from_air_temperature(air_temperature_mean)#[kPA]
    saturation_vapor_pressure = saturation_vapor_pressure_from_air_temperature(air_temperature)#[kPA]
    
    # ***REPLACED THIS WITH ACTUAL MEASTUREMENTS FROM TOWERS:
    # calculate relative humidity from water vapor pressure and saturation vapor pressure
    #     relative_humidity = water_vapor_pressure_mean / saturation_vapor_pressure
    # upper bound of relative humidity is one, results higher than one are capped at one
    # WE INSTEAD HAVE RH SO WE CALCULATE VPD FROM THIS:
    relative_humidity = filter_bad_values(RH,0.,1.)
    
    # floor saturation vapor pressure at 1
    if floor_saturation_vapor_pressure:
        saturation_vapor_pressure[saturation_vapor_pressure < 1] = 1
    water_vapor_pressure_mean = RH*saturation_vapor_pressure
    
    # calculate vapor pressure deficit from water vapor pressure
    vapor_pressure_deficit = saturation_vapor_pressure - water_vapor_pressure_mean # [kPa]
    
    # lower bound of vapor pressure deficit is zero, negative values replaced with nodata
    vapor_pressure_deficit[vapor_pressure_deficit < 0] = np.nan
    # AA['VPD_roll']=vapor_pressure_deficit
    # AA['VPD']= vapor_pressure_deficit
    
    # calculate relative surface wetness from relative humidity
    # AA['RH_roll']=relative_humidity
    relative_surface_wetness = relative_humidity ** 4
    relative_surface_wetness[air_temperature<=0]=0.; # FROZEN WATER
    # calculate slope of saturation to vapor pressure curve Pa/K
    delta = delta_from_air_temperature(air_temperature)
    # calculate delta / (delta + gamma)
    epsilon = delta / (delta + PSYCHROMETRIC_GAMMA)
    
    # calculate vegetation values
    if verbose:
        print('calculating vegetation values')
    
    # calculate fAPAR & fAPARmax from NDVI mean
    fAPAR = fAPAR_from_ndvi(ndvi_mean)
    fAPAR = enforce_boundaries(fAPAR,0.,1.)
    fAPARmax = np.nanmax(fAPAR)
    # AA['fAPAR']=fAPAR

    # calculate fIPAR from NDVI mean
    fIPAR = fIPAR_from_ndvi(ndvi_mean)
    # AA['fIPAR']=fIPAR

    # calculate green canopy fraction (fg) from fAPAR and fIPAR, constrained between zero and one
    green_canopy_fraction = enforce_boundaries(fAPAR / fIPAR, 0, 1)
    
    # calculate plant moisture constraint (fM) from fraction of photosynthetically active radiation,
    # constrained between zero and one
    plant_moisture_constraint = enforce_boundaries(fAPAR / fAPARmax, 0, 1)

    # calculate soil moisture constraint from relative humidity and vapor pressure deficit,
    # constrained between zero and one
    # soil_moisture_constraint = enforce_boundaries(AA.RH_roll.rolling(14,1).mean() ** (AA.VPD_roll.rolling(14,1).mean() / BETA), 0, 1) # github
    soil_moisture_constraint = enforce_boundaries(relative_humidity ** vapor_pressure_deficit, 0, 1)
    soil_moisture_constraint[air_temperature<0]=0.00
    # AA['soil_moisture_constraint']=soil_moisture_constraint
    
    # calculate SAVI from NDVI
    savi_mean = savi_from_ndvi(ndvi_mean)
    # AA['savi']=savi_mean

    # calculate optimum_temperature from: R𝑛 , Tair, SAVI, and VPD used to calculate Topt are all 2-week forward averages. (14*24*6 data)
    # optimum_temperature = Topt_fun(np.array(AA.NETRAD_day.rolling(30,1).mean()),np.array(air_temperature_mean.rolling(30,1).mean()),np.array(AA['savi'].rolling(30,1).mean()),np.array(AA['VPD'].rolling(30,1).mean()))
    # optimum_temperature = Topt_fun(np.array(AA.NETRAD.rolling(14*24*6,24*6).mean()),
    #                                np.array(AA.TA_C.rolling(14*24*6,24*6).mean()),
    #                                np.array(AA['savi'].rolling(14*24*6,24*6).mean()),
    #                                np.array(AA['VPD'].rolling(14*24*6,24*6).mean()))
    optimum_temperature = 23.5 #!!!!!

    if verbose:
        print('calculating plant optimum temperature')
        print(str(optimum_temperature)+ ' C')

    # calculate plant temperature constraint (fT) from optimal phenology
    # plant_temperature_constraint = fT_fun(air_temperature_mean, optimum_temperature)
    plant_temperature_constraint = fT_fun(air_temperature, optimum_temperature)
    # exp(-(((air_temperature_mean - optimum_temperature) / optimum_temperature) ** 2))

    # calculate leaf area index : now extract from MODIS dataset
    leaf_area_index = -log(1 - fIPAR) * (1 / KPAR)
    # leaf_area_index = np.array(AA.LAI);
    
    # soil evaporation
    if verbose:
        print('calculating soil evaporation')

    # caluclate net radiation of the soil from leaf area index
    soil_net_radiation = net_radiation * exp(-KRN * leaf_area_index)

    # take fractional vegetation cover from fraction of photosynthetically active radiation
    fractional_vegetation_cover = enforce_boundaries(fIPAR, 0, 1)

    # calculate instantaneous soil heat flux from net radiation and fractional vegetation cover
    soil_heat_flux = net_radiation * (0.05 + (1 - fractional_vegetation_cover) * 0.265)
    # change the above to METRIC G
    soil_heat_flux[soil_heat_flux < 0] = 0
    soil_heat_flux[soil_heat_flux > 0.35 * soil_net_radiation] = 0.35 * soil_net_radiation[soil_heat_flux > 0.35 * soil_net_radiation]

    # calculate soil evaporation (LEs) from relative surface wetness, soil moisture constraint,
    # priestley taylor coefficient, epsilon = delta / (delta + gamma), net radiation of the soil,
    # and soil heat flux
    soil_evaporation = (relative_surface_wetness + soil_moisture_constraint * (1 - relative_surface_wetness)) * PRIESTLEY_TAYLOR_ALPHA * epsilon * (soil_net_radiation - soil_heat_flux)
    #     soil_evaporation[numpy.isnan(soil_evaporation)] = 0
    soil_evaporation[soil_evaporation < 0] = np.nan

    # canopy transpiration
    if verbose:
        print('calculating canopy transpiration')

    # calculate net radiation of the canopy from net radiation of the soil
    canopy_net_radiation = net_radiation - soil_net_radiation

    # calculate canopy transpiration (LEc) from priestley taylor, relative surface wetness,
    # green canopy fraction, plant temperature constraint, plant moisture constraint,
    # epsilon = delta / (delta + gamma), and net radiation of the canopy
    canopy_transpiration = PRIESTLEY_TAYLOR_ALPHA * (1 - relative_surface_wetness) * green_canopy_fraction * plant_temperature_constraint * plant_moisture_constraint * epsilon * canopy_net_radiation
    canopy_transpiration[np.isnan(canopy_transpiration)] = 0
    canopy_transpiration[canopy_transpiration < 0] = 0

    # interception evaporation
    if verbose:
        print('calculating interception evaporation')
    # calculate interception evaporation (LEi) from relative surface wetness, priestley taylor coefficient,
    # epsilon = delta / (delta + gamma), and net radiation of the canopy
    interception_evaporation = relative_surface_wetness * PRIESTLEY_TAYLOR_ALPHA * epsilon * canopy_net_radiation
    interception_evaporation[interception_evaporation < 0] = 0
    
    # combined evapotranspiration
    if verbose:
        print('combining evapotranspiration')
    # combine soil evaporation (LEs), canopy transpiration (LEc), and interception evaporation (LEi)
    # into instantaneous evapotranspiration (LE)
    evapotranspiration = soil_evaporation + canopy_transpiration + interception_evaporation
    evapotranspiration[evapotranspiration > net_radiation] = net_radiation[evapotranspiration > net_radiation]
    evapotranspiration[np.isinf(evapotranspiration)] = np.nan
    evapotranspiration[evapotranspiration < 0] = np.nan

    # # daily evapotranspiration
    # if verbose:
    #     print('calculating daily evapotranspiration')
    # # calculate evaporative fraction (EF) from evapotranspiration, net radiation, and soil heat flux
    # evaporative_fraction = evapotranspiration / (net_radiation - soil_heat_flux)

    # # calculate daily evapotranspiration from daily net radiation and evaporative fraction
    # daily_evapotranspiration = daily_radiation * evaporative_fraction
    
    # potential evapotranspiration
    if verbose:
        print('calculating potential evapotranspiration')
    # calculate potential evapotranspiration (pET) from priestley taylor coefficient,
    # epsilon = delta / (delta + gamma), net radiation, and soil heat flux
    potential_evapotranspiration = PRIESTLEY_TAYLOR_ALPHA * epsilon * (net_radiation - soil_heat_flux)
    # potential_transpiration      = PRIESTLEY_TAYLOR_ALPHA * epsilon * (canopy_net_radiation)
    # potential_evaporation        = PRIESTLEY_TAYLOR_ALPHA * epsilon * (soil_net_radiation - soil_heat_flux)
    
    # # append new data to array
    # results = AA
    # results['evapotranspiration'] =           evapotranspiration    
    # # potential et added to results
    # results['potential_evapotranspiration']=  potential_evapotranspiration
    # # components evapotranspriation added for partitioning study
    # results['canopy_transpiration'] =         canopy_transpiration
    # results['interception_evaporation'] =     interception_evaporation
    # results['soil_evaporation'] =             soil_evaporation

    results = [evapotranspiration, canopy_transpiration, interception_evaporation, soil_evaporation, potential_evapotranspiration]
    return np.array(results)