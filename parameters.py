'''Sets user-defined parameters for calculating monthly mean IMB data, e.g. conductivity scheme, layers over which to calculate conductive flux'''

sfctemp_melt_threshold=-2.
layer_fcondbot_calc = [.5,.9]
layer_shu_calc =[0.,.5]
thin_snow_ice_threshold = 1.
thin_ice_threshold_for_fcondbot = layer_fcondbot_calc[1] * 1.
thin_ice_threshold_for_shu = layer_shu_calc[1] * 1.
layer_shu_label = '{:4.3f}_{:4.3f}'.format(*layer_shu_calc)
layer_fcondbot_label = '{:4.3f}_{:4.3f}'.format(*layer_fcondbot_calc)
snow_threshold_for_topmelt = 0.015
conductivity_formulation = 'maykut_untersteiner'
rhos = 330.
rhoi = 917.
