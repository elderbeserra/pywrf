import PyNGL_numpy.Ngl as ngl

def res_init(lat, lon, 
  lat_min=None,
  lat_max=None,
  lon_min=None,
  lon_max=None):
    """
    Assumes lat,lon are numpy arrays
    """
    if not lat_min:
        lat_min = lat.min()
    if not lat_max:
        lat_max = lat.max()
    if not lon_min:
        lon_min = lon.min()
    if not lon_max:
        lon_max = lon.max()
    res = ngl.Resources()
    res.mpProjection = 'Orthographic'
    res.mpCenterLonF = lon_min + (lon_max-lon_min)/2.
    res.mpCenterLatF = lat_min + (lat_max-lat_min)/2.
    res.mpLimitMode = 'LatLon'
    res.mpMinLatF = lat_min - 0.2
    res.mpMaxLatF = lat_max + 0.2
    res.mpMinLonF = lon_min - 0.2
    res.mpMaxLonF = lon_max + 0.2
    res.mpFillOn = True
    res.mpGridSpacingF = 2.
    res.vpXF = 0.1
    res.vpYF = 0.9
    res.vpWidthF = 0.7
    res.vpHeightF = 0.7
    return res

def make_land_gray(wks,res):
    ic = ngl.new_color(wks, 0.75, 0.75,0.75)
    res.mpFillColors = [0,-1,ic,-1]

def mslp_contour_levels(res):
    mnlvl = 980
    mxlvl = 1024
    spcng = 2
    ncn = (mxlvl - mnlvl) / spcng + 1
    res.cnLevelSelectionMode = 'ManualLevels'
    res.cnMinLevelValF = mnlvl
    res.cnMaxLevelValF = mxlvl
    res.cnLevelSpacingF = spcng

def sfc_pres_contour_levels(res):
    mnlvl = 840
    mxlvl = 1024
    spcng = 2
    ncn = (mxlvl - mnlvl) / spcng + 1
    res.cnLevelSelectionMode = 'ManualLevels'
    res.cnMinLevelValF = mnlvl
    res.cnMaxLevelValF = mxlvl
    res.cnLevelSpacingF = spcng  


