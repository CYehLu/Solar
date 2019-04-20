import datetime
import numpy as np

# === Abbreviation ===
# GMLS: Geom Mean Long Sum (deg)
# GMAS: Geom Mean Anom Sum (deg)
# EEO: Eccent Earth Orbit
# SEC: Sun Eq of Ctr
# STL: Sun True Long (deg)
# STA: Sun True Anom (deg)
# SRV: Sun Rad Vector (AUs)
# SAL: Sun App Long (deg)
# MOE: Mean Obliq Ecliptic (deg)
# OC: Obliq Corr (deg)
# SRA: Sun Rt Ascen (deg)
# SD: Sun Declin (deg)
# var_y: var y
# ET: Eq of Time (minutes)
# HAS: HA Time (minutes)
# SN: Solar Noon (LST)
# SrT: Sunrise Time (LST)
# SsT: Sunset Time (LST)
# SlD: Sunlight Duration (minutes)
# TST: True Solar Time (min)
# HA: Hour Angle (deg)
# SZA: Solar Zenith Angle (deg)
# SEA: Solar Elevation Angle (deg)
# AAR: Approx Atmospheric Refraction (deg)
# SEC: Solar Elevation corrected for atm refraction (deg)
# SAA: Solar Azimuth Angle (deg cw from N)
# ====================
def solar_angle(year, month, day, hour, minute, tz, lat, lon):
    """
    Calculate solar elevation & azimuth angle.
    Based on NOAA's Astronomical Algorithms, by Jean Meeus.
    Convert from EXCEL version: 
    https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
    
    Parameters:
    ----------
    year, month, day, hour, minute, tz, lat, lon:
        int. 'tz' is the time zone.
        
    Return:
    ------
    (Elevation_Angle, Azimuth_Angle)
    """
    date = datetime.datetime(year, month, day)
    reference = datetime.datetime(1899, 12, 31)
    julian_day = (date - reference).days + 2415018.5 + (hour*60+minute)/(24*60) - tz/24 + 1
    julian_century = (julian_day-2451545) / 36525

    GMLS = np.mod(280.46646 + julian_century*(36000.76983+julian_century*0.0003032), 360)
    GMAS = 357.52911 + julian_century * (35999.05029 - 0.0001537*julian_century)
    EEO = 0.016708634 - julian_century * (0.000042037+0.0000001267*julian_century)
    SEC = np.sin(np.deg2rad(GMAS)) * (1.914602 - julian_century*(0.004817+0.000014*julian_century))  \
          + np.sin(np.deg2rad(2*GMAS)) * (0.019993-0.000101*julian_century)  \
          + np.sin(np.deg2rad(3*GMAS))*0.000289
    STL = GMLS + SEC
    STA = GMAS + SEC
    SRV = 1.000001018 * (1-EEO**2) / (1 + EEO * np.cos(np.deg2rad(STA)))
    SAL = STL - 0.00569 - 0.00478 * np.sin(np.deg2rad(125.04-1934.136*julian_century))
    MOE = 23+(26+((21.448-julian_century*(46.815+julian_century*(0.00059-julian_century*0.001813))))/60)/60
    OC = MOE + 0.00256 * np.cos(np.deg2rad(125.04-1934.136*julian_century))
    SRA = np.rad2deg(
        np.arctan2(np.cos(np.deg2rad(OC))*np.sin(np.deg2rad(SAL)), np.cos(np.deg2rad(SAL)))
    )
    SD = np.rad2deg(np.arcsin(np.sin(np.deg2rad(OC)) * np.sin(np.deg2rad(SAL))))
    var_y = np.tan(np.deg2rad(OC/2)) ** 2
    ET = 4 * np.rad2deg(
        var_y*np.sin(2*np.deg2rad(GMLS)) - 2*EEO*np.sin(np.deg2rad(GMAS))
        + 4*EEO*var_y*np.sin(np.deg2rad(GMAS))*np.cos(2*np.deg2rad(GMLS))
        - 0.5*var_y**2*np.sin(4*np.deg2rad(GMLS))
        - 1.25*EEO**2*np.sin(2*np.deg2rad(GMAS))
    )
    HAS = np.rad2deg(
        np.arccos(np.cos(np.deg2rad(90.833))/(np.cos(np.deg2rad(lat))*np.cos(np.deg2rad(SD)))
                  - np.tan(np.deg2rad(lat))*np.tan(np.deg2rad(SD)))
    )
    SN = (720 - 4*lon - ET + tz*60) / 1440
    SrT = SN - HAS*4 / 1440
    SsT = SN + HAS*4 / 1440
    SlD = 8 * HAS
    TST = np.mod((hour*60 + minute)/(24*60) * 1440 + ET + 4*lon - 60*tz, 1440)

    if TST/4 < 0:
        HA = TST/4 + 180
    else:
        HA = TST/4 - 180

    SZA = np.rad2deg(
        np.arccos(np.sin(np.deg2rad(lat))*np.sin(np.deg2rad(SD))
                  + np.cos(np.deg2rad(lat))*np.cos(np.deg2rad(SD))*np.cos(np.deg2rad(HA)))
    )
    SEA = 90 - SZA

    radSEA = np.deg2rad(SEA)
    if SEA > 85:
        AAR = 0
    elif SEA > 5:
        AAR = (1 / np.tan(radSEA)) - (0.07 / np.tan(radSEA)**3) + (0.000086 / np.tan(radSEA)**5)
    elif SEA > -0.575:
        AAR = 1735 + SEA * (-518.2 + SEA*(103.4+SEA*(-12.79+SEA*0.711)))
    else:
        AAR = -20.772 / np.tan(np.deg2rad(SEA))
    AAR /= 3600

    SEC = SEA + AAR

    if HA > 0:
        temp = np.rad2deg(
            np.arccos(((np.sin(np.deg2rad(lat))*np.cos(np.deg2rad(SZA)))-np.sin(np.deg2rad(SD)))
                      / (np.cos(np.deg2rad(lat))*np.sin(np.deg2rad(SZA))))
        ) + 180
        SAA = np.mod(temp, 360)
    else:
        temp = 540 - np.rad2deg(
            np.arccos(((np.sin(np.deg2rad(lat))*np.cos(np.deg2rad(SZA)))-np.sin(np.deg2rad(SD)))
                      / (np.cos(np.deg2rad(lat))*np.sin(np.deg2rad(SZA))))
        )
        SAA = np.mod(temp, 360)
    
    return SEC, SAA