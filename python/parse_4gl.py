import xml.etree.ElementTree as ET
from astropy.coordinates import SkyCoord
import astropy.units as u

# Load the XML file
tree = ET.parse('/Users/richarddubois/Code/GLAST/tmp/gll_psc_v32.xml')
root = tree.getroot()

# Specify the source name you want to find
target_source_name = "4FGL J0240.5+6113"

# Use XPath to find the specific source by its name

sources = tree.findall("source")
for i_s, s in enumerate(sources):
    source_name = s.attrib["name"]
    if source_name == target_source_name:
        source_Energy_Flux100 = s.attrib["Energy_Flux100"]
        print(source_name, source_Energy_Flux100, i_s)

        spatial = s.find("spatialModel")
        parameters = spatial.findall("parameter")
        where_RADEC = {}
        where_gal = {}
        for p in parameters:
            name = p.attrib["name"]
            value = p.attrib["value"]
            where_RADEC[name] = float(value)

        print(where_RADEC)
        # Create a SkyCoord object in equatorial coordinates
        coords = SkyCoord(ra=where_RADEC["RA"] * u.deg, dec=where_RADEC["DEC"] * u.deg, frame='icrs')
        # Convert to Galactic coordinates
        galactic_coords = coords.galactic
        where_gal["lon"] = galactic_coords.l
        where_gal["lat"] = galactic_coords.b

        break

nearby_sources = []

for n_s, s in enumerate(sources):
    if n_s == i_s:  # bypass LS I
        continue

    source_name = s.attrib["name"]

    type = s.attrib["type"]
    if type != "PointSource":
        continue

    source_Energy_Flux100 = s.attrib["Energy_Flux100"]

    spatial = s.find("spatialModel")
    parameters = spatial.findall("parameter")
    where_RADEC_n = {}
    where_gal_n = {}
    for p in parameters:
        name = p.attrib["name"]
        value = p.attrib["value"]
        where_RADEC_n[name] = float(value)

    # Create a SkyCoord object in equatorial coordinates
    try:
        coords = SkyCoord(ra=where_RADEC_n["RA"] * u.deg, dec=where_RADEC_n["DEC"] * u.deg, frame='icrs')
    except KeyError:
        print("coords exception", source_name, "RA-DEC", where_RADEC_n)
        exit()
    # Convert to Galactic coordinates
    galactic_coords = coords.galactic
    where_gal_n["lon"] = galactic_coords.l
    where_gal_n["lat"] = galactic_coords.b

    lon_diff = abs(where_gal["lon"] - where_gal_n["lon"])
    lat_diff = abs(where_gal["lat"] - where_gal_n["lat"])

    if (lon_diff > 30.*u.deg and lon_diff < 35.*u.deg) and lat_diff < 5.*u.deg:
        n_info = [source_name, n_s, source_Energy_Flux100, where_RADEC_n]
        nearby_sources.append(n_info)

print(len(nearby_sources))
for res in nearby_sources:
    print(res)
