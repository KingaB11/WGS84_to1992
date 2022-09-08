
## Authors

- [@KingaB11](https://www.github.com/KingaB11)


## API Documentation

    Convert geographic coordinates from WGS84 to 1992 for the area of Poland

    Arguments:
        lon {float} -- longitude in WGS84
        lat {float} -- latitude in WGS84

    Returns:
        Xpuwg -- longitude in 1992
        Ypuwg -- latitude in 1992

## Usage/Examples

    In: 
    lon = np.array([51.7368, 53.0695, 52.3562])
    lat = np.array([14.3313, 14.2582, 14.5412])
    
    lat, lon = WGS84_to_1992(lon, lat)
    
    Out:
    lat = [177794.81116683 182483.48202851 196505.35939497] 
    lon = [440362.15520962 588750.26298626 508275.0857071]
## Acknowledgements

Originally written by Zbigniew Szymanski (reference: Zbigniew Szymanski, Stanislaw Jankowski, Jan Szczyrek, 
"Reconstruction of environment model by using radar vector field histograms.",
Photonics Applications in Astronomy, Communications, Industry, and 
High-Energy Physics Experiments 2012, Proc. of SPIE Vol. 8454, pp. 845422 - 1-8,
doi:10.1117/12.2001354) in MatLab

Presented version works with Python 3.x

