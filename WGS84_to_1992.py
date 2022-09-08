
''' 
Wersja oryginalna (MatLab):
Konwersja wspó³rzêdnych uk³adu WGS84 do uk³adu PUWG92
    Opis:
    konwersja wspolrzednych z ukladu WGS 84 do ukladu PUWG 1992
    Parametry:
    lon - szerokosc geograficzna wyrazona w stopniach
    lat - dlugosc geograficzna wyrazona w stopniach
    Wartosc zwracana
    Xpuwg - wspolrzedna X ukladu PUWG 1992 (UWAGA - wspolrzedna pionowa) 
            lub NaN, jesli parametry wejsciowe sa poza zakresem
    Ypuwg - wspolrzedna Y ukladu PUWG 1992 (UWAGA - wspolrzedna pozioma)
            lub NaN, jesli parametry wejsciowe sa poza zakresem
            
  Autor: Zbigniew Szymanski
  E-mail: z.szymanski@szymanski-net.eu
  Wersja: 1.0
  Data modyfikacji: 2012-11-17
  Uwagi: Oprogramowanie darmowe. Dozwolone jest wykorzystanie i modyfikacja 
        niniejszego oprogramowania do wlasnych celow pod warunkiem 
        pozostawienia wszystkich informacji z naglowka. W przypadku 
        wykorzystania niniejszego oprogramowania we wszelkich projektach
        naukowo-badawczych, rozwojowych, wdrozeniowych i dydaktycznych prosze
       o zacytowanie nastepujacego artykulu:
       
        Zbigniew Szymanski, Stanislaw Jankowski, Jan Szczyrek, 
        "Reconstruction of environment model by using radar vector field histograms.",
        Photonics Applications in Astronomy, Communications, Industry, and 
        High-Energy Physics Experiments 2012, Proc. of SPIE Vol. 8454, pp. 845422 - 1-8,
        doi:10.1117/12.2001354
       
 Literatura:
        Uriasz, J., “Wybrane odwzorowania kartograficzne”, Akademia Morska w Szczecinie,
        http://uriasz.am.szczecin.pl/naw_bezp/odwzorowania.html
        
        
    Wersja dla jêzyka Python:
    '''
import numpy as np

def WGS84_to_1992(lon:float, lat:float):
    """Convert geographic coordinates from WGS84 to 1992

    Arguments:
        lon {float} -- longitude in WGS84
        lat {float} -- latitude in WGS84

    Returns:
        Xpuwg -- longitude in WGS84
        Ypuwg -- latitude in WGS84
    """    
    #Parametry elipsoidy GRS-80
    e=0.0818191910428  #pierwszymimoœród elipsoidy
    R0=6367449.14577 #promieñ sfery Lagrange’a
    Snorm=2.0E-6   #parametr normuj¹cy
    xo=5760000; #parametr centruj¹cy
    a0=5765181.11148097
    a1=499800.81713800
    a2=-63.81145283
    a3=0.83537915
    a4=0.13046891
    a5=-0.00111138
    a6=-0.00010504

    #Parametry odwzorowania Gaussa-Kruegera dla uk³adu PUWG92
    L0_stopnie=19; #Pocz¹tek uk³adu wsp. PUWG92 (d³ugoœæ)
    m0=0.9993
    x0=-5300000
    y0= 500000

    #Zakres stosowalnosci metody
    Bmin=48*np.pi/180
    Bmax=56*np.pi/180
    dLmin=-6*np.pi/180
    dLmax=6*np.pi/180

    # Dane wejœciowe
    B=lon*np.pi/180
    dL_stopnie=lat-L0_stopnie
    dL=dL_stopnie*np.pi/180

    # etap I - elipsoida na kulê Lagrange’a
    U=1-e*np.sin(B)
    V=1+e*np.sin(B)
    K=(U/V)**(e/2)
    C=K*np.tan(B/2+np.pi/4)
    fi=2*np.arctan(C)-np.pi/2
    d_lambda=dL

    # etap II - kula na walec Merkatora
    p=np.sin(fi)
    q=np.cos(fi)*np.cos(d_lambda)
    r=1+np.cos(fi)*np.sin(d_lambda)
    s=1-np.cos(fi)*np.sin(d_lambda)
    XMERC=np.array(R0*np.arctan(p/q))
    YMERC=np.array(0.5*R0*np.log(r/s))
    

    # etap III - walec na p³aszczyznê
    Z=np.vectorize(complex)((XMERC-xo)*Snorm, YMERC*Snorm)
    Zgk=a0+Z*(a1+Z*(a2+Z*(a3+Z*(a4+Z*(a5+Z*a6)))))
    Xgk=np.real(Zgk)
    Ygk=np.imag(Zgk)

    # Przejœcie do uk³adu aplikacyjnego
    Xpuwg=m0*Xgk+x0
    Ypuwg=m0*Ygk+y0

    return Xpuwg, Ypuwg