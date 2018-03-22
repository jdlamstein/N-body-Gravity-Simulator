# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 03:05:45 2016

@author: joshualamstein
"""
import numpy as np
_all_=['x','y','z','vx','vy','vz','m']

# Data on celestial bodies in the Solar system found on nasa.gov. 
#==============================================================================
# 1988 - Jun - 20


solar = 1.988544*10**30 # in kg
mercury = 3.302e23
venus = 48.685e23
earth = 5.97219e24
mars = 6.4185e23
jup = 1898.13*10**24
saturn = 5.68319e26
neptune = 102.41e24
uranus = 86.8103e24
pluto= 1.307e22

planet = 3.2e24 # Planet refers to an additional planet I can add to the solar
# system to test if it is maximally packed. (The solar system is not
# maximally packed. There is room to add another planet and preserve the 
# stable orbit of all planets.)

dayToYear = 365.25; 


#==============================================================================
# Mercury
#==============================================================================
   
xMercury = 8.887985138765460E-02
yMercury = -4.426150338141062E-01
zMercury = -4.475716356484761E-02
vxMercury = 2.190877912081542E-02
vyMercury = 7.161568136528000E-03
vzMercury =  -1.425929443086507E-03
mMercury = mercury / solar

#==============================================================================
# Venus
#==============================================================================

xVenus = 4.043738093622098E-02 
yVenus = -7.239789211502183E-01
zVenus = -1.241560658530024E-02
vxVenus = 2.005742309538389E-02
vyVenus = 1.141448268256643E-03
vzVenus = -1.142174441569258E-03
mVenus = venus/solar
  
#==============================================================================
# Earth
#==============================================================================


xEarth =  -2.020844529756663E-02
yEarth =  -1.014332737790859E+00
zEarth =  -1.358267619371298E-05
#velocity in AU/day
vxEarth =  1.692836723212859E-02
vyEarth = -3.484006532982474E-04 
vzEarth =  6.028542314557626E-07
mEarth = earth/solar;

#target 1998
#  -3.466334931755365E-02 -1.013773181327570E+00  2.111689861662178E-04
#   1.692129864984556E-02 -5.252811129268817E-04  3.987686870581435E-07


#==============================================================================
# Mars
#==============================================================================

   
xMars = 7.462481663749645E-01
yMars = -1.181663652521456E+00
zMars = -4.321921404013512E-02
vxMars = 1.235610918162121E-02
vyMars = 8.680869489377649E-03 
vzMars = -1.220500608452554E-04
mMars = mars/solar


#==============================================================================
# Planet X
#==============================================================================

xPlanet = 0
yPlanet = 2.06
zPlanet = 0
vxPlanet = 1.235610918162121E-02
vyPlanet = 0
vzPlanet = 0
mPlanet = planet / solar


#==============================================================================
# Jupiter
#==============================================================================

xJup = 3.384805319103406E+00  
yJup =   3.658805636759595E+00
zJup =   -9.100441946210819E-02 
#velocity in AU/day
vxJup =  -5.634671617093230E-03
vyJup =  5.479180979634376E-03 
vzJup =    1.034981407898108E-04
mJup = jup/solar;

#==============================================================================
# Saturn
#==============================================================================

xSaturn = -1.083899692644216E-01
ySaturn = -1.003995196286016E+01
zSaturn =  1.793391553155583E-01
vxSaturn = 5.278410787728323E-03
vySaturn = -7.712342079566598E-05
vzSaturn =  -2.084447335785041E-04
mSaturn = saturn / solar

#==============================================================================
# Neptune
#==============================================================================

xNeptune = 4.675566709791660E+00
yNeptune = -2.985428200863175E+01
zNeptune = 5.070034142531887E-01
vxNeptune = 3.080716380724798E-03
vyNeptune = 5.030733458293977E-04
vzNeptune = -8.101711269674541E-05
mNeptune = neptune / solar

#==============================================================================
# Uranus
#==============================================================================

xUranus = -2.693448460292631E-01
yUranus = -1.927606446869220E+01
zUranus =  -6.808868692550485E-02
vxUranus =  3.903100242621723E-03
vyUranus = -2.380111092360100E-04
vzUranus =  -5.164025224695875E-05
mUranus = uranus / solar

#==============================================================================
# Pluto 
#==============================================================================

xPluto = -2.129074273328636E+01
yPluto = -1.896633337434039E+01
zPluto = 8.187955378677129E+00
vxPluto = 2.276295756013608E-03
vyPluto = -2.670481848836963E-03
vzPluto = -3.669545371032554E-04
mPluto = pluto / solar 


#target  -1.156541154581570E+01 -2.704864218000164E+01  6.239749761161465E+00
#   2.964408290188142E-03 -1.722224413824548E-03 -6.839434010481107E-04

#==============================================================================
# Sun 
#==============================================================================
       
xSun = -3.430031536367300E-03
ySun = 1.761881027012596E-03
zSun = 1.246691303879918E-05
vxSun =  3.433119412673547E-06
vySun =  -5.231300927361546E-06
vzSun = -2.972974735550750E-08
mSun = 1


#==============================================================================
# arrays
#==============================================================================


#Arrays for modeling the solar system
xSet = np.array([[xSun], [xMercury], [xVenus], [xEarth], [xMars], [xJup], [xSaturn], [xUranus], [xNeptune], [xPluto]])
ySet = np.array([[ySun], [yMercury], [yVenus], [yEarth], [yMars], [yJup], [ySaturn], [yUranus], [yNeptune], [yPluto]])
zSet = np.array([[zSun], [zMercury], [zVenus], [zEarth], [zMars], [zJup], [zSaturn], [zUranus], [zNeptune], [zPluto]])
vxSet = np.array([[vxSun], [vxMercury], [vxVenus], [vxEarth], [vxMars], [vxJup], [vxSaturn], [vxUranus], [vxNeptune], [vxPluto]])
vySet = np.array([[vySun], [vyMercury], [vyVenus], [vyEarth], [vyMars], [vyJup], [vySaturn], [vyUranus], [vyNeptune], [vyPluto]])
vzSet = np.array([[vzSun], [vzMercury], [vzVenus], [vzEarth], [vzMars], [vzJup], [vzSaturn], [vzUranus], [vzNeptune], [vzPluto]])
m = np.array([mSun, mMercury, mVenus, mEarth, mMars, mJup, mSaturn, mUranus, mNeptune, mPluto])

#Arrays for testing if the solar system is maximally packed. 
xPacked = np.array([[xSun], [xMercury], [xVenus], [xEarth], [xMars], [xPlanet],[xJup], [xSaturn], [xUranus], [xNeptune], [xPluto]])
yPacked = np.array([[ySun], [yMercury], [yVenus], [yEarth], [yMars], [yPlanet], [yJup], [ySaturn], [yUranus], [yNeptune], [yPluto]])
zPacked = np.array([[zSun], [zMercury], [zVenus], [zEarth], [zMars], [zPlanet],[zJup], [zSaturn], [zUranus], [zNeptune], [zPluto]])
vxPacked = np.array([[vxSun], [vxMercury], [vxVenus], [vxEarth], [vxMars], [vxPlanet],[vxJup], [vxSaturn], [vxUranus], [vxNeptune], [vxPluto]])
vyPacked = np.array([[vySun], [vyMercury], [vyVenus], [vyEarth], [vyMars], [vyPlanet],[vyJup], [vySaturn], [vyUranus], [vyNeptune], [vyPluto]])
vzPacked = np.array([[vzSun], [vzMercury], [vzVenus], [vzEarth], [vzMars], [vzPlanet],[vzJup], [vzSaturn], [vzUranus], [vzNeptune], [vzPluto]])
mPacked = np.array([mSun, mMercury, mVenus, mEarth, mMars, mPlanet, mJup, mSaturn, mUranus, mNeptune, mPluto])





