import pyfits
import pylab as pl
from datetime import date

def revise(file,
           freq_start=80e6,
           freq_end=300e6,
           steps=12,
           grid=(20,20),
           rect=(-2., -2, 2., 2.),
           obs=(215., 0., 0.),
           rsph = 25.,
           mode='Tbr',
           cnu=3.0,
           msun=1.0,
           units = 'K'):

    hdulist = pyfits.open(file, mode='update')

    print hdulist.ascardlist()[:]

#    prihdr = hdulist[0].header



    hdulist.flush()
    hdulist.close()

'''


    prihdr.update('OBJECT', 'Sun', 'Source Name')

    prihdr.update('TELESCOP', 'Ray Sim','')
    prihdr.update('INSTRUME', 'HART','')
    prihdr.update('OBSERVER', '','')
    prihdr.update('DATE-OBS', str(date.today()),'')

    prihdr.update('BSCALE', 1.0,'')
    prihdr.update('BZERO', 0.0,'')
    if(units=='K'):
        prihdr.update('BUNIT', 'K', 'Temperature Brightness')
    if(units=='Jy'):
        prihdr.update('BUNIT', 'JY', 'Temperature Brightness')

    prihdr.update('CTYPE1', 'RA---SIN', 'Coordinates along the ecliptic')
    prihdr.update('CRVAL1', 0.0, '') 
    prihdr.update('CDELT1', prihdr['CDELT1']*3600, 'Degrees/Pixel')
    prihdr.update('CRPIX1', (-rect[0])/(rect[2]-rect[0])*grid[0], '')


    prihdr.update('CTYPE2', 'DEC--SIN', 'Coordinates perpendicular to the ecliptic')
    prihdr.update('CRVAL2', 0.0, '')
    prihdr.update('CDELT2', prihdr['CDELT2']*3600, 'Degrees/Pixel')
    prihdr.update('CRPIX2', (-rect[1])/(rect[3]-rect[1])*grid[1], '')


    prihdr.update('CRVAL3', freq_start, 'Frequency in Hz')
    prihdr.update('CDELT3', (freq_end-freq_start)/(steps-1), 'Frequency per index')
    prihdr.update('CRPIX3', 1.0, 'Reference Freq index')

    prihdr.update('CDELT4', 1.0, '')
    prihdr.update('CRPIX4', 1.0, '')
    prihdr.update('CROTA4', 0.0, '')

    prihdr.add_comment('Scattering on')
    print prihdr['SIMPLE']

'''




        





revise(file='data/theta45_phi30_I.fits')


