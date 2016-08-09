import datetime
import yaml
from astropy.time import Time
from matplotlib import ticker

__author__ = 'william'

import ephem
import ConfigParser
import sys
import numpy as np
import matplotlib.pyplot as plt

try:
    from astroquery.simbad import Simbad
except:
    Simbad = None

plt.clf()  # Only for interactive mode.
fig = plt.figure(1, figsize=(16.5, 8))
ax = plt.axes()


# fig.autofmt_xdate()

def calc_ephem(object, aux_i, tag=None, duration=None, start=None):
    # Calculate altitudes for the night
    aux = []
    for dt in alt_dates:
        observer.date = dt
        object.compute(observer)
        aux.append(object.alt * 57.2957795)
    ax.plot(alt_dates, aux, label=aux_i)
    if duration is not None:
        print 'Duration: %3.2f' % duration
        if start is not None:
            print 'start at %f' % start
            # 1.1574074074074073e-05 = 1. day /24 hours /3600. secs
            t = np.linspace(t0, t0 + duration * 1.1574074074074073e-05, 20)
            aux2 = []
            for dt in t:
                observer.date = dt
                object.compute(observer)
                aux.append(object.alt * 57.2957795)
            # observer.date = (sunrise - sunset) / 2.
            plt.fill_between(t, aux2, color='black', alpha=.3)  # , color='black')
    if tag is not None:
        if start is None:
            aux_arg = np.argmax(aux)
            x, y = alt_dates[aux_arg], aux[aux_arg] + .5
        else:
            x, y = start, np.interp(start, alt_dates, aux)
        plt.plot(x, y, '.', color='black')
        plt.text(x, y, tag)
    if output_chimera:
        print "chimera-tel --slew --ra %s --dec %s && chimera-cam --object '%s' -t NN # mag = %s" % (
            object.ra, object.dec, aux_i, object.mag)
    if output_ingstaralt:
        print "%s %s %s" % (aux_i.replace(" ", "_"), object.ra, object.dec)
    if output_chimera_sched:
        print '%s %s J2000 OBJECT "%s" NREP*(FILTER:TEXP)' % (object.ra.__str__(), object.dec.__str__(), aux_i)


def format_date(date, pos):
    return ephem.date(date).__str__().split(' ')[1].rpartition(':')[0]


# Read Configuration
c = ConfigParser.ConfigParser()
c.read("config.ini")
output_chimera = c.getboolean("output", "output_chimera")
output_ingstaralt = c.getboolean("output", "output_ingstaralt")
output_chimera_sched = c.getboolean("output", "output_chimera_sched")

# Read Ephemeris file
objects_id = []
objects_ephem = []

with open("%s" % c.get("ephemeris", "data_file")) as f:
    for line in f:
        # try:
        if not line.startswith('#'):
            objects_id.append(line.split(',')[0])
            objects_ephem.append(line)
            # except:
            #     pass

# Compute Observatory
observer = ephem.Observer()
observer.lat = c.get("observatory", "latitude")
observer.lon = c.get("observatory", "longitude")
observer.elevation = c.getfloat("observatory", "elevation")
# Compute next sunset and next sunrise
sunset = observer.next_setting(ephem.Sun(), start=datetime.date.today())
observer.date = sunset
sunrise = observer.next_rising(ephem.Sun())
alt_dates = np.linspace(sunset, sunrise, 100)  # For the graph

# Compute Ephemeris
alt = []  # For object altitudes

# Moon
# Load ephem of the desired object
object = ephem.Moon()
# Calculate altitudes for the night
aux = []
for dt in alt_dates:
    observer.date = dt
    object.compute(observer)
    aux.append((object.alt / np.pi) * 180)
ax.plot(alt_dates, aux, '--', label="Moon")

try:
    f = open(c.get("ephemeris", "radec_file"))
    for l in f.readlines():
        if not l.startswith('#'):
            d = l.replace('\n', '').split(',')
            # Calculate coordinates
            object = ephem.FixedBody()
            ra, dec = d[1:3]
            if ':' in ra:
                object._ra = ephem.hours(ra)
            else:
                object._ra = ephem.degrees(ra)
            object._dec = ephem.degrees(dec)
            calc_ephem(object, d[0])
except:  # ConfigParser.NoOptionError:
    pass

t = 0
i_sched = 0
object = None
try:
    with open(c.get("ephemeris", "chimera-sched_file")) as f:
        data = yaml.load(f)
    for prog in data['programs']:
        if 'slewAt' in prog:
            x = ephem.Date(Time(float(prog['slewAt']), format='mjd').datetime)
            t0 = x if x > sunset else None
            print x, t0, sunset
        else:
            t0 = None
        for action in prog['actions']:
            if action['action'] == 'point':
                if object is not None:
                    if t0 is not None:
                        t0 += t * 1.1574074074074073e-05
                    calc_ephem(object, 'scheduler_%i' % i_sched, tag='%s' % (i_sched + 1), duration=t, start=t0)
                t = 0
                if 'ra' in action:
                    object = ephem.FixedBody()
                    object._ra = ephem.hours(action['ra'])
                    object._dec = ephem.degrees(action['dec'])
                else:
                    if Simbad is not None:
                        pass
                    else:
                        print 'Skipping %i. This script does not support name searches yet.' % i_sched
                i_sched += 1
            elif action['action'] == 'expose':
                t += action['exptime']
                try:
                    t += float(c.get("telescope_times", "exposure_overhead"))
                    t += float(c.get("telescope_times", "point_overhead"))
                except:
                    pass
            elif action['action'] == 'autofocus':
                try:
                    if 'exptime' in action:
                        t += float(c.get("telescope_times", "autofocus_fit"))
                    else:
                        t += float(c.get("telescope_times", "autofocus_model"))
                except:
                    pass
except ConfigParser.NoOptionError:
    pass

for id_object in sys.argv[1:]:
    aux_i = np.argwhere(np.array([a.startswith(id_object) for a in objects_id]))
    if len(aux_i) == 0:
        print 'Couldn\'t find object on ephem database: %s' % id_object
    else:
        # Load ephem of the desired object
        object = ephem.readdb(objects_ephem[aux_i[0]])
        calc_ephem(object, objects_id[aux_i])

ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_date))

ax.legend(loc='upper center', ncol=5, fancybox=True, framealpha=0.5)
ylim = ax.get_ylim()
ax.set_ylim(0, 110)
ax.set_xlim(sunset, sunrise)
ax.set_xlabel('UTC time')
ax.set_ylabel('Altitude (deg)')
plt.title("%s / %s @ latitude: %s - longitude: %s" % (
    sunset.__str__().split(' ')[0], c.get("observatory", "name"), observer.lat, observer.long))
plt.grid()

try:
    plt.savefig(c.get('plot', 'outfile'))
except ConfigParser.NoSectionError, ConfigParser.NoOptionError:
    pass
