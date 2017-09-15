import datetime

import yaml
from astropy.time import Time
from matplotlib import ticker

from common import load_configuration

__author__ = 'william'

import ephem
import ConfigParser
import numpy as np
import matplotlib.pyplot as plt

try:
    from astroquery.simbad import Simbad
except:
    Simbad = None

fig = plt.figure(1, figsize=(14, 8))
plt.clf()  # Only for interactive mode.
ax = plt.axes()


class ObjectTooLowException(Exception):
    pass


class ProgramExpiredException(Exception):
    pass


def read_tle(fname):
    f = open(fname, 'r')
    tlines = []
    for line in f:
        if line.strip() != '' and line[0] != '#':
            tlines.append(line)
    tle = ephem.readtle(tlines[0], tlines[1], tlines[2])
    f.close()
    return tle


def calc_ephem(object, aux_i, tag=None, duration=None, start=None):
    # Calculate altitudes for the night
    aux = []
    moondist = []
    for dt in alt_dates:
        observer.date = dt
        object.compute(observer)
        aux.append(object.alt * 57.2957795)
        # Calulate the moon separation
        moondist.append(ephem.separation(moon, object) * 57.2957795)

    line, = ax.plot(alt_dates, aux, label=aux_i)
    if duration is not None:
        if start is not None:
            # 1.1574074074074073e-05 = 1. day /24 hours /3600. secs
            t = np.linspace(t0, t0 + duration * 1.1574074074074073e-05, 20)
            aux2 = np.zeros_like(t)
            for i, dt in enumerate(t):
                observer.date = dt
                object.compute(observer)
                aux2[i] = object.alt * 57.2957795
            # observer.date = (sunrise - sunset) / 2.
            ylim = plt.ylim()

            plt.fill_between(t, aux2, ylim[0], color=line.get_color(), alpha=.3)  # , color='black')
            plt.ylim(ylim)

    # Fill in RED zones where moondist is too close
    try:
        moondist_max = c.getfloat("observatory", "max_moondist")
        mask = np.array(moondist) < moondist_max
        if mask.sum() > 0:
            print 'moon too close for %s' % tag
            ylim = plt.ylim()
            plt.fill_between(np.array(alt_dates)[mask], np.array(aux)[mask], ylim[0], color='red',
                             alpha=0.3)  # , alpha=.3)  # , color='black')
    except ConfigParser.NoSectionError:
        pass

    if tag is not None:
        if start is None:
            aux_arg = np.argmax(aux)
            x, y = alt_dates[aux_arg], aux[aux_arg] + .5
        else:
            x, y = start, np.interp(start, alt_dates, aux)
        plt.plot(x, y, '.', color='black')
        plt.text(x, y, tag)
    if output_chimera:
        while True:
            for exptime in [30]:  # (2, 5, 10, 15, 30, 60):
                expnumber = 5
                observer.date = ephem.now()
                object.compute(observer)
                moon.compute(observer)
                # cmd = "ssh lna chimera-tel --slew --ra %s --dec %s && ssh lna chimera-cam --binning 2x2 --object '%s' --output %s -t %i -n %i # mag = %s moondist = %s" % (
                #    object.ra, object.dec, aux_i, aux_i.replace(' ', '_'), exptime, expnumber, object.mag, ephem.separation(moon, object))
                # print cmd
                # os.system(cmd)
                # os.system("say done")
                # time.sleep(1)
                # for i in range(5):
                # print("\a")
                # sleep(.2)
    if output_ingstaralt:
        print "%s %s %s" % (aux_i.replace(" ", "_"), object.ra, object.dec)
    if output_chimera_sched:
        print '%s %s J2000 OBJECT "%s" NREP*(FILTER:TEXP)' % (object.ra.__str__(), object.dec.__str__(), aux_i)


def format_date(date, pos):
    return ephem.date(date).__str__().split(' ')[1].rpartition(':')[0]


c = load_configuration()
output_chimera = c.getboolean("output", "output_chimera")
output_ingstaralt = c.getboolean("output", "output_ingstaralt")
output_chimera_sched = c.getboolean("output", "output_chimera_sched")
try:
    tle_list = c.get("satellite_ephemeris", "tle_list").split(',')
except ConfigParser.NoSectionError, ConfigParser.NoOptionError:
    tle_list = None

# Read Ephemeris file
objects_id = []
objects_ephem = []

try:
    with open("%s" % c.get("ephemeris", "data_file")) as f:
        for line in f:
            if not line.startswith('#'):
                objects_id.append(line.split(',')[0])
                objects_ephem.append(line)
except ConfigParser.NoSectionError:
    pass
except TypeError:
    pass

# Compute Observatory
observer = ephem.Observer()
observer.lat = c.get("observatory", "latitude")
observer.lon = c.get("observatory", "longitude")
observer.elevation = c.getfloat("observatory", "elevation")
# Compute next sunset and next sunrise
time_offset = 1 if (0 <= datetime.datetime.now().hour < 8) else 0
sunset = observer.next_setting(ephem.Sun(), start=datetime.date.today()) - time_offset
observer.date = sunset
sunrise = observer.next_rising(ephem.Sun())
alt_dates = np.linspace(sunset - 0.015, sunrise + 0.015, 100)  # For the graph

# Compute Ephemeris
alt = []  # For object altitudes

# Moon
# Load ephem of the desired object
moon = ephem.Moon()
# Calculate altitudes for the night
aux = []
for dt in alt_dates:
    observer.date = dt
    moon.compute(observer)
    aux.append((moon.alt / np.pi) * 180)
ax.plot(alt_dates, aux, '--', label="Moon")

# Minimum alt
try:
    ax.hlines(c.getfloat("observatory", "min_alt"), ax.get_xlim()[0], ax.get_xlim()[1], color='purple')
except ConfigParser.NoSectionError:
    pass

# RA, DEC file
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

# Satellite "TLE" ephemeris files
if tle_list is not None:
    for tle in tle_list:
        tle_file = "%s/%s" % (c.get("satellite_ephemeris", "tle_dir"), tle)
        tle_data = read_tle("%s.tle" % tle_file)
        calc_ephem(tle_data, "tle_%s" % tle)

t = 0
i_sched = 0
object = None
now = ephem.now()
t0 = float(sunset) if now < sunset else (float(now) - time_offset)
t_start_night = t0
try:
    with open(c.get("ephemeris", "chimera-sched_file")) as f:
        data = yaml.load(f)
    for prog in data['programs']:

        try:

            if 'startAt' in prog and 'validFor' in prog:
                start_at = float(ephem.Date(Time(float(prog['startAt']), format='mjd').datetime))
                t0 = start_at
                if start_at + float(prog['validFor']) * ephem.second < t0:
                    raise ProgramExpiredException
                print "Start at set to ", start_at, " with t0 ", t0
            else:
                start_at = t0

            t = 0


            for action in prog['actions']:


                if action['action'] == 'point':

                    if 'ra' in action:
                        object = ephem.FixedBody()
                        object._ra = ephem.hours(action['ra'])
                        object._dec = ephem.degrees(action['dec'])
                    elif 'name' in action and not "alt" in action:
                        if Simbad is not None:
                            target = Simbad.query_object(action["name"])[0]
                            object = ephem.FixedBody()
                            object._ra = ephem.hours(str(target['RA']))
                            object._dec = ephem.degrees(str(target['DEC']))
                        else:
                            print 'Skipping %i. This script does not support name searches yet.' % i_sched

                    if 'name' in action:
                        object.name = action["name"]

                    if object is not None:
                        # Calulate t0 for min_alt

                        try:
                            min_alt = c.getfloat("observatory", "min_alt")
                            observer.date = start_at
                            object.compute(observer)
                            aux_t0 = start_at
                            while object.alt * 57.2957795 < min_alt:
                                observer.date = aux_t0
                                object.compute(observer)
                                aux_t0 += 1e-2

                            # If the next time object is above min_alt is after sunset, break the loop and go to next object
                            if aux_t0 < sunrise:
                                start_at = aux_t0
                            else:
                                title = 'scheduler_%i' % i_sched if len(object.name) == 0 else object.name + " (%i)" % (
                                    i_sched + 1)
                                calc_ephem(object, title)
                                i_sched += 1
                                raise ObjectTooLowException

                        except ConfigParser.NoSectionError:
                            print "WARN: no min_alt configured!"
                            pass

                        if start_at is not None:
                            start_at += (t * 1.1574074074074073e-05)


                elif action['action'] == 'expose':
                    t += action['exptime'] * action['frames']
                    try:
                        t += c.getfloat("telescope_times", "exposure_overhead") * action['frames']
                        t += c.getfloat("telescope_times", "point_overhead")
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
                elif action['action'] == 'pointverify':
                    t += float(c.get("telescope_times", "pverify_overhead"))
            if 'startAt' in prog:
                x = float(ephem.Date(Time(float(prog['startAt']), format='mjd').datetime))
                if x > t_start_night:
                    if start_at < x:
                        start_at = x
                else:
                    observer.date = ephem.Date(Time(float(prog['startAt']), format='mjd').datetime)
                    object.compute(observer)
                    aux_min_alt = object.alt
                    while object.alt < aux_min_alt:
                        observer.date = start_at
                        object.compute(observer)
                        start_at += 1e-2
                    print("start_at < night start for object %s. Correct it with:\n        startAt: %f" % (
                        object.name, Time(ephem.Date(start_at).datetime(), format="datetime").mjd))
            object_name = object.name if object is not None else "None"
            title = 'scheduler_%i' % i_sched if len(object_name) == 0 else object_name + " (%i)" % (i_sched + 1)
            if t > 0:
                calc_ephem(object, title, tag='%s' % (i_sched + 1), duration=t, start=start_at)
            print title, Time(ephem.Date(start_at).datetime(), format="datetime").mjd, t
            t0 = ephem.Date(ephem.Date(t0).datetime() + datetime.timedelta(seconds=t))
            i_sched += 1

        except ObjectTooLowException:
            print "Object %i is too low." % i_sched
        except ProgramExpiredException:
            print "Object %i EXPIRED." % i_sched


            # if not plot:
            #     calc_ephem(object, 'scheduler_%i' % i_sched, tag='%s' % (i_sched + 1), duration=t, start=t0)
            #     t0 += (t * 1.1574074074074073e-05)
except ConfigParser.NoOptionError:
    pass

# for id_object in sys.argv[1:]:
#     aux_i = np.argwhere(np.array([a.startswith(id_object) for a in objects_id]))
#     if len(aux_i) == 0:
#         print 'Couldn\'t find object on ephem database: %s' % id_object
#     else:
#         # Load ephem of the desired object
#         object = ephem.readdb(objects_ephem[aux_i[0]])
#         calc_ephem(object, objects_id[aux_i])

ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_date))

ax.legend(loc='upper center', ncol=5, fancybox=True, framealpha=0.5)
ylim = ax.get_ylim()
ax.set_ylim(10, 90)
ax.set_xlim(sunset - 0.015, sunrise + 0.015)
ax.set_xlabel('UTC time')
ax.set_ylabel('Altitude (deg)')
plt.title("%i / %s @ latitude: %s - longitude: %s" % (
    np.ceil(float(sunset.__str__().split(' ')[0])), c.get("observatory", "name"), observer.lat, observer.long))
plt.grid()
ymin, ymax = plt.ylim()
plt.vlines(sunset, ymin, ymax, color='gray')
plt.vlines(sunrise, ymin, ymax, color='gray')
plt.vlines(ephem.now(), ymin, ymax, color='black', linestyles='dashed')

try:
    plt.savefig(c.get('plot', 'outfile'))
except ConfigParser.NoSectionError, ConfigParser.NoOptionError:
    pass
