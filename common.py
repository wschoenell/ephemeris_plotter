import ConfigParser
import sys


# Read Configuration
def load_configuration():
    c = ConfigParser.ConfigParser()
    if len(sys.argv) > 1:
        c.read(sys.argv[1])
    else:
        c.read("config.ini")

    return c