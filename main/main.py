# constants/geometry
import sys
import argparse
import logging

from aileron import aileron


aileron_o = aileron()
logger = logging.getLogger()

def main(args):
    
    parser = argparse.ArgumentParser(description="SVV Program 2018 A32")
    parser.add_argument("-o", "--out", help="Specify output file", default=sys.stdout, type=argparse.FileType('w'))
    parser.add_argument("-v", "--verbose", help="Give more updates during runnning", action='count')
    
    arguments = parser.parse_args(args)
    logger.setLevel(arguments.verbose*10)

    context = (logger, aileron)
    



    
    pass


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))