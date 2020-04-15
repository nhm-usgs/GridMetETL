"""Console script for gridmetetl."""
from gridmetetl.etl import FpoNHM
import argparse
import sys
import datetime


def valid_date(s):
    try:
        return datetime.datetime.strptime(s, "%Y-%m-%d")
    except ValueError:
        msg = "Not a valid date: '{0}'.".format(s)
        raise argparse.ArgumentTypeError(msg)


def main():
    """Console script for gridmetetl."""
    numdays = None
    startdate = None
    enddate = None
    idir = None
    odir = None
    wght_file = None
    extract_type = None
    file_prefix = None
    gm_vars = None

    my_parser = argparse.ArgumentParser(prog='gridmet_etl',
                                        description='map gridded climate data to polygon using zonal area weighted mean')

    my_parser.add_argument('-t', '--extract_type', type=str,
                           help='extract method: (days) or (date)', metavar='extraction type',
                           default=None, required=True, choices=['days', 'date'])

    my_parser.add_argument('-p', '--period', type=valid_date,
                           help='option: start date and end date of retrieval (YYYY-MM-DD)',
                           metavar='(YYYY-MM-DD)',
                           default=None,
                           nargs=2)

    my_parser.add_argument('-d', '--days', type=int,
                           help='option: number of days to retrieve; if specified take precedence over -s & -e option',
                           metavar='numdays', default=None)

    my_parser.add_argument('-f', '--file_prefix', type=str,
                           help='option: prefix for output files',
                           metavar='output_file_prefix', default='')

    my_parser.add_argument('-i', '--inpath', type=str,
                           help='input_path (location of HRU shapefiles)', metavar='input_path',
                           default=None, required=True)

    my_parser.add_argument('-o', '--outpath', type=str,
                           help='Output path (location of netcdf output files by shapefile output)', metavar='output_path',
                           default=None, required=True)

    my_parser.add_argument('-w', '--weightsfile', type=str,
                           help='path/weight.csv - path/name of weight file', metavar='weight_file',
                           default=None, required=True)

    my_parser.add_argument('-v', '--variables', nargs='*', type=str,
                           help='over-ride default vars',
                           choices=['tmax', 'tmin', 'ppt', 'rhmax', 'rhmin', 'ws', 'srad'],
                           metavar='GridMet_Variables',
                           default=['tmax', 'tmin', 'ppt', 'rhmax', 'rhmin', 'ws'])

    args = my_parser.parse_args()

    if all(i is not None for i in [args.period, args.days]):
        my_parser.error('Either the --days or --period option must be specified not both')
    if all(i is None for i in [args.period, args.days]):
        my_parser.error('Either the --days or --period option must be specified')
    # if args.period is not None:
    #     startdate = args.period[0]
    #     enddate = args.period[1]
    if args.extract_type is not None:
        extract_type = args.extract_type
        if args.extract_type == 'days':
            if args.days is not None:
                numdays = args.days
            else:
                my_parser.error('if -t --extract_type == days then -d must be specified')
        elif args.extract_type == 'date':
            if args.period is not None:
                startdate = args.period[0]
                enddate = args.period[1]
                if startdate >= enddate:
                    my_parser.error('when using -p the first date must occur before the second')
            else:
                my_parser.error('if -t --extract_type == dates then -p must be specified')

    if args.outpath is not None:
        odir = args.outpath
    if args.inpath is not None:
        idir = args.inpath
    if args.weightsfile is not None:
        wght_file = args.weightsfile
    if args.file_prefix is not None:
        file_prefix = args.file_prefix
    if args.variables is not None:
        gm_vars = args.variables

    print('starting Script', flush=True)
    # numdays = 2
    fp = FpoNHM()
    print('instantiated', flush=True)
    # initialize(self, iptpath, optpath, weights_file, type=None, days=None, start_date=None, end_date=None)
    # ready = fp.initialize(idir, odir, wght_file, extract_type, numdays, startdate, enddate, file_prefix)
    try:
        ready = fp.initialize(gm_vars, idir, odir, wght_file, etype=extract_type, days=numdays,
                              start_date=startdate, end_date=enddate,
                              fileprefix=file_prefix)
        if ready:
            print('initalized\n', flush=True)
            print('running', flush=True)
            fp.run_weights()
            print('finished running', flush=True)
            fp.finalize()
            print('finalized', flush=True)
            return 0
        else:
            if extract_type == 'days':
                print('Gridmet not updated continue with numdays -1', flush=True)
                fp.setnumdays(numdays - 1)
                print('initalized\n', flush=True)
                print('running', flush=True)
                status = fp.run_weights()
                print('finished running', flush=True)
                fp.finalize()
                print('finalized', flush=True)
                return 0
            else:
                print('error: extract did not return period specified, Gridmet no updated', flush=True)
                return 1

    except NameError:
        print('An error: ')
        raise



if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
