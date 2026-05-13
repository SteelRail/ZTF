#!/usr/bin/env python

import os
import sys
import argparse
from time import time
from pathlib import Path
from ztfquery import query, io, buildurl


def load_query(ra, dec, size, max_seeing, root):
    ra, dec, max_seeing = round(ra, 3), round(dec, 3), round(max_seeing, 1)
    query_hash = f'ra{ra}_dec{dec}_size{size}_seeing{max_seeing}'
    query_csv = root / 'query' / f'{query_hash}.csv'
    if query_csv.exists():
        print(f'Query {query_hash} already exists. Loading from disk.')
        zquery = query.ZTFQuery.from_metafile(str(query_csv), format='csv')
    else:
        zquery = query.ZTFQuery()
        zquery.load_metadata(
            radec=[ra, dec],
            size=size / 7200,
            sql_query=f'seeing<{max_seeing}',
            kind='sci',
        )
        zquery.metatable.to_csv(query_csv)

    return zquery


def main(args):
    # ======= Init the file structure =======
    root = Path(os.path.expanduser(args.outdir))
    root.mkdir(exist_ok=True, parents=True)
    (root / 'query').mkdir(exist_ok=True, parents=True)

    # update local storage for ZTF
    os.environ['ZTFDATA'] = str(root)
    io.LOCALSOURCE = buildurl.LOCALSOURCE = str(root) + os.sep

    # ======= Query ZTF =======
    print('Querying ZTF...')
    zquery = load_query(args.ra, args.dec, args.size, args.max_seeing, root)
    print(f'{len(zquery.metatable)} exposures found', end='\n\n')

    # ======= Download data =======
    print('Downloading cutouts and PSF info...')
    start_time = time()
    for suffix in ['sciimg.fits', 'sciimgdaopsfcent.fits', 'sciimgdao.psf']:
        print(f'\nDownloading {suffix}...')

        kw = {}
        if suffix == 'sciimg.fits':
            kw = dict(cutouts=True, radec=[args.ra, args.dec], cutout_size=args.size)

        if suffix == 'sciimgdao.psf': # TODO: Should be a better way to deal with it.
            sys.stderr = open(os.devnull, 'w')  # suppress file testing warnings

        zquery.download_data(
            suffix,
            nprocess=args.nprocess,
            show_progress=True,
            overwrite=False,
            **kw
        )
        print('\nDone!')

    print(f'\nDone in [{time() - start_time:.2f}s]!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ra', type=float, default=83.633, help='RA of the cutout center in degrees')
    parser.add_argument('--dec', type=float, default=22.015, help='Dec of the cutout center in degrees')
    parser.add_argument('--size', type=int, default=64, help='Cutout half-size in arcseconds')
    parser.add_argument('--max-seeing', type=float, default=1.7, help='Maximum allowed seeing in arcseconds')
    parser.add_argument('--outdir', default='~/datasets/ZTF', help='Output directory for downloaded cutouts')
    parser.add_argument('--nprocess', type=int, default=12, help='Number of parallel processes to use for downloading')
    main(parser.parse_args())
