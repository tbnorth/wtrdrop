"""
wtrdrop.py - find landscape elevation about stream

for each stream cell
  generate circular region around it radius r
  generate binary layer of cells < h higher than stream cell in circle
  floodfill from stream cell in circle
  flooded cells are riparian for h, r for that cell

  there will be a non-linear jump in riparian extent if a stream
  "moved" as the circle slips round the end of a ridge

Terry N. Brown terrynbrown@gmail.com Mon Mar 18 12:24:13 EDT 2019
"""

import argparse

import numpy as np

from osgeo import gdal
from scipy.ndimage import label


def make_parser():

    parser = argparse.ArgumentParser(
        description="""Map riparian areas""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--radius",
        type=float,
        default=300.0,
        help="Radius (m) for stream cells to affect",
    )
    parser.add_argument(
        "--levels",
        type=float,
        default=[2.0],
        nargs='+',
        metavar='LEVEL',
        help="Height(s) (m) for riparian threshold",
    )
    parser.add_argument('elevation', type=str, help="Path to elevation grid")
    parser.add_argument('stream', type=str, help="Path to stream grid")

    return parser


def get_options(args=None):
    """
    get_options - use argparse to parse args, and return a
    argparse.Namespace, possibly with some changes / expansions /
    validatations.

    Client code should call this method with args as per sys.argv[1:],
    rather than calling make_parser() directly.

    Args:
        args ([str]): arguments to parse

    Returns:
        argparse.Namespace: options with modifications / validations
    """
    opt = make_parser().parse_args(args)

    # modifications / validations go here

    return opt


def make_kernel(opt, elev_ds):
    _, w, _, _, _, h = elev_ds.GetGeoTransform()
    szx = int(opt.radius / w)
    szy = int(opt.radius / abs(h))
    if szx % 2 == 0:
        szx += 1
    if szy % 2 == 0:
        szy += 1
    return np.ones((szx, szy))


def dbg(*s):
    return
    print(s)


def rip(opt, elev, ripar, kernel, row, col, level):
    kr = (kernel.shape[0] - 1) // 2
    kc = (kernel.shape[1] - 1) // 2
    r0 = row - kr
    c0 = col - kc
    dbg("-----------------", level, elev[row, col], row, col)
    window = elev[r0 : r0 + kernel.shape[0], c0 : c0 + kernel.shape[1]]
    dbg(window.astype(int))
    window = window - elev[row, col] < level
    dbg(window)
    blobs, blobs_n = label(window)
    dbg(blobs)
    stream_label = blobs[kr, kc]
    dbg(kr, kc, stream_label)
    blobs = blobs == stream_label
    dbg(blobs)
    ripar[r0 : r0 + 2 * kr + 1, c0 : c0 + 2 * kc + 1] += blobs * kernel


def do_cell(opt, elev, ripar, kernel, row, col):
    for level_i, level in enumerate(opt.levels):
        rip(opt, elev, ripar[level_i], kernel, row, col, level)


def main():
    opt = get_options()

    elev_ds = gdal.Open(opt.elevation)
    stream_ds = gdal.Open(opt.stream)

    elev = elev_ds.GetRasterBand(1).ReadAsArray()
    stream = stream_ds.GetRasterBand(1).ReadAsArray()

    ripar = [np.zeros_like(stream, dtype=float) for i in opt.levels]

    kernel = make_kernel(opt, elev_ds)

    done = 0
    todo = (stream == 1).sum()
    for row in range(stream.shape[0]):
        for col in range(stream.shape[1]):
            if stream[row, col] == 1:
                do_cell(opt, elev, ripar, kernel, row, col)
                done += 1
                if done % 1000 == 0:
                    print("%s/%s" % (done, todo))

    for level_i, level in enumerate(opt.levels):
        filename = "rip%s.tif" % level
        out = elev_ds.GetDriver().CreateCopy(filename, elev_ds)
        out.GetRasterBand(1).WriteArray(ripar[level_i])
        del out


if __name__ == "__main__":
    main()
