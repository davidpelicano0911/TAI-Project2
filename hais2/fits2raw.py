#!/usr/bin/env python3
"""
Convert a FITS image to raw 16-bit big-endian binary for use with HAIS.

Usage:
    python3 fits2raw.py <input.fits> <output.raw>

Output:
    Prints image dimensions for use with the compressor:
        ./compress <output.raw> <out.hais> <width> <height>
"""

import sys
import struct

def read_fits_header(f):
    """Parse FITS fixed-format header (80-char cards, 36 cards per 2880-byte block)."""
    header = {}
    while True:
        block = f.read(2880)
        if not block or len(block) < 2880:
            break
        for i in range(36):
            card = block[i*80:(i+1)*80].decode('ascii', errors='replace')
            key = card[:8].strip()
            if key == 'END':
                return header
            if '=' in card[8:10]:
                val_comment = card[10:]
                # strip inline comment
                val = val_comment.split('/')[0].strip()
                # strip FITS string quotes
                if val.startswith("'"):
                    val = val.strip("'").strip()
                header[key] = val
    return header

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.fits output.raw")
        sys.exit(1)

    fits_path = sys.argv[1]
    raw_path  = sys.argv[2]

    with open(fits_path, 'rb') as f:
        header = read_fits_header(f)

        bitpix = int(header.get('BITPIX', 16))
        naxis  = int(header.get('NAXIS',  2))
        naxis1 = int(header.get('NAXIS1', 0))   # width (columns)
        naxis2 = int(header.get('NAXIS2', 0))   # height (rows)
        bzero  = float(header.get('BZERO',  0.0))
        bscale = float(header.get('BSCALE', 1.0))

        print(f"BITPIX={bitpix}  NAXIS={naxis}  width={naxis1}  height={naxis2}")
        print(f"BZERO={bzero}  BSCALE={bscale}")

        if naxis < 2 or naxis1 == 0 or naxis2 == 0:
            print("ERROR: not a 2D image")
            sys.exit(1)

        npix = naxis1 * naxis2
        bytes_per_pix = abs(bitpix) // 8
        raw = f.read(npix * bytes_per_pix)
        if len(raw) < npix * bytes_per_pix:
            print(f"WARNING: expected {npix*bytes_per_pix} bytes, got {len(raw)}")

        # Unpack according to BITPIX (FITS is always big-endian)
        if bitpix == 16:
            fmt = f'>{npix}h'   # signed int16
        elif bitpix == -32:
            fmt = f'>{npix}f'   # float32
        elif bitpix == -64:
            fmt = f'>{npix}d'   # float64
        elif bitpix == 32:
            fmt = f'>{npix}i'   # signed int32
        elif bitpix == 8:
            fmt = f'>{npix}B'   # uint8
        else:
            print(f"ERROR: unsupported BITPIX={bitpix}")
            sys.exit(1)

        pixels = struct.unpack(fmt, raw[:npix * bytes_per_pix])

    # Apply BSCALE/BZERO and convert to uint16 big-endian
    # Common case: BITPIX=16, BZERO=32768 → unsigned 16-bit
    with open(raw_path, 'wb') as out:
        buf = bytearray(npix * 2)
        for i, p in enumerate(pixels):
            val = int(round(p * bscale + bzero))
            val = max(0, min(65535, val))
            buf[i*2]   = (val >> 8) & 0xFF
            buf[i*2+1] = val & 0xFF
        out.write(buf)

    print(f"\nWritten {naxis1}x{naxis2} uint16 big-endian → {raw_path}")
    print(f"\nTo compress:")
    print(f"  ./compress {raw_path} out.hais {naxis1} {naxis2}")
    print(f"To decompress:")
    print(f"  ./decompress out.hais recovered.raw")

if __name__ == '__main__':
    main()
