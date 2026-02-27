#!/usr/bin/env python

from pathlib import Path

from pyvips import Image

dir = Path("..") / "data" / "VisiumHD_CRC"
filename = dir / "Visium_HD_Human_Colon_Cancer_tissue_image.btf"
new_filename = dir / "Visium_HD_Human_Colon_Cancer_tissue_image_pyramid.btf"

Image.tiffload(filename).tiffsave(new_filename, pyramid=True)
