# coding: latin-1
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************

_bibtexStr = """

@article{Bonvin2015,
title = {Fast and sensitive rigid-body fitting into cryo-EM density maps with PowerFit},
journal = {AIMS Biophysics},
volume = {2},
number = {20150273},
pages = {73-87},
year = {2015},
note = "",
issn = {2377-9098},
doi = {http://dx.doi.org/10.3934/biophy.2015.2.73},
url = {http://www.aimspress.com/biophysics/article/281.html},
author = {Gydo C.P.van Zundert and Alexandre M.J.J. Bonvin}
keywords = {cross correlation|exhaustive search|GPU acceleration|Fast Fourier Transform|optimized rotation sets|trimming|resampling|biomolecular complexes},
}

"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
