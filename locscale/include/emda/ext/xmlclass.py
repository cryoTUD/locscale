"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np

# XML Class to help with Balbes data
class Xml:
    def __init__(self):
        self.map1path = None
        self.map2path = None
        self.mapres1 = None
        self.mapres2 = None
        self.res_arr = None
        self.fsc = None
        self.outmap = None
        self.fsccutoff1 = None
        self.fsccutoff2 = None

    def write_xml(self):
        from lxml import etree
        import datetime

        now = datetime.datetime.now()
        date_string = now.strftime("%Y-%m-%d")
        time_string = now.strftime("%H:%M:%S")

        root = etree.Element("EMDA")
        datetim = etree.SubElement(root, "Session")
        datetim.set("Date", date_string)
        datetim.set("Time", time_string)
        input = etree.SubElement(root, "Inputs")
        input.set("map1", self.map1path)
        input.set("map2", self.map2path)
        data = etree.SubElement(root, "Data")
        resol = etree.SubElement(data, "Res")
        cc = etree.SubElement(data, "Fsc")
        resol.text = str(self.res_arr)
        cc.text = str(self.fsc)
        mapres = etree.SubElement(root, "Resolution")
        cutoff1 = etree.SubElement(mapres, "cutoff1")
        cutoff2 = etree.SubElement(mapres, "cutoff2")
        cutoff1.set("fsc", str(self.fsccutoff1))
        cutoff1.set("resol", str(self.mapres1))
        cutoff2.set("fsc", str(self.fsccutoff2))
        cutoff2.set("resol", str(self.mapres2))
        outmap = etree.SubElement(root, "Mapout")
        outmap.set("Outmapfile", self.outmap)
        #
        strng = etree.tostring(root, pretty_print=True)
        with open("EMDA.xml", "wb") as file:
            file.write(strng)


if __name__ == "__main__":
    xml = Xml()
    xml.map1path = "path/to/map1"
    xml.map2path = "path/to/map2"
    xml.res_arr = np.array([42.00, 30.00, 23.33, 19.09, 16.15], dtype=float)
    xml.fsc = np.array([1.00, 1.00, 1.00, 1.00, 1.00], dtype=float)

    xml.write_xml()
