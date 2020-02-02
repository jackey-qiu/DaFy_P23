import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
from pyqtgraph.Qt import QtGui

#color_lib = {'C':(1,0,0,1),'O':(0,1,0,1),'Cu':(1,0,1,1)}
# color_lib = {'C':0xFFFFFF,'O':(0,1,0,1),'Cu':(1,0,1,1)}
def color_to_rgb(hex_str):
    rgb=[]
    for i in [0,2,4]:
        rgb.append(int(hex_str[i:(i+2)],16)/255.)
    rgb.append(1)
    return tuple(rgb)

color_lib = {
    "H": "FFFFFF",
    "HE": "D9FFFF",
    "LI": "CC80FF",
    "BE": "C2FF00",
    "B": "FFB5B5",
    "C": "909090",
    "N": "3050F8",
    "O": "FF0D0D",
    "F": "90E050",
    "NE": "B3E3F5",
    "NA": "AB5CF2",
    "MG": "8AFF00",
    "AL": "BFA6A6",
    "SI": "F0C8A0",
    "P": "FF8000",
    "S": "FFFF30",
    "CL": "1FF01F",
    "AR": "80D1E3",
    "K": "8F40D4",
    "CA": "3DFF00",
    "SC": "E6E6E6",
    "TI": "BFC2C7",
    "V": "A6A6AB",
    "CR": "8A99C7",
    "MN": "9C7AC7",
    "FE": "E06633",
    "CO": "F090A0",
    "NI": "50D050",
    "CU": "C88033",
    "ZN": "7D80B0",
    "GA": "C28F8F",
    "GE": "668F8F",
    "AS": "BD80E3",
    "SE": "FFA100",
    "BR": "A62929",
    "KR": "5CB8D1",
    "RB": "702EB0",
    "SR": "00FF00",
    "Y": "94FFFF",
    "ZR": "94E0E0",
    "NB": "73C2C9",
    "MO": "54B5B5",
    "TC": "3B9E9E",
    "RU": "248F8F",
    "RH": "0A7D8C",
    "PD": "006985",
    "AG": "C0C0C0",
    "CD": "FFD98F",
    "IN": "A67573",
    "SN": "668080",
    "SB": "9E63B5",
    "TE": "D47A00",
    "I": "940094",
    "XE": "429EB0",
    "CS": "57178F",
    "BA": "00C900",
    "LA": "70D4FF",
    "CE": "FFFFC7",
    "PR": "D9FFC7",
    "ND": "C7FFC7",
    "PM": "A3FFC7",
    "SM": "8FFFC7",
    "EU": "61FFC7",
    "GD": "45FFC7",
    "TB": "30FFC7",
    "DY": "1FFFC7",
    "HO": "00FF9C",
    "ER": "00E675",
    "TM": "00D452",
    "YB": "00BF38",
    "LU": "00AB24",
    "HF": "4DC2FF",
    "TA": "4DA6FF",
    "W": "2194D6",
    "RE": "267DAB",
    "OS": "266696",
    "IR": "175487",
    "PT": "D0D0E0",
    "AU": "FFD123",
    "HG": "B8B8D0",
    "TL": "A6544D",
    "PB": "575961",
    "BI": "9E4FB5",
    "PO": "AB5C00",
    "AT": "754F45",
    "RN": "428296",
    "FR": "420066",
    "RA": "007D00",
    "AC": "70ABFA",
    "TH": "00BAFF",
    "PA": "00A1FF",
    "U": "008FFF",
    "NP": "0080FF",
    "PU": "006BFF",
    "AM": "545CF2",
    "CM": "785CE3",
    "BK": "8A4FE3",
    "CF": "A136D4",
    "ES": "B31FD4",
    "FM": "B31FBA",
}

class GLViewWidget_cum(gl.GLViewWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

    def show_structure(self, xyz):
        # self.setCameraPosition(distance=55, azimuth=-90)
        self.setCameraPosition(distance=55,azimuth=0)
        self.setProjection()
        i=0
        if len(self.items)==0:
            xgrid = gl.GLGridItem()
            ygrid = gl.GLGridItem()
            zgrid = gl.GLGridItem()

            xgrid.setSize(5,5,5)
            ygrid.setSize(5,5,5)
            zgrid.setSize(5,5,5)

            xgrid.translate(0,0,3.615,True)
            zgrid.translate(0.,0,2*3.615,True)

            self.addItem(xgrid)
            self.addItem(ygrid)
            self.addItem(zgrid)
 
            # Rotate x and y grids to face the correct direction
            # xgrid.rotate(90, 0, 1, 0)
            # ygrid.rotate(90, 1, 0, 0)
 
            # Scale grids to the appropriate dimensions
            xgrid.scale(3.615, 3.615, 3.615)
            ygrid.scale(3.615, 3.615, 3.615)
            zgrid.scale(3.615, 3.615, 3.615)
            """
            g = gl.GLGridItem(size = QtGui.QVector3D(30,30,1))
            g.scale(1,1,1)
            # g.setSpacing(3.54,3.54,3.54)
            g.rotate(90,1,0,0)
            self.addItem(g)
            g1 = gl.GLGridItem(size = QtGui.QVector3D(30,30,1))
            g1.scale(1,1,1)
            g1.setSpacing(3.54,3.54,3.54)
            #g1.rotate(90,1,0,0)
            self.addItem(g1)
            """
            for each in xyz:
                e, x, y, z = each
                md = gl.MeshData.sphere(rows=10, cols=20)
                m1 = gl.GLMeshItem(meshdata=md, smooth=True, color=color_to_rgb(color_lib[e.upper()]), shader='shaded', glOptions='opaque')
                # print(dir(m1.metaObject()))
                m1.translate(x, y, z)
                m1.scale(0.5, 0.5, 0.5)
                self.addItem(m1)
        else:
            for each in xyz:
                _,x,y,z = xyz[i]
                #first item is grid net
                self.items[i+1+1+1].resetTransform()
                self.items[i+1+1+1].translate(x,y,z)
                self.items[i+1+1+1].scale(0.5, 0.5, 0.5)
                i += 1
        self.setProjection()

    def update_structure(self, xyz):
        for i in range(len(xyz)):
            _,x,y,z = xyz[i]
            #first item is grid net
            self.items[i+1+1+1].resetTransform()
            self.items[i+1+1+1].translate(x,y,z)
            self.items[i+1+1+1].scale(0.5, 0.5, 0.5)

    def setup_view(self):
        self.setCameraPosition(distance=15, azimuth=-90)
        g = gl.GLGridItem()
        g.scale(2,2,1)
        self.addItem(g)

        md = gl.MeshData.sphere(rows=10, cols=20)
        x = np.linspace(-8, 8, 6)

        m1 = gl.GLMeshItem(meshdata=md, smooth=True, color=(1, 0, 0), shader='balloon', glOptions='additive')
        m1.translate(x[0], 0, 0)
        m1.scale(1, 1, 2)
        self.addItem(m1)

        m2 = gl.GLMeshItem(meshdata=md, smooth=True, shader='normalColor', glOptions='opaque')
        m2.translate(x[1], 0, 0)
        m2.scale(1, 1, 2)
        self.addItem(m2)

        m3 = gl.GLMeshItem(meshdata=md, smooth=True, shader='viewNormalColor', glOptions='opaque')
        m3.translate(x[2], 0, 0)
        m3.scale(1, 1, 2)
        self.addItem(m3)

        m4 = gl.GLMeshItem(meshdata=md, smooth=True, shader='shaded', glOptions='opaque')
        m4.translate(x[3], 0, 0)
        m4.scale(1, 1, 2)
        self.addItem(m4)

        m5 = gl.GLMeshItem(meshdata=md, smooth=True, color=(1, 0, 0, 1), shader='edgeHilight', glOptions='opaque')
        m5.translate(x[4], 0, 0)
        m5.scale(1, 1, 2)
        self.addItem(m5)

        m6 = gl.GLMeshItem(meshdata=md, smooth=True, color=(1, 0, 0, 1), shader='heightColor', glOptions='opaque')
        m6.translate(x[5], 0, 0)
        m6.scale(1, 1, 2)
        self.addItem(m6)